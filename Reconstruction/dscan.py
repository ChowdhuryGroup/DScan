import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from tqdm import tqdm
import csv

_EPS = 1e-30  # small guard against division by zero


class DScanLoader:
    """Load raw DScan TSV output from dscanner.py."""

    def __init__(self, filepath, bg_filepath=None):
        data, self.meta = self._parse(filepath)
        self.wavelengths = data[0, 1:]   # nm,       (Nlambda,)
        self.fundamental = data[1, 1:]   #           (Nlambda,)
        self.positions   = data[2:, 0]   # mm,       (Nz,)
        self.dscan       = data[2:, 1:]  # (Nz, Nlambda)
        self.background  = None
        if bg_filepath:
            bg, _ = self._parse(bg_filepath)
            self.background = bg[1, :]     # (Nlambda,)

    @staticmethod
    def _parse(path):
        meta, rows = {}, []
        with open(path) as f:
            for line in csv.reader(f, delimiter='\t'):
                if not line:
                    continue
                if line[0].startswith('#'):
                    if len(line) >= 2:
                        meta[line[0][2:].strip()] = line[1].strip()
                else:
                    rows.append([float(v) for v in line])
        return np.array(rows), meta


class TraceProcessor:
    """Background subtract, resample to uniform freq grid, truncate around SH."""

    def __init__(self, loader, Nsamps=256):
        wl   = loader.wavelengths
        fund = loader.fundamental.copy()
        raw  = loader.dscan.copy()

        if loader.background is not None:
            fund -= loader.background
            raw  -= loader.background[None, :]

        # Wavelength → uniform frequency (PHz); Jacobian corrects spectral density
        raw_f = 300.0 / wl          # decreasing order, so reverse for interp
        df    = np.mean(np.abs(np.diff(raw_f)))
        uni_f = np.arange(raw_f.min(), raw_f.max(), df)
        jac   = 300.0 / uni_f ** 2  # |dlambda/df|

        def _rs(y):
            fi = interp1d(raw_f[::-1], y[::-1], kind='cubic',
                          bounds_error=False, fill_value=0.0)
            return jac * fi(uni_f)

        rs_fund  = _rs(fund)
        rs_trace = np.stack([_rs(raw[i]) for i in range(len(loader.positions))])

        # Central fundamental freq (weighted mean over inner 90%)
        lo, hi  = int(0.05 * len(uni_f)), int(0.95 * len(uni_f))
        f0_fund = (np.sum(uni_f[lo:hi] * rs_fund[lo:hi])
                   / np.sum(rs_fund[lo:hi]))
        f0_SH   = 2.0 * f0_fund

        # Truncate Nsamps around SH / fund centres; clamp to valid array range
        N  = Nsamps // 2
        iS = int(np.clip(np.argmin(np.abs(uni_f - f0_SH)),  N, len(uni_f) - N))
        iF = int(np.clip(np.argmin(np.abs(uni_f - f0_fund)), N, len(uni_f) - N))

        self.freqs      = uni_f[iS - N : iS + N] - f0_SH
        self.dscan_proc = rs_trace[:, iS - N : iS + N]
        self.fund_spec  = rs_fund [iF - N : iF + N]
        self.positions  = loader.positions
        self.f0SH, self.df = f0_SH, df

    def plot(self):
        _, ax = plt.subplots(1, 2, figsize=(12, 4))
        ax[0].pcolormesh(self.freqs, self.positions, self.dscan_proc, cmap='jet')
        ax[0].set(xlabel='SH freq offset (PHz)', ylabel='Position (mm)',
                  title='Processed D-scan')
        ax[1].plot(self.freqs, self.fund_spec)
        ax[1].set(xlabel='Freq offset (PHz)', title='Fundamental spectrum')
        plt.tight_layout(); plt.show()

    def save(self, path):
        dpos   = self.positions[1] - self.positions[0]
        N      = len(self.fund_spec)
        # Header row padded to match Nsamps columns
        header = np.zeros((1, N))
        header[0, :4] = [self.f0SH, self.df, self.positions[0], dpos]
        np.savetxt(path, np.vstack([header, self.fund_spec, self.dscan_proc]),
                   delimiter='\t')

    @classmethod
    def from_file(cls, path):
        """Load a previously saved processed file."""
        obj = cls.__new__(cls)
        dat = np.loadtxt(path, delimiter='\t')
        obj.f0SH, obj.df, Pi, dP = dat[0, :4]
        obj.fund_spec  = dat[1]
        obj.dscan_proc = dat[2:]
        N = len(obj.fund_spec)
        obj.freqs     = (np.arange(N) - N // 2) * obj.df
        obj.positions = Pi + np.arange(len(obj.dscan_proc)) * dP
        return obj


class DERetrieval:
    """Differential evolution D-scan retrieval for a grating compressor.

    Based on DE/BoR/1/bin with multi-resolution spectral phase (Escoto 2017).
    GVD in fs²/mm, TOD in fs³/mm, frequencies in PHz.
    """

    def __init__(self, proc, GVD, TOD, NP=80, N1=16,
                 iters_per_res=None, Cr=0.65):
        self.GVD = GVD;  self.TOD = TOD
        self.NP  = NP;   self.N1  = N1;  self.Cr = Cr
        self.iters_per_res = iters_per_res or [50, 100, 200, 200, 250]
        self._setup(proc)

    def _setup(self, proc):
        freqs  = proc.freqs
        Nsamps = len(freqs)

        # Recenter z at position-marginal peak
        z_ci = np.argmax(proc.dscan_proc.sum(axis=1))
        z    = proc.positions - proc.positions[z_ci]

        meas = proc.dscan_proc / proc.dscan_proc.max()
        meas = np.roll(meas, Nsamps // 2 - 1, axis=1)

        w = 2 * np.pi * freqs
        self.Espec      = np.sign(proc.fund_spec) * np.sqrt(np.abs(proc.fund_spec))
        self.dscan_meas = meas
        self.phase_mat  = np.exp(1j * z[None, :] * (
            0.5 * self.GVD * w[:, None] ** 2
            + (1 / 6) * self.TOD * w[:, None] ** 3))
        self.freqs     = freqs
        self.z         = z
        self.positions = proc.positions
        self.Nsamps    = Nsamps
        self.t         = (np.arange(Nsamps) - Nsamps // 2) / (2 * freqs.max())  # fs

    # ------------------------------------------------------------------
    def _g_batch(self, phases):
        """phases (Nsamps, K) → G values (K,).

        Vectorised over K candidates simultaneously.
        """
        E    = self.Espec[:, None] * np.exp(1j * phases)
        Emat = E[:, :, None] * self.phase_mat[:, None, :]   # (Nsamps, K, Nz)
        ft   = np.fft.fft(Emat, axis=0)
        sim  = np.abs(np.fft.ifft(ft ** 2, axis=0)) ** 2   # (Nsamps, K, Nz)
        sim  = sim.transpose(1, 2, 0)                        # (K, Nz, Nsamps)
        meas = self.dscan_meas                               # (Nz, Nsamps)
        mu   = np.sum(meas * sim, axis=1) / (np.sum(sim ** 2, axis=1) + _EPS)
        diff = meas - mu[:, None, :] * sim
        return np.sqrt(np.mean(diff ** 2, axis=(1, 2)))      # (K,)

    # ------------------------------------------------------------------
    def run(self):
        NP, N1, Nsamps, freqs = self.NP, self.N1, self.Nsamps, self.freqs
        NrInc = int(np.log2(N1))

        ind0 = np.arange(0, Nsamps, N1)
        X    = np.random.uniform(-np.pi, np.pi, (NP, len(ind0)))
        Xrs  = interp1d(freqs[ind0], X.T, kind='cubic', axis=0,
                        bounds_error=False, fill_value=0)(freqs)  # (Nsamps, NP)

        self.G_history = []

        for q in range(NrInc + 1):
            step = N1 // (2 ** q)
            ind  = np.arange(0, Nsamps, step)
            nc   = len(ind)
            X    = Xrs[ind].T                               # (NP, nc)

            for _ in tqdm(range(self.iters_per_res[q]),
                          desc=f"Res {q + 1}/{NrInc + 1}"):
                # --- Mutation: DE/BoR/1 ---
                V = np.empty_like(X)
                for i in range(NP):
                    oth  = np.delete(np.arange(NP), i)
                    abc  = np.random.choice(oth, 3, replace=False)
                    G3   = self._g_batch(Xrs[:, abc])
                    s    = np.argsort(G3)
                    V[i] = (X[abc[s[0]]]
                            + np.random.rand() * (X[abc[s[1]]] - X[abc[s[2]]]))

                # --- Crossover ---
                U    = X.copy()
                mask = np.random.rand(NP, nc) <= self.Cr
                mask[np.arange(NP), np.random.randint(0, nc, NP)] = True
                U[mask] = V[mask]
                Urs  = interp1d(freqs[ind], U.T, kind='cubic', axis=0,
                                bounds_error=False, fill_value=0)(freqs)

                # --- Selection ---
                GX = self._g_batch(Xrs)
                GU = self._g_batch(Urs)
                better = GU <= GX
                X[better]      = U[better]
                Xrs[:, better] = Urs[:, better]
                self.G_history.append(min(GX.min(), GU.min()))

        all_G = self._g_batch(Xrs)
        best  = np.argmin(all_G)
        print(f"Final G = {all_G[best]:.6f}")
        self._compute_results(Xrs, best)

    # ------------------------------------------------------------------
    def _compute_results(self, Xrs, best):
        freqs, Nsamps, z = self.freqs, self.Nsamps, self.z

        E    = self.Espec * np.exp(1j * Xrs[:, best])
        Emat = E[:, None] * self.phase_mat          # (Nsamps, Nz)
        ft   = np.fft.fft(Emat, axis=0)

        # Retrieved D-scan trace
        dsim = np.abs(np.fft.ifft(ft ** 2, axis=0)) ** 2
        dsim = np.roll(dsim.T, Nsamps // 2, axis=1)
        self.dscan_retrieved = dsim / dsim.max()

        # Temporal intensity; find peak z and t
        Imat = np.roll(np.abs(ft).T ** 2, Nsamps // 2, axis=1)  # (Nz, Nsamps)
        zi, ti = np.unravel_index(np.argmax(Imat), Imat.shape)
        t0 = self.t[ti]

        # Remove linear phase to centre pulse at t = 0
        Emat *= np.exp(-1j * 2 * np.pi * freqs[:, None] * t0)
        Emat *= np.abs(self.Espec).max() / np.abs(Emat).max()
        ft    = np.fft.fft(Emat, axis=0)
        Imat  = np.roll(np.abs(ft).T ** 2, Nsamps // 2, axis=1)

        bwl     = np.abs(np.fft.fftshift(np.fft.fft(self.Espec))) ** 2
        bwl_max = bwl.max()
        self.bwl_intensity  = bwl / bwl_max
        self.best_intensity = Imat[zi] / bwl_max
        self.best_z_idx     = zi

        # Spectral phase at best-compression z position
        phase = np.unwrap(Xrs[:, best])
        self.best_spec_phase = (
            phase + 2 * np.pi * freqs * (-t0)
            + z[zi] * (0.5 * self.GVD * (2 * np.pi * freqs) ** 2
                       + (1 / 6) * self.TOD * (2 * np.pi * freqs) ** 3)
        )
        self.best_spec_phase -= self.best_spec_phase[Nsamps // 2] - np.pi
        self.fund_spec_norm   = np.abs(self.Espec) ** 2 / np.abs(self.Espec).max() ** 2

        print(f"Best compression at z = {z[zi]:.3f} mm  "
              f"(stage pos = {self.positions[zi]:.3f} mm)")

    # ------------------------------------------------------------------
    def plot_results(self):
        freqs = self.freqs
        z     = self.z
        meas  = np.roll(self.dscan_meas, -(self.Nsamps // 2 - 1), axis=1)

        fig, axes = plt.subplots(2, 3, figsize=(15, 8))

        axes[0, 0].semilogy(self.G_history, '.b', ms=2)
        axes[0, 0].set(title='G history', xlabel='Iteration')

        axes[0, 1].pcolormesh(freqs, z, meas, cmap='jet')
        axes[0, 1].set(title='Measured trace',
                       xlabel='SH freq offset (PHz)', ylabel='z (mm)')

        axes[0, 2].pcolormesh(freqs, z, self.dscan_retrieved, cmap='jet')
        axes[0, 2].set(title='Retrieved trace',
                       xlabel='SH freq offset (PHz)', ylabel='z (mm)')

        axes[1, 0].plot(self.t, self.bwl_intensity, 'b', label='BW-limited')
        axes[1, 0].plot(self.t, self.best_intensity, 'r', label='Retrieved')
        axes[1, 0].set(title='Temporal pulse', xlabel='Time (fs)')
        axes[1, 0].legend()

        ax_ph = axes[1, 1].twinx()
        axes[1, 1].plot(freqs, self.best_spec_phase / (2 * np.pi), 'r', label='Phase')
        ax_ph.plot(freqs, self.fund_spec_norm, 'b', label='Spectrum')
        axes[1, 1].set(title='Spectral phase at best compression',
                       xlabel='Freq offset (PHz)', ylabel='Phase (cycles)')
        ax_ph.set_ylabel('Normalised intensity')
        axes[1, 1].legend(loc='upper left'); ax_ph.legend(loc='upper right')

        axes[1, 2].axis('off')
        plt.tight_layout()
        plt.show()
