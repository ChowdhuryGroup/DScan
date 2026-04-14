# DScan Post-Processing — Python

Pulse retrieval from D-scan measurements using differential evolution (Escoto 2017, DE/BoR/1/bin with multi-resolution spectral phase).
Designed for a **grating compressor** scan.

## Files

| File | Purpose |
|------|---------|
| `dscan.py` | Core classes |
| `run.py` | Example workflow |

## Classes

**`DScanLoader(filepath, bg_filepath=None)`**  
Reads the `.tsv` output of `Collection/dscanner.py`.  Parses metadata comments, wavelengths, fundamental spectrum, stage positions, and the D-scan matrix.  Pass the accompanying `_background.tsv` file to enable background subtraction.

**`TraceProcessor(loader, Nsamps=256)`**  
Background-subtracts, resamples to a uniform frequency grid (with Jacobian), and truncates `Nsamps` points around the second-harmonic centre frequency.  Call `.plot()` to inspect, `.save(path)` to cache, or load a cached file with `TraceProcessor.from_file(path)`.

**`DERetrieval(proc, GVD, TOD, NP=80, N1=16, iters_per_res=None, Cr=0.65)`**  
Runs multi-resolution differential evolution to retrieve the spectral phase.  Call `.run()` then `.plot_results()`.

| Parameter | Default | Notes |
|-----------|---------|-------|
| `GVD` | — | fs²/mm — from `GratingCompressorDispersionCalculator.m` |
| `TOD` | — | fs³/mm |
| `NP`  | 80 | Population size |
| `N1`  | 16 = 2⁴ | Coarsest phase sampling; sets number of resolution levels |
| `iters_per_res` | `[50,100,200,200,250]` | Iterations per resolution level |
| `Cr`  | 0.65 | Crossover rate |

## Dependencies

```
numpy  scipy  matplotlib  tqdm
```
