from dscan import DScanLoader, TraceProcessor, DERetrieval

# --- 1. Load raw data from dscanner.py output ---------------------------------
loader = DScanLoader(
    'data/scan.tsv',
    bg_filepath='data/scan_background.tsv',   # optional but recommended
)

# --- 2. Process trace ----------------------------------------------------------
proc = TraceProcessor(loader, Nsamps=256)
proc.plot()
# proc.save('data/scan_processed.tsv')        # save for later reuse

# Load previously saved processed file instead of running the two steps above:
# proc = TraceProcessor.from_file('data/scan_processed.tsv')

# --- 3. Retrieve spectral phase -----------------------------------------------
# GVD and TOD from GratingCompressorDispersionCalculator.m (fs²/mm, fs³/mm)
ret = DERetrieval(proc, GVD=-4154.4, TOD=2582.8, NP=80)
ret.run()
ret.plot_results()
