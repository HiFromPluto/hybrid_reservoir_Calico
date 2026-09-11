BSim writes one subfolder per run, named by the config `output.dir`.

Each run contains:
  voxels.csv           16 samples per window, semicolon-delimited
  window_summary.csv   one row per window
  results.csv          per-sample summary
  matlab_meta.txt      reshape recipe for MATLAB
  feature_contract.txt official vs excluded channels
  run_status.txt       row counts

This folder starts empty. Do not commit large CSVs.
