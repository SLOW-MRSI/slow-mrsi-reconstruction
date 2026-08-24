# SLOW-MRSI reconstruction

MATLAB scripts for reconstruction, draft metabolite-map generation, and spectrIm export of SLOW-MRSI time-resolved EPSI data.

## Workflow

1. Put the required Siemens raw-data `.dat`/Twix files in the working folder.
2. Open MATLAB in this repository folder.
3. Edit the user settings in `Run_01_recon.m` (especially `cfg.measID`).
4. Run the scripts in order:

```matlab
Run_01_recon
Run_02_map
Run_03_export
```

The scripts create `processedData/` and `data4spectrIm/` locally. These generated folders and raw MRI data are excluded from version control by `.gitignore`.

## Notes

- `Run_02_map.m` produces draft peak-integration maps for quality control and region-of-interest selection; it is not a replacement for quantitative spectral fitting.
- The code expects MATLAB with the toolboxes required by the functions used in the scripts. Parallel execution can be disabled in `Run_01_recon.m` if needed.

## Author

Dr. Guodong Weng, University of Bern.
