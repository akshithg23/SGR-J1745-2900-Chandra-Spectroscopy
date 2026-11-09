# SPA660 Project 8 — Chandra spectroscopy of SGR J1745-2900

This repository collects the project scripts for Chandra ObsID 14946. The public observation is available from the [Chandra Data Archive](https://cxc.cfa.harvard.edu/cdaftp/byobsid/6/14946/primary/).

## Files

- `reduce.sh`, `specextract_compat.py`, `source.reg`, and `background.reg` cover CIAO reprocessing and spectrum extraction.
- `inspect_events.py`, `repair_spa660_exposure.py`, and `reextract_spa660.py` check the event selections and support the alternate extraction steps used on SciServer.
- `fit_xspec.py` fits the extracted spectra with PyXspec and exports the fit summaries, model curves, and plots.
- `polish_ciao.sh`, `polish_xspec.py`, and `project_chart.sh` contain supplementary aperture, pileup, and spectral checks.
- `build_final_report.py` builds the report from the verified results archive. It shares layout functions with `build_report.py`, the earlier draft builder. Student identifiers were removed from the repository copies.

## Verify the reported fits

Use the HEASoft environment with PyXspec, Astropy, NumPy, and Matplotlib. Extract `SPA660_Results_verified.zip` from the project files into the repository root so `products_verified/` contains the spectra and response files. This archive is needed for numerical verification of the saved sessions. The `.xcm` files refer to their companion spectra and responses by relative path.

```bash
python fit_xspec.py --products products_verified --rea-crosscheck
```

See [the reduction notes](WORKFLOW.md) for the CIAO commands. The reduction notes describe the standard extraction route and include draft-stage prompts. The verified fits used the corrected ARF and extracted spectra in `SPA660_Results_verified.zip`.

To build the report from the saved results, use Python with ReportLab and the Windows Calibri fonts.

```text
python build_final_report.py SPA660_Results_verified.zip SPA660_Report.pdf --code-url https://github.com/akshithg23/SGR-J1745-2900-Chandra-Spectroscopy
```
