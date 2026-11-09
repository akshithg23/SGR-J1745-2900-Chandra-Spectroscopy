# SPA660 Project 8: SGR J1745−2900

This project analyses Chandra ACIS-S ObsID 14946 and compares it with Rea et al. (2013), ApJL 775 L34. Download the observation from the [Chandra Data Archive](https://cxc.cfa.harvard.edu/cda/) using ObsID 14946. The [reference paper](https://doi.org/10.1088/2041-8205/775/2/L34) describes the aperture and comparison fits in Table 1 and §2.3.

## Verify the reported fits

Download the [validated results archive (`SPA660_Results_verified.zip`)](https://drive.google.com/file/d/1ZkVJmYkV6vciIvEKzS94Lnvsm1iXcyZA/view?usp=drive_link) using the Google account authorized by the submitting group. The archive contains `products_verified/`, including the source and background event subsets, PI spectra, RMF, both ARFs, REEF fractions, saved XSPEC sessions and log, fit summary, and plot data. It does not contain the full reprocessed event file, aspect histogram or original pileup and aperture-check logs. SHA-256 of the audited original archive:

```text
ee06fccbb326d79a7544eb3c4151c0d583b804cd495b9c7e1b8fda4ba139aa83
```

The read-only check requires Python with Astropy and NumPy. After extracting the archive into a separate workspace, run:

```bash
python verify_products.py /path/to/workspace/products_verified
```

This command was run against the saved archive products. It checks FITS files, exposure, counts, response associations, ARF multiplication and saved fit artifacts. A live refit additionally requires HEASoft 6.36 / XSPEC 12.15.1 / PyXspec 2.1.5 with Astropy, NumPy and Matplotlib. The saved sessions reference their spectra and responses by relative path. Use the copy-and-fit instructions in [WORKFLOW.md](WORKFLOW.md) so the archived fits remain available for comparison.

## Reduce and inspect the observation

Use a Linux environment with CIAO 4.17.0, CALDB 4.12.4, `ciao-contrib`, and SAOImage DS9. The saved event files record CIAO/CALDB in their processing history and retain an older archive CALDB keyword. WSL2 is an option on Windows. See [WORKFLOW.md](WORKFLOW.md) for installation, exact commands, the standard route and the separately documented physical-coordinate replay. CIAO, DS9 and XSPEC replay commands require student execution.

- `reduce.sh`, `specextract_compat.py`, and the two region files cover standard reduction. `reduce.sh --preprocess-only` prepares the image and pileup diagnostic before trying extraction.
- `reextract_spa660.py` selects physical SKY events and checks the PI spectra. `make_reef_arf.py` documents the manual REEF correction as a replay step. `repair_spa660_exposure.py` is a legacy utility for the earlier `products_fast` workflow.
- `fit_xspec.py` fits the spectra and exports model summaries, curves and figures. `polish_ciao.sh`, `polish_xspec.py`, and `project_chart.sh` contain supplementary checks.
- `build_final_report.py` and `build_report.py` are historical report builders. They require ReportLab and Windows Calibri fonts. The earlier draft builder uses local draft inputs. The final builder contains prose for the archived run. Treat new fit products as evidence to review before editing the report.

The fitted 27 groups in the archive cover 2.263–6.687 keV. The reported 1–10 keV fluxes are model integrals and include extrapolation outside those groups. The saved ARF uses the REEF annular fraction once. A standard CIAO extraction with `arfcorr` is a separate calibration implementation and should have its own results directory.
