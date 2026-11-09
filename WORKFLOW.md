# SPA660 Project 8 analysis instructions

The Level 2 archive event file and the archive figures are included. The response-based results require CIAO with CALDB and HEASoft with PyXspec. Use a lab Linux environment or a [SciServer HEASARC container](https://heasarc.gsfc.nasa.gov/docs/sciserver/) with the HEASARC data volume mounted. Unzip this bundle in a writable directory and open a Terminal there.

## Task A · Reprocess and inspect the observation

```bash
conda activate ciao
bash reduce.sh
```

The script uses a writable copy of the complete observation from `/FTP/chandra/data/byobsid/6/14946` when present. Otherwise, these are the download and reprocessing commands. The current CIAO contributed script is named `download_chandra_obsid`.

On a rerun, the script reuses the completed event file recorded in `products/event_path.txt` if it is still present. On SciServer, the extraction launcher detects the Python 3.11 syntax error in `ciao_contrib/_tools/merging.py` and runs with a temporary quote-only fix when needed.

```bash
download_chandra_obsid 14946
mkdir -p products
punlearn chandra_repro
chandra_repro indir=14946 outdir=products/repro verbose=2 clobber=yes
```

`reduce.sh` finds the sole reprocessed event file, saves its relative path in `products/event_path.txt`, and generates the following image. The bins are half an ACIS sky pixel. Inspect the full image and the extraction regions in DS9.

```bash
evt=$(cat products/event_path.txt)
dmcopy "$evt[ccd_id=7,energy=1000:8000][bin x=4057.5:4097.5:0.5,y=4089.0:4129.0:0.5]" \
    products/field_1_8keV.fits clobber=yes
ds9 products/field_1_8keV.fits -scale log -zoom to fit \
    -regions load source.reg -regions load background.reg
```

In DS9, set the coordinate readout to **FK5** in **degrees** and check that the bright peak is centred within the annulus at RA 266.4173708°, Dec −29.0082889°. Cross-check Sgr A* to the northwest at approximately RA 266.416837°, Dec −29.007810°. Compare the image peak with the adopted position before adjusting any aperture. DS9 region properties show the centre and both radii. Review the same overlays in JS9 if running through SciServer without a DS9 display. The image axes in the report are sky offsets from the adopted position.

To assess pileup, the script creates a **native-pixel, full-band CCD 7 image** from the standard reprocessed events and runs:

```bash
dmcopy "$evt[ccd_id=7][bin x=3072.5:5120.5:1,y=3072.5:5120.5:1]" \
    products/ccd7_native.fits clobber=yes
pileup_map products/ccd7_native.fits products/pileup_map.fits clobber=yes
punlearn dmcoords
dmcoords "$evt" asol=none option=cel celfmt=deg ra=266.4173708 dec=-29.0082889
sx=$(pget dmcoords x)
sy=$(pget dmcoords y)
dmstat "products/pileup_map.fits[sky=circle($sx,$sy,1.52439)]" centroid=no
dmstat "products/pileup_map.fits[sky=annulus($sx,$sy,1.52439,4.06504)]" centroid=no
```

Inspect the map with the source annulus overlaid and read the two maximum values in `products/pileup_summary.txt`. The summary uses physical sky pixels and runs after spectral extraction. If the summary reports an error, inspect the map separately before interpreting pileup. The map units are counts per frame. CIAO gives rough guide values of 0.02, 0.1 and 0.2 counts/frame for approximately 1%, 5% and 10% pileup. Record the core and annular maxima. Rea et al. estimated about 5–15% pileup across their ACIS observations and chose a 0.75″ inner cut. The 2.0″ outer radius limits diffuse background and contamination from Sgr A*, about 2.4″ away. State whether the map supports excluding the core, and treat any appreciable annular signal as possible residual pileup.

## Task B · Extract and check the response

`source.reg` and `background.reg` contain the FK5 annuli and the six circular background exclusions. The script runs:

```bash
punlearn specextract
python specextract_compat.py infile="$evt[ccd_id=7][sky=region(source.reg)]" \
    bkgfile="$evt[ccd_id=7][sky=region(background.reg)]" outroot=products/magnetar \
    refcoord="266.4173708 -29.0082889" weight=no correctpsf=yes binarfcorr=0.25 \
    bkgresp=no grouptype=NUM_CTS binspec=50 bkg_grouptype=NONE clobber=yes verbose=2
```

Use `magnetar_grp.pi` with `magnetar_bkg.pi`, `magnetar.rmf` and **`magnetar.corr.arf`**. CIAO uses `dmgroup` through `specextract` for a minimum of 50 source counts per group. XSPEC ignores bad-quality bins. Check the usable groups at both band edges and record the number of noticed bins. Inspect the event header for the actual reprocessing CALDB version, and inspect `BACKSCAL` and `EXPOSURE` for the PHA files. The archive geometric ratio of 0.069 is a cross-check for the CIAO extraction.

`correctpsf=yes` invokes `arfcorr`. The energy-dependent annular fraction is the encircled PSF fraction within 2.0″ minus that within 0.75″ for a concentric annulus. The corrected effective area is `SPECRESP_corrected = SPECRESP_uncorrected × PSF_FRAC`. The multiplicative aperture factor is `1 / PSF_FRAC`; the corrected ARF already accounts for it during the fit. The fitting script checks the corrected ARF against the uncorrected ARF and writes `psf_fraction.csv` and the 4 keV fraction and factor to `fit_results.json`. Record **[f(4 keV) = ?; 1/f(4 keV) = ?; 1–10 keV fraction range = ?]**. The narrow annulus warrants checking the fraction with finer PSF sampling or a simulated PSF before quoting precise fluxes. A separate ARF is needed if the extraction geometry changes.

The included archive count-rate spectrum uses 0.25 keV bins and peaks in the 3.50–3.75 keV bin, with the strongest signal around 3–5 keV and a steep drop above 6 keV. The spectrum remains in detector space. The reprocessed **PI** spectrum plot is generated as `products/grouped_spectrum.png`. To regenerate the field and the fixed-width count-rate figure from the reprocessed events, run this in the HEASoft environment:

```bash
python inspect_events.py "$(cat products/event_path.txt)" --output products/event_check
```

## Task C · Fit the four cases in XSPEC

```bash
conda activate heasoft
python fit_xspec.py --rea-crosscheck
```

Use the Python interpreter that imports `xspec`, `astropy`, `numpy` and `matplotlib`. The script runs blackbody and power-law fits with a free column and with `nH=10.3` fixed, in XSPEC units of 10²² cm⁻². It sets `abund angr`, `xsect bcmc` and `TBABSVERSION 2`. **TBabs uses its own cross-sections.** The optional flag fits each continuum again with `phabs`, which honors `xsect bcmc` and matches the absorption prescription of Rea et al. more closely. Use the free-column fits for the measured columns. The frozen fits are a sensitivity check, with the same blackbody reference column imposed on both continua.

To inspect a fit interactively, launch XSPEC inside `products`. For the blackbody free-column fit the commands are:

```text
data 1:1 magnetar_grp.pi
backgrnd 1 magnetar_bkg.pi
response 1 magnetar.rmf
arf 1 magnetar.corr.arf
abund angr
xsect bcmc
xset TBABSVERSION 2
setplot energy
ignore bad
ignore **-1.0 10.0-**
statistic chi
weight standard
model tbabs*bbody & 10.3 & 0.9 & 2.0e-4
fit
error 1.0 1-3
show fit
cpd /xs
plot data delchi
plot eufspec
```

For the frozen run, use `newpar 1 10.3`, `freeze 1`, `fit` and `error 1.0 2-3`. For the power law, replace the model line with `model tbabs*powerlaw & 17.3 & 4.2 & 0.1` and run the free and fixed sequences again. These starting parameters initialize the fit. They are separate from the output measurements. `error 1.0` profiles the other free parameters at Δχ² = 1 for each one-parameter 68.3% interval. A fixed column has no fitted interval. Read the error status in `fit_results.json` and the fit log. Investigate pegged or unconverged intervals before quoting them.

The fitting script replaces the continuum amplitude with `cflux` for each flux measurement. For the absorbed flux it uses `cflux*tbabs*bbody` (or `cflux*tbabs*powerlaw`) and profiles parameter 3 with `error 1.0 3`. For the unabsorbed flux it uses `tbabs*cflux*bbody` (or `tbabs*cflux*powerlaw`) and profiles parameter 4 with `error 1.0 4`. Both use frozen 1 and 10 keV bounds, a frozen positive continuum normalization, and the baseline fit's free or frozen column. Convert each log-flux limit separately using 10 raised to that limit. Check that the two cflux fits recover the baseline χ² within numerical tolerance before using their results.

`L_X = (8.24 × 10⁴⁵ cm²) F_unabs` at 8.3 kpc. Apply the same factor to the flux limits with distance fixed. The bolometric `bbody` luminosity is `10³⁹ K (8.3/10)² erg s⁻¹`, where K is the XSPEC normalization. The apparent radius is `sqrt(L_bol / (4π σ_SB T⁴))` with `T = 1.16045 × 10⁷ K × (kT / keV)`. Report the radius from the blackbody case you discuss and compare `L_X / (4.9 × 10³³ erg s⁻¹)` for the same case. For a radius interval, profile the equivalent `bbodyrad` normalization or propagate joint kT and normalization samples. Their fitted errors are generally correlated.

`*_flux_density.csv` tabulates absorbed and intrinsic **F_E** in erg cm⁻² s⁻¹ keV⁻¹ and provides the basis for the 4 keV interpolations in the JSON. The `*_eufspec.csv` points are model-dependent unfolded estimates. The script distinguishes the maximum of F_E in the fitted band from the intrinsic blackbody maximum at 2.82144 kT. If the peak falls at a band edge, report a band-edge maximum. A pure power law has no isolated intrinsic spectral peak. The intrinsic blackbody peak in **E F_E** is 3.92069 kT if reporting energy output per logarithmic interval.

## Task D · Finish the report

Use the four `tbabs_*_fit.png` panels to check the residuals, particularly at higher energies. Compare χ² and the residual structure at equal degrees of freedom for the free blackbody and free power law, then repeat for the frozen pair. The two continua are non-nested. A smaller χ² ranks the fits on these data but does not by itself establish a significant preference. If both are poor or the difference is small, describe that outcome.

Use `rea_comparison.csv` to fill the reference table. The percentage difference is signed and uses the reference value as the denominator. Rea et al. report separate joint columns for blackbody and power law (10.3 and 17.3 in units of 10²² cm⁻²). Their 14946 temperature and index are epoch-specific. The reference intrinsic flux is derived from their quoted luminosity. Flux and luminosity comparisons for the power law use the published blackbody benchmark and are marked accordingly. Frozen columns are excluded from the percentage differences.

Fill the report placeholders from `fit_results.json`, `rea_comparison.csv` and `pileup_summary.txt`. Use the selected `*_fit.png` for Figure 3, and its `*_flux_density.png` for Figure 4. In the interpretation, state the observed preference or lack of preference, what the fitted column implies for soft X-ray absorption, the luminosity and spin-down ratio, and which fit values agree with Rea et al. within their uncertainties. A small apparent thermal radius supports localized emission if the blackbody is acceptable. Discuss one-epoch versus joint fitting, the absorption prescription, PSF and background systematics, residual pileup and covariance between absorption and temperature or index.

The PDF source is `build_report.py`. After filling the text and inserting the fit figures, rebuild on Windows with Calibri and ReportLab installed using `python build_report.py SPA660_Report.pdf`. Replace the archive figures and count summary only if you intend to report the reprocessed measurements.

The scripts and report placeholders are prepared. The CIAO and XSPEC jobs need to run in the configured analysis environment before using the fit results.

## Documentation

- [Rea et al. (2013), Table 1 and §2.3](https://arxiv.org/pdf/1307.6331)
- [CIAO reprocessing](https://cxc.harvard.edu/ciao/ahelp/chandra_repro.html), [pileup map](https://cxc.harvard.edu/ciao/ahelp/pileup_map.html) and [point-source spectroscopy](https://cxc.harvard.edu/ciao/threads/pointlike/)
- [XSPEC TBabs](https://heasarc.gsfc.nasa.gov/docs/software/xspec/manual/XSmodelTbabs.html), [phabs](https://heasarc.gsfc.nasa.gov/docs/software/xspec/manual/XSmodelPhabs.html), [cflux](https://heasarc.gsfc.nasa.gov/docs/software/xspec/manual/XSmodelCflux.html) and [confidence intervals](https://heasarc.gsfc.nasa.gov/docs/software/xspec/manual/node107.html)
