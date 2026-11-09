# Live demonstration and replay

The read-only archive checks ran on the saved data in Python. The CIAO, DS9, HEASoft and XSPEC commands below are student replay instructions. A fresh end-to-end run, the optional W-statistic fits, and installation on the demonstration laptop are **NOT VERIFIED**. Capture logs and keep the archived and fresh results separate.

## 1. Set up the laptop

Use Linux or a tested WSL2 Ubuntu environment. CIAO has no native Windows build. On Windows, an administrator can run `wsl --install -d Ubuntu` and restart, then configure the Linux user. Use WSLg for DS9 and XSPEC graphics where available. A native Windows DS9 can instead open the exported FITS image and region files. Test the GUI before the demonstration. [Microsoft WSL GUI guide](https://learn.microsoft.com/en-us/windows/wsl/tutorials/gui-apps)

Install Linux [Miniforge](https://github.com/conda-forge/miniforge#install) first if needed. Use separate environments to avoid mixing CIAO and HEASoft libraries. These candidate pinned Conda commands require package availability for the laptop architecture. If solving fails, use the official older release installers or a preconfigured Linux lab environment and record the actual versions. The reference CALDB version comes from the saved reprocessing history.

```bash
conda create -n spa-ciao -c https://cxc.cfa.harvard.edu/conda/ciao -c conda-forge \
  'ciao=4.17' pyciao ciao-contrib 'caldb_main=4.12.4' ds9 astropy matplotlib
conda create -n spa-xspec -c https://heasarc.gsfc.nasa.gov/FTP/software/conda/ \
  -c conda-forge 'heasoft=6.36' astropy numpy matplotlib
conda activate spa-ciao
ciaover -v
command -v ds9 dmcoords chandra_repro specextract download_chandra_obsid
python -c 'import astropy, numpy, matplotlib, psf'
conda activate spa-xspec
command -v xspec ftgrouppha
python -c 'from xspec import Xset; import astropy, numpy, matplotlib; print(Xset.version)'
```

Target CIAO 4.17.0, CALDB 4.12.4 and XSPEC 12.15.1 / PyXspec 2.1.5. A newer stack is a new replay configuration. Record its versions and compare the results rather than assuming equality. HEASoft 6.37+ can also require a separate XSPEC model-data download. [CIAO 4.17 install](https://cxc.cfa.harvard.edu/ciao4.17/download/conda.html) · [HEASoft install](https://heasarc.gsfc.nasa.gov/docs/software/conda.html) · [Version matrix](https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/versions.html) · [XSPEC model data](https://heasarc.gsfc.nasa.gov/docs/software/xspec/faq.html)

## 2. Retrieve, reprocess and inspect

Run from the repository root after applying the audit patch. Use a fresh working copy for reduction. In the commands below, `products/` is the CIAO working directory and `products_replay/` is the fresh physical-coordinate replay. Keep `products_verified/` for the original saved archive.

```bash
conda activate spa-ciao
download_chandra_obsid 14946
bash reduce.sh --preprocess-only
evt=$(cat products/event_path.txt)
dmlist "$evt[EVENTS]" header,clean
dmlist "$evt[GTI]" data
dmkeypar "$evt[EVENTS]" OBS_ID echo=yes
dmkeypar "$evt[EVENTS]" TIMEDEL echo=yes
dmcopy "$evt[ccd_id=7,energy=1000:8000][bin x=4057.5:4097.5:0.5,y=4089.0:4129.0:0.5]" \
  products/field_1_8keV.fits clobber=yes
ds9 products/field_1_8keV.fits -scale log -zoom to fit
```

The reprocessing script saves `products/reprocess.log`, `products/event_header.txt`, `products/ciao_version.txt`, and the path to a new ObsID 14946 event file. Show CCD 7, the GTI and the live time. The half-pixel sky bins are 0.246 arcsec. The separate saved report PNG uses 0.25 arcsec bins. Use the source position in ICRS and identify the adjacent Sgr A* candidate in DS9. Load the physical region overlays after step 3B to match the reported selection. Confirm the cutout WCS and the six background exclusions in the GUI. An image histogram for aperture counts uses reconstructed ENERGY, whereas the fit uses PI channels.

The preprocess script attempts `pileup_map` on a native-pixel full-energy CCD 7 image and records `products/pileup_summary.txt`. A failure is logged and can be diagnosed separately. The output is a diagnostic in counts/frame. Use the core and annulus maximum and keep the map and log with the event file, frame time and version. Avoid interpreting this as a directly fitted pileup fraction. [CIAO pileup guide](https://cxc.cfa.harvard.edu/ciao/ahelp/pileup_map.html)

## 3. Extract, inspect and record

Choose the matching observation aspect solution from the files below. Prefer the reprocessed solution where available. If there are multiple segments, supply a CIAO stack such as `@products/asol.lis` containing all matching segments of one processing version. Check each header and coverage. Do not combine duplicates of the same time interval.

```bash
find 14946 products/repro -name '*asol*.fits*' -print
read -r -p 'Matching aspect file or @stack: ' asol
test -s "${asol#@}"
punlearn asphist
asphist infile="$asol" outfile=products/asphist_ccd7.fits \
  evtfile="$evt[ccd_id=7]" clobber=no
dmlist 'products/asphist_ccd7.fits[ASPHIST][cols DURATION]' data
```

### 3A. Standard CIAO route for diagnosis

Inspect the two region expressions before running `specextract`. The files in the repository preserve the original FK5 labels. The published workflow actually used the numeric centres in the FITS ICRS sky WCS. The distinction can move events across the narrow aperture edges. A standard region extraction with `correctpsf=yes` applies CIAO's `arfcorr` and writes to `products/`.

```bash
dmlist "$evt[ccd_id=7][sky=region(source.reg)]" counts
dmlist "$evt[ccd_id=7][sky=region(background.reg)]" counts
bash reduce.sh products/asphist_ccd7.fits
```

If either count is zero, save the exact command, `dmlist` output, tool log, region file, and relevant headers. Check that `evt` is the full reprocessed event file, that it contains CCD 7, that region text is readable, and that the DS9 positions lie inside the image. Compare with the physical filters from 3B. A nonempty file on disk can still be an empty spectrum. The original failed run's command and intermediate PHA are missing, so the reason for its zero counts remains unresolved.

### 3B. Physical-coordinate replay of the reported method

This is the candidate replay of the alternative route. It uses the same numerical centre, aperture radii, masks, and method as the archived event selections. It regenerates its own point-source responses for the same observation before applying a single REEF correction. The original response headers show an unweighted `specextract` call with `correctpsf=no` and an explicit CCD 7 aspect histogram. Exact commands from that run were not archived. Compare new results to the saved archive instead of assuming identity.

```bash
python reextract_spa660.py select --select-only --output products_replay
dmlist "products_replay/source_events.fits[EVENTS]" counts
dmlist "products_replay/background_events.fits[EVENTS]" counts
ds9 products/field_1_8keV.fits -scale log -zoom to fit \
  -regions load products_replay/source_physical.reg \
  -regions load products_replay/background_physical.reg

punlearn dmextract
dmextract "products_replay/source_events.fits[bin pi]" \
  products_replay/magnetar.pi opt=pha1 clobber=no
punlearn dmextract
dmextract "products_replay/background_events.fits[bin pi]" \
  products_replay/magnetar_bkg.pi opt=pha1 clobber=no

dmlist "$evt[ccd_id=7][sky=region(products_replay/source_physical.reg)]" counts
dmlist "$evt[ccd_id=7][sky=region(products_replay/background_physical.reg)]" counts
find products/repro -name '*bpix*.fits*' -o -name '*msk*.fits*'
read -r -p 'Matching reprocessed bad-pixel file: ' bpix
read -r -p 'Matching observation mask file: ' msk
test -s "$bpix" && test -s "$msk"
punlearn specextract
python specextract_compat.py \
  infile="$evt[ccd_id=7][sky=region(products_replay/source_physical.reg)]" \
  asp=products/asphist_ccd7.fits badpixfile="$bpix" mskfile="$msk" \
  outroot=products_replay/point_response refcoord="266.4173708 -29.0082889" \
  weight=no correctpsf=no bkgresp=no grouptype=NONE clobber=no verbose=2 \
  2>&1 | tee products_replay/response.log
cp products_replay/point_response.arf products_replay/magnetar.arf
cp products_replay/point_response.rmf products_replay/magnetar.rmf

punlearn dmcoords
dmcoords "$evt" asol=none option=cel celfmt=deg ra=266.4173708 dec=-29.0082889
theta=$(pget dmcoords theta)
phi=$(pget dmcoords phi)
python make_reef_arf.py products_replay/magnetar.arf products_replay/magnetar.corr.arf \
  --theta-arcmin "$theta" --phi-deg "$phi" | tee products_replay/reef.log
python reextract_spa660.py verify --output products_replay
punlearn dmgroup
dmgroup infile=products_replay/magnetar.pi outfile=products_replay/magnetar_grp.pi \
  grouptype=NUM_CTS grouptypeval=50 xcolumn=CHANNEL ycolumn=COUNTS clobber=no
python inspect_events.py "$evt" --output products_replay/event_check

for name in magnetar.pi magnetar_grp.pi magnetar_bkg.pi; do
  for key in EXPOSURE BACKSCAL BACKFILE RESPFILE ANCRFILE; do
    printf '%s %s ' "$name" "$key"
    dmkeypar "products_replay/$name[SPECTRUM]" "$key" echo=yes
  done
done
dmlist products_replay/magnetar.rmf blocks
dmlist products_replay/magnetar.corr.arf cols
```

Stop if any command fails or the CIAO filters disagree with the selected event lists. Confirm the aspect histogram duration against the CCD 7 GTI-derived live time, source and background counts channel by channel, and the response source location. BACKSCAL is computed from the masked sky areas in this replay. Inspect the footprint, chip edge and bad pixels before relying on a geometric area. Different sky regions require new spectra and response validation, with a new aperture correction where appropriate. The uncorrected ARF should lack `PSF_FRAC`. The corrected ARF should contain exactly one correction. `pget dmcoords theta` returns arcminutes even when the printed screen label uses arcseconds. The helper records its CALDB file basename and angles. [asphist](https://cxc.cfa.harvard.edu/ciao/ahelp/asphist.html) · [specextract](https://cxc.cfa.harvard.edu/ciao4.17/ahelp/specextract.html) · [REEF API](https://cxc.cfa.harvard.edu/ciao/ahelp/psf.html) · [dmcoords](https://cxc.cfa.harvard.edu/ciao/ahelp/dmcoords.html)

The older `products_fast` path and `repair_spa660_exposure.py` are historical. The repair utility only checks band totals approximately and cannot rescue a zero-count spectrum. The physical replay above uses `--select-only` and no dependency on `products_fast`.

## 4. Verify and refit the saved products

Run the validator first. Then copy only the response and PHA inputs into a fresh run directory. Keep the original archive for side-by-side comparison. For the fresh CIAO extraction, pass `--products products_replay` instead.

```bash
conda activate spa-xspec
unzip SPA660_Results_verified.zip
python verify_products.py products_verified
mkdir fits_replay
cp products_verified/{magnetar.pi,magnetar_grp.pi,magnetar_bkg.pi,magnetar.rmf,magnetar.arf,magnetar.corr.arf} fits_replay/
python fit_xspec.py --products fits_replay --rea-crosscheck
```

Expect the replay directory to contain six `.xcm` sessions, `xspec.log`, `fit_results.json`, `psf_fraction.csv`, spectra/residual PNGs, unfolded/flux-density CSVs, and comparison files. Compare actual parameters, error statuses, χ²/dof, and retained energy coverage to the archived JSON and log. The saved-data check is an executed audit. All new fits in this step are **NOT VERIFIED** until the students run them.

For the live XSPEC display, change into `fits_replay`, launch `xspec`, and enter:

```text
@tbabs_bbody_free.xcm
show all
fit
error 1.0 1-3
setplot energy
cpd /xs
plot data delchi
flux 1 10
@tbabs_powerlaw_free.xcm
show all
fit
error 1.0 1-3
plot data delchi
```

Also restore `tbabs_bbody_fixed.xcm`, `tbabs_powerlaw_fixed.xcm`, `phabs_bbody_free.xcm`, and `phabs_powerlaw_free.xcm` for the complete comparison. The two fixed sessions intentionally share the blackbody reference column (10.3 in XSPEC units). For the fixed cases, use `error maximum 10.0 1.0 2-3` and interpret intervals conditional on fit quality. Use `cflux*absorption*continuum` for absorbed flux and `absorption*cflux*continuum` for intrinsic flux. The script fixes the continuum norm and integration bounds and profiles the log flux. Convert log flux to linear units and multiply the intrinsic flux and its interval limits by `4π(8.3 kpc in cm)²`. Record which interval is statistical and which uncertainties (distance, response, model) were held fixed. [XSPEC cflux](https://heasarc.gsfc.nasa.gov/docs/software/xspec/manual/XSmodelCflux.html)

## 5. Optional grouping and PSF sensitivity tests

The archived grouped χ² fit ignores the large group crossing 1 keV and the poor-quality final group. In the archived RMF, the 1–10 keV PI-centre channels are 69–685. Validate the channel endpoints for any new RMF. There are 325 zero-background channels out of 617 in the archived nominal band. To test use of the edges, first use the ungrouped PHA in a separate analysis and inspect the low-count behaviour. Grouping by a fixed minimum background count can help reduce the many empty-background bins in W-statistic. Predefine a criterion and report its sensitivity for both continua rather than choosing a grouping after inspecting fit preferences.

```bash
cd fits_replay
ftgrouppha infile=magnetar.pi backfile=magnetar_bkg.pi \
  outfile=magnetar_wstat.pi grouptype=bmin groupscale=1 \
  minchannel=69 maxchannel=685 clobber=no
xspec
```

```text
data 1:1 magnetar_wstat.pi
backgrnd 1 magnetar_bkg.pi
response 1 magnetar.rmf
arf 1 magnetar.corr.arf
statistic cstat
abund angr
xsect bcmc
xset TBABSVERSION 2
setplot energy
ignore bad
ignore **-1.0 10.0-**
model tbabs*bbody & 10.3 & 0.9 & 2.0e-4
fit
show all
plot data delchi
model tbabs*powerlaw & 17.3 & 4.2 & 0.1
fit
show all
plot data delchi
```

With a supplied Poisson background, XSPEC's `cstat` uses W-statistic. Check the retained edge groups and sparse background, including the final group, and compare actual confidence intervals and residuals. W-statistic per degree of freedom is not a χ² probability. This optional experiment is **NOT VERIFIED** and does not replace the archived fit. [Low-count spectra](https://heasarc.gsfc.nasa.gov/docs/software/xspec/manual/node395.html) · [ftgrouppha](https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/help/ftgrouppha.html)

Run the archived-response aperture sensitivity from the repo root with `python polish_xspec.py sensitivity --products products_verified` in the HEASoft environment. Preserve its outputs in `polish_checks/`. The optional `band` action changes the integrated flux interval while keeping the fitted groups. In CIAO, `bash polish_ciao.sh arfcorr` produces a shared-CALDB method cross-check, not independent PSF calibration. `project_chart.sh` is optional and takes matching ChaRT ray files plus an explicit aspect file. MARX, the rays, and the model for the ray simulation are additional inputs. These supplementary commands remain **NOT VERIFIED** here. A uniform effective-area multiplier is a normalization sensitivity test, not a measurement of PSF uncertainty.

## 6. Capture and submit the evidence

Save the full reprocessed event file, aspect solution identifiers and histogram, a CIAO version record, extraction commands and logs, masked region files, pileup map and log, ARF/RMF/PHA files, and XSPEC log and sessions from the actual demonstration. Record any differences before the students edit the report. Show the full 1–10 keV event counts separately from the retained fitted groups.

Share the original audited `SPA660_Results_verified.zip` with the examiner in the course submission system or a Drive folder restricted to the examiner. The original FITS history includes system-specific paths. Give the examiner access and test the link with the intended permissions. Add that authorized link in the submission notes or repository documentation according to the course policy. Keep the official Chandra download link for raw data. Submit one group PDF and confirm that its page 1 GitHub link opens the correct public repo.
