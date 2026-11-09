#!/usr/bin/env bash
# Run from SPA660 after conda activate ciao.
set -euo pipefail
action=${1:-pileup}
evt=$(cat products/event_path.txt)
test -f "$evt"
mkdir -p polish_checks

punlearn dmcoords
dmcoords "$evt" option=cel celfmt=deg ra=266.4173708 dec=-29.0082889
xsky=$(pget dmcoords x)
ysky=$(pget dmcoords y)
theta=$(pget dmcoords theta)
phi=$(pget dmcoords phi)
printf 'Source X=%s Y=%s theta=%s arcmin phi=%s deg\n' "$xsky" "$ysky" "$theta" "$phi"
export SPA_X="$xsky" SPA_Y="$ysky"
read -r xmin xmax ymin ymax rin rout < <(python - <<'PY'
import math, os
x, y = float(os.environ['SPA_X']), float(os.environ['SPA_Y'])
print(math.floor(x)-64, math.floor(x)+64, math.floor(y)-64, math.floor(y)+64,
      0.75/0.492, 2.0/0.492)
PY
)
ann="annulus($xsky,$ysky,$rin,$rout)"
core="circle($xsky,$ysky,$rin)"
printf 'physical\n%s\n' "$ann" > polish_checks/annulus_physical.reg
printf 'physical\n%s\n' "$core" > polish_checks/core_physical.reg

if [[ "$action" == pileup ]]; then
    # Full event energy range and native ACIS sky-pixel sampling are essential.
    dmcopy "$evt[ccd_id=7][bin x=$xmin:$xmax:1,y=$ymin:$ymax:1]" \
        polish_checks/ccd7_native.fits clobber=yes
    dmkeypar polish_checks/ccd7_native.fits LIVETIME echo=yes
    dmkeypar polish_checks/ccd7_native.fits TIMEDEL echo=yes
    punlearn pileup_map
    pileup_map infile=polish_checks/ccd7_native.fits \
        outfile=polish_checks/pileup_map.fits clobber=yes verbose=1
    dmstat "polish_checks/pileup_map.fits[sky=$core]" centroid=no
    coremax=$(pget dmstat out_max)
    dmstat "polish_checks/pileup_map.fits[sky=$ann]" centroid=no
    annmax=$(pget dmstat out_max)
    printf 'core_max_counts_per_frame=%s\nannulus_max_counts_per_frame=%s\n' \
        "$coremax" "$annmax" | tee polish_checks/pileup_summary.txt
elif [[ "$action" == arfcorr ]]; then
    test -f products_verified/magnetar.arf
    # Image cutout for speed, subpixels for the small annulus.
    dmcopy "$evt[ccd_id=7][bin x=$xmin:$xmax:0.25,y=$ymin:$ymax:0.25]" \
        polish_checks/psf_stamp.fits clobber=yes
    punlearn arfcorr
    arfcorr infile=polish_checks/psf_stamp.fits arf=none \
        outfile=polish_checks/arfcorr_4keV.fits region="$ann" \
        x="$xsky" y="$ysky" energy=4 nsubpix=5 nfracpix=8 \
        radlim=2 clobber=yes verbose=1 2>&1 | tee polish_checks/arfcorr_4keV.log
    arfcorr infile=polish_checks/psf_stamp.fits \
        arf=products_verified/magnetar.arf \
        outfile=polish_checks/magnetar.arfcorr.arf region="$ann" \
        x="$xsky" y="$ysky" energy=0 e_step=50 nsubpix=5 nfracpix=8 \
        radlim=2 clobber=yes verbose=1 2>&1 | tee polish_checks/arfcorr_energy.log
    python - <<'PY'
import numpy as np
from astropy.io import fits
a = fits.getdata('polish_checks/magnetar.arfcorr.arf', 'SPECRESP')
e = (a['ENERG_LO'] + a['ENERG_HI']) / 2
f = a['PSF_FRAC']
good = (e >= 1) & (e <= 10)
if not np.all(np.isfinite(f[good]) & (f[good] > 0) & (f[good] <= 1)):
    raise ValueError('Invalid arfcorr fractions. Review coordinate/image interpretation.')
print('ARF-interpolated f_PSF(4 keV) =', np.interp(4, e, f))
print('Saved REEF reference at 4 keV = 0.1077')
np.savetxt('polish_checks/arfcorr_fraction.csv', np.column_stack((e, f)),
           delimiter=',', header='energy_keV,annular_fraction', comments='')
PY
else
    printf 'Use bash polish_ciao.sh pileup or bash polish_ciao.sh arfcorr\n' >&2
    exit 2
fi
