#!/usr/bin/env bash
# Usage: bash project_chart.sh broad 'polish_checks/rays_broad/*.fits*'
# Supply fresh ChaRT realizations for each spectrum, matching the same aspect file.
set -euo pipefail
label=${1:?Supply broad or mono4}
pattern=${2:?Supply a quoted glob for the downloaded ChaRT FITS rayfiles}
case "$label" in broad|mono4) ;; *) exit 2 ;; esac
evt=$(cat products/event_path.txt)
test -f "$evt"
asol_candidates=(14946/primary/pcadf14946_000N001_asol1.fits*)
asol=${asol_candidates[0]}
test -f "$asol"
command -v marx
command -v marx2fits
find polish_checks -maxdepth 3 -type f -path "$pattern" | sort > "polish_checks/${label}_rays.lis"
test -s "polish_checks/${label}_rays.lis"
punlearn simulate_psf
simulate_psf infile="$evt" outroot="polish_checks/chart_${label}" \
    ra=266.4173708 dec=-29.0082889 spectrum=none \
    simulator=file rayfile="@polish_checks/${label}_rays.lis" projector=marx \
    asolfile="$asol" blur=0.25 pileup=no readout_streak=no \
    ideal=no extended=yes binsize=0.25 minsize=256 keepiter=yes verbose=1
python polish_xspec.py fraction --label "$label" \
    --events polish_checks/chart_"${label}"_i[0-9][0-9][0-9][0-9]_projrays.fits
