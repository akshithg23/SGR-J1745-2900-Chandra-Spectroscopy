#!/usr/bin/env bash
# Use an active CIAO/CALDB environment. An explicit aspect file or @stack is
# required for extraction when there is no unique reprocessed aspect solution.
set -euo pipefail
cd -- "$(dirname -- "$0")"

aspect=${1:-}
if (( $# > 1 )); then
    echo "Usage: bash reduce.sh [--preprocess-only | ASPECT_FILE_OR_@STACK]" >&2
    exit 2
fi
for tool in chandra_repro ciaover dmlist dmcopy punlearn; do
    command -v "$tool" >/dev/null || { echo "Activate CIAO first: missing $tool" >&2; exit 1; }
done

mkdir -p products
ciaover -v > products/ciao_version.txt
evt=""
if [[ -s products/event_path.txt ]]; then
    IFS= read -r evt < products/event_path.txt
fi
if [[ -z "$evt" || ! -s "$evt" ]]; then
    # chandra_repro needs a writable copy of the complete archive dataset.
    archive=/FTP/chandra/data/byobsid/6/14946
    if [[ ! -d 14946/primary || ! -d 14946/secondary ]]; then
        if [[ -d "$archive/primary" && -d "$archive/secondary" ]]; then
            mkdir -p 14946
            cp -RL -- "$archive/primary" "$archive/secondary" 14946/
        else
            download_chandra_obsid 14946
        fi
    fi

    punlearn chandra_repro
    chandra_repro indir=14946 outdir=products/repro verbose=2 clobber=yes \
        2>&1 | tee products/reprocess.log
    shopt -s nullglob
    events=(products/repro/*repro_evt2.fits)
    [[ ${#events[@]} -eq 1 ]] || { echo "Expected one reprocessed event file" >&2; exit 1; }
    evt=${events[0]}
    printf '%s\n' "$evt" > products/event_path.txt
else
    echo "Using the reprocessed event file recorded in products/event_path.txt."
fi
dmlist "$evt" header,clean > products/event_header.txt

punlearn dmcopy
dmcopy "$evt[ccd_id=7,energy=1000:8000][bin x=4057.5:4097.5:0.5,y=4089.0:4129.0:0.5]" \
    products/field_1_8keV.fits clobber=yes

# The pileup map is a diagnostic, independent of PI and response generation.
# Use native ACIS pixels and the full energy range.
pileup_diagnostic() {
    dmcopy "$evt[ccd_id=7][bin x=3072.5:5120.5:1,y=3072.5:5120.5:1]" \
        products/ccd7_native.fits clobber=yes || return 1
    punlearn pileup_map || return 1
    pileup_map products/ccd7_native.fits products/pileup_map.fits clobber=yes || return 1
    punlearn dmcoords || return 1
    dmcoords "$evt" asol=none option=cel celfmt=deg ra=266.4173708 dec=-29.0082889 || return 1
    local sx sy
    sx=$(pget dmcoords x) || return 1
    sy=$(pget dmcoords y) || return 1
    punlearn dmstat || return 1
    # 0.75 and 2.0 arcsec in 0.492 arcsec ACIS sky pixels.
    echo "Core counts per frame (maximum and mean)"
    dmstat "products/pileup_map.fits[sky=circle($sx,$sy,1.52439)]" centroid=no || return 1
    echo "Annulus counts per frame (maximum and mean)"
    dmstat "products/pileup_map.fits[sky=annulus($sx,$sy,1.52439,4.06504)]" centroid=no
}
if ! pileup_diagnostic 2>&1 | tee products/pileup_summary.txt; then
    echo "Pileup summary unavailable. Check products/pileup_summary.txt before interpreting pileup." >&2
fi
echo "Open products/field_1_8keV.fits in DS9 or JS9 with source.reg and background.reg."
[[ "$aspect" == "--preprocess-only" ]] && exit 0

if [[ -z "$aspect" ]]; then
    shopt -s nullglob
    candidates=(products/repro/*asol*.fits products/repro/*asol*.fits.gz)
    if (( ${#candidates[@]} != 1 )); then
        echo "Pass the matching aspect solution, @stack, or CCD 7 aspect histogram to reduce.sh." >&2
        exit 2
    fi
    aspect=${candidates[0]}
fi
[[ -s "${aspect#@}" ]] || { echo "Missing aspect input: $aspect" >&2; exit 1; }

# Standard CIAO route. Keep its outputs separate from the archived REEF run.
punlearn specextract
python specextract_compat.py infile="$evt[ccd_id=7][sky=region(source.reg)]" \
    bkgfile="$evt[ccd_id=7][sky=region(background.reg)]" outroot=products/magnetar \
    asp="$aspect" refcoord="266.4173708 -29.0082889" weight=no correctpsf=yes binarfcorr=0.25 \
    bkgresp=no grouptype=NUM_CTS binspec=50 bkg_grouptype=NONE clobber=yes verbose=2 \
    2>&1 | tee products/extraction.log

for product in magnetar_grp.pi magnetar_bkg.pi magnetar.rmf magnetar.arf magnetar.corr.arf; do
    test -s "products/$product" || { echo "Missing spectral product: products/$product" >&2; exit 1; }
done
for pha in magnetar_grp.pi magnetar_bkg.pi; do
    punlearn dmstat
    dmstat "products/$pha[SPECTRUM][cols COUNTS]" centroid=no
    counts=$(pget dmstat out_sum)
    awk -v counts="$counts" 'BEGIN {exit !(counts+0 > 0)}' ||
        { echo "Empty spectrum: products/$pha. See WORKFLOW.md before fitting." >&2; exit 1; }
done
