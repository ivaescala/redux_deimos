#!/usr/bin/env bash
# reduce_mask.sh <SEMESTER>/<MASK> <stage>      e.g. reduce_mask.sh 2022B/M32RA1 setup
#
# PypeIt reduction of one Keck/DEIMOS slitmask, following Geha et al. (2026).
# Stages, run in order and inspect between them:
#   setup    configuration census, .pypeit file(s) with deimos.par, calib groups per night
#   calib    calibrations only + rule-based wavelength-solution check (wavecal_qa.csv)
#   science  full reduction (reuses calibrations)
#   qa       slit report + sky-residual statistics vs. distance from M32
#   collate  exclude slits failing wavecal_qa, merge with any recovered slits, coadd 1D spectra
#
# Every path is absolute, so the script can be run from any folder.
# Requires: DEIMOS_RAW (koa_deimos_fetch.py --outdir) and DEIMOS_RDX (or DEIMOS_REDUX).
# The kit files (this script, deimos.par, assign_calib_by_night.py, skysub_qa.py,
# wavecal_qa.py, wavecal_exclude.py, merge_spec1d.py)
# must sit in the same folder. Re-running "setup" overwrites the .pypeit file(s).
set -euo pipefail

if [ $# -ne 2 ]; then
  sed -n '2,15p' "${BASH_SOURCE[0]}"; exit 1
fi
MASK=$1
STAGE=$2

# ---- locations ----------------------------------------------------------------
: "${DEIMOS_RAW:?DEIMOS_RAW is not set (source ~/.bashrc)}"
RDX=${DEIMOS_RDX:-${DEIMOS_REDUX:-}}
: "${RDX:?DEIMOS_RDX / DEIMOS_REDUX is not set (source ~/.bashrc)}"

KIT=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
PAR=${DEIMOS_PAR:-$KIT/deimos.par}
for f in "$PAR" "$KIT/assign_calib_by_night.py" "$KIT/skysub_qa.py" "$KIT/wavecal_qa.py" \
         "$KIT/wavecal_exclude.py" "$KIT/merge_spec1d.py"; do
  [ -f "$f" ] || { echo "reduce_mask.sh: missing $f" >&2; exit 1; }
done

RAWDIR=${DEIMOS_RAW%/}/$MASK
OUT=${RDX%/}/$MASK
[ -d "$RAWDIR" ] || { echo "reduce_mask.sh: no raw folder $RAWDIR" >&2; exit 1; }
mkdir -p "$OUT"
echo "mask $MASK | raw $RAWDIR | redux $OUT | stage $STAGE"

# one folder per instrument configuration: $OUT/keck_deimos_A, _B, ...
configs() {
  local d found=0
  for d in "$OUT"/keck_deimos_?; do
    [ -f "$d/$(basename "$d").pypeit" ] && { echo "$d"; found=1; }
  done
  [ $found -eq 1 ] || { echo "reduce_mask.sh: no .pypeit files in $OUT; run setup first" >&2; exit 1; }
}

# ---- stages -------------------------------------------------------------------
case "$STAGE" in
  setup)
    pypeit_setup -s keck_deimos -r "$RAWDIR" -d "$OUT"
    cat "$OUT"/setup_files/keck_deimos*.sorted
    pypeit_setup -s keck_deimos -r "$RAWDIR" -d "$OUT" -c all -p "$PAR"
    for CFG in $(configs); do
      PF=$CFG/$(basename "$CFG").pypeit
      grep -Eq "spec_method *= *skip" "$PF" \
        || { echo "reduce_mask.sh: parameter block missing from $PF" >&2; exit 1; }
      python "$KIT/assign_calib_by_night.py" "$PF"
      echo "--- data block of $PF"
      sed -n '/data read/,/data end/p' "$PF"
    done
    ;;
  calib)
    for CFG in $(configs); do
      run_pypeit "$CFG/$(basename "$CFG").pypeit" -r "$CFG" -c
      pypeit_chk_wavecalib "$CFG"/Calibrations/WaveCalib_*.fits > "$CFG/wavecalib_report.txt"
      python "$KIT/wavecal_qa.py" "$CFG"          # -> $CFG/wavecal_qa.csv, prints failing slits
      echo "Inspect slit edges interactively, e.g.:"
      echo "  pypeit_chk_edges $CFG/Calibrations/Edges_$(basename "$CFG" | cut -d_ -f3)_0_MSC01.fits.gz"
    done
    ;;
  science)
    for CFG in $(configs); do
      run_pypeit "$CFG/$(basename "$CFG").pypeit" -r "$CFG" -o
    done
    ;;
  qa)
    for CFG in $(configs); do
      for f in "$CFG"/Science/spec2d_*.fits; do
        echo "### $(basename "$f")"
        pypeit_parse_slits "$f"
      done > "$CFG/slit_report.txt"
      python "$KIT/skysub_qa.py" "$CFG"/Science/spec2d_*.fits --out "$CFG/skysub_qa.csv"
    done
    ;;
  collate)
    for CFG in $(configs); do
      # PypeIt extracts slits whose wavelength fit is wrong but exists; drop them by maskdef_id
      python "$KIT/wavecal_exclude.py" "$CFG"     # -> none rows in $CFG/selection.csv
      # always merge, so the exclusions apply whether or not slits were recovered (section 7)
      python "$KIT/merge_spec1d.py" --main "$CFG/Science" --test "$OUT/test_maskfrac/Science" \
          --select "$CFG/selection.csv" --out "$CFG/Science_merged"
      pypeit_collate_1d --spec1d_files "$CFG"/Science_merged/spec1d_*.fits \
          --tolerance 0.5 --wv_rms_thresh 0.4 --outdir "$CFG/collate1d"
    done
    ;;
  *)
    echo "reduce_mask.sh: unknown stage '$STAGE' (setup|calib|science|qa|collate)" >&2; exit 1 ;;
esac
