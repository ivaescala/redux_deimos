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
#   marz     (optional) Marz input file for flagging extragalactic objects
#   dmost    velocities and equivalent widths with dmost, from the collate outputs
#
# Every path is absolute, so the script can be run from any folder.
# Requires: DEIMOS_RAW (koa_deimos_fetch.py --outdir) and DEIMOS_RDX (or DEIMOS_REDUX).
# marz/dmost also use DMOST_DIR (patched dmost clone, default <kit>/../dmost),
# DMOST_DATA (templates/ and Other_data/, default <kit>/../dmost_data) and
# DMOST_FRESH=1 (archive earlier dmost outputs and start over).
# The kit files (this script, deimos.par, assign_calib_by_night.py, skysub_qa.py,
# wavecal_qa.py, wavecal_exclude.py, merge_spec1d.py, dmost_prep.py)
# must sit in the same folder. Re-running "setup" overwrites the .pypeit file(s).
set -euo pipefail

if [ $# -ne 2 ]; then
  sed -n '2,21p' "${BASH_SOURCE[0]}"; exit 1
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
         "$KIT/wavecal_exclude.py" "$KIT/merge_spec1d.py" "$KIT/dmost_prep.py"; do
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

# ---- dmost --------------------------------------------------------------------
# dmost expects one flat folder per mask under $DEIMOS_REDUX, and templates and raw
# frames under $DEIMOS_RAW. Both are set only for the dmost processes, to
#   DEIMOS_REDUX = $DEIMOS_RDX/<semester>/dmost/          (DMOST_REDUX)
#   DEIMOS_RAW   = $DMOST_REDUX/<name>/dmost_raw/         (built by dmost_prep.py)
# so the rest of this script and the PypeIt stages are unaffected.
dmost_env() {
  DMOST_DIR=${DMOST_DIR:-$(dirname "$KIT")/dmost}
  DMOST_DATA=${DMOST_DATA:-$(dirname "$KIT")/dmost_data}
  DMOST_REDUX=${RDX%/}/$(dirname "$MASK")/dmost
  [ -f "$DMOST_DIR/scripts/run_dmost.py" ] \
    || { echo "reduce_mask.sh: no dmost clone at $DMOST_DIR (set DMOST_DIR)" >&2; exit 1; }
  [ -f "$DMOST_DIR/PATCHED_BY_REDUX" ] \
    || { echo "reduce_mask.sh: dmost is unpatched; run: python $KIT/patch_dmost.py $DMOST_DIR" >&2; exit 1; }
  NCFG=$(configs | wc -l)
}
# dmost's name for a configuration: the mask name, plus _<X> if the mask has several
dmost_name() {
  local base; base=$(basename "$MASK")
  if [ "$NCFG" -eq 1 ]; then echo "$base"; else echo "${base}_$(basename "$1" | cut -d_ -f3)"; fi
}
# run a dmost script inside its working folder, with dmost's view of the paths
dmost_run() {  # dmost_run <working folder> <dmost DEIMOS_RAW> <script> [args]
  local w=$1 raw=$2; shift 2
  ( cd "$w" && DEIMOS_REDUX="$DMOST_REDUX/" DEIMOS_RAW="$raw" MPLBACKEND=Agg \
      PYTHONPATH="$DMOST_DIR${PYTHONPATH:+:$PYTHONPATH}" python "$@" )
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
  marz)
    dmost_env
    for CFG in $(configs); do
      NAME=$(dmost_name "$CFG")
      W=$DMOST_REDUX/$NAME
      python "$KIT/dmost_prep.py" --for marz --cfg "$CFG" --name "$NAME" \
          --workroot "$DMOST_REDUX" --rawdir "$RAWDIR"
      dmost_run "$W" "${DEIMOS_RAW%/}/" "$DMOST_DIR/scripts/run_collate1d_marz.py" --mask "$NAME"
      echo "Marz input: $DMOST_REDUX/marz_files/marz_${NAME}.fits"
      echo "Save the Marz results as $DMOST_REDUX/marz_files/marz_${NAME}_MG.mz before the dmost stage."
    done
    ;;
  dmost)
    dmost_env
    FRESH=""
    if [ "${DMOST_FRESH:-0}" = 1 ]; then FRESH=--fresh; fi
    for CFG in $(configs); do
      NAME=$(dmost_name "$CFG")
      W=$DMOST_REDUX/$NAME
      python "$KIT/dmost_prep.py" --for dmost --cfg "$CFG" --name "$NAME" --workroot "$DMOST_REDUX" \
          --rawdir "$RAWDIR" --data "$DMOST_DATA" $FRESH
      [ -f "$DMOST_REDUX/marz_files/marz_${NAME}_MG.mz" ] \
        || echo "No Marz results for $NAME: every object is treated as a star."
      dmost_run "$W" "$W/dmost_raw/" "$DMOST_DIR/scripts/run_dmost.py" "$NAME"
      dmost_run "$W" "$W/dmost_raw/" - "$W/dmost/dmost_$NAME.fits" <<'EOF'
import sys
import numpy as np
from dmost.core.dmost_utils import read_dmost
slits, mask = read_dmost(sys.argv[1])
star = slits['marz_flag'] < 2
good = star & (slits['dmost_v_err'] > 0)
coadd = good & (slits['coadd_flag'] == 1)
print(f"{mask['maskname'][0]}: {len(mask)} exposures, {len(slits)} objects, {star.sum()} stars, "
      f"{good.sum()} with velocities ({coadd.sum()} from coadded spectra)")
if good.any():
    print(f"median velocity error {np.median(slits['dmost_v_err'][good]):.1f} km/s")
EOF
      echo "dmost table: $W/dmost/dmost_$NAME.fits | log: $W/${NAME}_dmost.log | plots: $W/QA/"
    done
    ;;
  *)
    echo "reduce_mask.sh: unknown stage '$STAGE' (setup|calib|science|qa|collate|marz|dmost)" >&2; exit 1 ;;
esac
