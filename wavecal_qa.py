#!/usr/bin/env python
"""
wavecal_qa.py  <configuration folder, e.g. $DEIMOS_RDX/M32RA1/keck_deimos_A>

Rule-based check of every slit's wavelength solution, applied identically to all
masks so the rejection list is reproducible. A slit is "bad" if any check fails:

  - fitted lines do not run blue-to-red (minWave < Wave_cen < maxWave)
  - dispersion differs from the detector median by more than DISP_TOL
  - fewer than NLIN_MIN arc lines in the fit
  - fitted lines cover less than COV_MIN of the slit spectrum, or more than 100%
  - RMS above RMS_MAX pixels (the Geha et al. 2026 cut)
  - PypeIt itself flagged the slit (BADWVCALIB)

The RMS cut alone misses solutions that lock onto the wrong lines with a small
RMS (e.g. a dispersion several times too large), hence the other checks.
Alignment boxes are labelled "box"; slits with no solution at all, "fail".

Writes <folder>/wavecal_qa.csv (one row per slit, with maskdef_id for matching
to targets downstream) and prints every slit that is not "ok".
"""
import argparse
import glob
import os

import numpy as np
from astropy.table import Table, vstack

RMS_MAX = 0.4      # pixels
NLIN_MIN = 20      # arc lines used in the fit
COV_MIN = 60.0     # percent of the slit spectrum spanned by fitted lines
DISP_TOL = 0.05    # fractional deviation of dWave from the detector median

COLS = ["minWave", "Wave_cen", "maxWave", "dWave", "Nlin", "IDs_Wave_cov(%)",
        "measured_fwhm", "RMS"]


def classify(diag, flags):
    """Return (status, reason) lists for a wave_diagnostics table."""
    solved = np.asarray(diag["Nlin"]) > 0
    med_disp = np.median(np.asarray(diag["dWave"])[solved]) if solved.any() else np.nan
    status, reason = [], []
    for row, fl in zip(diag, flags):
        if "BOXSLIT" in fl:
            status.append("box"); reason.append(""); continue
        if row["Nlin"] == 0:
            status.append("fail"); reason.append("no solution"); continue
        why = []
        if not row["minWave"] < row["Wave_cen"] < row["maxWave"]:
            why.append("wavelength order")
        if abs(row["dWave"] / med_disp - 1) > DISP_TOL:
            why.append(f"dWave {row['dWave']:.3f} (median {med_disp:.3f})")
        if row["Nlin"] < NLIN_MIN:
            why.append(f"Nlin {row['Nlin']}")
        cov = row["IDs_Wave_cov(%)"]
        if not COV_MIN <= cov <= 100.5:
            why.append(f"coverage {cov:.0f}%")
        if row["RMS"] > RMS_MAX:
            why.append(f"RMS {row['RMS']:.2f}")
        if "BADWVCALIB" in fl:
            why.append("PypeIt BADWVCALIB")
        status.append("bad" if why else "ok")
        reason.append("; ".join(why))
    return status, reason


def main():
    from pypeit.slittrace import SlitTraceSet
    from pypeit.wavecalib import WaveCalib

    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    p.add_argument("config_dir")
    a = p.parse_args()
    caldir = os.path.join(a.config_dir, "Calibrations")

    tables = []
    for wfile in sorted(glob.glob(os.path.join(caldir, "WaveCalib_*.fits"))):
        key = os.path.basename(wfile)[len("WaveCalib_"):-len(".fits")]      # e.g. A_0_MSC01
        diag = WaveCalib.from_file(wfile, chk_version=False).wave_diagnostics(print_diag=False)

        sfile = os.path.join(caldir, f"Slits_{key}.fits.gz")
        slits = SlitTraceSet.from_file(sfile, chk_version=False) if os.path.exists(sfile) else None
        index = {} if slits is None else {int(s): i for i, s in enumerate(slits.spat_id)}
        flags, mids = [], []
        for sid in diag["SpatOrderID"]:
            i = index.get(int(sid))
            flags.append([] if i is None else slits.bitmask.flagged_bits(slits.mask[i]))
            mids.append(-1 if i is None or slits.maskdef_id is None else int(slits.maskdef_id[i]))

        status, reason = classify(diag, flags)
        t = Table()
        t["calib_key"] = [key] * len(diag)
        t["det"] = [key.split("_")[-1]] * len(diag)
        t["spat_id"] = diag["SpatOrderID"]
        t["maskdef_id"] = mids
        t["status"] = status
        t["reason"] = reason
        for c in COLS:
            t[c] = diag[c]
        tables.append(t)

    if not tables:
        raise SystemExit(f"No WaveCalib files in {caldir}")
    t = vstack(tables)
    out = os.path.join(a.config_dir, "wavecal_qa.csv")
    t.write(out, overwrite=True)

    print(f"{'det':<6}{'spat_id':>8}{'maskdef_id':>11}  {'status':<6} reason")
    for r in t[t["status"] != "ok"]:
        print(f"{r['det']:<6}{r['spat_id']:>8}{r['maskdef_id']:>11}  {r['status']:<6} {r['reason']}")
    n = {s: int((t["status"] == s).sum()) for s in ("ok", "bad", "fail", "box")}
    print(f"\n{len(t)} slits: {n['ok']} ok, {n['bad']} bad, {n['fail']} no solution, "
          f"{n['box']} alignment boxes.  Table -> {out}")


if __name__ == "__main__":
    main()
