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

The arc solutions are then tested against the science frames, which the checks
above cannot do: a solution can be clean for the arc and still wrong for the science
(M32RA2, 2022-09-28: science wavelengths ~+32.5 A too red).

  - sky:  every science frame is processed minimally (overscan, trim, mosaic); for
    each "ok" slit a sky spectrum is taken along the middle 60% of the slit, given
    the arc wavelengths, and cross-correlated over +-SKY_SHMAX A, blue and red half
    separately, with PypeIt's sky_single_mg.dat lines (vacuum, like the solutions).
    A frame is flagged OFFSET when the median shift exceeds SKY_TOL. The wide range
    matters: OH lines repeat at ~30-35 A, so a shifted spectrum still matches lines
    locally. This is the test that decides.
  - tilt: the raw grating-tilt encoder (G<n>TLTRAW) of every science frame against
    the arc of its calibration group, noted as TILT above TILT_TOL steps. For
    information only: M32RD2 has the same difference as RA2 (~197 steps) and no offset.

Writes <folder>/wavecal_sky.csv (shift per science frame, slit and half) and prints
one summary line per science frame. --no-sky skips the sky test (~25 s per mosaic).
"""
import argparse
import glob
import os
import re
import warnings

import numpy as np
from astropy.io import ascii, fits
from astropy.table import Table, vstack
from scipy.ndimage import median_filter
from scipy.signal import correlate

RMS_MAX = 0.4      # pixels
NLIN_MIN = 20      # arc lines used in the fit
COV_MIN = 60.0     # percent of the slit spectrum spanned by fitted lines
DISP_TOL = 0.05    # fractional deviation of dWave from the detector median

TILT_TOL = 100     # encoder steps between science and arc; informative only (RA2 198 steps +32.5 A, RD2 197 steps -0.2 A)
SKY_TOL = 2.0      # A, median sky-line shift of a science frame under the arc solution
SKY_SHMAX = 60.0   # A, search range of the sky cross-correlation
SKY_STEP = 0.05    # A
SKY_SPLIT = 7800.  # A, boundary between the blue and red halves

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


def frames(config_dir):
    """(path, filename, frametype, calib groups) of every frame in the .pypeit file."""
    from pypeit.inputfiles import PypeItFile
    pf = PypeItFile.from_file(os.path.join(config_dir, os.path.basename(os.path.normpath(config_dir)) + ".pypeit"))
    out = []
    for r in pf.data:
        path = next(p for p in pf.file_paths if os.path.exists(os.path.join(p, r["filename"])))
        cal = str(r["calib"]).strip()
        out.append((path, r["filename"], r["frametype"], None if cal == "all" else set(cal.split(","))))
    return out


def tilt_raw(file):
    h = fits.getheader(file)
    pos = h.get("GRATEPOS")
    return h.get(f"G{pos}TLTRAW") if pos else None


def sky_shift(wave, sky, lines):
    """Shift (A) of the sky lines in `sky` relative to `lines`, by cross-correlation."""
    o = np.argsort(wave)
    wave, sky = wave[o], sky[o]
    grid = np.arange(wave[0], wave[-1], SKY_STEP)
    sub = np.interp(grid, wave, sky)
    sub = np.clip(sub - median_filter(sub, int(5 / SKY_STEP) | 1), 0, None)
    if not np.any(sub > 0):
        return np.nan
    n = int(round(SKY_SHMAX / SKY_STEP))
    ext = np.arange(grid[0] - n * SKY_STEP, grid[-1] + (n + 0.5) * SKY_STEP, SKY_STEP)[: grid.size + 2 * n]
    model = np.zeros_like(ext)
    for lam in lines[(lines > ext[0] - 5) & (lines < ext[-1] + 5)]:
        model += np.exp(-0.5 * ((ext - lam) / 0.7) ** 2)
    cc = correlate(model, sub, mode="valid", method="fft")
    return -(int(np.argmax(cc)) - n) * SKY_STEP


def science_checks(config_dir, qa, do_sky):
    """Tilt and sky-line tests of every science frame against its arc solution."""
    from pypeit import dataPaths
    from pypeit.images.buildimage import buildimage_fromlist
    from pypeit.par.pypeitpar import FrameGroupPar, ProcessImagesPar
    from pypeit.slittrace import SlitTraceSet
    from pypeit.spectrographs.util import load_spectrograph
    from pypeit.wavecalib import WaveCalib

    caldir = os.path.join(config_dir, "Calibrations")
    fr = frames(config_dir)
    arcs = [(os.path.join(p, f), cal) for p, f, ft, cal in fr if "arc" in ft]
    sci = [(os.path.join(p, f), cal) for p, f, ft, cal in fr if "science" in ft]
    if not sci:
        print("\nNo science frames in the .pypeit file: tilt and sky tests skipped.")
        return
    lines = np.array(ascii.read(dataPaths.sky_spec.get_file_path("sky_single_mg.dat"))["Wave"], float)
    spec = load_spectrograph("keck_deimos")
    proc = FrameGroupPar(frametype="science", process=ProcessImagesPar(
        use_biasimage=False, use_darkimage=False, use_pixelflat=False, use_illumflat=False,
        use_specillum=False, use_pattern=False, mask_cr=False, use_overscan=True, apply_gain=True,
        shot_noise=False, noise_floor=0.0))
    ok = {(r["calib_key"], int(r["spat_id"])): int(r["maskdef_id"]) for r in qa if r["status"] == "ok"}

    rows, summary = [], []
    for file, cal in sci:
        exp = re.search(r"DE\.\d+\.(\d+)\.", file).group(1)
        group = [a for a, c in arcs if cal is None or c is None or cal & c]
        t_sci, t_arc = tilt_raw(file), (tilt_raw(group[0]) if group else None)
        dtilt = None if t_sci is None or t_arc is None else t_sci - t_arc
        shifts = []
        if do_sky:
            for wfile in sorted(glob.glob(os.path.join(caldir, "WaveCalib_*.fits"))):
                key = os.path.basename(wfile)[len("WaveCalib_"):-len(".fits")]     # e.g. A_0_MSC01
                cid, det = key.split("_")[1], key.split("_")[-1]
                if cal is not None and cid not in cal:
                    continue
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    img = buildimage_fromlist(spec, spec.allowed_mosaics[int(det[-1]) - 1], proc, [file]).image
                wv = WaveCalib.from_file(wfile, chk_version=False)
                slits = SlitTraceSet.from_file(os.path.join(caldir, f"Slits_{key}.fits.gz"), chk_version=False)
                left, right, _ = slits.select_edges()
                for k, sid in enumerate(slits.spat_id):
                    j = np.where(wv.spat_ids == sid)[0]
                    if (key, int(sid)) not in ok or not j.size or wv.wv_fits[j[0]] is None:
                        continue
                    wave = np.asarray(wv.wv_fits[j[0]].wave_soln, float)
                    lo = np.ceil(left[:, k] + 0.2 * (right[:, k] - left[:, k])).astype(int)
                    hi = np.floor(right[:, k] - 0.2 * (right[:, k] - left[:, k])).astype(int)
                    sky = np.array([np.median(img[r, max(lo[r], 0):max(hi[r], lo[r] + 1)])
                                    for r in range(img.shape[0])])
                    good = np.isfinite(sky) & (wave > 0)
                    for half, sel in (("blue", good & (wave < SKY_SPLIT)), ("red", good & (wave >= SKY_SPLIT))):
                        if sel.sum() > 300:
                            d = sky_shift(wave[sel], sky[sel], lines)
                            rows.append((exp, det, int(sid), ok[(key, int(sid))], half, d))
                            shifts.append(d)
        d = np.array(shifts, float)
        d = d[np.isfinite(d)]
        med = np.median(d) if d.size else np.nan
        agree = np.mean(np.abs(d - med) < SKY_TOL) * 100 if d.size else np.nan
        flags = []
        if dtilt is not None and abs(dtilt) > TILT_TOL:
            flags.append("TILT")
        if np.isfinite(med) and abs(med) > SKY_TOL:
            flags.append("OFFSET")
        summary.append((exp, t_arc, t_sci, dtilt, med, agree, d.size, " ".join(flags)))

    if rows:
        t = Table(rows=rows, names=("exp", "det", "spat_id", "maskdef_id", "half", "shift_A"))
        t.write(os.path.join(config_dir, "wavecal_sky.csv"), overwrite=True)
    print(f"\nScience frames vs arc solution (tilt tolerance {TILT_TOL} steps, offset tolerance {SKY_TOL} A):")
    print(f"{'exp':>6} {'arc tilt':>9} {'sci tilt':>9} {'diff':>6} {'sky shift':>10} {'agree':>6} {'N':>5}  flags")
    for exp, ta, ts, dt, med, ag, n, fl in summary:
        f = lambda v, fmt, unit="": format(v, fmt) + unit if v is not None and np.isfinite(v) else "-"
        print(f"{exp:>6} {f(ta, 'd'):>9} {f(ts, 'd'):>9} {f(dt, '+d'):>6} {f(med, '+.2f', ' A'):>10} "
              f"{f(ag, '.0f', '%'):>6} {n:>5}  {fl}")
    bad = [s for s in summary if "OFFSET" in s[7]]
    tilt = [s for s in summary if "TILT" in s[7]]
    if bad:
        print(f"\nWARNING: {len(bad)} of {len(summary)} science frames do not match their arc solution "
              f"({', '.join(s[0] for s in bad)}): their wavelengths will be offset by the sky shift "
              "shown. Do not run the science stage before deciding how to handle them.")
    elif do_sky:
        print(f"All {len(summary)} science frames match their arc solution"
              + (" (the TILT difference did not shift the wavelengths)." if tilt else "."))
    elif tilt:
        print(f"\nNOTE: {len(tilt)} science frames differ in grating tilt from their arc; the tilt alone "
              "does not decide (RD2: no offset, RA2: +32.5 A). Run without --no-sky to measure the offset.")


def main():
    from pypeit.slittrace import SlitTraceSet
    from pypeit.wavecalib import WaveCalib

    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    p.add_argument("config_dir")
    p.add_argument("--no-sky", action="store_true", help="skip the sky-line test of the science frames")
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

    science_checks(a.config_dir, t, not a.no_sky)


if __name__ == "__main__":
    main()
