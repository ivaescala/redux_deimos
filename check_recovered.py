#!/usr/bin/env python
"""
check_recovered.py  SPEC1D [SPEC1D ...]  [--maskdef ID [ID ...]]

For each target (design-matched object only) and exposure, print:
  fwhm    optimal-extraction profile FWHM (arcsec, as in the spec1d .txt)
  med_opt median OPT_COUNTS in 8400-8800 A
  med_box median BOX_COUNTS in 8400-8800 A (PypeIt boxcar, 3" wide by default, so it
          includes neighbours within 1.5")
  opt/box ratio of the two. The optimal flux should not exceed everything inside the
          boxcar; well above 1 means the optimal extraction is unreliable.
  sky_r   correlation between target and sky spectra after removing the continuum of
          each (51-pixel running median). This measures shape only: with bright sky
          lines even a small residual drives it towards -1 or +1.
  sky_k   slope of the same relation: the fraction of the sky-line flux left in the
          target spectrum (negative = over-subtracted, positive = under-subtracted).
          This is the size of the residual; -0.05 means 5% of the sky-line flux was
          removed in excess.

Without --maskdef, every design-matched target in the files is measured and a summary
of the distribution is printed at the end: run it on the main reduction's spec1d files
to get the reference values for normal slits.
"""
import argparse
import os
import re

import numpy as np
from scipy.ndimage import median_filter

WIN = (8400.0, 8800.0)


def highpass(y, n=51):
    return y - median_filter(y, size=n, mode="nearest")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("spec1d", nargs="+")
    ap.add_argument("--maskdef", nargs="+", type=int, default=None)
    a = ap.parse_args()

    from pypeit.specobjs import SpecObjs
    rows = []
    for fn in sorted(a.spec1d):
        m = re.search(r"DE\.\d+\.(\d+)", os.path.basename(fn))
        exp = m.group(1) if m else os.path.basename(fn)
        sobjs = SpecObjs.from_fitsfile(fn, chk_version=False)
        for i in range(len(sobjs)):
            o = sobjs[i]
            if str(o.MASKDEF_OBJNAME).strip() == "SERENDIP":
                continue
            if a.maskdef is not None and o.MASKDEF_ID not in a.maskdef:
                continue
            if o.OPT_WAVE is None or o.OPT_COUNTS is None:
                rows.append((o.MASKDEF_ID, exp, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan))
                continue
            w = o.OPT_WAVE
            sel = (w > WIN[0]) & (w < WIN[1]) & np.isfinite(o.OPT_COUNTS)
            opt = np.nanmedian(o.OPT_COUNTS[sel])
            box = np.nanmedian(o.BOX_COUNTS[sel]) if o.BOX_COUNTS is not None else np.nan
            r = k = np.nan
            if o.OPT_COUNTS_SKY is not None and sel.sum() > 100:
                t_hp = highpass(o.OPT_COUNTS[sel])
                s_hp = highpass(o.OPT_COUNTS_SKY[sel])
                r = np.corrcoef(t_hp, s_hp)[0, 1]
                k = np.sum(t_hp * s_hp) / np.sum(s_hp * s_hp)
            fw = o.FWHM * 0.1185 if o.FWHM is not None else np.nan
            rows.append((o.MASKDEF_ID, exp, fw, opt, box, opt / box if box else np.nan, r, k))

    print(f"{'maskdef':>8} {'exp':>6} {'fwhm':>5} {'med_opt':>8} {'med_box':>8} {'opt/box':>7} "
          f"{'sky_r':>6} {'sky_k':>6}")
    ids = a.maskdef if a.maskdef is not None else sorted({r[0] for r in rows})
    for mid in ids:
        sub = [r for r in rows if r[0] == mid]
        if not sub:
            print(f"{mid:>8}  no design-matched object in these spec1d files")
            continue
        for r in sub:
            print(f"{r[0]:>8} {r[1]:>6} {r[2]:5.2f} {r[3]:8.1f} {r[4]:8.1f} {r[5]:7.2f} {r[6]:6.2f} {r[7]:6.3f}")
    if a.maskdef is None and rows:
        arr = np.array([r[2:] for r in rows], float)
        print(f"\nsummary over {len(rows)} target spectra: 10th / 50th / 90th percentile")
        for j, name in enumerate(["fwhm", "med_opt", "med_box", "opt/box", "sky_r", "sky_k"]):
            v = arr[:, j][np.isfinite(arr[:, j])]
            if v.size:
                p10, p50, p90 = np.percentile(v, [10, 50, 90])
                print(f"  {name:>8}: {p10:8.3f} {p50:8.3f} {p90:8.3f}")


if __name__ == "__main__":
    main()
