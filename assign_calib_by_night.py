#!/usr/bin/env python
"""
assign_calib_by_night.py  <file.pypeit>

Give each UT night its own PypeIt calibration group, so science frames are
calibrated with flats/arcs from their own night. Prints a per-night census and
flags nights with science but fewer than 3 flats or no arc (Geha+2026 criteria).
Rewrites the .pypeit file in place (a .bak copy is kept).
"""
import shutil
import sys
from collections import Counter

from pypeit.inputfiles import PypeItFile

fn = sys.argv[1]
shutil.copy(fn, fn + ".bak")
pf = PypeItFile.from_file(fn, preserve_comments=True)
tbl = pf.data

night = [str(d)[:10] for d in tbl["dateobs"]]
gid = {n: str(i) for i, n in enumerate(sorted(set(night)))}
tbl["calib"] = [gid[n] for n in night]

print(f"{'UT night':<12}{'calib':>6}{'sci':>5}{'flat':>6}{'arc':>5}")
for n, g in gid.items():
    c = Counter()
    for ft, nn in zip(tbl["frametype"], night):
        if nn != n:
            continue
        ft = str(ft)
        c["sci"] += "science" in ft
        c["flat"] += "pixelflat" in ft
        c["arc"] += "arc" in ft
    flag = "  <-- incomplete calibs" if c["sci"] and (c["flat"] < 3 or c["arc"] < 1) else ""
    print(f"{n:<12}{g:>6}{c['sci']:>5}{c['flat']:>6}{c['arc']:>5}{flag}")

pf.write(fn)
