#!/usr/bin/env python
"""
koa_deimos_fetch.py

Find DEIMOS science frames around a target in a given Keck semester on the
Keck Observatory Archive (KOA), pull the matching slitmask calibrations
(flats + arcs taken through the same mask/grating on the same nights), download
the raw files, and sort them into one folder per slitmask, ready for PypeIt.

Defaults: M32, semester 2022B.

Output layout
-------------
  <outdir>/<slitmask>/   science + flats + arcs taken through that mask
                         (+ biases from the nights that mask was observed, with --biases)
  <outdir>/meta/         query tables and download manifest
  <outdir>/lev0/         PyKOA's staging area; emptied and removed after sorting

Biases are hard-linked (no extra disk space) into every mask folder observed
on the bias's night; a bias with no matching mask goes to <outdir>/bias/.

Typical use
-----------
  # 1. Inspect what exists (no download)
  python koa_deimos_fetch.py --outdir ~/deimos/2022B --list-only

  # 2. Download + sort (science + flats + arcs)
  python koa_deimos_fetch.py --outdir ~/deimos/2022B

  # 2b. Restrict to specific masks, add biases
  python koa_deimos_fetch.py --outdir ~/deimos/2022B --masks M32maskA M32maskB --biases

  # 3. Already downloaded into <outdir>/lev0? Sort without querying KOA again
  python koa_deimos_fetch.py --outdir ~/deimos/2022B --organize-only

  # Proprietary data: log in once, then reuse the cookie
  python koa_deimos_fetch.py --outdir ... --login --cookie ~/.koa_cookie.txt

Re-running is safe: files already anywhere under <outdir> are not downloaded
again, and truncated downloads are deleted so the next run re-fetches them.

Requires: pykoa, astropy  (pip install pykoa astropy)
"""

import argparse
import datetime as dt
import glob
import os
import re
import shutil
import sys
from collections import defaultdict

from astropy.io import fits
from astropy.table import Table, unique, vstack
from pykoa.koa import Koa

# M32, J2000 (deg)
M32_RA, M32_DEC = 10.674270, 40.865170

COLS = ("koaid, filehand, targname, koaimtyp, "
        "to_char(date_obs,'YYYY-MM-DD') as date_obs, ut, ra, dec, elaptime, "
        "slmsknam, gratenam, obsmode, semid, progpi, progtitl")

TYPE_LABEL = {"object": "sci", "flatlamp": "flat", "arclamp": "arc", "bias": "bias"}


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def strcol(tbl, name):
    """Column as a list of stripped strings (masked -> '')."""
    col = tbl[name]
    if hasattr(col, "filled"):
        col = col.filled("")
    return [str(v).strip() for v in col]


def run_adql(query, outpath, cookie=None):
    kw = dict(overwrite=True, format="ipac")
    if cookie:
        kw["cookiepath"] = cookie
    Koa.query_adql(query, outpath, **kw)
    if not os.path.exists(outpath):
        sys.exit(f"[error] KOA returned no table for query:\n{query}")
    try:
        return Table.read(outpath, format="ascii.ipac")
    except Exception as e:  # empty or error table
        print(f"[warn] could not parse {outpath} ({e}); treating as empty")
        return Table()


def sql_list(values):
    return ", ".join("'" + v.replace("'", "''") + "'" for v in values)


def near_night(date_str, nights, pad_days):
    try:
        d = dt.date.fromisoformat(date_str)
    except ValueError:
        return False
    return any(abs((d - n).days) <= pad_days for n in nights)


def mask_dirname(name):
    """Slitmask name -> safe folder name."""
    s = re.sub(r"[^A-Za-z0-9._+-]+", "_", name.strip()).strip("._")
    return s or "unknown_mask"


def summarize(sci):
    """Print one line per (mask, grating) configuration."""
    masks, grats = strcol(sci, "slmsknam"), strcol(sci, "gratenam")
    targs, dates = strcol(sci, "targname"), strcol(sci, "date_obs")
    semids, pis = strcol(sci, "semid"), strcol(sci, "progpi")
    exptime = [float(x) for x in sci["elaptime"]]

    groups = {}
    for i in range(len(sci)):
        g = groups.setdefault((masks[i], grats[i]), dict(
            targ=set(), nights=set(), n=0, t=0.0, semid=set(), pi=set()))
        g["targ"].add(targs[i]); g["nights"].add(dates[i])
        g["n"] += 1; g["t"] += exptime[i]
        g["semid"].add(semids[i]); g["pi"].add(pis[i])

    print(f"\n{'slitmask':<18}{'grating':<10}{'targname':<20}{'nexp':>5}"
          f"{'t_exp[ks]':>10}  {'nights (UT)':<34}{'semid / PI'}")
    print("-" * 125)
    for (m, gr), g in sorted(groups.items()):
        print(f"{m:<18}{gr:<10}{','.join(sorted(g['targ']))[:19]:<20}"
              f"{g['n']:>5}{g['t'] / 1e3:>10.1f}  "
              f"{','.join(sorted(g['nights']))[:33]:<34}"
              f"{','.join(sorted(g['semid']))} / {','.join(sorted(g['pi']))}")
    print()


def verify_fits(lev0dir):
    """Delete truncated/corrupt FITS so a re-run downloads them again."""
    bad = []
    for f in sorted(glob.glob(os.path.join(lev0dir, "D*.fits*"))):
        try:
            with fits.open(f, lazy_load_hdus=False) as hl:
                expected = sum(h.filebytes() for h in hl)
            if expected > os.path.getsize(f):
                raise IOError("short file")
        except Exception:
            bad.append(f)
            os.remove(f)
    return bad


def find_raw(outdir):
    """KOAID -> list of paths for every DEIMOS raw file under outdir (meta/ excluded)."""
    found = defaultdict(list)
    for root, dirs, files in os.walk(outdir):
        dirs[:] = [d for d in dirs if d != "meta"]
        for f in files:
            if f.startswith("DE.") and ".fits" in f:
                found[f].append(os.path.join(root, f))
    return found


def place(src, dest_dirs):
    """Move src into dest_dirs[0]; hard-link (or copy) into the rest."""
    for d in dest_dirs:
        os.makedirs(d, exist_ok=True)
        dst = os.path.join(d, os.path.basename(src))
        if os.path.exists(dst):
            continue
        try:
            os.link(src, dst)
        except OSError:
            shutil.copy2(src, dst)
    os.remove(src)


def organize(manifest, outdir):
    """Sort files from <outdir>/lev0 into <outdir>/<slitmask>/ using the manifest."""
    lev0 = os.path.join(outdir, "lev0")
    koaid, typ, mask, date = (strcol(manifest, c)
                              for c in ("koaid", "koaimtyp", "slmsknam", "date_obs"))

    sci_nights = defaultdict(set)          # mask -> UT dates with science
    for t, m, d in zip(typ, mask, date):
        if t == "object":
            sci_nights[m].add(d)

    for k, t, m, d in zip(koaid, typ, mask, date):
        src = os.path.join(lev0, k)
        if not os.path.exists(src):
            continue
        if t == "bias":
            names = sorted(mask_dirname(mm) for mm, n in sci_nights.items() if d in n) or ["bias"]
        else:
            names = [mask_dirname(m)]
        place(src, [os.path.join(outdir, n) for n in names])

    leftovers = glob.glob(os.path.join(lev0, "*"))
    if os.path.isdir(lev0) and not leftovers:
        os.rmdir(lev0)
    elif leftovers:
        print(f"[warn] {len(leftovers)} files in {lev0} are not in the manifest; left in place.")

    # per-folder summary
    lookup = dict(zip(koaid, typ))
    counts = defaultdict(lambda: defaultdict(int))
    for k, paths in find_raw(outdir).items():
        for p in paths:
            folder = os.path.relpath(os.path.dirname(p), outdir)
            counts[folder][TYPE_LABEL.get(lookup.get(k, ""), "other")] += 1
    print(f"\n{'folder':<28}{'sci':>6}{'flat':>6}{'arc':>6}{'bias':>6}{'other':>7}")
    print("-" * 59)
    for folder in sorted(counts):
        c = counts[folder]
        print(f"{folder:<28}{c['sci']:>6}{c['flat']:>6}{c['arc']:>6}{c['bias']:>6}{c['other']:>7}")


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def build_manifest(a, meta, sem):
    """Query KOA; return the merged table of science + calibration files (or None)."""
    # ---- 1. science frames -------------------------------------------------
    q_sci = (f"select {COLS} from koa_deimos where "
             f"contains(point('J2000',ra,dec), circle('J2000',{a.ra},{a.dec},{a.radius}))=1 "
             f"and koaimtyp = 'object' "
             f"and (semid like '{sem}%' or semid like '{sem.upper()}%')")
    sci = run_adql(q_sci, os.path.join(meta, "science_raw_query.tbl"), a.cookie)
    if len(sci) == 0:
        sys.exit("No DEIMOS science frames found. Try a larger --radius or check --semester.")

    keep = [k.startswith("DE.") for k in strcol(sci, "koaid")]   # drop FCS (DF.) frames
    if a.target_regex:
        rx = re.compile(a.target_regex, re.I)
        keep = [k and bool(rx.search(t)) for k, t in zip(keep, strcol(sci, "targname"))]
    if a.masks:
        keep = [k and m in a.masks for k, m in zip(keep, strcol(sci, "slmsknam"))]
    sci = sci[keep]
    if len(sci) == 0:
        sys.exit("All science frames removed by --target-regex/--masks.")

    print(f"{len(sci)} science frames within {a.radius} deg of "
          f"({a.ra:.5f}, {a.dec:+.5f}) in {a.semester}:")
    summarize(sci)
    sci.write(os.path.join(meta, "science.tbl"), format="ascii.ipac", overwrite=True)
    if a.list_only:
        print(f"Science table written to {meta}/science.tbl. Re-run without --list-only "
              "(optionally with --masks / --target-regex) to download.")
        return None

    # ---- 2. calibrations through the same masks ----------------------------
    masks, grats, dates = strcol(sci, "slmsknam"), strcol(sci, "gratenam"), strcol(sci, "date_obs")
    cfg = {}   # mask -> {'grat': set, 'nights': set}
    for m, g, d in zip(masks, grats, dates):
        c = cfg.setdefault(m, dict(grat=set(), nights=set()))
        c["grat"].add(g); c["nights"].add(dt.date.fromisoformat(d))

    q_cal = (f"select {COLS} from koa_deimos where "
             f"slmsknam in ({sql_list(cfg)}) "
             f"and koaimtyp in ('flatlamp','arclamp')")
    cal = run_adql(q_cal, os.path.join(meta, "calib_raw_query.tbl"), a.cookie)

    pieces = [sci]
    if len(cal):
        cm, cg, cd, ck = (strcol(cal, c) for c in ("slmsknam", "gratenam", "date_obs", "koaid"))
        ok = [k.startswith("DE.") and m in cfg and g in cfg[m]["grat"]
              and near_night(d, cfg[m]["nights"], a.pad_days)
              for m, g, d, k in zip(cm, cg, cd, ck)]
        cal = cal[ok]
        pieces.append(cal)
    print(f"{len(cal)} matching flats/arcs "
          f"({sum(t == 'flatlamp' for t in strcol(cal, 'koaimtyp')) if len(cal) else 0} flats).")

    # ---- 3. optional biases --------------------------------------------------
    if a.biases:
        all_nights = set().union(*(c["nights"] for c in cfg.values()))
        q_b = (f"select {COLS} from koa_deimos where koaimtyp = 'bias' "
               f"and (semid like '{sem}%' or semid like '{sem.upper()}%')")
        bia = run_adql(q_b, os.path.join(meta, "bias_raw_query.tbl"), a.cookie)
        if len(bia):
            ok = [k.startswith("DE.") and near_night(d, all_nights, 0)
                  for k, d in zip(strcol(bia, "koaid"), strcol(bia, "date_obs"))]
            bia = bia[ok]
            pieces.append(bia)
        print(f"{len(bia)} bias frames from the science nights.")

    allf = unique(vstack(pieces, metadata_conflicts="silent"), keys="koaid")
    allf["instrume"] = "DEIMOS"          # pykoa.download needs instrume/koaid/filehand
    allf.sort("koaid")
    return allf


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", required=True,
                   help="root folder; files are sorted into outdir/<slitmask>/")
    p.add_argument("--ra", type=float, default=M32_RA)
    p.add_argument("--dec", type=float, default=M32_DEC)
    p.add_argument("--radius", type=float, default=0.10,
                   help="cone radius on telescope pointing, deg (default 0.10)")
    p.add_argument("--semester", default="2022B")
    p.add_argument("--target-regex", default=None,
                   help="keep only science frames whose TARGNAME matches (case-insensitive)")
    p.add_argument("--masks", nargs="+", default=None,
                   help="keep only these slitmask names (SLMSKNAM)")
    p.add_argument("--biases", action="store_true",
                   help="also fetch bias frames from the science nights")
    p.add_argument("--pad-days", type=int, default=1,
                   help="calibs accepted within +/- this many UT days of a science night")
    p.add_argument("--list-only", action="store_true", help="query + summarize, no download")
    p.add_argument("--organize-only", action="store_true",
                   help="skip KOA; sort files already in outdir/lev0 using meta/download_manifest.tbl")
    p.add_argument("--cookie", default=None, help="KOA cookie file (proprietary data)")
    p.add_argument("--login", action="store_true", help="prompt for KOA login, write --cookie")
    a = p.parse_args()

    # Pre-flight: fail with a readable message if outdir can't be created/written.
    a.outdir = os.path.abspath(os.path.expanduser(a.outdir))
    probe = a.outdir
    while not os.path.exists(probe):
        probe = os.path.dirname(probe)
    if not os.access(probe, os.W_OK | os.X_OK):
        sys.exit(f"[error] cannot create/write {a.outdir}: no write permission on "
                 f"existing parent {probe}. Choose an --outdir under a directory you own "
                 f"(e.g. your scratch/project space).")

    os.makedirs(a.outdir, exist_ok=True)
    meta = os.path.join(a.outdir, "meta")
    os.makedirs(meta, exist_ok=True)
    lev0 = os.path.join(a.outdir, "lev0")
    manifest = os.path.join(meta, "download_manifest.tbl")

    if a.organize_only:
        if not os.path.exists(manifest):
            sys.exit(f"[error] {manifest} not found; run once without --organize-only.")
        allf = Table.read(manifest, format="ascii.ipac")
        bad = verify_fits(lev0)
    else:
        if a.login:
            if not a.cookie:
                sys.exit("--login needs --cookie <path>")
            Koa.login(a.cookie)

        allf = build_manifest(a, meta, a.semester.lower())
        if allf is None:          # --list-only
            return
        allf.write(manifest, format="ascii.ipac", overwrite=True)

        verify_fits(lev0)         # clear truncated leftovers from an interrupted run
        present = find_raw(a.outdir)
        todo = allf[[k not in present for k in strcol(allf, "koaid")]]
        print(f"\n{len(allf)} files in manifest ({manifest}); "
              f"{len(allf) - len(todo)} already on disk, {len(todo)} to download.")
        if len(todo):
            todo_path = os.path.join(meta, "download_todo.tbl")
            todo.write(todo_path, format="ascii.ipac", overwrite=True)
            dl_kw = dict(cookiepath=a.cookie) if a.cookie else {}
            Koa.download(todo_path, "ipac", a.outdir, **dl_kw)
        bad = verify_fits(lev0)

    organize(allf, a.outdir)

    if bad:
        print(f"\n{len(bad)} truncated downloads removed -- re-run the same command to re-fetch:")
        for b in bad:
            print("   ", os.path.basename(b))
    missing = set(strcol(allf, "koaid")) - set(find_raw(a.outdir))
    if missing:
        print(f"\n{len(missing)} files in the manifest are not on disk "
              "(proprietary without --cookie, or download errors). Re-run to retry.")


if __name__ == "__main__":
    main()
