#!/usr/bin/env python
"""
dmost_prep.py  --cfg CFG  --name NAME  --workroot DIR  --rawdir DIR  [--data DIR] [--for dmost|marz] [--fresh]

Build the flat working folder dmost expects, from the outputs of the `collate` stage
of reduce_mask.sh, without copying or modifying them. Called by the `marz` and
`dmost` stages; it can also be run by hand.

  CFG       configuration folder, e.g. $DEIMOS_RDX/2022B/M32RA1/keck_deimos_A
  NAME      mask name dmost uses in its file names (M32RA1, or M32RA1_B for a 2nd config)
  workroot  $DEIMOS_RDX/<semester>/dmost, passed to dmost as DEIMOS_REDUX
  rawdir    raw folder of the mask, $DEIMOS_RAW/<semester>/<mask>
  data      dmost data root holding templates/ and Other_data/ (only for --for dmost)

Result (W = workroot/NAME):

  W/Science             -> CFG/Science_merged       spec1d files after exclusions/recoveries
  W/collate1d           -> CFG/collate1d            0.5" coadds
  W/collate_report.dat  -> collate1d/collate_report.dat
  W/dmost/ QA/ emcee/                               dmost outputs
  W/dmost_raw/                                      passed to dmost as DEIMOS_RAW:
      rawdata_<year>/DE.*.fits -> raw science frames (dmost reads slit widths from them)
      templates, Other_data    -> DATA/templates, DATA/Other_data
  W/not_in_dmost.csv                                objects PypeIt extracted that dmost will not see
  workroot/marz_files/                              Marz input, and marz_<NAME>_MG.mz if you make one

Checks, before anything is written:
  - Science_merged has spec1d files and collate1d has collate_report.dat and coadds;
  - collate_report.dat refers only to files in Science_merged;
  - (--for dmost) every raw science frame is found, and the template and sky-line
    files dmost reads are present under DATA;
  - existing dmost outputs are not older than the collate outputs. dmost resumes from
    its own table and rebuilds its flexure-corrected coadds only when none exist, so
    stale outputs would silently mix two reductions. --fresh moves old outputs to
    W/archive_<time>/ and starts over.
"""
import argparse
import csv
import datetime
import glob
import os
import shutil
import sys
from pathlib import Path

OUTPUTS = ["dmost", "QA", "emcee", "Science_flex", "collate1d_flex"]   # + <NAME>_dmost.log


def fail(msg):
    sys.exit(f"dmost_prep.py: {msg}")


def link(target, linkpath, relative=True):
    """Create or update a symlink; never replace a real file or folder."""
    linkpath = Path(linkpath)
    dest = os.path.relpath(target, linkpath.parent) if relative else str(target)
    if linkpath.is_symlink():
        if os.readlink(linkpath) == dest:
            return
        linkpath.unlink()
    elif linkpath.exists():
        fail(f"{linkpath} exists and is not a symlink; move it away first")
    linkpath.symlink_to(dest)


def check_data(data):
    """Files dmost opens under DEIMOS_RAW (besides the raw frames). Returns a list of problems."""
    need = {
        "templates/pheonix/grid1/dmost*.fits": "stellar templates, S/N >= 25",
        "templates/pheonix/grid2/dmost*.fits": "stellar templates, 10 < S/N < 25",
        "templates/pheonix/grid3/dmost*.fits": "stellar templates, S/N <= 10",
        "templates/tellurics/telluric_0.02A*fits": "coarse telluric grid",
        "templates/fine_tellurics/telluric_0.02A_h2o_*fits": "fine telluric grid",
        "Other_data/sky_single_mg.dat": "sky emission lines for flexure",
    }
    return [f"  {data}/{pat}   ({what})" for pat, what in need.items() if not glob.glob(str(data / pat))]


def jyear(mjd):
    """Year exactly as dmost_create_maskfile.parse_year computes it."""
    from astropy.time import Time
    return Time(mjd, format="mjd").to_value("jyear", subfmt="str").split(".")[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cfg", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--workroot", required=True)
    ap.add_argument("--rawdir", required=True)
    ap.add_argument("--data", default=None)
    ap.add_argument("--for", dest="purpose", choices=["dmost", "marz"], default="dmost")
    ap.add_argument("--fresh", action="store_true", help="archive existing dmost outputs and start over")
    a = ap.parse_args()

    from astropy.io import ascii, fits

    cfg = Path(a.cfg).resolve()
    rawdir = Path(a.rawdir).resolve()
    workroot = Path(a.workroot).resolve()
    W = workroot / a.name

    # ---- collate outputs ----------------------------------------------------------
    merged = sorted((cfg / "Science_merged").glob("spec1d_*.fits"))
    report = cfg / "collate1d" / "collate_report.dat"
    coadds = list((cfg / "collate1d").glob("J*.fits"))
    if not merged:
        fail(f"no spec1d files in {cfg / 'Science_merged'}; run the collate stage first")
    if not report.is_file() or len(coadds) < 2:
        fail(f"no collate_report.dat or coadds in {cfg / 'collate1d'}; run the collate stage first")
    rep = ascii.read(report)
    merged_names = {p.name for p in merged}
    foreign = sorted(set(rep["spec1d_filename"]) - merged_names)
    if foreign:
        fail(f"{report} lists spec1d files that are not in Science_merged ({', '.join(foreign[:3])} ...);"
             " re-run the collate stage")

    # ---- raw frames and data files (dmost only) -----------------------------------
    raw_links, problems = [], []
    if a.purpose == "dmost":
        if not a.data:
            fail("--data (dmost data root with templates/ and Other_data/) is required for dmost")
        data = Path(a.data).expanduser().resolve()
        missing = check_data(data)
        if missing:
            problems.append("dmost data files missing (set DMOST_DATA; see README section 9):\n"
                            + "\n".join(missing))
        for f in merged:
            hdr = fits.getheader(f, 0)
            fn, mjd = hdr.get("FILENAME"), hdr.get("MJD")
            if not fn or mjd is None:
                problems.append(f"  {f.name}: no FILENAME or MJD in the primary header")
                continue
            src = rawdir / fn
            if not src.is_file():
                hits = sorted(rawdir.rglob(fn))
                src = hits[0] if hits else None
            if src is None:
                problems.append(f"  raw frame {fn} (for {f.name}) not found under {rawdir}")
                continue
            raw_links.append((f"rawdata_{jyear(mjd)}", src))
    if problems:
        fail("cannot run dmost:\n" + "\n".join(problems))

    # ---- stale or old outputs -----------------------------------------------------
    newest_input = max(p.stat().st_mtime for p in merged + [report])
    existing = [W / d for d in OUTPUTS if (W / d).exists()] + \
               [p for p in [W / f"{a.name}_dmost.log", W / "not_in_dmost.csv"] if p.exists()]
    table = W / "dmost" / f"dmost_{a.name}.fits"
    flexdir = W / "collate1d_flex"
    stale = [p for p in (table, flexdir) if p.exists() and p.stat().st_mtime < newest_input]
    if existing and (a.fresh or stale):
        if stale and not a.fresh:
            fail(f"dmost outputs in {W} predate the current collate outputs "
                 f"({', '.join(p.name for p in stale)}). Re-run with DMOST_FRESH=1 to archive them and start over.")
        arch = W / f"archive_{datetime.datetime.now():%Y%m%d-%H%M%S}"
        arch.mkdir(parents=True)
        for p in existing:
            shutil.move(str(p), str(arch / p.name))
        print(f"moved previous dmost outputs to {arch}")

    # ---- build the working folder -------------------------------------------------
    for d in ["dmost", "QA", "emcee"]:
        (W / d).mkdir(parents=True, exist_ok=True)
    (workroot / "marz_files").mkdir(parents=True, exist_ok=True)
    link(cfg / "Science_merged", W / "Science")
    link(cfg / "collate1d", W / "collate1d")
    link(W / "collate1d" / "collate_report.dat", W / "collate_report.dat")

    if a.purpose == "dmost":
        env = W / "dmost_raw"
        env.mkdir(exist_ok=True)
        link(data / "templates", env / "templates", relative=False)
        link(data / "Other_data", env / "Other_data", relative=False)
        for year, src in raw_links:
            (env / year).mkdir(exist_ok=True)
            link(src, env / year / src.name, relative=False)

    # ---- objects dmost will not see -----------------------------------------------
    # dmost builds its slit list from collate_report.dat. PypeIt 2.0.1 collate writes no
    # coadd, and no report row, for an object with a single spectrum (e.g. a slit
    # extracted in only one exposure), so such objects never reach dmost.
    in_report = {(r["spec1d_filename"], str(r["pypeit_name"]).strip()) for r in rep}
    dropped = []
    for f in merged:
        with fits.open(f) as hdul:
            for h in hdul[1:]:
                if not h.name.startswith("SPAT"):
                    continue
                if (f.name, h.name) not in in_report:
                    dropped.append({"spec1d": f.name, "pypeit_name": h.name,
                                    "maskdef_id": h.header.get("MASKDEF_ID", ""),
                                    "objname": h.header.get("MASKDEF_OBJNAME", "")})
    with open(W / "not_in_dmost.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=["spec1d", "pypeit_name", "maskdef_id", "objname"])
        wr.writeheader()
        wr.writerows(dropped)

    n_obj = len({str(x) for x in rep["filename"]})
    n_targ = len({d["maskdef_id"] for d in dropped if d["objname"] not in ("SERENDIP", "")})
    n_ser = sum(d["objname"] == "SERENDIP" for d in dropped)
    print(f"{a.name}: {len(merged)} exposures, {n_obj} coadded objects -> {W}")
    if a.purpose == "dmost":
        years = sorted({y for y, _ in raw_links})
        print(f"raw science frames linked into {', '.join(years)} under {W / 'dmost_raw'}")
    if dropped:
        print(f"NOTE: {len(dropped)} extracted spectra have no coadd and are not measured by dmost "
              f"({n_targ} targets, {n_ser} serendipitous spectra): see {W / 'not_in_dmost.csv'}")


if __name__ == "__main__":
    main()
