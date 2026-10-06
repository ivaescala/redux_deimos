#!/usr/bin/env python
"""
merge_spec1d.py  --main DIR  --test DIR  --select CSV  --out DIR

Build one spec1d file per exposure for pypeit_collate_1d, taking each slit's objects
(the target and any SERENDIP on that slit, which carry the slit's maskdef_id) from the
run named in the selection table:

  maskdef_id,exp,source[,reason]
  1201457,40287,test,recovered with max_mask_frac=0.97
  1201475,39022,none,wrong wavelength solution

  source = main   keep the main reduction's objects for this slit and exposure
           test   take them from the test reduction instead
           none   drop every object on this slit for this exposure

Slit-exposures not listed default to `main`. `exp` is the exposure number as in the
file names (e.g. 39022 in DE.20220927.39022.70). The reason column is ignored.

Duplicate objects (same NAME twice in one file, as in the main run's first exposure)
are reduced to one copy, and the number removed is reported.

The output files keep the main reduction's file names and primary-header keywords, and
a matching .txt summary is written next to each. The input files are not modified.
"""
import argparse
import csv
import glob
import os
import re


def exp_tag(fn):
    m = re.search(r"DE\.\d+\.(\d+)", os.path.basename(fn))
    return m.group(1) if m else None


def read_selection(path):
    sel = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            src = row["source"].strip().lower()
            if src not in ("main", "test", "none"):
                raise ValueError(f"bad source '{row['source']}' in {path}")
            key = (int(row["maskdef_id"]), row["exp"].strip())
            if key in sel and sel[key] != src:
                raise ValueError(f"conflicting entries for maskdef_id {key[0]}, exposure {key[1]}")
            sel[key] = src
    return sel


SKIP = {"SIMPLE", "BITPIX", "NAXIS", "EXTEND", "COMMENT", "HISTORY", "NSPEC", ""}


def subheader(hdr):
    out = {}
    for k in hdr.keys():
        if k in SKIP or k.startswith(("EXT", "CLBS_", "DMOD", "NAXIS")):
            continue
        out[k] = hdr[k]
    return out


def mdid(sobj):
    v = getattr(sobj, "MASKDEF_ID", None)
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--main", required=True, help="Science folder of the main reduction")
    ap.add_argument("--test", required=True, help="Science folder of the test reduction")
    ap.add_argument("--select", required=True, help="selection CSV")
    ap.add_argument("--out", required=True, help="output folder (must not be either input folder)")
    a = ap.parse_args()

    if os.path.abspath(a.out) in (os.path.abspath(a.main), os.path.abspath(a.test)):
        raise SystemExit("--out must be a new folder, not one of the input folders")
    from pypeit.specobjs import SpecObjs

    sel = read_selection(a.select)
    main_files = {exp_tag(f): f for f in sorted(glob.glob(os.path.join(a.main, "spec1d_*.fits")))}
    test_files = {exp_tag(f): f for f in sorted(glob.glob(os.path.join(a.test, "spec1d_*.fits")))}
    if not main_files:
        raise SystemExit(f"no spec1d files in {a.main}")
    unknown = {e for _, e in sel} - set(main_files)
    if unknown:
        raise SystemExit(f"selection lists exposures with no main spec1d file: {sorted(unknown)}")
    os.makedirs(a.out, exist_ok=True)
    runs_used = {}   # maskdef_id -> set of runs its kept objects came from

    for exp, fmain in main_files.items():
        smain = SpecObjs.from_fitsfile(fmain, chk_version=False)
        stest = SpecObjs.from_fitsfile(test_files[exp], chk_version=False) if exp in test_files else None
        need_test = any(src == "test" and e == exp for (_, e), src in sel.items())
        if need_test and stest is None:
            raise SystemExit(f"selection asks for test objects in exposure {exp}, but there is no test spec1d")

        out = SpecObjs(header=smain.header)
        out.calibs = smain.calibs
        read, kept, ndup, n_main, n_test, n_drop = set(), set(), 0, 0, 0, 0
        for i in range(len(smain)):
            o = smain[i]
            if o.NAME in read:
                ndup += 1
                continue
            read.add(o.NAME)
            if sel.get((mdid(o), exp), "main") == "main":
                out.add_sobj(o)
                kept.add(o.NAME)
                runs_used.setdefault(mdid(o), set()).add("main")
                n_main += 1
            else:
                n_drop += 1
        if stest is not None:
            tseen = set()
            for i in range(len(stest)):
                o = stest[i]
                if o.NAME in tseen:
                    continue
                tseen.add(o.NAME)
                if sel.get((mdid(o), exp)) == "test":
                    if o.NAME in kept:
                        raise SystemExit(f"exposure {exp}: object {o.NAME} would appear twice")
                    out.add_sobj(o)
                    kept.add(o.NAME)
                    runs_used.setdefault(mdid(o), set()).add("test")
                    n_test += 1

        # every slit selected from the test run should contribute at least one object
        for (mid, e), src in sel.items():
            if e == exp and src == "test" and not any(mdid(out[j]) == mid for j in range(len(out))):
                print(f"WARNING: exposure {exp}: no test-run object found for maskdef_id {mid}")

        fout = os.path.join(a.out, os.path.basename(fmain))
        out.write_to_fits(subheader(smain.header), fout, overwrite=True)
        out.write_info(fout.replace(".fits", ".txt"), "MultiSlit")
        print(f"{exp}: {n_main} objects from main, {n_test} from test, {n_drop} dropped, "
              f"{ndup} duplicates removed -> {fout}")

    mixed = sorted(m for m, r in runs_used.items() if m is not None and len(r) > 1)
    if mixed:
        print(f"NOTE: {len(mixed)} slit(s) combine exposures from both runs: {mixed}")
    else:
        print("Every slit takes all of its exposures from a single run.")


if __name__ == "__main__":
    main()
