#!/usr/bin/env python
"""
wavecal_exclude.py  <configuration folder, e.g. $DEIMOS_RDX/2022B/M32RA1/keck_deimos_A>

Turn the wavelength-solution check into exclusions for merge_spec1d.py: for every
slit that wavecal_qa.py marked "bad" or "fail", write one `none` row per science
exposure of that slit's calibration group into <folder>/selection.csv.

PypeIt only flags a DEIMOS slit BADWVCALIB (and skips it) when no wavelength fit
exists at all. With the default full_template method a fit above the RMS
threshold is kept, so slits with wrong solutions are extracted like any other.

Rows written here carry a reason starting with "wavecal_qa:" and are regenerated
on every run, so re-running after a new calib stage replaces them. Rows added by
hand (section 7 of the README) are kept, except that a hand-written `main` or
`test` row for a slit with a bad wavelength solution is replaced by `none`: the
test reduction uses the same calibrations, so its wavelengths are wrong too.
"""
import argparse
import csv
import os
import re

AUTO = "wavecal_qa:"
EXCLUDE = ("bad", "fail")
FIELDS = ["maskdef_id", "exp", "source", "reason"]


def exp_tag(fn):
    """Exposure number as merge_spec1d.py reads it from file names."""
    m = re.search(r"DE\.\d+\.(\d+)", os.path.basename(fn))
    return m.group(1) if m else None


def science_by_group(pypeit_file):
    """{calib group (str): [exposure tags]} from the data block of a .pypeit file."""
    with open(pypeit_file) as f:
        lines = f.read().splitlines()
    try:
        i0 = next(i for i, l in enumerate(lines) if l.strip().lower() == "data read")
        i1 = next(i for i, l in enumerate(lines) if l.strip().lower() == "data end")
    except StopIteration:
        raise SystemExit(f"no data block in {pypeit_file}")
    rows = [l for l in lines[i0 + 1:i1] if "|" in l]
    if not rows:
        raise SystemExit(f"empty data block in {pypeit_file}")
    cols = [c.strip() for c in rows[0].strip().strip("|").split("|")]
    for need in ("filename", "frametype", "calib"):
        if need not in cols:
            raise SystemExit(f"column '{need}' missing from the data block of {pypeit_file}")
    groups, every = {}, []
    for l in rows[1:]:
        if l.lstrip().startswith("#"):
            continue
        r = dict(zip(cols, (c.strip() for c in l.strip().strip("|").split("|"))))
        if "science" not in r["frametype"].split(","):
            continue
        tag = exp_tag(r["filename"])
        if tag is None:
            raise SystemExit(f"cannot read an exposure number from {r['filename']}")
        every.append(tag)
        if r["calib"].lower() == "all":
            groups.setdefault("all", []).append(tag)
        else:
            for g in r["calib"].split(","):
                groups.setdefault(g.strip(), []).append(tag)
    for tag in groups.pop("all", []):
        for g in groups:
            groups[g].append(tag)
    return groups, every


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    p.add_argument("config_dir")
    a = p.parse_args()
    cfg = os.path.abspath(a.config_dir)
    qa_file = os.path.join(cfg, "wavecal_qa.csv")
    pf = os.path.join(cfg, os.path.basename(cfg) + ".pypeit")
    sel_file = os.path.join(cfg, "selection.csv")
    for f in (qa_file, pf):
        if not os.path.exists(f):
            raise SystemExit(f"missing {f}")

    groups, every = science_by_group(pf)

    # exclusions from the wavelength check
    auto = {}
    with open(qa_file) as f:
        for r in csv.DictReader(f):
            if r["status"] not in EXCLUDE:
                continue
            mid = int(r["maskdef_id"])
            if mid < 0:
                print(f"WARNING: {r['calib_key']} spat_id {r['spat_id']} has no maskdef_id; "
                      f"exclude it by hand")
                continue
            grp = r["calib_key"].split("_")[1]           # A_0_MSC01 -> calibration group 0
            exps = groups.get(grp)
            if exps is None:
                raise SystemExit(f"calibration group {grp} ({r['calib_key']}) has no science "
                                 f"frames in {pf}")
            why = f"{AUTO} {r['status']} {r['calib_key']} spat {r['spat_id']}"
            if r.get("reason"):
                why += f" ({r['reason']})"
            for e in exps:
                auto[(mid, e)] = why

    # hand-written rows, minus earlier automatic ones
    manual = {}
    if os.path.exists(sel_file):
        with open(sel_file) as f:
            for r in csv.DictReader(f):
                if (r.get("reason") or "").startswith(AUTO):
                    continue
                manual[(int(r["maskdef_id"]), r["exp"].strip())] = r

    out, replaced = [], []
    for key, r in manual.items():
        if key in auto and r["source"].strip().lower() != "none":
            replaced.append((key, r["source"]))
            continue
        out.append({k: r.get(k, "") for k in FIELDS})
    for (mid, e), why in auto.items():
        if (mid, e) in manual and manual[(mid, e)]["source"].strip().lower() == "none":
            continue                                     # already excluded by hand
        out.append(dict(maskdef_id=mid, exp=e, source="none", reason=why))
    out.sort(key=lambda r: (int(r["maskdef_id"]), r["exp"]))

    with open(sel_file, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(out)

    slits = sorted({m for m, _ in auto})
    print(f"{len(slits)} slit(s) excluded for their wavelength solution, "
          f"{len(auto)} slit-exposure(s) over {len(every)} science exposure(s): {slits}")
    for (mid, e), src in replaced:
        print(f"NOTE: maskdef_id {mid}, exposure {e}: hand-written '{src}' row replaced by 'none'")
    print(f"{len(out)} row(s) -> {sel_file}")


if __name__ == "__main__":
    main()
