#!/usr/bin/env python
"""
patch_dmost.py  DMOST_DIR  [--toler 0.5] [--wv-rms 0.4]

Apply the changes dmost needs to run in the PypeIt 2.0.1 environment and on crowded
M31-system masks, to a clone of https://github.com/marlageha/dmost:

  numpy >= 2.4 (required by PypeIt 2.0.1)
    dmost_create_maskfile.py, run_collate1d_marz.py   np.in1d   -> np.isin
    dmost_telluric.py                                  np.float  -> float
    dmost_emcee.py, dmost_coadd_emcee.py               numpy.core.multiarray.interp -> numpy.interp

  collate tolerance
    dmost_flexure.py re-collates the flexure-corrected spec1d files with a hard-coded
    1.0" tolerance and no wavelength-RMS cut. Set both to the values the `collate`
    stage of reduce_mask.sh uses (0.5", 0.4 pix), so the coadds used for templates,
    low-S/N velocities and equivalent widths group objects the same way as the
    collate_report.dat dmost builds its slit table from.

  telluric file names
    dmost_telluric.parse_tfile() reads H2O and O2 as fields 3 and 5 of the full template
    path split on "_", so any underscore in a folder name (dmost_raw, dmost_data, a
    user name) shifts the fields. Parse the file name only.

  pyspherematch (imported by dmost, not on PyPI; used only to match Marz results)
    writes DMOST_DIR/pyspherematch.py, a stand-in built on astropy. reduce_mask.sh
    puts DMOST_DIR on PYTHONPATH so it is found.

Each change is checked: the script stops without writing anything if the code it
expects is not there (dmost changed upstream). Re-running is safe; already-applied
changes are skipped. A stamp file DMOST_DIR/PATCHED_BY_REDUX records what was applied;
reduce_mask.sh refuses to run dmost without it.
"""
import argparse
import datetime
import re
import subprocess
import sys
from pathlib import Path

SHIM = '''"""
Stand-in for pyspherematch, written by patch_dmost.py. dmost uses it only in
dmost_create_maskfile.add_marz():  m1, m2, d = spherematch(ra1, dec1, ra2, dec2, tol)

Each (ra1, dec1) is matched to its nnearest-th nearest (ra2, dec2); pairs closer than
tol are returned as (index into 1, index into 2, separation). Angles in degrees.
"""
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord


def spherematch(ra1, dec1, ra2, dec2, tol=None, nnearest=1):
    ra1, dec1, ra2, dec2 = (np.atleast_1d(np.asarray(x, dtype=float)) for x in (ra1, dec1, ra2, dec2))
    if ra1.size == 0 or ra2.size < nnearest:
        empty = np.array([], dtype=int)
        return empty, empty, np.array([], dtype=float)
    c1 = SkyCoord(ra1 * u.deg, dec1 * u.deg)
    c2 = SkyCoord(ra2 * u.deg, dec2 * u.deg)
    idx2, sep, _ = c1.match_to_catalog_sky(c2, nthneighbor=nnearest)
    sep = sep.deg
    keep = np.ones(ra1.size, dtype=bool) if tol is None else sep <= tol
    return np.nonzero(keep)[0], np.asarray(idx2)[keep], sep[keep]
'''


def literal(path, old, new):
    """Replace every occurrence of `old`; skip if only `new` is present."""
    def apply(text):
        if old in text:
            return text.replace(old, new), f"{text.count(old)} x {old!r} -> {new!r}"
        if new in text:
            return text, "already applied"
        raise LookupError(f"neither {old!r} nor {new!r} found")
    return path, apply


def regex(path, pattern, repl, desc):
    """Replace the single match of `pattern` (re-applied each run, so values can change)."""
    def apply(text):
        hits = re.findall(pattern, text)
        if len(hits) != 1:
            raise LookupError(f"expected one match of {pattern!r}, found {len(hits)}")
        new_text = re.sub(pattern, repl, text)
        return new_text, ("already applied" if new_text == text else desc)
    return path, apply


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dmost_dir", help="clone of github.com/marlageha/dmost")
    ap.add_argument("--toler", default="0.5", help="collate tolerance in arcsec (default 0.5, as the collate stage)")
    ap.add_argument("--wv-rms", default="0.4", help="collate wavelength-RMS threshold in pix (default 0.4)")
    a = ap.parse_args()

    root = Path(a.dmost_dir).expanduser().resolve()
    core = root / "dmost" / "core"
    if not (root / "scripts" / "run_dmost.py").is_file() or not core.is_dir():
        sys.exit(f"patch_dmost.py: {root} is not a dmost clone (no scripts/run_dmost.py)")
    float(a.toler), float(a.wv_rms)  # fail early on non-numbers

    interp_old = "from numpy.core.multiarray import interp as compiled_interp"
    interp_new = "from numpy import interp as compiled_interp"
    edits = [
        literal(core / "dmost_create_maskfile.py", "np.in1d(", "np.isin("),
        literal(root / "scripts" / "run_collate1d_marz.py", "np.in1d(", "np.isin("),
        literal(core / "dmost_telluric.py", "np.float(", "float("),
        literal(core / "dmost_telluric.py", "spl = tfile.split('_')", "spl = os.path.basename(tfile).split('_')"),
        literal(core / "dmost_emcee.py", interp_old, interp_new),
        literal(core / "dmost_coadd_emcee.py", interp_old, interp_new),
        regex(core / "dmost_flexure.py",
              r"(pypeit_collate_1d --spec1d_files \.\./Science_flex/spec1d_\*fits) --toler [^'\"]*(['\"])",
              rf"\1 --toler {a.toler} --wv_rms_thresh {a.wv_rms}\2",
              f"flexure re-collate: --toler {a.toler} --wv_rms_thresh {a.wv_rms}"),
    ]

    # dry run first: nothing is written unless every edit applies
    results, errors = {}, []
    for path, apply in edits:
        text = results.get(path, (None, None))[0] or path.read_text()
        try:
            new_text, msg = apply(text)
            results[path] = (new_text, results.get(path, (None, []))[1] + [msg])
        except LookupError as e:
            errors.append(f"{path.relative_to(root)}: {e}")
    if errors:
        sys.exit("patch_dmost.py: dmost code differs from what this script expects; nothing written.\n  "
                 + "\n  ".join(errors))

    for path, (new_text, msgs) in results.items():
        if new_text != path.read_text():
            path.write_text(new_text)
        for m in msgs:
            print(f"{path.relative_to(root)}: {m}")

    shim = root / "pyspherematch.py"
    if shim.exists() and shim.read_text() != SHIM and "patch_dmost.py" not in shim.read_text():
        sys.exit(f"patch_dmost.py: {shim} exists and was not written by this script; remove it first")
    shim.write_text(SHIM)
    print("pyspherematch.py: stand-in written")

    try:
        commit = subprocess.run(["git", "-C", str(root), "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    (root / "PATCHED_BY_REDUX").write_text(
        f"patched {datetime.datetime.now():%Y-%m-%d %H:%M} by patch_dmost.py\n"
        f"dmost commit {commit}\ntoler {a.toler}\nwv_rms_thresh {a.wv_rms}\n")
    print(f"done: dmost {commit} patched (collate tolerance {a.toler}\", RMS cut {a.wv_rms} pix)")


if __name__ == "__main__":
    main()
