#!/usr/bin/env python
"""
make_dmost_data.py  COMMAND  [--data DIR] [options]

Build the files dmost reads from DMOST_DATA (default: $DEIMOS_ROOT/dmost_data),
following dmost's notebooks prepare_pheonix_templates.ipynb and mk_telluric_grid.ipynb,
for the grids dmost's code requests.

  skylines          Other_data/sky_single_mg.dat, copied from PypeIt          [pypeit env]
  phoenix-download  PHOENIX HiRes spectra + wavelength file -> phoenix_raw/   [login node]
  phoenix-prepare   templates/pheonix/grid{1,2,3}/dmost_lte_*.fits            [pypeit env + ppxf]
  install-telfit    TelFit with LBLRTM, LNFL and the AER line list            [telfit env, login node]
  telluric          templates/tellurics/ and fine_tellurics/ with TelFit      [telfit env]
                      --test           one model (H2O 50%, O2 1.00), timed
                      --task I --ntasks N   every N-th model, starting at I
  telluric-slurm    write (and --submit) a SLURM array job running `telluric`
  check             what is complete and what is missing

Stellar templates (notebook cell 3: grid1 active, grid2 and grid3 commented out):
  grid1  S/N >= 25       [Fe/H] 0,-0.5,-1,-2,-3,-4  logg 1,3,5  Teff 2500-7000 (500 K), 8000
  grid2  10 < S/N < 25   [Fe/H] 0,-1,-2,-3,-4       logg 1,3,5  Teff 3000,3500,4000,5000,6000,7000,8000
  grid3  S/N <= 10       [Fe/H] 0,-1,-3             logg 2,5    Teff 3000,3500,4000,5000,6000,7000
  Each spectrum is converted to air, trimmed to 6000-9600 A and log-rebinned with
  ppxf.log_rebin. The notebook's call log_rebin([lmin, lmax], flux, 0.6) passes 0.6 as
  `oversample` in ppxf <= 8.0 (`velscale` only became the third argument in ppxf 9):
  0.6 x as many log pixels as input pixels. That is reproduced here with the explicit
  velscale it implies, using the actual PHOENIX wavelength array.

Telluric grids (notebook cells 2-5, TelFit at Mauna Kea, 600-960 nm, smoothed and
resampled to 0.02 A in air):
  tellurics/        H2O 5-100% (step 5) x O2 0.70-2.15 (step 0.05)              600 models
                    = notebook cell 4. The notebook names them new_telluric_*, which
                    dmost's glob (telluric_0.02A*) does not match; named telluric_* here.
  fine_tellurics/   the files dmost_telluric.final_telluric_values() copies:
                    H2O 2-100% (step 2) x O2 0.60-2.02 (step 0.02), and
                    H2O 5-100% (step 5) x O2 2.04-2.20 (step 0.02)              3780 models
  Models that appear in both grids are computed once. Not covered: an exposure whose
  fit gives H2O < 2.5% together with O2 > 2.02, which dmost rounds to h2o_0.
"""
import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path

C_KMS = 299792.458

# ---------------------------------------------------------------- stellar templates
PHOENIX_URL = "https://phoenix.astro.physik.uni-goettingen.de/data/HiResFITS/"
PHOENIX_LIB = "PHOENIX-ACES-AGSS-COND-2011"
WAVE_FILE = f"WAVE_{PHOENIX_LIB}.fits"
GRIDS = {
    "grid1": dict(z=["-0.0", "-0.5", "-1.0", "-2.0", "-3.0", "-4.0"], logg=["1.00", "3.00", "5.00"],
                  teff=[2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 8000]),
    "grid2": dict(z=["-0.0", "-1.0", "-2.0", "-3.0", "-4.0"], logg=["1.00", "3.00", "5.00"],
                  teff=[3000, 3500, 4000, 5000, 6000, 7000, 8000]),
    "grid3": dict(z=["-0.0", "-1.0", "-3.0"], logg=["2.00", "5.00"],
                  teff=[3000, 3500, 4000, 5000, 6000, 7000]),
}
OVERSAMPLE = 0.6          # third positional argument of log_rebin in the notebook
WMIN, WMAX = 6000., 9600.  # air, Angstrom

# ---------------------------------------------------------------- telluric templates
TELL_STEP = 0.02           # Angstrom
TELFIT_REPO = "https://github.com/freddavies/Telluric-Fitter-py3.git"
TELFIT_COMMIT = "aab37f5"


def data_root(arg):
    if arg:
        return Path(arg).expanduser().resolve()
    if os.environ.get("DMOST_DATA"):
        return Path(os.environ["DMOST_DATA"]).expanduser().resolve()
    if os.environ.get("DEIMOS_ROOT"):
        return Path(os.environ["DEIMOS_ROOT"]).expanduser().resolve() / "dmost_data"
    sys.exit("make_dmost_data.py: set DMOST_DATA or DEIMOS_ROOT, or pass --data")


def phoenix_name(teff, logg, z):
    return f"lte0{teff}-{logg}{z}.{PHOENIX_LIB}-HiRes.fits"


def phoenix_needed():
    """{(z, raw file name): [grids]} for every spectrum any grid uses."""
    need = {}
    for g, p in GRIDS.items():
        for z in p["z"]:
            for lg in p["logg"]:
                for t in p["teff"]:
                    need.setdefault((z, phoenix_name(t, lg, z)), []).append(g)
    return need


def telluric_models():
    """{(h2o, o2): [relative output paths]} for the coarse and fine grids."""
    def name(h, o):
        return f"telluric_0.02A_h2o_{int(h)}_o2_{o:2.2f}_.fits"
    models = {}
    for h in range(5, 105, 5):                                   # notebook cell 4
        for k in range(30):                                      # np.arange(0.7, 2.2, 0.05)
            o = round(0.70 + 0.05 * k, 2)
            models.setdefault((h, o), []).append(f"templates/tellurics/{name(h, o)}")
    for h in range(2, 102, 2):                                   # dmost final_telluric_values
        for k in range(30, 102):                                 # 0.60 ... 2.02
            o = round(0.02 * k, 2)
            models.setdefault((h, o), []).append(f"templates/fine_tellurics/{name(h, o)}")
    for h in range(5, 105, 5):
        for k in range(102, 111):                                # 2.04 ... 2.20
            o = round(0.02 * k, 2)
            models.setdefault((h, o), []).append(f"templates/fine_tellurics/{name(h, o)}")
    return dict(sorted(models.items()))


# ================================================================ skylines
def cmd_skylines(a):
    from pypeit import dataPaths
    src = Path(dataPaths.sky_spec.get_file_path("sky_single_mg.dat"))
    out = data_root(a.data) / "Other_data" / "sky_single_mg.dat"
    out.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, out)
    from astropy.io import ascii
    t = ascii.read(out)
    print(f"{out}: {len(t)} sky lines, {min(t['Wave']):.1f}-{max(t['Wave']):.1f} A (from {src})")


# ================================================================ PHOENIX
def fetch(url, out, tries=3):
    """Download url to out (via out.part). Returns 'ok', 'skip' or 'missing'."""
    if out.exists():
        return "skip"
    part = out.with_name(out.name + ".part")
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r, open(part, "wb") as fh:
                size = int(r.headers.get("Content-Length") or -1)
                shutil.copyfileobj(r, fh, 1 << 20)
            if size >= 0 and part.stat().st_size != size:
                raise OSError(f"incomplete: {part.stat().st_size} of {size} bytes")
            part.rename(out)
            return "ok"
        except urllib.error.HTTPError as e:
            if e.code == 404:
                part.unlink(missing_ok=True)
                return "missing"
            err = e
        except OSError as e:
            err = e
        time.sleep(5 * (i + 1))
    part.unlink(missing_ok=True)
    raise SystemExit(f"make_dmost_data.py: {url} failed after {tries} tries ({err}); re-run to resume")


def cmd_phoenix_download(a):
    raw = data_root(a.data) / "phoenix_raw"
    raw.mkdir(parents=True, exist_ok=True)
    print(f"{WAVE_FILE}: {fetch(PHOENIX_URL + WAVE_FILE, raw / WAVE_FILE)}")
    need = phoenix_needed()
    missing = []
    for i, (z, fn) in enumerate(need, 1):
        d = raw / f"Z{z}"
        d.mkdir(exist_ok=True)
        st = fetch(f"{PHOENIX_URL}{PHOENIX_LIB}/Z{z}/{fn}", d / fn)
        if st == "missing":
            missing.append(f"Z{z}/{fn}  ({', '.join(need[(z, fn)])})")
        if st != "skip" or i == len(need):
            print(f"[{i}/{len(need)}] Z{z}/{fn}: {st}")
    (raw / "missing.txt").write_text("".join(m + "\n" for m in missing))
    n_have = len(need) - len(missing)
    print(f"{n_have} of {len(need)} spectra in {raw}; {len(missing)} not in the PHOENIX library"
          + (f" (listed in {raw / 'missing.txt'})" if missing else ""))


def cmd_phoenix_prepare(a):
    import numpy as np
    from astropy.io import fits
    from astropy.table import Table
    from ppxf import ppxf_util

    root = data_root(a.data)
    raw = root / "phoenix_raw"
    wfile = raw / WAVE_FILE
    if not wfile.exists():
        sys.exit(f"make_dmost_data.py: {wfile} missing; run phoenix-download first")

    # notebook cell 7: vacuum -> air, trim to the DEIMOS range
    vwave = fits.getdata(wfile).astype(float)
    s2 = (1.e4 / vwave) ** 2
    wave = vwave / (1. + 0.0000834254 + (0.02406147 / (130. - s2)) + (0.00015998 / (38.9 - s2)))
    m = (wave >= WMIN) & (wave <= WMAX)
    wave_deimos = wave[m]
    n_in = wave_deimos.size

    # the notebook's log_rebin call assumes linear sampling between the two end points
    lin = np.linspace(wave_deimos[0], wave_deimos[-1], n_in)
    dev = np.max(np.abs(wave_deimos - lin))
    dl = np.diff(wave_deimos)
    print(f"PHOENIX grid in {WMIN:.0f}-{WMAX:.0f} A (air): {n_in} pixels, step {dl.min():.5f}-{dl.max():.5f} A, "
          f"max departure from linear sampling {dev:.4f} A")

    # velscale implied by log_rebin(lam_range, spec, oversample=0.6) in ppxf <= 8.0
    step = (wave_deimos[-1] - wave_deimos[0]) / (n_in - 1)
    edges = np.array([wave_deimos[0] - step / 2, wave_deimos[-1] + step / 2])
    velscale = C_KMS * np.diff(np.log(edges))[0] / int(n_in * OVERSAMPLE)
    print(f"log rebinning: oversample {OVERSAMPLE} -> velscale {velscale:.4f} km/s "
          f"({velscale / C_KMS * 8500:.4f} A per pixel at 8500 A)")

    done = skipped = 0
    for g in GRIDS:
        gdir = root / "templates" / "pheonix" / g
        gdir.mkdir(parents=True, exist_ok=True)
        for (z, fn), grids in phoenix_needed().items():
            if g not in grids:
                continue
            f = raw / f"Z{z}" / fn
            if not f.exists():
                continue
            with fits.open(f) as hdu:
                data = hdu[0].data.astype(float)
                hdr = hdu[0].header
            logg, teff, feh = hdr["PHXLOGG"], hdr["PHXTEFF"], hdr["PHXM_H"]
            iso = "{:0.6f}".format(0.02 * 10 ** feh)                      # notebook convert_zfeh
            outfile = gdir / "dmost_lte_{:0.0f}_{}_{}_.fits".format(teff, logg, feh)
            if outfile.exists() and not a.clobber:
                skipped += 1
                continue
            spec, ln_lam, vs = ppxf_util.log_rebin(wave_deimos, data[m], velscale=velscale)
            t = Table([[ln_lam], [spec], [logg], [teff], [feh], [iso], [vs]],
                      names=("wave", "flux", "logg", "teff", "feh", "Z", "vscale"))
            tmp = outfile.with_name(outfile.name + ".tmp")
            t.write(tmp, format="fits", overwrite=True)
            os.replace(tmp, outfile)
            done += 1
            if len(f"{g}/{outfile.name}") > 36:
                print(f"WARNING: {g}/{outfile.name} is longer than dmost's 36-character chi2_tfile column")
        print(f"{g}: {len(list(gdir.glob('dmost*.fits')))} templates in {gdir}")
    print(f"{done} written, {skipped} already present")


# ================================================================ TelFit
def patch_telfit_setup(text):
    """TELLURICMODELING from the environment; -std=legacy for GNU Fortran targets."""
    old = "TELLURICMODELING = '{}/.TelFit/'.format(os.environ['HOME'])"
    new = "TELLURICMODELING = os.environ.get('TELLURICMODELING', '{}/.TelFit/'.format(os.environ['HOME']))"
    if old in text:
        text = text.replace(old, new)
    elif new not in text:
        raise SystemExit("make_dmost_data.py: TelFit setup.py differs from the tested version (TELLURICMODELING)")
    anchor = "    #Build the executables\n"
    legacy = ("    # added by make_dmost_data.py: gfortran >= 8 needs -std=legacy for LNFL/LBLRTM\n"
              "    import re as _re\n"
              "    for _mk in ['lnfl/build/makefile.common', 'lblrtm/build/makefile.common']:\n"
              "        _p = TELLURICMODELING + _mk\n"
              "        _out, _target = [], ''\n"
              "        for _line in open(_p).read().splitlines(True):\n"
              "            _m = _re.match(r'^(\\w+):', _line)\n"
              "            _target = _m.group(1) if _m else _target\n"
              "            if 'GNU' in _target and 'FCFLAG=\"' in _line and '-std=legacy' not in _line:\n"
              "                _line = _line.replace('FCFLAG=\"', 'FCFLAG=\"-std=legacy ', 1)\n"
              "            _out.append(_line)\n"
              "        open(_p, 'w').write(''.join(_out))\n\n")
    if "added by make_dmost_data.py" not in text:
        if anchor not in text:
            raise SystemExit("make_dmost_data.py: TelFit setup.py differs from the tested version (build step)")
        text = text.replace(anchor, legacy + anchor, 1)
    return text


def telfit_root():
    return Path(os.environ.get("TELLURICMODELING", Path.home() / ".TelFit")).expanduser().resolve()


def cmd_install_telfit(a):
    src = Path(a.src).expanduser().resolve()
    if not (src / "setup.py").exists():
        subprocess.check_call(["git", "clone", TELFIT_REPO, str(src)])
        subprocess.check_call(["git", "-C", str(src), "checkout", "-q", TELFIT_COMMIT])
    setup = src / "setup.py"
    setup.write_text(patch_telfit_setup(setup.read_text()))
    root = telfit_root()
    env = dict(os.environ, TELLURICMODELING=str(root) + "/")
    print(f"building LBLRTM and LNFL in {root} (downloads ~470 MB from Zenodo; LNFL then builds TAPE3)")
    subprocess.check_call([sys.executable, "setup.py", "install"], cwd=src, env=env)
    if not (root / "rundir1").is_dir():
        sys.exit(f"make_dmost_data.py: TelFit installed, but {root}/rundir1 was not created")
    out = subprocess.run([sys.executable, "-c", "from telfit import Modeler"], cwd="/", env=env,
                         capture_output=True, text=True)
    if out.returncode:
        sys.exit("make_dmost_data.py: `from telfit import Modeler` failed:\n" + out.stderr[-2000:])
    print(f"TelFit installed; run directories in {root}. Next: make_dmost_data.py telluric --test")


class PrivateRundir:
    """Copy TelFit's rundir1 into a private folder, so concurrent jobs never share one.

    TelFit hands out run directories via lock files in the shared TELLURICMODELING root;
    on a cluster those locks are not reliable and stale ones block forever. TAPE3 and the
    lblrtm executable stay symlinks to the shared installation."""

    def __enter__(self):
        src = telfit_root() / "rundir1"
        if not src.is_dir():
            sys.exit(f"make_dmost_data.py: no TelFit run directory {src}; run install-telfit "
                     "(and set TELLURICMODELING if it is not ~/.TelFit)")
        self.root = Path(tempfile.mkdtemp(prefix="telfit_", dir=os.environ.get("TMPDIR")))
        shutil.copytree(src, self.root / "rundir1", symlinks=True,
                        ignore=shutil.ignore_patterns("OutputModels", "*.lock", "TAPE1[0-9]*"))
        (self.root / "rundir1" / "OutputModels").mkdir()
        return str(self.root) + "/"

    def __exit__(self, *exc):
        shutil.rmtree(self.root, ignore_errors=True)


def run_telfit(root, h2o, o2, wave_deimos):
    """mk_telluric_grid.ipynb run_telfit(), without the plot."""
    import numpy as np
    from astropy import convolution
    from telfit import Modeler

    alt, lat = 4.2, 19.8                       # hardwired to Mauna Kea
    wavestart, waveend = 600.0, 960.0          # nm
    modeler = Modeler(TelluricModelingDirRoot=root)
    model = modeler.MakeModel(humidity=h2o, o2=o2 * 1e5, lowfreq=1e7 / waveend,
                              highfreq=1e7 / wavestart, lat=lat, alt=alt)
    tell = model.toarray()
    wave = 10. * tell[:, 0]                    # nm -> A (air: vac2air=True by default)
    flux = tell[:, 1]
    bins = wave - np.roll(wave, 1)
    sig_res = TELL_STEP / np.median(bins)
    smooth_flux = convolution.convolve(flux, convolution.Gaussian1DKernel(sig_res))
    return np.interp(wave_deimos, wave, smooth_flux)


def write_telluric(root_data, rels, h2o, o2, wave_deimos, tflux):
    from astropy.table import Table
    t = Table([wave_deimos, tflux], names=("wave", "flux"))
    t.meta["h2o"] = h2o
    t.meta["o2"] = o2
    for rel in rels:
        out = root_data / rel
        out.parent.mkdir(parents=True, exist_ok=True)
        tmp = out.with_name(out.name + f".tmp{os.getpid()}")
        t.write(tmp, format="fits", overwrite=True)
        os.replace(tmp, out)


def cmd_telluric(a):
    import numpy as np
    root_data = data_root(a.data)
    wave_deimos = np.arange(6000, 9600, TELL_STEP)
    models = telluric_models()
    keys = list(models)
    if a.test:
        keys = [(50, 1.00)]
    else:
        if a.ntasks < 1 or not 0 <= a.task < a.ntasks:
            sys.exit("make_dmost_data.py: need 0 <= --task < --ntasks")
        keys = keys[a.task::a.ntasks]
    todo = [k for k in keys if a.test or not all((root_data / r).exists() for r in models[k])]
    print(f"{len(models)} models in total; this task: {len(keys)}, {len(todo)} still to compute")

    with PrivateRundir() as root:
        for i, (h2o, o2) in enumerate(todo, 1):
            t0 = time.time()
            tflux = run_telfit(root, float(h2o), o2, wave_deimos)
            write_telluric(root_data, models[(h2o, o2)], float(h2o), o2, wave_deimos, tflux)
            dt = time.time() - t0
            print(f"[{i}/{len(todo)}] H2O {h2o:3d}  O2 {o2:.2f}  {dt:6.1f} s  -> "
                  + ", ".join(models[(h2o, o2)]), flush=True)
            if a.test:
                print(f"transmission {tflux.min():.3f}-{tflux.max():.3f}; "
                      f"all {len(models)} models need ~{dt * len(models) / 3600:.0f} CPU-hours")


def cmd_telluric_slurm(a):
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from make_slurm import conda_base, slurm_account, slurm_partition, sh
    import getpass

    root_data = data_root(a.data)
    base = conda_base()
    conda_sh = base / "etc" / "profile.d" / "conda.sh"
    if not (base / "envs" / a.env / "bin" / "python").exists():
        sys.exit(f"make_dmost_data.py: no conda environment {a.env} in {base}; pass --env")
    partition = a.partition if a.partition is not None else slurm_partition()
    account = a.account if a.account is not None else slurm_account(getpass.getuser())
    mail = a.mail if a.mail is not None else sh(["git", "config", "--get", "user.email"])
    jobdir = root_data / "slurm"
    jobdir.mkdir(parents=True, exist_ok=True)
    sb = [f"#SBATCH --job-name=dmost_telluric", f"#SBATCH --array=0-{a.ntasks - 1}",
          "#SBATCH --nodes=1", "#SBATCH --ntasks=1", "#SBATCH --cpus-per-task=1",
          f"#SBATCH --mem={a.mem}", f"#SBATCH --time={a.time}", f"#SBATCH --output={jobdir}/%x-%A_%a.out"]
    if partition:
        sb.append(f"#SBATCH --partition={partition}")
    if account:
        sb.append(f"#SBATCH --account={account}")
    if mail:
        sb += [f"#SBATCH --mail-user={mail}", "#SBATCH --mail-type=END,FAIL,ARRAY_TASKS"]
    script = jobdir / "telluric_array.slurm"
    script.write_text(f"""#!/bin/bash
{chr(10).join(sb)}
#
# dmost telluric grid: {len(telluric_models())} TelFit models in {a.ntasks} tasks.
# Re-submitting is safe: finished models are skipped.
set -euo pipefail
export TELLURICMODELING={telfit_root()}/
source {conda_sh}
conda activate {a.env}
python {Path(__file__).resolve()} telluric --data {root_data} \\
    --task $SLURM_ARRAY_TASK_ID --ntasks $SLURM_ARRAY_TASK_COUNT
""")
    script.chmod(0o755)
    print(f"wrote {script}")
    if a.submit:
        r = subprocess.run(["sbatch", str(script)], capture_output=True, text=True)
        print("  " + (r.stdout.strip() or r.stderr.strip()))


# ================================================================ check
def cmd_check(a):
    root = data_root(a.data)
    ok = True
    sky = root / "Other_data" / "sky_single_mg.dat"
    print(f"sky lines        {'present' if sky.exists() else 'MISSING (run skylines)'}")
    ok &= sky.exists()

    raw = root / "phoenix_raw"
    absent = {l.split()[0].split("/")[1] for l in (raw / "missing.txt").read_text().splitlines()
              } if (raw / "missing.txt").exists() else set()
    for g, p in GRIDS.items():
        n_list = len(p["z"]) * len(p["logg"]) * len(p["teff"])
        n_absent = sum(1 for (z, fn), gs in phoenix_needed().items() if g in gs and fn in absent)
        files = glob.glob(str(root / "templates" / "pheonix" / g / "dmost*.fits"))
        need = n_list - n_absent
        state = "complete" if files and len(files) >= need else "INCOMPLETE"
        print(f"stellar {g}    {len(files):4d} of {need:4d} available in PHOENIX "
              f"({n_list} listed, {n_absent} not in the library)  {state}")
        ok &= state == "complete"

    models = telluric_models()
    for sub in ("tellurics", "fine_tellurics"):
        rels = [r for v in models.values() for r in v if f"/{sub}/" in r]
        miss = [r for r in rels if not (root / r).exists()]
        print(f"telluric {sub:15s} {len(rels) - len(miss):5d} of {len(rels):5d}"
              + (f"  missing e.g. {Path(miss[0]).name}" if miss else "  complete"))
        ok &= not miss
    print("all dmost data present" if ok else "incomplete")
    sys.exit(0 if ok else 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--data", help="dmost data root (default $DMOST_DATA or $DEIMOS_ROOT/dmost_data)")
    sub.add_parser("skylines", parents=[common])
    sub.add_parser("phoenix-download", parents=[common])
    p = sub.add_parser("phoenix-prepare", parents=[common])
    p.add_argument("--clobber", action="store_true", help="rewrite existing templates")
    p = sub.add_parser("install-telfit", parents=[common])
    p.add_argument("--src", default=os.path.expandvars("$DEIMOS_ROOT/telfit_src"),
                   help="where to clone TelFit (default $DEIMOS_ROOT/telfit_src)")
    p = sub.add_parser("telluric", parents=[common])
    p.add_argument("--test", action="store_true", help="compute one model and report the time")
    p.add_argument("--task", type=int, default=0)
    p.add_argument("--ntasks", type=int, default=1)
    p = sub.add_parser("telluric-slurm", parents=[common])
    p.add_argument("--ntasks", type=int, default=100, help="array size (default 100)")
    p.add_argument("--env", default="telfit", help="conda environment with TelFit (default telfit)")
    p.add_argument("--time", default="24:00:00")
    p.add_argument("--mem", default="8G")
    p.add_argument("--partition", default=None)
    p.add_argument("--account", default=None)
    p.add_argument("--mail", default=None)
    p.add_argument("--submit", action="store_true")
    sub.add_parser("check", parents=[common])
    a = ap.parse_args()
    {"skylines": cmd_skylines, "phoenix-download": cmd_phoenix_download,
     "phoenix-prepare": cmd_phoenix_prepare, "install-telfit": cmd_install_telfit,
     "telluric": cmd_telluric, "telluric-slurm": cmd_telluric_slurm, "check": cmd_check}[a.cmd](a)


if __name__ == "__main__":
    main()
