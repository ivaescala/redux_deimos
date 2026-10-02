#!/usr/bin/env python
"""
make_slurm.py  SEMESTER/MASK [SEMESTER/MASK ...]  [--stages calib science qa collate] [--submit] [options]

Write (and optionally submit) one SLURM batch script per mask that runs the
requested stages of reduce_mask.sh in order, stopping at the first failure.

Everything about the user and environment is read from the cluster session that
runs this script. No conda environment needs to be active: the script needs only the
Python standard library, and the batch job itself activates the PypeIt environment.
  user, home, host          getpass / os / socket
  DEIMOS_RAW, DEIMOS_RDX    environment (DEIMOS_REDUX accepted for DEIMOS_RDX)
  conda installation        `conda info --base`, CONDA_EXE, or ~/miniforge3 etc.
  conda environment         --env (default: pypeit), checked to contain run_pypeit
  kit folder                location of reduce_mask.sh on PATH (or this script's folder)
  partition                 SLURM default partition from `sinfo` (or --partition)
  account                   first SLURM association from `sacctmgr` (or --account)
  e-mail                    `git config user.email` (or --mail; omitted if neither)

The batch script exports these explicitly, so the job does not depend on ~/.bashrc.
Scripts and SLURM logs go to $DEIMOS_RDX/<SEMESTER>/<MASK>/slurm/.
"""
import argparse
import datetime
import getpass
import os
import shutil
import socket
import subprocess
import sys
from pathlib import Path

STAGES = ["setup", "calib", "science", "qa", "collate"]


def sh(cmd):
    """Run a command and return stripped stdout, or '' if it is missing or fails."""
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=30, check=True).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return ""


def slurm_partition():
    parts = sh(["sinfo", "-h", "-o", "%P"]).split()
    for p in parts:
        if p.endswith("*"):
            return p[:-1]
    return parts[0] if parts else ""


def slurm_timelimit(partition):
    out = sh(["sinfo", "-h", "-p", partition, "-o", "%l"]).split()
    return out[0] if out else "unknown"


def slurm_account(user):
    for line in sh(["sacctmgr", "-nP", "show", "assoc", f"user={user}", "format=account"]).splitlines():
        if line.strip():
            return line.strip()
    return ""


def conda_base():
    """Locate the conda installation without needing any environment to be active."""
    out = sh(["conda", "info", "--base"])
    if out:
        return Path(out.splitlines()[-1]).resolve()
    exe = os.environ.get("CONDA_EXE") or shutil.which("conda")
    if exe:
        return Path(exe).resolve().parents[1]
    for name in ("miniforge3", "mambaforge", "miniconda3", "anaconda3"):
        cand = Path.home() / name
        if (cand / "etc" / "profile.d" / "conda.sh").exists():
            return cand
    sys.exit("make_slurm.py: no conda installation found (tried `conda info --base`, "
             "CONDA_EXE, ~/miniforge3, ~/mambaforge, ~/miniconda3, ~/anaconda3)")


def conda_setup(env_name):
    """Return (conda.sh path, env name, env prefix) for the environment the job will activate."""
    base = conda_base()
    conda_sh = base / "etc" / "profile.d" / "conda.sh"
    prefix = base / "envs" / env_name
    if not conda_sh.exists():
        sys.exit(f"make_slurm.py: {conda_sh} not found")
    if not (prefix / "bin" / "run_pypeit").exists():
        sys.exit(f"make_slurm.py: no run_pypeit in {prefix}; pass --env <name of your PypeIt environment>")
    return conda_sh, env_name, prefix


def main():
    p = argparse.ArgumentParser(description="Write SLURM jobs that run reduce_mask.sh stages.")
    p.add_argument("masks", nargs="+", help="mask folder(s) relative to $DEIMOS_RAW, e.g. 2022B/M32RA1")
    p.add_argument("--stages", nargs="+", default=["calib", "science", "qa", "collate"],
                   choices=STAGES, help="stages to run, in pipeline order (default: all but setup)")
    p.add_argument("--partition", default=None)
    p.add_argument("--account", default=None)
    p.add_argument("--mail", default=None, help="address for END/FAIL e-mails")
    p.add_argument("--env", default="pypeit",
                   help="conda environment the job activates (default: pypeit)")
    p.add_argument("--cpus", type=int, default=2)
    p.add_argument("--mem", default="32G")
    p.add_argument("--time", default="24:00:00", help="wall time, HH:MM:SS or D-HH:MM:SS")
    p.add_argument("--force", action="store_true", help="allow 'setup' to overwrite existing .pypeit files")
    p.add_argument("--submit", action="store_true", help="submit with sbatch after writing")
    a = p.parse_args()

    stages = [s for s in STAGES if s in a.stages]          # enforce pipeline order
    user, home, host = getpass.getuser(), Path.home(), socket.gethostname()

    raw = os.environ.get("DEIMOS_RAW", "")
    rdx = os.environ.get("DEIMOS_RDX") or os.environ.get("DEIMOS_REDUX") or ""
    if not raw or not rdx:
        sys.exit("make_slurm.py: DEIMOS_RAW and DEIMOS_RDX must be set (source ~/.bashrc)")
    raw, rdx = Path(raw).resolve(), Path(rdx).resolve()

    found = shutil.which("reduce_mask.sh")
    kit = Path(found).resolve().parent if found else Path(__file__).resolve().parent
    if not (kit / "reduce_mask.sh").exists():
        sys.exit(f"make_slurm.py: reduce_mask.sh not found on PATH or in {kit}")

    conda_sh, env, prefix = conda_setup(a.env)
    partition = a.partition if a.partition is not None else slurm_partition()
    account = a.account if a.account is not None else slurm_account(user)
    mail = a.mail if a.mail is not None else sh(["git", "config", "--get", "user.email"])

    print(f"user {user} on {host} | env {env} ({prefix}) | kit {kit}")
    print(f"raw {raw} | redux {rdx}")
    print(f"partition {partition or '(cluster default)'}"
          + (f" [time limit {slurm_timelimit(partition)}]" if partition else "")
          + f" | account {account or '(none)'} | mail {mail or '(none)'}")
    print(f"stages: {' '.join(stages)} | {a.cpus} cpus, {a.mem}, {a.time}\n")

    for mask in a.masks:
        out = rdx / mask
        if not (raw / mask).is_dir():
            print(f"SKIP {mask}: no raw folder {raw / mask}")
            continue
        has_pypeit = any(out.glob("keck_deimos_?/keck_deimos_?.pypeit"))
        if "setup" in stages and has_pypeit and not a.force:
            print(f"SKIP {mask}: .pypeit files exist and 'setup' would overwrite them (use --force)")
            continue
        if "setup" not in stages and not has_pypeit:
            print(f"SKIP {mask}: no .pypeit files yet; include 'setup' in --stages")
            continue

        name = mask.strip("/").replace("/", "_")                 # 2022B/M32RA1 -> 2022B_M32RA1
        jobdir = out / "slurm"
        jobdir.mkdir(parents=True, exist_ok=True)
        tag = "-".join(stages)
        script = jobdir / f"{name}_{tag}.slurm"

        sb = [f"#SBATCH --job-name=pypeit_{name}",
              "#SBATCH --nodes=1",
              "#SBATCH --ntasks=1",
              f"#SBATCH --cpus-per-task={a.cpus}",
              f"#SBATCH --mem={a.mem}",
              f"#SBATCH --time={a.time}",
              f"#SBATCH --output={jobdir}/%x-%j.out"]
        if partition:
            sb.append(f"#SBATCH --partition={partition}")
        if account:
            sb.append(f"#SBATCH --account={account}")
        if mail:
            sb += [f"#SBATCH --mail-user={mail}", "#SBATCH --mail-type=END,FAIL"]

        body = f"""#!/bin/bash
{chr(10).join(sb)}
#
# Generated by make_slurm.py on {datetime.datetime.now():%Y-%m-%d %H:%M} by {user}@{host}
# Mask {mask}: stages {' '.join(stages)}

set -euo pipefail

export DEIMOS_RAW={raw}
export DEIMOS_RDX={rdx}
export DEIMOS_REDUX={rdx}/
export PATH={kit}:$PATH
export OMP_NUM_THREADS=${{SLURM_CPUS_PER_TASK:-1}}
export MKL_NUM_THREADS=$OMP_NUM_THREADS
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS

source {conda_sh}
conda activate {env}

echo "job $SLURM_JOB_ID on $(hostname) | $(python -c 'import pypeit; print("PypeIt", pypeit.__version__)')"
for stage in {' '.join(stages)}; do
  echo "===== $stage started $(date '+%F %T')"
  reduce_mask.sh {mask} "$stage"
  echo "===== $stage finished $(date '+%F %T')"
done
"""
        script.write_text(body)
        script.chmod(0o755)
        print(f"wrote {script}")
        if a.submit:
            res = subprocess.run(["sbatch", str(script)], capture_output=True, text=True)
            print("  " + (res.stdout.strip() or res.stderr.strip()))


if __name__ == "__main__":
    main()
