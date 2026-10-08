# Keck/DEIMOS slitmask reductions with PypeIt

Scripts for downloading archival Keck/DEIMOS multi-object spectroscopy from the
Keck Observatory Archive (KOA), reducing it uniformly with
[PypeIt](https://pypeit.readthedocs.io), and measuring velocities with
[dmost](https://github.com/marlageha/dmost), following the approach of
Geha et al. (2026) for Milky Way dwarf spheroidals, with additions for crowded,
high-surface-brightness fields in the M31 system (M32). The same workflow runs on
a laptop or on a SLURM cluster.

## Contents

| File | Purpose |
|---|---|
| `koa_deimos_fetch.py` | Query KOA, download science frames plus matching flats/arcs (optionally biases), sort them into one folder per slitmask |
| `reduce_mask.sh` | Run one stage (`setup`, `calib`, `science`, `qa`, `collate`, `marz`, `dmost`) for one mask |
| `deimos.par` | PypeIt parameter block appended to every `.pypeit` file |
| `assign_calib_by_night.py` | One calibration group per UT night; flags nights with incomplete calibrations |
| `wavecal_qa.py` | Rule-based check of every slit's wavelength solution |
| `wavecal_exclude.py` | Writes `none` rows to the selection table for slits that fail `wavecal_qa.py`, one per science exposure (run by the `collate` stage) |
| `skysub_qa.py` | Object-free sky-residual statistics and a Ca II triplet leakage test vs. distance from M32 |
| `make_slurm.py` | Write (and optionally submit) a SLURM job that runs `reduce_mask.sh` stages |
| `inspect_slits.py` | Spatial profiles and Ca II triplet spectra of chosen slits from existing spec2d files, including slits PypeIt rejected (section 7) |
| `check_recovered.py` | Per-exposure extraction and sky-line residual checks of target spectra in spec1d files (section 7) |
| `merge_spec1d.py` | Build one spec1d per exposure from the main and test reductions according to the selection table, dropping excluded slits and duplicate objects (run by the `collate` stage) |
| `patch_dmost.py` | One-time fixes to a dmost clone for PypeIt 2.0.1 / numpy ≥ 2.4 and this reduction's collate settings (section 9) |
| `dmost_prep.py` | Build dmost's working folder from the `collate` outputs and check its inputs (run by the `marz` and `dmost` stages) |

---

## 1. Installing PypeIt

All scripts were developed against **PypeIt 2.0.1**. Use the same version on every
machine.

### 1a. Local machine (macOS or Linux)

Install [Miniforge](https://github.com/conda-forge/miniforge) if you don't have a
conda installation, then create a dedicated environment:

```bash
conda create -n pypeit python=3.12
conda activate pypeit
pip install "pypeit==2.0.1" pykoa
```

PypeIt 2.0.1 installs its viewers (ginga, Qt) as core dependencies, so the
interactive tools (`pypeit_chk_edges`, `pypeit_show_2dspec`, …) work out of the box.

### 1b. Cluster

Install Miniforge in your home directory (no administrator rights needed):

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3
$HOME/miniforge3/bin/conda init bash      # then log out and back in
```

Create the environment. On clusters with an older system C library, the newest PyQt6
releases have no compatible binary wheel and pip tries to compile Qt from source,
which fails with a `qmake` error. Pin the last compatible release and forbid source
builds:

```bash
conda create -n pypeit python=3.11
conda activate pypeit
pip install "pypeit==2.0.1" "pyqt6==6.9.1" pykoa --only-binary=pyqt6,pyqt6-qt6,pyqt6-sip
```

Repeat the same pin and flag on any later `pip install --upgrade pypeit`.

**Compute nodes without internet access.** PypeIt downloads wavelength templates
on first use, and astropy looks up observatory locations online. Populate both
caches once, from a login node:

```bash
conda activate pypeit
pypeit_cache_github_data keck_deimos
python -c "from astropy.coordinates import EarthLocation; EarthLocation.of_site('Keck Observatory')"
```

The second command also serves dmost, which looks up Keck's location for the
heliocentric correction (section 9).

### 1c. Verify

```bash
run_pypeit -h
python -c "import pypeit; print(pypeit.__version__)"     # 2.0.1
python -c "import qtpy; print(qtpy.API_NAME)"            # PyQt6
```

---

## 2. Directory structure

Raw data, reduced data and this repository sit side by side under one root
(`$DEIMOS_ROOT`, e.g. `~/deimos`). Raw and reduced trees mirror each other as
`<semester>/<mask>`:

```
$DEIMOS_ROOT/
├── raw/                               $DEIMOS_RAW
│   └── 2022B/                         koa_deimos_fetch.py --outdir
│       ├── M32RA1/                    DE.*.fits: science, flats, arcs (+ biases)
│       ├── M32RA2/
│       ├── meta/                      KOA query tables, download manifest
│       └── bias/                      only biases with no matching mask night
├── redux/                             $DEIMOS_RDX
│   └── 2022B/
│       ├── M32RA1/                    same name as the raw mask folder
│       │   ├── setup_files/           keck_deimos_*.sorted (configuration census)
│       │   ├── slurm/                 batch scripts and job logs
│       │   ├── keck_deimos_A/         one folder per instrument configuration
│       │   │   ├── keck_deimos_A.pypeit
│       │   │   ├── keck_deimos_A.log
│       │   │   ├── Calibrations/      Edges_, Slits_, Arc_, Tilts_, WaveCalib_, Flat_ …
│       │   │   ├── Science/           spec2d_*.fits, spec1d_*.fits
│       │   │   ├── QA/
│       │   │   ├── collate1d/         coadded 1D spectra and collate_report.dat
│       │   │   ├── wavecal_qa.csv
│       │   │   ├── skysub_qa.csv
│       │   │   ├── slit_report.txt    slit flags per exposure (qa stage)
│       │   │   ├── inspect_slits/     PNGs from inspect_slits.py (section 7)
│       │   │   ├── selection.csv      test or none per exposure of slits inspected in section 7; wavecal exclusions added by collate
│       │   │   └── Science_merged/    merged spec1d files, input to collate
│       │   └── test_maskfrac/         re-reduction of rejected slits (section 7)
│       │       ├── keck_deimos_A.pypeit
│       │       ├── Calibrations -> ../keck_deimos_A/Calibrations
│       │       └── Science/
│       └── dmost/                     dmost working folders, one per mask (section 9)
│           ├── marz_files/
│           └── M32RA1/
├── pypeit/                            this repository (on $PATH)
├── dmost/                             patched dmost clone ($DMOST_DIR, section 9)
└── dmost_data/                        dmost templates and sky lines ($DMOST_DATA, section 9)
```

Set it up:

```bash
mkdir -p ~/deimos/raw ~/deimos/redux
git clone https://github.com/ivaescala/redux_deimos.git ~/deimos/pypeit
```

Raw data are never modified: the `.pypeit` files point to the raw frames by
absolute path, so a reduction can be deleted and redone without touching the
downloads.

---

## 3. Environment variables

Add this block to `~/.bashrc` (Linux, clusters) or `~/.zshrc` (macOS). Edit the
first line to match your root folder.

```bash
# --- DEIMOS / PypeIt reductions ---
export DEIMOS_ROOT="$HOME/deimos"
export DEIMOS_RAW="$DEIMOS_ROOT/raw"
export DEIMOS_RDX="$DEIMOS_ROOT/redux"
export PATH="$DEIMOS_ROOT/pypeit:$PATH"
# optional, only if dmost or its data are not at the default locations:
# export DMOST_DIR="$DEIMOS_ROOT/dmost"
# export DMOST_DATA="$DEIMOS_ROOT/dmost_data"
# --- end DEIMOS ---
```

The semester is given in the paths passed to the scripts, e.g.
`--outdir $DEIMOS_RAW/2022B` or `reduce_mask.sh 2022B/M32RA1 setup`.

Load and check:

```bash
source ~/.bashrc
echo "$DEIMOS_RAW  $DEIMOS_RDX"
which reduce_mask.sh                       # should be in $DEIMOS_ROOT/pypeit
```

Cluster notes:

- **Login shells** often read `~/.bash_profile` and not `~/.bashrc`. If the variables
  are empty after logging in, add `[ -f ~/.bashrc ] && . ~/.bashrc` to `~/.bash_profile`.
- **SLURM jobs** written by `make_slurm.py` export these variables explicitly, so
  they do not depend on `~/.bashrc` being read.

---

## 4. Downloading raw data from KOA

`koa_deimos_fetch.py` runs a cone search on the telescope pointing (default: M32,
0.10° radius) for science frames in one semester, then selects every flat and arc
taken through the same slitmask and grating within ±1 UT day of a science night.
Files are downloaded with PyKOA and sorted into `<outdir>/<slitmask>/`.


**1. Inspect what exists** (no download):

```bash
conda activate pypeit
python $DEIMOS_ROOT/pypeit/koa_deimos_fetch.py --outdir $DEIMOS_RAW/2022B \
    --semester 2022B --list-only
```

This prints one line per mask and grating, with target name, nights, exposure
time and PI. For other targets use `--ra`, `--dec` and `--radius`.

**2. Download and sort:**

```bash
python $DEIMOS_ROOT/pypeit/koa_deimos_fetch.py --outdir $DEIMOS_RAW/2022B \
    --semester 2022B --masks <mask1> <mask2>
```

Omit `--masks` to download every mask the query returns. Add `--biases` for bias
frames (PypeIt does not use them for DEIMOS by default). Files land in `lev0/`
first and are sorted into mask folders in one pass at the end, followed by a
per-folder table of science, flat, arc and bias counts.

**3. Re-runs are safe.** Files already on disk are skipped and truncated downloads
are deleted and re-fetched, so if the run is interrupted or reports missing files,
run the same command again with the same flags. To sort files already in `lev0/`
without querying KOA: `--organize-only`.

**Proprietary data:** log in once with `--login --cookie ~/.koa_cookie.txt`, then
pass `--cookie ~/.koa_cookie.txt` on later runs.

---

## 5. Running PypeIt

### 5a. Stages

`reduce_mask.sh <semester>/<mask> <stage>` runs one stage for one mask, reading
raw frames from `$DEIMOS_RAW/<semester>/<mask>` and writing to
`$DEIMOS_RDX/<semester>/<mask>`. It builds every path from the environment
variables, so it can be run from any folder. Run the stages in order and inspect
the output between them.

| Stage | What it does | Check before moving on |
|---|---|---|
| `setup` | Configuration census (`.sorted`), writes `keck_deimos_<X>.pypeit` with `deimos.par` appended, assigns one calibration group per night, prints the data block | Grating, number of configurations, frame types, exposure times; ≥3 flats and ≥1 arc per night |
| `calib` | Calibrations only; writes `wavecalib_report.txt` and `wavecal_qa.csv` | Slit edges; slits failing the wavelength checks (section 6) |
| `science` | Full reduction, reusing the calibrations | `Science/spec2d_*.fits` and `spec1d_*.fits`: one of each per exposure, each holding every slit on all four mosaics |
| `qa` | `skysub_qa.csv` and `slit_report.txt` (slit flags per exposure) | Sky residuals vs. radius (section 6); slits flagged `BADSKYSUB` (section 7) |
| `collate` | Adds `none` rows to `selection.csv` for slits failing `wavecal_qa.py`, merges the main and any test reduction into `Science_merged/`, and coadds 1D spectra across exposures (0.5″ matching, slits with wavelength RMS > 0.4 pix excluded) | `collate1d/`; the excluded and merged object counts printed per exposure |
| `marz` | Optional. Writes a Marz input file of the coadds, for flagging galaxies and QSOs (section 9) | Save the Marz results as `marz_<mask>_MG.mz` |
| `dmost` | Flexure, telluric, template and velocity fits, and Ca II triplet equivalent widths with dmost (section 9) | `not_in_dmost.csv`; plots in `QA/` |

Re-running `setup` overwrites the `.pypeit` file. After an interruption, re-run
the same stage: finished calibrations are reused, so only the step in progress is
repeated.

### 5b. Local machine

```bash
conda activate pypeit
bash reduce_mask.sh 2022B/M32RA1 setup
bash reduce_mask.sh 2022B/M32RA1 calib
bash reduce_mask.sh 2022B/M32RA1 science
bash reduce_mask.sh 2022B/M32RA1 qa
bash reduce_mask.sh 2022B/M32RA1 collate
bash reduce_mask.sh 2022B/M32RA1 marz      # optional
bash reduce_mask.sh 2022B/M32RA1 dmost
```

### 5c. Cluster: interactive sessions

Use a compute node, not the login node, for anything heavier than `setup` or
inspecting files:

```bash
srun --cpus-per-task=1 --mem=16G --time=02:00:00 --pty bash -l
conda activate pypeit
bash reduce_mask.sh 2022B/M32RA1 setup
```

### 5d. Cluster: batch jobs

`make_slurm.py` writes one batch script per mask that activates the PypeIt
environment, exports the paths, and runs the requested stages in order, stopping
at the first failure. It needs no environment active and reads everything from
the session: user, `DEIMOS_RAW`/`DEIMOS_RDX`, the conda installation, the default
partition (`sinfo`), your account (`sacctmgr`) and e-mail (`git config`).

```bash
python make_slurm.py 2022B/M32RA1                                  # write only; inspect the script
python make_slurm.py 2022B/M32RA1 --stages calib --submit          # write and submit
python make_slurm.py 2022B/M32RA1 2022B/M32RA2 --stages science qa collate --submit
python make_slurm.py 2022B/M32RA1 --stages dmost --submit
```

Defaults: environment `pypeit`, stages `calib science qa collate`, 1 CPU, 32 GB,
24 h. Override with `--env`, `--stages`, `--cpus`, `--mem`, `--time`,
`--partition`,  `--account`, and `--mail`. `setup` only
runs when requested, and is refused for masks that already have a `.pypeit`
file unless you add `--force`. Scripts and logs go to
`$DEIMOS_RDX/<semester>/<mask>/slurm/`.

Recommended sequence per mask, so the calibrations are checked before the
science stage:

1. `bash reduce_mask.sh <semester>/<mask> setup` in an interactive session (minutes), and check the data block.
2. `python make_slurm.py <semester>/<mask> --stages calib --submit`, then review `wavecal_qa.py`'s output in the job log.
3. `python make_slurm.py <semester>/<mask> --stages science qa --submit`, then check `slit_report.txt` for slits rejected for sky subtraction (section 7).
4. If slits need recovering, follow section 7 up to writing the selection table. Then `python make_slurm.py <semester>/<mask> --stages collate --submit`. Slits failing `wavecal_qa.py` are excluded in either case.
5. Optionally `bash reduce_mask.sh <semester>/<mask> marz` and flag extragalactic objects in Marz, then `python make_slurm.py <semester>/<mask> --stages dmost --submit` (section 9).

---

## 6. Quality checks

**Slit edges.** In `pypeit_chk_edges`, each slit's traced edges should follow the
illuminated slit along its full length, with no slit split in two or two slits
merged, and every slit labelled with a mask-design ID. Text-only alternative that
works over SSH: `pypeit_parse_slits <config folder>/Calibrations/Slits_A_0_MSC01.fits.gz`.

The viewers need a display. On a cluster either log in with `ssh -Y` (XQuartz
on macOS) and use `--mpl`, or copy the calibrations to your laptop and view them
there with the same PypeIt version:

```bash
rsync -av <user>@<cluster>:<redux path>/<semester>/<mask>/keck_deimos_A/Calibrations ~/qa/<mask>/
```

**Wavelength solutions.** `wavecal_qa.py` marks a slit bad if its lines don't run
blue to red, its dispersion differs by more than 5% from the detector median
(≈0.32 Å/pix for 1200G), it uses fewer than 20 arc lines, the lines span less
than 60% (or more than 100%) of the spectrum, its RMS exceeds 0.4 pix, or PypeIt
flagged it. Slits with no solution at all are `fail`; alignment boxes appear as
`box`.

PypeIt itself excludes few of these slits. With DEIMOS's default method, it flags
a slit `BADWVCALIB`, and skips it in the science stage, only when no fit exists. A fit above its RMS threshold (0.15 × arc FWHM) is logged as poor
but kept, so the slit is extracted with that solution. Collate's RMS cut then
drops those above 0.4 pix but misses solutions locked onto the wrong lines with a
small RMS. Bad slits therefore stay in the reduction and are excluded at the
`collate` stage: `wavecal_exclude.py` writes a `none` row to `selection.csv` for
every `bad` or `fail` slit and every science exposure of its calibration group,
and `merge_spec1d.py` drops those slits' objects (the target and any serendipitous
detections) before coadding. Only `Science_merged/` and `collate1d/` are
filtered; `Science/` still holds every extracted slit.

**Sky subtraction.** `skysub_qa.py` computes the residual
(science − sky − object model)/noise on object-free pixels of every slit, binned
by distance from M32. Noise-limited subtraction gives a median near 0 and a width
near 1 (Geha et al. 2026 find 0.94). `d_line` compares the residual in the three
Ca II triplet lines at M32's velocity with the neighbouring continuum: negative
values mean unresolved galaxy light was left in the spectra, positive values mean
over-subtraction. Use `--vsys` to test M31's local velocity instead. The test
skips slits with too few pixels outside the object masks, which on crowded
slits can be most of them; section 7 covers those.

**Slit flags.** `slit_report.txt` lists each slit's PypeIt flags for every
exposure. Alignment boxes carry `BOXSLIT`; slits without a wavelength solution
carry `BADWVCALIB` and related flags; slits rejected for sky subtraction carry
`BADSKYSUB` and have no extracted spectrum in that exposure (section 7).

---

## 7. Crowded slits: recovering slits rejected for sky subtraction

### 7a. Why slits are rejected

Before fitting the sky on a slit, PypeIt masks a region of ±1 FWHM around every
detected object, plus a few trimmed pixels at each slit edge. If more than 80% of
the slit is masked (`max_mask_frac`), PypeIt fits no sky, flags the slit
`BADSKYSUB`, and extracts nothing from it in that exposure, not even the forced
extraction of an undetected target. The log records each case as
`Giving up on global sky-subtraction`, with the masked percentage.

Crowded fields push slits over this threshold in three ways. Short slits start
with a large fraction masked, because the trimmed edge pixels count toward the
total. Every additional source on the slit, such as a serendipitous neighbour,
adds its own masked region. And a neighbour blended with the target makes object
finding fit a single, much wider profile, so the masked region around the target
can cover most of the slit. Because these detections and fits change slightly
between exposures, the same slit can be rejected in one exposure and extracted
in another.

PypeIt also leaves the wavelength image and tilts at zero on rejected slits and
flags their pixels as off-slit. `skysub_qa.py` therefore reports no usable pixels
for them, while `inspect_slits.py` rebuilds their wavelengths from the
calibration files.

### 7b. Additional quality checks

| Check | Tool | What it tells you |
|---|---|---|
| Which slits were rejected, in which exposures | `slit_report.txt` | `BADSKYSUB` per slit and exposure |
| What is on a rejected slit | `inspect_slits.py` | Whether a compact source sits at the design position, how far and how bright the nearest other source is, what share of the slit's light the target carries, and whether the Ca II triplet is visible |
| Whether recovered slits are sky-subtracted normally | `skysub_qa.py` on the test reduction | Residual statistics against the rest of the mask (only for slits with enough sky pixels) |
| Whether each extraction is reliable | `check_recovered.py` | Profile width, optimal vs. boxcar flux, and the fraction of sky-line flux left in the spectrum (`sky_k`), compared with the main reduction |

`inspect_slits.py` works on the existing spec2d files and needs no re-reduction.
For each slit it prints and plots:

- `t_off`: offset of the target peak from the design position (pix)
- `t_fwhm`: width of that peak (pix)
- `s_fwhm`: FWHM of a single star on the same mosaic, from the spec1d (pix)
- `n_off`, `n/t`: offset and relative height of the brightest other source
- `t_frac`: fraction of the slit's light within ±1 FWHM of the target

`check_recovered.py` prints, for each target and exposure, the optimal-extraction
FWHM, the median optimal and boxcar counts in 8400–8800 Å and their ratio, and two
measures of sky-line residuals in the target spectrum:

- `sky_r`: correlation between the continuum-subtracted target and sky spectra,
  showing whether residuals follow the sky-line pattern.
- `sky_k`: fraction of the sky-line flux left in the target spectrum (negative =
  over-subtracted, positive = under-subtracted, near 0 = well subtracted).

`check_recovered.py`, run without `--maskdef` on the main reduction, prints the
distribution of every quantity over all targets. Those percentiles are the
reference for judging recovered spectra. A recovered slit-exposure is accepted
when all of the following hold:

- the target is matched to its design entry, at a position consistent with the
  other exposures;
- the optimal-extraction FWHM is within the main reduction's range. Wider
  profiles absorb sky and neighbour light;
- the optimal flux does not exceed the boxcar flux by a large factor, and the
  boxcar flux is consistent between exposures. The boxcar is 3″ wide by default
  and already includes neighbours within 1.5″, so an optimal flux several times
  higher means the extraction failed;
- `sky_k` lies within the main reduction's 1st–99th percentile range.

### 7c. Workflow

Set the paths once:

```bash
CFG=$DEIMOS_RDX/2022B/M32RA1/keck_deimos_A
TEST=$DEIMOS_RDX/2022B/M32RA1/test_maskfrac
KIT=$DEIMOS_ROOT/pypeit
```

**1. Identify rejected slits.** After the `qa` stage, list slits flagged
`BADSKYSUB` in `slit_report.txt` in at least one exposure, leave out those that
fail `wavecal_qa.py` (they are excluded at `collate` regardless), and note the
remaining `maskdef_id`s.

**2. Inspect them in the existing reduction.**

```bash
python $KIT/inspect_slits.py $CFG/Science/spec2d_*.fits \
    --maskdef <maskdef_id> ... --outdir $CFG/inspect_slits
```

A slit is worth recovering when a compact source sits within ~2 pix of the
design position. Treat it as a blend, and leave it out, when the target peak is
two or more times the width of a single star. When a neighbour of comparable
brightness lies within ~1″, keep the slit but check in step 4 that the
extraction sits on the target and that the neighbour is extracted separately.

**3. Re-reduce the candidates in a separate folder.** Copy the `.pypeit` file
and link the existing calibrations:

```bash
mkdir -p $TEST
cp $CFG/keck_deimos_A.pypeit $TEST/
ln -s $CFG/Calibrations $TEST/Calibrations
```

In `$TEST/keck_deimos_A.pypeit`, add to the parameter block, under the existing
`[rdx]` and `[reduce]` headers if present:

```
[rdx]
    slitspatnum = <det>:<spat_id>, <det>:<spat_id>, ...     # e.g. MSC01:<spat_id>, one per candidate
[reduce]
    [[skysub]]
        max_mask_frac = 0.97
```

Then run:

```bash
run_pypeit $TEST/keck_deimos_A.pypeit -r $TEST
```

**4. Check the test reduction.**

```bash
for f in $TEST/Science/spec2d_*.fits; do echo "### $(basename "$f")"; pypeit_parse_slits "$f"; done > $TEST/slit_report.txt
python $KIT/skysub_qa.py $TEST/Science/spec2d_*.fits --out $TEST/skysub_qa.csv
cat $TEST/Science/spec1d_*.txt
python $KIT/check_recovered.py $TEST/Science/spec1d_*.fits > $TEST/check_test.txt
python $KIT/check_recovered.py $CFG/Science/spec1d_*.fits  > $CFG/check_main.txt
```

Apply the acceptance criteria of section 7b to each recovered slit-exposure.

**5. Write the selection table** `$CFG/selection.csv`. It decides, for the
slits inspected in step 2 only, which of their exposures reach `collate`:

```
maskdef_id,exp,source,reason
<maskdef_id>,<exp>,test,recovered with max_mask_frac=0.97
<maskdef_id>,<exp>,none,test sky_k -0.57
<maskdef_id>,<exp>,none,not recovered by the test run (BADSKYSUB at max_mask_frac=0.97)
```

`exp` is the exposure number in the file names (the field after the date in
`DE.<date>.<exp>.<n>`); `source` is `main`, `test` or `none`, and the reason
column is free text that `merge_spec1d.py` ignores. The rules:

- **Only slits inspected in step 2 get rows.** Every other slit keeps its
  main-run spectra unchanged and needs no row: unlisted slit-exposures default to
  `main`.
- **Slits judged in step 2 not worth re-reducing are `none` in every exposure**,
  including exposures the main run extracted. They take no spectra from either
  run.
- **Every exposure of a test-reduction slit gets a row**, with source `test` or
  `none`, never `main`, even for exposures the main run did extract. All of a
  slit's spectra then come from one processing.
- **`test`** when the recovered slit-exposure passes every acceptance criterion
  of section 7b.
- **`none`** when it fails any of them (failed extraction, `sky_k` outside the
  reference range, trace on a neighbour, inconsistent flux between exposures),
  when the test run did not extract it, or when the test extraction is a blend
  (FWHM two or more times a single star's).
- **No minimum number of exposures.** A slit with one or two accepted exposures
  is kept.

Slits with wrong wavelength solutions need no rows here: the `collate` stage adds
them (reason starting with `wavecal_qa:`).

**6. Merge and collate.**

```bash
bash $KIT/reduce_mask.sh 2022B/M32RA1 collate
```

This runs `wavecal_exclude.py`, then `merge_spec1d.py` with `$TEST/Science` as the
test reduction, then `pypeit_collate_1d` on `Science_merged/`. `merge_spec1d.py`
prints, per exposure, the objects taken from each reduction, the objects dropped
and the duplicates removed.

### 7d. Caveat

A subset run does not reproduce the main reduction exactly. With `slitspatnum`,
PypeIt fits no slitmask offset (`MaskOFF` 0.00 in `slit_report.txt`). Targets are
still matched to their design entries, but serendipitous objects get coordinates
shifted by the offset the main reduction applied. Object finding can also
differ: a slit-exposure that extracted normally in the main run can fail in the
subset run. Run the relaxed `max_mask_frac` only on the rejected slits (via
`slitspatnum`), never on the whole mask, so that slits extracted normally in the
main run are not put at risk.

---

## 8. Reduction choices

The parameter block reproduces the two overrides Geha et al. (2026) apply; all
other settings are PypeIt's `keck_deimos` defaults (four blue/red mosaics,
template-based wavelength calibration, overscan-only bias subtraction, slitmask
matching and forced extraction of undetected targets):

```
[calibrations]
    [[wavelengths]]
        refframe = observed      # heliocentric correction applied after velocity fitting
[flexure]
    spec_method = skip           # flexure measured downstream from sky lines
```

Their data-quality cuts carry over: 1200G grating, science exposures > 60 s,
≥3 flats and ≥1 arc per mask, wavelength RMS ≤ 0.4 pix.

Crowded masks add one optional override, applied only in a separate
re-reduction of rejected slits (section 7):

```
[reduce]
    [[skysub]]
        max_mask_frac = 0.97     # default 0.80
```

---

## 9. Velocities with dmost

[dmost](https://github.com/marlageha/dmost) (Geha et al. 2026) measures
line-of-sight velocities and Ca II triplet equivalent widths from the PypeIt
spectra. It replaces PypeIt's flexure and heliocentric corrections (switched off in
section 8) with its own: a linear flexure fit to sky lines per slit, a synthetic
telluric model per exposure, a PHOENIX stellar template per star, and an MCMC fit
of velocity and telluric shift in each exposure. Velocities are combined across
exposures as a weighted mean. When a single-exposure fit fails (typically at
S/N ≈ 2 per pixel), dmost fits the coadded spectrum instead.

The `marz` and `dmost` stages replace steps 1–4 of dmost's own README
(`run_dmost_planfiles`, `run_dmost_mask_setup`, its `run_pypeit` call and its
collate): those read the authors' mask list and write their own `.pypeit` file.
dmost starts from this reduction's `collate` outputs instead.

### 9a. Installing dmost

In the PypeIt environment:

```bash
conda activate pypeit
git clone https://github.com/marlageha/dmost.git $DEIMOS_ROOT/dmost
pip install emcee corner h5py astroplan
python $DEIMOS_ROOT/pypeit/patch_dmost.py $DEIMOS_ROOT/dmost
```

dmost is not pip-installed; the stages put the clone on `PYTHONPATH`.
`patch_dmost.py` makes the following change, including modifications for
compatibility with PypeIt 2.0.1:

| Change | Why |
|---|---|
| Flexure re-collate: `--toler 1.` → `--toler 0.5 --wv_rms_thresh 0.4` | dmost re-coadds the flexure-corrected spectra for its template, low-S/N and equivalent-width fits. At 1″ it can group a target with a neighbour that the 0.5″ `collate` stage keeps separate |

### 9b. Templates and sky lines

dmost reads three sets of files that are not in its repository. Put them under
`$DEIMOS_ROOT/dmost_data/` (or set `DMOST_DATA`):

```
dmost_data/
├── templates/
│   ├── pheonix/grid1/dmost*.fits      PHOENIX stellar templates, S/N ≥ 25
│   ├── pheonix/grid2/dmost*.fits      10 < S/N < 25
│   ├── pheonix/grid3/dmost*.fits      S/N ≤ 10
│   ├── tellurics/telluric_0.02A*.fits          coarse TelFit grid
│   └── fine_tellurics/telluric_0.02A_h2o_*.fits  fine TelFit grid
└── Other_data/sky_single_mg.dat       sky emission lines for the flexure fit
```

The notebooks in `dmost/notebooks/` show how the template grids were built
(PHOENIX spectra from Husser et al. 2013; telluric spectra with TelFit, which
needs LBLRTM). The sky-line list has no source in the repository. Ask the dmost
authors for these folders. The `dmost` stage lists whatever is missing and stops.

### 9c. Working folder and outputs

Run the `marz` and `dmost` stages as described in sections 5a, 5b and 5d. Both
first run `dmost_prep.py`, which builds dmost's working folder
`$DEIMOS_RDX/<semester>/dmost/<mask>/` (`<mask>_A`, `<mask>_B`, … for masks with
several instrument configurations) from symlinks:

| In the dmost folder | Points to |
|---|---|
| `Science/` | `keck_deimos_A/Science_merged/`: dmost sees only the merged spectra, so the wavelength exclusions and recovered slits carry over |
| `collate1d/`, `collate_report.dat` | `keck_deimos_A/collate1d/` and the report inside it |
| `dmost_raw/rawdata_<year>/` | the raw science frames, from which dmost reads slit widths |
| `dmost_raw/templates`, `dmost_raw/Other_data` | `$DMOST_DATA` |

**Marz (optional).** The `marz` stage writes `marz_files/marz_<mask>.fits` from the
coadds. Load it in [Marz](https://samreay.github.io/Marz/), confirm or flag
galaxies and QSOs, and save the results as `marz_files/marz_<mask>_MG.mz` (dmost
looks for this exact name). Without it every object is treated as a star.

Assign the Marz quality flag (QOP) as follows. dmost reads it as `marz_flag` and
treats QOP > 2 as extragalactic, 2 as bad and anything lower as a star:

| QOP | Use for | What dmost does |
|---|---|---|
| 3 or 4 | galaxies and QSOs (3 possible, 4 secure) | excluded from the stellar fits; velocity is cz from the Marz redshift |
| 2 | bad data | excluded from everything |
| 0 or 1 | stars (0 = left unflagged in Marz, 1 = flagged; dmost sets 0 to 1) | fitted as stars |

Keep Marz's automatic QOP assignment switched off (the default): it
gives confident stellar matches QOP 6, which dmost treats as extragalactic.

**Outputs**, in the dmost folder:

- `dmost/dmost_<mask>.fits`: the `mask` (per exposure) and `slits` (per object)
  tables. Final velocities are `dmost_v`, `dmost_v_err`; `coadd_flag` = 1 marks
  velocities from the coadded spectrum; equivalent widths are `cat`, `naI`, `mgI`.
- `<mask>_dmost.log`, plots in `QA/`, MCMC chains in `emcee/`.
- `not_in_dmost.csv`: extracted spectra dmost does not measure (below).

**Re-running.** dmost resumes from its own table and skips finished steps, so an
interrupted run continues where it stopped. After re-running `collate`, the old
outputs no longer match: `dmost_prep.py` detects this and stops. Run with
`DMOST_FRESH=1 bash reduce_mask.sh 2022B/M32RA1 dmost` to move the old outputs to
`archive_<time>/` and start over.

### 9d. Objects dmost does not measure

dmost builds its object list from `collate_report.dat`. `pypeit_collate_1d` writes
no coadd, and no report line, for an object with a single spectrum, so a slit
extracted in only one exposure (for example, rejected for sky subtraction in the
others) never reaches dmost, although dmost could fit that one exposure. These
spectra are listed in `not_in_dmost.csv` (spec1d file, PypeIt name, `maskdef_id`,
object name), and their number is printed when the stage starts. Slits extracted
in two or more exposures are unaffected.

---

## References

- PypeIt: Prochaska et al. 2020, JOSS, 5, 2308; Prochaska et al. 2020, Zenodo
- Geha et al. 2026, ApJ, 999, 140 — uniform PypeIt/dmost reduction of Milky Way
  dwarf spheroidals ([dmost](https://github.com/marlageha/dmost))
- Marz: Hinton et al. 2016, Astronomy and Computing, 15, 61
- PHOENIX stellar library: Husser et al. 2013, A&A, 553, A6
- Data from the Keck Observatory Archive; follow KOA's acknowledgment guidelines
