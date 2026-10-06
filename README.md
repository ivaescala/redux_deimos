# Keck/DEIMOS slitmask reductions with PypeIt

Scripts for downloading archival Keck/DEIMOS multi-object spectroscopy from the
Keck Observatory Archive (KOA) and reducing it uniformly with
[PypeIt](https://pypeit.readthedocs.io), following the approach of
Geha et al. (2026) for Milky Way dwarf spheroidals, with additions for crowded,
high-surface-brightness fields in the M31 system (M32). The same workflow runs on
a laptop or on a SLURM cluster.

## Contents

| File | Purpose |
|---|---|
| `koa_deimos_fetch.py` | Query KOA, download science frames plus matching flats/arcs (optionally biases), sort them into one folder per slitmask |
| `reduce_mask.sh` | Run one PypeIt stage (`setup`, `calib`, `science`, `qa`, `collate`) for one mask |
| `deimos.par` | PypeIt parameter block appended to every `.pypeit` file |
| `assign_calib_by_night.py` | One calibration group per UT night; flags nights with incomplete calibrations |
| `wavecal_qa.py` | Rule-based check of every slit's wavelength solution |
| `skysub_qa.py` | Object-free sky-residual statistics and a Ca II triplet leakage test vs. distance from M32 |
| `make_slurm.py` | Write (and optionally submit) a SLURM job that runs `reduce_mask.sh` stages |
| `inspect_slits.py` | Spatial profiles and Ca II triplet spectra of chosen slits from existing spec2d files, including slits PypeIt rejected (section 7) |
| `check_recovered.py` | Per-exposure extraction and sky-line residual checks of target spectra in spec1d files (section 7) |
| `merge_spec1d.py` | Build one spec1d per exposure from two reductions according to a selection table, removing duplicate objects (section 7) |

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
│       └── M32RA1/                    same name as the raw mask folder
│           ├── setup_files/           keck_deimos_*.sorted (configuration census)
│           ├── slurm/                 batch scripts and job logs
│           ├── keck_deimos_A/         one folder per instrument configuration
│           │   ├── keck_deimos_A.pypeit
│           │   ├── keck_deimos_A.log
│           │   ├── Calibrations/      Edges_, Slits_, Arc_, Tilts_, WaveCalib_, Flat_ …
│           │   ├── Science/           spec2d_*.fits, spec1d_*.fits
│           │   ├── QA/
│           │   ├── collate1d/         coadded 1D spectra
│           │   ├── wavecal_qa.csv
│           │   ├── skysub_qa.csv
│           │   ├── slit_report.txt    slit flags per exposure (qa stage)
│           │   ├── inspect_slits/     PNGs from inspect_slits.py (section 7)
│           │   ├── selection.csv      per slit-exposure source: main, test or none (section 7)
│           │   └── Science_merged/    merged spec1d files, input to collate (section 7)
│           └── test_maskfrac/         re-reduction of rejected slits (section 7)
│               ├── keck_deimos_A.pypeit
│               ├── Calibrations -> ../keck_deimos_A/Calibrations
│               └── Science/
└── pypeit/                            this repository (on $PATH)
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
export DEIMOS_REDUX="$DEIMOS_RDX/"        # same folder; dmost expects this name, with trailing slash
export PATH="$DEIMOS_ROOT/pypeit:$PATH"
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
| `collate` | Coadds 1D spectra across exposures (0.5″ matching, slits with wavelength RMS > 0.4 pix excluded) | `collate1d/`. If slits were recovered (section 7), collate `Science_merged/` instead |

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
```

### 5c. Cluster: interactive sessions

Use a compute node, not the login node, for anything heavier than `setup` or
inspecting files:

```bash
srun --cpus-per-task=2 --mem=16G --time=02:00:00 --pty bash -l
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
4. `python make_slurm.py <semester>/<mask> --stages collate --submit` if no slits need recovering; otherwise follow section 7 and collate the merged files.

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
flagged it. The RMS cut alone misses solutions locked onto the wrong lines with a
small RMS. Alignment boxes appear as `box`. Keep bad slits in the reduction and
exclude them downstream by `maskdef_id` (`wavecal_qa.csv`), for example as
`none` rows in the selection table of section 7. PypeIt does not flag most of
these slits itself, and collate's RMS cut misses a wrong solution with a small RMS.

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
`BADSKYSUB` in `slit_report.txt` that have a wavelength solution, and note their
`maskdef_id`s.

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

**5. Write the selection table** `$CFG/selection.csv`, one row per slit-exposure
that does not come from the main reduction:

```
maskdef_id,exp,source,reason
<maskdef_id>,<exp>,test,recovered with max_mask_frac=0.97
<maskdef_id>,<exp>,none,wrong wavelength solution (wavecal_qa)
```

`exp` is the exposure number in the file names (the field after the date in
`DE.<date>.<exp>.<n>`); `source` is `main`, `test` or `none`; unlisted
slit-exposures default to `main`. The rules:

- **One reduction per slit.** All exposures of a slit come from the same
  reduction.
- **Excluded slit-exposures:** wrong wavelength solutions, failed extractions,
  `sky_k` outside the reference range (main and test spectra alike), traces on a
  neighbour, and blends.

**6. Merge and collate.**

```bash
python $KIT/merge_spec1d.py --main $CFG/Science --test $TEST/Science \
    --select $CFG/selection.csv --out $CFG/Science_merged

pypeit_collate_1d --spec1d_files $CFG/Science_merged/spec1d_*.fits \
    --tolerance 0.5 --wv_rms_thresh 0.4 --outdir $CFG/collate1d
```

`merge_spec1d.py` prints, per exposure, the objects taken from each reduction,
the objects dropped and the duplicates removed.

### 7d. Caveat

A subset run does not reproduce the main reduction exactly. With `slitspatnum`,
PypeIt fits no slitmask offset (`MaskOFF` 0.00 in `slit_report.txt`). Targets are
still matched to their design entries, but serendipitous objects get coordinates
shifted by the offset the main reduction applied. Object finding can also
differ: a slit-exposure that extracted normally in the main run can fail in the
subset run. Take every slit-exposure the main reduction already extracted from
the main reduction.

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

## References

- PypeIt: Prochaska et al. 2020, JOSS, 5, 2308; Prochaska et al. 2020, Zenodo
- Geha et al. 2026, ApJ, 999, 140 — uniform PypeIt/dmost reduction of Milky Way
  dwarf spheroidals ([dmost](https://github.com/marlageha/dmost))
- Data from the Keck Observatory Archive; follow KOA's acknowledgment guidelines
