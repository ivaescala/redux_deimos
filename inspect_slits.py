#!/usr/bin/env python
"""
inspect_slits.py  SPEC2D [SPEC2D ...]  --maskdef ID [ID ...]  [--outdir DIR] [--calib DIR]

Look at what is on a set of slits using only the existing spec2d (and spec1d) files,
without re-running PypeIt. Intended for slits that PypeIt rejected (BADSKYSUB), which
have no sky model and no extraction, but whose processed science image is still stored.

For each maskdef_id it writes one PNG (one row per exposure):
  1. spatial profile: median of the science image along wavelength (7000-9000 A),
     which suppresses sky lines. Marked: the target's design position, and a model
     single star (Gaussian with the median FWHM of the stars PypeIt detected on the same
     mosaic in that exposure) centred on the target peak and scaled to its height.
     Light the model does not account for comes from other sources on the slit.
  2. 2D cutout around the Ca II triplet, sky-subtracted along lines of constant
     wavelength (see below)
  3. 1D spectrum in the Ca II triplet window: sum within +-3 px of the target peak,
     5-px smoothed; the last row averages the exposures. Dotted lines mark the triplet
     at rest (vacuum wavelengths, as PypeIt calibrates); the shaded bands cover
     -600 to 0 km/s, where M31 and M32 stars fall.

Printed per slit/exposure: slit length, design position, target offset from it, target
FWHM, the model-star FWHM, the brightest other source (offset and height relative to the
target), and the fraction of the slit's light within +-1 FWHM of the target.

Sky estimate: pixels farther than 2 FWHM from the target, fit with a cubic spline in
wavelength (knots every 0.3 A; iterative clipping against the residual scatter measured
as a function of wavelength, rejecting positive outliers harder), then evaluated at every pixel's own wavelength and
never extrapolated beyond the wavelengths it was fit over. Using each pixel's wavelength follows the
tilt of the sky lines across the slit. Light from neighbouring stars that survives the
clipping enters the sky and lowers the target's continuum. It is still a crude sky: use
the spectra to recognise Ca II triplet absorption, not to measure it.

PypeIt writes no wavelengths for slits it rejected (the waveimg is zero there), so for
those slits the wavelengths are rebuilt from the Tilts and WaveCalib calibration files
(default: the Calibrations folder next to the Science folder).
"""
import argparse
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CAT = (8500.36, 8544.44, 8664.52)        # Ca II triplet, vacuum wavelengths
CAT_WIN = (8450.0, 8700.0)
PROF_WIN = (7000.0, 9000.0)
VRANGE = (-600.0, 0.0)                   # km/s, shaded around each triplet line
C_KMS = 299792.458
# Pixel-level flags that mark genuinely bad data. OFFSLITS and EXTRACT are left out on
# purpose: PypeIt drops rejected slits from its slit image, so every pixel on a
# BADSKYSUB slit carries OFFSLITS even though the data there are fine.
PIXFLAGS = ["BPM", "CR", "SATURATION", "MINCOUNTS", "IS_NAN", "IVAR0", "IVAR_NAN", "BADSCALE"]
BAND = 3                                 # half-width (px) of the 1D sum
SKY_KNOT = 0.3                           # Angstrom, knot spacing of the sky spline
FWHM_DEFAULT = 6.0                       # px, if no spec1d stars are available


def rectify(img, ok, left, right):
    """Shift each spectral row so column k is k pixels from the left slit edge."""
    nspec = img.shape[0]
    nk = int(np.ceil(np.nanmedian(right - left))) + 1
    out = np.full((nspec, nk), np.nan)
    for j in range(nspec):
        x0 = int(np.round(left[j]))
        x1 = min(x0 + nk, img.shape[1])
        if x0 < 0 or x1 <= x0:
            continue
        out[j, :x1 - x0] = np.where(ok[j, x0:x1], img[j, x0:x1], np.nan)
        out[j, int(np.ceil(right[j] - x0)):] = np.nan
    return out


def target_stats(prof, expected, star_fwhm):
    """Target peak near the design position, its width, and the brightest other source."""
    out = dict(it=np.nan, off=np.nan, fwhm=np.nan, nbr_off=np.nan, nbr_ratio=np.nan,
               tfrac=np.nan, h=np.nan, base=np.nan)
    good = np.isfinite(prof)
    if good.sum() < 5 or expected is None:
        return out
    base = np.nanpercentile(prof, 5)
    p = np.where(good, prof - base, np.nan)
    k = np.arange(p.size)
    k0 = int(np.round(expected))
    win = (k >= k0 - 3) & (k <= k0 + 3) & good
    if not win.any():
        return out
    it = int(k[win][np.nanargmax(p[win])])
    h = p[it]
    # width: walk outwards while above half height and not climbing into a neighbour
    lo = it
    while lo > 0 and np.isfinite(p[lo - 1]) and p[lo - 1] >= h / 2 and p[lo - 1] <= 1.05 * p[lo]:
        lo -= 1
    hi = it
    while hi < p.size - 1 and np.isfinite(p[hi + 1]) and p[hi + 1] >= h / 2 and p[hi + 1] <= 1.05 * p[hi]:
        hi += 1
    away = good & (np.abs(k - it) > 2 * star_fwhm)
    if away.any():
        inb = int(k[away][np.nanargmax(p[away])])
        out.update(nbr_off=float(inb - it), nbr_ratio=float(p[inb] / h))
    pos = np.clip(np.nan_to_num(p), 0, None)
    near = np.abs(k - it) <= star_fwhm
    out.update(it=float(it), off=float(it - expected), fwhm=float(hi - lo + 1), h=float(h),
               base=float(base), tfrac=float(pos[near].sum() / pos.sum()) if pos.sum() > 0 else np.nan)
    return out


def local_scatter(ws, res, use, width=0.5):
    """Robust residual scatter as a function of wavelength (MAD in bins of `width` A)."""
    edges = np.arange(ws[0], ws[-1] + width, width)
    idx = np.searchsorted(ws, edges)
    cen, sig = [], []
    for i0, i1 in zip(idx[:-1], idx[1:]):
        r = res[i0:i1][use[i0:i1]]
        if r.size >= 5:
            cen.append(np.median(ws[i0:i1]))
            sig.append(1.4826 * np.median(np.abs(r - np.median(r))))
    if len(cen) < 2:
        s0 = 1.4826 * np.median(np.abs(res[use]))
        return np.full(ws.size, s0)
    return np.interp(ws, cen, sig)


def fit_sky(ws, fs):
    """Cubic spline in wavelength (knots every SKY_KNOT A), iteratively clipped against the
    local scatter of the residuals.

    The scatter is measured as a function of wavelength because the noise on bright sky
    lines is much larger than in the continuum: clipping on one global scatter would reject
    the bright side of every line and leave the fit too low there. The fit is unweighted:
    PypeIt's ivarraw is computed from the noisy counts, so weighting by it favours pixels
    that fluctuated low and biases the sky low. Positive outliers are clipped harder than
    negative ones, because light from other stars on the slit only adds to the sky pixels.
    Returns the spline and the wavelength range it was fit over.
    """
    from scipy.interpolate import LSQUnivariateSpline
    use = np.ones(ws.size, bool)
    spl, rng = None, (np.nan, np.nan)
    for _ in range(6):
        lo, hi = ws[use][0], ws[use][-1]
        t = np.arange(lo + SKY_KNOT, hi - SKY_KNOT, SKY_KNOT)
        try:
            spl = LSQUnivariateSpline(ws[use], fs[use], t, k=3)
        except ValueError:
            return None, rng
        rng = (lo, hi)
        res = fs - spl(ws)
        sig = local_scatter(ws, res, use)
        new = (res < 3.0 * sig) & (res > -5.0 * sig)
        if np.array_equal(new, use):
            break
        use = new
    return spl, rng


def sky_subtract(rect, ivar2d, wave2d, wrow, it, star_fwhm):
    """Sky as a function of wavelength from pixels away from the target, then subtracted
    from every pixel at that pixel's own wavelength (this follows the tilted sky lines).
    Pixels outside the wavelength range covered by the fit are left undefined rather than
    extrapolated."""
    empty = (np.array([]), np.zeros((0, rect.shape[1])))
    sel = (wrow > CAT_WIN[0] - 5) & (wrow < CAT_WIN[1] + 5)
    if not sel.any():
        return empty
    r, w, v = rect[sel], wave2d[sel], ivar2d[sel]
    cols = np.arange(rect.shape[1])
    skycol = np.abs(cols - it) > 2 * star_fwhm
    ws, fs, iv = w[:, skycol].ravel(), r[:, skycol].ravel(), v[:, skycol].ravel()
    m = np.isfinite(ws) & np.isfinite(fs) & np.isfinite(iv) & (ws > 0) & (iv > 0)
    if m.sum() < 200:
        return empty
    o = np.argsort(ws[m])
    spl, (lo, hi) = fit_sky(ws[m][o], fs[m][o])
    if spl is None:
        return empty
    inside = np.isfinite(w) & (w >= lo) & (w <= hi)
    sky = np.full(w.shape, np.nan)
    sky[inside] = spl(w[inside])
    keep = (wrow[sel] >= CAT_WIN[0]) & (wrow[sel] <= CAT_WIN[1])
    return wrow[sel][keep], (r - sky)[keep]


def smooth(y, n=5):
    """Boxcar smooth that leaves gaps (NaN) as gaps instead of filling them with zeros."""
    y = np.asarray(y, float)
    good = np.isfinite(y)
    num = np.convolve(np.where(good, y, 0.0), np.ones(n), mode="same")
    den = np.convolve(good.astype(float), np.ones(n), mode="same")
    out = np.full(y.shape, np.nan)
    ok = den >= (n + 1) // 2
    out[ok] = num[ok] / den[ok]
    return out


def star_fwhm_from_spec1d(fn, det, cache):
    """Median FWHM (px) of objects PypeIt detected (not forced) on this mosaic."""
    f1d = os.path.join(os.path.dirname(fn), os.path.basename(fn).replace("spec2d_", "spec1d_", 1))
    if f1d not in cache:
        try:
            from pypeit.specobjs import SpecObjs
            cache[f1d] = SpecObjs.from_fitsfile(f1d, chk_version=False)
        except Exception as err:
            print(f"WARNING: cannot read {os.path.basename(f1d)} ({err}); using FWHM {FWHM_DEFAULT} px")
            cache[f1d] = None
    sobjs = cache[f1d]
    if sobjs is None:
        return FWHM_DEFAULT
    fw = []
    for i in range(len(sobjs)):
        o = sobjs[i]
        if o.DET != det or getattr(o, "MASKDEF_EXTRACT", False):
            continue
        if o.FWHM is not None and np.isfinite(o.FWHM) and o.FWHM > 0:
            fw.append(o.FWHM)
    return float(np.median(fw)) if len(fw) >= 5 else FWHM_DEFAULT


class CalibCache:
    """Load Tilts/WaveCalib calibration files per detector, once."""

    def __init__(self, calib_dir):
        self.dir = calib_dir
        self.cache = {}

    def _find(self, kind, det):
        import glob
        hits = sorted(glob.glob(os.path.join(self.dir, f"{kind}_*_{det}.fits*")))
        if not hits:
            raise FileNotFoundError(f"no {kind}_*_{det}.fits in {self.dir} (set --calib)")
        return hits[0]

    def get(self, det):
        if det not in self.cache:
            from pypeit.wavetilts import WaveTilts
            from pypeit.wavecalib import WaveCalib
            self.cache[det] = (WaveTilts.from_file(self._find("Tilts", det), chk_version=False),
                               WaveCalib.from_file(self._find("WaveCalib", det), chk_version=False))
        return self.cache[det]


def slit_wave(s2d, i, flex, det, calib):
    """Wavelength image of slit i; rebuilt from calibrations if PypeIt left it at zero."""
    slits = s2d.slits
    spat = int(slits.spat_id[i])
    slitmask = slits.slit_img(flexure=flex, slitidx=[i])
    thismask = slitmask == spat
    if np.any(s2d.waveimg[thismask] > 0):
        return s2d.waveimg, False
    wtilts, wcal = calib.get(det)
    tilts = s2d.tilts
    if tilts is None or not np.any(tilts[thismask] > 0):
        tflex = (0.0 if flex is None else flex) - (0.0 if wtilts.spat_flexure is None else wtilts.spat_flexure)
        tilts = wtilts.fit2tiltimg(slitmask, flexure=tflex)
    k = np.where(np.asarray(wcal.spat_ids) == spat)[0]
    wave = np.zeros(s2d.waveimg.shape, dtype=float)
    if k.size and wcal.wv_fits[k[0]] is not None and wcal.wv_fits[k[0]].pypeitfit is not None:
        wave[thismask] = wcal.wv_fits[k[0]].pypeitfit.eval(tilts[thismask])
    return wave, True


def report_empty(s2d, ok, waveimg, left, right, spat):
    """If a slit yields no usable pixels in the profile window, say which cut removed them."""
    x = np.arange(s2d.sciimg.shape[1])[None, :]
    on = (x >= left[:, None]) & (x <= right[:, None])
    inwin = (waveimg > PROF_WIN[0]) & (waveimg < PROF_WIN[1])
    if np.sum(on & ok & inwin) >= 100:
        return
    print(f"  [spat {spat}] slit pixels {on.sum()}, pixel-flag/ivar ok {np.sum(on & ok)}, "
          f"wave > 0 {np.sum(on & (waveimg > 0))}, wave in {PROF_WIN} {np.sum(on & inwin)}, "
          f"wave range {np.nanmin(np.where(on, waveimg, np.nan)):.0f}-"
          f"{np.nanmax(np.where(on, waveimg, np.nan)):.0f}")


def slit_data(s2d, i, flex, det, calib):
    slits = s2d.slits
    left, right, _ = slits.select_edges(flexure=flex)
    left_init, _, _ = slits.select_edges(initial=True, flexure=flex)
    left, right = left[:, i], right[:, i]
    ok = np.logical_not(s2d.bpmmask.flagged(flag=PIXFLAGS)) & (s2d.ivarraw > 0) & np.isfinite(s2d.sciimg)
    rect = rectify(s2d.sciimg, ok, left, right)
    ivar2d = rectify(s2d.ivarraw, ok, left, right)
    waveimg, rebuilt = slit_wave(s2d, i, flex, det, calib)
    report_empty(s2d, ok, waveimg, left, right, int(slits.spat_id[i]))
    wave2d = rectify(waveimg, waveimg > 0, left, right)
    wrow = np.nanmedian(wave2d, axis=1)
    specmid = left.size // 2
    tab = s2d.maskdef_designtab
    oidx = np.where(tab["MASKDEF_ID"] == slits.maskdef_id[i])[0][0]
    expected = None
    if slits.maskdef_objpos is not None and slits.maskdef_offset is not None:
        expected = (slits.maskdef_objpos[oidx] + slits.maskdef_offset
                    + left_init[specmid, i] - np.round(left[specmid]))
    name = str(tab["OBJNAME"][oidx]) if "OBJNAME" in tab.colnames else ""
    return rect, ivar2d, wave2d, wrow, expected, name, rebuilt


def spatial_profile(rect, wrow):
    sel = (wrow > PROF_WIN[0]) & (wrow < PROF_WIN[1])
    return np.nanmedian(rect[sel], axis=0)


def shade_cat(axis):
    for c in CAT:
        axis.axvline(c, color="C3", ls=":", lw=0.8)
        axis.axvspan(c * (1 + VRANGE[0] / C_KMS), c * (1 + VRANGE[1] / C_KMS), color="C3", alpha=0.12)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("spec2d", nargs="+")
    ap.add_argument("--maskdef", nargs="+", type=int, required=True)
    ap.add_argument("--outdir", default="inspect_slits")
    ap.add_argument("--calib", default=None,
                    help="Calibrations folder (default: ../Calibrations relative to the spec2d files)")
    a = ap.parse_args()
    calib = CalibCache(a.calib or os.path.join(os.path.dirname(os.path.abspath(a.spec2d[0])),
                                               os.pardir, "Calibrations"))

    from pypeit.spec2dobj import AllSpec2DObj
    os.makedirs(a.outdir, exist_ok=True)
    s1d_cache = {}

    data = {m: [] for m in a.maskdef}
    for fn in sorted(a.spec2d):
        allspec = AllSpec2DObj.from_fits(fn, chk_version=False)
        mt = re.search(r"DE\.\d+\.(\d+)", os.path.basename(fn))
        tag = mt.group(1) if mt else os.path.basename(fn)
        for det in allspec.detectors:
            s2d = allspec[det]
            slits = s2d.slits
            if slits.maskdef_id is None or s2d.maskdef_designtab is None:
                continue
            if not any(m in data for m in slits.maskdef_id):
                continue
            flex = s2d.sci_spat_flexure
            sfw = star_fwhm_from_spec1d(fn, det, s1d_cache)
            for i, mid in enumerate(slits.maskdef_id):
                if mid not in data:
                    continue
                rect, ivar2d, wave2d, wrow, exp_pos, name, rebuilt = slit_data(s2d, i, flex, det, calib)
                prof = spatial_profile(rect, wrow)
                st = target_stats(prof, exp_pos, sfw)
                flags = ", ".join(slits.bitmask.flagged_bits(slits.mask[i])) or "None"
                if rebuilt:
                    flags += " (wave rebuilt)"
                it = st["it"] if np.isfinite(st["it"]) else (exp_pos if exp_pos is not None else np.nanargmax(prof))
                w, sub = sky_subtract(rect, ivar2d, wave2d, wrow, it, sfw)
                k0 = int(np.round(it))
                spec = np.sum(sub[:, max(k0 - BAND, 0):k0 + BAND + 1], axis=1) if w.size else np.array([])
                data[mid].append(dict(exp=tag, det=det, spat=int(slits.spat_id[i]), name=name, flags=flags,
                                      exp_pos=exp_pos, prof=prof, st=st, sfw=sfw, w=w, sub=sub, spec=spec))

    print(f"{'maskdef':>8} {'objname':>10} {'exp':>6} {'det':>5} {'spat':>5} {'len':>4} {'design':>6} "
          f"{'t_off':>5} {'t_fwhm':>6} {'s_fwhm':>6} {'n_off':>5} {'n/t':>5} {'t_frac':>6}  flags")
    for mid, rows in data.items():
        if not rows:
            print(f"{mid:>8}  not found in the spec2d files")
            continue
        n = len(rows)
        fig, ax = plt.subplots(n + 1, 3, figsize=(15, 3.2 * (n + 1)), squeeze=False)
        stack = []
        for r_i, d in enumerate(rows):
            st = d["st"]
            e = np.nan if d["exp_pos"] is None else d["exp_pos"]
            print(f"{mid:>8} {d['name']:>10} {d['exp']:>6} {d['det']:>5} {d['spat']:>5} {d['prof'].size:>4} "
                  f"{e:6.1f} {st['off']:5.1f} {st['fwhm']:6.0f} {d['sfw']:6.1f} {st['nbr_off']:5.0f} "
                  f"{st['nbr_ratio']:5.2f} {st['tfrac']:6.2f}  {d['flags']}")
            a0 = ax[r_i, 0]
            k = np.arange(d["prof"].size)
            base = st["base"] if np.isfinite(st["base"]) else np.nanpercentile(d["prof"], 5)
            a0.plot(k, d["prof"] - base, "k-", drawstyle="steps-mid", label=f"slit {d['spat']}")
            if np.isfinite(st["it"]):
                kk = np.linspace(0, k[-1], 400)
                sig = d["sfw"] / 2.3548
                a0.plot(kk, st["h"] * np.exp(-0.5 * ((kk - st["it"]) / sig) ** 2), "C1--",
                        label=f"single star, FWHM {d['sfw']:.1f} px")
            if d["exp_pos"] is not None:
                a0.axvline(d["exp_pos"], color="C0", ls=":", label="design position")
            a0.set_title(f"{mid} {d['name']}  {d['exp']} {d['det']}  [{d['flags']}]", fontsize=9)
            a0.set_xlabel("pixels from left edge")
            a0.legend(fontsize=7)
            if d["w"].size == 0:
                ax[r_i, 1].set_title("no sky estimate in the Ca II triplet window", fontsize=9)
                continue
            v = np.nanpercentile(d["sub"], [5, 95])
            ax[r_i, 1].imshow(d["sub"].T, aspect="auto", origin="lower", cmap="gray", vmin=v[0], vmax=v[1],
                              extent=[d["w"][0], d["w"][-1], 0, d["sub"].shape[1]])
            ax[r_i, 1].set_title("CaT region, sky subtracted by wavelength", fontsize=9)
            ax[r_i, 2].plot(d["w"], smooth(d["spec"]), "k-", lw=0.8)
            shade_cat(ax[r_i, 2])
            ax[r_i, 2].set_title(f"1D, +-{BAND} px of target peak, 5-px smooth", fontsize=9)
            stack.append((d["w"], d["spec"]))
        ax[n, 0].axis("off")
        ax[n, 1].axis("off")
        if stack:
            w0 = stack[0][0]
            mean = np.nanmean([np.interp(w0, w[np.argsort(w)], s[np.argsort(w)]) for w, s in stack], axis=0)
            ax[n, 2].plot(w0, smooth(mean), "k-", lw=0.8)
            shade_cat(ax[n, 2])
            ax[n, 2].set_title("1D, mean of exposures (shaded: -600 to 0 km/s)", fontsize=9)
        fig.tight_layout()
        fig.savefig(os.path.join(a.outdir, f"slit_{mid}.png"), dpi=110)
        plt.close(fig)
    print(f"\nPNG files in {a.outdir}/")


if __name__ == "__main__":
    main()
