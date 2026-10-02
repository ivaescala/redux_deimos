#!/usr/bin/env python
"""
skysub_qa.py  Science/spec2d_*.fits  [--out skysub_qa.csv] [--vsys -200] [options]

Does imperfectly subtracted unresolved galaxy light leak into the stellar spectra?

For every slit, the script builds the residual image
    chi = (science - sky model - object model) * sqrt(inverse variance)
on OBJECT-FREE pixels only (every extracted object, including forced extractions
and serendips, is masked to +/- MASK_FWHM x its spatial FWHM; slit edges are
trimmed). If the background model is right, chi is pure noise: median 0, width 1.

Per slit it reports:
  med_chi, std_chi          median / robust width of chi over the full spectrum
  med_chi_cat, std_chi_cat  the same in the Ca II triplet region (8450-8750 A)
  d_line                    median chi in the three CaT lines at the galaxy velocity
                            minus median chi in the neighbouring continuum bands
  d_line_err, d_line_sig    its uncertainty and significance

d_line is the direct test. Leftover galaxy light (under-subtraction) carries the
galaxy's CaT absorption: the lines sit lower than the continuum, so d_line < 0.
Over-subtraction gives d_line > 0. The signal should grow toward M32's centre.
Run with --vsys set to M31's local disk velocity to test M31 light separately.

Each slit is placed at its design position and binned by projected distance from
M32; the bin summaries pool d_line with inverse-variance weights.

Wavelengths: PypeIt spec2d wavelengths are vacuum and, with refframe = observed,
not heliocentric-corrected. The line windows (default +/-4 A, about +/-140 km/s)
absorb the heliocentric offset (< 30 km/s).
"""
import argparse
import os

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.stats import mad_std
from astropy.table import Table

M32 = SkyCoord(10.674270 * u.deg, 40.865170 * u.deg)
C_KMS = 299792.458
AIR2VAC = 1.000275                       # refractive index of air near 8500 A
CAT_AIR = np.array([8498.02, 8542.09, 8662.14])
# Armandroff & Zinn (1988) CaT continuum bands (air, rest frame)
CONT_AIR = np.array([[8474., 8484.], [8563., 8577.], [8619., 8642.], [8700., 8725.]])
CAT_REGION = (8450., 8750.)              # vacuum, observed
PLATESCALE = 0.1185                      # DEIMOS arcsec / unbinned pixel
MIN_PIX = 200                            # fewer object-free pixels -> no statistics


def windows(vsys, half_width):
    """Vacuum, Doppler-shifted CaT line windows and continuum bands."""
    z = 1 + vsys / C_KMS
    lines = CAT_AIR * AIR2VAC * z
    line_win = np.column_stack((lines - half_width, lines + half_width))
    cont_win = CONT_AIR * AIR2VAC * z
    return line_win, cont_win


def in_windows(wave, win):
    sel = np.zeros(wave.shape, dtype=bool)
    for lo, hi in win:
        sel |= (wave >= lo) & (wave <= hi)
    return sel


def slit_stats(chi, wave, good, line_win, cont_win):
    """Residual statistics for one slit from object-free pixels (good = True)."""
    nan = np.nan
    out = dict(n_free=int(good.sum()), med_chi=nan, std_chi=nan, med_chi_cat=nan,
               std_chi_cat=nan, d_line=nan, d_line_err=nan, d_line_sig=nan)
    if out["n_free"] < MIN_PIX:
        return out
    c, w = chi[good], wave[good]
    out["med_chi"], out["std_chi"] = np.median(c), mad_std(c)

    cat = (w >= CAT_REGION[0]) & (w <= CAT_REGION[1])
    if cat.sum() >= MIN_PIX:
        out["med_chi_cat"], out["std_chi_cat"] = np.median(c[cat]), mad_std(c[cat])

    li, co = in_windows(w, line_win), in_windows(w, cont_win)
    if li.sum() >= MIN_PIX // 4 and co.sum() >= MIN_PIX // 4:
        # standard error of a median: 1.2533 * sigma / sqrt(n)
        e_l = 1.2533 * mad_std(c[li]) / np.sqrt(li.sum())
        e_c = 1.2533 * mad_std(c[co]) / np.sqrt(co.sum())
        out["d_line"] = np.median(c[li]) - np.median(c[co])
        out["d_line_err"] = np.hypot(e_l, e_c)
        out["d_line_sig"] = out["d_line"] / out["d_line_err"]
    return out


def process_file(fn, line_win, cont_win, mask_fwhm, edge_trim):
    from pypeit.spec2dobj import AllSpec2DObj
    from pypeit.specobjs import SpecObjs

    f1d = os.path.join(os.path.dirname(fn), os.path.basename(fn).replace("spec2d_", "spec1d_", 1))
    sobjs = SpecObjs.from_fitsfile(f1d, chk_version=False) if os.path.exists(f1d) else None
    if sobjs is None:
        print(f"WARNING: no {os.path.basename(f1d)}; objects are not masked for {os.path.basename(fn)}")
    allspec = AllSpec2DObj.from_fits(fn, chk_version=False)

    rows = []
    for det in allspec.detectors:
        s2d = allspec[det]
        slits, tab = s2d.slits, s2d.maskdef_designtab
        if tab is None or slits.maskdef_id is None:
            continue
        flex = s2d.sci_spat_flexure
        left, right, _ = slits.select_edges(flexure=flex)
        slitimg = slits.slit_img(flexure=flex)
        chi_img = (s2d.sciimg - s2d.skymodel - s2d.objmodel) * np.sqrt(s2d.ivarmodel)
        okpix = (s2d.bpmmask.mask == 0) & (s2d.ivarmodel > 0) & (s2d.waveimg > 0)

        objs = [] if sobjs is None else [sobjs[i] for i in range(len(sobjs)) if sobjs[i].DET == det]
        fwhm_med = np.nanmedian([o.FWHM for o in objs]) if objs else 7.0

        for i, (spat_id, mid) in enumerate(zip(slits.spat_id, slits.maskdef_id)):
            bits = slits.bitmask.flagged_bits(slits.mask[i])
            if "BOXSLIT" in bits or "BADWVCALIB" in bits:
                continue
            hit = np.where(tab["MASKDEF_ID"] == mid)[0]
            if len(hit) == 0:
                continue
            x0 = max(int(np.floor(left[:, i].min())), 0)
            x1 = min(int(np.ceil(right[:, i].max())) + 1, slitimg.shape[1])
            x = np.arange(x0, x1)[None, :]
            box = np.s_[:, x0:x1]

            good = (slitimg[box] == spat_id) & okpix[box]
            good &= (x >= left[:, i, None] + edge_trim) & (x <= right[:, i, None] - edge_trim)
            on_slit = [o for o in objs if o.SLITID == spat_id]
            for o in on_slit:
                fw = o.FWHM if np.isfinite(o.FWHM) and o.FWHM > 0 else fwhm_med
                good &= np.abs(x - o.TRACE_SPAT[:, None]) > mask_fwhm * fw

            st = slit_stats(chi_img[box], s2d.waveimg[box], good, line_win, cont_win)
            c = SkyCoord(tab["SLITRA"][hit[0]] * u.deg, tab["SLITDEC"][hit[0]] * u.deg)
            length = np.median(right[:, i] - left[:, i]) * PLATESCALE * slits.binspat
            rows.append(dict(spec2d=os.path.basename(fn), det=det, spat_id=int(spat_id),
                             maskdef_id=int(mid), r_arcsec=c.separation(M32).arcsec,
                             slit_len_arcsec=length, n_obj=len(on_slit), **st))
    return rows


def summarize(t, edges):
    print(f"\n{'r [arcsec]':<12}{'nslit':>6}{'std_chi':>9}{'med_chi':>9}{'med_cat':>9}"
          f"{'d_line':>9}{'+/-':>7}{'sig':>6}")
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (t["r_arcsec"] >= lo) & (t["r_arcsec"] < hi) & np.isfinite(t["med_chi"])
        if not m.any():
            continue
        d, e = np.asarray(t["d_line"][m]), np.asarray(t["d_line_err"][m])
        ok = np.isfinite(d) & (e > 0)
        if ok.any():
            w = 1 / e[ok] ** 2
            dm, de = np.sum(w * d[ok]) / np.sum(w), 1 / np.sqrt(np.sum(w))
        else:
            dm, de = np.nan, np.nan
        print(f"{lo:>4.0f}-{hi:<7.0f}{m.sum():>6}{np.nanmedian(t['std_chi'][m]):>9.2f}"
              f"{np.nanmedian(t['med_chi'][m]):>9.2f}{np.nanmedian(t['med_chi_cat'][m]):>9.2f}"
              f"{dm:>9.3f}{de:>7.3f}{dm / de:>6.1f}")


def main():
    p = argparse.ArgumentParser(description="Object-free sky-residual QA vs distance from M32.")
    p.add_argument("spec2d", nargs="+")
    p.add_argument("--out", default="skysub_qa.csv")
    p.add_argument("--vsys", type=float, default=-200.0,
                   help="galaxy velocity for the CaT line windows, km/s (default: M32, -200)")
    p.add_argument("--line_hw", type=float, default=4.0, help="CaT line half-width, A (default 4)")
    p.add_argument("--mask_fwhm", type=float, default=2.0,
                   help="mask objects to +/- this many FWHM (default 2)")
    p.add_argument("--edge_trim", type=float, default=3.0, help="pixels trimmed at slit edges (default 3)")
    a = p.parse_args()

    line_win, cont_win = windows(a.vsys, a.line_hw)
    rows = []
    for fn in a.spec2d:
        rows += process_file(fn, line_win, cont_win, a.mask_fwhm, a.edge_trim)
    if not rows:
        raise SystemExit("No slits with mask-design information found.")
    t = Table(rows=rows)
    t.sort("r_arcsec")
    t.write(a.out, overwrite=True)

    summarize(t, np.array([0, 30, 60, 120, 240, 480, 1e4]))
    nolight = np.sum(t["n_free"] < MIN_PIX)
    print(f"\n{len(t)} slit-exposures; {nolight} with < {MIN_PIX} object-free pixels (no statistics)."
          f"\nd_line < 0: galaxy CaT absorption left in the spectra (under-subtraction);"
          f" > 0: over-subtraction.  Table -> {a.out}")


if __name__ == "__main__":
    main()
