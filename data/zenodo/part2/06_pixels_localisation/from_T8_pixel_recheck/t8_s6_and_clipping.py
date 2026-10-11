"""T8 check 3: T9 method in S6 (and S98) - per-event breakdown, leave-one-event-out, and the REVIEW_2 claim that
the 11x11 vs 21x21 difference is caused by the stamp edge clipping T9's 4-5.5 px annulus.
Writes only to OUT/_work/T8_pixel_recheck/out/."""
import json, os, sys
import numpy as np
W = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(W, "pkg_copy", "scripts"))
from common import ROOT, load
sys.path.insert(0, os.path.join(ROOT, "addenda/independent_verification_tic4206066_SEALED/scripts"))
from lcio import fast_detrend
from astropy.io import fits
sys.path.insert(0, W)
from t8_t9_decomposition_lib import t9_stat, EPH9   # noqa: E402

ODIR = os.path.join(W, "out")


def lcs(t, cube, clip11=False):
    ny, nx = cube.shape[1:]; med = np.nanmedian(cube, 0)
    cy, cx = np.unravel_index(np.nanargmax(med[ny//2-2:ny//2+3, nx//2-2:nx//2+3]), (5, 5)); cy += ny//2-2; cx += nx//2-2
    yy, xx = np.mgrid[0:ny, 0:nx]; rr = np.hypot(yy - cy, xx - cx)
    ann = (rr >= 4) & (rr <= 5.5)
    if clip11:   # keep only annulus pixels that lie inside the central 11x11 sub-array (the T9 stamp)
        inside = (np.abs(yy - ny//2) <= 5) & (np.abs(xx - nx//2) <= 5); ann &= inside
    bk = np.nanmedian(cube[:, ann], axis=1)
    excl = np.zeros(len(t), bool)
    for P, T0, dur in EPH9.values():
        excl |= np.abs((((t - T0) / P + 0.5) % 1 - 0.5) * P) < dur
    LC = {}
    for r in (1.0, 2.0, 3.0):
        ap = rr <= r; fl = np.nansum(cube[:, ap], axis=1) - bk * ap.sum(); fl /= np.nanmedian(fl)
        fd, k = fast_detrend(t, fl, 1.0, mask=excl); LC[r] = np.where(k, fd, np.nan)
    return LC, excl, int(ann.sum())


res = {}
for sec, name, p11 in (("S06", "S6_ffi", "tess-s0006-1-1_80.859827_-1.326964_11x11_astrocut.fits"),
                       ("S98", "S98_ffi", "tess-s0098-1-2_80.859827_-1.326964_11x11_astrocut.fits")):
    d = fits.open(os.path.join(ROOT, "addenda/independent_verification_tic4206066_SEALED/downloads/tesscut", p11))[1].data
    t = d["TIME"].astype(float); cube = d["FLUX"].astype(float); q = d["QUALITY"]
    g = (q == 0) & np.isfinite(t) & np.all(np.isfinite(cube.reshape(len(t), -1)), axis=1)
    L11, e11, n11 = lcs(t[g], cube[g])
    t21, F21, _, _, _ = load(name, raw=True)
    L21, e21, n21 = lcs(t21, F21)
    L21c, e21c, n21c = lcs(t21, F21, clip11=True)
    tc_, Fc, bkg, _, _ = load(name); Fraw = Fc + bkg[:, None, None]
    Lcl, ecl, _ = lcs(tc_, Fraw)
    r = {"own11x11_T9sel": t9_stat(t[g], L11, e11, EPH9), "own11x11_annulus_px": n11,
         "21x21_T9sel": t9_stat(t21, L21, e21, EPH9), "21x21_annulus_px": n21,
         "21x21_T9sel_annulus_clipped_to_11x11": t9_stat(t21, L21c, e21c, EPH9), "21x21_clipped_annulus_px": n21c,
         "21x21_clean": t9_stat(tc_, Lcl, ecl, EPH9)}
    # leave-one-event-out on the published configuration (own 11x11, T9 selection)
    base = r["own11x11_T9sel"]; eps = sorted(int(k) for k in base["per_event"])
    loo = {}
    for e in eps:
        keep = [x for x in eps if x != e]
        s = t9_stat(t[g], L11, e11, EPH9, keep_epochs=keep); loo[str(e)] = [round(s["extra_ppm"], 1), round(s["null_std_ppm"], 1)]
    r["own11x11_leave_one_out"] = loo
    res[sec] = r
    for k, v in r.items():
        if isinstance(v, dict) and "extra_ppm" in v:
            print(sec, k, round(v["extra_ppm"], 1), "+-", round(v["null_std_ppm"], 1), "null_mean", round(v["null_mean_ppm"], 1),
                  "white", round(v["white_err_ppm"], 1), "evscatter", v["event_scatter_err_ppm"] and round(v["event_scatter_err_ppm"], 1),
                  {e: (x["n_it"], round(x["extra_ppm"])) for e, x in v["per_event"].items()})
    print(sec, "LOO", loo)
json.dump(res, open(os.path.join(ODIR, "t8_s6_and_clipping.json"), "w"), indent=1)
