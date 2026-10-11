"""T8 check 2: decompose the T9-method S98 value (clean selection) and test what drives the difference with the
per-event method; T9 method on S32 (REVIEW_2 scratch numbers); error models for the T9 statistic.
T9's fast_detrend is imported read-only from the sealed verification (PYTHONDONTWRITEBYTECODE=1 must be set).
Writes only to OUT/_work/T8_pixel_recheck/out/."""
import json, os, sys
import numpy as np
W = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(W, "pkg_copy", "scripts"))
from common import ROOT, EPH, load, DATA
sys.path.insert(0, os.path.join(ROOT, "addenda/independent_verification_tic4206066_SEALED/scripts"))
from lcio import fast_detrend
from astropy.io import fits
ODIR = os.path.join(W, "out"); os.makedirs(ODIR, exist_ok=True)
EPH9 = {"A": (3.182786, 3998.7368, 1.86 / 24), "B": (11.132736, 4000.8682, 1.6 / 24)}   # t9_aperture_diff.py L7
EPH_AUD = {"A": (EPH["A"][0], EPH["A"][1], EPH["A"][2]), "B": (EPH["B"][0], EPH["B"][1], EPH["B"][2])}
CLEAN8 = [0, 1, 4, 5, 9, 10, 13, 14]   # epochs used by the per-event method (pkg out/after_decorrelated.json A.S98_ffi)


def t9_lcs(t, cube, eph):
    ny, nx = cube.shape[1:]; med = np.nanmedian(cube, 0)
    cy, cx = np.unravel_index(np.nanargmax(med[ny//2-2:ny//2+3, nx//2-2:nx//2+3]), (5, 5)); cy += ny//2-2; cx += nx//2-2
    yy, xx = np.mgrid[0:ny, 0:nx]; rr = np.hypot(yy - cy, xx - cx)
    ann = (rr >= 4) & (rr <= 5.5)
    bk = np.nanmedian(cube[:, ann], axis=1)
    excl = np.zeros(len(t), bool)
    for P, T0, dur in eph.values():
        excl |= np.abs((((t - T0) / P + 0.5) % 1 - 0.5) * P) < dur
    LC = {}
    for r in (1.0, 2.0, 3.0):
        ap = rr <= r; fl = np.nansum(cube[:, ap], axis=1) - bk * ap.sum(); fl /= np.nanmedian(fl)
        fd, k = fast_detrend(t, fl, 1.0, mask=excl); LC[r] = np.where(k, fd, np.nan)
    return LC, excl, int(ann.sum())


def t9_stat(t, LC, excl, eph, sig="A", a=1.0, b=2.0, seed=11, win_half=None, keep_epochs=None, ntrial=300):
    P, T0, dur = eph[sig]
    half = dur / 2 if win_half is None else win_half
    dl = LC[a] - LC[b]; ok = np.isfinite(dl); tt, dd, ee = t[ok], dl[ok], excl[ok]
    def sel(T):
        ph = ((tt - T) / P + 0.5) % 1 - 0.5; it = np.abs(ph * P) < half; oot = ~(np.abs(ph * P) < dur) & ~ee
        return ph, it, oot
    ph, it, oot = sel(T0)
    ep = np.round((tt - T0) / P).astype(int)
    if keep_epochs is not None:
        it = it & np.isin(ep, keep_epochs)
    base = np.median(dd[oot])
    d0 = base - np.mean(dd[it])
    per = {}
    for e in np.unique(ep[it]):
        m = it & (ep == e); per[int(e)] = {"n_it": int(m.sum()), "extra_ppm": float(1e6 * (base - dd[m].mean()))}
    rng = np.random.default_rng(seed)
    nl = []
    for _ in range(ntrial):
        T = T0 + rng.uniform(0.15, 0.85) * P
        _, it2, oot2 = sel(T)
        if keep_epochs is not None:
            ep2 = np.round((tt - T) / P).astype(int); it2 = it2 & np.isin(ep2, keep_epochs)
        nl.append(np.median(dd[oot2]) - np.mean(dd[it2]) if it2.sum() >= 3 else np.nan)
    nl = np.array(nl); nl = nl[np.isfinite(nl)]
    v = np.array([p["extra_ppm"] for p in per.values()])
    return {"extra_ppm": float(d0 * 1e6), "null_std_ppm": float(np.std(nl) * 1e6), "null_mean_ppm": float(np.mean(nl) * 1e6),
            "sigma_null": float(d0 / np.std(nl)), "n_it": int(it.sum()),
            "white_err_ppm": float(1e6 * np.std(dd[oot], ddof=1) / np.sqrt(it.sum())),
            "sigma_white": float(d0 * 1e6 / (1e6 * np.std(dd[oot], ddof=1) / np.sqrt(it.sum()))),
            "n_events": len(per), "event_mean_ppm": float(v.mean()) if len(v) else None,
            "event_scatter_err_ppm": float(v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else None,
            "per_event": per}


res = {}
# S98 clean selection, raw cube (background added back) as in s03b
tc_, Fc, bkg, _, _ = load("S98_ffi"); Fraw = Fc + bkg[:, None, None]
LC, excl, nann = t9_lcs(tc_, Fraw, EPH9)
res["S98_clean_T9eph"] = {
    "baseline_fullT14": t9_stat(tc_, LC, excl, EPH9),
    "restricted_to_clean8_events": t9_stat(tc_, LC, excl, EPH9, keep_epochs=CLEAN8),
    "inTransit_pm_T14over4": t9_stat(tc_, LC, excl, EPH9, win_half=EPH9["A"][2] / 4),
    "r2_minus_r3": t9_stat(tc_, LC, excl, EPH9, a=2.0, b=3.0),
    "B_r1_minus_r2": t9_stat(tc_, LC, excl, EPH9, sig="B"),
    "annulus_pixels": nann}
# audit ephemeris
LCa, excla, _ = t9_lcs(tc_, Fraw, EPH_AUD)
res["S98_clean_auditEph"] = {"baseline_fullT14": t9_stat(tc_, LCa, excla, EPH_AUD),
                             "restricted_to_clean8_events": t9_stat(tc_, LCa, excla, EPH_AUD, keep_epochs=CLEAN8)}
# S98 T9 selection (QUALITY == 0, 21x21) per-event breakdown
t0, F0, _, _, _ = load("S98_ffi", raw=True)
LC0, excl0, _ = t9_lcs(t0, F0, EPH9)
res["S98_T9sel_T9eph"] = {"baseline_fullT14": t9_stat(t0, LC0, excl0, EPH9),
                          "restricted_to_clean8_events": t9_stat(t0, LC0, excl0, EPH9, keep_epochs=CLEAN8)}
# T9's own 11x11 S98 cut-out, T9 selection: per-event breakdown of the published 226.6
p11 = os.path.join(ROOT, "addenda/independent_verification_tic4206066_SEALED/downloads/tesscut/tess-s0098-1-2_80.859827_-1.326964_11x11_astrocut.fits")
d = fits.open(p11)[1].data; t = d["TIME"].astype(float); cube = d["FLUX"].astype(float); q = d["QUALITY"]
g = (q == 0) & np.isfinite(t) & np.all(np.isfinite(cube.reshape(len(t), -1)), axis=1)
LC11, excl11, nann11 = t9_lcs(t[g], cube[g], EPH9)
res["S98_T9own11x11"] = {"baseline_fullT14_seed11": t9_stat(t[g], LC11, excl11, EPH9), "annulus_pixels_in_stamp": nann11}

# S32: T9 method on FFI (QUALITY==0 and clean) and 2-min TPF (QUALITY==0; FLUX already background-subtracted by SPOC)
for name, raw in (("S32_ffi", True), ("S32_ffi", False), ("S32_tpf2min", True)):
    if raw:
        tt, FF, _, _, _ = load(name, raw=True)
    else:
        tt, FF, bb, _, _ = load(name); FF = FF + bb[:, None, None]
    L, e, na = t9_lcs(tt, FF, EPH9)
    r = t9_stat(tt, L, e, EPH9); r.pop("per_event")
    res[f"{name}_{'Q0' if raw else 'clean'}_T9method"] = {**r, "annulus_pixels": na, "stamp_shape": list(FF.shape[1:])}

for k, v in res.items():
    for kk, vv in (v.items() if isinstance(v, dict) else []):
        if isinstance(vv, dict) and "extra_ppm" in vv:
            print(k, kk, {x: (round(y, 2) if isinstance(y, float) else y) for x, y in vv.items() if x != "per_event"})
    if "extra_ppm" in v:
        print(k, {x: (round(y, 2) if isinstance(y, float) else y) for x, y in v.items()})
json.dump(res, open(os.path.join(ODIR, "t8_t9_decomposition.json"), "w"), indent=1)
