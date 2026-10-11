"""Replay of the old T9 method (independent verification, scripts/t9_aperture_diff.py; its own detrender lcio.fast_detrend
imported read-only) on (a) its original cadence selection (QUALITY == 0) and (b) the clean selection of this package
(pubprep a2 bitmask, S98 stray-light and QLP 2^29 intervals removed). Same apertures (r <= 1, 2, 3 px), same annulus
background (4-5.5 px), same 300-trial random-phase null, T9's own ephemeris. Purpose: separate the effect of the
method from the effect of the cadence selection."""
import json, os, sys
import numpy as np
from common import ROOT, OUT, load
sys.path.insert(0, os.path.join(ROOT, "addenda/independent_verification_tic4206066_SEALED/scripts"))
from lcio import fast_detrend   # read-only import of the sealed T9 helper
EPH9 = {"A": (3.182786, 3998.7368, 1.86 / 24), "B": (11.132736, 4000.8682, 1.6 / 24)}   # t9_aperture_diff.py L7

def t9(t, cube, seed=11):
    rng = np.random.default_rng(seed)
    ny, nx = cube.shape[1:]; med = np.nanmedian(cube, 0)
    cy, cx = np.unravel_index(np.nanargmax(med[ny//2-2:ny//2+3, nx//2-2:nx//2+3]), (5, 5)); cy += ny//2-2; cx += nx//2-2
    yy, xx = np.mgrid[0:ny, 0:nx]; rr = np.hypot(yy - cy, xx - cx); bk = np.nanmedian(cube[:, (rr >= 4) & (rr <= 5.5)], axis=1)
    excl = np.zeros(len(t), bool)
    for P, T0, dur in EPH9.values():
        excl |= np.abs((((t - T0) / P + 0.5) % 1 - 0.5) * P) < dur
    LC = {}
    for r in (1.0, 2.0, 3.0):
        ap = rr <= r; fl = np.nansum(cube[:, ap], axis=1) - bk * ap.sum(); fl /= np.nanmedian(fl)
        fd, k = fast_detrend(t, fl, 1.0, mask=excl); LC[r] = np.where(k, fd, np.nan)
    rec = {}
    for s, (P, T0, dur) in EPH9.items():
        for a, b in ((1.0, 2.0), (2.0, 3.0), (1.0, 3.0)):
            dl = LC[a] - LC[b]; ok = np.isfinite(dl); tt, dd, ee = t[ok], dl[ok], excl[ok]
            def dep(T):
                ph = ((tt - T) / P + 0.5) % 1 - 0.5; it = np.abs(ph * P) < dur / 2; oot = ~(np.abs(ph * P) < dur) & ~ee
                return np.median(dd[oot]) - np.mean(dd[it]) if it.sum() >= 3 else np.nan
            d0 = dep(T0); nl = np.array([dep(T0 + rng.uniform(0.15, 0.85) * P) for _ in range(300)]); nl = nl[np.isfinite(nl)]
            rec[f"{s}_r{a}_minus_r{b}"] = {"extra_depth_ppm": float(d0 * 1e6), "null_std_ppm": float(np.std(nl) * 1e6), "sigma": float(d0 / np.std(nl))}
    return rec

res = {}
for name in ("S98_ffi", "S6_ffi"):
    t, F, _, _, _ = load(name, raw=True); res[f"{name}_T9_selection"] = t9(t, F)
    # clean selection: take the raw (not background-subtracted) cube so that T9's own annulus background applies
    tc_, Fc, bkg, _, _ = load(name); Fraw = Fc + bkg[:, None, None]
    res[f"{name}_clean_selection"] = t9(tc_, Fraw)
    for k in (f"{name}_T9_selection", f"{name}_clean_selection"):
        v = res[k]["A_r1.0_minus_r2.0"]; print(k, "A r1-r2:", round(v["extra_depth_ppm"]), "±", round(v["null_std_ppm"]), f"({v['sigma']:.2f} sigma)")
# exact replay on T9's own 11x11 cut-outs (must reproduce t9_aperture_diff.json)
from astropy.io import fits
T9_CUTS = {"S98": "addenda/independent_verification_tic4206066_SEALED/downloads/tesscut/tess-s0098-1-2_80.859827_-1.326964_11x11_astrocut.fits",
           "S06": "addenda/independent_verification_tic4206066_SEALED/downloads/tesscut/tess-s0006-1-1_80.859827_-1.326964_11x11_astrocut.fits"}
for sec, p in T9_CUTS.items():
    d = fits.open(os.path.join(ROOT, p))[1].data; t = d["TIME"].astype(float); cube = d["FLUX"].astype(float); q = d["QUALITY"]
    g = (q == 0) & np.isfinite(t) & np.all(np.isfinite(cube.reshape(len(t), -1)), axis=1)
    res[f"{sec}_T9_own_11x11_cutout"] = {"file": p, **t9(t[g], cube[g])}
    v = res[f"{sec}_T9_own_11x11_cutout"]["A_r1.0_minus_r2.0"]; print(sec, "T9 own 11x11 cut-out A r1-r2:", round(v["extra_depth_ppm"], 1), "±", round(v["null_std_ppm"], 1))
json.dump(res, open(os.path.join(OUT, "t9_method_replay.json"), "w"), indent=1)
