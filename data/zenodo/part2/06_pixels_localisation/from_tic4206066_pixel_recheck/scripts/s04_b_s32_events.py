"""Diagnosis of the apparent B/S32 'off-target' effect: the two S32 2-min B events, per event, with local
normalisation (median of the +-0.35 d out-of-transit flux); per-pixel depths, aperture depths, flux-weighted centroid
and POS_CORR shifts; momentum-dump distances from the project's T7 output."""
import json, os
import numpy as np
from astropy.io import fits
from common import ROOT, OUT, EPH, DATA, load, radii_grid, centroid

t, F, bk, c, med = load("S32_tpf2min", raw=True)   # QUALITY == 0, as for POS_CORR below (TPF FLUX is background-subtracted)
d = fits.open(os.path.join(ROOT, DATA["S32_tpf2min"]))[1].data
q = d["QUALITY"]; good = (q == 0) & np.isfinite(d["TIME"]) & np.all(np.isfinite(d["FLUX"].reshape(len(q), -1)), 1)
pc1 = d["POS_CORR1"][good].astype(float); pc2 = d["POS_CORR2"][good].astype(float)
P, T0, T14 = EPH["B"]; PA, T0A, T14A = EPH["A"]
rr = radii_grid(F.shape[1:], c); cx, cy = centroid(F, c)
t7 = json.load(open(os.path.join(ROOT, "addenda/independent_verification_tic4206066_SEALED/out/t7_momentum_dumps.json")))
out = {"centre_pix": c, "momentum_dumps_S32": t7["S32"], "events": {}}
def ws(lc, tt, tc):
    inn = np.abs(tt - tc) < 0.25 * T14; oot = (np.abs(tt - tc) > 0.75 * T14)
    nA_ = np.round((tt - T0A) / PA); oot &= np.abs(tt - (T0A + nA_ * PA)) > 0.75 * T14A   # mask Signal A transits
    base = np.nanmedian(lc[oot])
    return [float(1e6 * (1 - np.nanmean(lc[inn]) / base)), float(1e6 * np.nanstd(lc[oot] / base) / np.sqrt(inn.sum()))]
for n in np.unique(np.round((t - T0) / P)):
    tc = T0 + n * P; w = np.abs(t - tc) < 0.35
    if w.sum() < 50: continue
    tt, Fw = t[w], F[w]; inn = np.abs(tt - tc) < 0.25 * T14; oot = np.abs(tt - tc) > 0.75 * T14
    nA_ = np.round((tt - T0A) / PA); oot &= np.abs(tt - (T0A + nA_ * PA)) > 0.75 * T14A   # mask Signal A transits
    nA = np.round((tc - T0A) / PA)
    e = {"t_mid_BTJD": float(tc), "hours_after_nearest_A_mid": float((tc - (T0A + nA * PA)) * 24),
         "nearest_dump_h": float(min(abs(tc - x) for x in t7["S32"]["dump_times"]) * 24),
         "apertures": {str(r): ws(Fw[:, rr <= r].sum(1), tt, tc) for r in (0.5, 1.0, 1.5, 2.0, 3.0)},
         "centroid_shift_mpix": {"x": [float((cx[w][inn].mean()-cx[w][oot].mean())*1e3), float(cx[w][oot].std()/np.sqrt(inn.sum())*1e3)],
                                 "y": [float((cy[w][inn].mean()-cy[w][oot].mean())*1e3), float(cy[w][oot].std()/np.sqrt(inn.sum())*1e3)]},
         "poscorr_shift_mpix": [float((pc1[w][inn].mean()-pc1[w][oot].mean())*1e3), float((pc2[w][inn].mean()-pc2[w][oot].mean())*1e3)],
         "core_pixel_range_pct_over_window": float(100*(np.nanmax(Fw[:, c[0], c[1]])-np.nanmin(Fw[:, c[0], c[1]]))/np.nanmedian(Fw[:, c[0], c[1]])),
         "pixels": {f"{j-c[0]:+d},{i-c[1]:+d}": ws(Fw[:, j, i], tt, tc) for j in range(c[0]-2, c[0]+3) for i in range(c[1]-2, c[1]+3)}}
    out["events"][str(int(n))] = e
    print(int(n), {k: round(v[0]) for k, v in e["apertures"].items()}, "dump h", round(e["nearest_dump_h"], 1))
json.dump(out, open(os.path.join(OUT, "b_s32_events.json"), "w"), indent=1)
