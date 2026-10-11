"""Verifier: per-epoch coverage of Signal A with QUALITY==0 (S6, S32, S98) and gap structure near each transit."""
import os
import numpy as np
from astropy.io import fits
exec(open("v00_inspect.py").read().split("for s, fn in")[0])
T14 = 1.87 / 24
res = []
for s, fn in F.items():
    d = fits.open(os.path.join(D, fn))[1].data
    t = np.asarray(d["TIME"], float); q = np.asarray(d["QUALITY"]).astype(np.int64)
    col = "SYS_RM_FLUX" if s == 98 else "PDCSAP_FLUX"
    f = np.asarray(d[col], float)
    g = (q == 0) & np.isfinite(t) & np.isfinite(f)
    tg = t[g]
    cad = np.median(np.diff(tg))
    e0, e1 = int(np.floor((tg.min() - T0A) / PA)), int(np.ceil((tg.max() - T0A) / PA))
    for e in range(e0, e1 + 1):
        tc = T0A + e * PA
        if tc < tg.min() - 0.5 or tc > tg.max() + 0.5:
            continue
        inn = np.abs(tg - tc) < T14 / 2
        pre = (tg > tc - 0.4) & (tg < tc - T14 / 2 - 0.02)
        post = (tg < tc + 0.4) & (tg > tc + T14 / 2 + 0.02)
        dB = abs(tc - (T0B + np.round((tc - T0B) / PB) * PB)) * 24
        # nearest gap (>0.1 d) before/after the transit centre
        w = tg[np.abs(tg - tc) < 1.5]
        gaps = [(w[i], w[i + 1]) for i in range(len(w) - 1) if w[i + 1] - w[i] > 0.1]
        frac = inn.sum() * cad / T14
        # also: fraction of ALL cadences (any quality) near transit
        allin = (np.abs(t - tc) < T14 / 2).sum()
        res.append((s, e, tc, int(inn.sum()), round(frac, 2), int(pre.sum()), int(post.sum()), round(dB, 1), allin,
                    [(round(a, 3), round(b, 3)) for a, b in gaps]))
for r in res:
    print(r)
n_any = sum(r[3] > 0 for r in res); n_half = sum(r[4] >= 0.5 for r in res)
n_2s = sum((r[4] >= 0.5) and r[5] > 0 and r[6] > 0 for r in res)
print("any in-transit:", n_any, " >=50% T14:", n_half, " >=50% and baseline on both sides:", n_2s)
