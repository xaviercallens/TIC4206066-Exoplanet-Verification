"""Verifier: binned raw S98 SYS_RM_FLUX / DET_FLUX around epochs -1, 0, 1, 10 (30-min bins, ppm vs local median)."""
import os
import numpy as np
from astropy.io import fits
exec(open("v00_inspect.py").read().split("for s, fn in")[0])
d = fits.open(os.path.join(D, F[98]))[1].data
t = np.asarray(d["TIME"], float); q = np.asarray(d["QUALITY"]).astype(np.int64)
print("S98 time range all:", np.nanmin(t), np.nanmax(t))
g = (q == 0) & np.isfinite(t)
print("S98 q==0 range:", t[g].min(), t[g].max())
# segments of q==0 data
tg = t[g]; br = np.where(np.diff(tg) > 0.1)[0]
segs = list(zip(np.r_[tg[0], tg[br + 1]], np.r_[tg[br], tg[-1]]))
print("q==0 segments:", [(round(a, 3), round(b, 3)) for a, b in segs])
for col in ("SYS_RM_FLUX", "DET_FLUX"):
    f = np.asarray(d[col], float)
    for e in (-1, 1, 10):
        tc = T0A + e * PA
        m = g & np.isfinite(f) & (np.abs(t - tc) < 0.8)
        tt, ff = t[m], f[m] / np.nanmedian(f[m])
        edges = np.arange(tc - 0.8, tc + 0.8001, 1 / 24)
        idx = np.digitize(tt, edges)
        row = []
        for i in range(1, len(edges)):
            sel = idx == i
            row.append("%+5.0f" % ((ff[sel].mean() - 1) * 1e6) if sel.sum() > 3 else "   . ")
        print(col, e, "tc=%.3f" % tc, " ".join(row))
