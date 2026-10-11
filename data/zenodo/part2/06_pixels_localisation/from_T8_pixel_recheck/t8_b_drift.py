"""T8 check 6: recompute the SUMMARY's un-stored drift numbers for the two S32 B events (core-pixel range in 30-min
bins over +-8 h; binned flux-weighted centroid range). Same selection as s04/s05 (QUALITY == 0). Writes to out/."""
import json, os, sys
import numpy as np
W = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(W, "pkg_copy", "scripts"))
from common import EPH, load, centroid
t, F, bk, c, med = load("S32_tpf2min", raw=True); cx, cy = centroid(F, c)
P, T0, T14 = EPH["B"]; out = {}
for n in (-164, -162):
    tc = T0 + n * P; w = np.abs(t - tc) < 8 / 24; tt = t[w]; core = F[w][:, c[0], c[1]]
    nb = 15; k = len(tt) // nb * nb
    cb = core[:k].reshape(-1, nb).mean(1); xb = cx[w][:k].reshape(-1, nb).mean(1); yb = cy[w][:k].reshape(-1, nb).mean(1)
    out[str(n)] = {"core_range_pct_30min_bins": float(100 * (cb.max() - cb.min()) / np.median(cb)),
                   "core_range_pct_unbinned": float(100 * (core.max() - core.min()) / np.median(core)),
                   "centroid_range_mpix_30min_bins_x_y": [float(1e3 * (xb.max() - xb.min())), float(1e3 * (yb.max() - yb.min()))]}
print(out); json.dump(out, open(os.path.join(W, "out/t8_b_drift.json"), "w"), indent=1)
