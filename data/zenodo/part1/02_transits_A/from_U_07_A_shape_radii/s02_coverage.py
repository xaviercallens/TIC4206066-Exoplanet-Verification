"""Count the Signal-A transits that actually have data in the joint-fit windows (per variant)."""
import sys, json
import numpy as np
from common import *

VARIANT = sys.argv[1] if len(sys.argv) > 1 else "base"
HALF_T14 = T14_A_H / 24.0 / 2.0
rows = []
for s in (6, 32, 98):
    z = np.load(os.path.join(OUT, f"lc_{VARIANT}_S{s}.npz"))
    t, cad = z["t"], float(z["cad"])
    dtB, _ = phase_dt(t, P_B, T0_B)
    keepB = np.abs(dtB) > T14_B_H / 24.0          # drop +-1 T14_B around B transits
    e0, e1 = int(np.floor((t.min() - T0_A) / P_A)), int(np.ceil((t.max() - T0_A) / P_A))
    for e in range(e0, e1 + 1):
        tc = T0_A + e * P_A
        if tc < t.min() - 0.3 or tc > t.max() + 0.3:
            continue
        nexp_in = 2 * HALF_T14 / (cad / 86400.0)
        inn = (np.abs(t - tc) < HALF_T14) & keepB
        win = (np.abs(t - tc) < 0.3) & keepB
        pre = (t > tc - 0.3) & (t < tc - HALF_T14) & keepB
        post = (t > tc + HALF_T14) & (t < tc + 0.3) & keepB
        # B overlap
        dB = np.min(np.abs(tc - (T0_B + np.round((tc - T0_B) / P_B) * P_B)))
        rows.append(dict(sector=s, epoch=e, tc=round(tc, 4), n_in=int(inn.sum()), frac_in=round(float(inn.sum() / nexp_in), 3),
                         n_win=int(win.sum()), n_pre=int(pre.sum()), n_post=int(post.sum()), dist_to_B_h=round(float(dB * 24), 2)))
n_any = int(sum(r["n_in"] > 0 for r in rows))
n_half = int(sum(r["frac_in"] >= 0.5 for r in rows))
n_full = int(sum(r["frac_in"] >= 0.9 for r in rows))
for r in rows:
    print(r)
print(f"variant={VARIANT}: epochs with any in-transit point {n_any}; >=50% coverage {n_half}; >=90% coverage {n_full}")
json.dump(dict(rows=rows, n_any=n_any, n_half=n_half, n_full=n_full), open(os.path.join(OUT, f"coverage_{VARIANT}.json"), "w"), indent=1)
