"""Verifier diagnostic: per-window BIC of the local out-of-transit baseline order (0-3), +-0.45 d, same segment,
|dt| > 0.75 T14 (transit excluded), B excluded. sigma per window from the order-3 residual rms."""
import numpy as np, json
from vcommon import *
rows = []
for s in (6, 32, 98):
    t, f, cad, ss = load(s)
    dB = t - (T0B + np.round((t - T0B) / PB) * PB)
    for e in range(int(np.floor((t.min() - T0A) / PA)), int(np.ceil((t.max() - T0A) / PA)) + 1):
        tc = T0A + e * PA
        if (np.abs(t - tc) < T14A / 2).sum() * cad / 86400 < 0.5 * T14A:
            continue
        w = (np.abs(t - tc) < 0.45) & segment_of(t, tc) & (np.abs(dB) > T14B) & (np.abs(t - tc) > 0.75 * T14A)
        x, y = t[w] - tc, f[w]
        sig = np.std(y - np.polyval(np.polyfit(x, y, 3), x))
        bic = [float(np.sum((y - np.polyval(np.polyfit(x, y, o), x)) ** 2) / sig ** 2 + (o + 1) * np.log(len(x))) for o in range(4)]
        rows.append(dict(sector=s, epoch=e, n=int(w.sum()), best_order=int(np.argmin(bic)), dBIC_1_minus_2=bic[1] - bic[2]))
for r in rows:
    print(r)
bo = [r["best_order"] for r in rows]
print({o: bo.count(o) for o in range(4)})
json.dump(rows, open(os.path.join(OUT, "bic_baseline_order.json"), "w"), indent=1)
