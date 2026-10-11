"""Verifier: per-epoch breakdown of the deepest box in the secondary scan (phase ~0.368, width T14)."""
import numpy as np, json
from vcommon import *
PH = float(__import__("sys").argv[1]) if len(__import__("sys").argv) > 1 else 0.368
W = 1.8876 / 24
out = []
for s in (6, 32, 98):
    t, f, cad, ss = load(s)
    dA = t - (T0A + np.round((t - T0A) / PA) * PA); dB = t - (T0B + np.round((t - T0B) / PB) * PB)
    ok = (np.abs(dA) > 0.75 * T14A) & (np.abs(dB) > T14B)
    c0 = T0A + PH * PA
    for e in range(int(np.floor((t.min() - c0) / PA)), int(np.ceil((t.max() - c0) / PA)) + 1):
        tc = c0 + e * PA
        w = ok & (np.abs(t - tc) < 0.45) & segment_of(t, tc)
        inn = w & (np.abs(t - tc) < W / 2); oot = w & (np.abs(t - tc) > 0.75 * T14A)
        if inn.sum() * cad / 86400 < 0.5 * W or oot.sum() < 6:
            continue
        x = t - tc; c = np.polyfit(x[oot], f[oot], 2)
        r = f[inn] / np.polyval(c, x[inn]) - 1
        sig = np.std(f[oot] / np.polyval(c, x[oot]) - 1)
        # also linear-baseline value
        c1 = np.polyfit(x[oot], f[oot], 1); r1 = f[inn] / np.polyval(c1, x[inn]) - 1
        out.append(dict(sector=s, epoch=e, tc=round(tc, 3), n_in=int(inn.sum()), flux_quad_ppm=round(float(r.mean() * 1e6)),
                        flux_lin_ppm=round(float(r1.mean() * 1e6)), err_ppm=round(float(sig / np.sqrt(inn.sum()) * 1e6)),
                        span_pre=round(float(tc - t[w].min()), 2), span_post=round(float(t[w].max() - tc), 2)))
for o in out:
    print(o)
json.dump(out, open(os.path.join(OUT, f"dip_check_phase{PH}.json"), "w"), indent=1)
