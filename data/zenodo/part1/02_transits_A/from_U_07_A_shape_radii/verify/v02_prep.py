"""Verifier light curves: per-transit LOCAL polynomial baselines fitted to out-of-transit data of the SAME contiguous
segment (gaps > 0.1 d split), transit masked at |dt| < 0.75*T14 (T14=1.87 h), B transits removed (|dt_B| < T14_B).
usage: v02_prep.py ORDER HALFWIN COL98   e.g. 2 0.45 SYS_RM_FLUX  (ORDER -1: per-window BIC choice of 0-3)
Output: out/lc_o{ORDER}_{COL98}.npz with t, f (detrended), sector, epoch, per-epoch flags (two-sided baseline)."""
import sys, json
import numpy as np
from vcommon import *
ORDER = int(sys.argv[1]); HW = float(sys.argv[2]); COL98 = sys.argv[3]
T, Fl, S, E = [], [], [], []
info = []
for s in (6, 32, 98):
    t, f, cad, ss = load(s, COL98 if s == 98 else None)
    e0, e1 = int(np.floor((t.min() - T0A) / PA)), int(np.ceil((t.max() - T0A) / PA))
    for e in range(e0, e1 + 1):
        tc = T0A + e * PA
        inn = np.abs(t - tc) < T14A / 2
        if inn.sum() * cad / 86400 < 0.5 * T14A:
            continue
        dB = t - (T0B + np.round((t - T0B) / PB) * PB)
        w = (np.abs(t - tc) < HW) & segment_of(t, tc) & (np.abs(dB) > T14B)
        oot = w & (np.abs(t - tc) > 0.75 * T14A)
        pre = oot & (t < tc); post = oot & (t > tc)
        span_pre = (tc - t[pre].min()) if pre.any() else 0.0
        span_post = (t[post].max() - tc) if post.any() else 0.0
        x = t - tc
        if ORDER >= 0:
            c = np.polyfit(x[oot], f[oot], ORDER)
        else:   # BIC over orders 0-3, sigma from the order-3 residual rms
            sg = np.std(f[oot] - np.polyval(np.polyfit(x[oot], f[oot], 3), x[oot]))
            bics = [np.sum((f[oot] - np.polyval(np.polyfit(x[oot], f[oot], o), x[oot])) ** 2) / sg ** 2 + (o + 1) * np.log(oot.sum()) for o in range(4)]
            c = np.polyfit(x[oot], f[oot], int(np.argmin(bics)))
            c = np.r_[np.zeros(3 - (len(c) - 1)), c]
        base = np.polyval(c, x[w])
        fd = f[w] / base
        rms = np.std(f[oot] / np.polyval(c, x[oot]))
        two = (span_pre > 0.15) and (span_post > 0.15)
        info.append(dict(sector=s, epoch=e, tc=tc, n=int(w.sum()), n_in=int((w & inn).sum()), span_pre=span_pre,
                         span_post=span_post, two_sided=bool(two), rms_ppm=float(rms * 1e6), slope_ppm_per_d=float(c[-2] * 1e6) if len(c) > 1 else 0.0, order=int(ORDER if ORDER >= 0 else 3 - int(np.sum(c[:3 - 0] == 0)))))
        T.append(t[w]); Fl.append(fd); S.append(np.full(w.sum(), s)); E.append(np.full(w.sum(), e))
T, Fl, S, E = map(np.concatenate, (T, Fl, S, E))
tag = f"o{ORDER}_hw{HW}_{COL98}"
np.savez(os.path.join(OUT, f"lc_{tag}.npz"), t=T, f=Fl, s=S, e=E)
json.dump(info, open(os.path.join(OUT, f"prep_{tag}.json"), "w"), indent=1)
for r in info:
    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
print(tag, "epochs", len(info), "two-sided", sum(r["two_sided"] for r in info))
