"""Grid posterior for B's k, b, rho (circular, fixed ephemeris, fixed LD 0.32/0.31), local linear baselines per event,
A subtracted with a fixed 522 ppm template. (i) rho ~ N(2.51,0.24); (ii) rho log-uniform 0.3-40. Rp with R*=0.623+-0.019."""
import json, os
import numpy as np
from v01_independent import *

ev = []
for sec in (6, 32, 98):
    t, f, q, cad = load(sec)
    dA = dtph(t, PA, TA) * PA
    f = f + 522e-6 * unit(PMA, dA, cad)
    for n, tc in btimes(t):
        w = np.abs(t - tc) < 0.3
        tt = t[w] - tc; ff = f[w]
        oot = np.abs(tt) > 0.05
        c = np.polyfit(tt[oot], ff[oot], 1); s = 1.4826 * np.median(np.abs(ff[oot] - np.polyval(c, tt[oot])))
        ev.append((sec, n, tt, ff, cad, s))
print([(e[0], e[1], len(e[2]), round(e[5] * 1e6)) for e in ev])

def chi2(k, b, rho):
    pm = tmodel(k, b, rho, PB)
    tot = 0.0
    for sec, n, tt, ff, cad, s in ev:
        ss = max(1, int(np.ceil(cad * 1440))); ss += (ss % 2 == 0)
        m = batman.TransitModel(pm, tt, supersample_factor=ss, exp_time=cad) if ss > 1 else batman.TransitModel(pm, tt)
        mod = m.light_curve(pm)
        X = np.vstack([mod, mod * tt]).T   # multiplicative linear baseline
        c, *_ = np.linalg.lstsq(X, ff, rcond=None)
        tot += np.sum(((ff - X @ c) / s) ** 2)
    return tot

ks = np.arange(0.022, 0.0501, 0.0007); bs = np.arange(0.0, 1.06, 0.025)
rhos_p = np.array([1.91, 2.15, 2.39, 2.51, 2.63, 2.87, 3.11]); wp = np.exp(-0.5 * ((rhos_p - 2.51) / 0.24) ** 2)
rhos_f = 10 ** np.linspace(np.log10(0.3), np.log10(40), 22)
res = {}
for label, rhos, wr in (("prior", rhos_p, wp), ("free_logU", rhos_f, np.ones(len(rhos_f)))):
    L = np.zeros((len(rhos), len(bs), len(ks)))
    for i, r in enumerate(rhos):
        for j, b in enumerate(bs):
            for l, k in enumerate(ks):
                if b > 1 + k - 0.002: L[i, j, l] = np.inf; continue
                L[i, j, l] = chi2(k, b, r)
        print(label, i, flush=True)
    P = np.exp(-0.5 * (L - np.nanmin(L[np.isfinite(L)]))) * wr[:, None, None]
    P[~np.isfinite(L)] = 0
    pk = P.sum(axis=(0, 1)); pk /= pk.sum(); pb = P.sum(axis=(0, 2)); pb /= pb.sum(); pr = P.sum(axis=(1, 2)); pr /= pr.sum()
    # Rp samples: draw k from grid marginal (with jitter) times R*
    kk = rng.choice(ks, 200000, p=pk) + rng.uniform(-0.00035, 0.00035, 200000)
    Rs = rng.normal(0.623, 0.019, 200000); Rp = kk * Rs * 109.076
    q = lambda x, p: np.interp(p, np.cumsum(x) / np.sum(x), None) if False else None
    def pct(vals, prob, p):
        c = np.cumsum(prob); return float(np.interp(p, c / c[-1], vals))
    res[label] = dict(k=[pct(ks, pk, 0.16), pct(ks, pk, 0.5), pct(ks, pk, 0.84)], b=[pct(bs, pb, 0.16), pct(bs, pb, 0.5), pct(bs, pb, 0.84)],
                      rho=[pct(rhos, pr, 0.16), pct(rhos, pr, 0.5), pct(rhos, pr, 0.84)],
                      Rp=[float(np.percentile(Rp, 16)), float(np.median(Rp)), float(np.percentile(Rp, 84))], chi2min=float(np.nanmin(L[np.isfinite(L)])))
    print(label, json.dumps(res[label]), flush=True)
json.dump(res, open(os.path.join(H, "v03_radius.json"), "w"), indent=1)
