"""Step 7c: robustness of the S98-only false-alarm test of signal B (s07b) in the discovery (frozen-mask) data.
Three statistics at B's peak and as null maxima: BLS SDE, BLS log-likelihood power, depth S/N (astropy depth_snr).
Three nulls, all with A's and B's transits removed:
  (a) circular shift of each continuous segment (as s07b)          -> wrap points could add steps (conservative side)
  (b) block bootstrap: 0.5-d blocks re-ordered at random inside each segment (no single wrap point, keeps < 0.5 d noise)
  (c) single realisations: the B-free data unshifted, and the inverted light curve (2 - f, B kept) -> 'anti-transits'
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import lcio  # noqa
import s07b_s98_B_significance as sb  # noqa
from audit_common import save_json  # noqa

rng = np.random.default_rng(20260928)


def stats(t, f, periods, iB=None):
    err = np.full(len(t), 1.4826 * np.median(np.abs(f - np.median(f))))
    r = sb.BoxLeastSquares(t, f, err).power(periods, sb.DUR, objective="likelihood")
    pw = r.power
    sde = (pw - pw.mean()) / pw.std()
    snr = r.depth_snr
    if iB is None:
        return dict(SDE=float(sde.max()), power=float(pw.max()), depth_snr=float(np.nanmax(snr)))
    return dict(SDE=float(sde[iB]), power=float(pw[iB]), depth_snr=float(snr[iB]), P=float(periods[iB]), depth_ppm=float(r.depth[iB] * 1e6))


def block_boot(t, f, segs, block=0.5):
    fs = f.copy()
    for s in segs:
        ts = t[s]
        b = np.floor((ts - ts[0]) / block).astype(int)
        blocks = [s[b == k] for k in np.unique(b)]
        order = rng.permutation(len(blocks))
        vals = np.concatenate([f[blocks[k]] for k in order])
        fs[s] = vals[:len(s)] if len(vals) >= len(s) else np.resize(vals, len(s))
    return fs


t, f = sb.prep(lcio.load_product("QLP_S98_SYS", "frozen"))
base = t.max() - t.min()
periods = np.exp(np.arange(np.log(sb.PMIN), np.log(sb.PMAX), sb.DUR.min() / (base * 5)))
err = np.full(len(t), 1.4826 * np.median(np.abs(f - np.median(f))))
r0 = sb.BoxLeastSquares(t, f, err).power(periods, sb.DUR, objective="likelihood")
sde0 = (r0.power - r0.power.mean()) / r0.power.std()
nearB = np.abs(periods / sb.PB - 1) < 0.005
iB = int(np.flatnonzero(nearB)[np.argmax(sde0[nearB])])
obs = stats(t, f, periods, iB)
keep = ~lcio.in_transit(t, "B", 1.0, eph=sb.EPH_B)
tn, fn = t[keep], f[keep]
segs = sb.segments(tn)
single = dict(B_free_unshifted=stats(tn, fn, periods), inverted_with_B=stats(t, 2 - f, periods))
N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
res = {}
for name, maker in (("circular_shift", lambda: np.concatenate([np.roll(fn[s], rng.integers(1, len(s))) for s in segs])),
                    ("block_bootstrap_0.5d", lambda: block_boot(tn, fn, segs))):
    vals = [stats(tn, maker(), periods) for _ in range(N)]
    out = {}
    for k in ("SDE", "power", "depth_snr"):
        v = np.array([x[k] for x in vals])
        n_ge = int(np.sum(v >= obs[k]))
        out[k] = dict(obs=obs[k], null_median=float(np.median(v)), null_p95=float(np.percentile(v, 95)), null_max=float(v.max()),
                      n_ge=n_ge, fap=float((n_ge + 1) / (N + 1)))
    res[name] = out
    print(name, {k: (round(o["obs"], 2), round(o["null_median"], 2), round(o["null_p95"], 2), round(o["fap"], 3)) for k, o in out.items()}, flush=True)
print("observed B:", {k: round(v, 3) for k, v in obs.items()})
print("single realisations:", {k: {kk: round(vv, 2) for kk, vv in v.items()} for k, v in single.items()})
save_json(dict(doc=__doc__.strip().split("\n")[0], observed_B=obs, single_realisations=single, nulls=res, n_trials=N),
          "out/s07c_s98_B_null_variants.json")
