"""Step 7b: false-alarm probability of signal B in Sector 98 ALONE (open item U18).

Observed statistic: the BLS SDE of the strongest peak within 0.5 % of P_B, from the same astropy BoxLeastSquares set-up
as s07 (log-likelihood objective, durations 0.04-0.13 d, log-uniform grid), on S98 data detrended with ONLY signal A
masked and A's transits removed (B is treated exactly like noise would be).
Null: the same data with B's transits (+-1 T14) removed, split into continuous segments (gaps > 0.5 d), and the flux of
every segment circularly shifted by a random offset. This keeps the time sampling, the gaps, the noise amplitude and
the within-segment red noise, but destroys any periodicity across segments. For each null realisation the MAXIMUM SDE
over the full period range is recorded (look-elsewhere included). FAP = (n_null >= obs + 1) / (N + 1).
Caveat: circular shifts move scattered-light ramps inside a segment; the null keeps their amount, not their timing.
"""
import json
import os
import sys
import time

import numpy as np
from astropy.timeseries import BoxLeastSquares

sys.path.insert(0, os.path.dirname(__file__))
import lcio  # noqa
from audit_common import save_json  # noqa

DUR = np.array([0.04, 0.055, 0.07, 0.085, 0.10, 0.13])
EPH_CUR = json.load(open(os.path.join(lcio.AUD, "out", "ephem_current.json")))
EPH_B = {"B": dict(EPH_CUR["B"])}
PB = EPH_CUR["B"]["P"]
EPH_A_ONLY = {"A": dict(lcio.EPH["A"]), "B": dict(P=1e9, T0=-1e9, T14=1e-9)}
rng = np.random.default_rng(20260927)
PMIN, PMAX = 1.0, 20.0


def prep(lc):
    d = lcio.detrend(lc, window=0.75, eph=EPH_A_ONLY)
    keep = ~lcio.in_transit(d["t"], "A", 1.0)
    return d["t"][keep], d["fd"][keep]


def bls_sde(t, f, periods):
    err = np.full(len(t), 1.4826 * np.median(np.abs(f - np.median(f))))
    r = BoxLeastSquares(t, f, err).power(periods, DUR, objective="likelihood")
    pw = r.power
    sde = (pw - np.mean(pw)) / np.std(pw)
    return r, sde


def segments(t, gap=0.5):
    br = np.flatnonzero(np.diff(t) > gap) + 1
    return np.split(np.arange(len(t)), br)


def run(name, lc, ntrial):
    t, f = prep(lc)
    base = t.max() - t.min()
    periods = np.exp(np.arange(np.log(PMIN), np.log(PMAX), DUR.min() / (base * 5)))
    r, sde = bls_sde(t, f, periods)
    nearB = np.abs(periods / PB - 1) < 0.005
    iB = np.flatnonzero(nearB)[np.argmax(sde[nearB])]
    iM = int(np.argmax(sde))
    obs = float(sde[iB])
    keep = ~lcio.in_transit(t, "B", 1.0, eph=EPH_B)
    tn, fn = t[keep], f[keep]
    segs = segments(tn)
    null, t1 = [], time.time()
    for j in range(ntrial):
        fs = fn.copy()
        for s in segs:
            if len(s) > 2:
                fs[s] = np.roll(fn[s], rng.integers(1, len(s)))
        _, sn = bls_sde(tn, fs, periods)
        null.append(float(sn.max()))
        if j == 0:
            print(f"  {name}: one null BLS takes {time.time() - t1:.1f} s", flush=True)
    null = np.array(null)
    n_ge = int(np.sum(null >= obs))
    res = dict(name=name, n_points=int(len(t)), n_segments_null=len(segs), n_periods=len(periods), P_range=[PMIN, PMAX],
               obs_SDE_B=obs, P_at_B=float(periods[iB]), depth_B_ppm=float(r.depth[iB] * 1e6),
               obs_global_max_SDE=float(sde[iM]), P_at_global_max=float(periods[iM]), B_is_global_max=bool(iB == iM),
               null_n=ntrial, null_median=float(np.median(null)), null_p95=float(np.percentile(null, 95)),
               null_p99=float(np.percentile(null, 99)), null_max=float(null.max()), n_null_ge_obs=n_ge,
               fap=float((n_ge + 1) / (ntrial + 1)))
    print(f"{name:22s} B SDE {obs:.2f} at {periods[iB]:.4f} d (global max {sde[iM]:.2f} at {periods[iM]:.3f} d); null max-SDE median "
          f"{res['null_median']:.2f} p95 {res['null_p95']:.2f} p99 {res['null_p99']:.2f} max {res['null_max']:.2f}; N>=obs {n_ge}/{ntrial} "
          f"-> FAP {res['fap']:.3f}", flush=True)
    return res


if __name__ == "__main__":
    ntr = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    out = [run("QLP S98 frozen mask", lcio.load_product("QLP_S98_SYS", "frozen"), ntr),
           run("QLP S98 default mask", lcio.load_product("QLP_S98_SYS", "default"), ntr),
           run("own S98 no straylight", lcio.load_own("tesscut_s98", 2.0, drop_stray=True), ntr)]
    save_json(dict(method=__doc__.strip().split("\n")[0], P_B=PB, results=out), "out/s07b_s98_B_significance.json")
