"""Independent significance tests for signal B (11.13 d), written from scratch for this check.

 T1  Cross-sector link (S6 + S32, and S6 alone): statistic = maximum 1.6-h box S/N of the archival data at the times
     predicted from the Sector 98 search ephemeris, over a period scan of +-0.004 d (1601 steps). Null: random period in
     9-13.5 d (excluding +-0.1 d around P_B) and random epoch, same scan width in time, B's own archival transits removed
     from the null data. Signal A is masked in the detrending and removed (+-1 T14). FAP = (n >= obs + 1)/(N + 1); the
     number of distinct argmax windows is counted.
 T2  Rank test in S6 (distribution free): S/N of all independent 1.6-h boxes; are B's two predicted dips the two
     strongest? p = (K x 3)/C(N, 2) with K = positions allowed by the period scan for the first dip, 3 spacings.
 T3  Sector 98 alone: BLS (astropy) SDE of B in the QLP frozen-mask light curve and in own TESScut photometry
     (straylight removed), A removed. Null: B removed, flux circularly shifted within each continuous segment
     (gaps > 0.5 d), 300 trials, maximum SDE over 1-20 d.
Outputs: out/a3_signalB.json
"""
import json, os, sys, warnings
import numpy as np
from astropy.timeseries import BoxLeastSquares
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import a1_timing_depths as a1
from wotan import flatten

H = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(H, "out")
rng = np.random.default_rng(4206066)
P98, T98 = 11.132975671662091, 4000.8681717061445      # frozen S98 search values (candidate.json)
PB = a1.EPH["B"]["P"]; DUR = 1.6 / 24

def prep_A_only(lc):
    """detrend with only A masked; remove A transits; return t, f, weights (1/sigma^2 with per-light-curve sigma)"""
    t, f = lc["t"], lc["f"]; msk = a1.in_tr(t, "A", 1.5)
    _, tr = flatten(t[~msk], f[~msk], method="biweight", window_length=0.75, return_trend=True, break_tolerance=0.3)
    ok = np.isfinite(tr); fd = f / np.interp(t, t[~msk][ok], tr[ok]); s = 1.4826 * np.median(np.abs(fd - np.median(fd)))
    keep = np.isfinite(fd) & (fd < 1 + 4 * s) & ~a1.in_tr(t, "A", 1.0)
    return t[keep], fd[keep], np.full(keep.sum(), 1 / s ** 2)

def scan_max(t, f, w, T0, P0, dps, dur=DUR, nmin=4):
    """box S/N at the predicted times for every trial period; cumulative sums + searchsorted (t must be sorted)"""
    cw = np.concatenate([[0.0], np.cumsum(w)]); cwf = np.concatenate([[0.0], np.cumsum(w * f)]); cn = np.arange(len(t) + 1, dtype=float)
    P = P0 + dps; k0 = int(np.floor((t.min() - T0) / P.max())) - 1; k1 = int(np.ceil((t.max() - T0) / P.min())) + 1
    ks = np.arange(k0, k1 + 1); tc = T0 + ks[None, :] * P[:, None]                       # (ndp, ncyc)
    lo = np.searchsorted(t, tc - dur / 2); hi = np.searchsorted(t, tc + dur / 2)
    n = (cn[hi] - cn[lo]).sum(1); sw = (cw[hi] - cw[lo]).sum(1); swf = (cwf[hi] - cwf[lo]).sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        snr = np.where(n >= nmin, (1 - swf / sw) * np.sqrt(sw), -np.inf)
    j = int(np.argmax(snr))
    if not np.isfinite(snr[j]): return np.nan, None
    ph = ((t - T0) / P[j] + 0.5) % 1 - 0.5
    return float(snr[j]), (float(P[j]), np.abs(ph * P[j]) < dur / 2)

def link_test(name, lcs, ntrial=10000):
    T, F, W = zip(*[prep_A_only(lc) for lc in lcs]); t, f, w = np.concatenate(T), np.concatenate(F), np.concatenate(W)
    o = np.argsort(t); t, f, w = t[o], f[o], w[o]
    dps = np.linspace(-0.004, 0.004, 1601)
    obs, arg = scan_max(t, f, w, T98, P98, dps)
    keep = ~a1.in_tr(t, "B", 1.0); tn, fn, wn = t[keep], f[keep], w[keep]
    null = []; cfgs = set()
    while len(null) < ntrial:
        P0 = rng.uniform(9.0, 13.5)
        if abs(P0 - PB) < 0.1: continue
        T0 = T98 + rng.uniform(0, 1) * P0
        m, a = scan_max(tn, fn, wn, T0, P0, dps * (P0 / P98))
        if not np.isfinite(m): continue
        null.append(m)
        P, inn = a; cyc = np.unique(np.round((tn[inn] - T0) / P)); cfgs.add(tuple(int(round((T0 + c * P) / (DUR / 2))) for c in cyc))
    null = np.array(null); n_ge = int(np.sum(null >= obs))
    res = dict(name=name, obs=float(obs), P_at_max=float(arg[0]), null_n=len(null), null_median=float(np.median(null)), null_p99=float(np.percentile(null, 99)),
               null_max=float(null.max()), n_null_ge_obs=n_ge, fap=float((n_ge + 1) / (len(null) + 1)), n_distinct_windows=len(cfgs))
    print(f"T1 {name:12s}: obs {obs:.2f} at P {arg[0]:.6f}; null median {res['null_median']:.2f} p99 {res['null_p99']:.2f} max {res['null_max']:.2f}; FAP {res['fap']:.2e} ({n_ge}/{len(null)}); distinct windows {len(cfgs)}", flush=True)
    return res

def rank_test(lc, name):
    t, f, w = prep_A_only(lc); nmin = 2 if lc["cad"] > 0.01 else 6
    cen = np.arange(t.min() + DUR / 2, t.max() - DUR / 2, DUR / 4); snr = np.full(len(cen), np.nan)
    for i, c in enumerate(cen):
        m = np.abs(t - c) < DUR / 2
        if m.sum() >= nmin: sw = w[m].sum(); snr[i] = (1 - np.sum(w[m] * f[m]) / sw) * np.sqrt(sw)
    ok = np.isfinite(snr); cen, snr = cen[ok], snr[ok]
    ep = np.arange(np.floor((t.min() - a1.EPH["B"]["T0"]) / PB), np.ceil((t.max() - a1.EPH["B"]["T0"]) / PB) + 1); tb = a1.EPH["B"]["T0"] + ep * PB
    tb = tb[(tb > t.min()) & (tb < t.max())]
    nearB = np.zeros(len(cen), bool)
    for x in tb: nearB |= np.abs(cen - x) < DUR
    peaks = []
    for i in np.argsort(-snr):
        if nearB[i]: continue
        if all(abs(cen[i] - p[0]) > DUR for p in peaks): peaks.append((float(cen[i]), float(snr[i])))
        if len(peaks) >= 6: break
    bvals = []
    for x in tb:
        m = np.abs(cen - x) < DUR / 2
        if not m.any(): bvals.append(dict(t_pred=float(x), snr=None, rank=None)); continue
        j = np.flatnonzero(m)[np.argmax(snr[m])]
        bvals.append(dict(t_pred=float(x), snr=float(snr[j]), rank=int(np.sum(np.array([p[1] for p in peaks]) > snr[j])) + 1))
    N = int((t.max() - t.min()) / DUR); K = 2 * 0.004 * 227 / DUR
    p = K * 3 / (N * (N - 1) / 2) if len(bvals) == 2 and all(b["rank"] == 1 for b in bvals) else None
    print(f"T2 {name}: B boxes " + ", ".join(f"{b['snr']} (rank {b['rank']})" for b in bvals) + f"; top others " + ", ".join(f"{p_[1]:.1f}@{p_[0]:.1f}" for p_ in peaks[:4]) + f"; N={N} K={K:.0f} p_rank={p}", flush=True)
    return dict(name=name, B_boxes=bvals, top_other=peaks, N_boxes=N, K=K, p_rank=p, sd_all=float(np.std(snr[~nearB])))

def s98_test(name, lc, ntrial=300):
    t, f, _ = prep_A_only(lc); base = t.max() - t.min(); durs = np.array([0.04, 0.055, 0.07, 0.085, 0.10, 0.13])
    periods = np.exp(np.arange(np.log(1.0), np.log(20.0), durs.min() / (base * 5)))
    def sde(tt, ff):
        err = np.full(len(tt), 1.4826 * np.median(np.abs(ff - np.median(ff)))); r = BoxLeastSquares(tt, ff, err).power(periods, durs, objective="likelihood")
        return (r.power - r.power.mean()) / r.power.std(), r
    s, r = sde(t, f); nearB = np.abs(periods / PB - 1) < 0.005; iB = np.flatnonzero(nearB)[np.argmax(s[nearB])]; iM = int(np.argmax(s))
    obs = float(s[iB]); keep = ~a1.in_tr(t, "B", 1.0); tn, fn = t[keep], f[keep]
    segs = np.split(np.arange(len(tn)), np.flatnonzero(np.diff(tn) > 0.5) + 1)
    null = []
    for j in range(ntrial):
        fs = fn.copy()
        for sg in segs:
            if len(sg) > 2: fs[sg] = np.roll(fn[sg], rng.integers(1, len(sg)))
        null.append(float(sde(tn, fs)[0].max()))
    null = np.array(null); n_ge = int(np.sum(null >= obs))
    res = dict(name=name, obs_SDE_B=obs, P_at_B=float(periods[iB]), depth_B_ppm=float(r.depth[iB] * 1e6), global_max_SDE=float(s[iM]), P_global_max=float(periods[iM]),
               B_is_global_max=bool(iB == iM), null_n=ntrial, null_median=float(np.median(null)), null_p95=float(np.percentile(null, 95)), null_max=float(null.max()), fap=float((n_ge + 1) / (ntrial + 1)), n_segments=len(segs))
    print(f"T3 {name}: B SDE {obs:.2f} at {periods[iB]:.4f} d (global max {s[iM]:.2f} at {periods[iM]:.3f} d); null median {res['null_median']:.2f} p95 {res['null_p95']:.2f} max {null.max():.2f}; FAP {res['fap']:.3f}", flush=True)
    return res

if __name__ == "__main__":
    S6, S32 = a1.load_lc("S6"), a1.load_lc("S32"); S6o, S32o = a1.own_phot(6, 2.5), a1.own_phot(32, 2.0)
    out = dict(P98_search=P98, T98_search=T98, P_B_joint=PB, dur_h=DUR * 24)
    out["T1_link"] = [link_test("pipe S6+S32", [S6, S32]), link_test("pipe S6", [S6]), link_test("pipe S32", [S32]), link_test("own S6+S32", [S6o, S32o]), link_test("own S6", [S6o])]
    out["T2_rank"] = [rank_test(S6, "pipe S6"), rank_test(S6o, "own S6"), rank_test(S32, "pipe S32")]
    S98 = a1.load_lc("S98"); S98o = a1.own_phot(98, 2.0)
    out["T3_s98"] = [s98_test("QLP S98 frozen mask", S98), s98_test("own S98 (straylight removed)", S98o)]
    json.dump(out, open(os.path.join(OUT, "a3_signalB.json"), "w"), indent=1)
    print("done")
