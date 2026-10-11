"""Step 11b: Signal B cross-sector link — null test with (nearly) independent trials.

Problem found in s11 (audit bug C7): the null trials there shift only the S98 epoch (u*P98, u ~ U(0.1, 0.9)) at the
fixed S98 period. The dP scan (+-0.004 d) moves the predicted archival times by n*dP, i.e. +-0.91 d in S6 (n ~ -227)
and +-0.66 d in S32 (n ~ -163). Trials whose epochs differ by less than that window share the same data, so the
20 000 trials contain only ~5 (S6) to ~8 (S32) independent windows, and the null maximum is simply the largest feature
at any other phase (hence p95 = p99 = max in out/s11_nullB.json). An empirical FAP below ~1/(N_indep + 1) cannot be
claimed from that test.

Fix: randomise the period as well. P' ~ U(9.0, 13.5) d with |P' - P_B| > 0.1 d; T0' = T98 + u P', u ~ U(0, 1).
The dP grid is scaled by P'/P98 so that each null scan covers the same time window as the real scan. B's own archival
transits (+-1 T14 at the linked ephemeris) are removed from the data used by the null trials; A is handled exactly as in
s11 (masked in the detrending, its transits removed). The statistic is identical to s11 (maximum box S/N over the scan).
Reported: empirical FAP; the number of distinct argmax configurations (a direct count of independent null windows);
the per-position S/N distribution (to test the Gaussian noise model) and a Gaussian-model FAP for K independent box
positions per scan window.
"""
import json
import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
import lcio  # noqa
import s11_nullB as s11  # noqa  (same prep, same fast/brute scans; its __main__ block does not run on import)
from audit_common import save_json  # noqa

rng = np.random.default_rng(20260926)
EPH_CUR = json.load(open(os.path.join(lcio.AUD, "out", "ephem_current.json")))
EPH_B = {"B": dict(P=EPH_CUR["B"]["P"], T0=EPH_CUR["B"]["T0"], T14=EPH_CUR["B"]["T14"])}
PB = EPH_CUR["B"]["P"]
PLO, PHI = 9.0, 13.5
DUR = 1.6 / 24


def scan_grid(t, f, w, T0, dur, P0, grid):
    """s11.scan_stat on an arbitrary grid (same exact difference-array algorithm)."""
    old = s11.DPS
    try:
        s11.DPS = grid
        return s11.scan_stat(t, f, w, T0, dur, P0)
    finally:
        s11.DPS = old


def brute_grid(t, f, w, T0, dur, P0, grid):
    old = s11.DPS
    try:
        s11.DPS = grid
        return s11.scan_stat_brute(t, f, w, T0, dur, P0)
    finally:
        s11.DPS = old


def prepared(lcs):
    T, F, W = [], [], []
    for lc in lcs:
        t, f, w = s11.prep(lc)
        T.append(t), F.append(f), W.append(w)
    return np.concatenate(T), np.concatenate(F), np.concatenate(W)


def run(name, lcs, ntrial, n_brute=12):
    t, f, w = prepared(lcs)
    obs_scan = s11.scan_stat(t, f, w, s11.T98, DUR)
    obs = float(np.max(obs_scan))
    keep = ~lcio.in_transit(t, "B", 1.0, eph=EPH_B)           # null data: B's own archival transits removed
    tn, fn, wn = t[keep], f[keep], w[keep]
    # lever arm -> number of independent box positions inside one scan window (Gaussian-model look-elsewhere)
    ncyc = np.abs(np.round((t - s11.T98) / s11.P98))
    K = float(2 * 0.004 * ncyc.max() / DUR)
    Ps = []
    while len(Ps) < ntrial:
        p = rng.uniform(PLO, PHI)
        if abs(p - PB) > 0.1:
            Ps.append(p)
    Ps = np.array(Ps)
    us = rng.uniform(0, 1, ntrial)
    null = np.empty(ntrial)
    single = np.empty(ntrial)
    cfg = []
    maxdiff = 0.0
    for j in range(ntrial):
        P0, T0 = Ps[j], s11.T98 + us[j] * Ps[j]
        grid = s11.DPS * (P0 / s11.P98)
        sc = scan_grid(tn, fn, wn, T0, DUR, P0, grid)
        if j < n_brute:
            br = brute_grid(tn, fn, wn, T0, DUR, P0, grid)
            fin = np.isfinite(br) & np.isfinite(sc)
            assert np.array_equal(np.isfinite(br), np.isfinite(sc)), "finite masks differ"
            if fin.any():
                maxdiff = max(maxdiff, float(np.max(np.abs(br[fin] - sc[fin]))))
        fin = np.isfinite(sc)
        if not fin.any():
            null[j], single[j] = np.nan, np.nan
            cfg.append(None)
            continue
        i = int(np.argmax(np.where(fin, sc, -np.inf)))
        null[j] = sc[i]
        single[j] = sc[rng.choice(np.flatnonzero(fin))]
        Pm = P0 + grid[i]
        ph = ((tn - T0) / Pm + 0.5) % 1 - 0.5
        inn = np.abs(ph * Pm) < DUR / 2
        cyc = np.unique(np.round((tn[inn] - T0) / Pm))
        cfg.append(tuple(int(round((T0 + c * Pm) / (DUR / 2))) for c in cyc))
    assert maxdiff < 1e-6, f"fast scan != brute ({maxdiff})"
    ok = np.isfinite(null)
    null, single, Ps_ok = null[ok], single[ok], Ps[ok]
    cfg = [c for c in cfg if c is not None]
    n_ge = int(np.sum(null >= obs))
    n_distinct = len(set(cfg))
    sd = float(np.std(single, ddof=1))
    fap_gauss_unit = float(1 - stats.norm.cdf(obs) ** K) if obs < 37 else 0.0
    fap_gauss_sd = float(stats.norm.sf(obs / sd) * K)       # first-order union bound (tiny values)
    res = dict(name=name, n_points=int(len(t)), n_points_null=int(len(tn)), obs_max_snr=obs, null_n=int(ok.sum()),
               P_range=[PLO, PHI], null_median=float(np.median(null)), null_p95=float(np.percentile(null, 95)),
               null_p99=float(np.percentile(null, 99)), null_p999=float(np.percentile(null, 99.9)),
               null_max=float(null.max()), n_null_ge_obs=n_ge, fap_empirical=float((n_ge + 1) / (ok.sum() + 1)),
               n_distinct_argmax_configurations=n_distinct,
               fap_floor_from_distinct=float(1 / (n_distinct + 1)),
               single_position_snr=dict(mean=float(np.mean(single)), sd=sd, frac_gt3=float(np.mean(single > 3)),
                                        frac_gt4=float(np.mean(single > 4)), max=float(single.max()),
                                        gauss_expect_gt3=float(stats.norm.sf(3)), gauss_expect_gt4=float(stats.norm.sf(4))),
               K_independent_positions_per_window=K,
               fap_gaussian_unit_sd=fap_gauss_unit, fap_gaussian_measured_sd_union=fap_gauss_sd,
               brute_validation_maxdiff=maxdiff, n_brute_checked=n_brute)
    print(f"{name:22s} obs {obs:.2f}; null(P,T0 random) median {res['null_median']:.2f} p99 {res['null_p99']:.2f} "
          f"p99.9 {res['null_p999']:.2f} max {res['null_max']:.2f}; N>=obs {n_ge}/{ok.sum()} FAP_emp {res['fap_empirical']:.1e}; "
          f"distinct windows {n_distinct}; single-pos S/N sd {sd:.2f} (>3: {res['single_position_snr']['frac_gt3']:.4f}); "
          f"K {K:.1f}; Gaussian FAP {fap_gauss_sd:.1e}", flush=True)
    return res


if __name__ == "__main__":
    ntr = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    S6p, S32p = lcio.load_product("TESSSPOC_S6_PDC"), lcio.load_product("SPOC2min_S32_PDC")
    S6o, S32o = lcio.load_own("tesscut_s6", 2.5), lcio.load_own("tesscut_s32", 2.0)
    out = [run("pipe S6+S32", [S6p, S32p], ntr), run("pipe S6 only", [S6p], ntr), run("pipe S32 only", [S32p], ntr),
           run("own S6+S32", [S6o, S32o], ntr), run("own S6 only", [S6o], ntr)]
    save_json(dict(method=__doc__.strip().split("\n")[0], P98=s11.P98, T98=s11.T98, dur_h=DUR * 24, results=out),
              "out/s11b_nullB_periods.json")
