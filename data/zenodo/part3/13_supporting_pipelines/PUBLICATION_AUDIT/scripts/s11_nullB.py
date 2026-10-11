"""Step 11: Signal B cross-sector link and its null test, reproduced independently with many more trials.

Test statistic (as in deep_followup signal_11d.py): for trial period P' = P98 + dP, dP in [-0.004, +0.004] d
(1601 steps), the weighted box depth S/N of the archival data at the phases predicted from the S98 epoch; the
statistic is the MAXIMUM S/N over the dP scan. Null: the same scan with the S98 epoch shifted by u*P98,
u ~ U(0.1, 0.9) (20 000 trials). FAP = (N_null >= obs + 1)/(N + 1).
FAIRNESS: archival light curves are detrended with ONLY signal A masked (the detrending must not know where B is,
otherwise the real phase would be protected from the trend and favoured over null phases).
Variants: data set (pipeline / own photometry), sectors (S6+S32, S6, S32), box duration (1.3/1.6/2.0 h), S98 epoch
source (frozen BLS epoch / LEO trapezoid epoch), and the S/N at dP = 0 exactly (no scan, 'local').
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import lcio  # noqa
from audit_common import AUD, save_json  # noqa

P98 = 11.132975671662091          # frozen S98 search period (candidate.json)
T98 = 4000.8681717061445          # frozen S98 search epoch
T98_LEO = 4000.865724493904       # LEO-vetter trapezoid epoch (candidate.json leo_vetter.trap_epo)
EPH_A_ONLY = {"A": dict(lcio.EPH["A"]), "B": dict(P=1e9, T0=-1e9, T14=1e-9)}   # B unknown to the detrending
rng = np.random.default_rng(20260925)
DPS = np.linspace(-0.004, 0.004, 1601)


def prep(lc):
    d = lcio.detrend(lc, window=0.75, eph=EPH_A_ONLY)
    keep = ~lcio.in_transit(d["t"], "A", 1.0)          # remove A's transits (+-1 T14), as in the original test
    return d["t"][keep], d["fd"][keep], np.full(keep.sum(), 1 / d["sig"] ** 2)


def scan_stat_brute(t, f, w, T0, dur, P0=P98):
    """reference implementation: S/N for every dP (used to validate the fast version)."""
    snr = np.full(len(DPS), -np.inf)
    for i, dP in enumerate(DPS):
        P = P0 + dP
        inn = np.abs((((t - T0) / P + 0.5) % 1 - 0.5) * P) < dur / 2
        if inn.sum() < 5:
            continue
        sw = w[inn].sum()
        snr[i] = (1 - np.sum(w[inn] * f[inn]) / sw) * np.sqrt(sw)
    return snr


def scan_stat(t, f, w, T0, dur, P0=P98):
    """fast exact version: for fixed cycle number n the in-transit condition |D0 - n dP| < dur/2 is an interval in dP,
    so each point adds to a contiguous range of the dP grid (difference-array accumulation)."""
    n = np.round((t - T0) / P0)
    D0 = (t - T0) - n * P0
    nz = n != 0
    n, D0, fw, ww = n[nz], D0[nz], (f * w)[nz], w[nz]
    a = (D0 - dur / 2) / n
    b = (D0 + dur / 2) / n
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    # exact strict-inequality membership on the actual grid values (np.isclose must NOT be used here: its default
    # rtol=1e-5 silently shrank ranges -> caught by the brute-force validation, see BUG_REGRESSION_TESTS.md)
    i_lo = np.searchsorted(DPS, lo, side="right")          # first grid index with DPS > lo
    i_hi = np.searchsorted(DPS, hi, side="left") - 1       # last grid index with DPS < hi
    ok = (i_hi >= 0) & (i_lo <= len(DPS) - 1) & (i_hi >= i_lo)
    i_lo, i_hi = np.clip(i_lo[ok], 0, len(DPS) - 1), np.clip(i_hi[ok], 0, len(DPS) - 1)
    L = len(DPS) + 1
    S_w, S_fw, S_n = np.zeros(L), np.zeros(L), np.zeros(L)
    np.add.at(S_w, i_lo, ww[ok]); np.add.at(S_w, i_hi + 1, -ww[ok])
    np.add.at(S_fw, i_lo, fw[ok]); np.add.at(S_fw, i_hi + 1, -fw[ok])
    np.add.at(S_n, i_lo, 1); np.add.at(S_n, i_hi + 1, -1)
    sw, sfw, sn = np.cumsum(S_w)[:-1], np.cumsum(S_fw)[:-1], np.cumsum(S_n)[:-1]
    snr = np.full(len(DPS), -np.inf)
    g = sn >= 5
    snr[g] = (1 - sfw[g] / sw[g]) * np.sqrt(sw[g])
    return snr


def fast_null(t, f, w, dur, ntrial, P0=P98, T0=T98):
    out = np.empty(ntrial)
    us = rng.uniform(0.1, 0.9, ntrial)
    for j, u in enumerate(us):
        out[j] = np.max(scan_stat(t, f, w, T0 + u * P0, dur, P0))
    return out, us


def run(name, lcs, dur, T0, ntrial):
    T, F, W = [], [], []
    for lc in lcs:
        t, f, w = prep(lc)
        T.append(t)
        F.append(f)
        W.append(w)
    t, f, w = np.concatenate(T), np.concatenate(F), np.concatenate(W)
    obs_scan = scan_stat(t, f, w, T0, dur)
    ref = scan_stat_brute(t, f, w, T0, dur)
    fin = np.isfinite(ref) & np.isfinite(obs_scan)
    maxdiff = float(np.max(np.abs(ref[fin] - obs_scan[fin]))) if fin.any() else 0.0
    assert maxdiff < 1e-6 and np.array_equal(np.isfinite(ref), np.isfinite(obs_scan)), f"fast scan != brute ({maxdiff})"
    i = int(np.argmax(obs_scan))
    local = float(obs_scan[np.argmin(np.abs(DPS))])
    null, us = fast_null(t, f, w, dur, ntrial, T0=T0)
    n_ge = int(np.sum(null >= obs_scan[i]))
    res = dict(name=name, n_points=int(len(t)), dur_h=dur * 24, T0_s98=T0, obs_max_snr=float(obs_scan[i]),
               P_at_max=float(P98 + DPS[i]), local_snr_dP0=local, null_n=int(ntrial), null_median=float(np.median(null)),
               null_p95=float(np.percentile(null, 95)), null_p99=float(np.percentile(null, 99)),
               null_p999=float(np.percentile(null, 99.9)), null_max=float(null.max()), n_null_ge_obs=n_ge,
               fap=float((n_ge + 1) / (ntrial + 1)),
               z_equiv=float(np.sqrt(2) * __import__("scipy").special.erfcinv(2 * (n_ge + 1) / (ntrial + 1))))
    print(f"{name:40s} dur={dur*24:.1f}h obs max S/N {res['obs_max_snr']:.2f} at P={res['P_at_max']:.6f} (local {local:.2f}); "
          f"null median {res['null_median']:.2f} p99 {res['null_p99']:.2f} p99.9 {res['null_p999']:.2f} max {res['null_max']:.2f}; "
          f"N>=obs {n_ge}/{ntrial} -> FAP {res['fap']:.2e} (~{res['z_equiv']:.1f} sigma one-sided)", flush=True)
    return res


if __name__ == "__main__":
    ntr = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    S6p, S32p = lcio.load_product("TESSSPOC_S6_PDC"), lcio.load_product("SPOC2min_S32_PDC")
    S6o, S32o = lcio.load_own("tesscut_s6", 2.5), lcio.load_own("tesscut_s32", 2.0)
    out = []
    D16 = 1.6 / 24
    out.append(run("pipe S6+S32 (reproduces original)", [S6p, S32p], D16, T98, ntr))
    out.append(run("pipe S6 only", [S6p], D16, T98, ntr))
    out.append(run("pipe S32 only", [S32p], D16, T98, ntr))
    out.append(run("own S6+S32", [S6o, S32o], D16, T98, ntr))
    out.append(run("own S6 only", [S6o], D16, T98, ntr))
    out.append(run("pipe S6+S32 dur 1.3h", [S6p, S32p], 1.3 / 24, T98, ntr // 4))
    out.append(run("pipe S6+S32 dur 2.0h", [S6p, S32p], 2.0 / 24, T98, ntr // 4))
    out.append(run("pipe S6+S32 LEO epoch", [S6p, S32p], D16, T98_LEO, ntr // 4))
    save_json(dict(P98=P98, T98=T98, dP_grid=[float(DPS[0]), float(DPS[-1]), len(DPS)], results=out), "out/s11_nullB.json")
