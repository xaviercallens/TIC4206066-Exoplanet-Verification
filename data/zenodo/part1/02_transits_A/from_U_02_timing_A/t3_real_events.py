"""Audit 2, step 1: per-transit mid-times of signal A in S6 (TESS-SPOC 30 min PDCSAP), S32 (SPOC 2 min PDCSAP),
S98 (QLP 200 s SYS_RM_FLUX, scattered-light intervals removed), B masked.

Error models per event: (a) formal white (Delta chi2 = 1, chi2 normalised by the window residual variance),
(b) prayer bead (all cyclic residual shifts), (d) iid residual bootstrap; (e) time-averaging beta per sector.
(c) injection-recovery is in t4_injections.py.

Usage: python t3_real_events.py [variant]
variants: base (default), ld_lo, ld_hi, rmB_manu, doc98 (S98 with documented QLP bits only = adds the events in the
2^29/2^30 intervals), fulldetrend (whole-sector wotan detrending instead of local chunks), W030 (window +-0.3 d).
"""
import json, sys, time
import numpy as np
import pandas as pd
from wotan import flatten
import t2_common as C

variant = sys.argv[1] if len(sys.argv) > 1 else "base"
rng = np.random.default_rng(20260925)
cfg = dict(C.CFG)
tplkw = {}
s98_mode = "clean"
if variant == "ld_lo":
    tplkw = dict(u1=0.2604, u2=0.3657)      # Claret 2017 logg 5.0
elif variant == "ld_hi":
    tplkw = dict(u1=0.5586, u2=0.1295)      # Claret 2017 logg 4.0
elif variant == "rmB_manu":
    cfg["rmB"] = 1.0 * C.T14_B
elif variant == "doc98":
    s98_mode = "doc"
elif variant == "W030":
    cfg["W"] = 0.30
elif variant == "shape_manu_b0":
    tplkw = dict(b=0.0)

tpl = C.Template(**tplkw)
T14 = tpl.t14
print(f"variant {variant}: template T14 = {T14*24:.3f} h, a/R* = {tpl.pm.a:.3f}", flush=True)

rows = []
lcs = {s: C.load_lc(s, s98_mode=s98_mode if s == 98 else "clean") for s in (6, 32, 98)}
prov = {s: {k: lcs[s][k] for k in ("file", "col", "sha256", "url", "utc", "bitmask", "n_good", "n_total", "cad")} | {"meta": lcs[s]["meta"]} for s in lcs}

P, T0 = C.P_A0, C.T0_A0
for it in range(2):
    rows = []
    for s, lc in lcs.items():
        if variant == "fulldetrend":
            t, f = lc["t"], lc["f"]
            n_all, tcs_all = C.a_epochs(t, P, T0)
            m = C.near(t, tcs_all, cfg["maskA"] * T14) | C.near(t, C.b_times(t), cfg["maskB"] * C.T14_B)
            tr = flatten(t, f, method="biweight", window_length=cfg["win_len"], edge_cutoff=0, break_tolerance=0.375,
                         return_trend=True, mask=m)[1]
            ok = np.isfinite(tr)
            fd_full = f / np.interp(t, t[ok], tr[ok])
        n_all, tcs_all = C.a_epochs(lc["t"], P, T0)
        for n, tc in zip(n_all, tcs_all):
            if variant == "fulldetrend":
                w = np.abs(lc["t"] - tc) < cfg["W"]
                w &= ~C.near(lc["t"], C.b_times(lc["t"]), cfg["rmB"])
                tt, yy = lc["t"][w], fd_full[w]
                if len(tt) < 5:
                    continue
                oot = np.abs(tt - tc) > 0.75 * T14
                med = np.median(yy[oot]); sig = 1.4826 * np.median(np.abs(yy[oot] - med))
                k = yy < med + cfg["clip"] * sig; tt, yy = tt[k], yy[k]
                cov = float((np.abs(tt - tc) < T14 / 2).sum() * lc["cad"] / T14)
                prep = (tt, yy, cov)
            else:
                prep = C.prepare_window(lc, tc, tcs_all, T14, cfg=cfg)
            if prep is None:
                continue
            tt, yy, cov = prep
            if cov < cfg["cov_set"] or len(tt) < C.n_min(lc, cfg):
                continue
            r = C.measure(tt, yy, tc, lc["cad"], tpl, rng, n_boot=300, cfg=cfg)
            r.pop("chi_curve")
            r.update(sector=s, epoch=int(n), tc_pred=float(tc), coverage=cov,
                     tc_obs=float(tc + r["dt"] / 1440.0), tc_obs_alt=float(tc + r["dt_alt"] / 1440.0))
            rows.append(r)
    df = pd.DataFrame(rows)
    sel = (df.coverage >= cfg["cov_tim"]) & (df.snr >= cfg["snr_tim"]) & (~df.edge)
    # re-centre on a white-error linear ephemeris of the selected events (window placement only)
    e = df.loc[sel, "e_white"].values / 1440.0
    A = np.vstack([np.ones(sel.sum()), df.loc[sel, "epoch"].values]).T
    Wm = np.diag(1 / e ** 2)
    cvec = np.linalg.solve(A.T @ Wm @ A, A.T @ Wm @ df.loc[sel, "tc_obs"].values)
    print(f" iteration {it}: n={len(df)} selected={int(sel.sum())} -> T0={cvec[0]:.5f} P={cvec[1]:.8f}", flush=True)
    T0, P = float(cvec[0]), float(cvec[1])

df["sel"] = sel
# ---------------- beta (time averaging) per sector from out-of-transit residuals of the whole detrended sector
betas = {}
for s, lc in lcs.items():
    t, f = lc["t"], lc["f"]
    n_all, tcs_all = C.a_epochs(t, P, T0)
    m = C.near(t, tcs_all, cfg["maskA"] * T14) | C.near(t, C.b_times(t), cfg["maskB"] * C.T14_B)
    tr = flatten(t, f, method="biweight", window_length=cfg["win_len"], edge_cutoff=0, break_tolerance=0.375,
                 return_trend=True, mask=m)[1]
    ok = np.isfinite(tr) & ~m
    r = f[ok] / tr[ok] - 1.0
    med = np.median(r); sg = 1.4826 * np.median(np.abs(r - med)); k = r < med + 4 * sg
    taus = [10, 15, 20, 30, 40, 60, 90, 120, 180]
    b = C.binned_beta(t[ok][k], r[k], lc["cad"], taus)
    # adopted: median over 15-60 min (S32, S98), 60-120 min (S6: 30-min cadence cannot be binned below 60 min)
    rng_t = (55, 125) if s == 6 else (14, 61)
    bad = [v for kk, v in b.items() if rng_t[0] <= kk <= rng_t[1]]
    betas[s] = dict(curve=b, adopted=float(max(1.0, np.median(bad))), range_min=rng_t, rms_ppm=float(np.std(r[k]) * 1e6))
df["beta"] = df.sector.map({s: betas[s]["adopted"] for s in betas})
df["e_beta"] = df.e_white * df.beta
out = f"{C.OUT}/t3_events_{variant}.csv"
df.to_csv(out, index=False)
json.dump(dict(variant=variant, template=dict(T14_h=T14 * 24, a_rs=tpl.pm.a, inc=tpl.pm.inc, **C.SHAPE, **tplkw),
               cfg=cfg, ephem_window=dict(P=P, T0=T0), betas=betas, provenance=prov),
          open(f"{C.OUT}/t3_meta_{variant}.json", "w"), indent=1, default=str)
pd.set_option("display.width", 250)
print(df[["sector", "epoch", "coverage", "n", "dt", "dt_alt", "D", "snr", "e_white", "e_pb", "e_boot", "e_beta", "edge", "sel"]].round(3).to_string())
print({s: (betas[s]["adopted"], betas[s]["curve"]) for s in betas})
