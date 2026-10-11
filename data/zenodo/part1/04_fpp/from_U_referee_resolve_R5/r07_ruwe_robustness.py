"""R5: how much does the Gaia-constrained FPP depend on MOLUSC's RUWE test?  Re-derive each companion's RUWE rejection
probability with a more lenient 'observed' RUWE (1.4, the usual binarity threshold, and 2.0) instead of the target's 1.214,
using MOLUSC's own rule p = halfnorm.cdf(pred, loc=obs, scale=s).  The scale s is recovered per companion from MOLUSC's
output (p, pred at obs = 1.214).  Gaia-contrast rejections are kept as drawn by MOLUSC.  Then the r03 machinery
(TRICERATOPS default prior x survival, verifier grids) gives the FPP.  trice_env."""
import json, os
import numpy as np, pandas as pd
from scipy.special import erfinv
from triceratops.priors import lnprior_bound_TP, sample_q_companion

H = os.path.dirname(os.path.abspath(__file__))
V = os.path.abspath(os.path.join(H, "..", "..", "08_fpp", "verify"))
ALL = "/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/r5_cache/r5_all.csv"
M, plx = 0.606, 28.075; d_pc = 1000 / plx; OBS = 1.2142984867095947
df = pd.read_csv(ALL, usecols=["mass ratio", "semi-major axis(AU)", "Projected Separation(AU)", "Predicted RUWE", "RUWE Rejection Prob", "Gaia Rejected"])
m = df["mass ratio"].to_numpy() * M; a = df["semi-major axis(AU)"].to_numpy(); rho = df["Projected Separation(AU)"].to_numpy() / d_pc
pred = df["Predicted RUWE"].to_numpy(); p0 = df["RUWE Rejection Prob"].to_numpy(); rej_g = df["Gaia Rejected"].to_numpy().astype(bool)
ok = (p0 > 0) & (p0 < 1) & (pred > OBS)
s = np.full(len(df), np.nan); s[ok] = (pred[ok] - OBS) / (np.sqrt(2) * erfinv(p0[ok]))
s_med = np.nanmedian(s); s[~np.isfinite(s)] = s_med
applies = ~((p0 == 0) & (pred > OBS))           # MOLUSC sets p=0 outside its calibration box; keep that
rng = np.random.RandomState(3)
def p_rej(obs):
    from scipy.stats import halfnorm
    p = np.where(pred > obs, halfnorm.cdf(pred, loc=obs, scale=s), 0.0)
    return np.where(applies, p, 0.0)

a_edges = np.geomspace(3.0, 2.2 * d_pc, 31)
Fa = np.exp(lnprior_bound_TP(M, plx, a_edges / d_pc)); Fa[~np.isfinite(Fa)] = 0.0; dF = np.diff(Fa); f22 = Fa[-1]
m_edges = np.round(np.arange(0.10, 0.6061, 0.02), 3); m_edges[-1] = 0.6061; m_c = 0.5 * (m_edges[1:] + m_edges[:-1])
sel = (m >= 0.1) & (a >= a_edges[0]) & (a < a_edges[-1])
im = np.clip(np.digitize(m[sel], m_edges) - 1, 0, len(m_c) - 1); ia = np.clip(np.digitize(a[sel], a_edges) - 1, 0, len(dF) - 1)
n = np.zeros((len(m_c), len(dF))); np.add.at(n, (im, ia), 1)
def wm(surv_prob):
    sm = np.zeros_like(n); np.add.at(sm, (im, ia), surv_prob[sel]); return (np.where(n > 0, sm / np.maximum(n, 1), 1.0) @ dF)
rs = np.random.RandomState(11); m_def = sample_q_companion(rs.rand(1_000_000), M) * M
imd = np.clip(np.digitize(m_def, m_edges) - 1, 0, len(m_c) - 1)
out = {"recovered_scale_median": float(s_med), "check_p_recomputed_vs_molusc_at_obs": None, "FPP": {}}
pr_chk = p_rej(OBS); out["check_p_recomputed_vs_molusc_at_obs"] = float(np.nanmax(np.abs(pr_chk[ok] - p0[ok])))
for sig, gname in (("A", "grid_A_v_A_S32+S98_5min_N2e6.json"), ("B", "grid_B_v_B_S32+S98_5min.json")):
    g = json.load(open(os.path.join(V, "runs_v", gname))); rows = sorted(g["rows"], key=lambda r: r["m"])
    mg = np.array([r["m"] for r in rows]); lS = np.array([r["lnL_STP"] for r in rows]); lP = np.array([r["lnL_PTP"] for r in rows])
    ref = max(lS[np.isfinite(lS)].max(), lP.max())
    EL = lambda ms, lg: np.exp(np.interp(ms, mg, np.where(np.isfinite(lg), lg, -1e4)) - ref)
    runs = json.load(open(os.path.join(V, "runs_v", f"default_{sig}_v_{sig}_S32+S98_5min.json")))["runs"]
    base = np.log(np.mean(np.exp(np.array([r["lnZ"] for r in runs]) - ref), axis=0)) + ref
    def fpp(zS, zP):
        z = base.copy(); z[6] = zS; z[3] = zP; fin = np.isfinite(z); p = np.zeros_like(z); p[fin] = np.exp(z[fin] - z[fin].max()); p /= p.sum()
        return float(1 - p[0] - p[3] - p[9])
    LS, LP = EL(m_def, lS), EL(m_def, lP)
    f = lambda w: fpp(np.log((LS * w[imd]).mean()) + ref, np.log((LP * w[imd]).mean()) + ref)
    res = {"default": f(np.full(len(m_c), f22))}
    for obs in (OBS, 1.4, 2.0, 3.0):
        surv_ruwe = 1 - p_rej(obs)
        res[f"RUWE_obs{obs:.2f}+contrast"] = f(wm(surv_ruwe * (~rej_g)))
        res[f"RUWE_obs{obs:.2f}+contrast+overlum0.38<0.2\""] = f(wm(surv_ruwe * (~rej_g) * ~((rho < 0.2) & (m > 0.38))))
    res["contrast_only"] = f(wm((~rej_g).astype(float)))
    out["FPP"][sig] = res
json.dump(out, open(os.path.join(H, "out_r07_ruwe_robustness.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
