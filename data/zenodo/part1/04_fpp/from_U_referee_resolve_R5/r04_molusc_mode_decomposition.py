"""R5 (a): what actually lowers the FPP in TRICERATOPS' MOLUSC mode?  Factorial decomposition with the fresh MOLUSC
--all output (r01), which contains every generated companion and its Gaia/RUWE rejection flags (the audit and the verifier
only had the survivors file, so they could not switch the Gaia limits off while keeping everything else).

MOLUSC mode as coded in TRICERATOPS 1.1.0 (marginal_likelihoods.py l.946-951, l.1038): take companions from the file,
keep a(1-e) > 10 AU, clamp q >= 0.1/M, pad with q = 0 to N, prior 1 for every draw (Z = sum L / N).
Here the same arithmetic is applied, via the verifier's fixed-mass likelihood grids, to these populations:
   population = all generated (no Gaia limits) | survivors (RUWE + Gaia contrast)
   separation selection = none | a(1-e) > 10 AU (as coded) | TRICERATOPS default window (P > 2500 d with M = 1 Msun,
                          i.e. a >= 3.62 AU, and a <= 2.2" x d)
Run with trice_env.
"""
import json, os
import numpy as np, pandas as pd
from triceratops.priors import lnprior_bound_TP, sample_q_companion

H = os.path.dirname(os.path.abspath(__file__))
V = os.path.abspath(os.path.join(H, "..", "..", "08_fpp", "verify"))
ALL = "/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/r5_cache/r5_all.csv"
M, plx, N = 0.606, 28.075, 1_000_000
d_pc = 1000 / plx
df = pd.read_csv(ALL, usecols=["mass ratio", "period(days)", "semi-major axis(AU)", "eccentricity", "Full Rejected"])
q = df["mass ratio"].to_numpy(); a = df["semi-major axis(AU)"].to_numpy(); e = df["eccentricity"].to_numpy()
Pd = df["period(days)"].to_numpy(); surv = ~df["Full Rejected"].to_numpy().astype(bool)
m = np.maximum(q, 0.1 / M) * M                                  # TRICERATOPS clamp
a_in = (1.0 * (2500 / 365.25) ** 2) ** (1 / 3)                  # TRICERATOPS default inner limit (P > 2500 d, M = 1)
win = (a >= a_in) & (a <= 2.2 * d_pc)
peri = a * (1 - e) > 10
f22 = float(np.exp(lnprior_bound_TP(M, plx, np.array([2.2]))[0]))
rs = np.random.RandomState(11); m_def = sample_q_companion(rs.rand(1_000_000), M) * M
out = dict(n_generated=len(df), n_survivors=int(surv.sum()), a_inner_default_AU=a_in, a_outer_2p2_AU=2.2 * d_pc,
           counts={}, FPP={})
for sig, gname in (("A", "grid_A_v_A_S32+S98_5min_N2e6.json"), ("B", "grid_B_v_B_S32+S98_5min.json")):
    g = json.load(open(os.path.join(V, "runs_v", gname)))
    rows = sorted(g["rows"], key=lambda r: r["m"])
    mg = np.array([r["m"] for r in rows]); lS = np.array([r["lnL_STP"] for r in rows]); lP = np.array([r["lnL_PTP"] for r in rows])
    ref = max(lS[np.isfinite(lS)].max(), lP.max())
    EL = lambda ms, lg: np.exp(np.interp(ms, mg, np.where(np.isfinite(lg), lg, -1e4)) - ref)
    runs = json.load(open(os.path.join(V, "runs_v", f"default_{sig}_v_{sig}_S32+S98_5min.json")))["runs"]
    base = np.log(np.mean(np.exp(np.array([r["lnZ"] for r in runs]) - ref), axis=0)) + ref
    def fpp(zS, zP):
        z = base.copy(); z[6] = zS; z[3] = zP
        fin = np.isfinite(z); p = np.zeros_like(z); p[fin] = np.exp(z[fin] - z[fin].max()); p /= p.sum()
        return float(1 - p[0] - p[3] - p[9])
    LSall, LPall = EL(m, lS), EL(m, lP)
    def mol(sel):   # MOLUSC-mode arithmetic: sum over selected draws / N, prior 1
        return fpp(np.log(LSall[sel].sum() / N) + ref, np.log(LPall[sel].sum() / N) + ref)
    res = {}
    res["default_TRICERATOPS_prior"] = fpp(np.log(f22 * EL(m_def, lS).mean()) + ref, np.log(f22 * EL(m_def, lP).mean()) + ref)
    for pop, P_ in (("all_generated_noGaia", np.ones(len(df), bool)), ("survivors_RUWE+Gaia", surv)):
        for sname, S_ in (("no_sep_cut", np.ones(len(df), bool)), ("peri>10AU(as_coded)", peri), ("default_window_3.6-78AU", win)):
            k = P_ & S_
            res[f"{pop}|{sname}"] = mol(k)
            out["counts"][f"{pop}|{sname}"] = {"n": int(k.sum()), "n_M_gt_0.38": int((k & (m > 0.38)).sum())}
    out["FPP"][sig] = res
json.dump(out, open(os.path.join(H, "out_r04_molusc_decomposition.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
