"""R5 (b): constrained FPP where the overluminosity mass limit is applied only to Gaia-blended companions, and the Gaia
source-detection (contrast) and RUWE limits are applied at every separation through MOLUSC's own forward model.

Survival function S(m, a) = fraction of MOLUSC-generated companions (fresh run r01, 1e6, --all, age 5 Gyr, RUWE + Gaia
contrast, b26efbf) in a (mass, semi-major-axis) cell that
    (i) are not rejected by MOLUSC (RUWE test, Gaia contrast test), and
   (ii) are not [projected separation < rho_b  AND  m > M_lim]   (overluminosity, blended companions only).
This S(m, a) is applied to TRICERATOPS 1.1.0's DEFAULT companion prior (Moe & Di Stefano mass ratios, f_comp(<a/d) with
P > 2500 d, a <= 2.2" x d), i.e. no 10 AU periastron cut and no 100%-companion assumption (the two artefacts of
TRICERATOPS' MOLUSC mode).  Likelihoods: verifier's fixed-mass grids L_STP(m), L_PTP(m) (TRICERATOPS' own lnZ at fixed m).
Run with trice_env.
"""
import json, os
import numpy as np, pandas as pd
from triceratops.priors import lnprior_bound_TP, sample_q_companion

H = os.path.dirname(os.path.abspath(__file__))
V = os.path.abspath(os.path.join(H, "..", "..", "08_fpp", "verify"))
ALL = "/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/r5_cache/r5_all.csv"
M, plx = 0.606, 28.075
d_pc = 1000 / plx

cols = ["mass ratio", "semi-major axis(AU)", "Projected Separation(AU)", "DeltaG", "RUWE Rejected", "Gaia Rejected", "Full Rejected"]
df = pd.read_csv(ALL, usecols=cols)
m = df["mass ratio"].to_numpy() * M
a = df["semi-major axis(AU)"].to_numpy()
rho_proj = df["Projected Separation(AU)"].to_numpy() / d_pc          # arcsec
rej_all = df["Full Rejected"].to_numpy().astype(bool)
rej_ruwe = df["RUWE Rejected"].to_numpy().astype(bool)
rej_gaia = df["Gaia Rejected"].to_numpy().astype(bool)
n_gen = len(df)

# TRICERATOPS prior over a (a = rho * d convention), restricted to 2.2"
a_max = 2.2 * d_pc
a_edges = np.geomspace(3.0, a_max, 31)
Fa = np.exp(lnprior_bound_TP(M, plx, a_edges / d_pc)); Fa[~np.isfinite(Fa)] = 0.0
dF = np.diff(Fa)                                                       # prior mass per a bin
f22 = float(Fa[-1])
m_edges = np.round(np.arange(0.10, 0.6061, 0.02), 3); m_edges[-1] = 0.6061
m_cent = 0.5 * (m_edges[1:] + m_edges[:-1])
sel = (m >= 0.1) & (a >= a_edges[0]) & (a < a_edges[-1])
im = np.clip(np.digitize(m[sel], m_edges) - 1, 0, len(m_cent) - 1)
ia = np.clip(np.digitize(a[sel], a_edges) - 1, 0, len(dF) - 1)
ncell = np.zeros((len(m_cent), len(dF))); np.add.at(ncell, (im, ia), 1)

def S_table(M_lim, rho_b, use_molusc=True, only=None):
    rej = np.zeros(n_gen, bool)
    if use_molusc:
        rej |= {"ruwe": rej_ruwe, "gaia": rej_gaia}.get(only, rej_all)
    rej |= (rho_proj < rho_b) & (m > M_lim)
    s = np.zeros_like(ncell); np.add.at(s, (im, ia), (~rej[sel]).astype(float))
    with np.errstate(invalid="ignore", divide="ignore"):
        S = np.where(ncell > 0, s / ncell, 1.0)
    return S

def w_of_m(S):
    """prior weight (absolute) per mass bin centre: sum_a dF(a) S(m, a)"""
    return S @ dF

out = dict(molusc_all_file=ALL, n_generated=int(n_gen), n_survivors=int((~rej_all).sum()), f_comp_2p2=f22,
           cells_min_count=int(ncell[:, dF > 0].min()), cells_median_count=float(np.median(ncell)),
           prior_share_by_a={f"{a_edges[i]:.1f}-{a_edges[i+1]:.1f}AU": float(dF[i] / f22) for i in range(len(dF))})
# diagnostic: MOLUSC survival (RUWE+Gaia) of M>0.38 companions vs projected separation (TRICERATOPS a range)
diag = {}
for lo, hi in ((0.1, 0.2), (0.2, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 1.0), (1.0, 1.5), (1.5, 2.2)):
    k = (m > 0.38) & (rho_proj >= lo) & (rho_proj < hi) & (a >= 3.0) & (a < a_max)
    diag[f"{lo}-{hi}\""] = dict(n=int(k.sum()), surv_all=float((~rej_all[k]).mean()), surv_ruwe=float((~rej_ruwe[k]).mean()),
                                surv_gaia=float((~rej_gaia[k]).mean()))
out["survival_M_gt_0.38_by_projected_sep"] = diag

rs = np.random.RandomState(11); m_def = sample_q_companion(rs.rand(1_000_000), M) * M
imd = np.clip(np.digitize(m_def, m_edges) - 1, 0, len(m_cent) - 1)
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
    LS, LP = EL(m_def, lS), EL(m_def, lP)
    def fpp_w(wm):
        w = wm[imd]
        return fpp(np.log((LS * w).mean()) + ref, np.log((LP * w).mean()) + ref)
    res = dict(grid=gname, FPP_runs_mean=float(np.mean([r["FPP"] for r in runs])),
               FPP_default_binned_prior=fpp_w(np.full(len(m_cent), f22)))
    tab = {}
    tab["molusc_only(RUWE+Gaia)"] = fpp_w(w_of_m(S_table(9, 0.0)))
    tab["molusc_RUWE_only"] = fpp_w(w_of_m(S_table(9, 0.0, only="ruwe")))
    tab["molusc_Gaia_contrast_only"] = fpp_w(w_of_m(S_table(9, 0.0, only="gaia")))
    for M_lim in (0.38, 0.45, 0.54):
        for rho_b in (0.2, 0.3, 0.5, 1.0, 2.2):
            tab[f"Mlim{M_lim:.2f}_blended<{rho_b:.1f}\"_+molusc"] = fpp_w(w_of_m(S_table(M_lim, rho_b)))
            tab[f"Mlim{M_lim:.2f}_blended<{rho_b:.1f}\"_no_molusc"] = fpp_w(w_of_m(S_table(M_lim, rho_b, use_molusc=False)))
    res["FPP"] = tab
    # STP evidence share of M>0.38 after the pessimistic constraint (blended<0.2" + MOLUSC)
    wm = w_of_m(S_table(0.38, 0.2)); w = wm[imd]
    res["STP_share_M_gt_0.38_after_blended0.2_molusc"] = float((LS * w)[m_def > 0.38].sum() / (LS * w).sum())
    out[sig] = res
json.dump(out, open(os.path.join(H, "out_r03_constrained_molusc.json"), "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k not in ("A", "B", "prior_share_by_a")}, indent=1))
for s in ("A", "B"):
    r = out[s]; print(s, "runs=%.4f default(binned prior)=%.4f" % (r["FPP_runs_mean"], r["FPP_default_binned_prior"]),
                      "STP share M>0.38 after pessimistic constraint: %.2f" % r["STP_share_M_gt_0.38_after_blended0.2_molusc"])
    for k, v in r["FPP"].items():
        print("   %-44s %.4f" % (k, v))
