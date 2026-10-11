"""R5 (b): FPP with the companion-mass limit applied ONLY to companions that Gaia cannot resolve (referee R4(b)),
evaluated on the verifier's fixed-companion-mass likelihood grids (08_fpp/verify/runs_v/grid_*.json; TRICERATOPS 1.1.0
lnZ_STP / lnZ_PTP at fixed companion mass, own fold 5-min bins, Teff 4012) and the verifier's 8 default runs for the
other scenarios.  Run with trice_env (uses TRICERATOPS' own priors).

Companion population = TRICERATOPS 1.1.0 default:
  mass      m = q M_s, q ~ sample_q_companion (independent of separation)
  separation: cumulative prior f_comp(<rho) = exp(lnprior_bound_TP(M_s, plx, rho))  (Moe & Di Stefano 2017; only
              P > 2500 d counted; rho is converted to a = rho * d, i.e. TRICERATOPS' own convention); default rho_max = 2.2"
Constraints (hard limits, TRICERATOPS contrast-curve style, no orbital projection):
  overluminosity (Way et al. 2026 XP; limit M_lim): excludes m > M_lim only for rho < rho_b ("Gaia-blended")
  Gaia source detection (no Gaia DR3 source within 11"): MOLUSC gaia_contrast.txt (Wood et al. 2021), used either as a
      hard limit at a given recovery level or as a detection probability (MOLUSC's step rule); DeltaG(m) from MOLUSC/BHAC15
      5 Gyr (tabulated from the MOLUSC survivors file).
  Allowed weight for a companion of mass m: w(m) = E_rho[ 1{rho >= rho_lo(m)} * (1 - P_det(DeltaG(m), rho)) ] * f_comp(2.2")
FPP from run-averaged default evidences with lnZ_STP and lnZ_PTP replaced (same as verifier v06_combine.py).
"""
import json, os, sys
import numpy as np, pandas as pd
from triceratops.priors import lnprior_bound_TP, sample_q_companion

H = os.path.dirname(os.path.abspath(__file__))
V = os.path.abspath(os.path.join(H, "..", "..", "08_fpp", "verify"))
SCR = "/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad"
M, plx = 0.606, 28.075
d_pc = 1000 / plx

# ---- separation prior (TRICERATOPS) ----
rho = np.concatenate([np.linspace(0.001, 2.2, 22000)])
F = np.exp(lnprior_bound_TP(M, plx, rho)); F[~np.isfinite(F)] = 0.0
f22 = float(np.exp(lnprior_bound_TP(M, plx, np.array([2.2]))[0]))
dF = np.diff(np.concatenate([[0.0], F]))            # prior mass in each rho cell (upper edge = rho)
rho_c = rho                                          # use upper edge (conservative for the inner cut: cell counted as outside if rho>=rho_b)
sep_share = {f"<{r:.2f}\"": float(np.exp(lnprior_bound_TP(M, plx, np.array([r]))[0]) / f22) if r > 0.1 else 0.0 for r in (0.2, 0.3, 0.5, 0.7, 1.0, 1.5)}

# ---- DeltaG(m) from MOLUSC (BHAC15, 5 Gyr) ----
kept = pd.read_csv(os.path.join(SCR, "a08_cache", "molusc_run", "fresh5_kept.csv"), nrows=300000, usecols=["mass ratio", "DeltaG"])
mm = kept["mass ratio"].to_numpy() * M; dg = kept["DeltaG"].to_numpy(); ok = np.isfinite(dg)
o = np.argsort(mm[ok]); m_tab, dg_tab = mm[ok][o], dg[ok][o]
def DeltaG(m):
    return np.interp(m, m_tab, dg_tab)

# ---- Gaia recovery table (MOLUSC gaia_contrast.txt, Sep in mas) ----
gc = pd.read_csv(os.path.join(SCR, "a08_cache", "MOLUSC", "code", "gaia_contrast.txt"), sep=r"\s+")
g_sep = gc["Sep"].to_numpy() / 1000.0
rates = [0.2, 0.5, 0.7, 0.9]; cols = ["20%", "50%", "70%", "90%"]
def contrast_at(r, col):
    return np.interp(r, g_sep, gc[col].to_numpy(), left=0.0)
def p_det(dG, r, mode):
    """mode 'prob' = MOLUSC step rule (largest rate whose contrast limit >= dG); 'hardXX' = hard limit at rate XX."""
    if mode == "none":
        return np.zeros_like(r)
    if mode.startswith("hard"):
        col = mode[4:] + "%"
        return (dG <= contrast_at(r, col)).astype(float)
    p = np.zeros_like(r)
    for rate, col in zip(rates, cols):
        p = np.where(dG <= contrast_at(r, col), rate, p)
    return p

def weights(m_draws, M_lim, rho_b, gaia_mode):
    """Allowed prior weight per mass draw (absolute, includes f_comp)."""
    mu = np.unique(np.round(m_draws, 3))
    wtab = {}
    for m in mu:
        dG = DeltaG(m)
        surv = 1.0 - p_det(dG, rho_c, gaia_mode)
        if m > M_lim:
            surv = surv * (rho_c >= rho_b)
        wtab[m] = float(np.sum(dF * surv))
    return np.array([wtab[x] for x in np.round(m_draws, 3)])

rs = np.random.RandomState(11); m_def = sample_q_companion(rs.rand(1_000_000), M) * M
out = dict(f_comp_2p2=f22, sep_share_of_prior=sep_share, DeltaG_at={f"{x:.2f}": float(DeltaG(x)) for x in (0.2, 0.3, 0.38, 0.45, 0.54, 0.6)},
           note="rho_b = separation below which a companion is Gaia-blended (overluminosity applies); gaia_mode = treatment of Gaia source detection outside")
for sig, gname in (("A", "grid_A_v_A_S32+S98_5min_N2e6.json"), ("A_N4e5", "grid_A_v_A_S32+S98_5min.json"), ("B", "grid_B_v_B_S32+S98_5min.json")):
    s0 = sig[0]
    g = json.load(open(os.path.join(V, "runs_v", gname)))
    rows = sorted(g["rows"], key=lambda r: r["m"])
    mg = np.array([r["m"] for r in rows]); lS = np.array([r["lnL_STP"] for r in rows]); lP = np.array([r["lnL_PTP"] for r in rows])
    ref = max(lS[np.isfinite(lS)].max(), lP.max())
    def EL(ms, lg):
        lg2 = np.where(np.isfinite(lg), lg, -1e4); return np.exp(np.interp(ms, mg, lg2) - ref)
    runs = json.load(open(os.path.join(V, "runs_v", f"default_{s0}_v_{s0}_S32+S98_5min.json")))["runs"]
    base = np.log(np.mean(np.exp(np.array([r["lnZ"] for r in runs]) - ref), axis=0)) + ref
    def fpp(zS, zP):
        z = base.copy(); z[6] = zS; z[3] = zP
        fin = np.isfinite(z); p = np.zeros_like(z); p[fin] = np.exp(z[fin] - z[fin].max()); p /= p.sum()
        return float(1 - p[0] - p[3] - p[9])
    LS, LP = EL(m_def, lS), EL(m_def, lP)
    res = {"grid_file": gname, "grid_N": g["N"], "FPP_runs_mean": float(np.mean([r["FPP"] for r in runs]))}
    res["FPP_default_grid"] = fpp(np.log(f22 * LS.mean()) + ref, np.log(f22 * LP.mean()) + ref)
    tab = {}
    for M_lim in (0.38, 0.45, 0.54):
        for rho_b in (0.2, 0.3, 0.5, 0.7, 1.0, 2.2):
            for gm in ("none", "prob", "hard50", "hard20"):
                w = weights(m_def, M_lim, rho_b, gm)
                tab[f"Mlim{M_lim:.2f}_rhob{rho_b:.1f}_gaia-{gm}"] = fpp(np.log((LS * w).mean()) + ref, np.log((LP * w).mean()) + ref)
    for gm in ("prob", "hard50"):
        w = weights(m_def, 9.0, 0.0, gm)
        tab[f"no_overlum_gaia-{gm}"] = fpp(np.log((LS * w).mean()) + ref, np.log((LP * w).mean()) + ref)
    res["FPP"] = tab
    # STP evidence share by separation for M>0.38 (default prior): which rho the massive-companion evidence sits at
    out[sig] = res
json.dump(out, open(os.path.join(H, "out_r02_constrained_grid.json"), "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k not in ("A", "A_N4e5", "B")}, indent=1))
for s in ("A", "A_N4e5", "B"):
    r = out[s]; print(s, "default(grid)=%.4f runs=%.4f" % (r["FPP_default_grid"], r["FPP_runs_mean"]))
    for k, v in r["FPP"].items():
        print("   %-40s %.4f" % (k, v))
