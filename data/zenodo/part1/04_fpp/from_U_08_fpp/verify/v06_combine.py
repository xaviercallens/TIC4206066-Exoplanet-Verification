"""Verifier: combine the fixed-mass likelihood grids (v04) with companion-mass priors (v05 logic) and the other
scenarios' lnZ from own default runs (v03) to decompose the STP evidence and to evaluate counterfactual priors:
  default (TRICERATOPS Moe&DiStefano q prior x f_comp(<2.2"))
  default truncated at M<=0.38 / 0.45
  MOLUSC mode as coded (survivors with a(1-e)>10 AU, clamp 0.1 Msun, pad to N=1e6, prior 1)
  MOLUSC 'Gaia only' (all survivors, no periastron cut)            -> isolates the effect of the 10 AU cut
  MOLUSC 'kept' but massive (>0.38) removed                          -> overluminosity on top of MOLUSC
trice_env."""
import json, os, sys
import numpy as np, pandas as pd
from triceratops.priors import lnprior_bound_TP, sample_q_companion
H = os.path.dirname(os.path.abspath(__file__))
M, plx, N = 0.606, 28.075, 1_000_000
F = "/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/a08_cache/molusc_run/fresh5_kept.csv"
f_comp = float(np.exp(lnprior_bound_TP(M, plx, np.array([2.2]))[0]))
rs = np.random.RandomState(7); m_def = sample_q_companion(rs.rand(2_000_000), M) * M
d = pd.read_csv(F)
aa, ee = d["semi-major axis(AU)"].to_numpy(), d["eccentricity"].to_numpy()
m_mol = np.maximum(d["mass ratio"].to_numpy(), 0.1 / M) * M
peri = aa * (1 - ee) > 10
out = {}
TAG = sys.argv[1] if len(sys.argv) > 1 else ""
SIGS = sys.argv[2].split(",") if len(sys.argv) > 2 else ["A", "B"]
for sig in SIGS:
    g = json.load(open(os.path.join(H, "runs_v", f"grid_{sig}_v_{sig}_S32+S98_5min{TAG}.json")))
    rows = sorted(g["rows"], key=lambda r: r["m"])
    mg = np.array([r["m"] for r in rows])
    lS = np.array([r["lnL_STP"] for r in rows]); lP = np.array([r["lnL_PTP"] for r in rows])
    ref = max(lS[np.isfinite(lS)].max(), lP.max())
    def EL(ms, lg, weights=None):
        lg2 = np.where(np.isfinite(lg), lg, -1e4)
        return np.exp(np.interp(ms, mg, lg2) - ref)
    runs = json.load(open(os.path.join(H, "runs_v", f"default_{sig}_v_{sig}_S32+S98_5min.json")))["runs"]
    lnZ_runs = np.array([r["lnZ"] for r in runs])
    base = np.log(np.mean(np.exp(lnZ_runs - ref), axis=0)) + ref     # run-averaged evidences
    def fpp(lnZ_STP, lnZ_PTP):
        z = base.copy(); z[6] = lnZ_STP; z[3] = lnZ_PTP
        fin = np.isfinite(z); p = np.zeros_like(z); p[fin] = np.exp(z[fin] - z[fin].max()); p /= p.sum()
        return float(1 - p[0] - p[3] - p[9]), float(p[6])
    res = {"grid_N": g["N"], "grid": {f"{m:.3f}": [float(a), float(b)] for m, a, b in zip(mg, lS, lP)}}
    # validation: grid-predicted default evidences vs direct runs
    zS = np.log(f_comp * EL(m_def, lS).mean()) + ref; zP = np.log(f_comp * EL(m_def, lP).mean()) + ref
    res["check_lnZ_STP_grid_vs_runs"] = [float(zS), float(base[6])]
    res["check_lnZ_PTP_grid_vs_runs"] = [float(zP), float(base[3])]
    res["FPP_runs_mean"] = float(np.mean([r["FPP"] for r in runs]))
    res["FPP_default_from_grid"] = fpp(zS, zP)[0]
    wS = EL(m_def, lS)
    edges = [0.1, 0.2, 0.3, 0.38, 0.45, 0.54, 0.607]
    res["STP_evidence_share_by_mass_default"] = {f"{edges[i]:.2f}-{edges[i+1]:.2f}": float(wS[(m_def >= edges[i]) & (m_def < edges[i + 1])].sum() / wS.sum()) for i in range(len(edges) - 1)}
    res["STP_evidence_share_gt0.38_default"] = float(wS[m_def > 0.38].sum() / wS.sum())
    for mx in (0.30, 0.38, 0.45, 0.54):
        k = m_def <= mx
        res[f"FPP_default_trunc{mx:.2f}"] = fpp(np.log(f_comp * (EL(m_def, lS) * k).mean()) + ref,
                                              np.log(f_comp * (EL(m_def, lP) * k).mean()) + ref)[0]
    def mol(sel):
        return fpp(np.log(EL(m_mol[sel], lS).sum() / N) + ref, np.log(EL(m_mol[sel], lP).sum() / N) + ref)[0]
    res["FPP_molusc_as_coded"] = mol(peri)
    res["FPP_molusc_gaia_only_no_peri_cut"] = mol(np.ones(len(m_mol), bool))
    res["FPP_molusc_kept_and_M_le_0.38"] = mol(peri & (m_mol <= 0.38))
    res["FPP_molusc_gaia_only_and_M_le_0.38"] = mol(m_mol <= 0.38)
    wk = EL(m_mol[peri], lS); res["STP_evidence_share_gt0.38_molusc_kept"] = float(wk[m_mol[peri] > 0.38].sum() / wk.sum())
    wa = EL(m_mol, lS); res["STP_evidence_share_gt0.38_molusc_all"] = float(wa[m_mol > 0.38].sum() / wa.sum())
    # relative STP evidence vs default
    res["Z_STP_ratio_molusc_coded_over_default"] = float(EL(m_mol[peri], lS).sum() / N / (f_comp * EL(m_def, lS).mean()))
    res["Z_STP_ratio_molusc_gaia_only_over_default"] = float(EL(m_mol, lS).sum() / N / (f_comp * EL(m_def, lS).mean()))
    out[sig] = res
json.dump(out, open(os.path.join(H, f"v06_combine{TAG}.json"), "w"), indent=1)
for s, r in out.items():
    print(s, {k: v for k, v in r.items() if k != "grid"})
