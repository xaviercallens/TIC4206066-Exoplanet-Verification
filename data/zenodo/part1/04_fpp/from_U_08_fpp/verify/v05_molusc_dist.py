"""Verifier: companion-mass distributions entering the STP scenario in (i) TRICERATOPS default mode,
(ii) MOLUSC mode (survivors with a(1-e)>10 AU, clamped at 0.1 Msun, padded with q=0 to N), and
(iii) a counterfactual MOLUSC mode without the periastron cut (all survivors).
MOLUSC survivor file = the auditor's fresh MOLUSC b26efbf output (REUSED, sha256 recorded; not re-run here).
trice_env."""
import hashlib, json, os
import numpy as np, pandas as pd
from triceratops.priors import lnprior_bound_TP, sample_q_companion
H = os.path.dirname(os.path.abspath(__file__))
M, plx, N = 0.606, 28.075, 1_000_000
F = "/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/a08_cache/molusc_run/fresh5_kept.csv"
out = {"molusc_file": F, "molusc_sha256": hashlib.sha256(open(F, "rb").read()).hexdigest()}
f_comp = float(np.exp(lnprior_bound_TP(M, plx, np.array([2.2]))[0]))
out["f_comp_default_2.2arcsec"] = f_comp
rs = np.random.RandomState(7); m_def = sample_q_companion(rs.rand(4_000_000), M) * M
d = pd.read_csv(F)
aa, ee = d["semi-major axis(AU)"].to_numpy(), d["eccentricity"].to_numpy()
q = d["mass ratio"].to_numpy().copy()
m_all = np.maximum(q, 0.1 / M) * M
peri = aa * (1 - ee) > 10
out["n_survivors"] = int(len(d)); out["n_peri_gt_10AU"] = int(peri.sum())
out["survivor_period_range_d"] = [float(d["period(days)"].min()), float(d["period(days)"].max())]
out["survivor_q_range"] = [float(q.min()), float(q.max())]
edges = [0.0, 0.1001, 0.2, 0.3, 0.38, 0.45, 0.54, 0.607]
def hist(m, norm):
    return {f"{edges[i]:.2f}-{edges[i+1]:.2f}": float(((m >= edges[i]) & (m < edges[i + 1])).sum() / norm) for i in range(len(edges) - 1)}
# effective prior mass per companion-mass bin (probability per star)
out["eff_prior_default"] = {k: v * f_comp for k, v in hist(m_def, len(m_def)).items()}
out["eff_prior_molusc_kept"] = hist(m_all[peri], N)
out["eff_prior_molusc_all_survivors"] = hist(m_all, N)
for k in ("eff_prior_default", "eff_prior_molusc_kept", "eff_prior_molusc_all_survivors"):
    v = out[k]; out[k + "_total"] = float(sum(v.values()))
    out[k + "_gt0.38"] = float(sum(x for kk, x in v.items() if float(kk.split("-")[0]) >= 0.38))
    out[k + "_gt0.54"] = float(v["0.54-0.61"])
# where do the massive survivors sit? (orbits of survivors with M>0.38)
big = m_all > 0.38
out["survivors_M_gt_0.38"] = dict(n=int(big.sum()), frac_peri_gt_10=float(peri[big].mean()),
                                  a_quantiles_AU=[float(x) for x in np.quantile(aa[big], [0.05, 0.25, 0.5, 0.75, 0.95])],
                                  frac_a_lt_1AU=float((aa[big] < 1).mean()), frac_a_lt_10AU=float((aa[big] < 10).mean()),
                                  frac_P_lt_1yr=float((d["period(days)"].to_numpy()[big] < 365.25).mean()))
json.dump(out, open(os.path.join(H, "v05_molusc_dist.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
