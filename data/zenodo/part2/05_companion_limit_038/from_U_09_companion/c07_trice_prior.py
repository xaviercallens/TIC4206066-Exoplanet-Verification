"""Audit 9, step 9: evaluate TRICERATOPS 1.1.0 prior functions (imported from the installed package, not modified)
for this target: bound-companion prior within 2.2", companion-mass distribution, eccentricity prior."""
import json, numpy as np
from pathlib import Path
from triceratops.priors import lnprior_bound_TP, sample_q_companion, sample_ecc
from triceratops.funcs import flux_relation, stellar_relations
M_s, plx = 0.606, 28.075
out = {"lnprior_bound_TP_2.2as": float(lnprior_bound_TP(M_s, plx, np.array([2.2]))[0])}
out["f_bound_2.2as"] = float(np.exp(out["lnprior_bound_TP_2.2as"]))
np.random.seed(3); q = sample_q_companion(np.random.rand(1_000_000), M_s); Mc = q * M_s
out["Mc_range"] = [float(Mc.min()), float(Mc.max())]
out["frac_Mc_gt_0.38"] = float(np.mean(Mc > 0.38)); out["frac_Mc_gt_0.52"] = float(np.mean(Mc > 0.5232))
out["frac_Mc_lt_0.168"] = float(np.mean(Mc < 0.1677))
fr = flux_relation(Mc) / (flux_relation(Mc) + flux_relation(np.array([M_s])))
out["trice_fluxratio_at"] = {str(m): float(flux_relation(np.array([m]))[0] / (flux_relation(np.array([m]))[0] + flux_relation(np.array([M_s]))[0])) for m in [0.1, 0.2, 0.38, 0.5]}
R, T = stellar_relations(np.array([0.1, 0.2, 0.38, 0.5, 0.606]), np.full(5, 0.623), np.full(5, 4030.))
out["trice_radius_at_[0.1,0.2,0.38,0.5,0.606]"] = R.tolist()
e = sample_ecc(np.random.rand(1_000_000), planet=True, P_orb=3.18)
out["ecc_prior_mean"] = float(e.mean()); out["ecc_prior_frac_gt_0.2"] = float(np.mean(e > 0.2)); out["ecc_prior_frac_gt_0.46"] = float(np.mean(e > 0.46))
print(json.dumps(out, indent=1)); json.dump(out, open(Path(__file__).resolve().parent / "out/triceratops_priors.json", "w"), indent=1)
