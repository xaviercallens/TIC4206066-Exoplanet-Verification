"""Audit 8: compare the effective bound-companion prior in TRICERATOPS default mode vs MOLUSC mode (trice_env)."""
import json, os, sys
import numpy as np, pandas as pd
from triceratops.priors import lnprior_bound_TP, sample_q_companion
H = os.path.dirname(os.path.abspath(__file__))
M, plx, N = 0.606, 28.075, 1_000_000
out = {}
f_def = float(np.exp(lnprior_bound_TP(M, plx, np.array([2.2]))[0]))
out["default_f_comp_within_2.2arcsec"] = f_def
out["default_2.2arcsec_in_AU"] = 2.2 * 1000 / plx
np.random.seed(1); q = sample_q_companion(np.random.rand(N), M); m = q * M
out["default_q_prior"] = {"frac_mass_gt_0.38": float((m > 0.38).mean()), "frac_mass_gt_0.54": float((m > 0.54).mean()), "m_min": float(m.min()), "m_max": float(m.max())}
out["default_effective_prior_mass_gt_0.38"] = f_def * float((m > 0.38).mean())
files = {"prev_age5 (tic4206066_PUBLICATION_AUDIT, reused)": "/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/a08_cache/molusc_prev/molusc_age5_kept.csv",
         "fresh_age5 (this audit)": "/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/a08_cache/molusc_run/fresh5_kept.csv"}
for k, f in files.items():
    if not os.path.exists(f): continue
    d = pd.read_csv(f); a, e = d["semi-major axis(AU)"].values, d["eccentricity"].values
    kept = d[a * (1 - e) > 10]; qk = kept["mass ratio"].values.copy(); mk = np.maximum(qk, 0.1 / M) * M
    sep = kept["Projected Separation(AU)"].values
    out[k] = {"n_survivors": int(len(d)), "n_periastron_gt_10AU": int(len(kept)), "effective_prior_at_N1e6": len(kept) / N,
                  "effective_prior_at_N5e5": len(kept) / 5e5,
                  "frac_kept_q_below_0.1Msun_clamped": float((qk < 0.1 / M).mean()),
                  "frac_kept_mass_gt_0.38": float((mk > 0.38).mean()), "frac_kept_mass_gt_0.54": float((mk > 0.54).mean()),
                  "effective_prior_mass_gt_0.38_at_N1e6": float((mk > 0.38).sum() / N),
                  "frac_kept_projsep_gt_2.2arcsec": float((sep > 2.2 * 1000 / plx).mean()),
                  "frac_survivors_mass_gt_0.38": float((np.maximum(d["mass ratio"].values, 0.1 / M) * M > 0.38).mean())}
json.dump(out, open(os.path.join(H, "out_f03_priors.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
