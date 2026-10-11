"""Verifier check: recompute RUWE and Gaia-RV-gof percentiles from the reconciler's fresh TAP sample
(../data/ruwe_sample_fresh.csv, read-only), including the 'apparently single' subset for the RV statistics."""
import json, numpy as np, pandas as pd
from pathlib import Path
H = Path(__file__).resolve().parent
f = pd.read_csv(H.parent / "data/ruwe_sample_fresh.csv")
TID, TR, TG, TP = 3220388198192519424, 1.2142985, 2.998642, 4.0627176e-05
tgt = f[f.source_id == TID]
f = f[f.source_id != TID]
single = f[(f.non_single_star == 0) & (f.ipd_frac_multi_peak == 0) & (~f.duplicated_source.astype(bool))]
out = {"target_row": tgt.to_dict("records")}
for lab, s in [("all", f), ("single", single), ("single_ruwe_lt1p4", single[single.ruwe < 1.4]), ("all_ruwe_lt1p4", f[f.ruwe < 1.4])]:
    r = s.ruwe.values
    rv = s[np.isfinite(s.rv_renormalised_gof) & (s.rv_nb_transits >= 10)]
    out[lab] = dict(n=len(s), ruwe_pct_below=float(np.mean(r < TR) * 100), ruwe_median=float(np.median(r)),
                    n_rv=len(rv), gof_pct_below=float(np.mean(rv.rv_renormalised_gof < TG) * 100),
                    frac_p_le_target=float(np.mean(rv.rv_chisq_pvalue <= TP)),
                    frac_variable_rule=float(np.mean((rv.rv_chisq_pvalue < 0.01) & (rv.rv_renormalised_gof > 4))),
                    gof_median=float(np.median(rv.rv_renormalised_gof)))
json.dump(out, open(H / "v01_ruwe_rv_check.json", "w"), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
