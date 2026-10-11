"""Is the Gaia DR3 RV p-value of TIC 4206066 (rv_chisq_pvalue 4.06e-5, gof 3.0, 8 of 33 transits deblended) unusual among
comparable stars, and does deblending raise RV scatter? Same selection as T6 (G 10-11, BP-RP 1.6-1.9, parallax > 15), plus RV
columns. Read-only Gaia archive TAP query."""
import json, numpy as np
from astroquery.gaia import Gaia
q = ("SELECT source_id, phot_g_mean_mag, bp_rp, parallax, ruwe, non_single_star, ipd_frac_multi_peak, duplicated_source, "
     "rv_renormalised_gof, rv_chisq_pvalue, rv_nb_transits, rv_nb_deblended_transits, rv_amplitude_robust, rv_template_teff, "
     "radial_velocity_error, grvs_mag FROM gaiadr3.gaia_source WHERE phot_g_mean_mag BETWEEN 10.0 AND 11.0 "
     "AND bp_rp BETWEEN 1.60 AND 1.90 AND parallax > 15 AND rv_nb_transits >= 10")
t = Gaia.launch_job(q).get_results().to_pandas(); t.to_csv("rv01_sample.csv", index=False)
TID = 3220388198192519424
s = t[(t.source_id != TID) & (t.non_single_star == 0) & (t.ipd_frac_multi_peak == 0) & (~t.duplicated_source.astype(bool)) & (t.ruwe < 1.4)].copy()
s["fdeb"] = s.rv_nb_deblended_transits / s.rv_nb_transits
out = dict(n=len(s), target=t[t.source_id == TID].to_dict("records"))
for lab, m in [("no_deblend", s.fdeb == 0), ("some_deblend", s.fdeb > 0), ("deblend_ge_15pct", s.fdeb >= 0.15)]:
    x = s[m]; out[lab] = dict(n=int(len(x)), frac_p_le_target=float(np.mean(x.rv_chisq_pvalue <= 4.0627e-5)),
                             frac_p_lt_0p01=float(np.mean(x.rv_chisq_pvalue < 0.01)), median_gof=float(np.median(x.rv_renormalised_gof)),
                             frac_gof_ge_3=float(np.mean(x.rv_renormalised_gof >= 2.9986)),
                             median_amp=float(np.median(x.rv_amplitude_robust)), frac_amp_ge_target=float(np.mean(x.rv_amplitude_robust >= 3.2247)))
json.dump(out, open("rv01_deblend.json", "w"), indent=1, default=str); print(json.dumps(out, indent=1, default=str))
