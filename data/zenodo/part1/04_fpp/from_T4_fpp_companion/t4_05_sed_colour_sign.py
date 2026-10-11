"""C024: sign of the G-Ks SED residual. resid = -2.5 log10(obs/model) = m_obs - m_model (companion_check sed_residuals.py l.25)."""
import json, os
C = "/Users/paulsmith/space <3/tess-jev/addenda/companion_check_tic4206066/out/sed_residuals.json"
d = json.load(open(C)); f = d["filters"]; r = d["per_model"]["BMA"]["resid_mag"]
iG = [i for i, x in enumerate(f) if "G" in x and "BP" not in x and "RP" not in x][0]
iK = [i for i, x in enumerate(f) if "Ks" in x or x.endswith("_K")][0]
dGK = r[iG] - r[iK]   # (G-Ks)_obs - (G-Ks)_model
pred = {"M2V_0.44": 0.115, "M4V_0.23": 0.067, "M5V_0.16": 0.038}   # companion_check SUMMARY.md L59-L60
out = dict(filters=[f[iG], f[iK]], resid_G=r[iG], resid_Ks=r[iK], obs_minus_model_GmKs=dGK,
           meaning="positive = observed redder than single-star model",
           n_sigma_pred_minus_obs_err0p03={k: round((v - dGK) / 0.03, 2) for k, v in pred.items()})
json.dump(out, open(os.path.dirname(os.path.abspath(__file__)) + "/t4_05_sed_colour_sign.json", "w"), indent=1)
print(json.dumps(out, indent=1))
