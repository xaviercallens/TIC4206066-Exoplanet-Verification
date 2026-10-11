"""T2_A_shape r03: small derived-quantity checks (no fitting).
(1) C058: T14 implied by the pubprep timing-template parameters vs the audit posterior-median T14 and the T14 of the
    posterior-median parameters; (2) Rp from k with R* = 0.623 +- 0.019 (MC); (3) Teq / insolation conventions (C064);
(4) V-shape / phase-0.5 sign conventions read from the stored audit s10 JSON; (5) S6-drop changes in sigma units."""
import json, os
import numpy as np
ROOT = "/Users/paulsmith/space <3/tess-jev"
W = os.path.join(ROOT, "addenda/tic4206066_FINAL_RECONCILIATION/_work/T2_A_shape")
G = 6.67430e-8; MSUN = 1.98841e33; RSUN = 6.957e10; RHOSUN = MSUN / (4 / 3 * np.pi * RSUN ** 3); AU = 1.495978707e13
R_OVER_RE = 6.957e8 / 6.3781e6
def ars(rho, P): return (G * rho * RHOSUN * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
def t14(P, k, b, a):
    inc = np.arccos(b / a); return P / np.pi * np.arcsin(np.sqrt((1 + k) ** 2 - b ** 2) / (a * np.sin(inc)))
out = {}
P = 3.1827853
# (1) template parameters of pubprep a1_timing_depths.py L35
a_t = ars(2.507, P); out["C058_template_T14_h"] = float(t14(P, 0.02116, 0.35, a_t) * 24)
FP = json.load(open(os.path.join(ROOT, "addenda/tic4206066_PUBLICATION_AUDIT/FINAL_PARAMETERS.json")))["signal_A"]["shape_mann_prior"]
k, b, aR = FP["k"]["value"], FP["b"]["value"], FP["a_over_Rs"]["value"]
out["C058_T14_of_posterior_medians_h"] = float(t14(P, k, b, aR) * 24)
out["C058_audit_posterior_median_T14_h"] = FP["T14_h"]["value"]
out["C058_ars_template"] = float(a_t)
# (2) Rp(k) Monte Carlo
rng = np.random.default_rng(1)
Rs = rng.normal(0.623, 0.019, 400000)
for lab, kk, ek in (("k0.0212_0.0007", 0.0212, 0.0007), ("k0.021_0.001", 0.021, 0.001), ("k0.0204_0.0008", 0.0204, 0.0008),
                    ("k0.0204_0.0010", 0.0204, 0.0010), ("k0.0195", 0.0195, 0.0007), ("k0.02116", 0.02116, 0.0007)):
    rp = rng.normal(kk, ek, Rs.size) * Rs * R_OVER_RE
    out[f"Rp_{lab}"] = [float(np.median(rp)), float(np.std(rp))]
# LD-corrected sqrt(depth): mid-transit depth / k^2 for the template (quadratic LD, b)
def ldfac(u1, u2, b):
    mu = np.sqrt(1 - b ** 2); I = 1 - u1 * (1 - mu) - u2 * (1 - mu) ** 2; return I / (1 - u1 / 3 - u2 / 6)
out["LD_factor_b0.33_u0.27_0.29"] = float(ldfac(0.27, 0.29, 0.33))
out["LD_factor_b0.35_u0.316_0.313"] = float(ldfac(0.316, 0.313, 0.35))
for kk in (0.0204, 0.021, 0.0212):
    out[f"mid_depth_ppm_k{kk}_b0.33"] = float(kk ** 2 * ldfac(0.27, 0.29, 0.33) * 1e6)
# (3) Teq / insolation; a from Kepler with M* = 0.606
Ms, Rst = 0.606, 0.623
a_au = ((G * Ms * MSUN * (P * 86400) ** 2) / (4 * np.pi ** 2)) ** (1 / 3) / AU
out["a_AU"] = float(a_au); out["a_over_Rs_kepler"] = float(a_au * AU / (Rst * RSUN))
for Teff in (4012.0, 4020.0, 4030.0):
    out[f"Teq_A0_full_Teff{int(Teff)}"] = float(Teff * np.sqrt(Rst * RSUN / (2 * a_au * AU)))
    out[f"Teq_A0.3_full_Teff{int(Teff)}"] = float(Teff * 0.7 ** 0.25 * np.sqrt(Rst * RSUN / (2 * a_au * AU)))
    out[f"S_Searth_Teff{int(Teff)}"] = float(Rst ** 2 * (Teff / 5772.0) ** 4 / a_au ** 2)
cand = json.load(open(os.path.join(ROOT, "candidates/tic4206066_p3.182945/candidate.json")))["leo_vetter"]
out["LEO_Teq_reported"] = cand["Teq"]; out["LEO_aRs"] = cand["aRs"]
# (4) sign conventions in the stored audit s10 output (box() = 1 - weighted mean flux, s10_eb_tests.py L37-42)
s10 = json.load(open(os.path.join(ROOT, "addenda/tic4206066_PUBLICATION_AUDIT/out/s10_eb_tests.json")))["A_pipe"]
out["phase05_box_depth_ppm"] = s10["secondary_0p5"]["depth_ppm"]; out["phase05_flux_ppm"] = -s10["secondary_0p5"]["depth_ppm"]
out["phase05_err_emp_robust_ppm"] = s10["secondary_0p5"]["err_emp_ppm"]; out["phase05_upper3sigma_ppm"] = s10["secondary_0p5"]["upper_3sigma_ppm"]
out["P3_box_depths_ppm"] = {k: v["depth_ppm"] for k, v in s10["P_over_3_check"].items()}
out["vshape_dchi2_V_minus_U_rescaled"] = s10["shape"]["dchi2_V_minus_U"]; out["vshape_chi2_scale"] = s10["shape"]["chi2_scale"]
out["vshape_dchi2_V_minus_U_unscaled"] = s10["shape"]["dchi2_V_minus_U"] * s10["shape"]["chi2_scale"]
# (5) S6 drop, like-for-like, U/referee/resolve_R-S6/out/R-S6_summary.json
RS = json.load(open(os.path.join(ROOT, "addenda/opus55_final_audit_tic4206066/referee/resolve_R-S6/out/R-S6_summary.json")))
jf = RS["joint_fits"]; ch = {}
for tag_all, tag_sub in (("S6_32_98_raw_jit_ss615", "S32_98_raw_jit_ss615"), ("S6_32_98_bw_jit_ss615", "S32_98_bw_jit_ss615")):
    A_, B_ = jf[tag_all], jf[tag_sub]; row = {}
    for q in ("k", "b", "T14_h"):
        sa = 0.5 * (A_[q][2] - A_[q][0]); sb = 0.5 * (B_[q][2] - B_[q][0]); d = B_[q][1] - A_[q][1]
        row[q] = dict(all=A_[q][1], noS6=B_[q][1], diff=d, n_sigma_all=d / sa, n_sigma_noS6=d / sb)
    d = B_["Rp_mean_std"][0] - A_["Rp_mean_std"][0]; row["Rp"] = dict(all=A_["Rp_mean_std"][0], noS6=B_["Rp_mean_std"][0], n_sigma=d / A_["Rp_mean_std"][1])
    ch[tag_all] = row
c = RS["s6_checks"]
ch["phase05"] = dict(all=c["phase05"]["all3"], noS6=c["phase05"]["S32_S98"],
                     n_sigma_all=(c["phase05"]["S32_S98"]["depth_ppm"] - c["phase05"]["all3"]["depth_ppm"]) / c["phase05"]["all3"]["e_ppm"],
                     n_sigma_noS6=(c["phase05"]["S32_S98"]["depth_ppm"] - c["phase05"]["all3"]["depth_ppm"]) / c["phase05"]["S32_S98"]["e_ppm"])
ch["odd_even"] = c["odd_even"]; ch["vshape_dchi2"] = dict(all=c["vshape"]["all3"]["dchi2"], noS6=c["vshape"]["S32_S98"]["dchi2"])
ch["rho977_uniform"] = dict(all=jf["S6_32_98_raw_jit_ss615_rhouni"]["rho_q97_7"], noS6=jf["S32_98_raw_jit_ss615_rhouni"]["rho_q97_7"])
out["S6_drop"] = ch
json.dump(out, open(os.path.join(W, "out/r03_derived.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
