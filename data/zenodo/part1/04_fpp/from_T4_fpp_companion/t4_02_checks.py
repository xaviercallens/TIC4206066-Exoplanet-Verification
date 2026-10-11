"""Minimal re-derivations for topic T4 (read-only inputs; outputs only in this folder).
1) Gaia RV: sigma(V_j) from radial_velocity_error (Gaia DR3 datamodel, rv_method_used=1), mass-function limits.
2) Share of default FPP_A from companions > 0.38 Msun, from U/08_fpp/out_fpp_summary.csv (two definitions).
3) Default-FPP ranges, bin-width sensitivity, NFPP maxima, from t4_01_fpp_run_census.csv.
4) Overluminosity flux-ratio arithmetic (O + 2 sigma).
"""
import json, os, csv
import numpy as np
from scipy.optimize import brentq
H = os.path.dirname(os.path.abspath(__file__))
ROOT = "/Users/paulsmith/space <3/tess-jev"
U = ROOT + "/addenda/opus55_final_audit_tic4206066"
out = {}
# ---- 1) RV
G = 6.674e-11; Msun = 1.989e30; day = 86400.0
eps, N = 0.2206730842590332, 33   # U/09_companion/out/rv_limits.json radial_velocity_error, rv_nb_transits
sig_vj = np.sqrt((eps**2 - 0.11**2) * 2 * N / np.pi)
def msini(K_kms, P_d, M1):
    f = (K_kms * 1e3)**3 * P_d * day / (2 * np.pi * G) / Msun   # mass function (sin i = 1, e = 0)
    return brentq(lambda m2: m2**3 / (M1 + m2)**2 - f, 1e-6, 5)
def K_of(m2, P_d, M1):
    return (2 * np.pi * G / (P_d * day))**(1/3) * m2 * Msun / ((M1 + m2) * Msun)**(2/3) / 1e3
rv = {"sigma_Vj_kms": sig_vj}
for M1 in (0.606, 0.6062):
    for P in (3.1827856, 6.365571, 11.13274):
        for K in (3.2246542, 3.2, 1.16):
            rv[f"Msini_M1{M1}_P{P}_K{K}"] = round(msini(K, P, M1), 5)
        rv[f"K_for_0.08_M1{M1}_P{P}"] = round(K_of(0.08, P, M1), 2)
out["rv"] = rv
# ---- 2) share > 0.38 Msun
rows = list(csv.DictReader(open(U + "/08_fpp/out_fpp_summary.csv")))
share = {}
for r in rows:
    if r["signal"] != "A" or r["molusc"]: continue
    fpp = float(r["FPP_mean"]); cut = float(r["FPP_comp_mass_le_0.38_mean"])
    share[r["variant"]] = dict(FPP=fpp, FPP_le038=cut, share_by_FPP_drop=1 - cut / fpp, STP_wfrac_gt038=float(r["STP_wfrac_mass_gt_0.38"]))
out["share_gt038_A_default"] = share
s1 = [v["share_by_FPP_drop"] for v in share.values()]; s2 = [v["STP_wfrac_gt038"] for v in share.values()]
out["share_gt038_A_range"] = dict(by_FPP_drop=[min(s1), max(s1)], STP_wfrac=[min(s2), max(s2)], FPP_le038=[min(v["FPP_le038"] for v in share.values()), max(v["FPP_le038"] for v in share.values())])
B = {r["variant"]: dict(FPP=float(r["FPP_mean"]), FPP_le038=float(r["FPP_comp_mass_le_0.38_mean"]), STP_wfrac_gt038=float(r["STP_wfrac_mass_gt_0.38"])) for r in rows if r["signal"] == "B"}
out["B_variants"] = B
# ---- 3) census
cen = list(csv.DictReader(open(H + "/t4_01_fpp_run_census.csv")))
def g(v, src_sub):
    for r in cen:
        if r["variant"] == v and src_sub in r["source"]: return r
    return None
defaultA = {  # S32+S98 folds, default priors, 5-10 min bins
    "pubprep_ownfold_A (6 min)": g("pubprep_ownfold_A", "pubprep"),
    "audit A_base_S32+S98 (6 min)": g("A_base_S32+S98", "PUBLICATION_AUDIT"),
    "audit A_repro_S32+S98 (same fold, Teff 4080, exptime 200 s)": g("A_repro_S32+S98", "PUBLICATION_AUDIT"),
    "audit A_tic_S32+S98": g("A_tic_S32+S98", "PUBLICATION_AUDIT"),
    "U08 A_base (6 min)": g("A_base", "08_fpp/runs"),
    "U08 A_T0m": g("A_T0m", "08_fpp/runs"), "U08 A_T0p": g("A_T0p", "08_fpp/runs"),
    "U08 A_ap5x5": g("A_ap5x5", "08_fpp/runs"), "U08 A_sigformal": g("A_sigformal", "08_fpp/runs"),
    "U08 A_teff4012": g("A_teff4012", "08_fpp/runs"), "U08 A_tic": g("A_tic", "08_fpp/runs"),
    "U08 A_10min": g("A_10min", "08_fpp/runs"),
    "U08v default 5 min": g("default_A_v_A_S32+S98_5min", "runs_v"), "U08v default 10 min": g("default_A_v_A_S32+S98_10min", "runs_v"),
    "R5 default_teff4012 adopted fold": g("default_teff4012_A_adoptedfold", "resolve_R5"),
    "R5 default adopted fold (2-seed repro)": g("default_A_adoptedfold", "resolve_R5"),
}
tab = {k: dict(n=int(v["n_runs"]), mean=float(v["FPP_mean"]), min=float(v["FPP_min"]), max=float(v["FPP_max"]), NFPP_max=float(v["NFPP_max"])) for k, v in defaultA.items()}
out["A_default_S32S98"] = tab
means = [v["mean"] for v in tab.values()]; ns = [v["n"] for v in tab.values()]
out["A_default_S32S98_summary"] = dict(mean_range=[min(means), max(means)], run_counts=sorted(set(ns)),
    per_run_range=[min(v["min"] for v in tab.values()), max(v["max"] for v in tab.values())],
    NFPP_max=max(v["NFPP_max"] for v in tab.values()))
single = {k: g(k, "PUBLICATION_AUDIT") for k in ("A_base_S32only", "A_base_S6only", "A_base_S98only")}
single.update({"U08 " + k: g(k, "08_fpp/runs") for k in ("A_S32only", "A_S98only")})
out["A_default_single_sector"] = {k: dict(n=int(v["n_runs"]), mean=float(v["FPP_mean"]), NFPP_max=float(v["NFPP_max"])) for k, v in single.items()}
out["NFPP_max_all_default_A"] = max([float(r["NFPP_max"]) for r in cen if r.get("NFPP_max") not in ("", None) and r["signal"] == "A"])
out["binwidth"] = dict(U08_10_over_6=tab["U08 A_10min"]["mean"] / tab["U08 A_base (6 min)"]["mean"],
                       U08v_10_over_5=tab["U08v default 10 min"]["mean"] / tab["U08v default 5 min"]["mean"],
                       B_U08_10_over_6=float(g("B_10min", "08_fpp/runs")["FPP_mean"]) / float(g("B_base", "08_fpp/runs")["FPP_mean"]))
# B default ranges
Bdef = [r for r in cen if r["signal"] == "B" and r.get("FPP_mean") and not any(s in r["variant"] for s in ("molusc", "gaia", "ol0.38", "trunc"))]
out["B_default_all"] = {r["variant"] + " @ " + r["source"].split("/")[1]: (int(r["n_runs"]), round(float(r["FPP_mean"]), 4)) for r in Bdef}
# ---- 4) overluminosity arithmetic
for sig in (0.052, 0.0534, 0.068, 0.07, 0.098, 0.10):
    Ol = 0.003 + 2 * sig
    out.setdefault("overlum_fluxratio_2sig", {})[str(sig)] = dict(O_lim=Ol, fG_max=10**(0.4 * Ol) - 1, dG_comp=-2.5 * np.log10(10**(0.4 * Ol) - 1))
json.dump(out, open(H + "/t4_02_checks.json", "w"), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
