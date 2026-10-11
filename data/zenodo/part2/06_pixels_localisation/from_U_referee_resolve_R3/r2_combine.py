"""R3 resolution, step 2: combine per-sector positions under several error models; S6 inclusion test;
11.13 d position; distances of depth-capable Gaia/TIC stars from the fitted position.

Inputs (read-only): package a5_summary.json (manuscript numbers), auditor 03_pixels/summary.json and
poscorr_decorr_*.json (positions + injection covariances + jackknife), verifier 03_pixels/verify/v1_loc.json,
this folder's out/r1_boot.json (5000-resample event bootstrap, jackknife), 03_pixels/neighbours_vs_fitted_position.csv,
00_data TIC 8.2 cone (J2000 positions, for the 13_figure 'TIC-position' variant).
"""
import os, json
import numpy as np, pandas as pd
from scipy import stats

D = "/Users/paulsmith/space <3/tess-jev/addenda"
A = D + "/opus55_final_audit_tic4206066"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
a5 = json.load(open(D + "/rnaas_tic4206066_pubprep/analysis/out/a5_summary.json"))["localisation"]
aud = json.load(open(A + "/03_pixels/summary.json"))
pc32 = json.load(open(A + "/03_pixels/poscorr_decorr_S32_tpf.json"))
pc6 = json.load(open(A + "/03_pixels/poscorr_decorr_S6_tspoc_S6_tspoc_dumpx0.5h_S32_tspoc.json"))
r1 = {(r["key"], r["sig"], r["usepc"]): r for r in json.load(open(os.path.join(OUT, "r1_boot.json")))}
res = {}


def comb(items):
    """items: list of (mean[2], cov[2x2]); inverse-covariance weighted mean; chi2 vs target (0,0)."""
    W = sum(np.linalg.inv(np.array(c)) for _, c in items); C = np.linalg.inv(W)
    m = C @ sum(np.linalg.inv(np.array(c)) @ np.array(x) for x, c in items)
    chi2 = float(m @ W @ m); p = float(stats.chi2.sf(chi2, 2))
    scat = float(sum((np.array(x) - m) @ np.linalg.inv(np.array(c)) @ (np.array(x) - m) for x, c in items))
    return dict(mean=m.tolist(), sig=np.sqrt(np.diag(C)).tolist(), cov=C.tolist(), chi2_vs_target=chi2, p=p,
                sigma_equiv_2d=float(stats.norm.isf(p / 2)), z_E=float(m[0] / np.sqrt(C[0, 0])), z_N=float(m[1] / np.sqrt(C[1, 1])),
                chi2_between=scat, dof_between=2 * (len(items) - 1))


def diag(s): return [[s[0] ** 2, 0], [0, s[1] ** 2]]


def jk_se(jk):
    jk = np.array([j[1:3] for j in jk]); n = len(jk)
    return np.sqrt((n - 1) / n * ((jk - jk.mean(0)) ** 2).sum(0))


# ---- 1. package / manuscript reproduction -------------------------------------------------------------
LA = a5["A"]
ps = LA["per_set"]
pk = [(ps[k]["offset_arcsec"], diag(ps[k]["sigma_arcsec"])) for k in ("S32_tpf", "S98_ffi")]
res["package_S32tpf+S98"] = comb(pk)
res["package_S32tpf+S98_file"] = LA["combined_S32tpf_S98clean"]
res["package_S6+S32tpf+S98"] = comb(pk + [(ps["S6_ffi"]["offset_arcsec"], diag(ps["S6_ffi"]["sigma_arcsec"]))])
res["package_S6_point"] = dict(off=ps["S6_ffi"]["offset_arcsec"], sig=ps["S6_ffi"]["sigma_arcsec"], z=ps["S6_ffi"]["z"],
                               maha=ps["S6_ffi"]["mahalanobis"], dump_h=ps["S6_ffi"]["transit_to_dump_h"])

# ---- 2. corrected 3.18 d position: S32 (POS_CORR-decorrelated 2-min TPF) + S98 (clean, 9 events) ------
s32 = aud["sets"]["S32_tpf_poscorr_A"]; s98 = aud["sets"]["S98_clean29_A"]
b32 = r1[("S32_tpf", "A", True)]; b98 = r1[("S98_ffi", "A", False)]
inj32, inj98 = np.array(s32["inj_std"]), np.array(s98["inj_std"])
jk32_aud = jk_se(pc32["S32_tpf_poscorr_A"]["jackknife"])
ev32 = np.max([inj32, b32["jackknife_se"], b32["boot_robust_sd"], jk32_aud], 0)
ev98 = np.max([inj98, b98["jackknife_se"], b98["boot_robust_sd"]], 0)
res["per_set_errors_3.18d"] = dict(
    S32PC=dict(off_auditor=s32["real_off"], off_verifier_code=b32["off"], inj_sd=inj32.tolist(), jk_se_auditor=jk32_aud.tolist(),
               jk_se_this=b32["jackknife_se"], boot_sd=b32["boot_sd"], boot_robust_sd=b32["boot_robust_sd"], adopted_max=ev32.tolist()),
    S98=dict(off_auditor=s98["real_off"], off_verifier_code=b98["off"], inj_sd=inj98.tolist(), inj_bias=s98["inj_mean"],
             jk_se_this=b98["jackknife_se"], boot_sd=b98["boot_sd"], boot_robust_sd=b98["boot_robust_sd"], adopted_max=ev98.tolist()))
res["A_injcov_auditor"] = comb([(s32["real_off"], s32["inj_cov"]), (s98["real_off"], s98["inj_cov"])])
res["A_eventscatter_auditor_pos"] = comb([(s32["real_off"], diag(ev32)), (s98["real_off"], diag(ev98))])
res["A_eventscatter_verifier_pos"] = comb([(b32["off"], diag(ev32)), (b98["off"], diag(ev98))])
res["A_bootcov_verifier_pos"] = comb([(b32["off"], b32["boot_cov"]), (b98["off"], b98["boot_cov"])])
bias98 = np.array(s98["off_bias"]) if "off_bias" in s98 else np.array(s98["inj_mean"])
res["A_eventscatter_biascorr"] = comb([(np.array(s32["real_off"]) - np.array(s32["inj_mean"]), diag(ev32)),
                                       (np.array(s98["real_off"]) - bias98, diag(ev98))])
res["A_round_pm3"] = comb([([5.0, 0.0], diag([3.0, 3.0]))])
# raw (not decorrelated) S32 + S98 for reference
b32raw = r1[("S32_tpf", "A", False)]
res["A_raw_S32_eventscatter"] = comb([(b32raw["off"], diag(np.max([b32raw["jackknife_se"], b32raw["boot_robust_sd"]], 0))),
                                      (b98["off"], diag(ev98))])

# ---- 3. Sector 6 -------------------------------------------------------------------------------------
s6pc = aud["sets"]["S6_tspoc_poscorr_A"]
jk6 = pc6["S6_tspoc_poscorr_A"]["jackknife"]
ev6 = np.max([s6pc["inj_std"], jk_se(jk6)], 0)
res["S6"] = dict(raw_this=r1[("S6_ffi", "A", False)]["off"], raw_this_jk_se=r1[("S6_ffi", "A", False)]["jackknife_se"],
                 raw_this_jackknife=r1[("S6_ffi", "A", False)]["jackknife"],
                 dumpx05h_auditor=aud["sets"]["S6_dumpx0.5h_A"]["real_off"], pc_auditor=s6pc["real_off"], pc_inj_sd=s6pc["inj_std"],
                 pc_jackknife=jk6, pc_jk_se=jk_se(jk6).tolist(), package_point=ps["S6_ffi"]["offset_arcsec"])
res["A_with_S6pc_eventscatter"] = comb([(s32["real_off"], diag(ev32)), (s98["real_off"], diag(ev98)), (s6pc["real_off"], diag(ev6))])
res["A_with_S6pc_injcov"] = comb([(s32["real_off"], s32["inj_cov"]), (s98["real_off"], s98["inj_cov"]), (s6pc["real_off"], s6pc["inj_cov"])])

# ---- 4. 11.13 d --------------------------------------------------------------------------------------
LB = a5["B"]
res["B_package"] = LB["combined_S32tpf_S98clean"]
res["B_package_S32_point"] = LB["per_set"]["S32_tpf"]
sb98 = aud["sets"]["S98_clean29_B"]; bb98 = r1[("S98_ffi", "B", False)]
evb98 = np.max([sb98["inj_std"], bb98["jackknife_se"]], 0)
res["B_S98only_auditor"] = comb([(sb98["real_off"], diag(evb98))])
res["B_S98only_verifiercode"] = comb([(bb98["off"], diag(evb98))])
sb32 = aud["sets"]["S32_tpf_poscorr_B"]
res["B_S32PC+S98_auditor"] = comb([(sb98["real_off"], diag(evb98)), (sb32["real_off"], sb32["inj_cov"])])
res["B_S32_events"] = dict(auditor_PC_events=[j[0] for j in pc32["S32_tpf_poscorr_B"]["jackknife"]],
                           verifier_clean_events=r1[("S32_tpf", "B", True)]["tcs"],
                           verifier_single_clean_event_PC_off=r1[("S32_tpf", "B", True)]["off"])

# ---- 5. depth-capable stars vs fitted position ---------------------------------------------------------
nb = pd.read_csv(A + "/03_pixels/neighbours_vs_fitted_position.csv")
cap = nb[nb.dreqA_min <= 1].copy()
pos = {"manuscript(+8+-5,-8+-7)": ([8, -8], diag([5, 7])), "package_exact": (LA["combined_S32tpf_S98clean"]["offset_arcsec"], diag(LA["combined_S32tpf_S98clean"]["sigma_arcsec"])),
       "corrected_eventscatter": (res["A_eventscatter_auditor_pos"]["mean"], res["A_eventscatter_auditor_pos"]["cov"]),
       "corrected_rounded(+5+-3,0+-3)": ([5, 0], diag([3, 3]))}
dist = {}
for name, (m, C) in pos.items():
    Ci = np.linalg.inv(np.array(C)); m = np.array(m)
    d = [float(np.sqrt((np.array([r.dE, r.dN]) - m) @ Ci @ (np.array([r.dE, r.dN]) - m))) for r in cap.itertuples()]
    k = int(np.argmin(d))
    dist[name] = dict(min_maha=d[k], tic=int(cap.TIC.iloc[k]), sep=float(cap.sep.iloc[k]),
                      all=[(int(t), round(float(s), 1), round(x, 2)) for t, s, x in zip(cap.TIC, cap.sep, d)])
# TIC J2000 positions variant (as used by the 13_figure audit)
tic = pd.read_csv(A + "/00_data/catalogs/tic82_cone300as.csv")
t0 = tic[tic.ID == 4206066].iloc[0]
tic["dE"] = (tic.ra - t0.ra) * np.cos(np.radians(t0.dec)) * 3600; tic["dN"] = (tic.dec - t0.dec) * 3600
tc = tic[tic.ID.isin(cap.TIC.astype(int))]
for name in ("manuscript(+8+-5,-8+-7)", "corrected_rounded(+5+-3,0+-3)"):
    m, C = pos[name]; Ci = np.linalg.inv(np.array(C)); m = np.array(m)
    d = [float(np.sqrt((np.array([r.dE, r.dN]) - m) @ Ci @ (np.array([r.dE, r.dN]) - m))) for r in tc.itertuples()]
    k = int(np.argmin(d)); dist[name + "_TICJ2000"] = dict(min_maha=d[k], tic=int(tc.ID.iloc[k]),
                                                             sep=float(np.hypot(tc.dE.iloc[k], tc.dN.iloc[k])))
res["capable_star_distance"] = dist
res["offtarget_injection_recovery"] = {k: dict(true=v["true"], mean=v["mean"], std=v["std"], ratio=v["recovered_over_true"],
                                               dchi2_target_p5=v["dchi2_target_p5"]) for k, v in aud["shifts"].items()}

json.dump(res, open(os.path.join(OUT, "r2_combine.json"), "w"), indent=1, default=float)


def show(k):
    r = res[k]
    print(f"{k:34s} ({r['mean'][0]:+.2f}+-{r['sig'][0]:.2f}, {r['mean'][1]:+.2f}+-{r['sig'][1]:.2f})  chi2={r['chi2_vs_target']:.2f} "
          f"p={r['p']:.3g} ({r['sigma_equiv_2d']:.2f} sig 2D; zE={r['z_E']:.2f})  between={r['chi2_between']:.2f}/{r['dof_between']}")


for k in ("package_S32tpf+S98", "package_S6+S32tpf+S98", "A_injcov_auditor", "A_eventscatter_auditor_pos", "A_eventscatter_verifier_pos",
          "A_bootcov_verifier_pos", "A_eventscatter_biascorr", "A_round_pm3", "A_raw_S32_eventscatter", "A_with_S6pc_eventscatter",
          "A_with_S6pc_injcov", "B_S98only_auditor", "B_S98only_verifiercode", "B_S32PC+S98_auditor"):
    show(k)
print("per-set errors", json.dumps(res["per_set_errors_3.18d"], default=float))
print("S6", json.dumps({k: v for k, v in res["S6"].items() if k != "raw_this_jackknife"}, default=float))
print("B package", res["B_package"]); print("B S32 events", res["B_S32_events"])
for k, v in dist.items(): print("capable-star distance", k, {a: b for a, b in v.items() if a != "all"})
