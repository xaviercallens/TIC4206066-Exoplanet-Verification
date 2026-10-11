"""T3: robustness of the combined 3.18 d (A) position to the S98 and S32 per-set choices, and distances of the
depth-capable Gaia DR3 stars and of the 11" neighbour. Inputs are stored outputs only (read-only):
  U 03_pixels/summary.json (sets), U referee/resolve_R3 per-set adopted errors (r2_combine.log 'per-set errors'),
  audit out/s12c_clean_localisation.json (S98 clean A), pubprep analysis/out/a5_summary.json (S98_ffi A),
  this dir t3_s98_driver.json (audit recipe with U clean29 mask), U 03_pixels/neighbours_vs_fitted_position.csv.
Per-set offsets are REAL offsets (not injection-bias-corrected) except audit/pubprep (their stored values are
real minus injection mean, as they publish them)."""
import json, numpy as np, pandas as pd
from scipy.stats import chi2 as C2, norm
R = "/Users/paulsmith/space <3/tess-jev/addenda/"
U = R + "opus55_final_audit_tic4206066/"; W = R + "tic4206066_FINAL_RECONCILIATION/_work/T3_pixels_localisation/"
S = json.load(open(U + "03_pixels/summary.json"))["sets"]
AUD = json.load(open(R + "tic4206066_PUBLICATION_AUDIT/out/s12c_clean_localisation.json"))
PUB = json.load(open(R + "rnaas_tic4206066_pubprep/analysis/out/a5_summary.json"))["localisation"]
DRV = json.load(open(W + "t3_s98_driver.json"))
NB = pd.read_csv(U + "03_pixels/neighbours_vs_fitted_position.csv")

def diag(e): return np.diag(np.square(e))
# per-set entries: (offset E,N ; covariance)
sets = {
 # U S32 SPOC 2-min POS_CORR-decorrelated; error = R3 adopted max(inj SD, jackknife, robust boot) -> (3.08, 2.78)
 "S32PC_U": (S["S32_tpf_poscorr_A"]["real_off"], diag([3.0811, 2.7826])),
 "S32tspocPC_U": (S["S32_tspoc_poscorr_A"]["real_off"], np.array(S["S32_tspoc_poscorr_A"]["inj_cov"])),
 "S32raw_U": (S["S32_tpf_A"]["real_off"], np.array(S["S32_tpf_A"]["inj_cov"])),
 "S32raw_pubprep": (PUB["A"]["per_set"]["S32_tpf"]["offset_arcsec"], diag(PUB["A"]["per_set"]["S32_tpf"]["sigma_arcsec"])),
 # U S98 clean29 9 events, R3 adopted errors (5.27, 4.52)
 "S98_U_clean29": (S["S98_clean29_A"]["real_off"], diag([5.2686, 4.5166])),
 "S98_U_clean_deg2": (S["S98_clean_deg2_A"]["real_off"], np.array(S["S98_clean_deg2_A"]["inj_cov"])),
 "S98_U_clean14": (S["S98_clean_A"]["real_off"], np.array(S["S98_clean_A"]["inj_cov"])),
 "S98_audit_9ev": (AUD["clean_S98"]["A"]["real_sky_arcsec"], np.array(AUD["clean_S98"]["A"]["sky_cov_arcsec2"])),
 "S98_pubprep_10ev": (PUB["A"]["per_set"]["S98_ffi"]["offset_arcsec"], diag(PUB["A"]["per_set"]["S98_ffi"]["sigma_arcsec"])),
}
def comb(keys):
    Si = np.zeros((2, 2)); Sx = np.zeros(2); chi_b = 0
    for k in keys:
        m, C = np.array(sets[k][0], float), np.array(sets[k][1], float); Ci = np.linalg.inv(C); Si += Ci; Sx += Ci @ m
    Cc = np.linalg.inv(Si); mu = Cc @ Sx
    for k in keys:
        m, C = np.array(sets[k][0], float), np.array(sets[k][1], float); chi_b += float((m - mu) @ np.linalg.inv(C) @ (m - mu))
    c2 = float(mu @ Si @ mu); p = float(C2.sf(c2, 2))
    return dict(sets=keys, mean=mu.tolist(), sig=np.sqrt(np.diag(Cc)).tolist(), cov=Cc.tolist(), chi2_vs_target=c2, p=p,
                sigma_equiv_2d=float(norm.isf(p / 2)), z_E=float(mu[0] / np.sqrt(Cc[0, 0])), z_N=float(mu[1] / np.sqrt(Cc[1, 1])),
                chi2_between=chi_b, dof_between=2 * (len(keys) - 1))
cap = NB[NB.dreqA_min <= 1.0][["source_id", "TIC", "Tmag", "sep", "dE", "dN", "dreqA_min"]]
n11 = NB.iloc[0]
def dist(mu, C, row):
    d = np.array([row.dE, row.dN]) - np.array(mu); return float(np.sqrt(d @ np.linalg.solve(np.array(C), d)))
variants = {
 "U_headline_S32PC+S98U": ["S32PC_U", "S98_U_clean29"],
 "S32PC+S98audit9ev": ["S32PC_U", "S98_audit_9ev"],
 "S32PC+S98pubprep10ev": ["S32PC_U", "S98_pubprep_10ev"],
 "S32PC+S98U_clean_deg2": ["S32PC_U", "S98_U_clean_deg2"],
 "S32PC+S98U_clean14": ["S32PC_U", "S98_U_clean14"],
 "S32tspocPC+S98U": ["S32tspocPC_U", "S98_U_clean29"],
 "S32raw_pubprep+S98pubprep (package)": ["S32raw_pubprep", "S98_pubprep_10ev"],
 "S32raw_U+S98audit9ev": ["S32raw_U", "S98_audit_9ev"],
}
out = {"per_set": {k: dict(off=list(map(float, v[0])), sig=np.sqrt(np.diag(v[1])).tolist()) for k, v in sets.items()}, "variants": {}}
for name, keys in variants.items():
    r = comb(keys)
    ds = [(int(row.TIC), round(row.sep, 1), dist(r["mean"], r["cov"], row)) for row in cap.itertuples()]
    r["capable_star_mahalanobis"] = ds; r["min_capable"] = min(x[2] for x in ds)
    r["nearest_gaia_11as_mahalanobis"] = dist(r["mean"], r["cov"], n11)
    rC = np.diag([9.0, 9.0])   # rounded +-3" form
    r["min_capable_rounded_pm3"] = min(dist(r["mean"], rC, row) for row in cap.itertuples())
    out["variants"][name] = r
    print(f"{name:38s} ({r['mean'][0]:+5.2f}+-{r['sig'][0]:.2f}, {r['mean'][1]:+5.2f}+-{r['sig'][1]:.2f}) zE={r['z_E']:.2f} chi2={r['chi2_vs_target']:.2f} p={r['p']:.3f} "
          f"between={r['chi2_between']:.2f}/{r['dof_between']} minCap={r['min_capable']:.1f} (pm3: {r['min_capable_rounded_pm3']:.1f}) 11as={r['nearest_gaia_11as_mahalanobis']:.1f}")
# rounded printed forms
for lab, mu, sg in [("U printed (+5+-3,0+-3)", [5, 0], [3, 3]), ("pubprep printed (+8+-5,-8+-7)", [8, -8], [5, 7]), ("audit (+2.5,-7.5)+-5.7", [2.509, -7.534], [5.684, 5.692])]:
    C = np.diag(np.square(sg)); m = [dist(mu, C, row) for row in cap.itertuples()]
    out.setdefault("printed", {})[lab] = dict(min_capable=min(m), per_star=[(int(r.TIC), round(r.sep, 1), round(x, 2)) for r, x in zip(cap.itertuples(), m)],
                                            nearest_11as=dist(mu, C, n11), z_E=mu[0] / sg[0], p2d=float(C2.sf((mu[0] / sg[0]) ** 2 + (mu[1] / sg[1]) ** 2, 2)))
    print(lab, "min capable", round(min(m), 2), "11as", round(dist(mu, C, n11), 2), "zE", round(mu[0] / sg[0], 2), "p2d", round(out["printed"][lab]["p2d"], 3))
out["capable_list"] = cap.to_dict(orient="records"); out["nearest_gaia_11as"] = n11[["source_id", "TIC", "Tmag", "sep", "dE", "dN", "dreqA_min"]].to_dict()
json.dump(out, open(W + "t3_combine.json", "w"), indent=1, default=float)
