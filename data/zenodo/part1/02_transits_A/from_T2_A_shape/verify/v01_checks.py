"""T2_A_shape verifier v01: cheap read-back checks of the reconciler's numbers (no fitting).
(1) s17 all-sector weighted means per detrending; (2) SD/mean of the seven edge-excluded k fits and all-23 k fits;
(3) pubprep chi2 of constant (weighted/unweighted); (4) S6 SAP-PDCSAP paired per-event difference (s09 rows);
(5) S6 transit-to-dump distances with the manuscript ephemeris; (6) no-S6 values vs the revised manuscript values."""
import json, os, numpy as np
ROOT = "/Users/paulsmith/space <3/tess-jev"
A = os.path.join(ROOT, "addenda/tic4206066_PUBLICATION_AUDIT"); U = os.path.join(ROOT, "addenda/opus55_final_audit_tic4206066")
OUT = os.path.join(ROOT, "addenda/tic4206066_FINAL_RECONCILIATION/_work/T2_A_shape/verify/out/v01_checks.json")
o = {}
s17 = json.load(open(os.path.join(A, "out/s17_detrend_robust.json")))
m = {}
for k, v in s17.items():
    ps = v["A"]["per_sector"]; d = np.array([x["depth"] for x in ps.values()]); e = np.array([x.get("err", x.get("depth_err")) for x in ps.values()])
    w = 1 / e ** 2; m[k] = [float(np.sum(w * d) / w.sum()), float(1 / np.sqrt(w.sum()))]
o["s17_allsector_wmean_from_per_sector"] = m
o["s17_allsector_stored"] = {k: {kk: vv for kk, vv in v["A"].items() if kk != "per_sector"} for k, v in s17.items()}
kex = {"U07 mann_claret_jit_raw_excl-1_8": 0.02062, "U07 mann_claret_jit_base_excl-1_8": 0.02034, "verify o-1 98e-1": 0.02014,
       "verify o1 98e-1": 0.02100, "verify o2 SYS 98e-1": 0.01999, "verify o2 DET 98e-1": 0.01991, "R-S6 excl98e-1": 0.02070}
kin = {"U07 mann_claret_jit_raw": 0.02095, "U07 mann_claret_jit_base": 0.01965, "U07 det98": 0.02084, "U07 b29": 0.02068,
       "U07 beta_base": 0.02001, "verify o-1 none": 0.02010, "verify o2 SYS none": 0.01952, "verify o2 DET none": 0.01953,
       "audit FINAL": 0.02116, "R-S6 raw all3": 0.020958, "R-S6 bw all3": 0.020139}
for lab, dct in (("k_edge_excluded", kex), ("k_all23_mann_claret", kin)):
    a = np.array(list(dct.values())); o[lab] = dict(values=dct, mean=float(a.mean()), sd=float(a.std(ddof=1)), min=float(a.min()), max=float(a.max()))
# paired epoch -1 effect on k (same code, with vs without)
o["k_effect_of_excluding_epoch-1"] = {"U07 raw": 0.02062 - 0.02095, "U07 biweight base": 0.02034 - 0.01965,
                                      "verify o-1": 0.02014 - 0.02010, "verify o2 SYS": 0.01999 - 0.01952, "verify o2 DET": 0.01991 - 0.01953,
                                      "R-S6 raw": 0.02070 - 0.020958}
pp = [(519.29, 70.38), (483.96, 104.23), (574.27, 41.89)]
d = np.array([x[0] for x in pp]); e = np.array([x[1] for x in pp]); w = 1 / e ** 2; mu = np.sum(w * d) / w.sum()
o["pubprep_chi2_weighted"] = float(np.sum(((d - mu) / e) ** 2)); o["pubprep_chi2_unweighted_mean"] = float(np.sum(((d - d.mean()) / e) ** 2))
s09 = json.load(open(os.path.join(A, "out/s09_events.json")))["A"]
pdc = {r["epoch"]: r for r in s09["pipe"]["rows"] if r["sector"] == 6 and r.get("depth_ppm") is not None}
sap = {r["epoch"]: r for r in s09["others"]["rows"] if r["dataset"] == "TESSSPOC_S6_SAP" and r.get("depth_ppm") is not None}
diff = np.array([sap[k]["depth_ppm"] - pdc[k]["depth_ppm"] for k in sorted(pdc) if k in sap])
o["S6_SAP_minus_PDCSAP_paired"] = dict(per_event=diff.tolist(), mean=float(diff.mean()), sem=float(diff.std(ddof=1) / np.sqrt(len(diff))), n=len(diff))
# dump distances, manuscript ephemeris
P, T0 = 3.182785, 2460998.736 - 2457000
RS = json.load(open(os.path.join(U, "referee/resolve_R-S6/out/R-S6_summary.json")))
dumps = np.array([x["t"] for x in RS["s6_checks"]["S6_dumps"]["oot_steps_at_dumps"]])
dist = {}
for n in (-795, -794, -793, -791, -790, -789):
    tc = T0 + n * P; dist[n] = float(np.min(np.abs(dumps - tc)) * 24)
o["S6_transit_to_nearest_dump_h_manuscript_ephem"] = dist
# no-S6 (R-S6 joint fits) vs the revised manuscript values
jf = RS["joint_fits"]; ms = dict(k=(0.021, 0.001), b=(0.35, 0.12), T14_h=(1.87, 0.05))
cmp = {}
for tag in ("S32_98_raw_jit_ss615", "S32_98_bw_jit_ss615", "S32_98_raw_jit_ss615_excl98e-1"):
    cmp[tag] = {q: float((jf[tag][q][1] - ms[q][0]) / ms[q][1]) for q in ms}
o["noS6_vs_revised_manuscript_nsigma"] = cmp
json.dump(o, open(OUT, "w"), indent=1)
print(json.dumps(o, indent=1))
