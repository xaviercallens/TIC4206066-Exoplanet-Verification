"""T2_A_shape r01: recompute per-sector and combined template depths of signal A from the STORED per-event rows of
pubprep (a1_per_transit.csv) and audit (s09_events.json), with and without S98 epoch -1 (segment-edge transit),
and chi2 of a constant across sectors for every candidate per-sector set. Read-only inputs; output to OUT/_work."""
import csv, json, os
import numpy as np
ROOT = "/Users/paulsmith/space <3/tess-jev"
W = os.path.join(ROOT, "addenda/tic4206066_FINAL_RECONCILIATION/_work/T2_A_shape")
PP = os.path.join(ROOT, "addenda/rnaas_tic4206066_pubprep/analysis/out/a1_per_transit.csv")
AU = os.path.join(ROOT, "addenda/tic4206066_PUBLICATION_AUDIT/out/s09_events.json")

def wmean(d, e):
    d, e = np.asarray(d, float), np.asarray(e, float); w = 1 / e ** 2
    return float(np.sum(w * d) / w.sum()), float(1 / np.sqrt(w.sum()))

def chi2_const(vals):
    d = np.array([v[0] for v in vals]); e = np.array([v[1] for v in vals]); mu, _ = wmean(d, e)
    return float(np.sum(((d - mu) / e) ** 2)), len(d) - 1, mu

out = {}
# pubprep: rows A,A,pipe ; D = free-dt best-fit depth ; per-sector error = injection SD (same for every event of a sector)
rows = [r for r in csv.DictReader(open(PP)) if r["signal"] == "A" and r["template"] == "A" and r["photometry"] == "pipe"]
for excl in (None, (98, -1)):
    per = {}
    for s in (6, 32, 98):
        R = [r for r in rows if int(r["sector"]) == s and not (excl and int(r["sector"]) == excl[0] and int(r["epoch"]) == excl[1])]
        d = [float(r["D"]) * 1e6 for r in R]; e = [float(r["e_D_inj"]) * 1e6 for r in R]
        per[f"S{s}"] = list(wmean(d, e)) + [len(R)]
    c2, dof, mu = chi2_const([per[k] for k in per])
    comb = wmean([per[k][0] for k in per], [per[k][1] for k in per])
    out[f"pubprep_inj{'_excl98e-1' if excl else ''}"] = dict(per_sector=per, combined=comb, chi2_const=c2, dof=dof)
# pubprep: odd/even from the same rows
for excl in (None, (98, -1)):
    oe = {}
    for name, par in (("odd", 1), ("even", 0)):
        R = [r for r in rows if int(r["epoch"]) % 2 == par and not (excl and int(r["sector"]) == excl[0] and int(r["epoch"]) == excl[1])]
        oe[name] = list(wmean([float(r["D"]) * 1e6 for r in R], [float(r["e_D_inj"]) * 1e6 for r in R])) + [len(R)]
    oe["z"] = (oe["odd"][0] - oe["even"][0]) / np.hypot(oe["odd"][1], oe["even"][1])
    out[f"pubprep_oddeven_inj{'_excl98e-1' if excl else ''}"] = oe
# audit s09: fixed-tc depth, covariance error x beta
A = json.load(open(AU))["A"]["pipe"]["rows"]
A = [r for r in A if r.get("depth_ppm") is not None]
for excl in (None, (98, -1)):
    per = {}
    for s in (6, 32, 98):
        R = [r for r in A if r["sector"] == s and not (excl and r["sector"] == excl[0] and r["epoch"] == excl[1])]
        per[f"S{s}"] = list(wmean([r["depth_ppm"] for r in R], [r["depth_err_ppm"] for r in R])) + [len(R)]
    c2, dof, mu = chi2_const([per[k] for k in per])
    allw = wmean([r["depth_ppm"] for r in A if not (excl and r["sector"] == excl[0] and r["epoch"] == excl[1])],
                 [r["depth_err_ppm"] for r in A if not (excl and r["sector"] == excl[0] and r["epoch"] == excl[1])])
    out[f"audit_s09_beta{'_excl98e-1' if excl else ''}"] = dict(per_sector=per, combined_all_events=allw, chi2_const=c2, dof=dof)
    oe = {}
    for name, par in (("odd", 1), ("even", 0)):
        R = [r for r in A if r["epoch"] % 2 == par and not (excl and r["sector"] == excl[0] and r["epoch"] == excl[1])]
        oe[name] = list(wmean([r["depth_ppm"] for r in R], [r["depth_err_ppm"] for r in R])) + [len(R)]
    oe["z"] = (oe["odd"][0] - oe["even"][0]) / np.hypot(oe["odd"][1], oe["even"][1])
    out[f"audit_s09_oddeven_beta{'_excl98e-1' if excl else ''}"] = oe
# chi2 of a constant for the other published per-sector sets (values copied with their source keys)
sets = {
 "U07_raw_localbaseline_emp (07_A_shape_radii/out/depths_localbaseline_mann_claret_jit_raw.json S6/S32/S98 depth,sig_emp)": [(522.92, 51.58), (475.30, 77.19), (513.62, 44.05)],
 "U07_biweight_emp (07_A_shape_radii/out/tests_results.json depths.S*.depth_ppm,sig_emp_std_ppm)": [(536.00, 55.10), (461.19, 73.94), (403.05, 52.50)],
 "U07verify_localquad_emp (07_A_shape_radii/verify/out/tests_ms_template.json S*.depth_ppm,sig_emp_std_ppm)": [(513.71, 60.67), (444.99, 73.86), (468.38, 46.27)],
 "R-S6_raw_linear_emp (referee/resolve_R-S6/out/R-S6_summary.json s6_checks.depths.S*.D,eD_used)": [(523.31, 58.17), (480.48, 76.45), (516.70, 46.99)],
}
for k, v in sets.items():
    c2, dof, mu = chi2_const(v); out["chi2_" + k] = dict(chi2=c2, dof=dof, weighted_mean=mu, combined=wmean([a for a, b in v], [b for a, b in v]))
json.dump(out, open(os.path.join(W, "out/r01_depths.json"), "w"), indent=1)
for k, v in out.items(): print(k, json.dumps(v))
