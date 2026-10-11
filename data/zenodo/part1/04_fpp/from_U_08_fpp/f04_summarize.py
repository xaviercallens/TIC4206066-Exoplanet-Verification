"""Audit 8: aggregate runs/*/runs.json into out_fpp_summary.{json,csv} and a plot (audit_env or trice_env)."""
import json, os, glob, csv
import numpy as np
H = os.path.dirname(os.path.abspath(__file__))
rows = []; full = {}
def ms(x):
    x = np.array([v for v in x if v is not None and np.isfinite(v)], float)
    if len(x) == 0: return (None, None, None, None)
    return float(x.mean()), float(x.std(ddof=1)) if len(x) > 1 else None, float(x.min()), float(x.max())
for f in sorted(glob.glob(os.path.join(H, "runs", "*", "runs.json"))):
    d = json.load(open(f)); v = d["variant"]
    if v.startswith("test"): continue
    R = d["runs"]
    r = dict(variant=v, signal=d["signal"], fold=d["fold"], n_runs=len(R), N=d["args"]["N"], M=d["args"]["M"], R=d["args"]["R"], Teff=d["args"]["Teff"],
             ap=d["args"]["ap"], exptime_min=d["args"]["exptime_min"], sigma_ppm=round(d["sigma_used"] * 1e6, 1), molusc=os.path.basename(d["args"]["molusc"] or ""))
    for key, lab in (("FPP", "FPP"), ("NFPP", "NFPP")):
        r[lab + "_mean"], r[lab + "_sd"], r[lab + "_min"], r[lab + "_max"] = ms([x[key] for x in R])
    for s in ("TP", "PTP", "STP", "DTP", "SEBx2P", "BTP"):
        r["P_" + s] = ms([x["probs"].get(s, 0.0) for x in R])[0]
    for k in ("FPP_comp_mass_le_0.30", "FPP_comp_mass_le_0.38", "FPP_comp_mass_le_0.45", "FPP_comp_mass_le_0.54", "FPP_comp_mass_le_0.38_STP_e_lt_0.1", "FPP_STP_e_lt_0.1"):
        vals = [x.get(k) for x in R]
        r[k + "_mean"], r[k + "_sd"] = ms(vals)[:2]
    st = [x.get("STP_draws") or {} for x in R]
    for k in ("wfrac_mass_gt_0.38", "wfrac_mass_gt_0.54", "wfrac_ecc_lt_0.1", "wmean_mass", "wmean_ecc"):
        r["STP_" + k] = ms([s.get(k) for s in st])[0]
    mb = [s.get("wfrac_mass_bins") for s in st if s.get("wfrac_mass_bins")]
    if mb: r["STP_mass_bins"] = {k: float(np.mean([m[k] for m in mb])) for k in mb[0]}
    eb = [s.get("wfrac_ecc_bins") for s in st if s.get("wfrac_ecc_bins")]
    if eb: r["STP_ecc_bins"] = {k: float(np.mean([m[k] for m in eb])) for k in eb[0]}
    r["best_TP_Rp_mean"] = ms([x["best_TP"]["R_p"] for x in R])[0]
    r["sec_per_run"] = ms([x["seconds"] for x in R])[0]
    rows.append(r)
json.dump(rows, open(os.path.join(H, "out_fpp_summary.json"), "w"), indent=1)
keys = [k for k in rows[0] if not isinstance(rows[0][k], dict)]
with open(os.path.join(H, "out_fpp_summary.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore"); w.writeheader(); [w.writerow(r) for r in rows]
for r in rows:
    f = lambda x: "nan" if x is None else "%.4f" % x
    print("%-16s n=%2d FPP=%s±%s [%s-%s] NFPP=%.1e STP=%s PTP=%s | FPP(M<=.38)=%s FPP(M<=.45)=%s | STP w(M>.38)=%s w(e<.1)=%s" % (
        r["variant"], r["n_runs"], f(r["FPP_mean"]), f(r["FPP_sd"]), f(r["FPP_min"]), f(r["FPP_max"]), r["NFPP_mean"] or 0, f(r["P_STP"]), f(r["P_PTP"]),
        f(r["FPP_comp_mass_le_0.38_mean"]), f(r.get("FPP_comp_mass_le_0.45_mean")), f(r["STP_wfrac_mass_gt_0.38"]), f(r["STP_wfrac_ecc_lt_0.1"])))
