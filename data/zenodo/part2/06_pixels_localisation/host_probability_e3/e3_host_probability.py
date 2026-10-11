"""Per-star host probability for signals A and B from the PRF difference-image fits of the opus55 audit (U/03_pixels).
Relative likelihood (as the RL of TESS-Localize, Higgins & Bell 2023): L_i ~ exp(-dchi2_i / (2 s)), dchi2_i = chi2 increase when
the PRF source is forced to star i (U/03_pixels/neighbours_dchi2.json), s = error-scale calibration per data set from on-target
injections: s = p95(inj dchi2_target) / 5.991 (p95 of chi2 with 2 dof). Sectors combined: S32 SPOC 2-min pointing-decorrelated +
S98 TESScut clean (the headline sets; S6 excluded: momentum dumps). Flat prior over the target and the listed neighbours
(all Gaia DR3 stars within 15" plus every star able to produce the depth in some product, U/03_pixels/neighbours_vs_fitted_position.csv).
Variant 'depth-feasible': stars that would need an eclipse deeper than 100 % in every product (dreq_anyproduct > 1) get zero prior.
Caveat: a companion unresolved from the target (<~1-2") is indistinguishable from the target here."""
import json, os, numpy as np, csv
U3 = "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels"
D = json.load(open(os.path.join(U3, "neighbours_dchi2.json")))
nb = {int(r["source_id"]): r for r in csv.DictReader(open(os.path.join(U3, "neighbours_vs_fitted_position.csv")))}
out = {}
for sig, dreq_key in (("A", "dreqA_anyproduct"), ("B", "dreqB_min")):
    sets = [f"S32_tpf_PC_{sig}", f"S98_clean29_{sig}"]
    logL = {"target": 0.0}; info = {}
    for s in sets:
        d = D[s]; scale = max(1.0, d["inj_dchi2_target_p95"] / 5.991); info[s] = dict(scale=scale, dchi2_target=d["dchi2_target"])
        logL["target"] += -d["dchi2_target"] / (2 * scale)
        for n in d["neighbours"]:
            k = str(n["source_id"]); logL.setdefault(k, 0.0); logL[k] += -n["dchi2"] / (2 * scale)
    ids = list(logL); lv = np.array([logL[i] for i in ids]); p = np.exp(lv - lv.max()); p /= p.sum()
    feas = np.array([1.0 if i == "target" else float(float(nb[int(i)][dreq_key]) <= 1.0) for i in ids])
    pf = p * feas; pf /= pf.sum()
    rows = []
    for i, pp, ppf in zip(ids, p, pf):
        r = nb.get(int(i)) if i != "target" else None
        rows.append(dict(star=i, TIC=(r["TIC"] if r else "4206066"), Tmag=(float(r["Tmag"]) if r else 9.68),
                         sep_arcsec=(round(float(r["sep"]), 1) if r else 0.0),
                         required_eclipse_fraction=(round(float(r[dreq_key]), 3) if r else None),
                         P_position_only=float(pp), P_depth_feasible=float(ppf)))
    rows.sort(key=lambda x: -x["P_position_only"])
    out[sig] = dict(sets=info, stars=rows)
    print(f"--- signal {sig}")
    for r in rows: print(f"  {r['star']:>20} TIC {r['TIC']:>11} T={r['Tmag']:5.2f} sep={r['sep_arcsec']:5.1f}\" dreq={r['required_eclipse_fraction']}  P_pos={r['P_position_only']:.2e}  P_pos+depth={r['P_depth_feasible']:.2e}")
json.dump(out, open("e3_host_probability.json", "w"), indent=1)
