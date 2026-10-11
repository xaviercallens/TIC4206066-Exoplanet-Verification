"""R5: collect every FPP number used in the resolution into one JSON (out_r06_summary.json).
Sources: package runs (tic4206066_PUBLICATION_AUDIT/out/fpp/runs, rnaas pubprep out/fpp), audit 08 (out_fpp_summary.json),
verifier (08_fpp/verify/runs_v), this folder (r02-r05). Plain python3 (json + statistics only)."""
import json, os, glob, statistics as st
H = os.path.dirname(os.path.abspath(__file__))
ADD = os.path.abspath(os.path.join(H, "..", "..", ".."))
AUD08 = os.path.join(ADD, "opus55_final_audit_tic4206066", "08_fpp")
out = {}
pk = {}
for f in sorted(glob.glob(os.path.join(ADD, "tic4206066_PUBLICATION_AUDIT", "out", "fpp", "runs", "*", "summary.json"))):
    L = json.load(open(f))["log"]
    pk[os.path.basename(os.path.dirname(f))] = dict(mean=L["FPP_mean"], sd=L["FPP_std"], min=L["FPP_min"], max=L["FPP_max"], n=len(json.load(open(f))["runs"]),
                                                   Teff=L["star_params"]["Teff"], exptime_min=round(L["exptime_d"] * 1440, 2))
for s in "AB":
    j = json.load(open(os.path.join(ADD, "rnaas_tic4206066_pubprep", "analysis", "out", "fpp", f"fpp_ownfold_{s}.json")))
    pk[f"{s}_ownfold_pubprep"] = dict(mean=j["FPP_mean"], sd=j["FPP_sd"], min=j["FPP_min"], max=j["FPP_max"], n=len(j["runs"]), Teff=j["star"]["Teff"], exptime_min=6)
out["package_runs"] = pk
out["audit08_runs"] = {r["variant"]: dict(mean=r["FPP_mean"], sd=r["FPP_sd"], n=r["n_runs"], Teff=r["Teff"], exptime_min=r["exptime_min"])
                       for r in json.load(open(os.path.join(AUD08, "out_fpp_summary.json")))}
vr = {}
for f in sorted(glob.glob(os.path.join(AUD08, "verify", "runs_v", "*.json"))):
    j = json.load(open(f))
    if "runs" in j and j.get("FPP_mean") is not None:
        vr[os.path.basename(f)[:-5]] = dict(mean=j["FPP_mean"], sd=j.get("FPP_sd"), n=len(j["runs"]), Teff=j.get("Teff"), exptime_min=j.get("exptime_min"))
out["verifier_runs"] = vr
mine = {}
for f in sorted(glob.glob(os.path.join(H, "runs", "*.json"))):
    j = json.load(open(f)); fp = [r["FPP"] for r in j["runs"]]
    mine[os.path.basename(f)[:-5]] = dict(mean=st.mean(fp), sd=st.stdev(fp) if len(fp) > 1 else None, min=min(fp), max=max(fp), n=len(fp),
                                          seeds=[r["seed"] for r in j["runs"]], E_s=j.get("E_s"), NFPP_max=j.get("NFPP_max"))
out["r05_direct_adopted_fold"] = mine
# paired comparisons on identical seeds (default A/B from the package, same seeds 1000+i)
paired = {}
for s in "AB":
    base = {r["seed"]: r["FPP"] for r in json.load(open(os.path.join(ADD, "tic4206066_PUBLICATION_AUDIT", "out", "fpp", "runs", f"{s}_base_S32+S98", "summary.json")))["runs"]}
    for f in sorted(glob.glob(os.path.join(H, "runs", f"*_{s}_adoptedfold.json"))):
        j = json.load(open(f)); c = [(base[r["seed"]], r["FPP"]) for r in j["runs"] if r["seed"] in base]
        if c:
            ratios = [b / a for a, b in c]
            paired[os.path.basename(f)[:-5]] = dict(n=len(c), base_mean_same_seeds=st.mean(a for a, _ in c), mode_mean=st.mean(b for _, b in c),
                                                    ratio_mean=st.mean(ratios), ratio_sd=st.stdev(ratios) if len(ratios) > 1 else None)
out["paired_vs_package_default_same_seeds"] = paired
for k in ("out_r02_constrained_grid.json", "out_r03_constrained_molusc.json", "out_r04_molusc_decomposition.json"):
    p = os.path.join(H, k)
    if os.path.exists(p):
        out[k[:-5]] = json.load(open(p))
# ranges quoted in the recommended text
A_def = [pk["A_base_S32+S98"]["mean"], pk["A_repro_S32+S98"]["mean"], pk["A_ownfold_pubprep"]["mean"], out["audit08_runs"]["A_base"]["mean"],
         out["audit08_runs"]["A_10min"]["mean"], vr["default_A_v_A_S32+S98_5min"]["mean"], vr["default_A_v_A_S32+S98_10min"]["mean"]]
B_def = [pk["B_base_S32+S98"]["mean"], pk["B_prevfold_repro_S32+S98"]["mean"], pk["B_ownfold_pubprep"]["mean"], pk["B_tic_S32+S98"]["mean"],
         out["audit08_runs"]["B_base"]["mean"], out["audit08_runs"]["B_10min"]["mean"], vr["default_B_v_B_S32+S98_5min"]["mean"]]
out["ranges"] = dict(A_default_folds_and_bins=[min(A_def), max(A_def)], B_default_folds_bins_TICparams=[min(B_def), max(B_def)],
                     A_package_fold_means=[pk["A_repro_S32+S98"]["mean"], pk["A_base_S32+S98"]["mean"], pk["A_ownfold_pubprep"]["mean"]],
                     B_package_means_incl_TIC=[pk["B_tic_S32+S98"]["mean"], pk["B_prevfold_repro_S32+S98"]["mean"], pk["B_base_S32+S98"]["mean"], pk["B_ownfold_pubprep"]["mean"]])
json.dump(out, open(os.path.join(H, "out_r06_summary.json"), "w"), indent=1)
print(json.dumps({k: out[k] for k in ("r05_direct_adopted_fold", "paired_vs_package_default_same_seeds", "ranges")}, indent=1))
