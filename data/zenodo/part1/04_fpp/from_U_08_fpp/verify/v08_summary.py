"""Verifier: collect all verify outputs into v_summary.json (numbers only)."""
import glob, json, os
import numpy as np
H = os.path.dirname(os.path.abspath(__file__))
S = {"direct_runs": {}, "grids": {}}
for f in sorted(glob.glob(os.path.join(H, "runs_v", "*.json"))):
    d = json.load(open(f)); k = os.path.basename(f)[:-5]
    if "runs" in d:
        fp = np.array([r["FPP"] for r in d["runs"]])
        S["direct_runs"][k] = dict(n=len(fp), FPP_mean=float(fp.mean()), FPP_sd=float(fp.std(ddof=1)) if len(fp) > 1 else None,
                                   FPP_min=float(fp.min()), FPP_max=float(fp.max()),
                                   NFPP_max=float(max(r["NFPP"] for r in d["runs"])))
for f in sorted(glob.glob(os.path.join(H, "v06_combine*.json"))):
    d = json.load(open(f))
    S["grids"][os.path.basename(f)[:-5]] = {s: {k: v for k, v in r.items() if k != "grid"} for s, r in d.items()}
S["molusc_dist"] = json.load(open(os.path.join(H, "v05_molusc_dist.json")))
S["folds"] = json.load(open(os.path.join(H, "folds", "folds_meta.json")))
json.dump(S, open(os.path.join(H, "v_summary.json"), "w"), indent=1)
for k, v in S["direct_runs"].items():
    print("%-45s n=%d FPP=%.4f+-%s [%.4f-%.4f] NFPPmax=%.1e" % (k, v["n"], v["FPP_mean"], "%.4f" % v["FPP_sd"] if v["FPP_sd"] else "-", v["FPP_min"], v["FPP_max"], v["NFPP_max"]))
