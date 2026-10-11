"""Delta-chi2 (PRF fit of the mean difference image) at the target vs at each Gaia neighbour able to produce the
depth, for the headline data sets; calibrated against the on-target injection distribution of Delta-chi2(target)."""
import pandas as pd
from run_real import *
from poscorr_decorr import decorrelate
S = json.load(open(os.path.join(OUT, "summary.json")))["sets"]
nb = pd.read_csv(os.path.join(OUT, "neighbours_vs_fitted_position.csv"))
capA = set(nb[nb["dreqA_anyproduct"] <= 1]["source_id"]) | set(nb[nb["sep"] < 15]["source_id"])
capB = set(nb[nb["dreqB_min"] <= 1.2]["source_id"]) | set(nb[nb["sep"] < 15]["source_id"])
res = {}
for cfg, dec, injkey in [("S32_tpf", True, "S32_tpf_poscorr"), ("S98_clean29", False, "S98_clean29"),
                         ("S6_tspoc", True, "S6_tspoc_poscorr")]:
    name, key, var, kw = [c for c in CONFIGS if c[0] == cfg][0]
    ds = Dataset(key, var, **kw)
    if dec:
        decorrelate(ds)
    sf = ds.scene_fit()
    for sig, cap in (("A", capA), ("B", capB)):
        out, Dm, Ve, _ = analyse(ds, sig, sf, want_nb=True)
        inj = S.get(f"{injkey}_{sig}")
        rows = [dict(source_id=n["source_id"], dE=n["dE"], dN=n["dN"], Tmag=n["Tmag"], dchi2=n["dchi2"],
                     dchi2_minus_target=n["dchi2"] - out["dchi2_target"]) for n in out["neighbours"] if n["source_id"] in cap]
        res[f"{cfg}{'_PC' if dec else ''}_{sig}"] = dict(dchi2_target=out["dchi2_target"], best_fit_chi2=out["fit"]["chi2"],
                                                         dof=out["fit"]["dof"], inj_dchi2_target_p95=inj and inj["inj_dchi2_target_p95"],
                                                         neighbours=sorted(rows, key=lambda r: np.hypot(r["dE"], r["dN"])))
        print(cfg, "PC" if dec else "", sig, f"dchi2(target)={out['dchi2_target']:.1f} (on-target injections p95={inj and round(inj['inj_dchi2_target_p95'],1)})",
              " ".join(f"[{np.hypot(r['dE'],r['dN']):.0f}\" T{r['Tmag']:.1f}: {r['dchi2']:.0f}]" for r in sorted(rows, key=lambda r: np.hypot(r['dE'], r['dN']))))
# summed over sectors
for sig in "AB":
    keys = [k for k in res if k.endswith("_" + sig) and not k.startswith("S6")]
    tot = {}
    for k in keys:
        for r in res[k]["neighbours"]:
            tot.setdefault(r["source_id"], 0.0)
            tot[r["source_id"]] += r["dchi2_minus_target"]
    res[f"sum_S32PC_S98_{sig}"] = tot
    print(sig, "summed dchi2(neighbour)-dchi2(target), S32 PC + S98:", {k: round(v, 0) for k, v in tot.items()})
dump(res, "neighbours_dchi2.json")
