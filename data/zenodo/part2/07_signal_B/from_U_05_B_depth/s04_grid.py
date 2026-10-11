"""Sensitivity grid: per-sector B depths for every configuration in s_configs.CONFIGS.

For each configuration: per-sector depth D (fixed-shape template amplitude = mid-transit depth), white error,
injection-scaled error (kappa per sector from s03 injections of the SAME configuration if available, else from BASE),
S98-minus-archival difference in sigma, common depth and delta-chi2 (common vs per-sector, 2 dof).
Output: out/s04_grid.json, out/s04_grid.csv
"""
import csv, glob, sys
from scipy import stats
from b_run import *
from s_configs import CONFIGS, template_for


def load_kappa(name, mode="sector"):
    fn = os.path.join(OUT, f"s03_inj_{name}_{mode}.json")
    if not os.path.exists(fn):
        return None
    d = json.load(open(fn))
    return {int(s): (v["kappa"], v["kappa_bs_sd"], v["mean_bias_ppm"]) for s, v in d["sectors"].items()}


def evaluate(name, cfg, kap_base):
    tpl = template_for(cfg)
    out = run_config(cfg, tpl)
    kap_own = load_kappa(name)
    kap = kap_own or kap_base
    row = dict(config=name, kappa_source=(name if kap_own else "BASE") if kap else "none", template_T14_h=tpl.t14_d * 24,
               template=tpl.shape)
    dep_w, dep_k = {}, {}
    for s in (6, 32, 98):
        r = out[s]
        k = max(kap[s][0], 1.0) if kap else 1.0
        dep_w[s] = (r["D"], r["eD_white"]); dep_k[s] = (r["D"], r["eD_white"] * k)
        row[f"D{s}"] = r["D"]; row[f"eW{s}"] = r["eD_white"]; row[f"eK{s}"] = r["eD_white"] * k; row[f"kappa{s}"] = k
        row[f"chi2red{s}"] = r["chi2_red"]; row[f"events{s}"] = r["epochs"]; row[f"cov{s}"] = r["cov"]
        row[f"Dev{s}"] = r["D_event"]; row[f"eWev{s}"] = r["eD_event_white"]; row[f"prod{s}"] = r["product"]
        row[f"dt{s}"] = r["dt_min"]
    for tag, dep in (("W", dep_w), ("K", dep_k)):
        sm = s98_vs_archival(dep)
        row[f"arch_{tag}"] = sm["archival"][0]; row[f"earch_{tag}"] = sm["archival"][1]
        row[f"diff_{tag}"] = sm["diff"]; row[f"nsig_{tag}"] = sm["nsigma"]
        row[f"common_{tag}"] = sm["common"][0]; row[f"ecommon_{tag}"] = sm["common"][1]
        row[f"dchi2_{tag}"] = sm["chi2_common_3sectors"]; row[f"p_common_{tag}"] = float(stats.chi2.sf(sm["chi2_common_3sectors"], 2))
    return row


if __name__ == "__main__":
    kap_base = load_kappa("BASE")
    only = sys.argv[1:] if len(sys.argv) > 1 else None
    rows = []
    for name, cfg in CONFIGS.items():
        if only and name not in only:
            continue
        try:
            row = evaluate(name, cfg, kap_base)
        except Exception as e:
            print("FAILED", name, e, flush=True); continue
        rows.append(row)
        print(f"{name:16s} S6 {row['D6']:5.0f}±{row['eK6']:3.0f}  S32 {row['D32']:5.0f}±{row['eK32']:3.0f}  S98 {row['D98']:5.0f}±{row['eK98']:3.0f} | "
              f"arch {row['arch_K']:5.0f}±{row['earch_K']:3.0f} nsig W {row['nsig_W']:.2f} K {row['nsig_K']:.2f} | common {row['common_K']:.0f}±{row['ecommon_K']:.0f} "
              f"dchi2 W {row['dchi2_W']:.1f} K {row['dchi2_K']:.1f} (p {row['p_common_K']:.3f}) ev98 {np.round(row['Dev98']).astype(int).tolist()} kap {row['kappa_source']}", flush=True)
    tag = "" if not only else "_" + "_".join(only)[:60]
    jdump(rows, f"s04_grid{tag}.json")
    keys = ["config", "kappa_source", "D6", "eW6", "eK6", "D32", "eW32", "eK32", "D98", "eW98", "eK98", "arch_W", "earch_W", "nsig_W",
            "arch_K", "earch_K", "nsig_K", "common_K", "ecommon_K", "dchi2_W", "dchi2_K", "p_common_K", "prod6", "prod32", "prod98",
            "events6", "events32", "events98", "cov98", "Dev98", "template_T14_h"]
    with open(os.path.join(OUT, f"s04_grid{tag}.csv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore"); w.writeheader()
        for r in rows:
            w.writerow({k: (round(r[k], 3) if isinstance(r[k], float) else r[k]) for k in keys})
