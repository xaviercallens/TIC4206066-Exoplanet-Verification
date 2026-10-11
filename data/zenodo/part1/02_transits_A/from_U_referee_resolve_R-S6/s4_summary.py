"""R-S6: compile the key numbers of s1 (timing), s1b (aliases), s2 (joint fits) and s3 (S6 checks) into out/R-S6_summary.json."""
import json, glob
import numpy as np
from rs6_common import OUT, os

S = {}
t = json.load(open(os.path.join(OUT, "s1_timing_R60.json")))
eph = t["ephemerides"]
S["timing_R60"] = {k: {kk: eph[k][kk] for kk in ("n", "P", "eP", "T0", "eT0", "chi2", "dof")} for k in eph}
S["timing_R60_events"] = [{k: e[k] for k in ("sector", "epoch", "coverage", "snr", "dt", "sel", "n_pre", "n_post")} for e in t["events"]]
S["timing_R60_injections"] = t["injections"]
S["timing_R60_compare"] = t["compare"]
S["cycle_count"] = t["cycle_count"]
S["S6_OC_vs_S32S98_eph"] = t["S6_vs_S32S98"]
p144 = os.path.join(OUT, "s1_timing_R144.json")
if os.path.exists(p144):
    t2 = json.load(open(p144))
    S["timing_R144"] = {k: {kk: t2["ephemerides"][k][kk] for kk in ("n", "P", "eP", "T0", "eT0", "chi2", "dof")} for k in t2["ephemerides"]}
    S["timing_R144_compare"] = t2["compare"]
pa = os.path.join(OUT, "s1b_alias.json")
if os.path.exists(pa):
    S["alias_lightcurve_test"] = json.load(open(pa))["aliases"]
pm = os.path.join(OUT, "s1c_cycle_mc.json")
if os.path.exists(pm):
    S["cycle_count_mc"] = json.load(open(pm))
fits = {}
for f in sorted(glob.glob(os.path.join(OUT, "fit_*.json"))):
    d = json.load(open(f))
    row = dict(k=d["k"][1:4], k_mean_std=d["k_mean_std"], b=d["b"][1:4], T14_h=d["T14_h"][1:4], Rp_Re=d["Rp_Re"][1:4],
               Rp_mean_std=d["Rp_mean_std"], P=d["P"][1:4], T0=d["T0"][1:4], rho_q97_7=d["rho_q97_7"], n_over_tau=d["n_over_tau"],
               npts=d["npts"], supersample=d["supersample"], betas={k: v["beta"] for k, v in d["betas"].items()})
    for kk in ("k6", "dt6_min", "k6_minus_k", "z_k6", "z_dt6", "depth_ratio_k6sq_over_ksq"):
        if kk in d:
            row[kk] = d[kk]
    fits[d["tag"]] = row
S["joint_fits"] = fits
p3 = os.path.join(OUT, "s3_s6_checks.json")
if os.path.exists(p3):
    c = json.load(open(p3))
    S["s6_checks"] = {k: v for k, v in c.items() if k not in ("dumps_S6",)}
    if "S6_timing" in S["s6_checks"]:
        S["s6_checks"]["S6_timing"].pop("chi2_curve", None)
json.dump(S, open(os.path.join(OUT, "R-S6_summary.json"), "w"), indent=1, default=float)
print(json.dumps(S["joint_fits"], indent=0, default=float)[:6000])
