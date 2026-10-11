"""Audit 2, steps 2c-4: injection statistics, calibration of the error models, linear ephemeris under every error
model, Monte-Carlo chi2 / sigma_P under the empirical injection error distribution, leave-one-out, outliers,
effect of the dropped events.  Outputs out/t5_*.json/csv."""
import json
import numpy as np, pandas as pd
from scipy import stats
import t2_common as C

rng = np.random.default_rng(77)
ev = pd.read_csv(f"{C.OUT}/t3_events_base.csv")
meta = json.load(open(f"{C.OUT}/t3_meta_base.json"))
inj = {s: pd.read_csv(f"{C.OUT}/t4_inj_S{s}.csv") for s in (6, 32, 98)}
beta = {int(k): v["adopted"] for k, v in meta["betas"].items()}
res = {"betas": beta}


def rsig(x):
    x = np.asarray(x); return float(1.4826 * np.median(np.abs(x - np.median(x))))


def dstats(x, R=None):
    x = np.asarray(x, float)
    d = dict(n=int(len(x)), mean=float(np.mean(x)), sd=float(np.std(x, ddof=1)), robust_sigma=rsig(x),
             p2_5=float(np.percentile(x, 2.5)), p16=float(np.percentile(x, 16)), p50=float(np.median(x)),
             p84=float(np.percentile(x, 84)), p97_5=float(np.percentile(x, 97.5)),
             frac_abs_gt_30min=float(np.mean(np.abs(x) > 30)), frac_abs_gt_60min=float(np.mean(np.abs(x) > 60)),
             excess_kurtosis=float(stats.kurtosis(x)))
    rs = d["robust_sigma"]
    d["frac_gt_3robust"] = float(np.mean(np.abs(x - d["p50"]) > 3 * rs))
    d["frac_gt_3robust_gauss_expect"] = 0.0027
    if R is not None:
        d["frac_at_search_edge"] = float(np.mean(np.abs(x) >= R - 0.5))
    return d


# ------------------------------------------------------------------ 1. injection statistics
istat = {}
for s, d in inj.items():
    d["sel"] = (d.coverage >= 0.7) & (d.snr >= 2) & (~d.edge)
    d["sel_alt"] = (d.coverage >= 0.7) & (d.snr_alt >= 2) & (~d.edge_alt)
    istat[s] = dict(
        D_inj_ppm=float(d.D_inj.iloc[0] * 1e6), n_accepted_cov_ge_0p5=int(len(d)), n_selected=int(d.sel.sum()),
        frac_selected=float(d.sel.mean()),
        all_R144=dstats(d.dt, 144), selected_R144=dstats(d.dt[d.sel], 144),
        all_R60=dstats(d.dt_alt, 60), selected_R60=dstats(d.dt_alt[d.sel_alt], 60),
        cov_ge_0p7_R60_no_snr_cut=dstats(d.dt_alt[d.coverage >= 0.7], 60),
        rejected_by=dict(cov_lt_0p7=int((d.coverage < 0.7).sum()), snr_lt_2=int(((d.coverage >= 0.7) & (d.snr < 2)).sum()),
                         edge=int(((d.coverage >= 0.7) & (d.snr >= 2) & d.edge).sum())),
        median_e_white_sel=float(d.e_white[d.sel].median()), median_e_pb_sel=float(d.e_pb[d.sel].median()),
        median_snr_sel=float(d.snr[d.sel].median()))
res["injections"] = istat

# ------------------------------------------------------------------ 2. calibration of per-event sigma models
cal = {}
for s, d in inj.items():
    q = d[d.sel].copy()
    sd_const = q.dt.std(ddof=1)
    models = {"white": q.e_white, "beta_white": q.e_white * beta[s], "prayer_bead": q.e_pb, "bootstrap": q.e_boot,
              "const_sector_SD": pd.Series(sd_const, index=q.index), "const_sector_robust": pd.Series(rsig(q.dt), index=q.index)}
    cs = {}
    for k, sig in models.items():
        z = q.dt / sig
        cs[k] = dict(sd_z=float(np.std(z, ddof=1)), robust_sigma_z=rsig(z), frac_abs_z_gt2=float(np.mean(np.abs(z) > 2)),
                     frac_abs_z_gt3=float(np.mean(np.abs(z) > 3)), mean_z2=float(np.mean(z ** 2)),
                     median_z2=float(np.median(z ** 2)))
    # does |dt| scale with the per-event white sigma?  quartile bins of e_white
    bins = np.quantile(q.e_white, [0, .25, .5, .75, 1])
    bb = []
    for a, b in zip(bins[:-1], bins[1:]):
        m = (q.e_white >= a) & (q.e_white <= b)
        bb.append(dict(e_white_lo=float(a), e_white_hi=float(b), n=int(m.sum()), e_white_median=float(q.e_white[m].median()),
                       sd_dt=float(q.dt[m].std(ddof=1)), robust_dt=rsig(q.dt[m]), ratio_robust_to_white=float(rsig(q.dt[m]) / q.e_white[m].median()),
                       frac_abs_dt_gt_30=float(np.mean(np.abs(q.dt[m]) > 30))))
    # same in bins of snr and coverage
    bs = []
    for a, b in [(2, 3), (3, 4), (4, 5), (5, 20)]:
        m = (q.snr >= a) & (q.snr < b)
        if m.sum() > 10:
            bs.append(dict(snr_lo=a, snr_hi=b, n=int(m.sum()), sd_dt=float(q.dt[m].std(ddof=1)), robust_dt=rsig(q.dt[m]),
                           frac_abs_dt_gt_30=float(np.mean(np.abs(q.dt[m]) > 30))))
    bc = []
    for a, b in [(0.7, 0.9), (0.9, 2)]:
        m = (q.coverage >= a) & (q.coverage < b)
        if m.sum() > 10:
            bc.append(dict(cov_lo=a, cov_hi=b, n=int(m.sum()), sd_dt=float(q.dt[m].std(ddof=1)), robust_dt=rsig(q.dt[m])))
    # Gaussian log-likelihood comparison (per event sigma = a*e_white  vs  constant c), a and c by ML
    a_ml = np.sqrt(np.mean((q.dt / q.e_white) ** 2)); c_ml = np.sqrt(np.mean(q.dt ** 2))
    ll_scaled = float(np.sum(stats.norm.logpdf(q.dt, scale=a_ml * q.e_white)))
    ll_const = float(np.sum(stats.norm.logpdf(q.dt, scale=c_ml)))
    # Student-t (nu = 3) version, less dominated by the tails
    def tfit(scale_vec_fn):
        best = None
        for a in np.linspace(0.2, 10, 400):
            ll = np.sum(stats.t.logpdf(q.dt, df=3, scale=scale_vec_fn(a)))
            if best is None or ll > best[0]:
                best = (ll, a)
        return best
    t_scaled = tfit(lambda a: a * q.e_white); t_const = tfit(lambda a: a * 5.0)
    cal[s] = dict(z_stats=cs, by_e_white_quartile=bb, by_snr=bs, by_coverage=bc,
                  gauss_ML=dict(a_scaled=float(a_ml), ll_scaled=ll_scaled, c_const_min=float(c_ml), ll_const=ll_const),
                  student3_ML=dict(a_scaled=float(t_scaled[1]), ll_scaled=float(t_scaled[0]), c_const_min=float(t_const[1] * 5.0), ll_const=float(t_const[0])),
                  spearman_absdt_vs_ewhite=float(stats.spearmanr(np.abs(q.dt), q.e_white)[0]),
                  spearman_p=float(stats.spearmanr(np.abs(q.dt), q.e_white)[1]))
res["calibration"] = cal


# ------------------------------------------------------------------ 3. linear ephemeris under each error model
def linfit(n, tc, e):
    A = np.vstack([np.ones_like(n), n]).T; W = 1 / e ** 2
    cov = np.linalg.inv(A.T @ (A * W[:, None])); c = cov @ (A.T @ (W * tc))
    oc = tc - A @ c; chi2 = float(np.sum((oc / e) ** 2)); dof = len(n) - 2
    return dict(T0=float(c[0]), P=float(c[1]), eT0=float(np.sqrt(cov[0, 0])), eP=float(np.sqrt(cov[1, 1])),
                chi2=chi2, dof=dof, p_upper=float(stats.chi2.sf(chi2, dof)), p_lower=float(stats.chi2.cdf(chi2, dof)),
                eP_scaled=float(np.sqrt(cov[1, 1]) * np.sqrt(chi2 / dof)), oc_min=(oc * 1440).tolist(),
                rms_oc_min=float(np.sqrt(np.mean((oc * 1440) ** 2))))


S = ev[ev.sel].copy().reset_index(drop=True)
n = S.epoch.values.astype(float); tc = S.tc_obs.values
sd_sel = {s: istat[s]["selected_R144"]["sd"] for s in inj}
rob_sel = {s: istat[s]["selected_R144"]["robust_sigma"] for s in inj}
sd_all60 = {s: istat[s]["cov_ge_0p7_R60_no_snr_cut"]["sd"] for s in inj}
a_rob = {s: rsig(inj[s].dt[inj[s].sel] / inj[s].e_white[inj[s].sel]) for s in inj}
a_sd = {s: float(np.std(inj[s].dt[inj[s].sel] / inj[s].e_white[inj[s].sel], ddof=1)) for s in inj}
errmodels = {
    "white": S.e_white.values,
    "beta_white": S.e_beta.values,
    "prayer_bead": S.e_pb.values,
    "bootstrap": S.e_boot.values,
    "inj_SD_selected_R144 (per-sector const)": S.sector.map(sd_sel).values,
    "inj_robust_selected_R144 (per-sector const)": S.sector.map(rob_sel).values,
    "inj_SD_cov0.7_R60_no_snr_cut (manuscript-like)": S.sector.map(sd_all60).values,
    "inj_scaled_white_robust (a_s*e_white)": (S.e_white * S.sector.map(a_rob)).values,
    "inj_scaled_white_SD (a_s*e_white)": (S.e_white * S.sector.map(a_sd)).values,
}
eph = {}
for k, e in errmodels.items():
    eph[k] = linfit(n, tc, np.maximum(e, 0.25) / 1440)
    eph[k]["errors_min_by_event"] = [float(x) for x in e]
res["ephemeris_20"] = eph
res["events_20"] = S[["sector", "epoch", "tc_obs", "dt", "coverage", "snr", "D", "e_white", "e_beta", "e_pb", "e_boot"]].to_dict("records")
res["factors_a_robust"] = a_rob; res["factors_a_sd"] = a_sd

# ------------------------------------------------------------------ 4. Monte Carlo: chi2 distribution and true scatter of P
ADOPT = "inj_SD_selected_R144 (per-sector const)"
base = eph[ADOPT]
pools = {s: inj[s].dt[inj[s].sel].values for s in inj}
pools_ratio = {s: (inj[s].dt[inj[s].sel] / inj[s].e_white[inj[s].sel]).values for s in inj}
mc = {}
for label, mode in (("sector_pool_const_sigma", "pool"), ("ratio_pool_scaled_sigma", "ratio")):
    Ps, chis = [], []
    e_fit = S.sector.map(sd_sel).values if mode == "pool" else (S.e_white * S.sector.map(a_sd)).values
    for it in range(20000):
        if mode == "pool":
            err = np.array([rng.choice(pools[s]) for s in S.sector])
        else:
            err = np.array([rng.choice(pools_ratio[s]) for s in S.sector]) * S.e_white.values
        tsim = base["T0"] + base["P"] * n + err / 1440
        f = linfit(n, tsim, e_fit / 1440)
        Ps.append(f["P"] - base["P"]); chis.append(f["chi2"])
    Ps = np.array(Ps); chis = np.array(chis)
    obs = linfit(n, tc, e_fit / 1440)
    mc[label] = dict(sd_P=float(np.std(Ps)), robust_sd_P=rsig(Ps), formal_eP=float(obs["eP"]),
                     chi2_obs=obs["chi2"], chi2_mc_median=float(np.median(chis)), chi2_mc_mean=float(np.mean(chis)),
                     chi2_mc_p16=float(np.percentile(chis, 16)), chi2_mc_p84=float(np.percentile(chis, 84)),
                     frac_chi2_below_obs=float(np.mean(chis <= obs["chi2"])))
res["mc_empirical_errors"] = mc

# bootstrap over events (ephemeris level)
bP = []
for it in range(20000):
    idx = rng.integers(0, len(n), len(n))
    if len(np.unique(S.sector.values[idx])) < 2:
        continue
    bP.append(linfit(n[idx], tc[idx], S.sector.map(sd_sel).values[idx] / 1440)["P"])
res["event_bootstrap_sd_P"] = float(np.std(bP)); res["event_bootstrap_robust_sd_P"] = rsig(bP)

# ------------------------------------------------------------------ 5. leave-one-out, outliers, added events
e_ad = S.sector.map(sd_sel).values / 1440
loo = []
for i in range(len(n)):
    m = np.arange(len(n)) != i
    f = linfit(n[m], tc[m], e_ad[m])
    loo.append(dict(sector=int(S.sector[i]), epoch=int(S.epoch[i]), P=f["P"], dP_sigma=(f["P"] - base["P"]) / base["eP"],
                    chi2=f["chi2"], dof=f["dof"], oc_min=base["oc_min"][i], z=base["oc_min"][i] / (e_ad[i] * 1440)))
res["leave_one_out"] = loo
res["outliers_abs_z_gt_2"] = [r for r in loo if abs(r["z"]) > 2]
# add the events dropped by the selection (and the doc98 events)
alle = ev.copy()
add = {}
for label, rows in (("add_all_23_in_set", alle), ("add_cov_ge_0p7_only", alle[alle.coverage >= 0.7])):
    e = rows.sector.map(sd_sel).values / 1440
    add[label] = dict(n=int(len(rows)), **{k: v for k, v in linfit(rows.epoch.values.astype(float), rows.tc_obs.values, e).items() if k != "oc_min"})
d98 = pd.read_csv(f"{C.OUT}/t3_events_doc98.csv")
for label, rows in (("doc98_selected", d98[d98.sel]), ("doc98_all_in_set", d98)):
    e = rows.sector.map(sd_sel).values / 1440
    add[label] = dict(n=int(len(rows)), **{k: v for k, v in linfit(rows.epoch.values.astype(float), rows.tc_obs.values, e).items() if k != "oc_min"})
# sensitivity variants (same selection rule, adopted errors)
for v in ("ld_lo", "ld_hi", "rmB_manu", "fulldetrend", "W030"):
    d = pd.read_csv(f"{C.OUT}/t3_events_{v}.csv"); d = d[d.sel]
    e = d.sector.map(sd_sel).values / 1440
    add["variant_" + v] = dict(n=int(len(d)), **{k: vv for k, vv in linfit(d.epoch.values.astype(float), d.tc_obs.values, e).items() if k != "oc_min"})
# R = 60 min search (manuscript-like range) with R60 selection and R60 injection SDs
d = ev.copy(); d["sel60"] = (d.coverage >= 0.7) & (d.snr_alt >= 2) & (~d.edge_alt)
d = d[d.sel60]
sd60 = {s: istat[s]["selected_R60"]["sd"] for s in inj}
add["search_R60_sel60"] = dict(n=int(len(d)), epochs=[f"S{a}:{b}" for a, b in zip(d.sector, d.epoch)],
                               **{k: vv for k, vv in linfit(d.epoch.values.astype(float), d.tc_obs_alt.values, d.sector.map(sd60).values / 1440).items() if k != "oc_min"})
add["search_R60_sel60_white"] = {k: vv for k, vv in linfit(d.epoch.values.astype(float), d.tc_obs_alt.values, np.maximum(d.e_white_alt.values, .25) / 1440).items() if k != "oc_min"}
res["added_events_and_variants"] = add

json.dump(res, open(f"{C.OUT}/t5_results.json", "w"), indent=1, default=float)
# ------------------------------------------------------------------ console summary
for s in inj:
    I = istat[s]
    print(f"S{s}: Dinj {I['D_inj_ppm']:.0f} ppm  n={I['n_accepted_cov_ge_0p5']} sel={I['n_selected']}  "
          f"SD all144 {I['all_R144']['sd']:.1f} sel144 {I['selected_R144']['sd']:.1f} rob {I['selected_R144']['robust_sigma']:.1f} | "
          f"R60 cov0.7 nosnr SD {I['cov_ge_0p7_R60_no_snr_cut']['sd']:.1f} sel60 {I['selected_R60']['sd']:.1f} rob60 {I['selected_R60']['robust_sigma']:.1f} "
          f"edge144 {I['selected_R144']['frac_at_search_edge']:.3f} >60 {I['selected_R144']['frac_abs_gt_60min']:.3f} kurt {I['selected_R144']['excess_kurtosis']:.1f}")
    print("   z:", {k: (round(v['sd_z'], 2), round(v['robust_sigma_z'], 2), round(v['frac_abs_z_gt3'], 3)) for k, v in cal[s]["z_stats"].items()})
    print("   by e_white:", [(round(b['e_white_median'], 1), round(b['sd_dt'], 1), round(b['robust_dt'], 1)) for b in cal[s]["by_e_white_quartile"]])
    print("   by snr:", [(b['snr_lo'], b['n'], round(b['sd_dt'], 1), round(b['robust_dt'], 1)) for b in cal[s]["by_snr"]])
    print("   gaussML", cal[s]["gauss_ML"], "t3", cal[s]["student3_ML"])
for k, v in eph.items():
    print(f"{k:50s} P={v['P']:.8f}+-{v['eP']:.8f} T0={v['T0']:.5f}+-{v['eT0']:.5f} chi2={v['chi2']:.1f}/{v['dof']} p_up={v['p_upper']:.3g} p_low={v['p_lower']:.3g} rmsOC={v['rms_oc_min']:.1f}")
print("MC", json.dumps(mc, indent=0))
print("event bootstrap sd_P", res["event_bootstrap_sd_P"], res["event_bootstrap_robust_sd_P"])
for r in loo:
    print("LOO", r)
for k, v in add.items():
    print("ADD", k, v.get("n"), f"P={v['P']:.8f}+-{v['eP']:.8f} chi2={v['chi2']:.1f}/{v['dof']}", v.get("epochs", ""))
