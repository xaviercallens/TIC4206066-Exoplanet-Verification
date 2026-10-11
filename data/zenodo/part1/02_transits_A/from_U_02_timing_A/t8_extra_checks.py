"""Audit 2: extra checks.
(1) S6 injections near momentum dumps (all real S6 transits are 0.4-3.5 h from a dump) vs elsewhere.
(2) R60 (manuscript-like) configuration: Monte-Carlo chi2 distribution and true sigma_P with the empirical
    injection errors; where does chi2 = 12.6 (manuscript) / 11.6 (this audit) fall?
(3) Uncertainty of the per-sector injection SD (bootstrap over injections; subsampling to 150 injections as in
    the manuscript's analysis) -> how much the adopted sigma_P can move from the SD estimate alone.
(4) Sensitivity of sigma_P to the per-sector SD (manuscript SDs 13.8/24.6/10.1 min applied to our event set)."""
import json
import numpy as np, pandas as pd
from scipy import stats
import t2_common as C

rng = np.random.default_rng(8)
ev = pd.read_csv(f"{C.OUT}/t3_events_base.csv")
t5 = json.load(open(f"{C.OUT}/t5_results.json"))
inj = {s: pd.read_csv(f"{C.OUT}/t4_inj_S{s}.csv") for s in (6, 32, 98)}
for s, d in inj.items():
    d["sel"] = (d.coverage >= 0.7) & (d.snr >= 2) & (~d.edge)
    d["sel60"] = (d.coverage >= 0.7) & (d.snr_alt >= 2) & (~d.edge_alt)
out = {}
rs = lambda x: float(1.4826 * np.median(np.abs(np.asarray(x) - np.median(x))))

# (1) dumps from 00_data README (TESS-SPOC/TESScut QUALITY bit 6), BTJD
dumps = np.array([1468.3912, 1471.5161, 1474.6411, 1481.2868, 1484.4117, 1487.5366])
d = inj[6]; q = d[d.sel].copy()
q["ddump_h"] = np.min(np.abs(q.t_inj.values[:, None] - dumps[None, :]), axis=1) * 24
near_ = q.ddump_h < 4.0
out["S6_dump_proximity"] = dict(n_near=int(near_.sum()), sd_near=float(q.dt[near_].std()), rob_near=rs(q.dt[near_]),
                                n_far=int((~near_).sum()), sd_far=float(q.dt[~near_].std()), rob_far=rs(q.dt[~near_]),
                                real_event_dump_distance_h=[float(np.min(np.abs(t - dumps)) * 24) for t in ev[ev.sector == 6].tc_pred])

# (2) R60 MC
rows = ev[(ev.coverage >= 0.7) & (ev.snr_alt >= 2) & (~ev.edge_alt)]
n = rows.epoch.values.astype(float)
sd60 = {s: float(inj[s].dt_alt[inj[s].sel60].std(ddof=1)) for s in inj}
pools = {s: inj[s].dt_alt[inj[s].sel60].values for s in inj}
e = rows.sector.map(sd60).values / 1440
A = np.vstack([np.ones_like(n), n]).T
def lf(y, e):
    W = 1 / e ** 2; cov = np.linalg.inv(A.T @ (A * W[:, None])); c = cov @ (A.T @ (W * y))
    return c, float(np.sum(((y - A @ c) / e) ** 2)), np.sqrt(np.diag(cov))
c_obs, chi_obs, se = lf(rows.tc_obs_alt.values, e)
Ps, chis = [], []
for k in range(20000):
    err = np.array([rng.choice(pools[s]) for s in rows.sector]) / 1440
    c, ch, _ = lf(c_obs[0] + c_obs[1] * n + err, e)
    Ps.append(c[1] - c_obs[1]); chis.append(ch)
Ps, chis = np.array(Ps), np.array(chis)
out["R60_MC"] = dict(sd60=sd60, P=float(c_obs[1]), formal_eP=float(se[1]), T0=float(c_obs[0]), eT0=float(se[0]),
                     chi2_obs=chi_obs, mc_sd_P=float(Ps.std()), mc_robust_sd_P=rs(Ps),
                     chi2_mc_median=float(np.median(chis)), chi2_mc_p5=float(np.percentile(chis, 5)),
                     chi2_mc_p95=float(np.percentile(chis, 95)), frac_mc_chi2_le_obs=float(np.mean(chis <= chi_obs)),
                     frac_mc_chi2_le_12p6=float(np.mean(chis <= 12.6)),
                     gauss_p_lower_12p6_18=float(stats.chi2.cdf(12.6, 18)), gauss_p_lower_obs=float(stats.chi2.cdf(chi_obs, 18)))

# (3) SD uncertainty: bootstrap of the injection sample, and spread of the SD from random 150-injection subsets
sdu = {}
for s, d in inj.items():
    x = d.dt_alt[d.coverage >= 0.7].values          # manuscript-like: coverage cut only, R60
    bs = [np.std(rng.choice(x, len(x)), ddof=1) for _ in range(2000)]
    sub = [np.std(rng.choice(x, 150, replace=False), ddof=1) for _ in range(2000)]
    sdu[s] = dict(sd=float(np.std(x, ddof=1)), boot_se=float(np.std(bs)), sub150_p16=float(np.percentile(sub, 16)),
                  sub150_p50=float(np.median(sub)), sub150_p84=float(np.percentile(sub, 84)),
                  sub150_p2_5=float(np.percentile(sub, 2.5)), sub150_p97_5=float(np.percentile(sub, 97.5)))
out["SD_uncertainty_R60_cov0p7"] = sdu

# (4) manuscript's per-sector injection SDs applied to this audit's R60 event set
ms = {6: 13.834, 32: 24.627, 98: 10.066}   # values quoted in the package (a1_results via MANUSCRIPT_NUMBER_AUDIT/a1 csv): only used for this what-if
c, ch, se2 = lf(rows.tc_obs_alt.values, rows.sector.map(ms).values / 1440)
out["R60_with_package_SDs"] = dict(P=float(c[1]), eP=float(se2[1]), T0=float(c[0]), eT0=float(se2[0]), chi2=ch)
json.dump(out, open(f"{C.OUT}/t8_results.json", "w"), indent=1)
print(json.dumps(out, indent=1))
