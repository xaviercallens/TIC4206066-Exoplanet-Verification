"""Summarise injection-calibrated localisation, combine sectors, compare with neighbours (audit 03_pixels)."""
import glob
import pandas as pd
from common import *
from scipy.stats import chi2 as chi2dist, norm

real = json.load(open(os.path.join(OUT, "real_results.json")))


def load_inj(job):
    p = os.path.join(OUT, "inj", f"inj_{job}.json")
    if not os.path.exists(p):
        return None
    r = json.load(open(p))
    if r.get("partial", True):
        return None
    return r


def robust_cov(x):
    """Covariance from the MAD-scaled spread along the principal axes plus the plain covariance."""
    x = np.asarray(x)
    C = np.cov(x.T)
    med = np.median(x, axis=0)
    # robust: scale plain covariance so that the 68.3% Mahalanobis radius matches chi2(2) expectation
    d2 = np.einsum("ij,jk,ik->i", x - med, np.linalg.inv(C), x - med)
    q = np.quantile(d2, 0.683) / chi2dist.ppf(0.683, 2)
    return C, C * q, med


def sig_equiv(d2, dof=2):
    p = chi2dist.sf(d2, dof)
    return p, float(norm.isf(p / 2)) if p > 0 else np.inf


sets = {}
for job in ["S32_tpf_A", "S32_ffi_A", "S98_clean_A", "S98_all_A", "S6_ffi_A", "S6_step_A", "S6_dumpx0.5h_A",
            "S32_tpf_deg2_A", "S98_clean_deg2_A", "S98_clean29_A", "S98_clean29_deg2_A", "S98_clean29_W3_A",
            "S98_clean29_W6_A", "S32_tpf_W3_A", "S32_tpf_W6_A", "S6_dumpx1h_A",
            "S32_tpf_B", "S32_ffi_B", "S98_clean_B", "S98_all_B", "S6_ffi_B", "S98_clean29_B",
            "S32_tpf_poscorr_A", "S32_tpf_poscorr_B", "S32_tspoc_poscorr_A", "S32_tspoc_poscorr_B",
            "S6_tspoc_poscorr_A", "S6_tspoc_poscorr_B", "S6_tspoc_dumpx0.5h_poscorr_A"]:
    r = load_inj(job)
    if r is None:
        print("missing", job); continue
    tr = pd.DataFrame(r["trials"])
    X = tr[["off_E", "off_N"]].values
    C, Crob, med = robust_cov(X)
    mu = X.mean(0)
    xr = np.array([r["real"]["off_E"], r["real"]["off_N"]])
    # Mahalanobis distance of the real offset from the on-target injection distribution
    d2_inj = float((xr - mu) @ np.linalg.solve(C, xr - mu))
    d2_all = np.einsum("ij,jk,ik->i", X - mu, np.linalg.inv(C), X - mu)
    p_emp = float(np.mean(d2_all >= d2_inj))
    d2_zero = float(xr @ np.linalg.solve(C, xr))
    dchi = tr["dchi2_target"].values
    sets[job] = dict(n_trials=len(tr), n_events=r["n_events"], depth_injected_ppm=1e6 * r["depth_injected"],
                     depth_source=r["depth_source"],
                     real_off=xr, real_dchi2_target=r["real"]["dchi2_target"], real_snr=r["real"]["snr_target"],
                     inj_mean=mu, inj_median=med, inj_std=np.sqrt(np.diag(C)), inj_corr=C[0, 1] / np.sqrt(C[0, 0] * C[1, 1]),
                     inj_cov=C, inj_cov_robust=Crob, inj_std_robust=np.sqrt(np.diag(Crob)),
                     inj_p16_p84_E=np.percentile(X[:, 0], [16, 84]), inj_p16_p84_N=np.percentile(X[:, 1], [16, 84]),
                     d2_real_vs_injdist=d2_inj, p_emp_real_vs_injdist=p_emp,
                     d2_real_vs_target=d2_zero, p_sig_real_vs_target=sig_equiv(d2_zero),
                     inj_dchi2_target_median=float(np.median(dchi)), inj_dchi2_target_p95=float(np.percentile(dchi, 95)),
                     frac_inj_dchi2_ge_real=float(np.mean(dchi >= r["real"]["dchi2_target"])),
                     inj_snr_median=float(tr["snr_target"].median()))

# shifted injections (off-target, non-circular)
shifts = {}
for p in sorted(glob.glob(os.path.join(OUT, "inj", "inj_*shift*.json"))):
    r = json.load(open(p))
    tr = pd.DataFrame(r["trials"])
    X = tr[["off_E", "off_N"]].values
    shifts[r["job"]] = dict(true=r["true_offset_arcsec"], mean=X.mean(0), std=X.std(0, ddof=1),
                            recovered_over_true=float(np.dot(X.mean(0), r["true_offset_arcsec"]) /
                                                      np.dot(r["true_offset_arcsec"], r["true_offset_arcsec"])),
                            dchi2_target_median=float(tr["dchi2_target"].median()),
                            dchi2_target_p5=float(np.percentile(tr["dchi2_target"], 5)),
                            dchi2_truepos_median=float(tr["dchi2_truepos"].median()))


def combine(names, bias_correct=False, robust=False):
    Ws, Wx = np.zeros((2, 2)), np.zeros(2)
    for n in names:
        s = sets[n]
        C = np.array(s["inj_cov_robust"] if robust else s["inj_cov"])
        x = np.array(s["real_off"]) - (np.array(s["inj_mean"]) if bias_correct else 0)
        Wi = np.linalg.inv(C)
        Ws += Wi; Wx += Wi @ x
    Cc = np.linalg.inv(Ws)
    mu = Cc @ Wx
    d2 = float(mu @ Ws @ mu)
    # consistency between the inputs
    chi2_between = 0.0
    for n in names:
        s = sets[n]
        C = np.array(s["inj_cov_robust"] if robust else s["inj_cov"])
        x = np.array(s["real_off"]) - (np.array(s["inj_mean"]) if bias_correct else 0)
        chi2_between += float((x - mu) @ np.linalg.solve(C, x - mu))
    p, sg = sig_equiv(d2)
    return dict(inputs=names, bias_corrected=bias_correct, robust_cov=robust, mean=mu, sig=np.sqrt(np.diag(Cc)),
                corr=Cc[0, 1] / np.sqrt(Cc[0, 0] * Cc[1, 1]), cov=Cc, chi2_vs_target=d2, p_vs_target=p,
                sigma_equiv_vs_target=sg, chi2_between_inputs=chi2_between, dof_between=2 * (len(names) - 1))


comb = {}
for sig in "AB":
    for label, names in {"S32tpfPC+S98clean29": [f"S32_tpf_poscorr_{sig}", f"S98_clean29_{sig}"],
                         "S32tspocPC+S98clean29": [f"S32_tspoc_poscorr_{sig}", f"S98_clean29_{sig}"],
                         "S6PC+S32tpfPC+S98clean29": [f"S6_tspoc_poscorr_{sig}", f"S32_tpf_poscorr_{sig}", f"S98_clean29_{sig}"],
                         "S32tpfPC+S98clean": [f"S32_tpf_poscorr_{sig}", f"S98_clean_{sig}"],
                         "S32tpfPC+S98all": [f"S32_tpf_poscorr_{sig}", f"S98_all_{sig}"],
                         "S32tpf+S98clean29": [f"S32_tpf_{sig}", f"S98_clean29_{sig}"],
                         "S32ffi+S98clean29": [f"S32_ffi_{sig}", f"S98_clean29_{sig}"],
                         "S6+S32tpf+S98clean29": [f"S6_ffi_{sig}" if sig == "B" else "S6_dumpx0.5h_A",
                                                  f"S32_tpf_{sig}", f"S98_clean29_{sig}"],
                         "S6+S98clean29": [f"S6_ffi_B", f"S98_clean29_B"] if sig == "B" else None,
                         "W3:S32tpf+S98clean29": ["S32_tpf_W3_A", "S98_clean29_W3_A"] if sig == "A" else None,
                         "W6:S32tpf+S98clean29": ["S32_tpf_W6_A", "S98_clean29_W6_A"] if sig == "A" else None,
                         "deg2:S32tpf+S98clean29": ["S32_tpf_deg2_A", "S98_clean29_deg2_A"] if sig == "A" else None,
                         "S32tpf+S98clean": [f"S32_tpf_{sig}", f"S98_clean_{sig}"],
                         "S32ffi+S98clean": [f"S32_ffi_{sig}", f"S98_clean_{sig}"],
                         "S32tpf+S98all": [f"S32_tpf_{sig}", f"S98_all_{sig}"],
                         "S6+S32tpf+S98clean": [f"S6_ffi_{sig}" if sig == "B" else "S6_dumpx0.5h_A",
                                                f"S32_tpf_{sig}", f"S98_clean_{sig}"],
                         "S6step+S32tpf+S98clean": ["S6_step_A", "S32_tpf_A", "S98_clean_A"] if sig == "A" else None,
                         "deg2:S32tpf+S98clean": ["S32_tpf_deg2_A", "S98_clean_deg2_A"] if sig == "A" else None,
                         "S6+S98clean": [f"S6_ffi_B", f"S98_clean_B"] if sig == "B" else None,
                         }.items():
        if names is None or any(n not in sets for n in names):
            continue
        for bc in [False, True]:
            for rb in [False, True]:
                comb[f"{sig}:{label}:bias{int(bc)}:rob{int(rb)}"] = combine(names, bc, rb)

# neighbours: sigma-distance of every capable star from the combined position (A)
cen = pd.read_csv(os.path.join(OUT, "neighbour_census_per_product.csv"))
manus_products = ["S6_tspoc_tp", "S32_tpf", "QLP98_3x3_best"]
sub = cen[cen["product"].isin(manus_products)]
agg = sub.groupby("source_id").agg(TIC=("TIC", "first"), Tmag=("Tmag", "first"), sep=("sep", "mean"),
                                   dE=("dE", "mean"), dN=("dN", "mean"),
                                   dreqA_min=("dreq_A", "min"), dreqA_max=("dreq_A", "max"),
                                   dreqA_lo_min=("dreq_A_lo", "min"), dreqA_lo_max=("dreq_A_lo", "max"),
                                   dreqB_min=("dreq_B", "min"), dreqB_max=("dreq_B", "max"))
allp = cen.groupby("source_id").agg(dreqA_anyproduct=("dreq_A", "min"))
agg = agg.join(allp)
nb_out = {}
for key in ["A:S32tpfPC+S98clean29:bias0:rob0", "A:S32tpfPC+S98clean29:bias1:rob0", "B:S32tpfPC+S98clean29:bias0:rob0",
            "A:S32tpf+S98clean29:bias0:rob0", "A:S32tpf+S98clean29:bias1:rob0", "A:S32tpf+S98clean29:bias0:rob1",
            "A:S32tpf+S98clean:bias0:rob0", "A:S6+S32tpf+S98clean29:bias0:rob0",
            "B:S32tpf+S98clean29:bias0:rob0", "B:S6+S32tpf+S98clean29:bias0:rob0"]:
    if key not in comb:
        continue
    c = comb[key]
    mu, C = np.array(c["mean"]), np.array(c["cov"])
    col = f"z_{key.split(':')[0]}_{key.split(':')[1]}_b{key.split(':')[2][-1]}"
    agg[col] = [float(np.sqrt((np.array([e, n]) - mu) @ np.linalg.solve(C, np.array([e, n]) - mu)))
                for e, n in zip(agg["dE"], agg["dN"])]
agg = agg.sort_values("sep")
agg.to_csv(os.path.join(OUT, "neighbours_vs_fitted_position.csv"))
for crit, mask in {"A_nominal500_min<=1": agg["dreqA_min"] <= 1, "A_nominal500_min<=0.5": agg["dreqA_min"] <= 0.5,
                   "A_low400_min<=1": agg["dreqA_lo_min"] <= 1, "A_nominal500_max<=1": agg["dreqA_max"] <= 1,
                   "A_nominal500_anyproduct<=1": agg["dreqA_anyproduct"] <= 1,
                   "B_nominal1000_min<=1": agg["dreqB_min"] <= 1}.items():
    c = agg[mask]
    zc = [k for k in agg.columns if k.startswith("z_")]
    nb_out[crit] = dict(n=int(len(c)), min_sep=float(c["sep"].min()) if len(c) else None,
                        min_z={k: float(c[k].min()) for k in zc} if len(c) else None,
                        stars=c.reset_index()[["source_id", "TIC", "Tmag", "sep", "dE", "dN", "dreqA_min", "dreqA_max"] + zc]
                        .to_dict("records"))

out = dict(sets=sets, shifts=shifts, combined=comb, neighbours=nb_out)
dump(out, "summary.json")

# console digest
for k, s in sets.items():
    print(f"{k:18s} n={s['n_trials']} ev={s['n_events']:2d} dep={s['depth_injected_ppm']:5.0f} real=({s['real_off'][0]:+6.1f},{s['real_off'][1]:+6.1f}) "
          f"inj mean=({s['inj_mean'][0]:+5.1f},{s['inj_mean'][1]:+5.1f}) sd=({s['inj_std'][0]:4.1f},{s['inj_std'][1]:4.1f}) "
          f"rob=({s['inj_std_robust'][0]:4.1f},{s['inj_std_robust'][1]:4.1f}) r={s['inj_corr']:+.2f} "
          f"p_emp={s['p_emp_real_vs_injdist']:.3f} dchi2 real={s['real_dchi2_target']:.1f} inj95={s['inj_dchi2_target_p95']:.1f} "
          f"f(inj>=real)={s['frac_inj_dchi2_ge_real']:.3f}")
for k, s in shifts.items():
    print(k, np.round(s["true"], 1), np.round(s["mean"], 1), np.round(s["std"], 1), round(s["recovered_over_true"], 3),
          "dchi2_t med", round(s["dchi2_target_median"], 1))
for k, c in comb.items():
    print(f"{k:42s} ({c['mean'][0]:+5.1f} +- {c['sig'][0]:3.1f}, {c['mean'][1]:+5.1f} +- {c['sig'][1]:3.1f}) r={c['corr']:+.2f} "
          f"chi2={c['chi2_vs_target']:.2f} p={c['p_vs_target']:.3f} ({c['sigma_equiv_vs_target']:.1f} sig) between={c['chi2_between_inputs']:.1f}/{c['dof_between']}")
for k, v in nb_out.items():
    print(k, v["n"], v["min_sep"], v["min_z"])
