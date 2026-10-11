"""Summaries from a1/a2 outputs: combined PRF localisation (inverse-variance, injection covariances), B depth
comparison without the corrupted own-S32 event, the frozen-search period explanation, and the O-C plot of signal A.
Outputs: out/a5_summary.json, out/a1_oc.png"""
import csv, glob, json, os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(H, "out")
S = {}
# ---------------- localisation
res = {}
for f in sorted(glob.glob(os.path.join(OUT, "a2_localisation_*.json"))):
    if f.endswith("_all98.json"): continue
    res.update(json.load(open(f)))
loc = {}
for sig in ("A", "B"):
    rows = {}
    for k, r in res.items():
        if not k.endswith("_" + sig) or "free_fit" not in r: continue
        mu = np.array(r["sky"]["real_minus_injmean_arcsec"]); C = np.array(r["sky"]["cov_arcsec2"])
        rows[k[:-2]] = dict(mu=mu, C=C, n_epochs=r["n_epochs"], z=[r["injections"]["z_x"], r["injections"]["z_y"]], maha=r["injections"]["mahalanobis"],
                            frac_further=r["injections"]["frac_injections_further"], depth=r["depth_ppm_prf"], dump_h=r["transit_to_nearest_dump_h"],
                            raw=[r["free_fit"]["dRA_arcsec"], r["free_fit"]["dDec_arcsec"]])
    def comb(keys):
        sel = [rows[k] for k in keys if k in rows]
        W = sum(np.linalg.inv(x["C"]) for x in sel); Cc = np.linalg.inv(W); m = Cc @ sum(np.linalg.inv(x["C"]) @ x["mu"] for x in sel)
        chi = float(m @ W @ m); chi_scatter = float(sum((x["mu"] - m) @ np.linalg.inv(x["C"]) @ (x["mu"] - m) for x in sel))
        return dict(sets=[k for k in keys if k in rows], offset_arcsec=[float(m[0]), float(m[1])], sigma_arcsec=[float(np.sqrt(Cc[0, 0])), float(np.sqrt(Cc[1, 1]))],
                    chi2_offset_2dof=chi, chi2_scatter=chi_scatter, dof_scatter=2 * (len(sel) - 1))
    loc[sig] = dict(per_set={k: dict(n_epochs=v["n_epochs"], offset_arcsec=[float(v["mu"][0]), float(v["mu"][1])], sigma_arcsec=[float(np.sqrt(v["C"][0, 0])), float(np.sqrt(v["C"][1, 1]))],
                                     z=v["z"], mahalanobis=v["maha"], frac_further=v["frac_further"], depth_ppm=v["depth"], transit_to_dump_h=v["dump_h"]) for k, v in rows.items()},
                    combined_S32tpf_S98clean=comb(["S32_tpf", "S98_ffi"]), combined_all_independent=comb(["S6_ffi", "S32_tpf", "S98_ffi"]),
                    combined_S32ffi_S98clean=comb(["S32_ffi", "S98_ffi"]))
    for lab in ("combined_S32tpf_S98clean", "combined_all_independent"):
        c = loc[sig][lab]; print(f"{sig} {lab}: ({c['offset_arcsec'][0]:+.1f}, {c['offset_arcsec'][1]:+.1f})\" ± ({c['sigma_arcsec'][0]:.1f}, {c['sigma_arcsec'][1]:.1f})\" chi2={c['chi2_offset_2dof']:.2f} scatter chi2={c['chi2_scatter']:.1f}/{c['dof_scatter']}")
S["localisation"] = loc
# ---------------- B depths without the corrupted own-S32 event (timing fit at scan edge, negative depth)
rows = list(csv.DictReader(open(os.path.join(OUT, "a1_per_transit.csv"))))
def depth_cmp(sel):
    d = np.array([float(r["D"]) for r in sel]); e = np.array([float(r["e_D_inj"]) for r in sel]); sec = np.array([int(r["sector"]) for r in sel])
    def wm(m): return np.sum(d[m] / e[m] ** 2) / np.sum(1 / e[m] ** 2), 1 / np.sqrt(np.sum(1 / e[m] ** 2))
    mu, emu = wm(np.ones(len(d), bool)); chi_c = float(np.sum(((d - mu) / e) ** 2)); chi_s = 0.0
    for s in np.unique(sec): ms, _ = wm(sec == s); chi_s += float(np.sum(((d[sec == s] - ms) / e[sec == s]) ** 2))
    d98, e98 = wm(sec == 98); da, ea = wm(sec != 98); n = len(d); k = len(np.unique(sec))
    return dict(n=n, per_sector={int(s): [float(x * 1e6) for x in wm(sec == s)] for s in np.unique(sec)}, common_ppm=[mu * 1e6, emu * 1e6], S98_ppm=[d98 * 1e6, e98 * 1e6],
                archival_ppm=[da * 1e6, ea * 1e6], diff_sigma=float((d98 - da) / np.hypot(e98, ea)), chi2_common=chi_c, chi2_per_sector=chi_s, dchi2=chi_c - chi_s,
                dBIC_common_minus_sep=float((chi_c + np.log(n)) - (chi_s + k * np.log(n))))
S["B_depths"] = {}
for tpl in ("B_mann", "B_free"):
    for ph in ("pipe", "own"):
        sel = [r for r in rows if r["signal"] == "B" and r["template"] == tpl and r["photometry"] == ph and r["edge"] != "True" and float(r["snr_white"]) > 0]
        S["B_depths"][f"{tpl}__{ph}"] = depth_cmp(sel); c = S["B_depths"][f"{tpl}__{ph}"]
        print(f"B {tpl} {ph}: n={c['n']} per-sector {{{', '.join(f'{k}: {v[0]:.0f}±{v[1]:.0f}' for k, v in c['per_sector'].items())}}} S98-arch = {c['diff_sigma']:.2f} sigma; dchi2(common-sep)={c['dchi2']:.2f} dBIC={c['dBIC_common_minus_sep']:.2f}")
# ---------------- signal A O-C plot and ephemeris explanation
A = json.load(open(os.path.join(OUT, "a1_results.json")))["signals"]["A"]
selA = [r for r in rows if r["signal"] == "A" and r["photometry"] == "pipe" and float(r["coverage"]) >= 0.7 and float(r["snr_white"]) >= 2 and r["edge"] != "True"]
eph = A["A__pipe"]["ephemeris"]["e_dt_inj_min"]
fig, ax = plt.subplots(figsize=(6, 3))
col = {6: "#1f4e9c", 32: "#b2182b", 98: "#1b7837"}
for r in selA:
    n = int(r["epoch"]); tc = float(r["tc_pred"]) + float(r["dt_min"]) / 1440; oc = (tc - (eph["T0"] + n * eph["P"])) * 1440
    ax.errorbar(n, oc, float(r["e_dt_inj_min"]), fmt="o", ms=3.5, color=col[int(r["sector"])], lw=0.8, capsize=0)
    ax.errorbar(n, oc, float(r["e_dt_white_min"]), fmt="none", ecolor="k", lw=1.6, capsize=0, alpha=0.6)
ax.axhline(0, color="k", lw=0.6); ax.set_xlabel("epoch (0 = Sector 98)"); ax.set_ylabel("O$-$C (min)")
ax.set_title(f"Signal A: P = {eph['P']:.7f} ± {eph['eP']:.7f} d, χ² = {eph['chi2']:.1f}/{eph['dof']} (injection errors; black: white-noise errors)", fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "a1_oc.png"), dpi=150)
# frozen-search period 3.182945 d: S98-only precision
s98 = A["A__pipe"]["S98_only"]
S["period_A"] = dict(final_linear_inj=dict(P=eph["P"], eP=eph["eP"], T0=eph["T0"], eT0=eph["eT0"], chi2=eph["chi2"], dof=eph["dof"], n=eph["n"], rms_oc_min=eph["rms_oc_min"]),
                     own_photometry=A["A__own"]["ephemeris"]["e_dt_inj_min"], white_errors=A["A__pipe"]["ephemeris"]["e_dt_white_min"], perm_errors=A["A__pipe"]["ephemeris"]["e_dt_perm_min"],
                     S98_only_linear=dict(P=s98["P"], eP=s98["eP"], n=s98["n"]), frozen_search_P=3.182945191903069,
                     frozen_search_minus_final_in_S98only_sigma=(3.182945191903069 - eph["P"]) / s98["eP"],
                     S98_only_minus_final_sigma=(s98["P"] - eph["P"]) / s98["eP"], cycle_offset_chi2=s98["chi2_vs_cycle_offset"],
                     alias_spacing_S32_d=eph["P"] / 569, alias_spacing_S6_d=eph["P"] / 792,
                     note="3.182945 d is the Sector 98 BLS grid period of the frozen search (single-sector, sigma ~4e-4 d); the S98-only linear ephemeris gives %.6f ± %.6f d" % (s98["P"], s98["eP"]))
print("period A:", {k: v for k, v in S["period_A"].items() if k in ("final_linear_inj", "S98_only_linear", "frozen_search_minus_final_in_S98only_sigma", "S98_only_minus_final_sigma", "alias_spacing_S32_d")})
json.dump(S, open(os.path.join(OUT, "a5_summary.json"), "w"), indent=1, default=float)
