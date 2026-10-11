"""Audit 9, steps 4,5,6,10: combine constraints on an unresolved bound companion (mass x separation), and the
eccentric-orbit fraction in the TRICERATOPS 1.1.0 planet eccentricity prior (Beta(0.867, 3.030), Kipping 2013;
uniform omega) that satisfies the photo-eccentric requirement rho_app <= rho_app,max for hosts below 0.38 Msun.
Astrometry (approximate, own model): photocentre semi-major axis alpha = plx * a * |q/(1+q) - beta_G|
(beta_G = G-band flux fraction of the companion); along-scan rms ~ 0.58 alpha (random orientation) for P < T=1038 d,
~ 0.58 alpha * 1.47 (T/P)^2 for P > T (residual curvature after a linear fit); excluded if rms > observed
astrometric_excess_noise (0.093 mas; all excess attributed to the companion = conservative)."""
import json
from pathlib import Path
import numpy as np
from scipy.stats import beta as betad
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import sys; sys.path.insert(0, str(Path(__file__).resolve().parent))
from c00_relations import companion_props

HERE = Path(__file__).resolve().parent
prof = json.load(open(HERE / "out/profile_A.json")); profB = json.load(open(HERE / "out/profile_B.json"))
rv = json.load(open(HERE / "out/rv_limits.json"))
PLX, AEN, M1, TBASE = 28.075, 0.0935, 0.606, 1038.0
MG1 = 7.758   # Way+2026 predicted single-star M_G (catalogue row; not re-downloaded here)
res = {}
# (5) overluminosity -> mass: O_lim -> companion M_G -> mass via Mamajek sequence (own relation)
Ms = np.linspace(0.08, 0.606, 2000); pr = companion_props(Ms)
for lab, O in [("2sig_sigma0.068", 0.003 + 2 * 0.0682), ("2sig_sigma0.07", 0.003 + 2 * 0.07), ("3sig", 0.003 + 3 * 0.0682)]:
    f = 10 ** (0.4 * O) - 1; MGc = MG1 - 2.5 * np.log10(f)
    res[f"overlum_{lab}"] = dict(O_lim=O, flux_ratio=f, MG_c_min=MGc, M_max=float(np.interp(MGc, pr["MG"][::-1], Ms[::-1])))
# (6) circular-orbit minimum host mass from the profile (interpolate dchi2 = 4 and 9)
rows = [r for r in prof["companions"] if r["M"] >= 0.1]
Mm = np.array([r["M"] for r in rows]); dc = np.array([r["dchi2_circ"] for r in rows]); o = np.argsort(Mm)
Mm, dc = Mm[o], dc[o]
def cross(th):
    i = np.where((dc[:-1] > th) & (dc[1:] <= th))[0][-1]
    return float(np.interp(th, [dc[i + 1], dc[i]], [Mm[i + 1], Mm[i]]))
res["circ_min_host_mass_dchi2_4"] = cross(4.0); res["circ_min_host_mass_dchi2_9"] = cross(9.0)
de = np.array([r["dchi2_ecc"] for r in rows])[o]
def cross_e(th):
    i = np.where((de[:-1] > th) & (de[1:] <= th))[0]
    return float(np.interp(th, [de[i[-1] + 1], de[i[-1]]], [Mm[i[-1] + 1], Mm[i[-1]]])) if len(i) else None
res["any_orbit_min_host_mass_dchi2_4"] = cross_e(4.0); res["any_orbit_min_host_mass_dchi2_9"] = cross_e(9.0)
res["rho_app_max_target_host_dchi2_4"] = [r for r in prof["companions"] if r["M"] == 0.606][0]["rho_app_max_2sig"]
# eccentric fraction under the TRICERATOPS prior
rng = np.random.default_rng(1); e = betad.rvs(0.867, 3.030, size=400000, random_state=rng); w = rng.uniform(0, 2 * np.pi, e.size)
g3 = ((1 + e * np.sin(w)) / np.sqrt(1 - e ** 2)) ** 3
tab = []
for r in rows:
    if r["M"] > 0.40: continue
    ok = r["rho"] * g3 <= r["rho_app_max_2sig"]
    # weight by transit probability factor (1+e sin w)/(1-e^2) as TRICERATOPS does
    wt = (1 + e * np.sin(w)) / (1 - e ** 2)
    tab.append(dict(M=r["M"], rho=r["rho"], rho_app_max=r["rho_app_max_2sig"], e_min=r["e_min_w270"],
                    frac_prior_ok=float(ok.mean()), frac_transit_weighted_ok=float((wt * ok).sum() / wt.sum()),
                    dchi2_ecc=r["dchi2_ecc"], f_c=r["f_c"]))
res["ecc_prior_fraction"] = tab
res["signalB_circ_dchi2_range_M0.14-0.606"] = [float(min(r["dchi2_circ"] for r in profB["companions"] if r["M"] >= 0.14)),
                                               float(max(r["dchi2_circ"] for r in profB["companions"] if r["M"] >= 0.14))]
# (4) astrometric limit grid and (3) RV limit grid
Pk = np.array([r["P_d"] for r in rv["limits_vs_period"]]); Kk = np.array([r["K95_sd_kms"] for r in rv["limits_vs_period"]])
agrid = np.logspace(-1.3, np.log10(78), 120); mgrid = np.linspace(0.08, 0.60, 105)
ast_ex = np.zeros((mgrid.size, agrid.size), bool); rv_ex = np.zeros_like(ast_ex); rv_prob = np.zeros(ast_ex.shape)
cosi = rng.uniform(0, 1, 3000); sini = np.sqrt(1 - cosi ** 2)
prm = companion_props(mgrid); FG1 = 10 ** (-0.4 * MG1)
for i, m in enumerate(mgrid):
    bG = 10 ** (-0.4 * prm["MG"][i]) / (FG1 + 10 ** (-0.4 * prm["MG"][i])); q = m / M1
    for j, a in enumerate(agrid):
        P = 365.25 * np.sqrt(a ** 3 / (M1 + m))
        alpha = PLX * a * abs(q / (1 + q) - bG)
        rms = 0.58 * alpha * (1.0 if P <= TBASE else min(1.0, 1.47 * (TBASE / P) ** 2))
        ast_ex[i, j] = rms > AEN
        K = (2 * np.pi * 6.674e-11 / (P * 86400)) ** (1 / 3) * m * 1.989e30 / ((M1 + m) * 1.989e30) ** (2 / 3) / 1e3
        Klim = float(np.exp(np.interp(np.log(P), np.log(Pk), np.log(Kk)))) if P <= Pk[-1] else np.inf
        rv_prob[i, j] = np.mean(K * sini > Klim); rv_ex[i, j] = rv_prob[i, j] > 0.95
ol = res["overlum_2sig_sigma0.068"]["M_max"]; shape_lo = res["any_orbit_min_host_mass_dchi2_4"] or 0.08
allowed = (~ast_ex) & (~rv_ex) & (mgrid[:, None] <= ol) & (mgrid[:, None] >= shape_lo)
res["allowed_summary"] = []
for m in [0.10, 0.15, 0.2, 0.25, 0.3, 0.35, 0.38]:
    i = np.argmin(abs(mgrid - m)); aa = agrid[allowed[i]]
    res["allowed_summary"].append(dict(M=float(mgrid[i]), a_min_au=float(aa.min()) if aa.size else None,
                                       a_max_au=float(aa.max()) if aa.size else None,
                                       rv_95_excluded_below_au=float(agrid[rv_ex[i]].max()) if rv_ex[i].any() else None,
                                       astro_excluded_au=[float(agrid[ast_ex[i]].min()), float(agrid[ast_ex[i]].max())] if ast_ex[i].any() else None))
json.dump(res, open(HERE / "out/companion_space.json", "w"), indent=1, default=float)
import csv
with open(HERE / "out/companion_space_table.csv", "w", newline="") as fh:
    wr = csv.writer(fh); wr.writerow(["M_c_Msun", "a_allowed_min_au", "a_allowed_max_au", "circular_orbit_allowed_A", "e_min_A_(w=270)", "frac_TRICERATOPS_ecc_prior_ok_A", "circular_allowed_B"])
    for s in res["allowed_summary"]:
        t = min(tab, key=lambda r: abs(r["M"] - s["M"])) if s["M"] <= 0.40 else None
        wr.writerow([round(s["M"], 3), s["a_min_au"] and round(s["a_min_au"], 2), s["a_max_au"] and round(s["a_max_au"], 1),
                     "no" if s["M"] < res["circ_min_host_mass_dchi2_4"] else "yes",
                     t and round(t["e_min"], 2), t and round(t["frac_prior_ok"], 3), "yes" if s["M"] >= 0.14 else "marginal"])
print(json.dumps({k: v for k, v in res.items() if k not in ("ecc_prior_fraction",)}, indent=1, default=float))
for r in tab: print(r)
# plot
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
ax[0].contourf(agrid, mgrid, rv_ex.astype(float), levels=[0.5, 1.5], colors=["#c9d7ea"])
ax[0].contourf(agrid, mgrid, ast_ex.astype(float), levels=[0.5, 1.5], colors=["#f2d1a8"], alpha=0.8)
ax[0].contourf(agrid, mgrid, allowed.astype(float), levels=[0.5, 1.5], colors=["#9bd19b"])
ax[0].axhline(ol, color="k", ls="--"); ax[0].axhline(shape_lo, color="k", ls=":")
ax[0].axhline(res["circ_min_host_mass_dchi2_4"], color="#b03030", ls="-.")
ax[0].set_xscale("log"); ax[0].set_xlabel("companion semi-major axis (au); 2.2\" = 78 au"); ax[0].set_ylabel("companion mass (Msun)")
ax[0].set_title("green: allowed; blue: Gaia RV (95%); orange: astrometric excess noise\n-- overluminosity 2σ, : transit shape (any e), -. circular-orbit minimum (A)", fontsize=8)
ax[1].plot([r["M"] for r in tab], [r["e_min"] for r in tab], "o-", label="e_min (ω=270°), signal A")
ax[1].plot([r["M"] for r in tab], [r["frac_prior_ok"] for r in tab], "s-", label="fraction of TRICERATOPS e-prior draws allowed")
ax[1].set_xlabel("host (companion) mass (Msun)"); ax[1].legend(fontsize=8); ax[1].set_ylim(0, 1)
fig.tight_layout(); fig.savefig(HERE / "plots/companion_space.png", dpi=130)
