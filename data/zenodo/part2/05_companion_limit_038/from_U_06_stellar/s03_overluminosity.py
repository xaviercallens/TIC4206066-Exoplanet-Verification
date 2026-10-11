"""Audit 6, step 7: Gaia overluminosity and the unresolved-companion mass limit.

(a) Way, Lepine, Gagne & Medan 2026, AJ 171, 252 (doi 10.3847/1538-3881/ae48ae; arXiv:2602.20096), Zenodo
    10.5281/zenodo.18500082 v1, overluminous_catalog.csv (md5 d27fb6e782fc2c46eb21b80707a90731 = Zenodo checksum,
    kept in the session scratchpad). Columns: absG_dr (measured, dereddened M_G), pabsG_dr (predicted single-star M_G
    from XP RP coefficients), g_rp_dr / pg_rp_dr, ol (flag, threshold 0.376 mag = 2.5 log10(2)/2).
    O = pabsG_dr - absG_dr (positive = brighter than predicted), as defined in the paper.
(b) Independent CMD check that does not use the Way et al. model: local main sequence M_G(BP-RP) from the Gaia DR3
    100-pc sample queried here (data/gaia_100pc_bprp_1.54_1.94.ecsv) and the Pecaut & Mamajek dwarf table.
(c) Companion mass limit: O_lim -> G flux ratio f = 10^(0.4 O_lim) - 1 -> companion M_G = M_G,single - 2.5 log10 f
    -> mass via (i) Pecaut & Mamajek M_G-Msun, (ii) PM M_G -> M_Ks then Mann+2019 Table 6.
Usage: python s03_overluminosity.py <dir with Way et al. csv>
"""
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.table import Table

HERE = Path(__file__).resolve().parent
DL = Path(sys.argv[1])
SID = 3220388198192519424
res = {}

# ---------------- Pecaut & Mamajek table ----------------
cols, rows = None, []
for line in open(HERE / "refs/EEM_dwarf_UBVIJHK_colors_Teff.txt"):
    if line.startswith("#SpT") and cols is None:
        cols = line[1:].split()
    elif cols and line.strip() and not line.startswith("#") and line.split()[0].endswith("V") and line[0] in "OBAFGKML":
        parts = line.split()
        if len(parts) >= len(cols):
            rows.append(parts[:len(cols)])
    elif cols and rows and line.strip() == "":
        break
pm = pd.DataFrame(rows, columns=cols).replace("...", np.nan)
for c in cols[1:-1]:
    pm[c] = pd.to_numeric(pm[c], errors="coerce")
pm.to_csv(HERE / "out/pecaut_mamajek_parsed.csv", index=False)
mg = pm.dropna(subset=["M_G", "Msun"]).sort_values("M_G")
mass_PM = lambda x: float(np.interp(x, mg.M_G, mg.Msun, left=np.nan, right=np.nan))
mk_tab = pm.dropna(subset=["M_G", "M_Ks"]).sort_values("M_G")
MKs_of_MG = lambda x: float(np.interp(x, mk_tab.M_G, mk_tab.M_Ks, left=np.nan, right=np.nan))
a6 = np.array([-0.642, -0.208, -8.43e-4, 7.87e-3, 1.42e-4, -2.13e-4])
mass_M19 = lambda mk: float(10 ** sum(a6[i] * (mk - 7.5) ** i for i in range(6)))

# ---------------- (a) Way et al. catalogue ----------------
ol = pd.read_csv(DL / "overluminous_catalog.csv", index_col=0)
row = ol[ol.source_id == SID]
assert len(row) == 1
row = row.iloc[0]
row.to_frame().T.to_csv(HERE / "out/way2026_target_row.csv", index=False)
O_all = (ol.pabsG_dr - ol.absG_dr).values
O_obs = float(row.pabsG_dr - row.absG_dr)
g = Table.read(HERE / "data/gaia_dr3_source.ecsv")[0]
MG_meas_mine = float(g["phot_g_mean_mag"]) + 5 * np.log10(float(g["parallax"]) / 100)
G_RP_mine = float(g["phot_g_mean_mag"] - g["phot_rp_mean_mag"])
res["way2026_row"] = {k: (row[k].item() if hasattr(row[k], "item") else row[k]) for k in row.index}
res["way2026_row"]["O_obs"] = O_obs
res["my_MG_from_gaia_dr3_no_extinction"] = MG_meas_mine
res["my_G_RP"] = G_RP_mine
res["catalogue_size"] = int(len(ol))
res["catalogue_frac_ol"] = float(ol.ol.mean())


def rsig(x):
    x = np.asarray(x)
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def core_sigma(x, n_iter=10, k=2.5):
    """iterative k-sigma clipped standard deviation (clips the binary tail)."""
    x = np.asarray(x)
    m, s = np.median(x), rsig(x)
    for _ in range(n_iter):
        sel = np.abs(x - m) < k * s
        m, s = np.mean(x[sel]), np.std(x[sel])
    return float(m), float(s)


variants = {}
for name, dM, dC in [("dMG0.25_dGRP0.05", 0.25, 0.05), ("dMG0.10_dGRP0.02", 0.10, 0.02), ("dMG0.50_dGRP0.10", 0.50, 0.10)]:
    sel = (np.abs(ol.pabsG_dr - row.pabsG_dr) < dM) & (np.abs(ol.pg_rp_dr - row.pg_rp_dr) < dC)
    x = O_all[sel.values]
    xs = O_all[(sel & ~ol.ol).values]
    m_c, s_c = core_sigma(x)
    p16, p50, p84 = np.percentile(x, [15.865, 50, 84.135])
    variants[name] = {"n": int(sel.sum()), "median": float(np.median(x)), "robust_sigma_MAD": rsig(x),
                      "half_16_84": float((p84 - p16) / 2), "clipped_mean": m_c, "clipped_sigma": s_c,
                      "not_flagged_robust_sigma": rsig(xs), "not_flagged_std": float(np.std(xs)),
                      "frac_O_gt_target": float(np.mean(x > O_obs)), "frac_flagged": float(ol.ol[sel].mean()),
                      "target_z_vs_median_MADsigma": float((O_obs - np.median(x)) / rsig(x))}
res["local_O_distribution"] = variants

# optional: clean-astrometry subset through the 100-pc Gaia sample (RUWE < 1.2, no NSS, no multi-peak IPD)
samp_fn = HERE / "data/gaia_100pc_bprp_1.54_1.94.ecsv"
if samp_fn.exists():
    s = Table.read(samp_fn).to_pandas()
    m = ol.merge(s, on="source_id", how="inner")
    m["O"] = m.pabsG_dr - m.absG_dr
    near = (np.abs(m.pabsG_dr - row.pabsG_dr) < 0.25) & (np.abs(m.pg_rp_dr - row.pg_rp_dr) < 0.05)
    clean = near & (m.ruwe < 1.2) & (m.non_single_star == 0) & (m.ipd_frac_multi_peak == 0) & (m.duplicated_source == False)
    rvok = clean & m.rv_chisq_pvalue.notna() & (m.rv_chisq_pvalue > 0.01)
    for nm, sel in [("near_all_in_sample", near), ("near_clean_astrometry", clean), ("near_clean_and_rv_constant", rvok)]:
        x = m.O[sel].values
        mc, sc = core_sigma(x) if len(x) > 20 else (np.nan, np.nan)
        variants[nm] = {"n": int(sel.sum()), "median": float(np.median(x)), "robust_sigma_MAD": rsig(x),
                        "clipped_mean": mc, "clipped_sigma": sc, "frac_O_gt_target": float(np.mean(x > O_obs))}

# ---------------- (c) companion mass limits for several sigma choices ----------------
MG1 = float(row.pabsG_dr)  # single-star M_G predicted from the spectrum
lim = {}
for sname, sig in [("sigma0.068_author", 0.068), ("sigma_MAD_local", variants["dMG0.25_dGRP0.05"]["robust_sigma_MAD"]),
                   ("sigma_clipped_local", variants["dMG0.25_dGRP0.05"]["clipped_sigma"]),
                   ("sigma_0.10_conservative", 0.10), ("sigma_0.05_optimistic", 0.05)]:
    for k in (1, 2, 3):
        O_lim = O_obs + k * sig
        f = 10 ** (0.4 * O_lim) - 1
        MG2 = MG1 - 2.5 * np.log10(f)
        mk2 = MKs_of_MG(MG2)
        lim[f"{sname}_{k}sigma"] = {"O_lim": O_lim, "flux_ratio_G": f, "MG2_min": MG2, "mass_PM": mass_PM(MG2),
                                    "MKs2_PM": mk2, "mass_Mann2019_via_PM_MKs": mass_M19(mk2)}
# paper threshold (what the ol flag itself excludes)
f_thr = 10 ** (0.4 * 0.376) - 1
MG2_thr = MG1 - 2.5 * np.log10(f_thr)
lim["paper_threshold_0.376"] = {"flux_ratio_G": f_thr, "MG2_min": MG2_thr, "mass_PM": mass_PM(MG2_thr),
                                "mass_Mann2019_via_PM_MKs": mass_M19(MKs_of_MG(MG2_thr))}
res["companion_limits"] = lim
res["primary_mass_PM_at_predicted_MG"] = mass_PM(MG1)

# flux ratio of companions of given mass (PM table) in G, Ks and TESS-like band, for the 'allowed companions' statement
comp = []
for _, r in pm.dropna(subset=["M_G", "Msun", "M_Ks"]).iterrows():
    if r.Msun > 0.62 or r.Msun < 0.075:
        continue
    fG = 10 ** (-0.4 * (r.M_G - MG1))
    fK = 10 ** (-0.4 * (r.M_Ks - 4.938))
    comp.append({"SpT": r.SpT, "Msun": r.Msun, "M_G": r.M_G, "M_Ks": r.M_Ks, "flux_ratio_G": fG,
                 "O_induced_mag": 2.5 * np.log10(1 + fG), "flux_ratio_Ks": fK})
res["companion_grid_PM"] = comp

# ---------------- (b) independent CMD check ----------------
cmd = {}
bprp_t = float(g["bp_rp"])
# Pecaut & Mamajek dwarf locus
pmc = pm.dropna(subset=["Bp-Rp", "M_G"]).sort_values("Bp-Rp")
MG_PM_at_colour = float(np.interp(bprp_t, pmc["Bp-Rp"], pmc["M_G"]))
cmd["PM_locus_MG_at_bprp"] = MG_PM_at_colour
cmd["overluminosity_vs_PM"] = MG_PM_at_colour - MG_meas_mine
if samp_fn.exists():
    s = Table.read(samp_fn).to_pandas()
    s["MG"] = s.phot_g_mean_mag + 5 * np.log10(s.parallax / 100)
    s["bpc"] = s.phot_bp_mean_mag - s.phot_rp_mean_mag
    # Riello+2021 (A&A 649, A3) corrected excess factor C* = C - f(BP-RP), Table 2 coefficients for 0.5 <= BP-RP < 4
    # (checked at https://arxiv.org/html/2012.01916v1). A loose fixed cut |C*| < CCUT removes blended photometry.
    c = s.bp_rp
    s["Cstar"] = s.phot_bp_rp_excess_factor - (1.162004 + 0.011464 * c + 0.049255 * c**2 - 0.005879 * c**3)
    CCUT = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
    cmd["Cstar_cut"] = CCUT
    base = (np.abs(s.Cstar) < CCUT) & (s.parallax_error / s.parallax < 0.05)
    clean = base & (s.ruwe < 1.2) & (s.non_single_star == 0) & (s.ipd_frac_multi_peak == 0) & (s.duplicated_source == False)
    for nm, sel in [("all_quality", base), ("clean_astrometry", clean)]:
        for w in (0.01, 0.02, 0.04):
            win = sel & (np.abs(s.bp_rp - bprp_t) < w)
            dM = s.MG[win].values
            mc, sc = core_sigma(dM) if len(dM) > 20 else (np.nan, np.nan)
            cmd[f"{nm}_win{w}"] = {"n": int(win.sum()), "median_MG": float(np.median(dM)), "clipped_mean_MG": mc,
                                   "robust_sigma_MAD": rsig(dM), "clipped_sigma": sc,
                                   "overluminosity_vs_median": float(np.median(dM) - MG_meas_mine),
                                   "overluminosity_vs_clipped_mean": mc - MG_meas_mine,
                                   "frac_brighter_than_target": float(np.mean(dM < MG_meas_mine))}
    # linear local fit M_G(BP-RP) within +-0.1 for clean stars, residual of the target
    win = clean & (np.abs(s.bp_rp - bprp_t) < 0.10)
    x, y = s.bp_rp[win].values, s.MG[win].values
    keep = np.ones_like(x, bool)
    for _ in range(10):
        p = np.polyfit(x[keep], y[keep], 2)
        r = y - np.polyval(p, x)
        sc = np.std(r[keep]); keep = np.abs(r - np.mean(r[keep])) < 2.5 * sc
    cmd["local_quadratic_fit_clean_pm0.10"] = {"n": int(win.sum()), "n_kept": int(keep.sum()), "coeffs": p.tolist(),
                                               "locus_MG_at_target": float(np.polyval(p, bprp_t)),
                                               "clipped_sigma": float(sc),
                                               "overluminosity": float(np.polyval(p, bprp_t) - MG_meas_mine)}
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    ax[0].scatter(s.bp_rp[base], s.MG[base], s=2, c="#4a6fa5", alpha=0.4, rasterized=True, label="100 pc, quality cuts")
    ax[0].scatter(s.bp_rp[clean], s.MG[clean], s=2, c="#1b2f4f", alpha=0.4, rasterized=True, label="+ RUWE<1.2, no NSS")
    xx = np.linspace(bprp_t - 0.1, bprp_t + 0.1, 50)
    ax[0].plot(xx, np.polyval(p, xx), c="#d97706", lw=2, label="clipped quadratic locus")
    ax[0].plot(pmc["Bp-Rp"], pmc["M_G"], "k--", lw=1, label="Pecaut & Mamajek")
    ax[0].scatter([bprp_t], [MG_meas_mine], marker="*", s=120, c="#c2410c", edgecolor="k", zorder=5, label="TIC 4206066")
    ax[0].set_xlim(1.55, 1.93); ax[0].set_ylim(9.2, 6.4); ax[0].set_xlabel("BP - RP"); ax[0].set_ylabel("M_G")
    ax[0].legend(fontsize=7, markerscale=1)
    sel = (np.abs(ol.pabsG_dr - row.pabsG_dr) < 0.25) & (np.abs(ol.pg_rp_dr - row.pg_rp_dr) < 0.05)
    ax[1].hist(O_all[sel.values], bins=np.linspace(-0.5, 1.1, 161), color="#4a6fa5")
    ax[1].axvline(O_obs, c="#c2410c", lw=2, label=f"target O = {O_obs:+.3f}")
    ax[1].axvline(0.376, c="k", ls="--", label="Way+2026 threshold 0.376")
    ax[1].set_xlabel("O = predicted - measured M_G (Way+2026)"); ax[1].set_ylabel(f"N (n={int(sel.sum())})")
    ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(HERE / "out/s03_cmd_overluminosity.png", dpi=130)
res["cmd"] = cmd
json.dump(res, open(HERE / "out/s03_overluminosity.json", "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
