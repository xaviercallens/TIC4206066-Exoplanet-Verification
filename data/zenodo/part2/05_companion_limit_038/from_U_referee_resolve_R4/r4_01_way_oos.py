"""R4 (a),(c): Was TIC 4206066 a training star of the Way et al. (2026) regressor, and what companion-mass limit
follows from an out-of-sample scatter?

(c) Way et al. 2026 (arXiv:2602.20096, Appendix A) define the Regression Sample, from which all training/testing sets
    are drawn, with 'AND gs.ruwe < 1.2'. The Full Sample (Appendix B, Table 2) drops the RUWE cut. A star with
    RUWE >= 1.2 was therefore never a training star: its Table-2 prediction is out-of-sample.
    Empirical check of in-sample shrinkage: compare the local O = pabsG - absG distribution of stars with
    RUWE in [1.0,1.2) (training-eligible) and [1.2,1.4) (never trained), same CMD window as the target.
(a) Companion-mass limit from the out-of-sample scatter, flux-sum conversion, masses from the Pecaut & Mamajek (2013)
    dwarf table (version 2022.04.16; M_G column is Gaia DR2 G), parsed here independently from the raw file.
Inputs (read-only): Way+2026 Zenodo overluminous_catalog.csv (md5 d27fb6e7... = Zenodo), 06_stellar Gaia 100-pc sample,
00_data fresh Gaia DR3 row. audit_env python.
"""
import json, os, csv
import numpy as np, pandas as pd
from astropy.table import Table

H = os.path.dirname(os.path.abspath(__file__))
AUD = os.path.abspath(os.path.join(H, "..", ".."))
WAY = "/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/06_stellar_dl/overluminous_catalog.csv"
SID = 3220388198192519424
res = {}

# fresh Gaia DR3 row (00_data)
rows = list(csv.DictReader(open(os.path.join(AUD, "00_data/catalogs/gaiadr3_gaia_source_cone240as.csv"))))
t = min(rows, key=lambda r: float(r["dist_arcsec_j2016pos"]))
assert int(t["source_id"]) == SID
ruwe_t = float(t["ruwe"])
res["target_ruwe_gaia_dr3"] = ruwe_t
res["regression_sample_cut"] = "Way+2026 Appendix A ADQL: ... AND gs.ruwe < 1.2 (Full Sample, Appendix B: no RUWE cut)"
res["target_in_regression_sample"] = bool(ruwe_t < 1.2)

ol = pd.read_csv(WAY, index_col=0)
row = ol[ol.source_id == SID].iloc[0]
O_t = float(row.pabsG_dr - row.absG_dr)
res["way_row"] = {"absG_dr": float(row.absG_dr), "pabsG_dr": float(row.pabsG_dr), "O": O_t, "ol": bool(row.ol)}

s = Table.read(os.path.join(AUD, "06_stellar/data/gaia_100pc_bprp_1.54_1.94.ecsv")).to_pandas()
m = ol.merge(s[["source_id", "ruwe", "ipd_frac_multi_peak", "non_single_star", "duplicated_source"]], on="source_id")
m["O"] = m.pabsG_dr - m.absG_dr
near = (np.abs(m.pabsG_dr - row.pabsG_dr) < 0.25) & (np.abs(m.pg_rp_dr - row.pg_rp_dr) < 0.05)


def rsig(x):
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def clip(x, k=2.5, it=10):
    mu, sd = np.median(x), rsig(x)
    for _ in range(it):
        sel = np.abs(x - mu) < k * sd
        mu, sd = x[sel].mean(), x[sel].std()
    return float(mu), float(sd)


def lower_tail(x, q=0.02275):
    # empirical 2-sigma-equivalent distance of the lower tail below the median (O below median = underluminous side)
    return float(np.median(x) - np.quantile(x, q))


stats = {}
for nm, sel in [("ruwe_1.0-1.2_training_eligible", near & (m.ruwe >= 1.0) & (m.ruwe < 1.2)),
                ("ruwe_1.2-1.4_never_trained", near & (m.ruwe >= 1.2) & (m.ruwe < 1.4)),
                ("ruwe_1.2-1.3_never_trained", near & (m.ruwe >= 1.2) & (m.ruwe < 1.3)),
                ("ruwe_lt1.2_all", near & (m.ruwe < 1.2)),
                ("ruwe_ge1.2_all", near & (m.ruwe >= 1.2))]:
    x = m.O[sel].to_numpy()
    mu, sd = clip(x)
    stats[nm] = {"n": int(sel.sum()), "median": float(np.median(x)), "sigma_MAD": rsig(x), "clipped_mean": mu,
                 "clipped_sigma": sd, "half_16_84": float((np.percentile(x, 84.135) - np.percentile(x, 15.865)) / 2),
                 "lower_tail_2sig_equiv": lower_tail(x), "frac_flagged_ol": float(m.ol[sel].mean()),
                 "frac_O_gt_target": float((x > O_t).mean())}
res["local_O_by_ruwe"] = stats

# Pecaut & Mamajek table, own parse
cols, tab = None, []
for line in open(os.path.join(AUD, "06_stellar/refs/EEM_dwarf_UBVIJHK_colors_Teff.txt")):
    if line.startswith("#SpT") and cols is None:
        cols = line[1:].split()
        continue
    if cols and line[:1] in "OBAFGKML" and line.split()[0].endswith("V"):
        p = line.split()
        tab.append(dict(zip(cols, p)))
    elif cols and tab and not line.strip():
        break
pm = pd.DataFrame(tab)
for c in ("M_G", "Msun", "M_Ks", "Mv", "V-Ic"):
    pm[c] = pd.to_numeric(pm[c].replace("...", np.nan), errors="coerce")
g = pm.dropna(subset=["M_G", "Msun"]).sort_values("M_G")
mass_of_MG = lambda x: float(np.interp(x, g.M_G, g.Msun))
MG_of_mass = lambda M: float(np.interp(-M, -g.Msun.to_numpy(), g.M_G.to_numpy()))
pmI = pm.dropna(subset=["Mv", "V-Ic", "Msun"]).copy(); pmI["M_I"] = pmI.Mv - pmI["V-Ic"]; pmI = pmI.sort_values("Msun")
MI_of_mass = lambda M: float(np.interp(M, pmI.Msun, pmI.M_I))

MG1 = float(row.pabsG_dr)


def limit(sig, k=2.0, O0=O_t):
    Ol = O0 + k * sig
    f = 10 ** (0.4 * Ol) - 1
    MG2 = MG1 - 2.5 * np.log10(f)
    return {"sigma": sig, "O_lim": Ol, "flux_ratio_G": f, "MG2": MG2, "mass_PM": mass_of_MG(MG2)}


lims = {}
for nm, sg in [("author_0.068", 0.068), ("clipped_0.052", 0.052), ("conservative_0.10", 0.10)]:
    lims[nm] = limit(sg)
for nm, st in stats.items():
    lims[f"{nm}_MAD"] = limit(st["sigma_MAD"])
    lims[f"{nm}_lowtail"] = limit(st["lower_tail_2sig_equiv"] / 2)
lims["paper_flag_threshold_0.376"] = {"flux_ratio_G": 10 ** (0.4 * 0.376) - 1,
                                      "mass_PM": mass_of_MG(MG1 - 2.5 * np.log10(10 ** (0.4 * 0.376) - 1))}
res["companion_mass_limits_2sigma"] = lims

# contrasts of companions of given mass relative to the primary (for the speckle calibration and MOLUSC regimes)
M1 = 0.606
MI1 = MI_of_mass(M1)
res["contrast_vs_mass"] = {f"{M:.2f}": {"dG": MG_of_mass(M) - MG1, "dI": MI_of_mass(M) - MI1} for M in (0.15, 0.2, 0.3, 0.35, 0.38, 0.42, 0.45, 0.5, 0.52)}
res["notes"] = ("dG uses the Way-predicted single-star M_G of the primary (7.758); dI uses the PM table at 0.606 Msun "
                "for the primary. PM M_G is Gaia DR2 G (table header).")
json.dump(res, open(os.path.join(H, "out", "r4_01_way_oos.json"), "w"), indent=1)
print(json.dumps(res, indent=1))
