"""Overluminosity check for TIC 4206066 (Gaia DR3 3220388198192519424) against Way et al. 2026 (AJ, arXiv:2602.20096).

Catalog: Zenodo 10.5281/zenodo.18500082, overluminous_catalog.csv (md5 d27fb6e782fc2c46eb21b80707a90731) and
catalog_of_systems_100pc_lite.csv (md5 3f1eb18bc2d027d6125647aa6ccef666). Overluminosity O = predicted - measured
M_G (from Gaia XP RP coefficients, random-forest model trained on single K/M dwarfs); paper threshold 0.376 mag
(= 2.5 log10(2) / 2, half the equal-mass offset).

Limit on a companion: O_lim = O_obs + 2 sigma_local, sigma_local = robust scatter (1.4826 MAD) of O for stars with a
similar predicted M_G and colour. Flux ratio in G: f = 10**(0.4 O_lim) - 1 (same flux-sum logic the paper uses for
the equal-mass offset 2.5 log10 2). Companion M_G = M_G,single - 2.5 log10 f -> mass by interpolating the
Pecaut & Mamajek dwarf table (v2022.04.16, copy in ../deep_followup_tic4206066/archival_checks/mamajek_table.txt).
"""
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DL = Path(sys.argv[1])  # folder with the two downloaded Zenodo CSVs
SID = 3220388198192519424

ol = pd.read_csv(DL / "overluminous_catalog.csv", index_col=0)
sysc = pd.read_csv(DL / "catalog_of_systems_100pc_lite.csv", low_memory=False)
row = ol[ol.source_id == SID].iloc[0]
srow = sysc[sysc.source_id_1 == SID]
row.to_frame().T.to_csv(HERE / "data/way2026_overluminous_row.csv", index=False)
srow.to_csv(HERE / "data/way2026_systems_row.csv", index=False)

O_all = ol.pabsG_dr - ol.absG_dr
O_obs = float(row.pabsG_dr - row.absG_dr)

def robust_sigma(x):
    x = np.asarray(x); return 1.4826 * np.median(np.abs(x - np.median(x)))

# neighbours in predicted M_G and colour (similar spectrum)
sel = (np.abs(ol.pabsG_dr - row.pabsG_dr) < 0.25) & (np.abs(ol.pg_rp_dr - row.pg_rp_dr) < 0.05)
O_loc = O_all[sel]
sig_loc = robust_sigma(O_loc)
med_loc = float(np.median(O_loc))
frac_brighter = float(np.mean(O_loc > O_obs))

# Mamajek dwarf table
cols = None; tab = []
for line in open(HERE.parent / "deep_followup_tic4206066/archival_checks/mamajek_table.txt"):
    if line.startswith("#SpT") and cols is None:
        cols = line[1:].split()
    elif cols and line[:1] in "OBAFGKML" and line.split()[0].endswith("V") and len(line.split()) >= len(cols):
        tab.append(line.split()[:len(cols)])
    elif cols and tab and line.strip() == "":
        break
mj = pd.DataFrame(tab, columns=cols).replace("...", np.nan)
mj = mj[["SpT", "M_G", "Msun", "Bp-Rp", "M_Ks"]].dropna(subset=["M_G", "Msun"])
mj[["M_G", "Msun"]] = mj[["M_G", "Msun"]].astype(float)
mj = mj.sort_values("M_G")
mass_of = lambda mg: float(np.interp(mg, mj.M_G, mj.Msun, left=np.nan, right=np.nan))

MG1 = float(row.pabsG_dr)  # single-star M_G predicted from the spectrum
res = {"source_id": SID, "measured_MG": float(row.absG_dr), "predicted_MG": MG1,
       "measured_G_RP": float(row.g_rp_dr), "predicted_G_RP": float(row.pg_rp_dr),
       "O_obs_mag": O_obs, "flag_overluminous": bool(row.ol), "in_banyan_sigma": bool(row.in_banyan_sigma),
       "systems_catalog": srow.iloc[0][["catalogs_1", "ol_1", "n_components"]].astype(str).to_dict(),
       "local_bin": {"n": int(sel.sum()), "median_O": med_loc, "robust_sigma": sig_loc,
                     "frac_with_larger_O": frac_brighter, "frac_flagged_ol": float(ol.ol[sel].mean())},
       "paper_threshold_mag": 0.376, "limits": {}}
for k in (2, 3):
    O_lim = O_obs + k * sig_loc
    f = 10 ** (0.4 * O_lim) - 1
    MG2 = MG1 - 2.5 * np.log10(f)
    res["limits"][f"{k}sigma"] = {"O_lim_mag": O_lim, "flux_ratio_G_max": f, "companion_MG_min": MG2,
                                   "companion_mass_max_Msun": mass_of(MG2),
                                   "mass_ratio_max": mass_of(MG2) / 0.606}
# paper threshold itself
f_thr = 10 ** (0.4 * 0.376) - 1
res["limits"]["paper_threshold"] = {"flux_ratio_G_max": f_thr, "companion_MG_min": MG1 - 2.5 * np.log10(f_thr),
                                    "companion_mass_max_Msun": mass_of(MG1 - 2.5 * np.log10(f_thr))}
res["mamajek_primary_mass_at_MG1"] = mass_of(MG1)
json.dump(res, open(HERE / "out/overluminosity.json", "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
bins = np.linspace(-0.6, 1.2, 181)
ax[0].hist(O_loc, bins=bins, color="#3b6ea5")
ax[0].axvline(O_obs, color="#c0392b", lw=2, label=f"TIC 4206066: {O_obs:+.3f}")
ax[0].axvline(0.376, color="k", ls="--", label="paper threshold 0.376")
ax[0].axvline(0.753, color="k", ls=":", label="equal-mass 0.753")
ax[0].set_xlabel("overluminosity  O = predicted − measured M_G (mag)"); ax[0].set_ylabel(f"stars (n={sel.sum()})")
ax[0].legend(fontsize=8); ax[0].set_title("Stars with the same XP spectrum class")
ax[1].scatter(ol.g_rp_dr, ol.absG_dr, s=0.2, c=np.where(ol.ol, "#e67e22", "#3b6ea5"), rasterized=True)
ax[1].scatter([row.g_rp_dr], [row.absG_dr], s=60, c="#c0392b", edgecolor="k", zorder=5)
ax[1].set_xlim(0.7, 1.1); ax[1].set_ylim(9.5, 6.0)
ax[1].set_xlabel("G − RP"); ax[1].set_ylabel("M_G"); ax[1].set_title("Way+2026 sample (orange = overluminous)")
fig.tight_layout(); fig.savefig(HERE / "plots/overluminosity.png", dpi=130)
