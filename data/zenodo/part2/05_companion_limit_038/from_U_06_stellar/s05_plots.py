"""Audit 6: diagnostic plot of the observed SED (Gaia XP + 2MASS + AllWISE) and of the IR colours against the
empirical 100-pc locus. Reads only data/ and out/ of this audit."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.table import Table
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sed = json.load(open(HERE / "out/s02_sed.json"))
ph = sed["photometry"]
xp = pd.read_csv(HERE / "data/gaia_dr3_xp_sampled.csv")
c = 2.99792458e8
zp = {"J": (1235.0, 1594.0), "H": (1662.0, 1024.0), "Ks": (2159.0, 666.7), "W1": (3352.6, 309.540),
      "W2": (4602.8, 171.787), "W3": (11560.8, 31.674), "W4": (22088.3, 8.363)}  # W3/W4 for display only


def flam(m, l, f0):
    return f0 * 1e-26 * 10 ** (-0.4 * m) * c / (l * 1e-9) ** 2 * 1e-9


fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
ax[0].plot(xp.wavelength / 1000, xp.wavelength * xp.flux, c="#1b2f4f", lw=1.2, label="Gaia DR3 XP sampled")
for b, (l, f0) in zp.items():
    m, e = ph[b]
    y = l * flam(m, l, f0)
    ax[0].errorbar(l / 1000, y, yerr=y * 0.4 * np.log(10) * e, fmt="o", c="#c2410c", ms=5)
    ax[0].annotate(b, (l / 1000, y), textcoords="offset points", xytext=(4, 6), fontsize=8)
# expected photospheric W1-W4 from the target Ks and the 100-pc locus median Ks-Wn colour (same BP-RP)
irl = sed["ir_excess_empirical_locus"]
for b in ("W1", "W2", "W3", "W4"):
    l, f0 = zp[b]
    m_exp = ph["Ks"][0] - irl[f"Ks-{b}"]["locus_median"]
    ax[0].plot(l / 1000, l * flam(m_exp, l, f0), marker="s", mfc="none", mec="k", ms=9, ls="none",
               label="expected from Ks + 100-pc locus colour" if b == "W1" else None)
ax[0].set_xscale("log"); ax[0].set_yscale("log")
ax[0].set_xlabel("wavelength (micron)"); ax[0].set_ylabel(r"$\lambda F_\lambda$ (W m$^{-2}$)")
fb = sed["fbol_integration"]
ax[0].set_title(f"F_bol = {fb['fbol_MC']['median']*1e3:.3e} erg/s/cm2, L = {fb['L_MC_Lsun']['median']:.4f} Lsun", fontsize=9)
ax[0].legend(fontsize=8)

s = Table.read(HERE / "data/gaia_100pc_bprp_1.64_1.84_2mass_allwise.ecsv").to_pandas()
g = Table.read(HERE / "data/gaia_dr3_source.ecsv")[0]
ok = (np.abs(s.bp_rp - float(g["bp_rp"])) < 0.03) & (s.tm_qual.astype(str).str[2] == "A") & (s.w_sep < 1.0) & \
     (s.cc_flags.astype(str) == "0000") & (s.ext_flag == 0) & (s.ruwe < 1.4)
ir = sed["ir_excess_empirical_locus"]
for i, (col, (a, b), qi) in enumerate([("Ks-W3", ("ks_m", "w3mpro"), 2), ("Ks-W4", ("ks_m", "w4mpro"), 3)]):
    sel = ok & s.w_qual.astype(str).str[qi].isin(["A"])
    x = (s[a] - s[b])[sel]
    ax[1].hist(x, bins=np.linspace(-0.3, 0.8, 111), histtype="step", lw=1.5, color=["#1b2f4f", "#4a6fa5"][i],
               label=f"{col}: 100-pc locus (n={int(sel.sum())})", density=True)
    ax[1].axvline(ir[col]["target"], color=["#c2410c", "#d97706"][i], lw=2,
                  label=f"target {col} = {ir[col]['target']:.3f} +- {ir[col]['target_err']:.3f}")
ax[1].set_xlabel("colour (mag)"); ax[1].set_ylabel("density"); ax[1].legend(fontsize=7)
ax[1].set_title("Mid-IR colours vs stars with the same BP-RP (+-0.03)", fontsize=9)
fig.tight_layout(); fig.savefig(HERE / "out/s05_sed_ir.png", dpi=130)
print("ok")
