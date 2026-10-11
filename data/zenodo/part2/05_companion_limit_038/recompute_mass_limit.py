"""Standalone re-derivation of the ~0.38 Msun companion limit from the Gaia XP overluminosity (runs from this package).

Inputs (all inside the package, paths relative to this file):
  from_companion_check/data/way2026_overluminous_row.csv   target row of the Way et al. (2026) catalogue (absG_dr, pabsG_dr)
  ../01_inputs/reference_tables/mamajek_table.txt           Pecaut & Mamajek (2013) dwarf table, version 2022.04.16
Equations (same as from_companion_check/overluminosity.py):
  O      = pabsG_dr - absG_dr                  (predicted single-star M_G minus measured M_G)
  O_lim  = O + k * sigma                       (k = 2 or 3; sigma = scatter of O among comparable stars)
  f_G    = 10**(0.4 * O_lim) - 1               (largest companion/primary flux ratio in G)
  M_G,2  = pabsG_dr - 2.5 log10(f_G)           (faintest-allowed companion absolute magnitude)
  M_2    = interp(M_G,2) in the Pecaut & Mamajek M_G -> mass table
Only numpy and pandas are needed.
"""
import os
import numpy as np
import pandas as pd

H = os.path.dirname(os.path.abspath(__file__))
row = pd.read_csv(os.path.join(H, "from_companion_check/data/way2026_overluminous_row.csv")).iloc[0]
MG1 = float(row.pabsG_dr); O = float(row.pabsG_dr - row.absG_dr)

cols, tab = None, []
for line in open(os.path.join(H, "../01_inputs/reference_tables/mamajek_table.txt")):
    if line.startswith("#SpT") and cols is None:
        cols = line[1:].split()
    elif cols and line[:1] in "OBAFGKML" and line.split()[0].endswith("V") and len(line.split()) >= len(cols):
        tab.append(line.split()[:len(cols)])
    elif cols and tab and line.strip() == "":
        break
mj = pd.DataFrame(tab, columns=cols).replace("...", np.nan)[["SpT", "M_G", "Msun"]].dropna()
mj[["M_G", "Msun"]] = mj[["M_G", "Msun"]].astype(float); mj = mj.sort_values("M_G")
mass_of = lambda mg: float(np.interp(mg, mj.M_G, mj.Msun, left=np.nan, right=np.nan))

print(f"target: measured M_G {row.absG_dr}, predicted single-star M_G {MG1}, O = {O:+.3f} mag")
print(f"{'sigma':>7} {'k':>2} {'O_lim':>7} {'f_G':>7} {'M_G,2':>7} {'M_2 (Msun)':>10}")
for sigma, label in ((0.0682, "1.4826 MAD of 5205 comparable stars (from_companion_check/out/overluminosity.json)"),
                     (0.052, "clipped core (from_U_referee_resolve_R4/out/r4_01_way_oos.json)"),
                     (0.0534, "RUWE 1.0-1.2 training-eligible MAD (r4_01_way_oos.json)"),
                     (0.0541, "RUWE 1.2-1.3 never-trained MAD (r4_01_way_oos.json)"),
                     (0.10, "conservative (r4_01_way_oos.json)")):
    for k in (2, 3):
        Ol = O + k * sigma; f = 10 ** (0.4 * Ol) - 1; MG2 = MG1 - 2.5 * np.log10(f)
        print(f"{sigma:7.4f} {k:2d} {Ol:7.4f} {f:7.4f} {MG2:7.3f} {mass_of(MG2):10.3f}   {label if k == 2 else ''}")
f_thr = 10 ** (0.4 * 0.376) - 1
print(f"Way et al. flag threshold O = 0.376 mag -> f_G {f_thr:.3f} -> M_2 {mass_of(MG1 - 2.5*np.log10(f_thr)):.3f} Msun")
