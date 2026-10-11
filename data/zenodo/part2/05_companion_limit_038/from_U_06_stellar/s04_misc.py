"""Audit 6, small checks: TESS magnitude (Stassun+2019 AJ 158, 138 eq. 1, checked at arxiv.org/html/1905.10694v1),
host-density threshold behind the manuscript's '0.54 Msun' statement (Mann+2015 R(MK) x Mann+2019 M(MK), and the
Pecaut & Mamajek table), spectral-type bracket from the PM table, and the distance-modulus / extinction bookkeeping."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.table import Table

HERE = Path(__file__).resolve().parent
res = {}
g = Table.read(HERE / "data/gaia_dr3_source.ecsv")[0]
tic = Table.read(HERE / "data/tic82.ecsv")[0]


def tmag(G, c):
    return G - 0.00522555 * c**3 + 0.0891337 * c**2 - 0.633923 * c + 0.0324473


res["Tmag"] = {"TIC82_catalogue": float(tic["Tmag"]), "TIC82_e_Tmag": float(tic["e_Tmag"]),
               "stassun2019_from_DR2_phot_in_TIC": float(tmag(float(tic["GAIAmag"]), float(tic["gaiabp"]) - float(tic["gaiarp"]))),
               "stassun2019_from_DR3_phot": float(tmag(float(g["phot_g_mean_mag"]), float(g["bp_rp"])))}

# host density threshold rho = 3.5 rho_sun (manuscript: density-free fit, 97.7 % upper limit)
a6 = np.array([-0.642, -0.208, -8.43e-4, 7.87e-3, 1.42e-4, -2.13e-4])
M19 = lambda mk: 10 ** sum(a6[i] * (mk - 7.5) ** i for i in range(6))
R15 = lambda mk: 1.9515 - 0.3520 * mk + 0.01680 * mk**2
mk = np.linspace(4.6, 9.0, 4401)
rho = M19(mk) / R15(mk) ** 3
i = np.argmin(np.abs(rho - 3.5))
res["host_density_threshold_mann"] = {"rho_limit": 3.5, "MK_at_limit": float(mk[i]), "M_at_limit": float(M19(mk[i])),
                                      "R_at_limit": float(R15(mk[i]))}
pm = pd.read_csv(HERE / "out/pecaut_mamajek_parsed.csv")
sub = pm[pm.SpT.isin(["K7V", "K8V", "K9V", "M0V", "M0.5V", "M1V", "M1.5V", "M2V"])][["SpT", "Teff", "Msun", "R_Rsun", "Bp-Rp", "M_G", "M_Ks"]].copy()
sub["rho_rhosun"] = sub.Msun / sub.R_Rsun**3
res["PM_density_by_SpT"] = sub.to_dict(orient="records")
# spectral type bracket
bprp = float(g["bp_rp"])
MG = float(g["phot_g_mean_mag"]) + 5 * np.log10(float(g["parallax"]) / 100)
res["spt_bracket"] = {"bp_rp": bprp, "M_G": MG,
                      "PM_rows": pm[pm.SpT.isin(["K6V", "K7V", "K8V", "K9V", "M0V"])][["SpT", "Teff", "Bp-Rp", "M_G", "M_Ks", "Msun"]].to_dict(orient="records")}
json.dump(res, open(HERE / "out/s04_misc.json", "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
