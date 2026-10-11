"""Step 4: stellar parameters recomputed from primary catalogue values and the ORIGINAL relations.

Relations and where their coefficients were checked (files in data/refs/):
  * Mann et al. 2015 ApJ 804, 64 Table 1 + Erratum ApJ 819, 87 (2016) Table 1:
      R = 1.9515 - 0.3520 MK + 0.01680 MK^2, scatter 2.89 %  (valid 4.6 < MK < 9.8, K7-M7, 2700-4100 K)
      [Fe/H] form: (1.9305 - 0.3466 MK + 0.01647 MK^2)(1 + 0.04458 [Fe/H]), 2.70 %
      NOTE: the ORIGINAL published Table 1 printed f = "-0.04458" (press error, corrected in the erratum);
      arXiv v3 source and the erratum both give +0.04458.
      Teff(V-J) and BC relations from erratum Tables 2-3.
  * Mann et al. 2019 ApJ 871, 63 Table 6 (n = 5, zp = 7.5, sigma_e = 0.020); total relation error Table 7
      (2.4 % at MK = 5.0). Second route: official posterior code/data (github awmann/M_-M_K-, commit ffe5f19).
  * Stassun et al. 2019 AJ 158, 138 eq. (T from G, BP-RP), valid -0.2 < BP-RP < 3.5.
  * Lindegren et al. 2021 parallax zero point (gaiadr3-zeropoint 0.1.0).
"""
import csv
import json
import os
import sys

import numpy as np
from astropy.io import fits

sys.path.insert(0, os.path.dirname(__file__))
from audit_common import AUD, save_json  # noqa

rng = np.random.default_rng(20260925)
N = 1_600_000

g = list(csv.DictReader(open(os.path.join(AUD, "data/gaia/gaia_dr3_source.csv"))))[0]
plx, eplx = float(g["parallax"]), float(g["parallax_error"])
G, BP, RP = float(g["phot_g_mean_mag"]), float(g["phot_bp_mean_mag"]), float(g["phot_rp_mean_mag"])
bprp = float(g["bp_rp"])

# 2MASS (II/246) row
def vizrow(fn):
    L = [l.rstrip("\n") for l in open(os.path.join(AUD, "data/catalogs", fn)) if l.strip() and not l.startswith("#")]
    hdr = [h.strip() for h in L[0].split("|")]
    vals = [v.strip() for v in L[3].split("|")]
    return dict(zip(hdr, vals))

tm = vizrow("2mass_psc.tsv")
J, eJ, H, eH, K, eK = (float(tm[k]) for k in ["Jmag", "e_Jmag", "Hmag", "e_Hmag", "Kmag", "e_Kmag"])
ap = vizrow("apass_dr9.tsv")
V_apass, eV_apass = float(ap["Vmag"]), float(ap["e_Vmag"])
tic = vizrow("tic82_vizier.tsv")

from zero_point import zpt  # noqa
zpt.load_tables()
zp = float(zpt.get_zpt(np.array([G]), np.array([float(g["nu_eff_used_in_astrometry"])]), np.array([np.nan]),
                       np.array([float(g["ecl_lat"])]), np.array([int(g["astrometric_params_solved"])]))[0])
plx_c = plx - zp

out = {"inputs": dict(parallax_mas=plx, parallax_err=eplx, lindegren2021_zp_mas=zp, parallax_corrected=plx_c,
                      G=G, BP=BP, RP=RP, bp_rp=bprp, J=J, eJ=eJ, H=H, eH=eH, Ks=K, eKs=eK, V_APASS=V_apass,
                      eV_APASS=eV_apass, ruwe=float(g["ruwe"]))}

d_inv = 1000.0 / plx_c
out["distance_pc_1_over_plx_corrected"] = d_inv
out["distance_pc_1_over_plx_uncorrected"] = 1000.0 / plx

# Monte Carlo of observables
Ks = rng.normal(K, eK, N)
pl = rng.normal(plx_c, np.hypot(eplx, 0.010), N)  # add 10 uas systematic floor (Lindegren+2021 ZP residual ~ few uas)
MK = Ks + 5 * np.log10(pl / 100.0)
out["MK"] = dict(median=float(np.median(MK)), std=float(np.std(MK)),
                 point=float(K + 5 * np.log10(plx_c / 100.0)), point_uncorrected_plx=float(K + 5 * np.log10(plx / 100.0)))

# --- radius, Mann+2015 (MK only) ---
R0 = 1.9515 - 0.3520 * MK + 0.01680 * MK**2
R = R0 * (1 + rng.normal(0, 0.0289, N))
# [Fe/H] variants
def R_feh(feh):
    return (1.9305 - 0.3466 * MK + 0.01647 * MK**2) * (1 + 0.04458 * feh) * (1 + rng.normal(0, 0.0270, N))

R_m027 = R_feh(-0.27)
R_0 = R_feh(0.0)
R_pub_wrong_sign = (1.9305 - 0.3466 * MK + 0.01647 * MK**2) * (1 - 0.04458 * (-0.27))  # press-error sign, for reference

# --- mass, Mann+2019 n=5 ---
x = MK - 7.5
a = [-0.642, -0.208, -8.43e-4, 7.87e-3, 1.42e-4, -2.13e-4]
logM = sum(a[i] * x**i for i in range(6))
M_table = 10**logM
# relation error: Table 7 gives 2.4 % total at MK = 5.0 (coefficient + sigma_e); use 2.4 %
M = M_table * (1 + rng.normal(0, 0.024, N))
# route 2: official posterior of coefficients (Mk-M_7_trim.fits) exactly as in mk_mass.py (sigma_e from posterior)
post = fits.open(os.path.join(AUD, "data/refs/mann2019_repo/Mk-M_7_trim.fits"))[0].data
npost = post.shape[0]
idx = rng.integers(0, npost, N)
pa = post[idx, :6].T
sige = np.median(np.exp(post[:, 6]))
xk = MK - 7.5
M_post = 10 ** (pa[0] + pa[1] * xk + pa[2] * xk**2 + pa[3] * xk**3 + pa[4] * xk**4 + pa[5] * xk**5)
M_post = M_post * (1 + sige * rng.normal(size=N))
# Mann+2015 semi-empirical quartic (for comparison only)
MK0 = np.median(MK)
M_2015 = 0.5858 + 0.3872 * MK0 - 0.1217 * MK0**2 + 0.0106 * MK0**3 - 2.7262e-4 * MK0**4

rho = M / R**3
Msun_g, Rsun_cm, Gcgs = 1.98841e33, 6.957e10, 6.6743e-8
rho_sun = Msun_g / (4 / 3 * np.pi * Rsun_cm**3)
logg = np.log10(Gcgs * M * Msun_g / (R * Rsun_cm) ** 2)


def summ(arr):
    p = np.percentile(arr, [2.275, 15.865, 50, 84.135, 97.725])
    return dict(median=float(p[2]), minus1=float(p[2] - p[1]), plus1=float(p[3] - p[2]), mean=float(np.mean(arr)),
                std=float(np.std(arr)), p2_5=float(p[0]), p97_5=float(p[4]))

out["radius_mann2015_MK"] = summ(R)
out["radius_mann2015_MK_FeH_-0.27"] = summ(R_m027)
out["radius_mann2015_MK_FeH_0.0"] = summ(R_0)
out["radius_published_table_wrong_sign_FeH_-0.27_point"] = float(np.median(R_pub_wrong_sign))
out["mass_mann2019_table6_n5"] = summ(M)
out["mass_mann2019_table6_point"] = float(np.median(M_table))
out["mass_mann2019_official_posterior"] = summ(M_post)
out["mass_mann2019_sigma_e_posterior_median"] = float(sige)
out["mass_mann2015_quartic_point"] = float(M_2015)
out["density_rho_sun"] = summ(rho)
out["density_g_cm3"] = summ(rho * rho_sun)
out["rho_sun_g_cm3"] = rho_sun
out["logg_cgs"] = summ(logg)

# --- Teff routes ---
VJ = rng.normal(V_apass, eV_apass, N) - rng.normal(J, eJ, N)
def teff_vj(X, feh=None):
    if feh is None:
        return 3500 * (2.840 - 1.3453 * X + 0.3906 * X**2 - 0.0546 * X**3 + 0.002913 * X**4)
    return 3500 * (2.515 - 1.054 * X + 0.2965 * X**2 - 0.04150 * X**3 + 0.002245 * X**4 + 0.05262 * feh)

T_vj = teff_vj(VJ) + rng.normal(0, np.hypot(55, 60), N)  # relation scatter + recommended 60 K (erratum Tab.2 note)
T_vj_feh = teff_vj(VJ, -0.27) + rng.normal(0, np.hypot(42, 60), N)
out["V_minus_J"] = summ(VJ)
out["teff_mann2015_VJ"] = summ(T_vj)
out["teff_mann2015_VJ_FeH_-0.27"] = summ(T_vj_feh)

# --- bolometric luminosity routes (Mann+2015 erratum Table 3, m_bol,sun = -26.8167 -> M_bol,sun = 4.7554) ---
Mbol_sun = -26.8167 + 5 * np.log10(648000 / np.pi) + 5  # M = m + 5 log10(206265 AU/pc) + 5
BCK = 1.421 + 0.6084 * VJ - 0.09655 * VJ**2 - 6.263e-3 * VJ**3 + rng.normal(0, 0.036, N)
BCV = 0.5817 - 0.4168 * VJ - 0.08165 * VJ**2 - 4.084e-3 * VJ**3 + rng.normal(0, 0.016, N)
dm = 5 * np.log10(100.0 / pl)
L_K = 10 ** (-0.4 * (Ks + BCK - dm - Mbol_sun))
Vs = rng.normal(V_apass, eV_apass, N)
L_V = 10 ** (-0.4 * (Vs + BCV - dm - Mbol_sun))
out["Mbol_sun_used"] = Mbol_sun
out["L_from_BC_K"] = summ(L_K)
out["L_from_BC_V"] = summ(L_V)
# Stefan-Boltzmann R check from (L, Teff)
Tsun = 5772.0
R_SB = np.sqrt(L_K) / (T_vj / Tsun) ** 2
out["radius_stefan_boltzmann_LK_TVJ"] = summ(R_SB)
out["L_from_R_mann_and_Teff_4010"] = float(np.median(R) ** 2 * (4010 / Tsun) ** 4)

# --- TESS magnitude from Gaia (Stassun+2019 eq.) ---
c = bprp
Tmag = G - 0.00522555 * c**3 + 0.0891337 * c**2 - 0.633923 * c + 0.0324473
out["Tmag_stassun2019"] = Tmag
import pandas as pd  # noqa
_t = pd.read_csv(os.path.join(AUD, "data/catalogs/tic82_cone5arcmin.csv"))
_r = _t[_t["ID"] == 4206066].iloc[0]
out["TIC82"] = {k: (None if pd.isna(_r[k]) else (float(_r[k]) if not isinstance(_r[k], str) else _r[k]))
                for k in ["Tmag", "e_Tmag", "Teff", "e_Teff", "logg", "e_logg", "rad", "e_rad", "mass", "e_mass", "rho",
                          "e_rho", "lum", "e_lum", "d", "e_d", "contratio", "numcont", "Vmag", "Kmag", "TESSflag",
                          "SPFlag", "LUMCLASS", "MH", "plx", "ebv", "version"] if k in _t.columns}
out["Tmag_TIC82"] = out["TIC82"].get("Tmag")
out["gaia_dr2_flame"] = {"teff_val": 4069.3333, "radius_val": 0.61791587, "lum_val": 0.09432926,
                         "source": "data/gaia/gaia_dr2_source.csv"}
out["gaia_dr3_gspspec"] = {"teff": 4182.0, "logg": 5.01, "mh": -0.27, "flags_gspspec": "00000001000009999...",
                           "note": "GSP-Spec extrapolation flag set (8th char = 1); logg 5.01 inconsistent with M/R"}
save_json(out, "out/s04_stellar.json")
print(json.dumps({k: out[k] for k in out if k != "TIC82"}, indent=1, default=str)[:6000])
print(json.dumps(out["TIC82"]))
