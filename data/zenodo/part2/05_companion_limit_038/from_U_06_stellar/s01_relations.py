"""Audit 6, step 2-5: distance, M_Ks, Mann+2015 radius, Mann+2019 mass, stellar density (Monte Carlo).

Inputs are only the freshly queried catalogue rows in data/ (Gaia DR3 gaia_source, Bailer-Jones+2021, 2MASS PSC).
Relations are typed in from the primary sources (refs/):
  Mann et al. 2015, ApJ 804, 64, Table 1 (arXiv:1501.01635v3 source, refs/m15/main.tex l.1684) and the 2016 erratum
    (ApJ 819, 87; same coefficients, the erratum reprints Tables 1-3 because the journal version had press errors):
    R = 1.9515 - 0.3520 MK + 0.01680 MK^2, sigma = 2.89 % (fractional scatter in R)
    R = (1.9305 - 0.3466 MK + 0.01647 MK^2)(1 + 0.04458 [Fe/H]), sigma = 2.70 %
    validity (Sect. 5, l.251): 4.6 < MK < 9.8, 2700 < Teff < 4100 K, K7-M7, -0.6 < [Fe/H] < 0.5
  Mann et al. 2019, ApJ 871, 63, Table 6 (refs/m19/Masses_1.tex l.858), n = 5, zp = 7.5:
    log10 M = sum a_i (MK - 7.5)^i, a = [-0.642, -0.208, -8.43e-4, 7.87e-3, 1.42e-4, -2.13e-4], sigma_e = 0.020
    Table 7: total relation error 2.4 % at MK = 4.5 and 5.0 (coefficients + sigma_e), excludes input errors.
    Official posterior code github.com/awmann/M_-M_K- (commit ffe5f19), resources/Mk-M_7_trim.fits;
    validity (Sect. 8 l.1220): 4.0 < MK < 11.0, safer 4.5 < MK < 10.5.
  Lindegren et al. 2021, A&A 649, A4 parallax zero point (gaiadr3-zeropoint 0.1.0, zpt.get_zpt).
"""
import json
from pathlib import Path
import numpy as np
from astropy.table import Table
from astropy.io import fits

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"; OUT.mkdir(exist_ok=True)
rng = np.random.default_rng(6_2026_09_25)
N = 1_600_000  # same length as the official trimmed posterior

g = Table.read(HERE / "data/gaia_dr3_source.ecsv")[0]
bj = Table.read(HERE / "data/bailerjones2021.ecsv")[0]
tm = Table.read(HERE / "data/2mass_psc.ecsv")[0]
plx, eplx = float(g["parallax"]), float(g["parallax_error"])
Ks, eKs = float(tm["Kmag"]), float(tm["e_Kmag"])
res = {"inputs": {"gaia_dr3_source_id": int(g["source_id"]), "parallax_mas": plx, "parallax_error_mas": eplx,
                  "G": float(g["phot_g_mean_mag"]), "BP": float(g["phot_bp_mean_mag"]), "RP": float(g["phot_rp_mean_mag"]),
                  "bp_rp": float(g["bp_rp"]), "ruwe": float(g["ruwe"]), "non_single_star": int(g["non_single_star"]),
                  "Ks_2MASS": Ks, "e_Ks": eKs, "Ks_qual": str(tm["Qflg"]), "2MASS": str(tm["2MASS"])}}

# ---------------- (2) distance ----------------
from zero_point import zpt
zpt.load_tables()
# arrays of shape (1,) (the package's scalar path fails under numpy 2 / NEP 50); pseudocolour is not used for 5-p
# solutions (astrometric_params_solved = 31), so a placeholder inside the valid range is passed
zp = float(zpt.get_zpt(np.array([float(g["phot_g_mean_mag"])]), np.array([float(g["nu_eff_used_in_astrometry"])]),
                       np.array([1.5]), np.array([float(g["ecl_lat"])]),
                       np.array([int(g["astrometric_params_solved"])]))[0])
plx_c = plx - zp
res["distance"] = {
    "lindegren2021_zp_mas": zp, "parallax_corrected_mas": plx_c,
    "d_1_over_plx_pc": 1000 / plx, "d_1_over_plx_err_pc": 1000 / plx**2 * eplx,
    "d_1_over_plx_zpcorr_pc": 1000 / plx_c, "d_1_over_plx_zpcorr_err_pc": 1000 / plx_c**2 * eplx,
    "bailer_jones2021_geo_pc": [float(bj["r_med_geo"]), float(bj["r_lo_geo"]), float(bj["r_hi_geo"])],
    "bailer_jones2021_photogeo_pc": [float(bj["r_med_photogeo"]), float(bj["r_lo_photogeo"]), float(bj["r_hi_photogeo"])],
    "bj_flag": str(bj["flag"]),
    "note": "BJ21 uses the Lindegren+2021 zero point internally (Bailer-Jones et al. 2021 AJ 161, 147, Sect. 2.3)"}

# ---------------- (3) M_Ks ----------------
def mk_of(k, p_mas):
    return k + 5 * np.log10(p_mas / 100.0)

# parallax error: catalogue error; a second variant inflates it by 1.3 (El-Badry+2021 wide-binary calibration for
# bright stars) and adds 0.01 mas systematic (Lindegren+2021 ZP residual). Neither matters: distance term ~0.002 mag.
res["MK"] = {"MK_uncorrected_plx": mk_of(Ks, plx), "MK_zpcorr_plx": mk_of(Ks, plx_c),
             "MK_bj_geo": Ks - 5 * np.log10(float(bj["r_med_geo"]) / 10),
             "sigma_from_Ks": eKs, "sigma_from_plx": 5 / np.log(10) * eplx / plx,
             "sigma_total": float(np.hypot(eKs, 5 / np.log(10) * np.hypot(1.3 * eplx, 0.010) / plx))}
Ks_s = rng.normal(Ks, eKs, N)
pl_s = rng.normal(plx_c, np.hypot(1.3 * eplx, 0.010), N)
MK = mk_of(Ks_s, pl_s)

# extinction at 35.6 pc: see REPORT (Local Bubble); A_Ks = 0.114 A_V (Cardelli+1989 R_V=3.1) -> < 0.002 mag for A_V < 0.02


def summ(a):
    p = np.percentile(a, [2.275, 15.865, 50, 84.135, 97.725])
    return {"median": float(p[2]), "minus": float(p[2] - p[1]), "plus": float(p[3] - p[2]), "std": float(np.std(a)),
            "mean": float(np.mean(a)), "p2.3": float(p[0]), "p97.7": float(p[4])}


res["MK"]["MC"] = summ(MK)

# ---------------- (4a) Mann+2015 radius ----------------
def R15(mk):
    return 1.9515 - 0.3520 * mk + 0.01680 * mk**2


def R15feh(mk, feh):
    return (1.9305 - 0.3466 * mk + 0.01647 * mk**2) * (1 + 0.04458 * feh)


MK0 = res["MK"]["MK_zpcorr_plx"]
R_point = R15(MK0)
R_mc = R15(MK) * (1 + rng.normal(0, 0.0289, N))
R_mc_noscatter = R15(MK)
res["radius_mann2015"] = {
    "point_at_MK": R_point, "MC_scatter_plus_inputs": summ(R_mc), "MC_inputs_only": summ(R_mc_noscatter),
    "analytic_total_sigma": float(np.hypot(0.0289 * R_point, abs(-0.3520 + 2 * 0.01680 * MK0) * res["MK"]["sigma_total"])),
    "dR_dMK": -0.3520 + 2 * 0.01680 * MK0}
fehs = {"ARIADNE_-0.08+-0.18": (-0.08, 0.18), "GaiaGSPSpec_-0.27_flagged": (-0.27, 0.10), "solar_0.0": (0.0, 0.0)}
for k, (f, ef) in fehs.items():
    Rf = R15feh(MK, rng.normal(f, ef, N) if ef > 0 else f) * (1 + rng.normal(0, 0.0270, N))
    res["radius_mann2015"][f"with_FeH_{k}"] = summ(Rf)
# published (press-error) sign of f, for reference only
res["radius_mann2015"]["FeH_-0.27_wrong_sign_point"] = float((1.9305 - 0.3466 * MK0 + 0.01647 * MK0**2) * (1 - 0.04458 * -0.27))

# ---------------- (4b) Mann+2019 mass ----------------
a6 = np.array([-0.642, -0.208, -8.43e-4, 7.87e-3, 1.42e-4, -2.13e-4])


def M19(mk, a=a6):
    x = mk - 7.5
    return 10 ** sum(a[i] * x**i for i in range(6))


M_point = M19(MK0)
M_mc_tab7 = M19(MK) * (1 + rng.normal(0, 0.024, N))   # Table 7: 2.4 % total relation error at MK 4.5-5.0
M_mc_sige = M19(MK) * (1 + rng.normal(0, 0.020, N))   # sigma_e only (under-estimate: no coefficient error)
# official code route (re-implemented from mk_mass.py; posterior read from the cloned repository)
post = fits.open(HERE / "refs/mann2019_repo/resources/Mk-M_7_trim.fits")[0].data
npost = post.shape[0]
xk = MK[:npost] - 7.5
M_off = 10 ** (post[:, 0] + post[:, 1] * xk + post[:, 2] * xk**2 + post[:, 3] * xk**3 + post[:, 4] * xk**4 + post[:, 5] * xk**5)
sige_med = float(np.median(np.exp(post[:, 6])))
M_off = M_off * (1 + sige_med * rng.normal(size=npost))
# posterior-median relation at MK0 without input errors (how far the best-fit Table 6 sits from the posterior median)
x0 = MK0 - 7.5
M_post_curve = 10 ** (post[:, 0] + post[:, 1] * x0 + post[:, 2] * x0**2 + post[:, 3] * x0**3 + post[:, 4] * x0**4 + post[:, 5] * x0**5)
res["mass_mann2019"] = {
    "table6_point_at_MK": M_point, "table6_plus_table7_2.4pct_MC": summ(M_mc_tab7), "table6_plus_sigma_e_only_MC": summ(M_mc_sige),
    "official_posterior_MC": summ(M_off), "official_sigma_e_median": sige_med, "n_posterior": int(npost),
    "posterior_relation_at_MK0_no_input_errors": summ(M_post_curve),
    "table7_check": {str(m): float(M19(m)) for m in (4.5, 5.0, 5.5)},
    "mann2015_semi_empirical_quartic_point": float(0.5858 + 0.3872 * MK0 - 0.1217 * MK0**2 + 0.0106 * MK0**3 - 2.7262e-4 * MK0**4)}

# ---------------- (5) density ----------------
Rn = R_mc[:npost]
rho_tab = M_mc_tab7 / R_mc**3
rho_off = M_off / Rn**3
res["density_rho_sun"] = {"point": M_point / R_point**3, "MC_table6_table7_x_mann2015": summ(rho_tab),
                          "MC_official_posterior_x_mann2015": summ(rho_off)}
Msun, Rsun, Gc = 1.98840987e33, 6.957e10, 6.67430e-8  # IAU 2015 nominal (GM_sun / G), cgs
rho_sun = Msun / (4 / 3 * np.pi * Rsun**3)
res["density_rho_sun"]["rho_sun_cgs"] = rho_sun
res["density_rho_sun"]["MC_table6_cgs"] = summ(rho_tab * rho_sun)
res["logg_cgs_MC"] = summ(np.log10(Gc * M_mc_tab7 * Msun / (R_mc * Rsun) ** 2))

# ---------------- validity ----------------
res["validity"] = {"MK": MK0, "mann2015_MK_range": [4.6, 9.8], "mann2015_Teff_range": [2700, 4100],
                   "mann2015_FeH_range": [-0.6, 0.5], "mann2019_MK_range": [4.0, 11.0], "mann2019_safer_MK_range": [4.5, 10.5],
                   "MK_minus_mann2015_bright_limit": MK0 - 4.6,
                   "frac_MC_MK_below_4.6": float(np.mean(MK < 4.6))}

# ---------------- manuscript-rounding checks ----------------
res["rounded"] = {"d": round(1000 / plx_c, 1), "MK": round(MK0, 2), "R15": [round(res["radius_mann2015"]["MC_scatter_plus_inputs"]["median"], 3),
                  round(res["radius_mann2015"]["MC_scatter_plus_inputs"]["std"], 3)],
                  "M19_table": [round(res["mass_mann2019"]["table6_plus_table7_2.4pct_MC"]["median"], 3),
                                round(res["mass_mann2019"]["table6_plus_table7_2.4pct_MC"]["std"], 3)],
                  "M19_official": [round(res["mass_mann2019"]["official_posterior_MC"]["median"], 3),
                                   round(res["mass_mann2019"]["official_posterior_MC"]["std"], 3)],
                  "rho": [round(res["density_rho_sun"]["MC_table6_table7_x_mann2015"][k], 2) for k in ("median", "plus", "minus")],
                  "rho_official": [round(res["density_rho_sun"]["MC_official_posterior_x_mann2015"][k], 2) for k in ("median", "plus", "minus")]}
json.dump(res, open(OUT / "s01_relations.json", "w"), indent=1)
print(json.dumps(res, indent=1))
