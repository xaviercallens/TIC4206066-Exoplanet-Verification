"""Audit 6, step 6: independent SED check (F_bol, Teff, L, R) and infrared excess; check of the ARIADNE output files.

Independent of ARIADNE:
  Teff  - Mann+2015 Table 2 colour relations (V-J, r-J, with [Fe/H] or J-H), 60 K added in quadrature as Mann
          recommend (Table 2 note a); Casagrande+2021 (MNRAS 507, 2684) Gaia DR3 + 2MASS IRFM relations via the
          authors' official code colte.py (github.com/casaluca/colte, commit 6a2bc70); Pecaut & Mamajek table
          interpolation (BP-RP, V-Ks).
  F_bol - (a) Mann+2015 Table 3 bolometric corrections BC_V, BC_J, BC_K (V-J) and BC_r (r-J) with m_bol,sun =
          -26.8167 (their Sect. 6); (b) Pecaut & Mamajek BC_V interpolated in V-Ks, M_bol,sun = 4.74;
          (c) direct integration of the Gaia DR3 XP sampled spectrum (336-1020 nm, W m^-2 nm^-1) plus 2MASS/WISE
          photometry (Cohen+2003 and Wright+2010/Jarrett+2011 zero points) with a Rayleigh-Jeans tail beyond W2 and a
          scaled 4000 K blackbody (x0.5 +- 0.5) below 336 nm.
  R     - Stefan-Boltzmann with Teff_sun = 5772 K (IAU 2015 B3), Monte Carlo.
IR excess - Ks-W1..W4 and W1-W3/W4 of the target against an empirical photospheric locus: Gaia DR3 stars within
          100 pc with BP-RP within +-0.03 of the target (Gaia-2MASS/AllWISE best-neighbour cross-matches).
"""
import json, sys, pickle
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.table import Table

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
rng = np.random.default_rng(20260925)
N = 400_000
res = {}

g = Table.read(HERE / "data/gaia_dr3_source.ecsv")[0]
tm = Table.read(HERE / "data/2mass_psc.ecsv")[0]
aw = Table.read(HERE / "data/allwise.ecsv")[0]
cw = Table.read(HERE / "data/catwise2020.ecsv")[0]
ap = Table.read(HERE / "data/apass9.ecsv")[0]
ty = Table.read(HERE / "data/tycho2.ecsv")[0]
rel = json.load(open(OUT / "s01_relations.json"))
d_pc = rel["distance"]["d_1_over_plx_zpcorr_pc"]

ph = {"G": (float(g["phot_g_mean_mag"]), 0.0028), "BP": (float(g["phot_bp_mean_mag"]), 0.0029),
      "RP": (float(g["phot_rp_mean_mag"]), 0.0038),
      "V_APASS": (float(ap["Vmag"]), float(ap["e_Vmag"])), "B_APASS": (float(ap["Bmag"]), float(ap["e_Bmag"])),
      "g_APASS": (float(ap["g'mag"]), float(ap["e_g'mag"])), "r_APASS": (float(ap["r'mag"]), float(ap["e_r'mag"])),
      "i_APASS": (float(ap["i'mag"]), float(ap["e_i'mag"])),
      "BT": (float(ty["BTmag"]), float(ty["e_BTmag"])), "VT": (float(ty["VTmag"]), float(ty["e_VTmag"])),
      "J": (float(tm["Jmag"]), float(tm["e_Jmag"])), "H": (float(tm["Hmag"]), float(tm["e_Hmag"])),
      "Ks": (float(tm["Kmag"]), float(tm["e_Kmag"])),
      "W1": (float(aw["W1mag"]), float(aw["e_W1mag"])), "W2": (float(aw["W2mag"]), float(aw["e_W2mag"])),
      "W3": (float(aw["W3mag"]), float(aw["e_W3mag"])), "W4": (float(aw["W4mag"]), float(aw["e_W4mag"])),
      "W1_CatWISE": (float(cw["W1mproPM"]), float(cw["e_W1mproPM"])), "W2_CatWISE": (float(cw["W2mproPM"]), float(cw["e_W2mproPM"]))}
res["photometry"] = ph
res["allwise_flags"] = {"ccf": str(aw["ccf"]), "ex": int(aw["ex"]), "qph": str(aw["qph"]), "snr": [float(aw[f"snr{i}"]) for i in range(1, 5)],
                        "sat_frac": [float(aw[f"sat{i}"]) for i in range(1, 5)], "nb": int(aw["nb"])}
S = {k: rng.normal(v[0], v[1], N) for k, v in ph.items()}


def summ(a):
    a = np.asarray(a); a = a[np.isfinite(a)]
    p = np.percentile(a, [15.865, 50, 84.135])
    return {"median": float(p[1]), "minus": float(p[1] - p[0]), "plus": float(p[2] - p[1]), "std": float(np.std(a))}


# ---------------- Teff ----------------
T = {}
VJ, rJ, JH = S["V_APASS"] - S["J"], S["r_APASS"] - S["J"], S["J"] - S["H"]
feh_s = rng.normal(-0.08, 0.18, N)   # ARIADNE (loose); Mann relations are only mildly [Fe/H] dependent here
add60 = lambda n: rng.normal(0, 60, n)
poly = lambda X, c: sum(ci * X**i for i, ci in enumerate(c))
T["mann15_VJ"] = 3500 * poly(VJ, [2.840, -1.3453, 0.3906, -0.0546, 0.002913]) + rng.normal(0, 55, N) + add60(N)
T["mann15_VJ_feh"] = 3500 * (poly(VJ, [2.515, -1.054, 0.2965, -0.04150, 0.002245]) + 0.05262 * feh_s) + rng.normal(0, 42, N) + add60(N)
T["mann15_VJ_JH"] = 3500 * (poly(VJ, [2.769, -1.421, 0.4284, -0.06133, 0.003310]) + 0.1333 * JH + 0.05416 * JH**2) + rng.normal(0, 48, N) + add60(N)
T["mann15_rJ"] = 3500 * poly(rJ, [2.445, -1.2578, 0.4340, -0.0720, 0.004502]) + rng.normal(0, 58, N) + add60(N)
T["mann15_rJ_feh"] = 3500 * (poly(rJ, [2.532, -1.319, 0.4449, -0.07151, 0.004333]) + 0.05629 * feh_s) + rng.normal(0, 47, N) + add60(N)
T["mann15_rJ_JH"] = 3500 * (poly(rJ, [2.151, -1.092, 0.3767, -0.06292, 0.003950]) + 0.1697 * JH + 0.03106 * JH**2) + rng.normal(0, 52, N) + add60(N)
# Mann+2015 BP-RP relation used pre-launch synthetic Gaia passbands; shown for reference only
BPRP = S["BP"] - S["RP"]
T["mann15_BPRP_prelaunch_passbands_REFERENCE_ONLY"] = 3500 * poly(BPRP, [3.245, -2.4309, 1.043, -0.2127, 0.01649]) + rng.normal(0, 52, N) + add60(N)
# Pecaut & Mamajek interpolation
pm = pd.read_csv(OUT / "pecaut_mamajek_parsed.csv")
kk = pm[pm.SpT.isin(["K5V", "K6V", "K7V", "K8V", "K9V", "M0V", "M0.5V", "M1V"])]
T["PM_BPRP"] = np.interp(BPRP, kk["Bp-Rp"], kk["Teff"])
T["PM_VKs"] = np.interp(S["V_APASS"] - S["Ks"], kk["V-Ks"], kk["Teff"])
T["PM_GRP"] = np.interp(S["G"] - S["RP"], kk["G-Rp"], kk["Teff"])
res["teff"] = {k: summ(v) for k, v in T.items()}
res["colours"] = {"V-J": float(np.median(VJ)), "r-J": float(np.median(rJ)), "J-H": float(np.median(JH)), "BP-RP": float(np.median(BPRP)),
                  "V-Ks": float(ph["V_APASS"][0] - ph["Ks"][0]), "G-Ks": ph["G"][0] - ph["Ks"][0], "G-RP": ph["G"][0] - ph["RP"][0],
                  "BP-Ks": ph["BP"][0] - ph["Ks"][0]}

# Casagrande+2021 via official code (DR3, dwarfs); logg from Mann M, R; [Fe/H] -0.08 +- 0.2; E(B-V) = 0 +- 0.005
sys.path.insert(0, str(HERE / "refs/colte_repo"))
try:
    from colte import colte
    import os
    cwd = os.getcwd(); os.chdir(OUT)
    np.random.seed(20260925)
    colte(["TIC4206066"], [4.63], [-0.08], [ph["G"][0]], [ph["BP"][0]], [ph["RP"][0]], [ph["J"][0]], [ph["H"][0]], [ph["Ks"][0]],
          [0.0], DR3=True, MC=True, trials=20000, elogg=[0.05], efeh=[0.2], egg=[0.0028], ebp=[0.0029], erp=[0.0038],
          ej2=[ph["J"][1]], eh2=[ph["H"][1]], ek2=[ph["Ks"][1]], eebv=[0.005], outfile="s02_colte_casagrande2021.csv")
    os.chdir(cwd)
    txt = open(OUT / "s02_colte_casagrande2021.csv").read().splitlines()
    hdr = [h.strip().lstrip("#") for h in txt[0].split(",")]
    vals = [v.strip() for v in txt[1].split(",")]
    res["teff"]["casagrande2021_colte"] = dict(zip(hdr, vals))
except Exception as e:  # noqa
    res["teff"]["casagrande2021_colte"] = f"failed: {e!r}"

# ---------------- F_bol and L ----------------
L = {}
AU_pc = 1 / 206264.806
dist_s = 1000 / rng.normal(rel["distance"]["parallax_corrected_mas"], 0.025, N)
DM = 5 * np.log10(dist_s / 10)


def L_mann(mbol):  # Mann+2015: m_bol,sun = -26.8167 -> L/Lsun = 10^(-0.4 (m_bol + 26.8167)) (d/AU)^2
    return 10 ** (-0.4 * (mbol + 26.8167)) * (dist_s / AU_pc) ** 2


BC = {}
BC["V_VJ"] = poly(VJ, [0.5817, -0.4168, -0.08165, 4.084e-3]) + rng.normal(0, 0.016, N)
BC["J_VJ"] = poly(VJ, [0.8694, 0.3667, -0.02920]) + rng.normal(0, 0.016, N)
BC["K_VJ"] = poly(VJ, [1.421, 0.6084, -0.09655, 6.263e-3]) + rng.normal(0, 0.036, N)
BC["r_rJ"] = poly(rJ, [0.8958, -0.5081, -0.07387, 3.999e-3]) + rng.normal(0, 0.016, N)
BC["J_rJ"] = poly(rJ, [0.8790, 0.5068, -0.07791, 4.338e-3]) + rng.normal(0, 0.016, N)
BC["K_rJ"] = poly(rJ, [1.719, 0.5236, -0.09085, 6.735e-3]) + rng.normal(0, 0.036, N)
mags = {"V_VJ": S["V_APASS"], "J_VJ": S["J"], "K_VJ": S["Ks"], "r_rJ": S["r_APASS"], "J_rJ": S["J"], "K_rJ": S["Ks"]}
for k in BC:
    L[f"mann15_BC{k}"] = L_mann(mags[k] + BC[k])
res["BC_values"] = {k: summ(v) for k, v in BC.items()}
# Pecaut & Mamajek BC_V(V-Ks), IAU scale M_bol,sun = 4.74
bcv_pm = np.interp(S["V_APASS"] - S["Ks"], kk["V-Ks"], kk["BCv"])
L["PM_BCv_VKs"] = 10 ** (-0.4 * (S["V_APASS"] + bcv_pm - DM - 4.74))
res["BC_values"]["PM_BCv_VKs"] = summ(bcv_pm)

# (c) direct integration: Gaia XP sampled spectrum + photometry
xp = pd.read_csv(HERE / "data/gaia_dr3_xp_sampled.csv")
lam_nm = xp.wavelength.values; fl = xp.flux.values  # W m^-2 nm^-1
eflx = xp.flux_error.values
c_ms = 2.99792458e8
# F_lambda of photometric points (W m^-2 nm^-1) from F_nu zero points (Jy): F_lam = F_nu c / lam^2
zp = {"J": (1235.0, 1594.0), "H": (1662.0, 1024.0), "Ks": (2159.0, 666.7), "W1": (3352.6, 309.540), "W2": (4602.8, 171.787)}


def flam(mag, lam_nm_, fnu0):
    fnu = fnu0 * 1e-26 * 10 ** (-0.4 * mag)             # W m^-2 Hz^-1
    return fnu * c_ms / (lam_nm_ * 1e-9) ** 2 * 1e-9     # W m^-2 nm^-1


def fbol_once(m, xp_scale=1.0, uv_fac=0.5):
    # XP part (trapezoid)
    f_xp = np.trapezoid(fl * xp_scale, lam_nm)
    # UV part below 336 nm: 4000 K blackbody scaled to the 336-360 nm XP flux, times uv_fac
    from astropy.modeling.physical_models import BlackBody
    import astropy.units as u
    bb = BlackBody(temperature=4000 * u.K)
    lg = np.linspace(150, 336, 400)
    bl = bb(lg * u.nm).to_value(u.W / u.m**2 / u.Hz / u.sr) * c_ms / (lg * 1e-9) ** 2
    ref = (lam_nm >= 336) & (lam_nm <= 360)
    blr = bb(lam_nm[ref] * u.nm).to_value(u.W / u.m**2 / u.Hz / u.sr) * c_ms / (lam_nm[ref] * 1e-9) ** 2
    scale = np.mean(fl[ref] * xp_scale) / np.mean(blr)
    f_uv = uv_fac * np.trapezoid(bl * scale, lg)
    # IR part: XP red end (mean 1000-1020 nm) -> J -> H -> Ks -> W1 -> W2, log-log linear interpolation
    red = lam_nm >= 1000
    pts_l = [float(np.mean(lam_nm[red]))] + [zp[b][0] for b in ("J", "H", "Ks", "W1", "W2")]
    pts_f = [float(np.mean(fl[red] * xp_scale))] + [flam(m[b], zp[b][0], zp[b][1]) for b in ("J", "H", "Ks", "W1", "W2")]
    lg2 = np.logspace(np.log10(1020), np.log10(zp["W2"][0]), 2000)
    fi = 10 ** np.interp(np.log10(lg2), np.log10(pts_l), np.log10(pts_f))
    f_ir = np.trapezoid(fi, lg2)
    f_rj = pts_f[-1] * zp["W2"][0] / 3.0   # int_{l2}^inf F(l2)(l2/l)^4 dl = F(l2) l2 / 3
    return f_xp, f_uv, f_ir, f_rj


base = fbol_once({b: ph[b][0] for b in ("J", "H", "Ks", "W1", "W2")})
fb0 = sum(base)
# MC over photometric errors, an XP absolute-calibration uncertainty of 1 % (assumption, see REPORT), UV factor 0-1
M = 3000
fb = np.empty(M)
for i in range(M):
    m = {b: rng.normal(ph[b][0], ph[b][1]) for b in ("J", "H", "Ks", "W1", "W2")}
    fb[i] = sum(fbol_once(m, xp_scale=rng.normal(1, 0.01), uv_fac=rng.uniform(0, 1)))
# add 2 % interpolation systematic (molecular bands between the NIR points; assumption) in quadrature
fb = fb * (1 + rng.normal(0, 0.02, M))
Lsun_W = 3.828e26  # IAU 2015 B3 nominal
pc_m = 3.0856775814913673e16
d_mc = 1000 / rng.normal(rel["distance"]["parallax_corrected_mas"], 0.025, M)
L_int = 4 * np.pi * (d_mc * pc_m) ** 2 * fb / Lsun_W
res["fbol_integration"] = {"parts_W_m2": {"XP_336_1020nm": base[0], "UV_below_336nm_half_BB": base[1],
                                          "IR_1020nm_to_W2": base[2], "RJ_tail_beyond_W2": base[3]},
                           "fbol_W_m2": fb0, "fbol_erg_s_cm2": fb0 * 1e3, "frac_XP": base[0] / fb0,
                           "fbol_MC": summ(fb), "L_MC_Lsun": summ(L_int)}
L["xp_plus_photometry_integration"] = np.resize(L_int, N)
res["luminosity_Lsun"] = {k: summ(v) for k, v in L.items()}

# ---------------- radius from Stefan-Boltzmann ----------------
Tsun = 5772.0
Rres = {}
combos = [("mann15_BCK_VJ", "mann15_VJ"), ("mann15_BCV_VJ", "mann15_VJ"), ("mann15_BCK_VJ", "mann15_VJ_JH"),
          ("xp_plus_photometry_integration", "mann15_VJ"), ("xp_plus_photometry_integration", "mann15_rJ_JH"),
          ("PM_BCv_VKs", "PM_VKs"), ("mann15_BCK_rJ", "mann15_rJ")]
for lk, tk in combos:
    Rres[f"L:{lk}|T:{tk}"] = summ(np.sqrt(L[lk]) * (Tsun / T[tk]) ** 2)
# with ARIADNE Teff and the XP-integrated L
Rres["L:xp_integration|T:4012+-43 (ARIADNE)"] = summ(np.sqrt(L["xp_plus_photometry_integration"]) * (Tsun / rng.normal(4012, 43, N)) ** 2)
res["radius_stefan_boltzmann"] = Rres

# ---------------- ARIADNE output files ----------------
AR = HERE.parent.parent / "companion_check_tic4206066/sed_out"


def read_dat(fn):
    out = {}
    for line in open(fn):
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        out[p[0]] = p[1:]
    return out


res["ariadne_files"] = {"best_fit_average.dat": read_dat(AR / "best_fit_average.dat"),
                        "best_fit_sample.dat": read_dat(AR / "best_fit_sample.dat"),
                        "model_probabilities.dat": open(AR / "model_probabilities.dat").read().split("\n")}


class _Stub:
    def __init__(self, *a, **k): pass
    def __setstate__(self, st): self.__dict__["_state"] = st


class _U(pickle.Unpickler):
    def find_class(self, mod, name):
        try:
            return super().find_class(mod, name)
        except Exception:
            return type(name, (_Stub,), {})


try:
    bma = _U(open(AR / "BMA.pkl", "rb")).load()
    w = {k: float(v) for k, v in bma["weights"].items()}
    ari = {"weights": w, "sum_w2": float(sum(v**2 for v in w.values()))}
    for par in ("teff", "rad", "lum"):
        ari[par] = {"weighted_average_arrays_pct": summ(bma["weighted_average"][par]),
                    "weighted_samples_arrays_pct": summ(bma["weighted_samples"][par]),
                    "per_model": {gname: summ(o[par]) for gname, o in bma["originals"].items() if par in o},
                    "dat_average_best_unc": [float(bma["best_fit_averaged"][par])] + [float(x) for x in bma["uncertainties_averaged"][par]],
                    "dat_samples_best_unc": [float(bma["best_fit_samples"][par])] + [float(x) for x in bma["uncertainties_samples"][par]]}
        # a proper BMA mixture: draw from each model posterior in proportion to its weight
        pool = []
        for gname, o in bma["originals"].items():
            k = int(round(w[gname] * 200000))
            pool.append(rng.choice(np.asarray(o[par]), k))
        ari[par]["mixture_redrawn_pct"] = summ(np.concatenate(pool))
    res["ariadne_bma_pickle"] = ari
except Exception as e:  # noqa
    res["ariadne_bma_pickle"] = f"failed: {e!r}"

# ---------------- infrared excess vs an empirical locus ----------------
irf = HERE / "data/gaia_100pc_bprp_1.64_1.84_2mass_allwise.ecsv"
if irf.exists():
    s = Table.read(irf).to_pandas()
    bprp_t = float(g["bp_rp"])
    ok = (np.abs(s.bp_rp - bprp_t) < 0.03) & (s.tm_qual.astype(str).str[2] == "A") & (s.w_sep < 1.0) & \
         (s.cc_flags.astype(str) == "0000") & (s.ext_flag == 0) & (s.ruwe < 1.4)
    ir = {}
    tgt = {"Ks-W1": ph["Ks"][0] - ph["W1"][0], "Ks-W2": ph["Ks"][0] - ph["W2"][0], "Ks-W3": ph["Ks"][0] - ph["W3"][0],
           "Ks-W4": ph["Ks"][0] - ph["W4"][0], "W1-W2": ph["W1"][0] - ph["W2"][0], "W2-W3": ph["W2"][0] - ph["W3"][0],
           "J-Ks": ph["J"][0] - ph["Ks"][0], "H-Ks": ph["H"][0] - ph["Ks"][0]}
    etgt = {"Ks-W1": np.hypot(ph["Ks"][1], ph["W1"][1]), "Ks-W2": np.hypot(ph["Ks"][1], ph["W2"][1]),
            "Ks-W3": np.hypot(ph["Ks"][1], ph["W3"][1]), "Ks-W4": np.hypot(ph["Ks"][1], ph["W4"][1]),
            "W1-W2": np.hypot(ph["W1"][1], ph["W2"][1]), "W2-W3": np.hypot(ph["W2"][1], ph["W3"][1]),
            "J-Ks": np.hypot(ph["J"][1], ph["Ks"][1]), "H-Ks": np.hypot(ph["H"][1], ph["Ks"][1])}
    colmap = {"Ks-W1": ("ks_m", "w1mpro"), "Ks-W2": ("ks_m", "w2mpro"), "Ks-W3": ("ks_m", "w3mpro"), "Ks-W4": ("ks_m", "w4mpro"),
              "W1-W2": ("w1mpro", "w2mpro"), "W2-W3": ("w2mpro", "w3mpro"), "J-Ks": ("j_m", "ks_m"), "H-Ks": ("h_m", "ks_m")}
    qual_idx = {"w1mpro": 0, "w2mpro": 1, "w3mpro": 2, "w4mpro": 3}
    for cname, (a, b) in colmap.items():
        sel = ok.copy()
        for band in (a, b):
            if band in qual_idx:
                sel &= s.w_qual.astype(str).str[qual_idx[band]].isin(["A"])
        # like-for-like W1/W2 saturation regime: compare with stars of similar W1 brightness
        if "w1mpro" in (a, b) or "w2mpro" in (a, b):
            sel &= (s.w1mpro > 7.0) & (s.w1mpro < 8.3)
        x = (s[a] - s[b])[sel].values
        x = x[np.isfinite(x)]
        med = float(np.median(x)); rs = float(1.4826 * np.median(np.abs(x - med)))
        ex = tgt[cname] - med
        sig = float(np.hypot(rs, 0))  # locus scatter already includes typical photometric errors
        ir[cname] = {"target": tgt[cname], "target_err": float(etgt[cname]), "locus_median": med, "locus_robust_sigma": rs,
                     "n": int(len(x)), "excess_mag": ex, "excess_sigma_units": ex / np.hypot(rs, 0.0) if rs > 0 else np.nan,
                     "excess_z_target_err_only": ex / etgt[cname],
                     "max_flux_excess_2sigma_pct": float(100 * (10 ** (0.4 * max(ex + 2 * max(rs, etgt[cname]), 0)) - 1))}
    res["ir_excess_empirical_locus"] = ir
    res["ir_sample_selection"] = "100 pc, parallax_over_error > 20, |BP-RP - 1.739| < 0.03, 2MASS Ks qual A, AllWISE best neighbour < 1 arcsec, cc_flags 0000, ext_flag 0, RUWE < 1.4, band ph_qual A; W1/W2 colours restricted to 7.0 < W1 < 8.3 (same saturation regime as the target, W1 = 7.605)"
json.dump(res, open(OUT / "s02_sed.json", "w"), indent=1, default=float)
print(json.dumps({k: res[k] for k in res if k not in ("photometry",)}, indent=1, default=float)[:20000])
