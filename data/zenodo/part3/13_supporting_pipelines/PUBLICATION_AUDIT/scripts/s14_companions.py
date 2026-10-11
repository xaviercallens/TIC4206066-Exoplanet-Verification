"""Step 14: resolved-neighbour census and hidden-companion constraints.
(1) Gaia DR3 sources within 180'': TESS magnitude (Stassun+2019; missing colour -> G-0.430 with 0.6 mag and the
    BRIGHTER 1-sigma value used), flux ratio, eclipse depth required to produce each signal with (a) full inclusion,
    (b) PRF-integrated fluxes in the actual apertures (S6 TESS-SPOC optimal 18 px, S32 SPOC optimal 18 px, S98 QLP 3x3,
    and an r = 2 px circle), official TESS PRF (TESS_PRF) at each sector's position.
(2) RUWE of the target vs. the RUWE distribution of Gaia DR3 stars with similar G, BP-RP and distance (TAP query).
(3) Photometric density constraint on the host of signal A -> which dwarf spectral types could host it on a circular
    orbit (Pecaut & Mamajek table v2022.04.16), vs. the Gaia-XP overluminosity limit (Way et al. 2026 row).
(4) WISE excess check (K-W1..W4 vs dwarf sequence).
(5) Gaia DR3 RV scatter -> approximate K limit and companion mass limit vs period.
"""
import csv
import glob
import json
import os
import sys
import urllib.parse

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

sys.path.insert(0, os.path.dirname(__file__))
import lcio  # noqa
import prf_tools as pt  # noqa
from audit_common import AUD, OLD_COMP, fetch, save_json  # noqa

TARGET = "3220388198192519424"
DEPTH = {"A": 480e-6, "B": 1000e-6}
out = {}

# ---------------- (1) neighbour census with aperture-integrated PRF fluxes
aps = json.load(open(os.path.join(AUD, "out", "fpp", "inputs", "apertures.json")))
cubes = {6: "data/tesscut/tess-s0006-*_astrocut.fits", 32: "data/tesscut/tess-s0032-*_astrocut.fits", 98: "data/tesscut/tess-s0098-*_astrocut.fits"}
cam_ccd = {6: (1, 1), 32: (1, 2), 98: (1, 2)}
src = pt.gaia_sources(180)
rows = []
per_sector = {}
for sec, pat in cubes.items():
    h = fits.open(glob.glob(os.path.join(AUD, pat))[0])
    w = WCS(h[2].header)
    c0, r0 = h[1].header["1CRV4P"], h[1].header["2CRV4P"]
    tmid = float(np.nanmedian(h[1].data["TIME"]))
    sc = pt.scene(src, w, tmid, int(TARGET), Tlim=22)
    tg = [s for s in sc if s["is_target"]][0]
    prf = pt.get_prf(cam_ccd[sec][0], cam_ccd[sec][1], sec, c0 + tg["x"], r0 + tg["y"])
    shape = (21, 21)
    ap_abs = np.array(aps[str(sec)]["pixels"])            # absolute CCD (col, row)
    apmask = np.zeros(shape, bool)
    for col, row in ap_abs:
        x, y = col - c0, row - r0
        if 0 <= x < 21 and 0 <= y < 21:
            apmask[int(y), int(x)] = True
    yy, xx = np.mgrid[:21, :21]
    circ = np.hypot(xx - tg["x"], yy - tg["y"]) <= 2.0
    mT = prf.locate(tg["x"], tg["y"], shape)
    per_sector[sec] = dict(target_frac_in_aperture=float(mT[apmask].sum()), target_frac_in_r2=float(mT[circ].sum()),
                           aperture_npix=int(apmask.sum()), aperture_pixels_inside_stamp=int(apmask.sum()) == len(ap_abs))
    for s in sc:
        if s["is_target"]:
            continue
        mN = prf.locate(s["x"], s["y"], shape) * s["rel_flux"]
        fr_bright = s["rel_flux"] * (10 ** (0.4 * s["T_err"]) if s["T_flag"] != "ok" else 1.0)
        rec = dict(source_id=s["source_id"], sector=sec, sep=s["sep"], T=s["T"], T_err=s["T_err"], T_flag=s["T_flag"],
                   rel_flux=s["rel_flux"], rel_flux_bright=fr_bright,
                   frac_in_ap=float(mN[apmask].sum() / s["rel_flux"]) if s["rel_flux"] > 0 else 0.0)
        for sig, d in DEPTH.items():
            rec[f"req_{sig}_full"] = float(d * (1 + fr_bright) / fr_bright)
            fN = mN[apmask].sum() * (fr_bright / s["rel_flux"])
            rec[f"req_{sig}_aperture"] = float(d * (mT[apmask].sum() + fN) / fN) if fN > 0 else np.inf
            fN2 = mN[circ].sum() * (fr_bright / s["rel_flux"])
            rec[f"req_{sig}_r2"] = float(d * (mT[circ].sum() + fN2) / fN2) if fN2 > 0 else np.inf
        rows.append(rec)
out["neighbour_census"] = dict(n_sources_T_lt_22_within_180=len({r["source_id"] for r in rows}), per_sector=per_sector)
for sig in ("A", "B"):
    cap_full = sorted({(r["source_id"], round(r["sep"], 1)) for r in rows if r[f"req_{sig}_full"] <= 1.0}, key=lambda x: x[1])
    cap_ap = sorted({(r["source_id"], round(r["sep"], 1), r["sector"]) for r in rows if r[f"req_{sig}_aperture"] <= 1.0}, key=lambda x: x[1])
    out["neighbour_census"][sig] = dict(capable_full_inclusion=cap_full, capable_in_actual_aperture=cap_ap,
                                        nearest_capable_full_inclusion_arcsec=cap_full[0][1] if cap_full else None)
    print(sig, "capable (full inclusion):", cap_full[:8])
    print(sig, "capable in actual apertures (source, sep, sector):", cap_ap[:12])
save_json(rows, "out/s14_neighbour_census_rows.json")
close = sorted({(r["source_id"], round(r["sep"], 1), round(r["T"], 2), r["T_flag"], round(r["req_A_full"], 2), round(r["req_B_full"], 2))
                for r in rows if r["sep"] < 60}, key=lambda x: x[1])
out["neighbours_within_60"] = close
for c in close:
    print("  ", c)

# ---------------- (2) RUWE percentile among comparable stars
q = ("SELECT ruwe, phot_g_mean_mag, bp_rp, parallax FROM gaiadr3.gaia_source WHERE phot_g_mean_mag BETWEEN 10.0 AND 11.0 "
     "AND bp_rp BETWEEN 1.60 AND 1.90 AND parallax > 15 AND astrometric_params_solved = 31 AND ruwe IS NOT NULL")
data = urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q}).encode()
fetch("https://gea.esac.esa.int/tap-server/tap/sync", "data/gaia/ruwe_comparison_sample.csv", note="RUWE comparison sample", data=data)
open(os.path.join(AUD, "data/gaia/ruwe_comparison_sample.adql"), "w").write(q + "\n")
R = np.array([float(r["ruwe"]) for r in csv.DictReader(open(os.path.join(AUD, "data/gaia/ruwe_comparison_sample.csv")))])
out["ruwe"] = dict(target=1.2142985, n_comparison=int(len(R)), percentile_of_target=float(np.mean(R < 1.2142985) * 100),
                   comparison_median=float(np.median(R)), comparison_p84=float(np.percentile(R, 84)), comparison_p95=float(np.percentile(R, 95)),
                   frac_comparison_gt_1p4=float(np.mean(R > 1.4)),
                   query=q)
print("RUWE:", out["ruwe"])

# ---------------- (3) density-host constraint vs overluminosity
pm = []
for line in open(os.path.join(AUD, "data/refs/EEM_dwarf_UBVIJHK_colors_Teff.txt")):
    if not line.strip() or line.startswith("#"):
        continue
    p = line.split()
    if not p[0].endswith("V") or len(p) < 30:
        continue
    try:
        spt, teff, R_, MG, Msun = p[0], float(p[1]), float(p[6]), float(p[13]) if p[13] != "..." else np.nan, float(p[30])  # column 30 = Msun
    except (ValueError, IndexError):
        continue
    pm.append((spt, teff, R_, MG, Msun, Msun / R_ ** 3))
    if spt == "M9.5V":
        break
fA = json.load(open(os.path.join(AUD, "out/fits/fit_A_pipe_free_frozen_w0.75.json")))
rho_up_logprior = fA["rho"][4]
rho_up_uniform = fA["rho_uniform_prior_reweighted"][4]
allowed = [(s, M, round(r, 2)) for s, T, R_, MG, M, r in pm if r <= max(rho_up_logprior, rho_up_uniform) and 0.08 <= M <= 0.8]
way = {}
for fn in glob.glob(os.path.join(OLD_COMP, "data", "way2026_*row.csv")):
    way[os.path.basename(fn)] = list(csv.DictReader(open(fn)))
out["host_density_constraint_A"] = dict(rho_97p7_logprior=rho_up_logprior, rho_97p7_uniformprior=rho_up_uniform,
                                        dwarf_types_allowed_circular=allowed,
                                        min_host_mass_allowed=min(M for s, M, r in allowed) if allowed else None,
                                        way2026_rows=way)
print("host density constraint A: allowed", allowed[:3], "... min mass", out["host_density_constraint_A"]["min_host_mass_allowed"])

# ---------------- (4) WISE excess
tm = {"J": 8.527, "H": 7.932, "Ks": 7.694}
wise = {"W1": (7.605, 0.031), "W2": (7.626, 0.020), "W3": (7.560, 0.019), "W4": (7.451, 0.133)}
k8 = [x for x in open(os.path.join(AUD, "data/refs/EEM_dwarf_UBVIJHK_colors_Teff.txt")) if x.startswith("K8V ")][0].split()
m0 = [x for x in open(os.path.join(AUD, "data/refs/EEM_dwarf_UBVIJHK_colors_Teff.txt")) if x.startswith("M0V ")][0].split()
out["wise"] = dict(Ks_minus_W1=tm["Ks"] - wise["W1"][0], Ks_minus_W2=tm["Ks"] - wise["W2"][0], Ks_minus_W3=tm["Ks"] - wise["W3"][0],
                   Ks_minus_W4=tm["Ks"] - wise["W4"][0], W4_err=wise["W4"][1], PM_K8V_Ks_minus_W1=float(k8[23]), PM_M0V_Ks_minus_W1=float(m0[23]),
                   note="PM table gives W1-W3/W1-W4 only for earlier types; for K7-M0 dwarfs photospheric Ks-W3 ~0.1-0.2 and Ks-W4 ~0.2-0.3 mag")
print("WISE:", out["wise"])

# ---------------- (5) Gaia RV constraint
g = list(csv.DictReader(open(os.path.join(AUD, "data/gaia/gaia_dr3_source.csv"))))[0]
N = int(g["rv_nb_transits"])
erv = float(g["radial_velocity_error"])
sig_epoch = np.sqrt(max(erv ** 2 - 0.11 ** 2, 0) / (np.pi / (2 * N)))   # Katz+2023 eq. for rv_method_used = 1
amp = float(g["rv_amplitude_robust"])
Ms = 0.606
periods = [3.18, 11.13, 30, 100, 300, 1000]
Klim = amp  # conservative: a full peak-to-peak amplitude as large as the robust amplitude -> K = amp/2; use K_lim = amp


def mc_min(P, Klim_kms):
    # solve K = (2 pi G / P)^(1/3) Mc sin i / (M* + Mc)^(2/3), sin i = 1
    Gc, Msun_kg = 6.674e-11, 1.989e30
    Ps = P * 86400
    lo, hi = 1e-4, 2.0
    for _ in range(100):
        m = 0.5 * (lo + hi)
        K = (2 * np.pi * Gc / Ps) ** (1 / 3) * m * Msun_kg / ((Ms + m) * Msun_kg) ** (2 / 3) / 1e3
        lo, hi = (m, hi) if K < Klim_kms else (lo, m)
    return 0.5 * (lo + hi)


out["gaia_rv"] = dict(rv=float(g["radial_velocity"]), rv_error=erv, n_transits=N, rv_method_used=int(g["rv_method_used"]),
                      epoch_scatter_kms_from_error=float(sig_epoch), rv_amplitude_robust=amp, rv_chisq_pvalue=float(g["rv_chisq_pvalue"]),
                      rv_renormalised_gof=float(g["rv_renormalised_gof"]),
                      katz2023_variability_criteria="rv_chisq_pvalue <= 0.01 AND rv_renormalised_gof > 4 AND rv_nb_transits >= 10 -> NOT met (gof 3.0)",
                      Klim_kms_used=Klim, min_companion_mass_Msun_sin_i_1={str(P): round(mc_min(P, Klim), 3) for P in periods},
                      note="K_lim = rv_amplitude_robust (peak-to-peak) is conservative; companions above these masses at these periods would "
                           "produce larger RV spreads than observed (approximate, sin i = 1, circular)")
print("Gaia RV:", out["gaia_rv"])
save_json(out, "out/s14_companions.json")
