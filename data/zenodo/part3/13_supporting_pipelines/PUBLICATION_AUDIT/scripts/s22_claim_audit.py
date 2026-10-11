"""Step 22: CLAIM_AUDIT.csv - one row per candidate manuscript claim (RNAAS draft + deep follow-up + FPP report + companion
check + archival checks + priority re-check). recalculated values are read from FINAL_PARAMETERS.json / out/*.json where
possible; the FPP rows are filled from out/fpp/runs/*/summary.json at run time."""
import csv
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from audit_common import AUD  # noqa

F = json.load(open(os.path.join(AUD, "FINAL_PARAMETERS.json")))
A, B = F["signal_A"], F["signal_B"]


def fpp(sig, var, fold="S32+S98"):
    p = os.path.join(AUD, "out", "fpp", "runs", f"{sig}_{var}_{fold}", "summary.json")
    if not os.path.exists(p):
        return "not run"
    L = json.load(open(p))["log"]
    if "FPP_mean" not in L:
        return "incomplete"
    return f"{L['FPP_mean']:.4f} ± {L['FPP_std']:.4f} (sd of {len(json.load(open(p))['runs'])} runs; range {L['FPP_min']:.4f}-{L['FPP_max']:.4f}); NFPP {L['NFPP_mean']:.1e}"


def v(x, fmt="{:.4g}"):
    return f"{fmt.format(x['value'])} +{fmt.format(x['plus'])}/-{fmt.format(x['minus'])}"


R = []
add = lambda *a: R.append(dict(zip(["id", "claim", "source_data", "source_code", "recalculated_value", "previous_value", "difference",
                                    "uncertainty", "status", "allowed_wording"], a)))
S = F["star"]
ad = S["adopted"]
add("S1", "Distance", "Gaia DR3 parallax (+Lindegren+2021 ZP)", "s04_stellar.py", f"{S['astrometry']['distance_pc']:.2f} pc (ZP-corrected); 35.62 pc uncorrected; Bailer-Jones 35.58",
    "35.62 ± 0.02 pc", "0.04 pc (ZP)", "±0.03 pc", "VERIFIED", "at 35.6 pc (Gaia DR3)")
add("S2", "TESS magnitude T = 9.68", "TIC 8.2; Gaia G, BP-RP", "s04_stellar.py", f"{S['photometry']['T_TIC82']} (TIC); {S['photometry']['T_stassun_from_gaia']:.3f} (Stassun eq.)", "9.68", "0", "0.006", "VERIFIED", "T = 9.68")
add("S3", "Spectral class 'K dwarf' / K7-M0", "SIMBAD/LAMOST/Lepine&Gaidos; M_Ks vs Pecaut&Mamajek", "s04, s14", "K7V-M0V", "K6 / K7 / dM1; 'K dwarf'", "-", "-", "VERIFIED",
    "a late-K (K7-M0) dwarf")
add("S4", "M_Ks = 4.94", "2MASS Ks, Gaia parallax", "s04", f"{S['MK']['median']:.3f} ± {S['MK']['std']:.3f}", "4.936", "+0.002 (ZP)", "0.029", "VERIFIED", "M_Ks = 4.94")
add("S5", "R* = 0.623 ± 0.019 R_sun (Mann+2015)", "M_Ks; Mann+2015 Table 1 / erratum 2016", "s04", f"{ad['R_Rsun'][0]:.3f} ± {ad['R_Rsun'][1]:.3f}", "0.623 ± 0.019", "0", "3.0 %",
    "VERIFIED", "R* = 0.623 ± 0.019 R_sun (Mann et al. 2015, 2016 erratum); SED 0.621 ± 0.015 agrees")
add("S6", "M* = 0.606 ± 0.013 M_sun (Mann+2019)", "M_Ks; Mann+2019 Table 6/7", "s04", f"{ad['M_Msun'][0]:.3f} ± {ad['M_Msun'][1]:.3f}", "0.606 ± 0.013", "error +0.002",
    "2.5 %", "CORRECTED", "M* = 0.606 ± 0.015 M_sun (Mann et al. 2019)")
add("S7", "rho* = 2.51 ± 0.23 rho_sun", "R*, M*", "s04", f"{ad['rho_rhosun']['value']:.2f} +{ad['rho_rhosun']['plus']:.2f}/-{ad['rho_rhosun']['minus']:.2f}", "2.51 ± 0.23", "0",
    "9 %", "VERIFIED", "rho* = 2.51 +0.25/-0.22 rho_sun")
add("S8", "Teff 4010 K (TIC) / adopted 4080 ± 100 K", "TIC, Mann V-J, ARIADNE, Gaia DR2, GSP-Spec", "s04, s20", f"{ad['Teff_K'][0]} ± {ad['Teff_K'][1]:.0f} K", "4080 ± 100 (adopted)",
    "-50 K", "90 K", "CORRECTED", "Teff = 4030 ± 90 K (mean of TIC, Mann V-J, SED, Gaia DR2; GSP-Spec excluded as flagged)")
add("S9", "[M/H] = -0.27 (Gaia GSP-Spec)", "Gaia DR3 AP table", "s03", "-0.27 with extrapolation flag; SED -0.08 ± 0.18", "-0.27", "-", "unknown", "UNRESOLVED",
    "metallicity poorly constrained (-0.3 to 0)")
add("S10", "RV 10.49 ± 0.22 km/s; not RV-variable", "Gaia DR3", "s14", "10.49 ± 0.22; p = 4.1e-5, gof 3.0 (< 4), N = 33; epoch scatter 0.88 km/s", "same", "0", "-",
    "VERIFIED", "Gaia DR3 RVs do not meet the variability criterion (gof 3.0 < 4), although the chi2 p-value is small (4e-5)")
add("S11", "RUWE 1.21, no NSS solution, ipd_frac_multi_peak 0", "Gaia DR3", "s03, s14", f"RUWE 1.214 = {F['star']['binarity_and_companions']['ruwe']['percentile_of_target']:.0f}th percentile of 828 comparable stars",
    "1.21 'normal'", "-", "-", "VERIFIED", "RUWE 1.21 (typical for comparable nearby stars); no Gaia NSS solution")
add("S12", "No comoving star within 1 deg", "Gaia DR3 TAP", "s03", "0 stars with plx 24-32 mas within 1 deg (other than target)", "none", "0", "-", "VERIFIED", "no comoving Gaia DR3 source within 1 deg")
add("S13", "Rotation period unknown", "TESS SAP, ASAS-SN", "s15", "not detected; catalogue ~8 d not confirmed", "UNKNOWN (deep follow-up); ~8 d in 4 catalogues (priority re-check)", "-", "-",
    "VERIFIED", "no rotation period is detected in TESS or ASAS-SN photometry")
add("S14", "Quadratic LD u1 = 0.351, u2 = 0.289", "Claret 2017 table15", "s08 transit_tools.claret_ld", "0.33, 0.30 at 4030 K, log g 4.63 (table values verified)", "0.351, 0.289 (4080 K)",
    "-0.02, +0.01", "prior sigma 0.1", "CORRECTED", "Gaussian priors (sigma 0.1) centred on Claret (2017) PHOENIX-COND values")
add("S15", "X-ray source (new)", "eRASS1 (Merloni+2024)", "s15", "1eRASS J052326.6-011938, 4.8'' away, DET_LIKE 8.0, log Lx/Lbol ~ -4.8 to -4.9", "not reported", "new", "±0.2 dex+",
    "VERIFIED", "a faint eRASS1 X-ray source (detection likelihood 8) coincides with the star")
add("A1", "Signal A detected in the S98 QLP BLS search, MES 10.6", "frozen search candidate.json", "LEO-vetter 1.2.0 main.py", "MES = dep/sqrt(sig_w^2/n_in + sig_r^2/N_tr) = 10.6 (S98 only)",
    "10.6", "0", "-", "VERIFIED", "signal-to-pink-noise ratio 10.6 in Sector 98 (LEO-vetter MES)")
add("A2", "Signal A recurs in S6 and S32", "SPOC/TESS-SPOC/TESScut", "s07, s09", "independent BLS per sector: S6 3.1826 d (SDE 8.7), S32 3.1825 d (SDE 7.6); combined blind BLS 3.18279 d", "yes",
    "-", "-", "VERIFIED", "recovered independently in each of Sectors 6, 32 and 98")
add("A3", "23 transits (6, 8, 9 per sector)", "pipeline set", "s09", "6 (S6), 8 (S32), 9 (S98 clean) = 23 with >=50 % coverage", "23 (6, 8, 9)", "0", "-", "VERIFIED", "23 transits")
add("A4", "Depth 509 ± 28 ppm", "pipeline set", "s09, s17, s08", f"mid-transit (template) 503 ± 29 ppm; k^2 = {v(A['shape_mann_prior']['depth_k2_ppm'],'{:.0f}')} ppm; box (1.87 h) 470 ± 32", "509 ± 28", "-6 ppm",
    "6 %", "VERIFIED", "transit depth ~500 ppm (k^2 = 448 ± 30 ppm); state which definition")
add("A5", "Rp = 1.44 +0.06/-0.07 R_earth", "joint fit, Mann rho prior", "s08", f"{v(A['shape_mann_prior']['Rp_Rearth'],'{:.3f}')} (Mann prior); free rho {v(A['free_density']['Rp_Rearth'],'{:.2f}')}; free LD 1.50 +0.09/-0.08",
    "1.44 +0.06/-0.07", "0", "4.5 %", "VERIFIED", "Rp = 1.44 ± 0.06 R_earth (1.46-1.50 without the density prior or with free limb darkening)")
add("A6", "P = 3.1827863 ± 0.0000028 d", "joint fit", "s08", f"{A['ephemeris']['P_d']['value']:.7f} ± {A['ephemeris']['P_d']['plus']:.7f}", "3.1827863 ± 0.0000028", "-1.0e-6 (0.4 sigma)", "2.7e-6",
    "CORRECTED", "P = 3.1827853 ± 0.0000027 d")
add("A7", "T0 = 2460998.7368 ± 0.0012 BJD_TDB", "joint fit", "s08", f"{A['ephemeris']['T0_BJD_TDB']:.5f} ± 0.00112", "2460998.7368 ± 0.0012", "-0.00001", "0.0011", "VERIFIED",
    "T0 = 2460998.7368 ± 0.0011 BJD_TDB")
add("A8", "rho* from transit (free) = 2.7 +0.5/-1.0 rho_sun, consistent with catalogue", "joint fit, free rho", "s08", f"uniform prior {v(A['free_density']['rho_uniformprior'],'{:.2f}')}; ln-uniform {v(A['free_density']['rho_logprior'],'{:.2f}')}",
    "2.68 +0.51/-1.03", "0.07", "prior-dependent lower side", "VERIFIED", "consistent with the catalogue density; the lower bound depends on the prior, the 97.7 % upper limit is 3.5 rho_sun")
add("A9", "rho* > 3.7 rho_sun disfavoured (Delta chi2 > 4)", "profile likelihood (white noise)", "s08 MCMC", "97.7 % upper limit 3.47 (ln prior) / 3.55 (uniform)", "3.7", "-0.2", "-", "VERIFIED",
    "rho* < 3.5 rho_sun (97.7 %, circular orbit)")
add("A10", "b = 0.36 +0.09/-0.13", "joint fit, Mann prior", "s08", v(A["shape_mann_prior"]["b"], "{:.2f}"), "0.36 +0.09/-0.13", "-0.01", "-", "VERIFIED", "b = 0.35 +0.10/-0.14 (with the stellar-density prior; prior-driven)")
add("A11", "T14 = 1.86 ± 0.05 h", "joint fit", "s08", v(A["shape_mann_prior"]["T14_h"], "{:.3f}"), "1.86 ± 0.05", "+0.01", "0.05", "VERIFIED", "T14 = 1.87 ± 0.05 h")
add("A12", "a = 0.036 AU", "P, M*", "s20", v(A["derived"]["a_AU"], "{:.4f}"), "0.0358", "0", "-", "VERIFIED", "a = 0.0358 AU")
add("A13", "Teq 820 K, S 75 S_earth", "R*, Teff, a", "s20", f"Teq {v(A['derived']['Teq_K_A0_full_redistribution'],'{:.0f}')} K; S {v(A['derived']['S_Searth'],'{:.0f}')}", "820 ± 13 K; 75 ± 5",
    "-10 K; -3", "±22 K; ±8", "CORRECTED", "Teq ~810 K (A_B = 0, full redistribution); S ~72 S_earth")
add("A14", "Per-sector depths 508 ± 58, 463 ± 47, 550 ± 44 ppm", "pipeline set", "s09, s17", "S6 521 ± 53, S32 476 ± 52, S98 511 ± 47 (per-event template, red-noise errors); 483-558 across detrendings",
    "508/463/550", "<= 40 ppm", "~50 ppm", "VERIFIED", "per-sector depths agree within 1 sigma")
add("A15", "Odd/even 550 ± 39 vs 463 ± 41 ppm (1.5 sigma)", "per-event depths", "s09, s10, s10b", "per-event: 536 ± 40 vs 466 ± 43 (1.2 sigma, pipeline); 530 vs 530 (0.0 sigma, own, matched events); box test 2.2 sigma with opposite sign in the independent photometry",
    "1.5 sigma", "-0.3 sigma", "-", "CORRECTED", "odd and even depths agree (per-event fits within 1.2 sigma; box comparisons give +2.2 and -2.1 sigma in two photometries, i.e. opposite signs)")
add("A16", "Flux at phase 0.5 = -22 ± 24 ppm", "pipeline set", "s10", "-20 ± 32 ppm (empirical red-noise error); 3 sigma upper limit 95 ppm", "-22 ± 24", "0", "32 (was 24, white)", "CORRECTED",
    "no secondary eclipse (-20 ± 32 ppm at phase 0.5; < 95 ppm at 3 sigma)")
add("A17", "Timing scatter up to ~20 min, chi2 = 120 for 21 dof", "per-event times", "s09", "chi2 47.5/19 (pipeline, red-noise-inflated errors), rms O-C 9.5 min; own set 42.6/22", "120/21 (RNAAS) vs 114/21 (report)",
    "inconsistent before", "-", "CORRECTED", "O-C rms ~10 min, consistent with underestimated single-transit timing errors; no TTV claim")
add("A18", "Median single-transit S/N 3.8", "per-event fits", "s09", "3.5 (pipeline set, red-noise errors)", "3.8", "-0.3", "-", "CORRECTED", "median single-transit S/N ~3.5")
add("A19", "No Gaia DR3 source within 43'' can produce the depth", "Gaia DR3 census + PRF apertures", "s14", "full inclusion: nearest capable 43.5'' (needs 45 %); in real apertures also a T=13.06 star at 83.8'' could in S32 (28 %)",
    "43'' (five stars T 16.3-17.3 at 43-56'')", "new star at 83.8''", "-", "CORRECTED", "no resolved star closer than 43'' can produce it; all capable stars (43-84'') are excluded by the pixel localisation")
add("A20", "Difference images from all three sectors peak on the target", "TESScut, TPF", "s12, s12c, s12d", "S32 2-min and S98 clean: consistent with target; combined (+2.5, -7.5)'' ± 5.7''; S6 compromised by momentum dumps",
    "all three on target", "S6 unusable", "±5.7''", "CORRECTED", "difference images in S32 and S98 place the source within ~6'' (1 sigma) of the target")
add("A21", "Offsets 1.6'' (S6, 0.2 sigma), 10.5'' (S32, 3.8 sigma), 19'' (S98, 1.6 sigma) with formal errors", "shifted-OOT fits + bootstrap", "prf_shift_test.py", "superseded by PRF fits with injection-calibrated errors",
    "as quoted", "-", "-", "SUPERSEDED", "do not quote")
add("A22", "POS_CORR1 moves -0.16 ± 0.04 mpx in transit; residual 3 ± 3''", "SPOC 2-min", "s12c", "MOM_CENTR1 -0.24, POS_CORR1 -0.17, residual -0.06 mpx; empirical null sd 0.14/0.14/0.10 mpx -> z -1.7/-1.1/-0.9",
    "4 sigma pointing shift", "not significant", "-", "CORRECTED", "no significant in-transit centroid or pointing shift")
add("A23", "FPP = 0.031 ± 0.003 (default)", "TRICERATOPS 1.1.0", "s13", fpp("A", "repro") + " | exptime 6 min: " + fpp("A", "base"), "0.031 ± 0.003", "0", "run-to-run sd 0.003", "VERIFIED",
    "TRICERATOPS FPP = 0.031 ± 0.003 (default assumptions, 2.2'' limit)")
add("A24", "NFPP < 2e-12", "TRICERATOPS", "s13", "1.5-2.4e-12 (0 after removing TIC 4206077)", "<2e-12", "0", "-", "VERIFIED",
    "NFPP ~ 1e-12 (note: TRICERATOPS's Gaussian PSF ignores the 83.8'' T=13.1 star; it is excluded by the pixel data)")
add("A25", "FPP = 0.005 ± 0.003 with MOLUSC", "TRICERATOPS + MOLUSC", "s13", fpp("A", "molusc5") + " | age 1 Gyr: " + fpp("A", "molusc1"), "0.0052 ± 0.0034", "see value", "-", "VERIFIED",
    "with MOLUSC companion limits FPP is lower, but this mode drops companions within 10 au and the companion prior (tool choice)")
add("A26", "FPP scenario: nearly all is STP", "TRICERATOPS scenario tables", "s13", "STP 0.031 of 0.031; EB scenarios ~1e-7", "same", "0", "-", "VERIFIED", "almost the entire FPP is a planet on an unresolved bound companion")
add("A27", "EB on the target excluded", "Gaia RV", "s14", "K < 3.2 km/s -> companion M sin i < 0.016 M_sun at 3.18 d, < 0.02 at 6.37 d", "excluded (K >> 3 km/s)", "-", "-", "VERIFIED", "an eclipsing stellar companion of the target is excluded by the Gaia RVs")
add("A28", "U-shaped (not V)", "folded light curve", "s10", "trapezoid: V-shape disfavoured by Delta chi2 = 31.8 (pipeline) / 31.4 (own); ingress fraction < 0.18 (2 sigma)", "V not seen", "-", "-", "VERIFIED", "the transit is flat-bottomed")
add("A29", "Planet on a bound companion (STP) not excluded", "density + overluminosity", "s14", "circular-orbit STP excluded (host must be >= 0.54 M_sun but O = +0.003 mag); eccentric STP open", "not excluded", "partly closed", "-",
    "CORRECTED", "a planet orbiting a bound companion is disfavoured unless the orbit is highly eccentric")
add("A30", "No SPOC TCE because S32 MES 6.7 < 7.1", "estimate", "s09", "S32-only S/N 6.7-9.1 depending on estimator", "6.7 < 7.1", "-", "-", "UNRESOLVED", "no SPOC TCE was issued (reason not established)")
add("A31", "Multiplicity boost not applied (Guerrero+2021 not valid for QLP-only)", "Guerrero+2021 text", "-", "verified in the paper text (factor ~20; ~54 < 6 R_earth; 'do not apply to ... QLP')", "same", "0", "-", "VERIFIED",
    "we do not apply a multiplicity boost")
add("B1", "Signal B detected in S98, MES 11.0", "frozen search; QLP S98 (three masks), TESScut S98", "LEO-vetter; s07, s07b, s07c, s11c",
    "LEO MES 11.0 (S98 only) reproduced, but it is not a calibrated significance. Calibrated BLS test (null: segment-scrambled S98 data, B removed): "
    "FAP 0.33 in the discovery (frozen-mask) data, 0.008 in own photometry, ~1 with the QLP default mask. The B-free S98 data contain a peak as strong as B "
    "(SDE 6.48 vs 6.09). Single S98 B events (box S/N 7.3, 7.9) are matched by a non-B dip of S/N 8.7 at BTJD 4017.68", "11.0", "0", "FAP 0.008-1", "CORRECTED",
    "a candidate signal found in Sector 98 (LEO-vetter signal-to-pink-noise ratio 11.0); Sector 98 alone does not establish it (calibrated false-alarm probability 0.008-0.33 depending on the photometry, ~1 when the scattered-light intervals are kept)")
add("B2", "S98: two complete transits (936 ± 109, 1062 ± 115 ppm) and one partial", "QLP frozen mask; TESScut", "s09, s09b", "clean: 1102 ± 132, 1245 ± 139, partial 939 ± 242 (QLP); TESScut shows epoch 4 complete (917 ± 196) and epoch 2 (1198 ± 198) in a QLP-flagged, SPOC-unflagged interval",
    "2 + partial", "+1 transit (epoch 2)", "~12 %", "CORRECTED", "Sector 98 contains three clean transits and a fourth in a scattered-light-affected interval")
add("B3", "BLS of S6+S32 alone does not return B as strongest", "S6, S32", "s07", "S6 alone: B/2 at 5.557 d is the 2nd peak after A; S32: B not in top peaks", "same", "-", "-", "VERIFIED", "archival sectors alone do not reveal B")
add("B4", "Link S/N 9.9 at P = 11.13274 d exceeds all 300 epoch-shifted trials (max 7.5)", "S6+S32", "s11, s11b, s11c",
    "statistic reproduced (10.4 at 11.132736 d). The epoch-shifted null at the fixed S98 period is invalid: its trials overlap (16 distinct windows per 1500 trials in S6; "
    "p99 = max). Period+epoch-randomised null (20,000 trials, >2,600 distinct windows): S6+S32 FAP 1.2e-3 (24/20,000; tail set by S32 edge artefacts); "
    "S6 alone 0/18,652 (null max 5.3 vs 7.6); S32 alone FAP 0.03. Rank test: the two strongest 1.6-h dips of all of S6 (outside A) are the two predicted B events "
    "(pipeline 6.2, 4.3 vs next 3.2; own photometry 5.6, 6.8 vs 3.7): distribution-free p = 1.6e-3. Gaussian noise model (S6 per-position S/N sd 1.02): ~1e-12",
    "FAP < 1/300 (null max 7.5)", "FAP ~1e-3 (distribution-free), not < 5e-5", "model-dependent (1e-3 to 1e-12)", "CORRECTED",
    "in Sector 6 the two strongest dips of the whole sector fall at the times predicted by the Sector 98 ephemeris (distribution-free false-alarm probability ~2e-3, about 3 sigma); Sector 32 alone is not significant")
add("B5", "Joint P = 11.132743 ± 0.000013 d, T0 = 2461000.8661 ± 0.0016", "joint fit", "s08", f"P {B['ephemeris']['P_d']['value']:.7f} ± {B['ephemeris']['P_d']['plus']:.7f}; T0 {B['ephemeris']['T0_BJD_TDB']:.5f} ± 0.0016",
    "11.132743 / 2461000.8661", "+2.5e-6 / -0.0002", "1.3e-5 / 0.0016", "VERIFIED", "P = 11.132745 ± 0.000013 d, T0 = 2461000.8659 ± 0.0016 BJD_TDB")
add("B6", "Archival depths 731 ± 122, 598 ± 104, 525 ± 131, 728 ± 161 ppm", "box depths", "s09", "template: S6 924 ± 153, 796 ± 134; S32 720 ± 157, 1055 ± 204 (pipeline)", "box values", "definition", "-",
    "CORRECTED", "quote per-event depths with one stated definition")
add("B7", "S98 and archival depths differ by 3.1 sigma", "per-event fits", "s09, s17", "2.2-2.7 sigma across six detrendings / two photometries", "3.1 sigma (RNAAS), 2.4 sigma (report)", "-0.4 to -0.9 sigma", "-",
    "CORRECTED", "the Sector 98 depth is larger than the archival depth at the 2-3 sigma level")
add("B8", "Phase 0.5 flux 85 ± 45 ppm rules out P = 5.57 d", "all sectors", "s10", "61 ± 62 ppm (empirical error); a 5.57-d period would need ~1000 ppm", "85 ± 45", "-", "62", "VERIFIED", "no dip at phase 0.5, excluding half the period")
add("B9", "Odd/even agree within 1.2 sigma", "per-event", "s10, s10b", "per-event 1.5 sigma (2 odd vs 5 even events); box 0.3 sigma", "1.2 sigma", "-", "-", "VERIFIED", "odd and even depths agree within 1.5 sigma")
add("B10", "Rp ~ 2.0 R_earth (68 % 1.9-2.2)", "joint fit", "s08, s08b", f"{v(B['shape_mann_prior']['Rp_Rearth'],'{:.2f}')} with the target's density (b = 0.87); free rho {v(B['free_density']['Rp_Rearth'],'{:.2f}')}; free LD 2.05 +0.26/-0.16",
    "2.0 (1.9-2.2)", "+0.25 (prior)", "prior-dependent", "CORRECTED", "Rp = 2.0-2.4 R_earth depending on the unconstrained impact parameter (2.25 ± 0.12 if B orbits the target)")
add("B11", "Fit does not constrain b or rho*", "joint fit", "s08", f"free rho {v(B['free_density']['rho_logprior'],'{:.1f}')} (ln prior) vs {v(B['free_density']['rho_uniformprior'],'{:.1f}')} (uniform)", "same", "-", "-", "VERIFIED",
    "the light curve does not constrain the impact parameter or the stellar density")
add("B12", "a = 0.083 AU, Teq ~540 K, S ~14", "P, M*, R*, Teff", "s20", f"a {v(B['derived']['a_AU'],'{:.4f}')}; Teq {v(B['derived']['Teq_K_A0_full_redistribution'],'{:.0f}')}; S {v(B['derived']['S_Searth'],'{:.1f}')}", "0.083 / 540 / 14", "-7 K", "-",
    "VERIFIED", "a = 0.083 AU, Teq ~530 K, S ~13.5 S_earth")
add("B13", "FPP = 0.146 ± 0.011 (default)", "TRICERATOPS 1.1.0", "s13, s13c", fpp("B", "repro") + " | exptime 6 min: " + fpp("B", "base") + " | TIC params: " + fpp("B", "tic")
    + " | earlier fold + earlier configuration: " + fpp("B", "prevfold_repro") + " (reproduces the earlier run; max per-seed difference 7e-5) | earlier fold + audit configuration: " + fpp("B", "prevfold_base"),
    "0.146 ± 0.011", "+0.03 (fold only)", "sd 0.011; fold- and host-dependent", "CORRECTED",
    "TRICERATOPS FPP 0.13-0.18 under default assumptions (fold- and stellar-parameter-dependent); not validated; the FPP assumes the signal is astrophysical")
add("B14", "FPP 0.05 ± 0.02 with MOLUSC", "TRICERATOPS + MOLUSC", "s13", fpp("B", "molusc5") + " | age 1 Gyr: " + fpp("B", "molusc1"), "0.0505 ± 0.0165", "see value", "-", "VERIFIED",
    "lower with MOLUSC limits, subject to the same tool caveats as A25")
add("B15", "B pixel offsets within 1.3 sigma of injections", "TESScut, TPF", "s12c", "PRF fits: combined (+1.0, +1.4)'' ± (2.3, 3.0)''; without S6 (+6.1, +1.6)'' ± (3.6, 4.0)''", "0.5/5.8/5.6'' within 1.3 sigma", "method changed", "±3-4''",
    "CORRECTED", "the Sector 6, 32 and 98 difference images of B are consistent with the target")
add("B16", "Period ratio 11.1327/3.1828 = 3.498 (0.06 % from 7:2)", "joint ephemerides", "s20", f"{F['period_ratio_B_over_A']['value']:.5f} ({F['period_ratio_B_over_A']['fractional_offset_from_7_2']*100:.3f} %)", "3.498 (0.06 %)", "0", "-", "VERIFIED",
    "period ratio 3.4977, 0.06 % inside the 7:2 commensurability")
add("O1", "TESS Sector 110 (2026-10-31 to 11-26) will observe the star", "HEASARC + MIT pointing, tess-point", "s16", "camera 1 CCD 2 (col 702.8, row 1577.0); 8-9 A and 3 B transits inside the sector dates", "S110, cam1/CCD2, 8 A, 3 B (±23 min)",
    "B timing ±2.7 min if the link holds", "-", "VERIFIED", "TESS Sector 110 (2026 Oct 31 - Nov 26) covers the star on camera 1")
add("O2", "No prior public report found", "sealed priority re-check + refresh", "s18", "negative in all searched sources at 04:46-05:25Z; TOI/CTOI/TCE/NASA re-checked 06:30Z", "same", "0", "-", "VERIFIED",
    "no prior public report was identified in the searched sources (as of 2026-09-25 UTC)")
add("O3", "Instrumental origin excluded (three sectors, two CCDs, several pipelines)", "all", "s07, s07b, s09, s11b, s11c, s17",
    "A recovered independently in each sector, in SPOC, TESS-SPOC, QLP and own photometry. B: S98 detection mask-dependent; S6 support ~3 sigma distribution-free (two photometries); S32 not significant",
    "excluded", "weaker for B", "-", "CORRECTED",
    "A persists across three sectors, two CCDs and independent reductions of the pixels; B is supported by Sector 6 at about 3 sigma")
add("O4", "Two-planet candidate system", "all", "-", "A: strong candidate; B: low-significance candidate (S98 alone FAP 0.008-1; S6 rank test 1.6e-3; overall 2.3-3.9 sigma; density unconstrained; FPP 0.13-0.18)", "two-planet candidate system", "-", "-",
    "CORRECTED", "one transiting planet candidate (3.18 d) and a second, low-significance candidate signal (11.13 d) that TESS Sector 110 can test")
add("B17", "(implied) B is a secure detection", "S98, S6, S32", "s07b, s11b, s11c, s11d", "indicative overall FAP 5e-5 to 1e-2 (Fisher: S98 test x S6 rank test x S32 null), i.e. 3.9 to 2.3 sigma depending on the S98 product",
    "FAP < 1/300 (link) plus S98 detection", "large", "model/product dependent", "CORRECTED", "B is significant at about 3 sigma overall; confirmation requires new data (Sector 110)")
keys = ["id", "claim", "source_data", "source_code", "recalculated_value", "previous_value", "difference", "uncertainty", "status", "allowed_wording"]
with open(os.path.join(AUD, "CLAIM_AUDIT.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=keys)
    w.writeheader()
    w.writerows(R)
from collections import Counter  # noqa
print(len(R), "claims;", dict(Counter(r["status"] for r in R)))
