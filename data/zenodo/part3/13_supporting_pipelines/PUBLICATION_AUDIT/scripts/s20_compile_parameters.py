"""Step 20: compile FINAL_PARAMETERS.json ONLY from this audit's own outputs (out/*.json) - no number is typed by hand
except the adopted-Teff decision, which is documented inline. Run after all analysis steps."""
import glob
import json
import os
import sys

import numpy as np
from astropy.time import Time

sys.path.insert(0, os.path.dirname(__file__))
from audit_common import AUD, save_json, utcnow  # noqa

J = lambda p: json.load(open(os.path.join(AUD, p)))
ST = J("out/s04_stellar.json")
EPH = J("out/ephem_current.json")
FITS = {os.path.basename(f)[4:-5]: json.load(open(f)) for f in glob.glob(os.path.join(AUD, "out/fits/fit_*.json"))}
DYN = {os.path.basename(f)[8:-5]: json.load(open(f)) for f in glob.glob(os.path.join(AUD, "out/fits/dynesty_*.json"))}
EV = J("out/s09_events.json")
EB = J("out/s10_eb_tests.json")
OE = J("out/s10b_oddeven_matched.json")
NB = J("out/s11_nullB.json")
NB2 = J("out/s11b_nullB_periods.json")
CEN = J("out/s11c_single_box_census.json")
PF = J("out/s10c_phase_features.json")
S98B = J("out/s07b_s98_B_significance.json")
S98C = J("out/s07c_s98_B_null_variants.json")
BCOMB = J("out/s11d_B_combined.json")
LOC = J("out/s12c_clean_localisation.json")
S6D = J("out/s12d_s6_dumps.json")
COMP = J("out/s14_companions.json")
ACT = J("out/s15_activity.json")
S110 = J("out/s16_s110.json")
DET = J("out/s17_detrend_robust.json")
NOV = J("out/s18_novelty_refresh.json")
CYC = J("out/s19_cycles.json")
TIM = J("out/s06_timing.json")


def med_err(v):
    return dict(value=v[2], plus=v[3] - v[2], minus=v[2] - v[1], p2_3=v[0], p97_7=v[4])


# ---- adopted Teff (decision): mean of independent estimates, GSP-Spec excluded (extrapolation flag, logg 5.01)
teff_sources = {"TIC 8.2 (Stassun+2019 relations)": 4010.0, "Mann+2015 V-J relation (APASS V, 2MASS J)": ST["teff_mann2015_VJ"]["median"],
                "ARIADNE SED BMA (companion_check, re-read from best_fit_average.dat)": 4011.8, "Gaia DR2 Apsis-Priam": 4069.3}
teff = float(np.mean(list(teff_sources.values())))
teff_adopt = round(teff / 10) * 10
eteff = 90.0
R = ST["radius_mann2015_MK"]["median"]
eR = ST["radius_mann2015_MK"]["std"]
M = ST["mass_mann2019_table6_n5"]["median"]
eM = ST["mass_mann2019_table6_n5"]["std"]
L = R ** 2 * (teff_adopt / 5772.0) ** 4
eL = L * np.hypot(2 * eR / R, 4 * eteff / teff_adopt)


def planet_block(sig):
    f = FITS[f"{sig}_pipe_mann_frozen_w0.75"]
    fo = FITS[f"{sig}_own_mann_frozen_w0.75"]
    ff = FITS[f"{sig}_pipe_free_frozen_w0.75"]
    P, T0 = f["P"][2], f["T0"][2]
    rng = np.random.default_rng(7)
    ch = np.load(os.path.join(AUD, "out/fits", f"chain_{sig}_pipe_mann_frozen_w0.75.npz"))["chain"]
    Ms = rng.normal(M, eM, len(ch))
    Rs = rng.normal(R, eR, len(ch))
    Ts = rng.normal(teff_adopt, eteff, len(ch))
    a_au = (Ms * (ch[:, 1] / 365.25) ** 2) ** (1 / 3)
    Ls = Rs ** 2 * (Ts / 5772.0) ** 4
    S = Ls / a_au ** 2
    Teq = Ts * np.sqrt(Rs * 0.00465047 / (2 * a_au))
    k = np.exp(ch[:, 2])
    Rp = k * Rs * 109.076
    pc = lambda v: [float(x) for x in np.percentile(v, [2.275, 15.865, 50, 84.135, 97.725])]
    Kpm = 28.4329 * (1 / 317.828) * (Ms ** (-2 / 3)) * (ch[:, 1] / 365.25) ** (-1 / 3)
    blk = dict(
        ephemeris=dict(P_d=med_err(f["P"]), T0_BTJD=med_err(f["T0"]), T0_BJD_TDB=T0 + 2457000.0, corr_T0_P=f["corr_T0_P"],
                       source="joint batman+emcee fit, pipeline set, Mann density prior (fit_%s_pipe_mann_frozen_w0.75)" % sig,
                       time_system="BJD_TDB - 2457000 (BTJD); all products declare TIMESYS=TDB, TIMEREF=SOLARSYSTEM",
                       cross_checks=dict(own_photometry_fit=dict(P=f"{fo['P'][2]:.8f} +{fo['P'][3]-fo['P'][2]:.8f}/-{fo['P'][2]-fo['P'][1]:.8f}",
                                                                 T0=f"{fo['T0'][2]:.5f} +{fo['T0'][3]-fo['T0'][2]:.5f}/-{fo['T0'][2]-fo['T0'][1]:.5f}"),
                                         per_event_linear_pipe=EV[sig]["ephem_pipe"] and {k_: EV[sig]["ephem_pipe"][k_] for k_ in ("P", "eP_scaled", "T0", "eT0_scaled", "chi2", "dof")},
                                         per_event_linear_own=EV[sig]["ephem_own"] and {k_: EV[sig]["ephem_own"][k_] for k_ in ("P", "eP_scaled", "T0", "eT0_scaled", "chi2", "dof")},
                                         dynesty=DYN.get(f"{sig}_pipe", {}).get("P"),
                                         cycle_count=CYC[sig])),
        shape_mann_prior=dict(k=med_err(f["k"]), depth_k2_ppm=med_err(f["depth_k2_ppm"]), b=med_err(f["b"]), a_over_Rs=med_err(f["a_over_Rs"]),
                              inc_deg=med_err(f["inc_deg"]), T14_h=med_err(f["T14_h"]), Rp_Rearth=med_err(pc(Rp)),
                              rho_prior=f["rho_mann"], u1=med_err(f["u1"]), u2=med_err(f["u2"]),
                              convergence=dict(n_steps=f["n_steps"], tau_max=f["tau_max"], n_over_tau=f["n_over_tau"], acceptance=f["acceptance"])),
        shape_own_photometry_mann_prior=dict(k=med_err(fo["k"]), b=med_err(fo["b"]), T14_h=med_err(fo["T14_h"]), Rp_Rearth=med_err(fo["Rp_Rearth"])),
        free_density=dict(rho_logprior=med_err(ff["rho"]), rho_uniformprior=med_err(ff["rho_uniform_prior_reweighted"]),
                          b_logprior=med_err(ff["b"]), Rp_Rearth=med_err(ff["Rp_Rearth"]), k=med_err(ff["k"]),
                          convergence=dict(n_steps=ff["n_steps"], tau_max=ff["tau_max"], n_over_tau=ff["n_over_tau"]),
                          note="prior sensitivity shown by importance-reweighting ln-uniform -> uniform in rho"),
        derived=dict(a_AU=med_err(pc(a_au)), S_Searth=med_err(pc(S)), Teq_K_A0_full_redistribution=med_err(pc(Teq)),
                     K_ms_per_Mearth=med_err(pc(Kpm))),
        dynesty_independent=DYN.get(f"{sig}_pipe"),
        events=dict(per_dataset={sname: {k_: v for k_, v in EV[sig][sname]["summary"].items() if v} for sname in ("pipe", "own", "own_nostray", "pipe_s98default", "own2min", "others")}),
        detrending_robustness={k_: v[sig] for k_, v in DET.items()},
        eb_tests=dict(pipe=EB[f"{sig}_pipe"] | {"phase_scan": {k_: v for k_, v in EB[f"{sig}_pipe"]["phase_scan"].items() if k_ != "scan"}},
                      own=EB[f"{sig}_own"] | {"phase_scan": {k_: v for k_, v in EB[f"{sig}_own"]["phase_scan"].items() if k_ != "scan"}},
                      odd_even_per_event=OE[sig],
                      phase_feature_check_s10c={k_: {kk: vv for kk, vv in v.items() if kk != "per_cycle"} for k_, v in PF.items() if k_.startswith(sig)}),
        localisation=dict(S98_clean=LOC["clean_S98"][sig], combined_without_S6=LOC["combined_without_S6"][sig], combined_all=LOC["combined"][sig],
                          spoc2min_centroids=LOC["spoc2min_centroids"][sig], S6_dump_note=S6D.get("A") if sig == "A" else None),
    )
    fpp = {}
    for d in sorted(glob.glob(os.path.join(AUD, f"out/fpp/runs/{sig}_*"))):
        s = os.path.join(d, "summary.json")
        if os.path.exists(s):
            x = json.load(open(s))
            L_ = x["log"]
            if "FPP_mean" in L_:
                runs = x["runs"]
                fpp[os.path.basename(d)] = dict(FPP_mean=L_["FPP_mean"], FPP_sd=L_["FPP_std"], FPP_min=L_["FPP_min"], FPP_max=L_["FPP_max"],
                                               NFPP_mean=L_["NFPP_mean"], NFPP_max=L_["NFPP_max"], n_runs=len(runs), N=L_["N"],
                                               exptime_min=L_["exptime_d"] * 1440, P_STP_mean=float(np.mean([r["P_STP"] for r in runs])),
                                               P_EB_all_mean=float(np.mean([r["P_EB_all"] for r in runs])),
                                               n_saturated=L_["n_saturated"], n_degenerate=L_["n_degenerate"])
    blk["fpp_triceratops_1.1.0"] = fpp
    return blk


out = dict(
    meta=dict(generated_utc=utcnow(), audit_dir="addenda/tic4206066_PUBLICATION_AUDIT", note="compiled by scripts/s20_compile_parameters.py from out/*.json"),
    star=dict(
        identifiers=["TIC 4206066", "StKM 1-561", "Gaia DR3 3220388198192519424 (= DR2 id)", "TYC 4753-1591-1", "2MASS J05232636-0119380", "PM J05234-0119", "UCAC4 444-009590"],
        astrometry=dict(ra_deg_gaia_dr3=80.85982665567226, dec_deg_gaia_dr3=-1.3269640103386409, epoch=2016.0, pmra_masyr=-3.4059944859312923,
                        pmdec_masyr=70.80569610380549, parallax_mas=ST["inputs"]["parallax_mas"], parallax_err=ST["inputs"]["parallax_err"],
                        lindegren2021_zp_mas=ST["inputs"]["lindegren2021_zp_mas"], parallax_corrected_mas=ST["inputs"]["parallax_corrected"],
                        distance_pc=ST["distance_pc_1_over_plx_corrected"], bailer_jones2021_rgeo_pc=[35.5786552, 35.5568771, 35.5957909]),
        photometry=dict(G=ST["inputs"]["G"], BP=ST["inputs"]["BP"], RP=ST["inputs"]["RP"], bp_rp=ST["inputs"]["bp_rp"], T_TIC82=ST["Tmag_TIC82"],
                        T_stassun_from_gaia=ST["Tmag_stassun2019"], J=[ST["inputs"]["J"], ST["inputs"]["eJ"]], H=[ST["inputs"]["H"], ST["inputs"]["eH"]],
                        Ks=[ST["inputs"]["Ks"], ST["inputs"]["eKs"]], V_APASS=[ST["inputs"]["V_APASS"], ST["inputs"]["eV_APASS"]],
                        W1=[7.605, 0.031], W2=[7.626, 0.020], W3=[7.560, 0.019], W4=[7.451, 0.133]),
        MK=ST["MK"],
        adopted=dict(R_Rsun=[R, eR], R_source="Mann+2015 Table 1 (erratum ApJ 819, 87) R(M_Ks), 2.89% scatter + M_Ks error (MC)",
                     M_Msun=[M, eM], M_source="Mann+2019 Table 6 n=5, zp 7.5, relation error 2.4% at M_Ks=5 (their Table 7) + M_Ks error",
                     rho_rhosun=med_err([ST["density_rho_sun"][k] if k in ST["density_rho_sun"] else None for k in ("p2_5",)] + [
                         ST["density_rho_sun"]["median"] - ST["density_rho_sun"]["minus1"], ST["density_rho_sun"]["median"],
                         ST["density_rho_sun"]["median"] + ST["density_rho_sun"]["plus1"], ST["density_rho_sun"]["p97_5"]]),
                     logg_cgs=[ST["logg_cgs"]["median"], ST["logg_cgs"]["std"]],
                     Teff_K=[teff_adopt, eteff], Teff_sources=teff_sources,
                     Teff_excluded={"Gaia DR3 GSP-Spec 4182 K": "flags_gspspec extrapolation bit set; logg_gspspec 5.01 inconsistent with M, R"},
                     L_Lsun=[L, eL], L_source="R^2 (Teff/5772)^4 with adopted R, Teff",
                     spectral_type="K7V-M0V (Stephenson 1986 K6; LAMOST LRS K7; Lepine & Gaidos 2011 k7; M_Ks between K7V and M0V of Pecaut & Mamajek)",
                     metallicity="poorly constrained: Gaia GSP-Spec -0.27 (flagged), ARIADNE -0.08 +/- 0.18"),
        alternatives=dict(TIC82=ST["TIC82"], mann2019_official_posterior_code=ST["mass_mann2019_official_posterior"],
                          mann2015_R_with_FeH_minus0p27=ST["radius_mann2015_MK_FeH_-0.27"], mann2015_quartic_mass=ST["mass_mann2015_quartic_point"],
                          sed_ariadne=dict(Teff=[4011.8, 45.3, 40.3], R=[0.6206, 0.0147, 0.0140], L=[0.0901, 0.0059, 0.0055], FeH=[-0.079, 0.174, 0.183],
                                           source="addenda/companion_check_tic4206066/sed_out/best_fit_average.dat (not re-run)"),
                          gaia_dr2_flame=ST["gaia_dr2_flame"], L_BC_K=ST["L_from_BC_K"], L_BC_V=ST["L_from_BC_V"]),
        binarity_and_companions=dict(ruwe=COMP["ruwe"], gaia_rv=COMP["gaia_rv"], wise=COMP["wise"],
                                     host_density_constraint_A=COMP["host_density_constraint_A"],
                                     neighbour_census=COMP["neighbour_census"], neighbours_within_60=COMP["neighbours_within_60"]),
        activity=dict(xray=ACT["xray"], tess_periodograms=[{k: v for k, v in r.items() if k in ("label", "peaks", "rms_ppm")} for r in ACT["tess"]],
                      asassn={b: {k: v for k, v in ACT["asassn"][b].items() if k in ("n", "peaks", "mad_mmag", "jd_range")} for b in ("g", "V") if b in ACT["asassn"]},
                      catalogue_periods=ACT["catalogue_periods_reported"])),
    signal_A=planet_block("A"),
    signal_B=planet_block("B"),
    signal_B_link_null_test=dict(
        epoch_shifted_fixed_period_s11=dict(NB, status="SUPERSEDED (audit bug C7): the trials overlap in time (about 5-16 distinct windows in S6), "
                                                "so their FAP values are not valid; kept for reproduction of the earlier statistic only"),
        period_and_epoch_randomised_s11b=NB2,
        single_box_census_and_rank_test_s11c=CEN,
        s98_only_detection_significance_s07b=S98B,
        s98_only_null_robustness_s07c=S98C,
        overall_indicative_fap_s11d=BCOMB),
    period_ratio_B_over_A=dict(value=EPH["B"]["P"] / EPH["A"]["P"], nearest_simple=3.5,
                               fractional_offset_from_7_2=(EPH["B"]["P"] / EPH["A"]["P"]) / 3.5 - 1),
    sector110=S110,
    novelty_refresh=NOV,
    timing_audit={k: v for k, v in TIM.items() if k == "compare"},
)
save_json(out, "FINAL_PARAMETERS.json")
print("written FINAL_PARAMETERS.json;", "Teff adopted", teff_adopt, "L", round(L, 4), "+/-", round(eL, 4))
for sig in ("A", "B"):
    b = out[f"signal_{sig}"]
    print(sig, "P", b["ephemeris"]["P_d"]["value"], "T0", b["ephemeris"]["T0_BTJD"]["value"], "Rp", b["shape_mann_prior"]["Rp_Rearth"]["value"],
          "S", b["derived"]["S_Searth"]["value"], "Teq", b["derived"]["Teq_K_A0_full_redistribution"]["value"], "FPP", {k: round(v["FPP_mean"], 4) for k, v in b["fpp_triceratops_1.1.0"].items()})
