# Audit 02_timing_A — REPORT

> Rendered verbatim from the agent's structured result (`_structured/audit__02_timing_A.json`); the sub-agent harness blocked agents from writing .md report files. Scripts, outputs and logs are in this folder.


**Verdict:** PASS_WITH_CHANGES

## Summary

The injection-calibrated timing errors are justified, not inflated. I rebuilt the measurement from the raw MAST light curves with my own code; every input was sha256-checked against 00_data/DATA_MANIFEST.json. Data: TESS-SPOC S6 PDCSAP; SPOC S32 2-min PDCSAP; QLP S98 SYS_RM_FLUX with the 2^29/2^30 scattered-light intervals removed (they match the HEASARC S98 note); B masked.

Method: fixed-shape batman template (k=0.0212, b=0.35, rho=2.51, Claret 2017 limb darkening). Local wotan detrending, then a least-squares mid-time scan over +-144 min, with the +-60-min sub-range (the package's range) evaluated on the same chi2 curve. Error models: white, prayer bead (Winn et al. 2008), iid bootstrap, time-averaging beta (Pont et al. 2006; Winn et al. 2008 eq. 2), and 2000 injections per sector. The injections go into the raw flux at the real per-sector depth, are measured identically, and mirror the real coverage/S/N/edge selection. No selected injection hits the search boundary.

Results:
- With the manuscript's +-60-min range: P = 3.1827850 +- 0.0000056 d, T0 = 2460998.7361 +- 0.0029, chi2 = 11.6/18. With an untruncated +-144-min search: P = 3.1827879 +- 0.0000063 d, chi2 = 26.6/18. The difference comes entirely from one S/N-2.5 event, S32 E-573, whose egress sits next to B's transit.
- The injection error distributions are heavy-tailed (robust sigma 8-11 min, SD 13-20 min, excess kurtosis 5-11).
- A Monte Carlo that draws errors from these empirical distributions shows the true scatter of the fitted period equals the quoted sigma_P (5.6e-6 and 6.3e-6 d). The manuscript's chi2 = 12.6 lies at the 30th percentile of the expected chi2 distribution, so it is not suspiciously low.
- Prayer bead, event bootstrap and the white fit rescaled to chi2_red = 1 all give sigma_P = 5.2-7.1e-6 d. The global-fit formal sigma_P (3.2-3.6e-6) is too small.
- White errors underestimate the injection scatter by 1.45x (S6), 6.4x (S32) and 3.5x (S98). My white-error fit reproduces the old chi2 problem (120.4/18).
- A per-sector constant sigma is defensible: per-event white sigma does not predict the scatter in S32 or S98.
- TTV: no evidence (permutation FAP 0.19-0.21). The 95% limit on a sinusoidal semi-amplitude is about 10-13 min for TTV periods under 64 d and 16-23 min for 64-640 d; longer periods are unconstrained.

Changes needed (all MINOR):
1. Quote chi2 as about 12 or 'consistent with the errors', not 12.6.
2. Increase sigma_T0 to about 0.003 d.
3. Require a baseline on both sides of each transit: S98 E-1 has no post-egress baseline.

The package's explanation of the old chi2 = 120 as 'grid quantisation' is incomplete (NOTE).

Harness notes: the harness blocked writing REPORT.md, so its full content is in this result. The other deliverables (7 scripts, CSV/JSON outputs, plots, logs) are in "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/02_timing_A/"; key outputs are out/t5_results.json, out/t6_results.json, out/t8_results.json, out/t9_results.json, out/t4_inj_S*.csv, out/t3_events_*.csv and out/t7_*.png. The ld_lo/ld_hi variants wrote their CSVs but not their meta JSON (a duplicate-keyword error at the final dump); their results were used only through the CSVs.

## Numeric rows

| quantity | manuscript | recomputed | status | method / note |
|---|---|---|---|---|
| P_A (linear ephemeris, 20 transits) | 3.1827856 +- 0.0000058 d | 3.1827850 +- 0.0000056 d (+-60-min search, injection SD); 3.1827879 +- 0.0000063 d (untruncated +-144-min search); MC true SD 5.64e-6 / 6.28e-6 | MATCH | own per-event template fits + WLS; errors from 2000 mirrored injections/sector; MC with empirical errors value within 0.1-0.4 sigma; sigma_P confirmed by MC, event bootstrap (5.7e-6), prayer bead (5.2e-6), rescaled white (5.5e-6) |
| sigma_P justified vs inflated | 0.0000058 d (injection-calibrated) | true scatter of P-hat under empirical injection errors 5.6e-6 (R60) / 6.3e-6 (R144); global-fit formal 3.2e-6, beta-scaled 3.6e-6 (too small) | MATCH | Monte Carlo drawing per-event errors from selected-injection pools (20,000 sets) not inflated |
| T0_A | 2460998.7375 +- 0.0023 BJD_TDB | 2460998.7361 +- 0.0029 (R60); 2460998.7354 +- 0.0031 (R144); white 2460998.7363 +- 0.0009 | MINOR_DIFF | same WLS fits, reference epoch 0 (S98) value within 0.5 sigma; error ~20% too small because the package's S98 injection SD (10.1 min, 150 injections) is low vs mirrored 12.0-13.2 min |
| number of transits in timing fit | 20 | 20 (6 S6 + 6 S32 + 8 S98 vs package 5 + 6 + 9) | MINOR_DIFF | same selection rule (coverage>=0.7, S/N>=2, not at scan edge) my set includes S6 E-795 and excludes S98 E-1 (no post-egress baseline, fit unstable) |
| chi2 of linear ephemeris | 12.6 for 18 dof | 11.6/18 (R60); 26.6/18 (R144; S32 E-573 at -82 min, z=-4.2); without E-573 8.7/17 | MINOR_DIFF | WLS with per-sector injection SD 12.6 sits at the 30th percentile of the MC chi2 distribution (not low); 26.6 at the 83rd (p_emp 0.17) |
| previous draft chi2 with white errors | 120/21 (draft) | 120.4/18 with formal white errors; beta-scaled 84.3/18 | MATCH | white Delta-chi2=1 errors cause: white sigma too small by 1.45x/6.4x/3.5x (S6/S32/S98), not grid quantisation |
| per-sector injection timing sigma | 13.8 / 24.6 / 10.1 min (package a1; not in manuscript text) | SD 14.4 / 19.8 / 13.2 (R60, coverage cut only); 12.9 / 14.7 / 12.0 (R60 selected); 14.7 / 19.8 / 12.9 (R144 selected); robust 10.7 / 8.7 / 7.8 | MINOR_DIFF | 2000 mirrored injections/sector package used 150 injections without S/N/edge mirroring |
| rms O-C (package claim) | 9 min (package) | 9.9 min (R60); 21.1 min (R144, dominated by E-573) | MATCH | residuals of the WLS fit  |
| joint-fit formal sigma_P (package, not in manuscript) | 0.0000027 (not used) | global fit formal 0.0000032, beta-scaled 0.0000036; both too small vs injection-calibrated 0.0000056-63 | MATCH | joint least squares over 20 windows correctly not used |
| joint fit transit count | 23 transits | 23 events with coverage>=0.5 (6+8+9) in my set | MATCH | event set definition composition differs slightly from the package (S6 E-795 in, S98 E-1 in) |
| TTV | no TTV statement (draft sentence dropped) | no evidence: permutation FAP 0.19-0.21, empirical-noise FAP 0.25-0.66; 95% semi-amplitude limit ~10-13 min (P_TTV<64 d), ~16-23 min (64-640 d) | MATCH | sinusoid periodogram (4040 frequencies), permutation/empirical FAP, Neyman upper limits  |

## Findings

### F1 [MINOR] chi2 = 12.6/18 depends on the +-60 min search range through one low-S/N event

With an untruncated +-144-min search, S32 E-573 (BTJD 2175.001; S/N 2.5; egress next to B's transit) moves to -82 min (z = -4.2). Then chi2 = 26.6/18 (Gaussian p = 0.086, empirical p = 0.17) and P = 3.1827879 +- 0.0000063 (+0.4 sigma). The conclusions do not change, but three significant figures overstate the precision of chi2.

*Evidence:* out/t5_results.json ephemeris_20 and leave_one_out; added_events_and_variants.search_R60_sel60 (chi2 11.6); logs/t5.log; out/t7_oc.png

*Location:* Sec. 2, first sentence

*Proposed fix:* Write 'consistent with these errors (chi2 ~ 12 for 18 dof)' or state the +-60-min search range.

### F2 [MINOR] sigma_T0 = 0.0023 d is about 20% too small

sigma_T0 comes mostly from S98. The package's S98 injection SD (10.1 min from 150 injections) is below the 2.5th percentile of 150-injection subsets of a mirrored 2000-injection run (95% range 10.7-15.6 min; full SD 12.0-13.2).

*Evidence:* out/t8_results.json SD_uncertainty_R60_cov0p7 and R60_MC (T0 3998.73611 +- 0.00287); R60 fit with the package SDs gives eT0 = 0.00244

*Location:* Sec. 2

*Proposed fix:* T0 = 2460998.736 +- 0.003 BJD_TDB

### F3 [MINOR] One of the manuscript's 20 events (S98 E-1) has no post-egress baseline

At BTJD 3995.554 the data stop about 5 min after egress (2^30 flags, then the downlink gap), so the time is set by the baseline slope. The package measured +19.5 min with depth 1093 ppm (twice the true depth); my fit prefers an anti-transit. The coverage criterion does not catch one-sided baselines. Swapping this event for S6 E-795 leaves P and chi2 essentially unchanged; adding it to my set raises chi2 to 56.9/21.

*Evidence:* logs/compare_with_package.log (S98 E-1: package +19.5 min vs mine -75.7 min, S/N -4.6); out/t3_events_base.png; t5 add_all_23_in_set

*Location:* Sec. 2 ('from 20 transits')

*Proposed fix:* Require a baseline on both sides of each transit (e.g. >= 0.5 T14 of data before and after); the count stays 20.

### F4 [NOTE] The old chi2 = 120 comes from non-Gaussian timing errors, not grid quantisation

At per-event S/N 3-6, white errors underestimate the injection scatter (z SD 1.45 / 6.4 / 3.5 for S6/S32/S98). Correlated noise explains little (beta = 1.00 / 1.26 / 1.10); the excess is fit-failure tails. My white-error fit reproduces the problem: 120.4/18.

*Evidence:* out/t5_results.json calibration.z_stats; betas; ephemeris_20.white

*Location:* package CHANGES_FROM_CURRENT_DRAFT.md (not manuscript)

*Proposed fix:* Correct the explanation in the change log.

### F5 [NOTE] Package injections: only 150 per sector, selection not mirrored, heavy tails

Mirrored injections (>= 1693 selected per sector) give different per-sector SDs but the same sigma_P (5.6 vs 5.8e-6). The distributions are heavy-tailed: robust sigma 8-11 min, SD 13-20 min, excess kurtosis 5-11, and 0.4-2.1% of events beyond 60 min. The SD, not the robust sigma, is the right Gaussian-equivalent error for WLS: the robust core alone gives chi2 = 116/18.

*Evidence:* out/t4_inj_S{6,32,98}.csv; out/t5_results.json injections; out/t7_injections.png

*Location:* Sec. 2 method

*Proposed fix:* none required

### F6 [NOTE] A per-sector constant sigma is adequate

Per-event white sigma is a poor predictor of the recovered-time scatter in S32 and S98: Gaussian log-likelihood favours the constant by 296 (S32) and 172 (S98); only S6 favours a scaled per-event sigma (+71). The Student-t likelihood gives the same ranking.

*Evidence:* out/t5_results.json calibration.gauss_ML, student3_ML, by_e_white_quartile

*Location:* Sec. 2 method

*Proposed fix:* none

### F7 [NOTE] Global-fit formal sigma_P is too small

The global joint fit gives sigma_P = 3.2e-6 (formal) or 3.6e-6 (beta-scaled), and the audit's joint MCMC gave 2.7e-6. All are about 1.7x smaller than the injection-calibrated value because they miss the error tails. The manuscript correctly does not use them.

*Evidence:* out/t6_results.json global_fit_20

*Location:* Sec. 2

*Proposed fix:* none

### F8 [NOTE] S6 momentum-dump proximity does not degrade timing

Real S6 transits lie 0.4-3.5 h from dumps, and the generic injections cannot sample that. Dump-mirrored injections (A divided out) give SD 12.5 min (robust 11.9) vs 13.2 (9.9) in a random control, with no bias.

*Evidence:* out/t9_results.json

*Location:* Sec. 2

*Proposed fix:* none

### F9 [NOTE] TTV: no evidence; the manuscript makes no TTV implication

Linear-ephemeris chi2 is acceptable (p = 0.74-0.87 for R60; 0.086-0.17 for R144). Best sinusoid: permutation FAP 0.19-0.21, empirical-noise FAP 0.25-0.66. 95% semi-amplitude upper limit: median ~10-13 min for P_TTV < 64 d and ~16-23 min for 64-640 d; longer periods are unconstrained. The abstract and Sec. 2 contain no TTV claim.

*Evidence:* out/t6_results.json ttv; out/t6_ttv_upper_limits_*.csv; out/t7_ttv.png

*Location:* abstract, Sec. 2

*Proposed fix:* Optional: add 'no evidence for timing variations (95% limit ~15 min for TTV periods below ~2 yr)'.

## Recommended wording

Transit times from a fixed-shape batman template, with per-sector errors (12-15 min) calibrated by injecting the template at random times, give P = 3.182785 +- 0.000006 d and T0 = 2460998.736 +- 0.003 BJD_TDB from 20 transits; the residuals are consistent with these errors (chi2 ~ 12 for 18 degrees of freedom), with no evidence for timing variations (95% limit ~15 min for TTV periods below ~2 yr).

## Files

- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/02_timing_A/t2_common.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/02_timing_A/t3_real_events.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/02_timing_A/t3b_plot_events.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/02_timing_A/t4_injections.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/02_timing_A/t5_calibration_ephemeris.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/02_timing_A/t6_global_ttv.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/02_timing_A/t7_plots.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/02_timing_A/t8_extra_checks.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/02_timing_A/t9_s6_dump_injections.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/02_timing_A/out/t3_events_base.csv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/02_timing_A/out/t4_inj_S6.csv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/02_timing_A/out/t4_inj_S32.csv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/02_timing_A/out/t4_inj_S98.csv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/02_timing_A/out/t5_results.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/02_timing_A/out/t6_results.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/02_timing_A/out/t8_results.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/02_timing_A/out/t9_results.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/02_timing_A/out/t7_injections.png
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/02_timing_A/out/t7_oc.png
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/02_timing_A/out/t7_ttv.png
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/02_timing_A/logs/compare_with_package.log