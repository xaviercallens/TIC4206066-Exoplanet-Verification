import json
W = "/Users/paulsmith/space <3/tess-jev/addenda/tic4206066_FINAL_RECONCILIATION/_work/T3_pixels_localisation/"
UA = "addenda/opus55_final_audit_tic4206066/"; PP = "addenda/rnaas_tic4206066_pubprep/"; AU = "addenda/tic4206066_PUBLICATION_AUDIT/"
WK = "addenda/tic4206066_FINAL_RECONCILIATION/_work/T3_pixels_localisation/"
F = []
def f(id, group, q, v, u, units, status, ev, code, rerun, sup, mf, notes):
    F.append(dict(id=id, group=group, quantity=q, value=v, uncertainty=u, units=units, status=status, evidence=ev, code=code, rerun=rerun, supersedes=sup, manuscript_form=mf, notes=notes))
f("A.offset", "PIXEL", "3.18 d source position relative to the target (dE, dN), S32 SPOC 2-min POS_CORR-decorrelated + S98 TESScut clean (9 events), inverse-variance combination",
  "(+5.32, +0.20)", "(2.66, 2.37) per axis; per-set error = max(injection SD, jackknife SE, robust bootstrap SD)", "arcsec", "FINAL",
  UA+"referee/resolve_R3/logs/r2_combine.log line 'A_eventscatter_auditor_pos' ; "+UA+"proposed_revision/figure_src/figure_values_v2.json key panel_b.combined.A.exact ; "+WK+"t3_combine.json key variants.U_headline_S32PC+S98U",
  UA+"03_pixels/pixloc.py, run_real.py, poscorr_decorr.py, summarize.py; "+UA+"referee/resolve_R3/r2_combine.py", WK+"t3_combine.json ; "+WK+"t3_s98_driver.json",
  "pubprep (+7.90, -8.20) +- (5.27, 6.84), printed (+8 +- 5, -8 +- 7) [a5_summary.json localisation.A.combined_S32tpf_S98clean]; audit (+2.51, -7.53) +- (5.68, 5.69) [FINAL_PARAMETERS.json signal_A.localisation.combined_without_S6]; U 03 injection-only (+5.47 +- 2.23, +0.17 +- 2.35) [summary.json combined A:S32tpfPC+S98clean29:bias0:rob0]",
  "$(+5\\pm3, 0\\pm3)''$",
  "Driver of the old -8\": raw S32 carries no weight (pubprep S32 sigma 17-22\", U raw inj cov 8.7-11\") so the S98 set alone set the Dec. After POS_CORR decorrelation S32 gives (+5.21, +0.39) +- (3.08, 2.78) (S/N 3.5 -> 7.7; independently reproduced by the U verifier (+5.13, +0.20)). With decorrelated S32, the combined position is (+4.5..+6.3, -2.2..+0.2)\" for every S98 recipe tried (U, audit 9-event, pubprep 10-event, U quadratic baseline, U 14-event 2048-only) and (+4.99, -2.01) with the S32 TESS-SPOC FFI decorrelated instead (own rerun t3_combine.json). Rounded form covers all variants.")
f("A.S32_per_set", "PIXEL", "S32 SPOC 2-min 3.18 d offset: raw vs POS_CORR-decorrelated", "decorrelated (+5.21, +0.39); raw (+9.12, -8.31) (U) / (+13.46, -2.25) bias-corrected (pubprep)", "decorrelated (3.08, 2.78); raw pubprep (17.19, 21.55)", "arcsec", "FINAL",
  UA+"03_pixels/summary.json keys sets.S32_tpf_poscorr_A.real_off, sets.S32_tpf_A.real_off ; "+UA+"03_pixels/verify/v1_loc.json (usepc true: off [5.13, 0.20]) ; "+PP+"analysis/out/a5_summary.json key localisation.A.per_set.S32_tpf",
  UA+"03_pixels/poscorr_decorr.py; "+UA+"03_pixels/verify/v1_loc.py", "", "audit z (0.746, -0.101) and pubprep z (+0.774, -0.191) of the raw S32 set (C077) - valid for raw data, superseded by the decorrelated set", "not in Note (only 'Sector 32, pointing-decorrelated')",
  "Raw-S32 Dec -8\" is driven by event 2175.001 (jackknife drop -> -2.1\"; U/referee/resolve_R3/logs/r5_s98_gap.log). The U claim of a 2.7-sigma transit-correlated POS_CORR2 shift is overstated: 2.0 sigma by U's own pointing.json, ~1.1 sigma unweighted (verify/v2_poscorr.json); the audit's empirical-null POS_CORR2 z is 0.30 (FINAL_PARAMETERS.json signal_A.localisation.spoc2min_centroids). Describe it as pointing jitter inflating the noise, not as a significant transit-correlated pointing shift.")
f("A.S98_per_set", "PIXEL", "S98 TESScut clean 3.18 d offset (per-set; enters the combination and Fig. 1b)", "U (+5.65, -0.32) [9 events]; audit recipe (+1.17, -8.41) [same 9 events]; pubprep (+7.03, -9.29) [10 events]; audit recipe with U's clean29 mask (+9.31, -10.41) raw; U quadratic baseline 14 events (+7.31, -7.22)",
  "U adopted (5.27, 4.52); audit inj (5.98, 5.93); pubprep (5.60, 7.29); audit-recipe jackknife swings Dec from +2.8 to -17 (one drop gives -56)", "arcsec", "UNRESOLVED",
  UA+"03_pixels/diagnose_s98.json keys S98_clean29_A, S98_clean29_priorlike_A ; "+AU+"out/s12c_clean_localisation.json key clean_S98.A ; "+PP+"analysis/out/a5_summary.json key localisation.A.per_set.S98_ffi ; "+WK+"t3_s98_driver.json keys AUD_frozen, AUD_clean29, U_clean29, U_frozenmask ; "+UA+"03_pixels/summary.json key sets.S98_clean_deg2_A",
  AU+"scripts/s12_pixels.py, s12c_clean_localisation.py (imported read-only); "+UA+"03_pixels/pixloc.py", WK+"t3_s98_driver.py -> t3_s98_driver.json/.log",
  "", "not in Note (Fig. 1b plots U value (+5.6 +- 5.3, -0.3 +- 4.5))",
  "Identical pixel data (t3_cube_compare.json). The audit recipe reproduces its stored real_px (0.05, -0.45) exactly. With the same 9 events the Dec is recipe-dependent at the 8-10\" level: U's continuous +-4 h box-regression recipe gives -0.3 and is jackknife-stable (N range -1.5..+2.0); the audit's windowed-baseline recipe gives -8.4 and is jackknife-unstable. So U's attribution of the -8\" to the S98 gap-edge event explains pubprep's 10-event value (U priorlike: -8.1 -> -1.6 when 3995.555 is dropped) but NOT the audit's -8.4, which already excludes that event. Does not affect A.offset (see there).")
f("A.S98_epochs", "PIXEL", "Number of S98 3.18 d events usable for difference imaging", "9 (BTJD 3998.737 ... 4043.296); event 3995.555 excluded", "", "events", "FINAL",
  WK+"t3_gap3995.json (TESScut and QLP gap 3995.608-3995.821, QUALITY 32 dump at 3995.821) ; "+UA+"referee/resolve_R3/logs/r5_s98_gap.log ; "+PP+"analysis/a2_prf_localisation.py L178-181 (one-sided events accepted with constant baseline) and L145 (fake/injection sets require both sides) ; "+AU+"out/s12c_clean_localisation.json clean_S98.A.epochs",
  PP+"analysis/a2_prf_localisation.py; "+AU+"scripts/s12_pixels.py epochs_for", WK+"t3_gap3995.py -> t3_gap3995.json", "pubprep 10 epochs [a5 localisation.A.per_set.S98_ffi.n_epochs]",
  "not in Note",
  "Egress of 3995.555 is at 3995.594; a 5.1-h data gap follows (3995.608-3995.821, ending with a momentum dump), so 0 cadences fall in the 1.40-4.68 h post-transit window. Pubprep used a one-sided constant baseline for it, but its injections never contain one-sided events, so that event is outside its own error calibration. Audit and U both use 9 events.")
f("A.offset_significance", "PIXEL", "Offset of the 3.18 d source from the target", "z_E = 2.00; chi2 = 4.01 (2 dof), p = 0.135 (1.49 sigma 2-D equivalent)", "across error models p = 0.049-0.166 (R3); across S98 recipes with decorrelated S32 p = 0.020-0.224, z_E 1.66-2.58; rounded (+5+-3, 0+-3): z_E 1.67, p = 0.25", "sigma", "FINAL",
  UA+"proposed_revision/figure_src/figure_values_v2.json key panel_b.combined.A.exact (z_E 1.999663721290341, p 0.1349557197582979) ; "+UA+"referee/resolve_R3/logs/r2_combine.log ; "+WK+"t3_combine.json variants.*.z_E/p",
  UA+"referee/resolve_R3/r2_combine.py", WK+"t3_combine.json", "pubprep 'consistent at 2 sigma' (chi2 5.28); audit chi2 2.09",
  "consistent with the target within $\\approx2\\sigma$",
  "A consistent ~5\" eastward offset appears in S6 (decorrelated/dump-excised), S32 (decorrelated) and S98; S98 injections show a +2.0\" E method bias (summary.json sets.S98_clean29_A.inj_mean). Not significant.")
f("A.capable_stars", "PIXEL", "Gaia DR3 stars able to produce the 3.18 d depth (required eclipse depth <= 1 in the real S6/S32/S98 apertures, delta = 500 ppm)",
  "5 stars: TIC 4206077 = Gaia DR3 3220388296975320448 (T 17.02) at 43.1\"; TIC 4206064 at 50.5\"; TIC 4206073 at 53.9\"; TIC 4206059 at 56.3\"; TIC 4206057 (T 13.06) at 84.1\" (PRF wings)", "separations at sector epochs; pubprep 43.5/50.3/54.0/56.0/83.8\"", "arcsec", "FINAL",
  UA+"03_pixels/neighbours_vs_fitted_position.csv rows with dreqA_min<=1 ; "+UA+"03_pixels/neighbour_census_summary.json capable ; "+PP+"MANUSCRIPT_NUMBER_AUDIT.csv L50",
  UA+"03_pixels/neighbours.py", "", "", "$\\ge43''$",
  "A per-product 'any aperture' census also lists TIC 712986149 (51.0\"), 4206069 (91.3\"), 4206056 (97.7\") - all farther than 43\". Only the T 13.06 star at 84\" qualifies in every aperture or with a 0.5 threshold.")
f("A.capable_sigma", "PIXEL", "Minimum Mahalanobis distance of depth-capable stars from the fitted 3.18 d position",
  "14.3 (rounded +-3\" errors; TIC 4206077); 17.6-18.0 (exact event-scatter covariance)", "14.2 (sector-epoch positions, figure_values_v2) to 14.8 (TIC J2000 positions); >= 15.5 for every decorrelated-S32 variant", "sigma (Mahalanobis)", "FINAL",
  UA+"referee/resolve_R3/logs/r2_combine.log 'capable-star distance corrected_rounded(+5+-3,0+-3)' ; "+UA+"proposed_revision/figure_src/figure_values_v2.json panel_b.neighbours[0].sigma_distance_drawn_3arcsec ; "+WK+"t3_combine.json printed / variants",
  UA+"referee/resolve_R3/r2_combine.py", WK+"t3_combine.json",
  ">=7 sigma (pubprep; 7.18 from (+8+-5,-8+-7), 6.99 with TIC J2000 positions); >=8 sigma (audit; 8.05 from (+2.5,-7.5)+-5.7); '>10 sigma' (R3 recommendation); all from superseded positions except R3",
  "($>10\\sigma$) [formal minimum 14$\\sigma$]",
  "Gaussian extrapolation to 14 sigma is not validated: per-set estimators have heavy tails (audit-recipe jackknife -56\"; U S98 clean29 deg2 variant (+32.9, -16.5)). Independent support: sources injected ~41\" off target are recovered at 0.99-1.00 of their offset with target dchi2 >= 150 (resolve__R3.json), and PRF-shaped injections on the six nearest capable stars reproduce the real images in 0 of 150 (audit PUBLICATION_TRUTH.md L157-158). '>10 sigma' is the defensible printed form.")
f("A.nearest_gaia_neighbour", "PIXEL", "Nearest Gaia DR3 source of any brightness", "Gaia DR3 3220388296975929984 (TIC 712986249), 11.06\" at J2016.0 (10.6-10.9\" at the TESS epochs), G 18.75, T 18.13 (TIC 8.2) / 18.15 (pubprep)",
  "needs a 120-123% eclipse at delta = 500 ppm (max ~409-410 ppm); 4.1 sigma (+-3\") / 5.0 sigma (exact) from the fitted position", "arcsec", "FINAL",
  UA+"referee/resolve_R3/logs/r4_panel_d_spec.log nearest_other_source ; "+UA+"03_pixels/neighbours_vs_fitted_position.csv row 1 ; "+PP+"MANUSCRIPT_NUMBER_AUDIT.csv L50 ; "+WK+"t3_combine.json printed.*.nearest_11as",
  UA+"03_pixels/neighbours.py; "+UA+"referee/resolve_R3/r3_gaia_contrast.py", WK+"t3_combine.json", "", "no other source within $10''$",
  "Capability depends on the adopted depth: at the conservative-low delta = 400 ppm the required eclipse is 0.90 in the most favourable product (neighbour_census_summary.json capable 'dreq_A_lo_min<=1.0'), i.e. marginally capable. The Note's '>=43\"' claim holds for delta ~>420 ppm; the depth itself is another topic (T8).")
f("A.S6_diffimg", "PIXEL", "Usability of Sector 6 difference images for the 3.18 d signal", "not usable: all 6 transits in the data lie 0.41, 0.76, 1.80, 2.15, 3.19, 3.54 h from a momentum dump (7th transit in the mid-sector gap); two dumps fall in transit, four in the baseline windows",
  "dump times resolved to +-0.25 h (30-min FFI flags); raw S6 (-8.7, -34.7)\", dump-excised (+6.3, -1.3)\", decorrelated (+1.0, -1.9)\" with jackknife SE (10.6, 32.8)\"", "h", "FINAL",
  UA+"03_pixels/s6_dumps.json S6_ffi.transits ; "+UA+"03_pixels/verify/v3_s6.json ; "+UA+"referee/resolve_R3/logs/r2_combine.log 'S6' ; "+PP+"analysis/out/a5_summary.json localisation.A.per_set.S6_ffi.transit_to_dump_h",
  UA+"03_pixels/s6_dumps.py; "+UA+"03_pixels/verify/v3_s6.py", "", "C079 per-set S6 values (pubprep (+8.19, -2.39), audit z (2.39, -0.64), '7\" E 2.4 sigma') - none adopted",
  "Its Sector 6 difference images are unusable because every transit falls 0.4--3.5\\,h from a reaction-wheel desaturation.", "Dumps every 3.125 d ~ P_A, so the phase drifts only 1.39 h per dump cycle.")
f("A.centroids", "PIXEL", "Flux-weighted / SPOC centroid shifts during 3.18 d transits (supporting only)", "SPOC S32 MOM_CENTR1 z -1.66, MOM_CENTR2 z +0.12 (empirical null, 2000 sets); raw-pixel S98 -0.52 / -0.33 sigma, S6 -1.23 / -0.93 sigma", "", "sigma", "FINAL",
  AU+"FINAL_PARAMETERS.json signal_A.localisation.spoc2min_centroids.columns ; addenda/independent_verification_tic4206066_SEALED/out/t10_centroid.json S98.A_x/A_y, S06.A_x/A_y",
  AU+"scripts/s12c_clean_localisation.py; verif scripts/t10_centroid.py", "", "", "not in Note", "A 43.5\" source would move the centroid by ~1 mpix (verif report T10); centroid errors up to 2.6 mpix in P make them non-decisive (pixel_recheck REVIEW_1 L111).")
f("A.aperture_anomaly_bearing", "PIXEL", "Does the S98 core-pixel depth excess (T9, 227 +- 60 ppm) bear on localisation?", "No bearing on the >=43\" exclusion; a deeper core is the opposite of an off-target signature; the aperture test cannot distinguish a 0.5-px (~10\") offset", "", "", "FINAL",
  "addenda/tic4206066_pixel_recheck/REVIEW_1.md L104-106 ; addenda/independent_verification_tic4206066_SEALED/verification_report.md L30", "", "", "", "not in Note", "Depth anomaly itself handled by T8.")
f("B.offset", "PIXEL", "11.13 d source position relative to the target, S32 decorrelated + S98 clean29", "(+1.14, -1.22); S98 alone (+1.19, -1.02)", "(2.20, 2.93); S98 alone (2.29, 3.12)", "arcsec", "FINAL",
  UA+"referee/resolve_R3/logs/r2_combine.log 'B_S32PC+S98_auditor', 'B_S98only_auditor' ; "+UA+"proposed_revision/figure_src/figure_values_v2.json panel_b.combined.B",
  UA+"03_pixels/pixloc.py; "+UA+"referee/resolve_R3/r2_combine.py", WK+"t3_s98_B_driver.json",
  "pubprep (+1.41, +0.09) +- (3.27, 4.43), printed (+1+-3, +1+-4) - Dec '+1' is a transcription mix with the audit value; audit w/o S6 (+6.12, +1.56) +- (3.65, 3.99); audit all three (+0.96, +1.38) +- (2.32, 3.03); pubprep all three (-0.04, +0.38)",
  "$(+1\\pm3, -1\\pm3)''$",
  "All versions are within 2.2 sigma of the target. The audit's S98 B (+8.33, +0.10) reproduces in my rerun but drops to (+5.23, -0.24) with U's mask and its jackknife spans +2.1..+13.5\" E (t3_s98_B_driver.json); U's value is jackknife-stable (E +0.3..+2.0). S32 has one clean B event (2197.361); 2175.096 lies 2.3 h after a 3.18 d transit and during a strong pointing drift (pixel_recheck SUMMARY L73, L85).")
f("B.bright_neighbour_84as", "PIXEL", "T = 13.06 star (TIC 4206057, Gaia DR3 3220385243255181184) at 84\" for signal B", "can produce B only through PRF wings (required eclipse 0.19-0.61 depending on aperture); 27.5 sigma (+-3\") / >=37 sigma (U injection cov) from the B position", "", "arcsec / sigma", "FINAL",
  UA+"03_pixels/neighbour_census_summary.json capable 'dreq_B_min<=1.0' ; "+UA+"_structured/audit__03_pixels.json numeric row 'Gaia DR3 stars able to produce the 11.13 d depth' ; "+AU+"PUBLICATION_TRUTH.md L255-256", UA+"03_pixels/neighbours.py", "", "audit '59% eclipse, excluded by the localisation' (83.8\")", "not in Note", "")
f("PRF.method", "PIXEL", "Method facts stated in the Note", "TESS_PRF (Bell & Higgins 2022, ascl:2207.008) with MAST start_s0004 PRF files (valid S6/S32/S98); TESScut = Astrocut (Brasseur et al. 2019, ascl:1905.007); POS_CORR decorrelation applied to S32 SPOC 2-min only (TESScut POS_CORR is zero); errors = max(injection SD, jackknife SE, robust bootstrap SD) from 300 model-free injections into the real pixels", "", "", "FINAL",
  UA+"03_pixels/g1_prf_check.json ; "+UA+"_structured/resolve__R3.json analysis (S98 POS_CORR zero) ; https://ascl.net/2207.008 and https://ascl.net/1905.007 (accessed 2026-09-25T15:35Z)",
  UA+"03_pixels/g1_prf_check.py, inject.py", "", "", "SPOC (Sector 32, pointing-decorrelated) and TESScut \\citep[Sector 98;][]{Brasseur2019} difference images, fitted with the TESS pixel response function \\citep{Bell2022}, with errors from injections and transit-to-transit scatter", "No known-planet control of the PRF pipeline exists (audit UNRESOLVED_ISSUES.md U12).")
f("verif.T3_neighbour_count", "DATA", "Number of QLP S98 neighbour stars in the verification common-mode test (C133)", "220", "", "stars", "FINAL",
  "addenda/independent_verification_tic4206066_SEALED/out/t3_neighbours.json key n_neighbours_used ; out/t3_neighbours.csv (220 data rows)", "verif scripts/t3_neighbours.py L34", "", "219 (verification_report.md L24, typo)", "not in Note", "0 neighbours reach the target S/N either way.")

C = []
def c(id, q, mids, claims, trace, diag, rerun, res, reason, status):
    C.append(dict(id=id, quantity=q, master_ids=mids, claims=claims, trace=trace, diagnosis=diag, rerun=rerun, resolution=res, reason=reason, status=status))
def cl(s, v, e, m): return dict(source=s, value=v, evidence=e, method=m)
c("X-T3-01", "A headline offset", "C076",
  [cl("pubprep", "(+7.90, -8.20) +- (5.27, 6.84)", PP+"analysis/out/a5_summary.json localisation.A.combined_S32tpf_S98clean", "raw S32 TPF + S98 10 events, windowed baseline, injection covariances"),
   cl("audit", "(+2.51, -7.53) +- (5.68, 5.69)", AU+"FINAL_PARAMETERS.json signal_A.localisation.combined_without_S6", "raw S32 TPF + S98 9 events (FROZEN QLP mask), windowed baseline"),
   cl("U", "(+5.32 +- 2.66, +0.20 +- 2.37)", UA+"referee/resolve_R3/logs/r2_combine.log A_eventscatter_auditor_pos", "S32 POS_CORR-decorrelated + S98 9 events, box regression, max(inj, jackknife, bootstrap) errors")],
  "pubprep a2_prf_localisation.py + a5_summaries.py on TESScut/SPOC TPF; audit s12_pixels.py/s12c on same pixels; U pixloc/run_real/poscorr_decorr on identical pixel data (t3_cube_compare.json)",
  "flawed input weighting: raw S32 is uninformative (sigma 14-22\") so pubprep/audit headlines are essentially their S98 values, whose Dec is recipe-dependent; decorrelated S32 (sigma ~3\") dominates the Dec",
  WK+"t3_s98_driver.py (reproduces audit S98 real_px (0.05,-0.45) and U (+5.65,-0.32)); "+WK+"t3_combine.py: decorrelated S32 + any S98 recipe -> (+4.5..+6.3, -2.2..+0.2)",
  "(+5 +- 3, 0 +- 3)\"", "U's result uses the only well-constrained S32 set and is insensitive to the S98 recipe; the two older values rest on S98 alone", "RESOLVED")
c("X-T3-02", "A S32 2-min per-set z", "C077",
  [cl("audit", "z (0.746, -0.101)", AU+"out/s12c_clean_localisation.json per_set.A[S32_tpf2min].z", "raw"), cl("pubprep", "z (+0.774, -0.191)", PP+"analysis/out/a5_summary.json localisation.A.per_set.S32_tpf.z", "raw"),
   cl("U", "decorrelated (+5.21, +0.39) +- (3.08, 2.78)", UA+"03_pixels/summary.json sets.S32_tpf_poscorr_A", "POS_CORR decorrelation")],
  "same raw SPOC TPF, different codes", "both raw z values agree (<0.1) and are valid for raw data; superseded by the decorrelated set", "none needed (U verifier reproduces decorrelated value independently, verify/v1_loc.json)",
  "use decorrelated S32 (+5.2 +- 3.1, +0.4 +- 2.8)\"", "higher S/N (7.7 vs 3.5) and independent reproduction", "RESOLVED")
c("X-T3-03", "A S98 clean per-set offset and epoch count", "C078",
  [cl("audit", "z (0.07, -1.45); sky (+1.17, -8.41), 9 epochs", AU+"FINAL_PARAMETERS.json signal_A.localisation.S98_clean", "windowed baseline, FROZEN QLP mask"),
   cl("pubprep", "z (+1.05, -1.41); (+7.03, -9.29), 10 epochs", PP+"analysis/out/a5_summary.json localisation.A.per_set.S98_ffi", "windowed baseline, one-sided event allowed"),
   cl("U", "(+5.65, -0.32), 9 epochs", UA+"03_pixels/diagnose_s98.json S98_clean29_A", "box regression +-4 h"),
   cl("T3 rerun", "audit recipe + U mask (+9.31, -10.41) raw; U recipe + audit mask (+4.91, +3.46)", WK+"t3_s98_driver.json", "cross-swap")],
  "identical S98 pixels; epoch difference = BTJD 3995.555 (pubprep only)", "epochs: pubprep included an event with no post-transit baseline (gap 3995.608-3995.821) outside its injection calibration -> 9 is correct. Offset: recipe (baseline window/estimator) dependence, not mask; audit recipe is jackknife-unstable",
  WK+"t3_s98_driver.json, t3_gap3995.json", "epochs FINAL = 9; per-set S98 Dec ranges -10..0\" by recipe (UNRESOLVED per set); no effect on the headline (X-T3-01)",
  "cannot be settled without a common injection+jackknife calibration of both recipes; not needed for the Note", "UNRESOLVED")
c("X-T3-04", "A Sector 6 localisation", "C079",
  [cl("pubprep", "(+8.19, -2.39), z (+3.02, -0.58); text (+6.6, -3.4)", PP+"analysis/out/a5_summary.json localisation.A.per_set.S6_ffi", "raw"), cl("audit", "z (2.39, -0.64); with dump step z (2.51, -2.48)", AU+"FINAL_PARAMETERS.json signal_A.localisation.S6_dump_note", "raw / step"),
   cl("U", "raw (-8.7, -34.7); dump-excised (+6.3, -1.3); decorrelated (+1.0, -1.9) jackknife SE (10.6, 32.8)", UA+"referee/resolve_R3/logs/r2_combine.log S6", "three treatments")],
  "S6 TESScut 30-min; every transit 0.41-3.54 h from a dump", "values depend entirely on dump treatment; none reliable", "none needed", "S6 not used; no S6 value quoted or plotted", "all packages agree S6 is unusable; U shows no treatment is stable", "RESOLVED")
c("X-T3-05", "A combined offset including S6", "C080",
  [cl("audit", "(+6.32, -3.93) +- (2.68, 3.06), chi2 7.34", AU+"FINAL_PARAMETERS.json signal_A.localisation.combined_all", "raw S6"), cl("pubprep", "(+8.34, -3.69) +- (2.39, 3.15), chi2 14.55", PP+"analysis/out/a5_summary.json localisation.A.combined_all_independent", "raw S6"),
   cl("U", "(+5.06 +- 2.58, +0.19 +- 2.36), p 0.145 with decorrelated S6 and event errors", UA+"referee/resolve_R3/logs/r2_combine.log A_with_S6pc_eventscatter", "decorrelated S6")],
  "all include an unusable S6 set", "not a headline quantity; both older values carry the corrupted S6 point with injection-only errors", "none needed", "not adopted; Note excludes S6", "S6 excluded (X-T3-04)", "RESOLVED")
c("X-T3-06", "A S32 10-min FFI offset", "C081",
  [cl("pubprep", "(+9.75, +2.65) +- (10.44, 17.11)", PP+"analysis/out/a5_summary.json localisation.A.per_set.S32_ffi", "real minus injection mean"), cl("audit", "free-fit total 7.42\" (dRA 7.13\")", AU+"out/s12_pixels_S32_ffi.json S32_ffi_A.free_fit", "free fit, all data"),
   cl("U", "raw TESScut (+9.01, -8.24); TESS-SPOC FFI decorrelated (+4.88, -3.03) +- (3.37, 3.58)", UA+"03_pixels/summary.json sets.S32_ffi_A, S32_tspoc_poscorr_A", "raw / decorrelated")],
  "same S32 FFI photons", "different quantities (bias-corrected component vs total free-fit offset); raw values share the S32 pointing noise", "none needed",
  "not in the headline; decorrelated S32 FFI agrees with the adopted position", "consistent within errors", "NOT_A_CONFLICT")
c("X-T3-07", "sigma-distance of depth-capable stars", "C037",
  [cl("pubprep", ">=7 sigma (7.18; 6.99 with TIC J2000)", PP+"manuscript.tex L39 ; "+UA+"referee/resolve_R3/logs/r2_combine.log", "from (+8+-5,-8+-7)"), cl("audit", ">=8 sigma", AU+"PUBLICATION_TRUTH.md L156", "from (+2.5,-7.5)+-5.7 (T3 recompute 8.05)"),
   cl("U manuscript", ">=14 sigma", UA+"proposed_revision/manuscript.tex L41", "from (+5+-3, 0+-3), 14.2-14.3"), cl("U R3", ">10 sigma", UA+"_structured/resolve__R3.json recommended_manuscript_change", "conservative")],
  "Mahalanobis distances from different fitted positions", "stale positions for pubprep/audit", WK+"t3_combine.json printed", "formal minimum 14.2-14.3 sigma; print '>10 sigma'", "position superseded (X-T3-01); Gaussian tails not validated beyond ~10 sigma", "RESOLVED")
c("X-T3-08", "B combined S32+S98", "C119",
  [cl("pubprep", "(+1.41, +0.09) +- (3.27, 4.43); printed (+1+-3, +1+-4)", PP+"analysis/out/a5_summary.json localisation.B.combined_S32tpf_S98clean", "raw S32 2 events"), cl("audit", "(+6.12, +1.56) +- (3.65, 3.99)", AU+"FINAL_PARAMETERS.json signal_B.localisation.combined_without_S6", "audit recipe"),
   cl("U", "(+1.14 +- 2.20, -1.22 +- 2.93)", UA+"referee/resolve_R3/logs/r2_combine.log B_S32PC+S98_auditor", "decorrelated S32 + S98")],
  "audit S98 B (+8.33, +0.10) drives its value", "recipe/mask sensitivity of the audit S98 estimator (clean29 mask -> +5.2; jackknife +2.1..+13.5)", WK+"t3_s98_B_driver.json", "(+1 +- 3, -1 +- 3)\"", "U value jackknife-stable and reproduced by the verifier code (S98 (+1.16, -0.84)); all values consistent with target", "RESOLVED")
c("X-T3-09", "B combined all three sectors", "C120",
  [cl("audit", "(+0.96, +1.38) +- (2.32, 3.03)", AU+"FINAL_PARAMETERS.json signal_B.localisation.combined_all", "incl. S6"), cl("pubprep", "(-0.04, +0.38) +- (2.21, 3.21); METHOD_AUDIT text (+0.2, +0.6)", PP+"analysis/out/a5_summary.json localisation.B.combined_all_independent", "incl. S6")],
  "different S98/S32 inputs; pubprep text vs a5 mismatch is a stale text value", "not adopted (S6 excluded from Note)", "none needed", "not used; all consistent with target", "Note quotes S32+S98 (or S98) only", "RESOLVED")
c("X-T3-10", "B per-set z S98 clean", "C121",
  [cl("audit", "(+2.14, -0.19)", AU+"PUBLICATION_TRUTH.md L253", "audit recipe"), cl("pubprep", "(+0.76, -0.12)", PP+"analysis/out/a5_summary.json localisation.B.per_set.S98_ffi.z", "pubprep recipe"), cl("U", "(+1.19, -1.02) raw, inj mean (+1.07, -1.08)", UA+"03_pixels/summary.json sets.S98_clean29_B", "U recipe")],
  "same 3 events", "audit estimator mask- and event-sensitive (T3 rerun)", WK+"t3_s98_B_driver.json", "U S98 B adopted", "stable estimator; all within ~2 sigma", "RESOLVED")
c("X-T3-11", "B per-set z S6", "C122", [cl("audit", "(-0.84, +0.34)", AU+"PUBLICATION_TRUTH.md L253", ""), cl("pubprep", "(-0.55, +0.41)", PP+"analysis/out/a5_summary.json localisation.B.per_set.S6_ffi.z", "")],
  "independent codes, 2 events", "noise-level difference, both on target", "none needed", "not used in Note", "S6 excluded", "COSMETIC")
c("X-T3-12", "B per-set z S32 2-min", "C123", [cl("audit", "(-0.77, +1.5)", AU+"PUBLICATION_TRUTH.md L253", "2 events"), cl("pubprep", "(-0.99, +1.22)", PP+"analysis/out/a5_summary.json localisation.B.per_set.S32_tpf.z", "2 events"), cl("U verifier", "1 clean event (2197.361)", UA+"referee/resolve_R3/logs/r5_s98_gap.log", "acceptance rules")],
  "event 2175.096 compromised (adjacent 3.18 d transit, pointing drift)", "both values rest on a compromised event", "none needed", "S32 B point not plotted or quoted", "only one clean event", "COSMETIC")
c("X-T3-13", "verification T3 neighbour count", "C133", [cl("verif", "220", "addenda/independent_verification_tic4206066_SEALED/out/t3_neighbours.json n_neighbours_used", "snr_A notna count"), cl("verif report", "219", "addenda/independent_verification_tic4206066_SEALED/verification_report.md L24", "text")],
  "CSV has 220 data rows; target not included", "transcription error in report", "none needed (row count checked)", "220", "JSON/CSV outrank narrative", "COSMETIC")
c("X-T3-14", "U internal: sigma wording for capable stars", "",
  [cl("U REQUIRED_CHANGES", "'>=43\" / >=7 sigma' checked, not changed", UA+"REQUIRED_CHANGES.md L60", ""), cl("U manuscript", ">=14 sigma", UA+"proposed_revision/manuscript.tex L41", ""), cl("U R3", ">10 sigma", UA+"_structured/resolve__R3.json", ""), cl("U figure values", "text '>= 43\" (>= 7 sigma)'", UA+"proposed_revision/figure_src/figure_values_v2.json manuscript_consistency.b_capable_stars.text", "stale")],
  "pass-2 edit changed the text after REQUIRED_CHANGES and figure consistency block were written", "stale documents", "none needed", "print '>10 sigma'", "see X-T3-07", "RESOLVED")
c("X-T3-15", "U internal: cause of the old -8\" Dec", "",
  [cl("U REQUIRED_CHANGES B1", "S32 pointing / S98 gap-edge artefact", UA+"REQUIRED_CHANGES.md L26", ""), cl("U verifier", "POS_CORR2 correlation 2.0 sigma (not 2.7), ~1.1 unweighted", UA+"03_pixels/verify/VERIFY_REPORT.md 'Issues the auditor missed'", ""), cl("T3 rerun", "audit S98 gives -8.4\" without the gap-edge event", WK+"t3_s98_driver.json AUD_frozen", "")],
  "", "attribution partly wrong: gap-edge event explains pubprep's S98 only; audit's -8\" comes from its S98 estimator; S32 raw effect is jitter noise plus event 2175.001", WK+"t3_s98_driver.py", "old -8\" = S98-dominated combination (raw S32 uninformative) with a recipe-dependent S98 Dec", "rerun evidence", "RESOLVED")
c("X-T3-16", "U internal: S6 'unusable' vs 'unreliable'; error size +-2 vs +-3", "",
  [cl("U 03 auditor", "'unreliable'; (+5+-2, 0+-2)", UA+"03_pixels/REPORT.md F1, F3", ""), cl("U verifier/R3/manuscript", "'unusable'; (+5+-3, 0+-3)", UA+"_structured/resolve__R3.json ; "+UA+"proposed_revision/manuscript.tex L41", "")],
  "", "injection-only errors understate event scatter (S98 E 5.3 vs 3.3)", "none needed", "'unusable'; +-3\"", "no S6 treatment is stable; event-scatter errors are larger", "RESOLVED")
c("X-T3-17", "U vs pixel_recheck on localisation", "",
  [cl("P", "localisation statements at pubprep l.39 'NONE' (checked vs (+8+-5,-8+-7), (+1+-3,+1+-4))", "addenda/tic4206066_pixel_recheck/REVIEW_1.md L124", "no new localisation"), cl("U", "(+5+-3, 0+-3), (+1+-3, -1+-3)", UA+"proposed_revision/manuscript.tex L37, L41", "")],
  "", "P does not localise; its check is against the superseded pubprep numbers; P confirms 0.4-3.5 h and that the aperture test cannot resolve 0.5-px offsets; P and U agree S32 B event 2175.096 (epoch -164) is compromised", "none needed", "no conflict in substance", "different quantities", "NOT_A_CONFLICT")

MF = [
 dict(quote="place the 3.18\\,d source at $(+5\\pm3, 0\\pm3)''$ from the target ($\\approx2\\sigma$)", problem="'(≈2σ)' is ambiguous (reads like a detection significance). The value is z_E = 2.0 with exact errors (1.7 with the printed ±3\"), 2-D p = 0.135; the source is consistent with the target.",
      fix="place the 3.18\\,d source at $(+5\\pm3, 0\\pm3)''$ from the target, consistent with it within $\\approx2\\sigma$", severity="MINOR"),
 dict(quote="All Gaia DR3 stars able to produce the depth lie $\\ge43''$ ($\\ge14\\sigma$) away.", problem="14σ is a Gaussian extrapolation of ±3\" errors (formal minimum 14.2–14.3σ, TIC 4206077); estimator tails are not Gaussian. Injection tests support >10σ-level exclusion.",
      fix="All Gaia DR3 stars able to produce the depth lie $\\ge43''$ ($>10\\sigma$) away.", severity="MINOR"),
 dict(quote="no non-single-star solution or other source within $11''$", problem="The nearest Gaia DR3 source is at 11.06\" at J2016.0 and 10.6–10.9\" at the TESS epochs, so '11″' is borderline.",
      fix="no non-single-star solution or other source within $10''$", severity="MINOR"),
 dict(quote="and Gaia DR3 stars able to produce the 3.18\\,d depth (crosses).", problem="One capable star, TIC 4206057 (T = 13.06, 84\"), is not plotted (figure_values_v2.json panel_b.neighbours[4].plotted = false), so the caption implies all capable stars are shown.",
      fix="and Gaia DR3 stars able to produce the 3.18\\,d depth (crosses; a fifth, $T=13.1$ at $84''$, lies outside the panel).", severity="MINOR"),
 dict(quote="SPOC (Sector 32, pointing-decorrelated) and TESScut \\citep[Sector 98;][]{Brasseur2019} difference images, fitted with the TESS pixel response function \\citep{Bell2022}, with errors from injections and transit-to-transit scatter", problem="Matches the method (decorrelation applied to S32 only; errors = max of injection SD, jackknife, bootstrap).", fix="none", severity="NONE"),
 dict(quote="Its Sector 6 difference images are unusable because every transit falls 0.4--3.5\\,h from a reaction-wheel desaturation.", problem="Confirmed: 0.41–3.54 h for all 6 transits in the data; no S6 treatment is stable.", fix="none", severity="NONE"),
 dict(quote="offset $(+1\\pm3, -1\\pm3)''$", problem="Matches (+1.14 ± 2.20, −1.22 ± 2.93)\" (±3 on E is conservative rounding of 2.2).", fix="none", severity="NONE"),
]
UNK = [
 dict(item="S98 per-set 3.18 d Dec: recipe dependence (-10 to 0\" for the same 9 events)", why="Windowed-baseline (audit/pubprep) and continuous box-regression (U) estimators disagree; audit estimator is jackknife-unstable; no common calibration of both recipes exists",
      would_resolve="Run both recipes with the same model-free injections and jackknife on the same 9 events; compare the injection-calibrated spread with the observed difference", blocks_submission=False),
 dict(item="Known-planet control of the PRF localisation pipeline (audit U12)", why="Only injections into this star's own pixels were used (including off-target shifted sources recovered at 0.99–1.02 of the true offset)",
      would_resolve="Run the same pipeline on a confirmed planet host on the same camera/CCD in S32/S98", blocks_submission=False),
 dict(item="Sources within ~5–10\" of the target and the ~5\" eastward offset (≈2σ)", why="Pixel fits cannot exclude a source within ~5–10\" (U 03 F6); Gaia misses depth-capable companions only inside ~2.7–4.5\" (resolve_R3); the eastward offset is not significant and partly method bias (+2\" E in S98 injections)",
      would_resolve="High-resolution imaging (already named in the Note)", blocks_submission=False),
 dict(item="Whether the 11\" Gaia neighbour (T 18.13) can produce the depth", why="Required eclipse 1.20–1.23 at δ = 500 ppm but ~0.90 at δ = 400 ppm; depends on the adopted dilution-corrected depth (topic T8)",
      would_resolve="Fix the A depth (T8); it is also 4–5σ from the fitted position", blocks_submission=False),
]
FW = [W + x for x in ("t3_cube_compare.py", "t3_cube_compare.json", "t3_s98_driver.py", "t3_s98_driver.json", "t3_s98_driver.log", "t3_s98_B_driver.py", "t3_s98_B_driver.json", "t3_s98_B_driver.log", "t3_combine.py", "t3_combine.json", "t3_combine.log", "t3_gap3995.py", "t3_gap3995.json", "t3_write_result.py", "result.json")]
SUM = ("FINAL: the 3.18 d source lies at (+5.32 ± 2.66, +0.20 ± 2.37)\", printed (+5 ± 3, 0 ± 3)\", consistent with the target within ≈2σ (z_E 2.0, 2-D p = 0.135). "
 "The old -8\" Dec in pubprep (+7.9, -8.2) and audit (+2.5, -7.5) came from S98 alone: raw S32 has 14–22\" errors and carries no weight. "
 "Once S32 is pointing-decorrelated (independently reproduced), every S98 recipe gives a combined (+4.5..+6.3, -2.2..+0.2)\". "
 "My rerun on identical pixels shows the S98 per-set Dec depends on the recipe (audit code: -8.4\" for the same 9 events, jackknife-unstable; U code: -0.3\", stable). "
 "So U's 'gap-edge artefact' explains pubprep but not the audit; the S98 per-set value stays UNRESOLVED, but it does not affect the headline. "
 "Nine S98 events are correct: 3995.555 is followed by a 5.1-h data gap. "
 "Depth-capable Gaia stars lie ≥43\" away (TIC 4206077 nearest), 14.2σ formally; I recommend printing '>10σ'. "
 "S6 is unusable (0.41–3.54 h from dumps; no stable treatment). "
 "11.13 d: (+1 ± 3, −1 ± 3)\", on target in all versions; pubprep's '+1 ± 4' Dec was a transcription error. "
 "Manuscript: four MINOR wording fixes (≈2σ phrasing, >10σ, 'within 10\"', caption for the off-panel 84\" star). "
 "Open, non-blocking: S98 recipe dependence, no known-planet control, sub-10\" blends, capability of the 11\" star at low depth.")
R = dict(topic="T3_pixels_localisation", summary=SUM, facts=F, conflicts=C, manuscript_findings=MF, unknowns=UNK, files_written=FW)
json.dump(R, open(W + "result.json", "w"), indent=1, ensure_ascii=False)
print(len(SUM.split()), "words;", len(F), "facts;", len(C), "conflicts")
