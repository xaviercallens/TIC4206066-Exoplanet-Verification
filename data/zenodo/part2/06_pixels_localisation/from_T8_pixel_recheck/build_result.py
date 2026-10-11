"""Assemble result.json for topic T8_pixel_recheck (values copied from the rerun outputs listed in each evidence field)."""
import json, os
W = os.path.dirname(os.path.abspath(__file__))
P = "addenda/tic4206066_pixel_recheck"
K = "addenda/tic4206066_FINAL_RECONCILIATION/_work/T8_pixel_recheck"
V = "addenda/independent_verification_tic4206066_SEALED"
U = "addenda/opus55_final_audit_tic4206066"
M = "addenda/tic4206066_MASTER_EVIDENCE"
MS = f"{U}/proposed_revision/manuscript.tex"


def F(id_, group, q, value, unc, units, status, ev, code, rerun, sup, mform, notes):
    return dict(id=id_, group=group, quantity=q, value=value, uncertainty=unc, units=units, status=status, evidence=ev,
                code=code, rerun=rerun, supersedes=sup, manuscript_form=mform, notes=notes)


facts = [
 F("PIX.reproduction", "DATA", "Reproducibility of the pixel re-check package P (v3) from its own scripts",
   "all 4 science JSONs and 3 PNGs byte-identical; inputs.json differs only in the recorded interpreter path", "", "", "FINAL",
   f"{K}/rerun_log.txt ; {K}/pkg_copy/out/*.json vs {P}/out/*.json (cmp) ; {P}/SHA256SUMS (verified OK before copy)",
   f"{P}/scripts/s01..s05 (copy with absolute ROOT in {K}/pkg_copy/scripts/common.py)", f"{K}/pkg_copy/out/", "",
   "not in Note", "s06_hashes not run. The verif __pycache__ listing was unchanged after the s03b import."),
 F("PIX.v3_fixes_review", "PIXEL", "Re-review of the post-REVIEW_2 fixes (gap/dump event rule, epoch 8 removal, s04 A-masking caveat)",
   "Fixes correct and applied. S98: epoch 8 (BTJD 4024.199) is removed by the gap rule: a 5.44-h gap from +2.55 h to +7.99 h after mid-transit contains a momentum dump at +7.87 h. Epoch -1 (BTJD 3995.554; 5.22-h gap; dump 6.4 h away) is also removed by the gap rule. Epoch 6 has <70% in-transit coverage. Epochs 3, 12 and 15 have no clean data in the window. The dump rule removes no additional S98 or S32 event; it removes all 5 fittable S6 events (0.43-3.53 h from dumps). Removing epoch -1 is immaterial (+22.3 +- 42.9 -> +19.5 +- 48.6 ppm). s04 now masks A transits in the B local baseline, and SUMMARY labels that row 'unstable under drift'.",
   "", "ppm", "FINAL", f"{K}/out/t8_a_events_and_nulls.json keys S98_ffi.bookkeeping, S6_ffi.bookkeeping ; {K}/out/t8_s98_rule_sensitivity.json keys v3_8, plus_epoch_-1 ; {K}/logs_epoch8.txt ; {P}/scripts/s03_after_decorrelated.py L26-L28 ; {P}/scripts/s04_b_s32_events.py L19, L26",
   f"{K}/t8_a_events_and_nulls.py ; {K}/t8_s98_rule_sensitivity.py", f"{K}/out/t8_a_events_and_nulls.json", "REVIEW_2 table row C4 'Not fixed' (it described v2)", "not in Note", ""),
 F("PIX.B_S32_offtarget", "SIGNAL_B", "Signal B, S32 2-min TPF: apparent off-target (depth grows with aperture) pattern of the first pass",
   "Methodological false positive (confirmed). Per-event quadratic-in-time fits give flat depths with or without centroid terms. Epoch -164: r1 951 raw / 924 +- 209 corrected; r3 895 raw / 949 +- 151 corrected. Epoch -162: r1 1071 / 1061 +- 256; r3 1016 / 999 +- 176 ppm. Combined r1-r3 = +8 +- 120 ppm (white errors). The on-target PRF prediction is r1/r3 = 1.02; a source at TIC 4206077 (43'') predicts r1/r3 = 0.12.",
   "white errors are lower bounds (S98 red-noise factor ~1.6 in verif T5); only 2 events", "ppm", "FINAL",
   f"{P}/out/after_decorrelated.json B.S32_tpf2min.events[*].raw / .corrected ; {K}/out/t8_b_and_paired.json key B_diff_1.0-3.0_combined_white ; {K}/out/t8_prf_offtarget_signature.json S32_tpf2min.*.depth_ratio_vs_r2 ; {K}/out/t8_b_drift.json (epoch -164: core range 3.2% in 30-min bins, y centroid range 38 mpix; -162: 0.8%, 4 mpix)",
   f"{P}/scripts/s03_after_decorrelated.py ; {K}/t8_b_and_paired.py ; {K}/t8_prf_offtarget_signature.py", f"{K}/out/t8_b_and_paired.json",
   "first-pass 'B missing from the core' (P/out/before_global_detrend.json)", "not in Note",
   "The first-pass pattern comes from the global 0.75-d running median (s02), and from a constant local-median baseline (s04), neither of which follows the ~38-mpix y drift in epoch -164. Because the time-only (raw) fit also flattens the pattern, the flattening is not an artefact of centroid regressors absorbing a real blend signal (REVIEW_1 C6)."),
 F("PIX.B_S32_firstpass", "SIGNAL_B", "First-pass (flawed) B depths vs aperture, S32 2-min, global detrend",
   "r0.5 -1619 +- 580; r1 336 +- 239; r2 740 +- 165; r3 1075 +- 155 ppm; r1-r3 -739 +- 159 (4.6 sigma); r2-r3 -336 +- 75 (4.5 sigma); r0.5-r2 -2358 +- 472 (5.0 sigma)",
   "random-phase null SD (400 trials, seed 1)", "ppm", "FINAL", f"{P}/out/before_global_detrend.json S32_tpf2min.apertures.*.B and .diffs.*.B ; {K}/out/t8_b_and_paired.json first_pass_*",
   f"{P}/scripts/s02_before_global_detrend.py", f"{K}/pkg_copy/out/before_global_detrend.json (identical)", "", "not in Note",
   "SUMMARY's 'core' in the first pass is r = 0.5 px (single pixel); its per-event table compares r = 1 px (5 pixels) with r = 3 px. These are different apertures."),
 F("PIX.B_S32_perevent_depths", "SIGNAL_B", "B per-event corrected depths r=1 / r=3 px, S32 2-min",
   "epoch -164: 924 +- 209 / 949 +- 151; epoch -162: 1061 +- 256 / 999 +- 176", "formal white errors", "ppm", "FINAL",
   f"{P}/out/after_decorrelated.json B.S32_tpf2min.events[0,1].corrected['1.0','3.0']", f"{P}/scripts/s03_after_decorrelated.py", f"{K}/pkg_copy/out/after_decorrelated.json", "", "not in Note",
   "Matches SUMMARY and the user's numbers exactly."),
 F("PIX.B_S32_r2_minus_r3", "SIGNAL_B", "B r2-r3 aperture difference, S32 2-min, per-event method",
   "-81 +- 59 ppm (white error of the difference light curve; 1.4 sigma); per event -61 +- 81, -103 +- 85",
   "the JSON's +-21 (err_empirical) is the scatter of N = 2 values and is not a usable error", "ppm", "FINAL",
   f"{K}/out/t8_b_and_paired.json key B_diff_2.0-3.0_combined_white ; {P}/out/after_decorrelated.json B.S32_tpf2min.summary.corrected.diff_2.0-3.0 = -82.0 +- 20.6",
   f"{K}/t8_b_and_paired.py", f"{K}/out/t8_b_and_paired.json", f"{P}/out/after_decorrelated.json diff_2.0-3.0 err_empirical 20.6 (nominal 4 sigma)", "not in Note",
   "On-target PRF prediction +15 ppm; a 43'' source predicts -728 ppm. No off-target signature."),
 F("PIX.A_S98_T9method", "SIGNAL_A", "T9 S98 r<=1 minus r<=2 px extra depth, T9's own 11x11 cut-out and QUALITY==0 selection (published T9 value)",
   "226.6 ppm (reproduced exactly); leave-epoch-8-out: +4.0 +- 70.7 ppm", "random-phase null SD 60.1 ppm (3.77 sigma)", "ppm", "FINAL",
   f"{V}/out/t9_aperture_diff.json S98.A_r1.0_minus_r2.0 ; {P}/out/t9_method_replay.json S98_T9_own_11x11_cutout.A_r1.0_minus_r2.0 ; {K}/out/t8_s6_and_clipping.json S98.own11x11_T9sel, S98.own11x11_leave_one_out['8']",
   f"{V}/scripts/t9_aperture_diff.py ; {K}/t8_s6_and_clipping.py", f"{K}/out/t8_s6_and_clipping.json",
   "verif T9 'unexplained, open' and M C092 UNRESOLVED", "not in Note",
   "The whole excess comes from epoch 8: 29 in-transit points at +3819 ppm (8x the transit depth) out of 497. Every other leave-one-out value lies between 175 and 315 ppm. Mechanism (not proven): the epoch-8 transit ends 1.6 h before a 5.4-h data gap that contains a momentum dump, so T9's 1-d running median, segmented at gaps > 0.1 d, is one-sided there."),
 F("PIX.A_S98_T9method_clean", "SIGNAL_A", "T9 method, S98 21x21, clean cadence selection (P s03b)",
   "199.8 ppm over 10 events (epochs -1,0,1,4,5,8,9,10,13,14); epoch 8 contributes +336 ppm of it", "null SD 73.5 (2.72 sigma); white 51.0; event scatter 405", "ppm", "FINAL",
   f"{P}/out/t9_method_replay.json S98_ffi_clean_selection.A_r1.0_minus_r2.0 ; {K}/out/t8_t9_decomposition.json S98_clean_T9eph.baseline_fullT14",
   f"{P}/scripts/s03b_t9_method_replay.py ; {K}/t8_t9_decomposition.py", f"{K}/out/t8_t9_decomposition.json", "", "not in Note",
   "SUMMARY calls this 'the same clean data' as the per-event method. It shares the cadences but not the events: it includes epochs 8 and -1."),
 F("PIX.A_S98_T9method_same8", "SIGNAL_A", "T9 method restricted to the 8 events used by the per-event method",
   "-73.9 ppm", "null SD 105.7; event scatter 54.2; white 56.6", "ppm", "FINAL",
   f"{K}/out/t8_t9_decomposition.json S98_clean_T9eph.restricted_to_clean8_events (audit ephemeris: -74.9 +- 105.1)",
   f"{K}/t8_t9_decomposition.py", f"{K}/out/t8_t9_decomposition.json", "", "not in Note", ""),
 F("PIX.A_S98_perevent", "SIGNAL_A", "Per-event local-correction method, S98 r<=1 minus r<=2 px, 8 clean events",
   "+19.5 ppm (median +6.2)", "event scatter 48.6; white 56.8; random-phase null 75.4 (null mean +7.0; 41% of null trials >= observed)", "ppm", "FINAL",
   f"{P}/out/after_decorrelated.json A.S98_ffi.summary.corrected.diff_1.0-2.0 ; {K}/out/t8_a_events_and_nulls.json S98_ffi.summary_diff_1.0-2.0, S98_ffi.random_phase_null",
   f"{P}/scripts/s03_after_decorrelated.py ; {K}/t8_a_events_and_nulls.py", f"{K}/out/t8_a_events_and_nulls.json",
   "v2 value +230 +- 211 (N=10, epoch 8 +2098 ppm included) in REVIEW_2", "not in Note",
   "With the T9-like full-T14 in-transit window the value is -1.8 +- 50.3 ppm, so the in-transit window does not drive the difference."),
 F("PIX.A_S98_method_difference", "SIGNAL_A", "Difference between the T9 method and the per-event method in S98",
   "Not a method effect. On their own selections: 199.8 - 19.5 = 180 ppm (about 2.0 sigma if the errors 73.5 and 48.6 are treated as independent). On identical events (same 8): paired difference -93 +- 63 ppm (1.5 sigma).",
   "paired error from event scatter", "ppm", "FINAL",
   f"{K}/out/t8_b_and_paired.json S98_paired_T9method_minus_perevent_same8 ; {K}/out/t8_t9_decomposition.json S98_clean_T9eph",
   f"{K}/t8_t9_decomposition.py ; {K}/t8_b_and_paired.py", f"{K}/out/t8_b_and_paired.json",
   "SUMMARY/REVIEW_2 'two methods differ by ~180 ppm (~2 sigma), cause not identified'; user (d)", "not in Note",
   "The ~180 ppm comes from including epoch 8 (and, with the opposite sign, epoch -1) in the T9 cadence selection."),
 F("PIX.A_S98_error_models", "SIGNAL_A", "Significance of the S98 core excess under different error models",
   "T9 published 226.6: null 3.77 sigma (3.0 sigma after subtracting the null mean +45.5), white 5.0 sigma, event scatter 0.80 sigma (283 ppm), without epoch 8 0.06 sigma. T9 clean 199.8: null 2.72, white 3.92, event scatter 0.49 sigma. Per-event 19.5: event scatter 0.40, white 0.34, null 0.26 sigma.",
   "", "sigma", "FINAL", f"{K}/out/t8_s6_and_clipping.json S98.* ; {K}/out/t8_t9_decomposition.json ; {K}/out/t8_a_events_and_nulls.json S98_ffi",
   f"{K}/t8_t9_decomposition.py ; {K}/t8_a_events_and_nulls.py ; {K}/t8_s6_and_clipping.py", f"{K}/out/", "", "not in Note",
   "The error model matters only for the epoch-8-contaminated statistic. The white and random-phase errors assume no single-event outlier."),
 F("PIX.A_S32", "SIGNAL_A", "S32 A aperture difference r<=1 minus r<=2 px",
   "per-event TPF -61.9 (event scatter 60.5; white 49.0; null 85.9); per-event FFI -36.2 +- 65.2; T9 method FFI QUALITY==0 +84.7 +- 135.1, FFI clean +81.9 +- 163.3, TPF +108.9 +- 117.9 (null SD)",
   "", "ppm", "FINAL",
   f"{P}/out/after_decorrelated.json A.S32_tpf2min / A.S32_ffi summary.corrected.diff_1.0-2.0 ; {K}/out/t8_a_events_and_nulls.json S32_tpf2min.random_phase_null ; {K}/out/t8_t9_decomposition.json S32_*_T9method",
   f"{K}/t8_t9_decomposition.py ; {K}/t8_a_events_and_nulls.py", f"{K}/out/t8_t9_decomposition.json",
   "REVIEW_2 scratch numbers +85 +- 135, +82 +- 163, +109 +- 118 (now reproduced exactly)", "not in Note",
   "No significant aperture dependence under any method (all < 1.1 sigma). With a full-T14 window the per-event TPF value is -87 +- 46 (event scatter)."),
 F("PIX.A_S6", "SIGNAL_A", "S6 A aperture difference r<=1 minus r<=2 px",
   "T9 method: 149.5 +- 64.8 (11x11, re-seeded null; T9 json 65.5); 232.4 +- 52.4 (21x21, T9 selection); 228.8 +- 50.3 (21x21, clean). Leave-epoch-(-794)-out on 11x11: 25.5 +- 64.3; on 21x21 clean without -794: ~101 ppm (from the stored per-event values). Per-event method: no usable event (all 0.43-3.53 h from dumps); without the dump rule +102 +- 127 (4 events).",
   "event-scatter errors for the T9 method about 129-130 ppm", "ppm", "FINAL",
   f"{K}/out/t8_s6_and_clipping.json S06.* ; {P}/out/t9_method_replay.json S6_ffi_* / S06_T9_own_11x11_cutout ; {P}/out/after_decorrelated.json A_sensitivity_no_dump_rule.S6_ffi",
   f"{K}/t8_s6_and_clipping.py", f"{K}/out/t8_s6_and_clipping.json", "", "not in Note",
   "REVIEW_2's annulus-clipping explanation is reproduced exactly: a 21x21 cube with the annulus clipped to the 11x11 stamp gives 149.5 +- 64.8 (S6) and 226.6 +- 60.1 (S98). S6 transits all lie near momentum dumps and cannot test the effect."),
 F("PIX.A_combined_S32_S98", "SIGNAL_A", "Per-event r<=1 minus r<=2 px, S32 TPF + S98 combined (independent photons)",
   "-12.4 +- 37.9 ppm (package, event-scatter weights); -15.9 +- 56.7 ppm with random-phase-null errors (my combination of 19.5 +- 75.4 and -61.9 +- 85.9)",
   "", "ppm", "FINAL", f"{P}/out/after_decorrelated.json A_combined.corrected.diff_1.0-2.0 ; {K}/out/t8_a_events_and_nulls.json *.random_phase_null.std_of_trial_means",
   f"{P}/scripts/s03_after_decorrelated.py", f"{K}/out/t8_a_events_and_nulls.json", "REVIEW_2 v2 combined -15 +- 53", "not in Note", "A 200 ppm core excess is 3.8 sigma above the conservative combination."),
 F("PIX.A_aperture_status", "SIGNAL_A", "Status of the Signal A central/aperture effect",
   "Not confirmed. The only >3 sigma evidence (T9 S98, 3.77 sigma) comes from one event (epoch 8, next to a 5.4-h gap with a momentum dump). Without it, S98 gives +4 +- 71 (T9 method) or +20 +- 75 (per-event method). A T9-sized excess (~200 ppm) is disfavoured at about 3.8 sigma by S98+S32 (-16 +- 57 ppm); effects below ~100 ppm are not excluded. S6 (dump-affected) retains 100-230 ppm with the T9 method and cannot be tested cleanly.",
   "", "", "FINAL", "see PIX.A_S98_T9method, PIX.A_S98_perevent, PIX.A_S32, PIX.A_S6, PIX.A_combined_S32_S98",
   "", f"{K}/out/", "user summary (b) 'not confirmed and not excluded'; P SUMMARY conclusion 1", "not in Note", ""),
 F("PIX.PRF_offtarget_signature", "FALSE_POSITIVE", "Aperture signature of an off-target source in P's apertures (official TESS PRF, Gaia DR3/TIC 8.2 scene)",
   "Source at TIC 4206077 (42.8-43.3''): depth ratio r1/r2 = 0.43 (S98), 0.22 (S32 TPF), 0.31 (S6); r3/r2 = 1.31, 1.81, 1.77. Predicted r1-r2 for a 470-ppm r2 depth: -268 (S98), -364 (S32 TPF), -323 (S6) ppm. On-target signal: r1/r2 = 1.002-1.004 (+1 to +2 ppm).",
   "model; PRF at the nearest MAST grid position", "ppm", "FINAL", f"{K}/out/t8_prf_offtarget_signature.json",
   f"{K}/t8_prf_offtarget_signature.py (uses {U}/03_pixels/common.py PRFModel read-only, PRF files {U}/00_data/prf/start_s0004)", f"{K}/out/t8_prf_offtarget_signature.json", "", "not in Note",
   f"Confirms that an off-target source makes the core SHALLOWER, so a 'core deeper' excess is not an off-target signature (also QLP DRN004 footnote, {U}/referee/resolve_R-DIL/QLP_DRN004.txt L59). The measured per-event values disfavour the 43'' star independently: S98 +20 +- 75 vs -268; S32 -62 +- 86 vs -364."),
 F("PIX.U_localisation_epoch8", "PIXEL", "Effect of the anomalous S98 epoch 8 on U's S98 A localisation",
   "U's S98 clean29 set (9 events) includes epoch 8 (BTJD 4024.20). U's own jackknife without it gives (+7.0, +0.3)'' against the full-set (+5.6, -0.3)'' (S98 alone).",
   "S98-only errors about +-3.3, +-4.5''", "arcsec", "FINAL",
   f"{U}/03_pixels/diagnose_s98.json S98_clean29_A.jackknife (dropped_tc 4024.200) ; {U}/03_pixels/real_results.json S98_clean29.A.event_tc, off_E, off_N",
   f"{U}/03_pixels/diagnose_s98.py", "", "", "not in Note",
   "Not material: a 1.4'' shift, and in the direction away from the target-consistent solution, not towards it."),
 F("PIX.manuscript_impact", "MANUSCRIPT_META", "Does the pixel re-check change any claim in U/proposed_revision/manuscript.tex?",
   "No. The Note makes no aperture-dependence, core-pixel or T9 claim. Its localisation, S6-dump, B-offset and >=43'' neighbour statements are consistent with P and with the reruns here. Do not promote this into the Note.",
   "", "", "FINAL", f"{MS} L22, L37, L41 ; {K}/out/*",
   "", "", "P SUMMARY 'Manuscript impact' (assessed the pubprep manuscript, not U's revision)", "no change",
   "P used the audit ephemeris (not U's final one) and box depths, not model fits. Swapping T9's ephemeris for the audit one moves the T9-method S98 value by 7 ppm (199.8 -> 206.9), and the box windows are +-28 min, so this does not matter. An optional record for the Zenodo or limitations notes: the verif T9 3.8 sigma S98 aperture excess traces to one event (epoch 8)."),
]

conflicts = [
 dict(id="X-T8-01", quantity="S98 Signal A core (r<=1) minus r<=2 px aperture depth excess", master_ids="C092",
      claims=[dict(source="verif", value="226.6 +- 60.1 ppm (3.77 sigma), 'unexplained'", evidence=f"{V}/out/t9_aperture_diff.json S98.A_r1.0_minus_r2.0 ; {V}/verification_report.md L30", method="T9: 11x11 TESScut, annulus bkg, global 1-d running median, full-T14 box, random-phase null"),
              dict(source="M", value="C092 UNRESOLVED: T9 226.6 +- 60.1 vs deep S32 TPF 316 (r0.7) rising to 452 (r2.0) ppm", evidence=f"{M}/CONTRADICTIONS.md L711-L721", method="compilation"),
              dict(source="P", value="per-event +20 +- 49 (8 events); T9 method on clean data 200 +- 73; 'not confirmed, not excluded'", evidence=f"{P}/SUMMARY.md Results/Conclusions ; {P}/out/after_decorrelated.json ; {P}/out/t9_method_replay.json", method="per-event local fit with centroid/background terms; T9 replay"),
              dict(source="other", value="T8 rerun: T9 published config without epoch 8 +4.0 +- 70.7; T9 method on the same 8 events -73.9 +- 105.7; per-event +19.5 +- 75.4 (null)", evidence=f"{K}/out/t8_s6_and_clipping.json S98.own11x11_leave_one_out ; {K}/out/t8_t9_decomposition.json ; {K}/out/t8_a_events_and_nulls.json", method="leave-one-out and event-restricted T9 replay; random-phase null for the per-event statistic")],
      trace="T9 used QUALITY==0 11x11 cut-outs (15 S98 events incl. epoch 8 and the scattered-light epochs 2,3,11,12). P s03 used the pubprep a2 bitmask plus QLP 2^29 and the gap/dump rules (8 events). The deep follow-up numbers are single-fold aperture depths in the S32 2-min TPF with no errors quoted (deep_followup_report.md L116-L117).",
      diagnosis="flawed statistic: one outlier event (epoch 8, next to a 5.4-h data gap with a momentum dump) produces the whole T9 excess; the deep S32 trend is a ~1 sigma, non-monotonic difference of a different quantity",
      rerun=f"{K}/t8_s6_and_clipping.py -> {K}/out/t8_s6_and_clipping.json (LOO epoch 8: +4.0 +- 70.7); {K}/t8_t9_decomposition.py -> {K}/out/t8_t9_decomposition.json",
      resolution="No significant Signal A aperture dependence in S98 or S32. The T9 3.77 sigma value is an artefact of S98 epoch 8. S32 TPF r1-r2 = -62 +- 60 ppm (the deep follow-up's 374 vs 452 ppm is the same ~1 sigma effect, and r3 < r2, so it is not an off-target pattern). S6 cannot test it.",
      reason="A clean rerun from raw pixels with leave-one-out outranks the single aggregated T9 number. The excess vanishes when one event with a documented gap and dump is removed, and it is absent in T9's own method on the clean events.",
      status="RESOLVED"),
 dict(id="X-T8-02", quantity="Cause of the ~180 ppm S98 difference between the T9 method and the per-event method", master_ids="",
      claims=[dict(source="P", value="'In S98 the two methods differ by about 180 ppm, about 2 sigma, and the cause is not identified'", evidence=f"{P}/SUMMARY.md Conclusions 1", method="comparison of s03b clean-selection T9 value with the s03 per-event value"),
              dict(source="P", value="REVIEW_2: 'two methods differ by about 2 sigma (178 ppm) ... candidates: in-transit window, detrend'", evidence=f"{P}/REVIEW_2.md section 4", method="reviewer scratch"),
              dict(source="other", value="user summary (d): 'S98 retains a method-dependent ~200-230 ppm difference whose significance depends on the error model'", evidence="task statement", method="summary"),
              dict(source="other", value="T8: event selection, not method. Same 8 events: T9 -73.9 vs per-event +19.5, paired -93 +- 63 (1.5 sigma). A full-T14 window in the per-event method gives -1.8 +- 50.3.", evidence=f"{K}/out/t8_b_and_paired.json S98_paired_T9method_minus_perevent_same8 ; {K}/out/t8_a_events_and_nulls.json S98_ffi.fullT14_window", method="event-restricted T9 replay + paired differences")],
      trace="s03b clean selection keeps all clean cadences, i.e. 10 A events including epoch 8 (+3822 ppm x 29 points) and epoch -1 (-761 ppm). s03 removes both by the gap rule.",
      diagnosis="definition difference: the two numbers are computed over different event sets",
      rerun=f"{K}/t8_t9_decomposition.py ; {K}/t8_b_and_paired.py",
      resolution="The methods agree within 1.5 sigma on identical events. No method-dependent S98 excess remains; the 180 ppm is epoch 8.",
      reason="Direct test on identical events. The candidate causes (window, background, ephemeris) were tested and each moves the value by <10-50 ppm.", status="RESOLVED"),
 dict(id="X-T8-03", quantity="Wording 'T9's own method still finds about 200 ppm on the same clean data'", master_ids="",
      claims=[dict(source="P", value="SUMMARY Conclusions 1; REVIEW_1 section C ('T9's own method on the same clean cadences still gives +207 +- 72')", evidence=f"{P}/SUMMARY.md ; {P}/REVIEW_1.md L85-L87", method="T9 replay on clean cadences"),
              dict(source="other", value="same cadences but 10 vs 8 events; on the same events -74 ppm", evidence=f"{K}/out/t8_t9_decomposition.json S98_clean_T9eph", method="event-restricted replay")],
      trace="as X-T8-02", diagnosis="overstatement", rerun=f"{K}/t8_t9_decomposition.py",
      resolution="Replace with: 'T9's method gives ~200 ppm on the clean cadences only because they still contain epoch 8; on the same 8 events it gives -74 +- 106 ppm.'",
      reason="rerun", status="RESOLVED"),
 dict(id="X-T8-04", quantity="S6 T9-method value and its error", master_ids="",
      claims=[dict(source="verif", value="149.5 +- 65.5 ppm (2.3 sigma)", evidence=f"{V}/out/t9_aperture_diff.json S06.A_r1.0_minus_r2.0", method="T9, one rng stream across S98 then S06"),
              dict(source="P", value="149.5 +- 64.8 (11x11 replay); 232 +- 52 (21x21); 229 +- 50 clean", evidence=f"{P}/out/t9_method_replay.json", method="s03b re-seeded per sector"),
              dict(source="other", value="clipped-annulus 21x21 = 149.5 +- 64.8 exactly; LOO epoch -794: 25.5 +- 64.3", evidence=f"{K}/out/t8_s6_and_clipping.json S06", method="rerun")],
      trace="Same pixels (the 11x11 cut-out is the centre of the 21x21). In the 11x11 stamp, T9's 4-5.5 px annulus is clipped by the edge (47 of 52 pixels in S6).",
      diagnosis="rounding/rng seeding (65.5 vs 64.8); definition difference (annulus pixels) for 149.5 vs 232; the value is dominated by one dump-adjacent event",
      rerun=f"{K}/t8_s6_and_clipping.py", resolution="S6 T9-method 149.5-232 ppm depending on the background annulus; not a usable test (all transits 0.4-3.5 h from dumps).",
      reason="exact reproduction of both numbers from one cube", status="RESOLVED"),
 dict(id="X-T8-05", quantity="B S32 r2 minus r3 aperture difference", master_ids="",
      claims=[dict(source="P", value="-82.0 +- 20.6 ppm (nominal 4.0 sigma) in JSON; not quoted in SUMMARY", evidence=f"{P}/out/after_decorrelated.json B.S32_tpf2min.summary.corrected.diff_2.0-3.0", method="error = SD of N=2 event values / sqrt(2)"),
              dict(source="other", value="-80.9 +- 58.7 ppm (1.4 sigma)", evidence=f"{K}/out/t8_b_and_paired.json B_diff_2.0-3.0_combined_white", method="white error of each event's difference light curve, inverse-variance")],
      trace="same fits; only the error differs", diagnosis="flawed error estimate (N=2 scatter)", rerun=f"{K}/t8_b_and_paired.py",
      resolution="-81 +- 59 ppm; not significant; on-target PRF predicts +15 ppm, a 43'' source predicts -728 ppm", reason="an error from N=2 scatter has 1 degree of freedom", status="RESOLVED"),
 dict(id="X-T8-06", quantity="B in-transit centroid shifts quoted in SUMMARY", master_ids="",
      claims=[dict(source="P", value="SUMMARY: x -0.8 +- 0.8, y +1.4 +- 2.6 mpix (-164); x -0.1 +- 0.6, y -0.2 +- 0.8 (-162)", evidence=f"{P}/SUMMARY.md Conclusions 2", method="text"),
              dict(source="P", value="v3 JSON: -164 x -0.56 +- 0.88, y +0.35 +- 2.80; -162 x -0.08 +- 0.57, y -0.30 +- 0.79", evidence=f"{P}/out/b_s32_events.json events.*.centroid_shift_mpix", method="s04 with A masked")],
      trace="SUMMARY numbers predate the v3 A-masking in s04", diagnosis="stale data in SUMMARY text", rerun="none needed (JSON reproduced byte-identically)",
      resolution="use the JSON values; the conclusion (centroids not decisive, a 43.5'' source ~1 mpix) is unchanged", reason="the JSON is produced by the current script", status="COSMETIC"),
 dict(id="X-T8-07", quantity="REVIEW_1/REVIEW_2 numbers vs v3", master_ids="",
      claims=[dict(source="P", value="REVIEW_1: S98 15 events +105 +- 180, 9 clean +84 +- 60, T9 clean +207 +- 72, S6 5 events; REVIEW_2: S98 +230 +- 211 (n=10), combined -15 +- 53, 'without epoch 8 +22 +- 43 (N=9)', C4 'not fixed'", evidence=f"{P}/REVIEW_1.md L78-L87 ; {P}/REVIEW_2.md sections 2-3", method="v1/v2 reviews"),
              dict(source="P", value="v3: S98 +19.5 +- 48.6 (n=8); combined -12.4 +- 37.9; T9 clean 199.8 +- 73.5; S6 0 events (4 without dump rule, +102 +- 127); s04 masks A", evidence=f"{P}/out/after_decorrelated.json ; {P}/out/t9_method_replay.json", method="v3 scripts")],
      trace="the reviews ran on v1/v2", diagnosis="stale data (earlier versions)",
      rerun=f"{K}/t8_s98_rule_sensitivity.py reproduces REVIEW_2's +22.3 +- 42.9 (N=9, epoch -1 kept) and v2 +229.9 +- 211.1 (N=10)",
      resolution="Still valid in v3: S32 TPF -62 +- 60, S32 FFI -36 +- 65, S6 no-dump-rule +102 +- 127, all T9 replay numbers, the annulus-clipping explanation, and the T9-on-S32 scratch numbers. Superseded: S98 per-event values, combined values, C3/C4 status, 'cause not identified'.",
      reason="v3 outputs reproduce byte-identically", status="RESOLVED"),
 dict(id="X-T8-08", quantity="Manuscript used for P's impact assessment", master_ids="",
      claims=[dict(source="P", value="assessed pubprep manuscript.tex: A (+8+-5, -8+-7)'', B (+1+-3, +1+-4)''", evidence=f"{P}/SUMMARY.md Manuscript impact ; {P}/REVIEW_1.md section E", method="text review"),
              dict(source="U", value="proposed revision: A (+5+-3, 0+-3)'', B (+1+-3, -1+-3)''", evidence=f"{MS} L37, L41", method="U 03_pixels")],
      trace="P never opened U", diagnosis="stale reference document", rerun="impact re-assessed on U's revision (PIX.manuscript_impact)",
      resolution="No change needed in U's revision either", reason="the Note makes no aperture claim, and P does not constrain positions at the 3'' level", status="RESOLVED"),
 dict(id="X-T8-09", quantity="Ephemeris used by P vs U final", master_ids="",
      claims=[dict(source="P", value="audit ephemeris P=3.182785293657263 d, T0=3998.7367930 BTJD", evidence=f"{P}/out/inputs.json ephemeris_BTJD.A", method="audit fit"),
              dict(source="U", value="P=3.182785+-0.000006 d, T0=2460998.736+-0.003 BJD_TDB", evidence=f"{MS} L33", method="U timing")],
      trace="the T0 values agree to 0.0003 d; P uses +-T14/4 boxes", diagnosis="negligible difference",
      rerun=f"{K}/out/t8_t9_decomposition.json S98_clean_auditEph vs S98_clean_T9eph: 206.9 vs 199.8 ppm (7 ppm)",
      resolution="no effect on any conclusion", reason="rerun", status="NOT_A_CONFLICT"),
 dict(id="X-T8-10", quantity="User summary (b): A central/aperture effect 'NOT confirmed and NOT excluded'", master_ids="C092",
      claims=[dict(source="other", value="not confirmed, not excluded", evidence="task statement", method="summary"),
              dict(source="other", value="not confirmed; the only >3 sigma evidence is epoch 8; a ~200 ppm excess is disfavoured at ~3.8 sigma by S98+S32 (-16 +- 57 ppm); effects <~100 ppm are not excluded; S6 is untestable", evidence=f"{K}/out/t8_s6_and_clipping.json ; {K}/out/t8_a_events_and_nulls.json", method="reruns")],
      trace="see X-T8-01", diagnosis="the summary predates the event-level decomposition", rerun="as X-T8-01",
      resolution="'Not confirmed. The T9 3.8 sigma value is caused by one event; the T9-sized excess is disfavoured, and smaller effects are not tested.'", reason="rerun evidence", status="RESOLVED"),
 dict(id="X-T8-11", quantity="User summary (a), (c), (e)", master_ids="",
      claims=[dict(source="other", value="(a) B/S32 ~4 sigma off-target = methodological false positive; (c) S32 no significant aperture dependence; (e) do not promote to the Note", evidence="task statement", method="summary"),
              dict(source="other", value="(a) confirmed (first pass 4.5-5.0 sigma; per-event flat, also without centroid terms); (c) confirmed under every method (<1.1 sigma); (e) confirmed, no claim changes", evidence=f"{K}/out/t8_b_and_paired.json ; {K}/out/t8_t9_decomposition.json ; {K}/out/t8_a_events_and_nulls.json", method="reruns")],
      trace="P outputs + T8 reruns", diagnosis="no conflict", rerun="as cited", resolution="confirmed as stated (the first-pass significance is 4.5-5.0 sigma, not ~4)", reason="", status="NOT_A_CONFLICT"),
 dict(id="X-T8-12", quantity="inputs.json byte difference on rerun", master_ids="",
      claims=[dict(source="P", value="interpreter '.../tic4206066_pixel_recheck/scripts/../../../.venv/bin/python'", evidence=f"{P}/out/inputs.json software.interpreter", method="s01"),
              dict(source="other", value="interpreter '/Users/paulsmith/space <3/tess-jev/.venv/bin/python'", evidence=f"{K}/pkg_copy/out/inputs.json", method="s01 in the copy")],
      trace="the same interpreter called by a different path", diagnosis="cosmetic", rerun="none needed", resolution="all science outputs are identical", reason="", status="COSMETIC"),
]

mf = [
 dict(quote="place the 3.18\\,d source at $(+5\\pm3, 0\\pm3)''$ from the target ($\\approx2\\sigma$)", problem="None from the pixel re-check. U's S98 set includes S98 epoch 8, the event with the pixel-level anomaly, but U's own jackknife without it moves S98 only to (+7.0, +0.3)'' (diagnose_s98.json). P does not constrain positions at the 3'' level.", fix="no change", severity="NONE"),
 dict(quote="All Gaia DR3 stars able to produce the depth lie $\\ge43''$ ($\\ge14\\sigma$) away.", problem="None. The PRF shows that a source at 43'' would make the r<=1 px depth shallower by 268-364 ppm. The measured per-event values (+20 +- 75 ppm in S98, -62 +- 86 in S32) independently disfavour it.", fix="no change", severity="NONE"),
 dict(quote="offset $(+1\\pm3, -1\\pm3)''$", problem="None. B's per-event S32 depths are flat across apertures (r1/r3 ~ 1.0, against 0.12 for a 43'' source). The first-pass off-target pattern was a detrending artefact.", fix="no change", severity="NONE"),
 dict(quote="Its Sector 6 difference images are unusable because every transit falls 0.4--3.5\\,h from a reaction-wheel desaturation.", problem="None. Rerun dump distances are 0.43-3.53 h for the fittable S6 A events.", fix="no change", severity="NONE"),
 dict(quote="Vetting, pixel-level localization, Gaia astrometry and photometry, radial velocities, and overluminosity", problem="None. The pixel re-check leaves no open aperture anomaly that undermines 'pixel-level localization'.", fix="no change", severity="NONE"),
]

unknowns = [
 dict(item="Physical cause of the S98 epoch-8 pixel anomaly (+3800 ppm r1-r2 in the T9 method, +2098 ppm per-event)", why="It coincides with a 5.4-h data gap that contains a momentum dump 7.9 h after mid-transit. The one-sided detrend is a plausible mechanism, but it was not isolated.",
      would_resolve="Inspect the epoch-8 pixel time series and the pointing/background around the gap. Not needed for the Note.", blocks_submission=False),
 dict(item="Residual S6 T9-method aperture excess (~100 ppm without epoch -794; 150-230 ppm with it)", why="All S6 transits lie 0.4-3.5 h from momentum dumps, and the per-event method has no usable S6 event.",
      would_resolve="Future sectors without dumps near transits (e.g. S110) or a full pixel-level decorrelation of S6", blocks_submission=False),
 dict(item="B aperture behaviour rests on 2 S32 events with white (lower-bound) errors", why="Only two S32 2-min B events exist.", would_resolve="S98/S110 B events at 2-min or 200-s cadence", blocks_submission=False),
]

summary = ("P reproduces byte-identically; only inputs.json differs, in the interpreter path. The v3 fixes are correct. "
 "(a) Confirmed: the B/S32 first-pass off-target pattern (r1-r3 -739 +- 159, 4.6 sigma) comes from global detrending and a constant baseline during a ~38 mpix drift. Per-event fits, with or without centroid terms, give flat depths (924/949, 1061/999 ppm). "
 "(b) Stronger than stated: T9's 3.77 sigma S98 excess comes entirely from one event, epoch 8, next to a 5.4-h gap containing a momentum dump. Without it the value is +4 +- 71 ppm. On the same 8 clean events T9's method gives -74 and the per-event method +20 +- 49/75. A ~200 ppm excess is disfavoured at ~3.8 sigma by S98+S32; effects below ~100 ppm are untested. "
 "(c) Confirmed: S32 is <1.1 sigma under every method (T9-method S32 numbers reproduced exactly). "
 "(d) Refuted as a method effect: the ~180 ppm S98 difference is event selection (epoch 8). The paired difference on identical events is -93 +- 63. The error model matters only for the contaminated statistic: white 5.0, null 3.8, event scatter 0.8 sigma. "
 "S6 T9-method values (149.5/232) follow from annulus clipping and one dump-adjacent event; S6 is untestable. "
 "The PRF confirms that a 43'' source makes the core shallower (r1-r2 -268 to -364 ppm), so 'core deeper' is not an off-target signature. "
 "(e) Confirmed: no claim in U's manuscript changes. M C092 is resolved.")

files = []
for dp, dn, fn in os.walk(W):
    dn[:] = [d for d in dn if d != ".cache"]
    for f in fn:
        files.append(os.path.relpath(os.path.join(dp, f), "/Users/paulsmith/space <3/tess-jev"))
files.append(f"{K}/result.json")
files = sorted(set(files))
res = dict(topic="T8_pixel_recheck", summary=summary, facts=facts, conflicts=conflicts, manuscript_findings=mf, unknowns=unknowns, files_written=files)
json.dump(res, open(os.path.join(W, "result.json"), "w"), indent=1, ensure_ascii=False)
print(len(summary.split()), "words;", len(facts), "facts;", len(conflicts), "conflicts;", len(files), "files")
