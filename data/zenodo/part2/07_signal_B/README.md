# 07_signal_B: the 11.13 d signal, and Sector 110 predictions

## Summary tables (start here)

| file | content |
|---|---|
| `ephemeris_B.json` | Printed ephemeris P = 11.13274 ± 0.00003 d, T0 = 2461000.865 ± 0.003 BJD_TDB, valid only if the archival dips are real. Two fits are included. The 4-event fit (S6 E-227, E-226; S98 E0, E1) gives P = 11.1327392 ± 0.0000149 and T0 = BTJD 4000.86552 ± 0.00158, χ² 1.45 for 2 dof. The 7-event fit adds S32 and S98 E4. The printed errors are conservative, about twice the formal ones. |
| `B_events_at_P_harmonic_audit.csv` | Box S/N at every predicted P-epoch in S6, S32 and S98, as written by the frozen harmonic audit at 03:08:22 UTC, before the recovery test was defined. S6 dips present (S/N 4.5, 6.1); S32 events "contradicting" (S/N 2.3, 2.6; they abut 3.18 d transits); S98 E0, E1 clean; E4 partial; E-1, E2, E3 masked. |
| `CHRONOLOGY_file_times.csv` | Creation and modification times (UTC), re-read read-only from the author's file system when this package was built, for the files that fix the post hoc order. |
| `S110_predictions.csv` | Predicted TESS Sector 110 transit times for both signals: BJD_TDB, UTC as seen by TESS (light-time corrected), 1σ in minutes, ingress/egress, and status against the published orbit table. Source: `from_T7_data_future/s110_predict.json`. |
| `FACTS_IN_THIS_SECTION.csv` | All `B.*` and `FUTURE.*` facts, with evidence mapped to package paths. |

## Significance: the archival recovery test (fact `B.significance`)

- **Test.** Link the Sector 98 ephemeris to Sectors 6+32 and take the maximum box S/N over P within ±0.004 d of the S98 period.
- **Null.** The same scan at a random period in 9-13.5 d (excluding ±0.1 d around P_B) and a random epoch, on the S6+S32 data with B's own archival transits removed. Signal A is masked. FAP = (n ≥ obs + 1)/(N + 1).
- **Result.** 19 of 10,000 null trials reach the observed S/N 10.03, so FAP = 2 × 10^-3, i.e. 2.88σ one-sided. The Note prints "≈3σ (false-alarm probability ≈2×10^-3)".
  - Script: `../13_supporting_pipelines/pubprep/analysis/a3_signalB_tests.py`.
  - Output: `../13_supporting_pipelines/pubprep/analysis/out/a3_signalB.json`, key `T1_link[0]`.
- **Other implementations:**
  - 1.25 × 10^-3: audit, `../13_supporting_pipelines/PUBLICATION_AUDIT/out/s11b_nullB_periods.json`;
  - 0 of 50,000: U, `from_U_04_B_existence/out/b04b_link_extra.json`, with the verifier in `verify/out/v02_link.json`;
  - σ conversions: `from_T5_signal_B/t5_rerun.json` → `sigma_one_sided`.
- **Single sectors.** Sector 6 alone ≈ 2 × 10^-3 (rank test; `a3_signalB.json` `T2_rank[0]`, `from_U_04_B_existence/out/b05_s6_rank.json`). Sector 32 alone 0.02-0.05 (`T1_link[2]`), hence "marginally in Sector 32".
- **Sector 98 alone** does not establish the signal. Its FAP ranges from 2 × 10^-4 to ~1 depending on the flux column and statistic (`from_U_04_B_existence/out/b03_null_*.json`).
- **Trials.** The FAP is per candidate. It does not include the trials of the 126,246-star search.

## Post hoc chronology (fact `B.post_hoc`)

All times are 2026-09-25 UTC (`CHRONOLOGY_file_times.csv`; also `from_T5_signal_B/t5_rerun.json` → `chronology_utc` and `from_T5_signal_B/verify/v_chron_a03.json`):

1. **02:54:49** The unanchored multi-sector search did not lock on B (`addenda/pool_multisector/multisector_results.json`).
2. **03:07:17** The pool summary recorded "not confirmed".
3. **03:08:22** The frozen harmonic audit wrote the S98-anchored S6 dips to file (`data/harmonic/per_candidate/tic4206066_p11.132976.json`).
4. **03:26:34** Only then was the anchored ±0.004 d scan created (`from_deep_followup_original_posthoc_test/signal_11d.py`).
5. **03:29:42** Its null test was written (`signal_11d_null.py`, `signal_11d_null.json`).

The recovery test was therefore defined after the unanchored search failed and after the S6 dips were on file. Hence "a recovery test defined after inspecting the archival data, where an unanchored search had failed".

The pool files and the harmonic-audit files belong to the author's wider TESS search and are **not** included; only their times are listed.

## Depths, archival events, Sector 110

- **Per-sector depths** depend on reduction choices (facts `B.depth_S6`, `B.depth_S32`, `B.depth_S98`):
  - S6: 830-920 ppm without the single SPOC bit-512 cadence;
  - S32: about 600-1000 ppm;
  - S98: 1100-1180 ppm with QLP products.

  Common depth ≈1000 ppm (`from_U_05_B_depth/out/s04_grid.csv`, `s07_summary.json`; `verify/v02_depths.json`). Rp ≈ 2.2 R⊕ if on the target (`../13_supporting_pipelines/PUBLICATION_AUDIT/FINAL_PARAMETERS.json` → `signal_B.shape_mann_prior`; `from_U_05_B_depth/out/s02_joint_prior_*.json`).
- **Archival events** (`B.archival_recovery`):
  - S6 E-227 is clean.
  - S6 E-226 rests on about 2 unflagged 30-min cadences.
  - S32 E-164 and E-162 start 0.16-0.5 h after 3.18 d transits end (`from_T5_signal_B/t5_rerun.json` → `S32_B_vs_A`; `from_U_05_B_depth/verify/v01_info.json`).
- **Sector 98 events** (`B.S98_transits`): two clean transits (E0, E1). E4 has about 50% effective coverage. The others fall in QLP-flagged intervals (`from_U_04_B_existence/out/b01_s98_events.json`, `logs/b01.log`).
- **Sector 110** (`FUTURE.*`): HEASARC gives 2026 Oct 31-Nov 26. The MIT orbit table (`from_T7_data_future/web/TESS_orbit_times.csv`, with HTTP headers and fetch log) was used to predict transit times:
  - 3.18 d: 8 mid-times from the sector start to Nov 26 00:00, 9 if data continue past Nov 26 ~14:48 UTC; one (Nov 13 ~21:15) falls in the orbit 231/232 downlink gap.
  - 11.13 d: 3 mid-times.

  Script `from_T7_data_future/s110_predict.py`, with an independent check in `verify/s110_indep.py`. The Sector 98 cadence and quality-bit counts are in `from_T7_data_future/s98_quality_bits.json`.

## Source folders

| folder | content |
|---|---|
| `from_U_04_B_existence/` | Sector 98 events and search (`b01`–`b03`), link test and nulls (`b04`, `b04b`), S6 rank (`b05`), ephemeris (`b08`), diagnostics; independent verifier (`verify/`). |
| `from_U_05_B_depth/` | Depth grid over reductions (`s04_grid`), joint fits (`s02`, chains not included), injections (`s03`), flux-scale check (`s05`), Claret LD (`s00`); verifier. |
| `from_T5_signal_B/` | Reconciliation rerun: ephemerides, σ conversions, S110 rows, S32 A/B gaps, chronology; light-time check (`t5_ltt`); verifier. |
| `from_T7_data_future/` | S110 predictions, S98 quality bits, QLP headers, orbit table (`web/`), TASOC check; verifier. |
| `from_deep_followup_original_posthoc_test/` | The original anchored scan and null (`signal_11d*.py`/`.json`), kept as the chronology evidence. **Superseded** by the pubprep/U tests for the printed FAP. |

**Not included.** MCMC chains (`s02_chain_*.npz`) and own-photometry arrays (`ownphot_*.npz`); see `../12_facts/EXCLUDED_FILES.csv`.
