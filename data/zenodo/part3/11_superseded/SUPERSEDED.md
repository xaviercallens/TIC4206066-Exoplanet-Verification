# Superseded results (relevant to numbers printed in the Note)

Earlier stages of this project produced values that were later corrected or replaced. This table lists the superseded values that concern a number printed in the Note. **Do not quote the old values.**

Columns:
- *Old value*: the superseded value.
- *Where it appeared*: the earlier stage that produced it.
- *Replaced by*: the final value, with its `FINAL_FACTS.csv` fact id in brackets.
- *Why*: the reason for the change.

Stage names:
- **pubprep**: the author's earlier publication-prep package (`../13_supporting_pipelines/pubprep`).
- **audit**: the independent publication audit (`../13_supporting_pipelines/PUBLICATION_AUDIT`).
- **U**: the opus55 final audit (all `from_U_*` folders).
- **M**: the master-evidence compilation.
- **P**: the pixel re-check.
- **pool**: the multi-sector pool pass of the wider search.

Source: `FINAL_TRUTH.md` §12 (final reconciliation, 2026-09-25), which merges `CONFLICT_RESOLUTION.md` and `tic4206066_MASTER_EVIDENCE/SUPERSEDED_RESULTS.md`; neither file is included here. 20 rows that concern wording, metadata or values the Note does not print are left out; they are listed at the end.

| Quantity | Old value | Where it appeared | Replaced by [fact id] | Why |
|---|---|---|---|---|
| P_A | 3.1827856 ± 0.0000058 d | pubprep (20-set) | 3.1827850 ± 0.0000057 [A.P] | S98 E-1 and S32 E-573 lack post-egress baseline; 19-set |
| P_A | 3.1827853 ± 0.0000027 | audit joint fit | [A.P] | formal errors ~1.7x too small [A.P_jointfit_formal] |
| P_A | 3.18281 d | pool SUMMARY | not used [DATA.pool_P_A] | stale pass-1 value (wrong phase) |
| T0_A | 2460998.7375 ± 0.0023 | pubprep | 2460998.73599 ± 0.00288 [A.T0] | entirely due to S98 E-1 [A.T0_shift_cause] |
| N timing | 20 / 21 transits | pubprep, U 01, audit | 19 [A.N_timing] | event selection |
| Timing chi2 | 12.6/18; 47.5/19 | pubprep; audit | 9.5/17 [A.timing_chi2] | set; red-noise error model too small [A.timing_error_model] |
| Timing chi2 (printed) | "χ²/dof ≈ 9/17" | submission candidate before pass 2 | "≈ 9.5/17" [A.timing_chi2] | 9.51 does not round to 9 (R2P2-02) |
| Rp/R\* | 0.0212 ± 0.0007; 0.021 ± 0.001 | pubprep/audit; U revision | 0.0204 ± 0.0010 [A.k] | centre of the reduction range |
| Rp_A | 1.44 ± 0.06 R⊕ | pubprep / M | 1.39 ± 0.08, print 1.4 ± 0.1 [A.Rp] | follows k |
| Depth_A | 551 ± 34 ppm | pubprep | 498 ± 29 [A.depth_combined] | epoch -1 and free-dt bias |
| Depth S98 | 574 ± 42 ppm | pubprep | 496 ± 48 [A.depth_S98] | epoch -1 excluded |
| Depth range | 483-558 ppm | audit | 456.9-558.2 [A.depth_six_detrendings] | wrong low end |
| rho_free | "< 3.5 for both priors" | audit | ≲3.6 (3.4-3.7) [A.rho_free] | stored uniform value 3.55 |
| Phase 0.5 | "flux −20 ± 32 ppm" | pubprep | depth −20.4 ± 31.7 [A.phase05] | sign |
| S6 effect | "<0.5σ" | U revision | ≲1σ (max 1.01) [A.S6_effect] | like-for-like refits |
| A offset | (+8 ± 5, −8 ± 7)″ | pubprep | (+5 ± 3, 0 ± 3)″ [A.offset] | raw S32 uninformative; decorrelated S32 |
| A offset | (+2.51, −7.53) ± 5.7″ | audit | [A.offset] | same |
| Capable stars | ≥43″ ≥7σ / ≥8σ / ≥14σ | pubprep / audit / U | ≳43″ (>10σ) [A.capable_stars, A.capable_sigma] | stale positions; 42.78″ at S98 |
| Nearest source | "within 11″" | U revision | "within 10″" [A.nearest_gaia_neighbour] | 10.35″ at the S98 epoch |
| S98 A epochs | 10 | pubprep | 9 [A.S98_epochs] | event without post-transit baseline |
| B offset | (+1 ± 3, +1 ± 4)″ | pubprep | (+1 ± 3, −1 ± 3)″ [B.offset] | transcription mix |
| T9 core excess | 226.6 ppm "unexplained" | verif / M C092 | not confirmed [PIX.A_aperture_status] | one event (epoch 8) |
| B/S32 off-target | "~4σ" | author summary of P | 4.5-5.0σ methodological false positive [PIX.B_S32_firstpass, PIX.B_S32_offtarget] | global detrend + drift |
| B r2−r3 | ± 21 ppm (nominal 4σ) | P JSON | −81 ± 59 ppm [PIX.B_S32_r2_minus_r3] | scatter of 2 values |
| FPP_A | 0.031 ± 0.003 (0.031-0.034) | pubprep / M | 0.03-0.04 [FPP.A_default] | one fold; mislabelled variant |
| FPP_A constrained | MOLUSC 0.005 ± 0.003 | pubprep / audit / U | 0.002-0.007 [FPP.A_constrained] | scales as 1/N [FPP.A_molusc_mode] |
| FPP share | "~95%" | U08 | about 90% (84-96%) [FPP.A_share_above_0.38] | A_base only |
| NFPP | ≤1e-7 in every variant | U08 | <1e-6 [FPP.A_NFPP] | audit folds reach 4.6e-7 |
| RV limit | K < 3.2 km/s, M sin i < 0.016 / 0.024 / 0.025 | pubprep / M / audit | K < 1.16 km/s, < 0.0058 / 0.0088 Msun [COMP.Gaia_RV] | 3.2 is peak-to-peak noise |
| Circular host | ≥0.54 Msun | pubprep / audit / M | ≈0.52 [COMP.host_mass_circular] | table interpolation |
| Eccentricity | e ≳ 0.6 near apoastron | M / audit | e ≳ 0.2-0.6 [COMP.host_lower_bound_and_ecc] | single value, lowest host only |
| Overluminosity | +0.003 ± 0.07 (Way et al.); "≲20 au ~ 0.6″" | pubprep; U06 | author scatter; ≲0.2″ [COMP.XP_overluminosity] | not in Way; angular limit |
| Parallax | 28.108 ± 0.019 | M | 28.1083 ± 0.0176 [STAR.parallax] | error not reproduced |
| P_B | 11.132976; 11.132339 | S98 grid; pool | 11.13274 ± 0.00003 [B.P] | grid peak; pool not locked |
| P_B error | ± 1.3e-5 | audit joint | ± 0.00003 [B.P] | white-noise emcee |
| B significance | "~3σ overall"; "0 of 10,000"; epoch-shifted nulls | pubprep / audit / deep | FAP ≈2e-3, post hoc [B.significance] | flawed / edge-inflated nulls |
| B S6-only | "0 of 10,000", "0 of 18,652", "<~7e-4" | pubprep / audit / U auditor | ≈2e-3 [B.significance_S6only] | excluded region in null; floor ≠ ceiling |
| B S98-only | 0.35 | pubprep | 2e-4 to ~1 [B.S98_alone] | one flux column + statistic |
| B S98 events | "three in clean intervals" | M / pubprep | two clean [B.S98_transits] | E4 half-flagged |
| Rp_B | "≈2 R⊕" (headline) | pubprep | ≈2.2 if on target [B.Rp_on_target] | free-density value |
| S110 count | "eight or nine ... and three" (no gap) | U revision | with "one of them in a scheduled data gap" [FUTURE.S110_counts] | official orbit table |
| Shortlist | "preregistered" | pubprep / U C5 | post hoc exploratory cohort [DATA.search.shortlist] | amendment 01 |

## Package-specific notes

| Item | Old | Where | Replaced by | Why |
|---|---|---|---|---|
| Signal B significance test | anchored ±0.004 d scan and null written 2026-09-25 03:26-03:29 UTC | `../07_signal_B/from_deep_followup_original_posthoc_test/` | FAP ≈ 2×10⁻³ (pubprep null with B windows removed) [B.significance] | kept only as chronology evidence; its null was superseded (edge-inflated / flawed nulls) |
| B, Sector 32 centroid numbers and the "~4σ off-target" first pass | as written in `SUMMARY.md` of the pixel re-check | `../06_pixels_localisation/from_tic4206066_pixel_recheck/SUMMARY.md` | methodological false positive; use the JSON outputs and the T8 verification [PIX.B_S32_offtarget, PIX.B_S32_firstpass] | SUMMARY's centroid numbers predate its v3 fixes |
| Radius / Teff quoted in the companion check | R 0.621 ± 0.015 Rsun, Teff 4012 ± 43 K, planet 1.44 R⊕ | `../05_companion_limit_038/from_companion_check/SUMMARY.md` | R 0.62 ± 0.02, Teff 4020 ± 60 K, Rp 1.4 ± 0.1 R⊕ [STAR.R_ARIADNE, STAR.Teff, A.Rp] | `.dat` errors are artefacts; BMA mixture adopted; k revised |
| Archival FAP of the 3.18 d signal (printable sentence of the first run) | "<0.7% in each of S6 and S32, <0.03% jointly" (N_eff from a heuristic count: 164 / 195 / 3658) | `../03_archival_test_A/a_fap_results.json` → `printable`, `OPTIMISED.*.fap_conservative` | "≲1% per sector, ≲0.1% jointly" [A.archival_FAP] | pass-3 reviewer measured N_eff directly (`../03_archival_test_A/reviewer_neff/`, `PASS3_NOTE.json`): S6 bound 1.0-1.3%, joint 0.006-0.06% |
| MOLUSC-mode FPP | 0.005 ± 0.003 | pubprep / audit / U runs `*molusc*` in `../04_fpp` and `../13_supporting_pipelines` | not a constraint; post hoc cut 0.002-0.007 [FPP.A_constrained, FPP.A_molusc_mode] | scales as 1/N |

Many copied REPORT.md and result files record the values that were current at their stage. Where they disagree with `../12_facts/FINAL_FACTS.csv`, FINAL_FACTS is authoritative.

## Rows left out (not numbers printed in the Note)

The ARIADNE rows (Teff, R_ARIADNE) moved here in the pass-3 revision: that clause is no longer printed.

- TTV
- Rp_A
- Depth consistency
- FPP_B
- Teff
- R_ARIADNE
- L
- Rotation
- Spectral type
- B depth diff
- B odd/even
- Teq_B
- P_B/P_A
- S110 A times
- S110 B times
- ExoFOP counts
- QLP citation
- ChatGPT models
- AI statement status
- Word count
