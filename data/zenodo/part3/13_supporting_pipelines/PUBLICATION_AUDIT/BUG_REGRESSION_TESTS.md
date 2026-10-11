# Bugs found and regression tests — TIC 4206066 publication audit

Scope: every script that produced a number meant for the manuscript (deep follow-up, companion check, FPP, RNAAS draft),
plus this audit's own code. Every bug has an executable test. Run both files with the audit environment:

```
audit_env/bin/python tests/test_regressions.py        # 15 tests  -> tests/test_regressions_result.json
audit_env/bin/python tests/test_prf_regressions.py    #  4 tests  -> tests/test_prf_regressions_result.json
```

Status at 2026-09-25: **19 / 19 pass.**

## A. Bugs and method errors in the earlier TIC 4206066 work

| # | Where | What was wrong | Effect on published-candidate numbers | Test |
|---|---|---|---|---|
| A1 | `deep_followup_tic4206066/lc_analysis.py` (`detrend`, `windows`, phase scan, odd/even) | Signal B was not masked when signal A was detrended, fitted, phase-scanned and split into odd/even. Two S32 B transits fall 2.0–2.4 h after A transits, inside A's ±0.35 d fit windows. | Small. With B masked, A gives P = 3.1827853 ± 0.0000027 d (was 3.1827863 ± 0.0000028; 0.4σ) and Rp = 1.44 ± 0.06 R⊕ (unchanged). The quoted phase-scan maximum (+89 ppm at A-phase 0.17) was B's S98 transit at 4012.00 BTJD, which falls at A-phase 0.167. With B masked the scan maximum is \|z\| = 2.2 (37 independent phases). | T1, T2 |
| A2 | `prf_shift_test.py`, `pixel_analysis.py` | Localisation errors came from a bootstrap over 6–9 transits (and white noise for the centroids). That underestimates the scatter about 3× for the S32 2-min data. | The claim "S32 2-min offset 10.5″ E at 3.8σ" is not supported. With errors from real-pixel injections the S32 offset is 0.75σ. The SPOC MOM_CENTR1 shift is −0.24 mpx with an empirical null of ±0.14 mpx (z = −1.7), not ±0.05. The pointing explanation (POS_CORR) is unnecessary. | T13 |
| A3 | `archival_checks/pixels_11d_injection.py` | Injections were made as A × (out-of-transit image) and fitted with the same image template (the TIC 158955499 "circular injection" pattern, empirical-template version). | They test noise only, not the model. Replaced here by model-free injection into real pixels plus official-PRF fits. B's pixel result stays on target. | PRF R3 |
| A4 | Pixel localisation, S6 | Not a coding error, but a missed systematic: in S6 the reaction-wheel desaturations came every 3.125 d, and **every one of the six S6 A transits lies 0.4–3.5 h from a dump** (P_A = 3.183 d). Pointing jumps enter the difference images. | The S6 A offset (7″ E, 2.4σ) cannot be calibrated: dump-matched injections are impossible because every dump carries an A transit. S6 is excluded from the headline localisation of A. | T10 |
| A5 | S98 localisation | Scattered-light intervals (QLP bit 2^29, SPOC straylight) were used in the difference images. The official S98 notes (HEASARC sector98_summary) state Earth scattered light in camera 1 at the start of each orbit. | The S98 A offset of ~15″ (0/300 injections as far) disappears in clean intervals: z = (0.1, −1.5). | T12 |
| A6 | `fpp/run_fpp.py` | TRICERATOPS `exptime` = 200 s for folds binned to 6 min. | Negligible for A (FPP 0.0304 vs 0.0312 with 6 min). Tested for B too. | T11 |
| A7 | `fpp/prepare_inputs.py` | B's fold used the S98 BLS epoch (4000.86817), not the joint-fit epoch (4000.86587, 3.3 min earlier). The calc_depths depth (783 ppm) sits below the fold's own depth. | B's FPP changes from 0.146 to 0.18 with the audit's fold (same TRICERATOPS version, same N). Test run (`s13c`): the earlier fold with the earlier configuration reproduces 0.1459 exactly (largest per-seed difference 7 × 10⁻⁵); with the audit configuration it gives 0.147. The fold alone causes the change. B's FPP is sensitive to fold details at the ~25 % level. | `s13c` (reproduction run) |
| A8 | TRICERATOPS 1.1.0 itself (documented behaviour, not a bug we can fix) | `calc_depths` uses a Gaussian PSF (σ = 0.75 px). The real TESS PRF puts 1–3.5 % of the flux of the T = 13.06 star 83.8″ away (Gaia DR3 3220385243255181184) into the 18-pixel SPOC apertures. In S32 that star could host A with a 28 % eclipse; TRICERATOPS does not list it as a possible host. | NFPP ignores that star. It is excluded by the pixel localisation (> 10σ), so the conclusion is unchanged, but NFPP ≈ 0 must not be quoted as evidence against it. | reported |
| A9 | MOLUSC b26efbf CLI | `star_jitter` is never set by the command-line parser (crash). | Same one-line patch as the deep follow-up (`out/molusc/PATCH.diff`). Survivors reproduced: 429,613 vs 429,780 of 10^6. | reported |
| A10 | `deep_followup_report.md` / RNAAS draft | TGLC S6 depth quoted as 190 ± 171 ppm (PSF flux column). | The TGLC aperture-flux column gives 982 ± 101 ppm. S6 depths depend strongly on the pipeline (521–982 ppm); only the SPOC/TESS-SPOC S6 value (521 ± 53) is used for science. | reported |
| A11 | RNAAS draft internal inconsistency | Timing χ² quoted as 120/21 (RNAAS) vs 114/21 (deep report). S98-vs-archival B depth difference quoted as 3.1σ (RNAAS) vs 2.4σ (report). | Audit values: A timing χ² = 47.5/19 with red-noise-inflated per-event errors (pipeline set); B depth difference 2.2–2.7σ across six detrendings (per-event fits with red-noise errors). | reported |
| A12 | Deep follow-up stellar mass error | Mann+2019 σ_e (2.0 %) used alone; Table 7 of Mann+2019 gives 2.4 % total at M_Ks = 5 (coefficient + intrinsic scatter). | M = 0.606 ± 0.015 M☉ (was ± 0.013). | T4 (same script) |
| A14 | `deep_followup_tic4206066` B link null test (300 trials) | The null shifted only the S98 epoch at the fixed S98 period. The ±0.004 d period scan moves the predicted archival times by ±0.9 d (S6) and ±0.66 d (S32), so trials whose epochs differ by less than that re-use the same data. The trials contain only ~5–16 distinct windows (S6), and the null maximum is just the largest feature at another phase. | "Exceeds all 300 trials" does not give FAP < 1/300. With period and epoch both randomised, the S6+S32 link has FAP 1.2 × 10⁻³. S6 alone: 0 of 18,652 trials reach it, and its two dips are the two strongest in the sector (distribution-free p = 1.6 × 10⁻³). S32 alone: FAP 0.03. | T14 |
| A13 | Teff adopted 4080 ± 100 K | Pulled up by Gaia GSP-Spec 4182 K, which carries the GSP-Spec extrapolation flag and logg 5.01 (inconsistent with M, R). | Adopted 4030 ± 90 K. Changes Teq(A) 820 → 810 K and S(A) 75 → 72 S⊕. LD coefficients change by < 0.03. | reported |

## B. Regression tests carried over from the TIC 158955499 PRF audit (tests/test_prf_regressions.py)

| # | TIC 158955499 mistake | Test here | Result |
|---|---|---|---|
| R1 | MAST PRF files read as one continuous 9×-supersampled image instead of 9×9 interleaved 13×13 PRFs | TESS_PRF sub-PRF centroids must track the requested intra-pixel position (< 0.05 px). The un-reversed interleave `prf[i::9, j::9]` must fail. | TESS_PRF max error 0.028 px; wrong reading 0.89 px → PASS |
| R2 | Host star duplicated into its own neighbour list (TIC J2000 vs Gaia 2016 positions, proper motion 70.8 mas/yr) | Scene built by Gaia source_id: exactly one target, no non-target source within 2″ | PASS |
| R3 | Circular injection (inject and fit with the same PRF model) | Injection function takes no PRF and multiplies real background-subtracted pixels exactly | PASS |
| R4 | Neighbour eclipse depth computed with full flux inclusion instead of the aperture actually used | PRF-integrated aperture requirement ≥ full-inclusion requirement for all real neighbours ≤ 60″ | PASS (e.g. the 43.5″ star needs 133 % in an r = 2 px aperture vs 45 % with full inclusion) |

Also checked (not bugs): TPF FLUX is background-subtracted, TESScut FFI FLUX is not (TESScut docs). The Stassun+2019
G→T relation degrades for BP−RP > 2; every neighbour within 60″ has BP−RP < 2 or is flagged. Gaia epoch photometry
does not exist for this star (the TIC 158955499 BTJD = t − 1802.5 conversion is not needed).

## C. Bugs found in this audit's own code and fixed before any result was used

| # | Bug | How it was caught | Test |
|---|---|---|---|
| C1 | `np.isclose` default `rtol = 1e-5` shrank the dP ranges in the fast null-test scan (S/N wrong by up to 0.23) | Built-in assertion against the brute-force scan | T3 |
| C2 | Bolometric zero point: M_bol,☉ computed with the sign of the "+5" wrong (−5.24 instead of 4.755) | Implausible L_bol ~1e-5 L☉ | T8 |
| C3 | Pecaut & Mamajek table read from the wrong column (z−Y instead of Msun) | Empty result list | T7 |
| C4 | Jacobian sign of the Gaussian ρ★ prior when sampling ln ρ (−ln ρ instead of +ln ρ) | Code review before the first fit | T9 |
| C5 | zsh does not word-split `$var`: the first batch of 8 fits started with merged arguments | Empty process list; logs removed, jobs relaunched | — |
| C6 | Frequency grid of the BLS search ~25× denser than needed (runtime only) | Runtime; replaced by a log-uniform grid with d ln P = dur_min/(baseline·5) | — |
| C7 | **This audit repeated bug A14 with 20,000 trials and first reported "FAP < 5 × 10⁻⁵" for the B link (S6 alone and S6+S32).** More correlated trials add no information. | Late check of the null percentiles: in s11 the S6-only p95 = p99 = max = 3.01. That is only possible when most trials share one window. | T14 (distinct windows: 16 per 1500 epoch-only trials vs 874 with the period randomised). Corrected by `s11b_nullB_periods.py` (period + epoch randomised) and `s11c_single_box_census.py` (rank test) |
| C8 | `s10_eb_tests.py` draws its random boxes at phases 0.1–0.9, but its phase scan covers 0.05–0.95. The "fraction of random boxes" therefore does not describe phases < 0.1 or > 0.9. | Own-photometry B scan maximum at phase 0.925 | T15: σ re-derived over 0.05–0.95 (`s10c`) agrees within 3 %, so no significance changes. Both out-of-transit features (B own \|z\| = 7.7 at 0.925; A own \|z\| = 3.8 at 0.27) come from one 4000-ppm, < 2-h dip at BTJD 2185.4 in the own S32 TESScut photometry, absent from SPOC PDCSAP. Leaving that stretch out removes them. |
| C10 | `s07_search.py` has no `__main__` guard. An exploratory `import` of it (09:19 CEST) re-ran the whole search and rewrote `out/s07_search.json`. | File time changed | No effect on results: BLS is deterministic, and all 1,102 high-precision values of the earlier printout are present unchanged in the rewritten file. |
| C9 | Wording in the first PUBLICATION_TRUTH draft: "no random box among 3000 is as deep" (§4.4). The random boxes sample only ~37 (A) or ~145 (B) independent phases of the fold. | Same check as C7 | Wording fixed. The empirical S/N (14.8 for A) stays as it is: its σ comes from ~37 independent phases, so the σ itself is uncertain by ~12 %. |
