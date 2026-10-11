# Verification of Audit 05_B_depth (signal B depth and radius)

Date: 2026-09-25. Verifier code is new, in this directory. It reads only the raw FITS files in `../../00_data/mast`. I did not import any code from the auditor or from earlier analyses. I read the auditor's scripts (`b_lib.py`, `b_run.py`, `s_configs.py`, `s02_joint_fit.py`) to look for bugs, and I read the earlier `rnaas_tic4206066_pubprep/analysis/a1_timing_depths.py` only to confirm how the A mask was applied.

## What I did

- **`v01_independent.py`: bookkeeping checks.**
  - Timing of signal A relative to signal B in S32.
  - The fraction of B in-transit cadences removed by a ±T14_A mask.
  - The S6 bit-512 cadences.
  - QLP S98 flag bits per B epoch.
- **`v02_depths.py`: per-sector depths with a method deliberately different from the auditor's.**
  - No wotan. The detrending is either none (a local polynomial per event) or my own running median of 0.75 d with the transits masked.
  - The polynomial baseline is order 1 or 2.
  - Signal A is handled four ways: a fixed 522 ppm subtraction, a free amplitude, a ±T14_A mask as in the manuscript, or a ±0.6 T14_A mask.
  - Errors are either white or "null-fit" errors: the standard deviation of the same fit at 300 random out-of-transit epoch sets in the same light curve.
  - Output: `v02_depths.json`, `v02.log`.
- **`v03_radius_grid.py`: grid posterior over (k, b, ρ*).** It uses a circular orbit, fixed limb darkening (0.32/0.31), per-event multiplicative linear baselines and A subtracted. It is run with ρ* ~ N(2.51, 0.24) and with ρ* log-uniform, then R* = 0.623 ± 0.019. Output: `v03_radius.json`, `v03.log`.

## Key independent results

**Bookkeeping (`v01_info.json`)**

- **S32 timing.** B minus A mid-time is 2.263 h (epoch −164) and 1.926 h (epoch −162). A ±T14_A mask removes 25% and 46% of the B in-transit cadences. The earlier code `a1_timing_depths.py` does use `~in_tr(t, other, 1.0)` with `|dt| < 1.0*T14`, so the mask half-width is the full T14.
- **S6 bit 512.** Exactly one cadence has it: BTJD 1484.870, 0.10 h from the predicted B mid-time, flux −2123 ppm.
- **S98 flags.**
  - Epoch −1 has bits 2048, 2^29 and 2^30 on all its in-transit cadences.
  - Epoch 2 has 2^29 on all of them (2^30 on 20%).
  - Epoch 3 has 2048 and 2^29 on all of them.
  - Epoch 4 has 2^30 on 48% (it is partial).
  - Epochs 0 and 1 are clean.

**Per-sector depths.** I ran 20 configurations; errors are null-fit errors.

| Config | S6 | S32 | S98 | nσ (S98 − archival) | common | Δχ² |
|---|---|---|---|---|---|---|
| local p1, A fixed | 866 ± 115 | 823 ± 164 | 1165 ± 89 | 2.32 | 1004 ± 67 | 5.4 |
| local p1, A ±T14 mask | 866 | 904 | 1165 | 2.13 | 1018 | 4.6 |
| local p2, A fixed | 915 ± 117 | 931 ± 145 | 1183 ± 84 | 1.93 | 1039 | 3.7 |
| local p2, A ±T14 mask | 915 | 970 | 1183 | 1.81 | 1047 | 3.4 |
| median p1, A fixed | 865 ± 113 | 778 ± 167 | 1151 ± 83 | 2.36 | 993 | 5.8 |
| median p1, A free | 865 | 736 | 1151 | 2.51 | 980 | 6.8 |
| median p2, A ±T14 mask | 921 | 990 | 1183 | 1.71 | 1060 | 3.1 |
| any, bit 512 kept | 1128–1191 | – | – | 0.64–1.00 | 1084–1137 | 2.0–4.0 |

Across all my configurations nσ is 0.6–2.5 (1.7–2.5 without the bit-512 cadence), Δχ² is 2.0–6.8, and the common depth is 980–1137 ppm.

**Radius (grid)**

- With the density prior: k = 0.0332 ± 0.0013, b = 0.878 ± 0.017 and **Rp = 2.28 ± 0.11 R⊕**.
- With ρ* log-uniform: b = 0.62 +0.28/−0.43, ρ* = 10 +11/−9 ρ☉ and **Rp = 2.05 +0.29/−0.14 R⊕**.

Both agree with the manuscript (2.25 ± 0.12; ≈2.0) and with the auditor (2.23; 2.02).

## Verdicts on the auditor's findings

| ID | Auditor | Verdict | Corrected severity | Reason |
|---|---|---|---|---|
| F1 S32 biased by A mask | MAJOR | PARTIAL | MINOR | The mask geometry is exact (25%/46%). In my fits, masking raises S32 by only 40–80 ppm at a fixed ephemeris. With the whole transit, S32 is 736–942 ppm (a quadratic baseline gives about 930; the auditor's own LOCALPOLY2 gives 934). The auditor's "710–770 for any full-transit A treatment" is overstated. The manuscript's 1003 ± 206 is about 1σ above the full-transit values, and my ±T14 mask with a quadratic baseline reaches 970–990. It is a disclosure and robustness issue, not a wrong number that changes a conclusion. |
| F2 1.6σ / Δχ² = 3.1 not robust | MAJOR | PARTIAL | MINOR | This is confirmed as choice-dependent: my central value is 2.0–2.4σ and my range 0.6–2.5σ. But 1.6σ lies inside the defensible range. The manuscript-like choices reproduce it (auditor 1.5σ; mine 1.7–1.8σ). The qualitative statement ("not significantly different; common depth not rejected at conventional levels") holds in every configuration (p ≥ 0.03). Rewording to "≈2σ (1–2.5σ depending on reduction)" is advisable but does not change any conclusion. |
| F3 S6 hinges on the bit-512 cadence | MINOR | CONFIRMED | MINOR | S6 moves from 865–921 to 1128–1191 ppm, and nσ drops to 0.6–1.0. Dropping the cadence is defensible: it lies about 8σ below the model in PDCSAP, and SPOC flags it. It is still deep (−1.5 to −1.7 ppt) in the FFI products, so it should be disclosed. |
| F4 S98 flagged-epoch count and partial error | MINOR | CONFIRMED | MINOR | Three epochs (−1, 2, 3) are flagged, not one. Only −1 and 3 carry the straylight bit 2048; epoch 2 is flagged only by the QLP bit 2^29. The partial epoch 4 (about 50% coverage) should carry about √2 × 161 ≈ 225 ppm. |
| F5 abstract "≈2 R⊕" | MINOR | CONFIRMED | MINOR | The on-target radius is 2.25–2.28 R⊕ and the free-density radius is 2.05. "≈2.2 R⊕" is better, since the other planet is given as 1.44 ± 0.06. Low priority. |
| Row: S6 depth | MINOR_DIFF | CONFIRMED | MINOR | This is the same issue as F3. The number itself (866 vs 838) matches within its error. |
| Row: S32 depth | DISCREPANT | PARTIAL | MINOR | See F1. The range is 736–990 ppm, and the manuscript value is at the high end, reached with its own masking choice. |
| Row: S98 per-event | MINOR_DIFF | CONFIRMED | MINOR | Only the partial event's error is wrong. |
| Row: S98 excess | DISCREPANT | PARTIAL | MINOR | See F2. |
| Row: Δχ² | MINOR_DIFF | CONFIRMED | MINOR | I get 2.0–6.8; 3.1 is inside that range. |
| Row: abstract | MINOR_DIFF | CONFIRMED | MINOR | See F5. |
| Row: previous 3.1σ | MINOR_DIFF | CONFIRMED | NOTE | My white-error nσ is 1.9–3.0; the draft value is superseded anyway. |

## Spot-checks of MATCH conclusions (circularity and errors)

1. **Rp with the density prior = 2.23 ± 0.12.**
   - The auditor's `s02_joint_fit.py` loads raw files through its own loaders and subtracts A with the manuscript's template values. It does not read any earlier derived tables; a grep for pubprep, deep_followup, PUBLICATION_AUDIT, SEALED or a1_results in its scripts finds nothing.
   - My independent grid gives 2.28 ± 0.11. **Confirmed.**
2. **Rp with the density free ≈ 2.0.** My grid gives 2.05 +0.29/−0.14, and b and ρ* are unconstrained. **Confirmed.**
3. **Common depth 996 ± 68 and S98 1099 ± 93.** I get a common depth of 980–1060 and S98 = 1151–1183 ± 85 across my configurations. **Confirmed within 1σ.** My S98 is systematically 50–80 ppm deeper than the manuscript's, which is within the errors.

## Issues the auditor missed or overstated

- **MINOR: the auditor's recommended S32 value is itself choice-dependent.** "770 ± 160" is a single number, but the S32 depth depends on the baseline order as much as on the A mask (p1 ≈ 780–820, p2 ≈ 930 with A subtracted). The manuscript should quote a range (≈750–950 ppm) or state its choice, rather than replace one choice-dependent number with another.
- **MINOR: the auditor's wording shifts the problem.** Its recommended "S98 is deeper by ≈2σ (1–3σ)" is fine, but the "3.3σ" upper end comes only from the auditor's own S98 aperture photometry, which gives depths 15% deeper than any archival product. The upper end should be phrased as about 2.5σ for the archive products.
- **NOTE: the Sector 6 error.** The manuscript's S6 error (114) matches my null-fit error (113–122). The auditor's κ-inflated 143 is somewhat conservative.

## Overall

The auditor's arithmetic, masks, time system (BTJD throughout) and flag bookkeeping are correct, and there is no circularity. The substantive points are real: S32 and S6 depend on choices, the S98 excess is ≈2σ rather than a firm 1.6σ, three S98 epochs are flagged, and the partial event's error is understated. But the two MAJOR findings are overstated. The manuscript's numbers are reproducible with its stated choices and fall inside the defensible range, and no conclusion of the paper changes. I rate them MINOR (wording and disclosure).
