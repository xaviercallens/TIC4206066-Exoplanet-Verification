# REVIEW_2: re-review of `tic4206066_pixel_recheck` v2

Reviewer: independent (Claude), 2026-09-25. I did not open `opus55_final_audit_tic4206066/` and edited nothing except this file.
I re-ran `run_all.sh` once. Scratch checks are in `…/scratchpad/review/v2/q1.py … q4.py`.

**Bottom line.**
- The package reproduces and the numbers are right.
- One key sentence in the Signal A conclusion is wrong: "S98 central value about 200–230 ppm with both methods". The per-event value comes from a single bad event.
- The two methods really do disagree, so the verdict "not confirmed, not excluded" still stands.
- Signal B wording is fine.
- The Note needs no change.

## 1. Reproducibility
- Before the re-run, `SHA256SUMS` was saved to scratch. `sh scripts/run_all.sh` then gave a **byte-identical** `SHA256SUMS` for all 20 entries.
  That covers all 5 `out/*.json`, the PNGs and `REVIEW_1.md`.
- This `REVIEW_2.md` will enter `SHA256SUMS` on the next run.

## 2. Numbers in SUMMARY.md vs out/*.json
All match:
- **`t9_method_replay.json`** → `*.A_r1.0_minus_r2.0`: 226.6 ± 60.1; 227 ± 66; 200 ± 73 (2.72σ); 149.5 ± 64.8; 232 ± 52; 229 ± 50 (4.55σ).
- **`after_decorrelated.json`** → `A.<set>.summary.corrected.diff_1.0-2.0`:
  - S98 +230 ± 211, n = 10;
  - S6 +102 ± 127, n = 4;
  - S32 TPF −62 ± 60, n = 8;
  - S32 FFI −36 ± 65;
  - `A_combined.corrected.diff_1.0-2.0` −15 ± 53.
- **B table:** unchanged from v1 and matches `b_s32_events.json`.
- **"3000–42000":** matches.
- **Not stored in any JSON:** "3.3 % / 0.9 %" and "about 40 mpix". They agree with my own recomputation (REVIEW_1, `r6.py`).
- **The explanation of the null σ, 64.8 vs 65.5**, is correct. T9 used one random stream across S98 and then S06, while s03b re-seeds it for each sector.

## 3. Status of REVIEW_1 findings, and new issues

| REVIEW_1 item | v2 status |
|---|---|
| C1: no sky subtraction; background regressor all zeros | **Fixed.** `common.py:68-73` subtracts a faint-pixel median, as in pubprep `a2_prf_localisation.py:90`, so the regressor is now a real per-cadence sky level. TPF handling is right: FLUX is already background-subtracted, the regressor is the sum of `FLUX_BKG`, and no TESScut time shift is applied. The shift signs match a2 L40/L84. |
| C2: S98 scattered light included | **Fixed** (`common.py:56-61`). |
| B4: T9 mis-described | **Fixed.** T9's quantity r1 − r2 is now computed. |
| **C3: S98 epoch 8 (gap and momentum dump inside the out-of-transit window)** | **Not fixed.** The 70 % rule (`s03:23-24`) tests only in-transit coverage. Epoch 8 has n_in = 14 of about 17 expected, so it passes. **It is the whole S98 per-event signal.** Per-event r1 − r2 is +2098 ppm for epoch 8 and between −132 and +212 for the other 9 events. Without epoch 8: **+22 ± 43 ppm (N = 9)**. Median of all 10: +70 ± 87 (`q1.py`). |
| C4: s04 local baseline does not mask A | **Not fixed** (`s04:18,24`, unchanged; `b_s32_events.json` hash is identical to v1). The "local baseline" row in the B table is still unstable: −164 at r = 1 goes from 1097 to −164 ppm once A is masked. The conclusion uses the s03 row, so this is minor. The row should either carry a caveat or be dropped. |
| s04/s05 `raw=True` for the TPF | **OK.** TPF FLUX is already background-subtracted, and `QUALITY == 0` keeps the cadences aligned with `POS_CORR` (`s04:11`). The only difference from s03 is the cadence mask. |
| **S6: 11×11 vs 21×21 ("cause not identified")** | **Cause found.** The pixel data are identical: same cadences, and the 11×11 is the central sub-array of the 21×21 with a maximum pixel difference of 0. In the 11×11 stamp, T9's 4–5.5 px background annulus is cut off by the stamp edge: 47 of 52 pixels in S6, 42 of 52 in S98. Clipping the annulus in the 21×21 reproduces **149.5 ± 64.8 (S6) and 226.6 ± 60.1 (S98) exactly** (`q3.py`). In S6, T9's result therefore depends on which background pixels are used, at the 80-ppm level. |
| Minor | The comment at `s03:56` still says "single-pixel (r=0.5 minus r=2.0)" but the code now combines both quantities. Cosmetic. |

Extra checks, all in scratch:
- **T9 method on clean data with the package's faint-pixel background:** S98 +191 ± 76, S6 +227 ± 52. T9's excess does not depend on how the sky is treated.
- **T9 method on S32** (T9 never tested S32):
  - FFI: +85 ± 135 with `QUALITY == 0`, +82 ± 163 with the clean selection;
  - TPF: +109 ± 118.

  S32 cannot tell the two methods apart.

## 4. Wording

**Signal A.** "Not confirmed and not excluded" is justified. But the key sentence is **not** supported:
> "In S98, however, the correction does not move the central value: about 200–230 ppm with both methods."

The per-event +230 ppm is epoch 8 alone. The clean per-event S98 value is +22 ± 43 ppm. The honest picture:
- **T9's method** gives about 190–230 ppm in clean S98 and S6, and this is robust to sky treatment and cadence selection.
- **The per-event local method** gives about zero: S98 +22 ± 43 without epoch 8, S32 TPF −62 ± 60, S6 +102 ± 127.
- **In S98 the two methods differ by about 2σ** (178 ppm, combined error about 85). What drives the difference is unknown. The in-transit window (full T14 vs ±T14/4) and the detrend (global 1-d vs a local 9-parameter fit) are the candidates.

Also:
- "not significant in any sector" is true only for the per-event method. T9's method gives 4.4–4.6σ in S6 on the 21×21 cut-out, although S6 is the sector hit by momentum dumps. This should be stated.
- Suggested fix: remove epoch 8, or require the out-of-transit data to have no gap or momentum dump. Then say "the per-event method finds no excess in any sector. T9's method finds about 200 ppm in S6 and S98 but is not decisive in S32. The discrepancy between the methods is unresolved."
- Replace "cause not identified" for S6 with the annulus-clipping explanation.

**Signal B.** Justified as written. It is limited to 2 events, says the centroid test is not decisive, notes the −164 drift, and cites pubprep's own exclusion of −164.
Only change: add a caveat to, or remove, the "local baseline" row (C4).

## 5. The Note (`rnaas_tic4206066_pubprep/manuscript.tex`)
Nothing in v2 changes REVIEW_1 §E:
- The Note makes no statement about aperture dependence, the core pixel, T9, or B's behaviour at pixel level.
- Its statements that can be checked here agree with the package. These are l.39 ("every transit falls 0.4--3.5\,h from a reaction-wheel desaturation"), the B S32 depth at l.35 ($1003\pm206$\,ppm vs per-event 949–999), and the localisations at l.39.
- The method disagreement above concerns an effect where the small aperture looks deeper, not an off-target signature. It stays an optional item for `UNRESOLVED_LIMITATIONS.md`.

Manuscript change required: NO
