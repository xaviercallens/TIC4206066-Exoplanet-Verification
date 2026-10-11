# TIC 4206066 — pixel re-check (aperture dependence of the transit depth)

Status: exploratory re-check, 2026-09-25. It does not change and was not merged into `../tic4206066_MASTER_EVIDENCE/`
or either RNAAS package. Nothing was submitted.

Version 3. Two independent reviews shaped it, and their findings were verified and fixed:
- [REVIEW_1.md](REVIEW_1.md): the TESScut cubes had no background subtracted; T9's quantity was misdescribed; S98 scattered-light intervals were not removed.
- [REVIEW_2.md](REVIEW_2.md): one S98 event with a data gap and a momentum dump in its window (epoch 8) carried the whole per-event S98 value; the B local baseline did not mask A transits.

The fixes made after REVIEW_2 have not been re-reviewed.

## Question
The independent verification (T9) reported that in S98 Signal A is 226.6 ± 60.1 ppm (3.77σ) deeper in an r ≤ 1 px aperture than in an r ≤ 2 px aperture. The r ≤ 1 px aperture is 5 pixels, not a single pixel. The comparison used an annulus background and a random-phase null.
- Sources: [t9_aperture_diff.json](../independent_verification_tic4206066_SEALED/out/t9_aperture_diff.json) key `S98.A_r1.0_minus_r2.0` and [t9_aperture_diff.py](../independent_verification_tic4206066_SEALED/scripts/t9_aperture_diff.py) L18-L19.
- In S6 T9 gave 149.5 ± 65.5 ppm (2.3σ).

QLP's 1×1/3×3 depth ratio points the other way:
- 0.46 for A ([p3.182945/candidate.json](../../candidates/tic4206066_p3.182945/candidate.json));
- 0.54 for B ([p11.132976/candidate.json](../../candidates/tic4206066_p11.132976/candidate.json)), key `metrics.depth_ratio_qlp_det_1x1_to_3x3`.

A first pass in this session also reported an apparent 4σ "off-target" pattern for Signal B in the S32 2-min data.

## Inputs and processing
Inputs, all existing project files, are listed with sha256 in [`out/inputs.json`](out/inputs.json):
- the TESScut 21×21 cut-outs for S6, S32 and S98;
- the SPOC S32 2-min TPF;
- T9's own 11×11 cut-outs;
- the QLP S98 light curve, used for its scattered-light flags.

Ephemeris: [`ephem_current.json`](../tic4206066_PUBLICATION_AUDIT/out/ephem_current.json) (publication audit):
- A: P = 3.182785293657263 d, T0 = 3998.7367930333035 BTJD;
- B: P = 11.132745380039134 d, T0 = 4000.8658739613074 BTJD.

Cadence selection and background are as in pubprep `a2_prf_localisation.py` (L39-L91):
- SPOC default bitmask; in S98 also bit 2048 and the QLP 2^29 scattered-light intervals;
- the TESScut cubes are not background-subtracted (header `BACKAPP = 0`, `FLUX_BKG = 0`), so a per-cadence background from faint pixels more than 5 px from the target is subtracted;
- the TESScut time offsets are applied.

Scripts:
- [`s02`](scripts/s02_before_global_detrend.py), **first pass (flawed, kept for provenance)**: QUALITY = 0 only, no background subtraction, global 0.75-d running-median detrend.
- [`s03`](scripts/s03_after_decorrelated.py), **per-event local correction**: a ±0.35-d out-of-transit fit with a quadratic in time plus flux-weighted centroid terms (x, y, x², y², xy) and the background level, following the centroid-decorrelation idea of Vanderburg & Johnson (2014).
  - Dropped events: less than 70 % in-transit coverage, a data gap > 2.4 h in the ±0.35-d window, or a momentum dump in the window ([T7](../independent_verification_tic4206066_SEALED/out/t7_momentum_dumps.json)).
  - Every S6 A transit lies 0.4–3.5 h from a dump, so S6 has no usable events under this rule; it is kept only as a sensitivity case without the dump rule.
  - Aperture differences are taken per event; the error is the scatter across events.
  - `raw` repeats the fit without the centroid and background terms.
- [`s03b`](scripts/s03b_t9_method_replay.py), **T9 replay**: T9's own detrender, imported read-only, and T9's apertures, background and null, run on (a) T9's selection and (b) the clean selection above.
- [`s04`](scripts/s04_b_s32_events.py): the two S32 2-min B events, one by one.

## Results — Signal A, quantity of T9 (r ≤ 1 px minus r ≤ 2 px)

| data | T9 method, T9 selection | T9 method, clean selection | per-event local correction (clean; gap and dump rule) |
|---|---|---|---|
| S98 | 226.6 ± 60.1 (T9's own 11×11 cut-out, reproduced exactly); 227 ± 66 (21×21) | **200 ± 73 ppm (2.7σ)** | **+20 ± 49 ppm, 8 events** |
| S6 | 149.5 ± 64.8 (11×11); 232 ± 52 (21×21) | 229 ± 50 ppm (4.6σ) | no usable events (all near dumps); without the dump rule +102 ± 127 ppm, 4 events |
| S32 2-min TPF (not tested by T9) | — | — | **−62 ± 60 ppm, 8 events** |
| S32 FFI (same photons as the TPF) | — | — | −36 ± 65 ppm, 8 events |

Sources: [`out/t9_method_replay.json`](out/t9_method_replay.json) and [`out/after_decorrelated.json`](out/after_decorrelated.json).

Notes:
- The S6 difference between the 11×11 and 21×21 cut-outs comes from the 11×11 stamp: its edge clips T9's 4–5.5 px background annulus. The pixel data are identical ([REVIEW_2.md](REVIEW_2.md)).
- Before the gap and dump rule, the per-event S98 value was +230 ± 211 ppm. All of it came from epoch 8, whose own r1 − r2 is +2098 ppm.

## Results — Signal B, S32 2-min ([`out/b_s32_events.json`](out/b_s32_events.json), [plot](plots/B_S32_events_core_vs_wide.png))

| | epoch −164 | epoch −162 |
|---|---|---|
| nearest momentum dump ([T7](../independent_verification_tic4206066_SEALED/out/t7_momentum_dumps.json)) | 105.5 h | 92.9 h |
| depth, r = 1 / r = 3 px, local median baseline, A masked (unstable under drift) | −164 / 921 ppm | 957 / 1052 ppm |
| depth, r = 1 / r = 3 px, per-event correction | 924 ± 209 / 949 ± 151 ppm | 1061 ± 256 / 999 ± 176 ppm |
| core-pixel flux range over ±8 h, 30-min bins | 3.3 % | 0.9 % |

- Epoch −164 has a strong pointing drift: the flux-weighted centroid moves by about 40 mpix. The same event is excluded from pubprep's own photometry as "corrupted by a centroid excursion" ([METHOD_AUDIT.md](../rnaas_tic4206066_pubprep/METHOD_AUDIT.md) L58).
- The first-pass global detrend did not follow this drift. It produced the apparent "B missing from the core" pattern: core −1619 ± 580 ppm, r = 3 px 1075 ± 155 ppm, in [`out/before_global_detrend.json`](out/before_global_detrend.json).

## Conclusions (cautious)
1. **Signal A.** The aperture dependence **is not confirmed after local, per-event correction for pointing drift**:
   - S98: +20 ± 49 ppm (8 clean events);
   - S32 2-min: −62 ± 60 ppm;
   - no usable S6 events.

   It is **not excluded**, because T9's own global-detrend method still finds about 200 ppm on the same clean data: S98 200 ± 73 ppm, S6 229 ± 50 ppm. In S98 the two methods differ by about 180 ppm, about 2σ, and the cause is not identified.

   So S98 is no longer the limit. Its per-event precision (±49 ppm) matches T9's (±60–73 ppm). What is open is the disagreement between the two methods. The first-pass claim that the pointing correction "removes" the anomaly was premature; the version-2 statement that the per-event S98 value was "+230" was driven by one bad event.
2. **Signal B.** The apparent off-target effect in S32 is a **methodological false positive**. It was caused by global detrending during a strong pointing drift in epoch −164. Per event, B has the same depth in small and wide apertures, and no momentum dump lies within 90 h.
   - This rests on only 2 events.
   - The per-event corrected depths are used; the local-median baseline is unstable for epoch −164 because of the drift.
   - The in-transit centroid shifts (x −0.8 ± 0.8, y +1.4 ± 2.6 mpix; x −0.1 ± 0.6, y −0.2 ± 0.8 mpix) are not by themselves decisive: a source at 43.5″ would move the centroid by about 1 mpix (T10).
   - The statement is about the method, not a new localisation of B.
3. **The first pass was wrong for three reasons:** no background subtraction of the TESScut cubes, no scattered-light cut in S98, and a global detrend. Its S98 output, with every pixel "brightening" by 3000–42000 ppm, is unphysical and not used.
4. **Limitations.**
   - This is a box-depth test, not a transit-model fit.
   - The audit ephemeris is used; the pubprep ephemeris differs by minutes, which was not tested.
   - Full pixel-level decorrelation was not done.

## Manuscript impact
According to both reviews ([REVIEW_1.md](REVIEW_1.md), [REVIEW_2.md](REVIEW_2.md)), the current Note ([manuscript.tex](../rnaas_tic4206066_pubprep/manuscript.tex)) makes no statement about aperture dependence or the core pixel, so no change is required. Optionally, the open T9 question could be listed in pubprep `UNRESOLVED_LIMITATIONS.md`. That is not done here.

## Reproduce
`sh scripts/run_all.sh`. It uses the project's `.venv` read-only: Python 3.12.13, numpy 2.5.3, astropy 8.0.1, scipy 1.18.1, matplotlib 3.11.2. The full freeze is in `env/pip_freeze.txt`. Checksums are in `SHA256SUMS`.
