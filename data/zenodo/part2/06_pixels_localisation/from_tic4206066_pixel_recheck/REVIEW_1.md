# REVIEW_1 — independent review of `tic4206066_pixel_recheck`

Reviewer: independent (Claude), 2026-09-25. Scope: package files, the T9 source, and the RNAAS pubprep files only.
`addenda/opus55_final_audit_tic4206066/` was not opened. Nothing in the package was edited. `run_all.sh` was re-run once.
Scratch checks (not part of the package) are in
`/private/tmp/claude-501/-Users-paulsmith-space--3/7d7db2d2-c496-4e2c-8bef-565716340671/scratchpad/review/r1.py … r10.py`.

**Bottom line.** The package reproduces exactly and its tables match its JSON. The B/S32 retraction holds.
The S98 part of the Signal A conclusion is weaker than stated, for three reasons:
- the FFI cut-outs are not background-subtracted;
- S98 scattered-light windows are included;
- the comparison with T9 uses a different aperture.

Once S98 is restricted to clean events, it is precise enough to be in 2–3σ tension with T9. However, T9's own method still finds
the excess in the same clean data. The anomaly therefore stays open, but it is now a disagreement between methods, not a lack of precision.
None of this touches the Note.

---

## A. Reproducibility

- Before the re-run, `SHA256SUMS` was copied to scratch. `sh scripts/run_all.sh` then ran in about 15 s with the project `.venv`
  (Python 3.12.13, numpy 2.5.3, astropy 8.0.1).
- The new `SHA256SUMS` is **byte-identical** to the saved one. This covers all 17 entries, including all 4 `out/*.json` and all 3 PNGs.
- Console output equals the SUMMARY tables.
- Note: `s06_hashes.py` hashes every file in the package. A future re-run will therefore add this `REVIEW_1.md` to `SHA256SUMS`.

## B. Numbers (SUMMARY.md vs out/*.json)

All table values match the JSON after rounding:
- Signal A table: `after_decorrelated.json` → `A.<set>.summary.corrected.diff_*`, combined result under `A_combined_single_pixel_minus_r2.corrected`.
- "before" values: `before_global_detrend.json` → `apertures`.
- B table: `b_s32_events.json` and `after_decorrelated.json` → `B.S32_tpf2min.events`.
- The S98 raw scatter of ±769 ppm: `summary.raw.diff_0.5-2.0.err_empirical`.

The ephemeris matches `ephem_current.json`. The S6 dump distances (0.4–3.5 h for the 5 events used) match `t7_momentum_dumps.json`.

Mismatches or imprecisions:
1. **"every pixel 'brightening' by 3000–40000 ppm"**. `before_global_detrend.json` → `S98_ffi.pixels.*.A` spans −3041 to −41827 ppm, so the range is 3000–42000. Trivial.
2. **"core pixel's flux changes by 4.9% (−164) and 3.8% (−162)"**. These match `core_pixel_range_pct_over_window` (4.89, 3.77),
   but that key is the max − min of the unbinned 2-min points. For −162 the per-point scatter is 0.61%, so the range is mostly noise:
   - 30-min binned range: 3.3% for −164, **0.9% for −162**;
   - binned centroid range: y ≈ 39 mpix for −164, **3.6 mpix for −162** (scratch `r6.py`).

   The "about 40 mpix" is correct for −164 but is not stored in any JSON. "30–50 times the transit depth" holds only for −164, at about 35×.
3. The QLP 1×1/3×3 ratio 0.54 for B is cited to `candidates/tic4206066_p3.182945/candidate.json`. It actually lives in
   `candidates/tic4206066_p11.132976/candidate.json` (`metrics.depth_ratio_qlp_det_1x1_to_3x3` = 0.544). The value is correct.
4. **T9 is mis-described.** SUMMARY §Question calls T9 the "single brightest pixel". T9's key `S98.A_r1.0_minus_r2.0` uses
   `ap = rr <= r` (`independent_verification_tic4206066_SEALED/scripts/t9_aperture_diff.py:19`), so r = 1.0 is a **5-pixel plus-shaped aperture**,
   not one pixel. The mislabel starts in `verification_report.md:30` ("1-px core aperture"). The package's "single pixel" is r = 0.5, a different quantity.
   The like-for-like package number, r1 − r2, is not in `out/`. From the events it is **+105 ± 180 ppm** for S98 corrected (scratch `r1.py`).
   So "S98 central value +233 ppm is close to [227]" compares two different apertures.

## C. Code review

Checked and correct:
- **Signs:** depth = 1 − mean; differences are small minus large; centroid shift is in − out.
- **Other-signal masking** in the s03 out-of-transit window (`s03:21`).
- **No A/B in-transit overlap in any event used.** Only S32 A −573/−566 have a B transit 1.9–2.3 h later, and it is masked.
- **Event counts:**
  - S32 has 8 A events and 2 B events; the third B transit, −163, falls in the data gap.
  - S6 uses 5 of the 7 A windows in range.
- **Combination:** S32 FFI is correctly left out of the combined estimate because it shares photons with the TPF.

| # | issue | evidence | changes a conclusion? |
|---|---|---|---|
| C1 | **FFI cut-outs are not background-subtracted, and the "background" regressor is all zeros.** S6, S32 FFI and S98 have `BACKAPP = 0` and `FLUX_BKG ≡ 0`. The annulus sky is 120–135 e⁻/s/pix, and the S98 95th percentile is 463. `lstsq` silently drops the zero column. SUMMARY §Methods ("…and the background") and "adding the centroid and background terms…" are therefore wrong for all FFI sets. Only the S32 TPF has a real background, and its FLUX is already background-subtracted. Unsubtracted sky dilutes wide apertures more than the core, which biases the result towards "core deeper". | `common.py:40,45`; `s03:23-24`; FITS headers (scratch `r3.py`) | **Yes, for S98** (see C2/C3) |
| C2 | **S98 scattered-light windows are included.** TESScut `QUALITY == 0` does not carry QLP bit 2^29. The Note and all project analyses remove those intervals (`manuscript.tex:27`; `MANUSCRIPT_NUMBER_AUDIT.csv` row 16: 9 clean S98 A transits). Here 5 of the 15 S98 events used lie in flagged windows: epochs 2, 3, 6, 11, 12; epoch 7 has no data. At epochs 2, 3, 11 and 12 the r = 2 flux is 9–46% above its median. **This is the cause of the S98 "before" failure** that SUMMARY calls unidentified. Those four windows give −26,000 to −56,000 ppm at r = 2, while all other windows give +190 to +1200 ppm. The effect grows with aperture area, as additive sky should (scratch `r2.py`, `r10.py`). | `common.py:41` | **Yes, for S98** |
| C3 | **No rejection of events at gaps or edges.** S98 epoch 8 (BTJD 4024.20) straddles a 0.23-d gap next to a momentum dump. Its single-pixel depth is 5993 ppm, and it dominates the S98 scatter. | `after_decorrelated.json` → `A.S98_ffi.events[epoch=8]` | Yes, for S98 |
| C4 | **s04 "local baseline" does not mask A.** Both B events lie 1.9–2.3 h after an A mid-transit, inside B's out-of-transit window. Masking A moves −164 r = 1 from **1097 to −164 ppm**. The local-baseline row is therefore unstable under the drift. The conclusion rests on the s03 "per-event correction" row, which does mask A. | `s04:18,24` (scratch `r9.py`) | No, but the local-baseline row should be flagged as unreliable |
| C5 | **Formal errors are white-noise only.** They use `std(res[oot])/sqrt(n_in)`, with no degrees-of-freedom correction for the 9 parameters and no variance from extrapolating the fit into the transit. In S6 this is 9 parameters on about 25 out-of-transit points. Empirical errors are 2–5× the formal ones (S98 r = 0.5: 434 vs 94; S32 TPF: 316 vs 152). The B "use the per-event formal errors" advice therefore gives lower bounds. | `s03:32`, `s04:20` | No: the B r1/r3 differences stay ≪ 1σ |
| C6 | **Centroid terms can absorb a real in-transit centroid shift.** They are fitted out of transit and applied in transit, so a blend's centroid shift could be partly removed. | `s03:24` | No: `raw` (no centroid terms) gives the same verdict, combined +357 ± 300 |
| C7 | **S6 relaxed window.** ±0.45 T14 with 30-min integrations reaches about ±65 min, beyond T14/2 ≈ 56 min, and `min_in = 2`. Absolute depths are diluted, and the fit has 9 parameters on about 25 points. S6 is also the sector the project calls unusable for A. | `s03:51` | No: S6 has little weight. Dropping it gives +196 ± 229 |
| C8 | **Combined estimate.** It uses inverse-variance weights from empirical errors with N = 5/8/15, whose own uncertainty is 25–35%. It combines the single-pixel − r2 difference, which is not the T9 quantity. | `s03:55-60` | No |

**Effect of C1–C3 on S98** (scratch `r4.py`, `r7.py`, `r8.py`; r1 − r2 is the T9 quantity):

| S98 variant | N | r1 − r2 (ppm) | single pixel − r2 |
|---|---|---|---|
| package as is | 15 | +105 ± 180 | +233 ± 389 |
| package method, 9 clean events (no QLP 2^29, no epoch 8) | 9 | **+84 ± 60** | +222 ± 104 |
| same + annulus sky subtraction, 6–9 px | 9 | +12 ± 42 | −5 ± 144 |
| same + annulus sky subtraction, 4–5.5 px (T9-like) | 9 | +59 ± 62 | +147 ± 112 |
| **T9 method**, reproduced exactly (226.6 ± 60.1), all cadences | — | +227 ± 60 | — |
| **T9 method, scattered-light cadences removed** | — | **+207 ± 72** | — |

With clean data, the per-event method is about as precise as T9, not 6× worse. It sits 1.7–2.9σ below T9's value, but T9's own method on the same clean cadences still gives +207 ± 72.
The anomaly is therefore **method-dependent and unresolved**. The two methods use different in-transit windows: T9 uses the full T14 and a global 1-d detrend,
while the package uses ±T14/4 and a per-event local fit. The background choice alone moves the result by about 50 ppm.

## D. Method and wording

**(1) Signal A.** "Not confirmed, not excluded" is the right bottom line. The supporting text is not:
- **The lead sentence suggests pointing explained the T9 excess.** It reads "After local, per-event correction for pointing drift, the aperture dependence … is not confirmed".
  However, T10 found no pointing excursion at S98 transit times (`verification_report.md:31`), and in S98 the correction did not move the value (raw +180 → corrected +233).
  The trend the correction did remove is the package's **own first-pass** trend in S6 and S32, not T9's. These are two separate anomalies and should be described separately.
- **"S98 is not precise enough (±389)"** is an artefact of C1–C3. It should say that with the project's clean-interval convention S98 gives r1 − r2 ≈ +10 to +85 ± 40–60 ppm,
  against T9's 227 ± 60 and against T9's own method on the same clean data, 207 ± 72. Open.
- **The single-pixel vs r ≤ 1 (5-pixel) mislabel** (B4) should be corrected.
- **The empirical errors** are reasonable in principle, since the formal ones are 2–5× too small, but mean ± SEM is not robust here:
  - N = 5 in S6;
  - S98 is dominated by outliers: median-based single pixel − r2 is +212 ± 166, against a mean of +233 ± 389.
  - Report robust values and apply event quality cuts.
- **PRF consistency.** The A position (+8 ± 5, −8 ± 7)″ (`manuscript.tex:39`; audit (+2.5, −7.5)″ ± 5.7″ in `PUBLICATION_TRUTH.md:147-148`)
  is an offset of about 0.5 px. The aperture test has ±90 ppm on r1 − r3 in S32, so it cannot tell a 0.5-px offset from the target. It neither supports nor contradicts the localisation.
  A "core deeper" pattern is also the opposite of the usual off-target signature (`verification_report.md:30`).

**(2) Signal B.** "Methodological false positive" is justified. The test only needs to re-examine the same two events that produced the first-pass claim,
and per event r = 1 and r = 3 agree after correction (924 vs 949 and 1061 vs 999 ppm). Wording fixes:
- "caused by global detrending during strong pointing drift": strong drift is shown only for −164, and −162 drifts about 0.9% (B2). Say "attributable to…, demonstrated for epoch −164".
- **"no centroid shift during transit"** is not a constraining test. The errors are up to 2.6 mpix, while T10 puts a 43.5″ blend at about 1 mpix.
  The support comes from aperture-depth agreement and the PRF position (+1 ± 3, +1 ± 4)″, not from the centroids.
- **Cite the prior flag on −164.** `rnaas_tic4206066_pubprep/METHOD_AUDIT.md:58` already excluded own-photometry S32 B event −164 as "corrupted by a centroid excursion".
  That is consistent with this package and should be cited.
- **Use the s03 row, not the local-baseline row** (C4).

## E. Manuscript impact (`rnaas_tic4206066_pubprep/manuscript.tex`, not submitted)

The Note contains **no** aperture-dependence, core-pixel-depth or T9 statement, and never mentions the B/S32 off-target effect. Lines checked:
- l.20 "pixel-level localization, Gaia astrometry, … leave a planet on a bound companion … as the main untested false-positive scenario" — **NONE.**
  The package does not test localisation. A "core deeper" excess is not an off-target signature.
- l.31 "Per-sector depths agree within errors" — **NONE.** Package box depths (±T14/4, about 430 ppm in S32) are not comparable, and SUMMARY says so.
- l.35 "Per-sector depths are $838\pm114$, $1003\pm206$, and $1099\pm93$\,ppm" — **NONE.** The package's S32 B depths (r = 3 corrected 949 ± 151 and 999 ± 176) are consistent with 1003.
- l.39 "Sectors 32 and 98 place the source at $(+8\pm5, -8\pm7)''$ … Sector 6 is unusable because every transit falls 0.4--3.5\,h from a reaction-wheel desaturation … The 11.13\,d signal lies at $(+1\pm3, +1\pm4)''$" — **NONE.** The package confirms the 0.4–3.5 h figure and does not constrain the positions.
- `MANUSCRIPT_NUMBER_AUDIT.csv`: no row is affected.
- **OPTIONAL, outside the Note:** `UNRESOLVED_LIMITATIONS.md` has no entry for the open T9 S98 aperture excess. One line could record:
  - the excess is method-dependent: T9 method +207 ± 72 on clean data, per-event method +12 to +84 ± 42–62;
  - it is unresolved;
  - it does not point to an off-target source.

  Putting it in the Note is not needed.

**Manuscript change required: NO.** The Note makes no statement about aperture dependence, core-pixel depth or B's pixel behaviour.
Its localisation, S6 and B-depth statements agree with this package. The package's own S98 conclusion needs rewording (D1) before it is cited anywhere.
