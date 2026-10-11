# Verification of Audit 07 (Signal A shape, EB tests, planet radii)

Verifier: adversarial re-check of `07_A_shape_radii` (the auditor's structured result; no REPORT.md was written there).
Everything below was recomputed with new code in this folder (`v*.py`). I read the auditor's scripts to look for bugs. I read the earlier analyses (`tic4206066_PUBLICATION_AUDIT`) only to find where the manuscript's numbers came from. No derived product from either one was used as an input to my calculations.

## Summary

The auditor's central finding holds, and it is somewhat stronger than stated. The radius depends on the detrending, and the manuscript's value sits at the **high end** of the reasonable range.

The main driver is a single edge transit, Sector 98 epoch -1 (BTJD 3995.555). It lies at the end of a 0.45-d data segment that starts after the scattered-light interval. That segment sits on a ramp of about 4000 ppm/d and has no baseline after egress. Its fitted depth ranges from -109 to +1285 ppm depending on the detrending. The detrending style used for the manuscript gives about 1285 ppm, 2.5 times the mean depth: the trend is interpolated under a +-1.5 T14 mask after removing the masked points.

My local-baseline fits without that transit give:
- k = 0.0200-0.0201 (BIC-chosen or quadratic baselines);
- k = 0.0210 (linear baselines, which BIC rejects in 3 of 23 windows);
- Rp = 1.36-1.43 Re.

The manuscript quotes k = 0.0212 +- 0.0007 and 1.44 +- 0.06 Re.

Other findings:
- The mislabelled "1.50 Re with free LD" is confirmed.
- The phase-0.5 sign error is confirmed.
- The rho limit changes by only +-0.1-0.2 with the prior or detrending. I downgrade this to NOTE.
- The auditor's red-noise V-shape delta-chi2 of 21 is specific to its own biweight light curves. With my local baselines, beta is 1.0-1.17 and delta-chi2 is 31.5-35 after rescaling, so the manuscript's 32 stands. Downgraded to NOTE.
- No circularity was found in the auditor's work: it used its own code and its own VizieR query for Claret, and it takes k_B from the earlier chains by design.

## Data and provenance

- Light curves, read from `00_data/mast` (sha256 matches `DATA_MANIFEST.json`):
  - S6 TESS-SPOC PDCSAP: `ac7736ea...882e`
  - S32 SPOC 2-min PDCSAP: `764dad34...b0e7`
  - S98 QLP SYS_RM_FLUX / DET_FLUX: `a0341366...188e`
  - QUALITY == 0 in all three. For S98 this drops bits 2^11 (straylight), 2^29 and 2^30.
  - The MIT QLP page (`00_data/docs/qlp/mit_qlp.html.txt`) says bits 16-30 (1-indexed) are "unused", yet 10,676 S98 cadences carry 2^29. This confirms F7.
- Claret (2017) TESS quadratic LD coefficients, queried by me from VizieR J/A+A/600/A30:
  - `tableab` (Type q, PHOENIX), sha256 `900ebba3...`, retrieved 2026-09-25T12:29:46Z;
  - `table25` (ATLAS Z=0, xi=2), sha256 `d3d430ee...`, retrieved 12:30:19Z;
  - bilinear interpolation at Teff 4012 K, logg 4.631 gives LSM u = (0.316, 0.313) for PHOENIX and (0.375, 0.278) for ATLAS. This matches the auditor.
- Python: `audit_env`. No packages were installed. All jobs ran with at most 2 processes and OMP_NUM_THREADS=1.

## Method (independent of the auditor's)

- `v02_prep.py` builds a per-transit **local** baseline:
  - uses out-of-transit data (|dt| > 0.75 T14) within +-0.45 d, from the **same contiguous segment only** (the light curve is split at gaps > 0.1 d);
  - removes the 11.13-d signal's transits;
  - the polynomial order is 1, 2, or chosen per window by BIC (`o-1`);
  - the auditor instead used a biweight, or linear baselines profiled over +-0.30 d that ignore gaps.
- `v03_fit.py` fits batman + emcee (stretch move, 40 walkers, 15,000 steps, N/tau 27-108):
  - supersampling x15 for S6, x3 for S32 and x5 for S98;
  - per-sector jitter;
  - density prior: split-normal Mann, uniform, or log-uniform;
  - limb darkening: Claret prior N(0.316, 0.1) and N(0.313, 0.1), or free Kipping q1, q2.
- `v04_tests.py` measures template-scaled depths with local quadratic baselines fitted simultaneously. Empirical sigma comes from about 220 randomly shifted ephemerides. It also measures the phase-0.5 box.
- `v05_vshape.py` tests V shape: supersampled trapezoid vs triangle, with per-sector white sigma and with time-averaging beta.
- `v06_epoch_m1.py` compares four detrendings for each S98 epoch.
- `v07_secscan.py` scans for a secondary eclipse, requiring a baseline of at least 0.15 d on both sides.
- `v08_bic_order.py` picks the baseline order by BIC.
- `v00*.py` handles headers, flags, coverage and segments.

## Results: manuscript vs recomputed

| quantity | manuscript | auditor | verifier (this folder) |
|---|---|---|---|
| transits in fit | 23 | 23 (>=50% T14), 24 any | 23 with >=50% T14 (S6 6, S32 8, S98 9); 24 with any in-transit point; 21 with a baseline on both sides |
| k (Mann prior, Claret LD) | 0.0212 +- 0.0007 | 0.02095 raw; 0.0197 biweight; 0.0206 / 0.0203 excluding S98 -1 and 8 | all 23: 0.0195 (quadratic), 0.0201 (BIC order), 0.0195 (DET_FLUX); without S98 -1: 0.0200 (quadratic), **0.0201 +- 0.0007 (BIC)**, 0.0210 (linear), 0.0199 (DET_FLUX) |
| b | 0.35 +0.10/-0.14 | 0.33 +0.09/-0.12 | 0.31 +0.12/-0.17 |
| T14 | 1.87 +- 0.05 h | 1.88 +- 0.04 h | 1.886 +- 0.060 h |
| Rp | 1.44 +- 0.06 Re | 1.42 +- 0.06 (raw), 1.33 (biweight) | 1.36-1.37 +- 0.06 (quadratic / BIC / DET), 1.43 (linear); all 23 with quadratic: 1.33 |
| Rp, free LD | 1.50 Re | 1.35-1.45 | 1.367 (Mann prior), 1.385 (uniform rho + free LD); +0.01 relative to the Claret-prior fits |
| rho* 97.7%, density free | < 3.5 | 3.40-3.72 | 3.51 (log-uniform), 3.58 (uniform), 3.62 (uniform, DET_FLUX), 3.69 (uniform + free LD) |
| per-sector depths (template scaling, empirical sigma) | agree | chi2 0.24/2 and 2.9/2 | 514+-61, 445+-74, 468+-47 ppm; chi2 = 0.59/2 |
| odd / even | 536+-40 / 466+-43 | 539+-47 / 449+-51 | 488+-47 / 446+-52 (0.6 sigma) |
| phase 0.5 (flux, negative = dip) | "-20 +- 32" | +12 +- 31 | -9 +- 30 (empirical); source value is a **depth** of -20.4, i.e. flux +20 |
| secondary scan (width T14) | - | min 2.4-2.9 sigma | min -61 ppm at phase 0.59 (1.7 sigma) |
| V-shape delta-chi2 | 32 | 34.2 white / 21.2 beta | 35.4 white / 31.5 beta (all 23); 38.1 / 35.0 without S98 -1; beta = 1.00, 1.14, 1.05 |
| Signal B Rp | 2.25 +- 0.12; ~2.0 free density; abstract ~2 | 2.246 / 2.011 | 0.03305 x 0.623 x 109.076 = 2.246; 0.0296 gives 2.011 (arithmetic verified; k_B comes from the earlier chains) |

Per-epoch S98 depths (ppm, manuscript template; `out/s98_epoch_detrend_diag.json`). Detrending (a) mimics the manuscript: trend interpolated under a +-1.5 T14 mask. (b) is a biweight with the wotan mask. (c) is a local quadratic within the same segment. (d) is a raw linear fit over +-0.3 d.

| epoch | (a) | (b) | (c) | (d) |
|---|---|---|---|---|
| -1 | **1285** | **-109** | **78** | **905** |
| 0 | 503 | 527 | 502 | 480 |
| 1 | 555 | 505 | 490 | 671 |
| 4 | 385 | 403 | 351 | 386 |
| 5 | 549 | 566 | 563 | 520 |
| 9 | 460 | 451 | 443 | 446 |
| 10 | 449 | 453 | 413 | 449 |
| 13 | 541 | 530 | 538 | 534 |
| 14 | 419 | 401 | 406 | 397 |

Every epoch except -1 agrees to within about 60 ppm across methods; epoch 1 under the linear baseline, which BIC rejects, is the one exception.

## Verdicts on the auditor's findings

- **F1** (k and Rp depend on detrending): **CONFIRMED, MINOR** (upper end).
  - The spread is real.
  - The auditor's "reproduction" of 0.02095 includes epoch -1 at 898 ppm.
  - Without epoch -1, the reasonable methods give k = 0.0199-0.0210 (mean about 0.0203).
  - Recommended: k = 0.0203 +- 0.0010 and Rp = 1.38 +- 0.08 Re. In the abstract, "~1.4 Re".
  - State that the S98 edge transit (epoch -1) is excluded, or that it is poorly constrained.
- **F2** ("1.50 Re with free LD"): **CONFIRMED, MINOR**.
  - The source (`dynesty_A_pipe.json`) had free LD **and** a free, log-uniform density (median b = 0.62, median rho = 1.59).
  - Freeing LD alone moves Rp by only about +0.01 Re.
- **F3** (rho limit not robust): **PARTIAL, NOTE**. The limit is 3.4-3.7 across priors and detrendings, so 3.5 lies in the range. Either value still implies a host above about 0.52-0.54 Msun, which the 0.38 Msun overluminosity limit already excludes, so no conclusion changes. "<~3.6" would be the safer wording.
- **F4** (phase-0.5 sign): **CONFIRMED, MINOR**. The earlier `box()` returns 1 - weighted mean flux, so -20.4 is a depth (flux +20). Every estimate is consistent with zero.
- **F5** (V-shape delta-chi2 without red noise): **PARTIAL, NOTE**.
  - The value of 21 depends on the auditor's biweight light curves: beta for S98 is 1.62, inflated in part by the badly detrended epoch -1.
  - With local baselines, delta-chi2 is 31.5-35 after beta rescaling.
- **F6** (abstract "~2 Re" vs 2.25): **CONFIRMED, MINOR** (wording; low priority).
- **F7** (the 2^29 bit): **CONFIRMED, NOTE**, from the MIT QLP page ("Bits 16-30 are unused").
- **F9** (odd/even errors): **CONFIRMED, NOTE**. The empirical errors are 47 and 52 ppm, against the quoted 40 and 43.

## Missed issues

1. **MINOR: the manuscript's own k is probably inflated by S98 epoch -1.** The effect is worth about +0.0007 in k and +0.05 Re in Rp. The same epoch is odd, so it also inflates the odd depth of 536 ppm. Recommend excluding it (and quoting "22 transits" for the shape fit), or showing robustness to it.
2. **NOTE: baseline order.** BIC prefers a constant or linear baseline in 20 of 23 windows and a quadratic in 3: S32 -573 and -569, and S98 1. Using linear everywhere raises k by about 0.001 (fit o1 vs o-1).
3. **NOTE: "depth ~500 ppm".** Without the edge transit, the mid-transit model depth is about 460-475 ppm. "~500" is acceptable only as a round number.
4. **NOTE: S98 photometry method.** Since Sector 94, QLP light curves are made with TGLC photometry (MIT QLP news, 2025-09-03; Han & Brandt 2023). Consider citing it next to Huang et al. (2020).
5. **NOTE: secondary-eclipse scan.** A scan without a two-sided baseline requirement finds a spurious "4.2 sigma" dip at phase 0.368, caused by the 0.18-d segment BTJD 4003.05-4003.23. Any scan quoted later should use two-sided baselines.

## Recommended wording (Sec. 2)

"A joint emcee fit of 22 transits (circular orbit; Claret 2017 limb-darkening priors; Mann-relation density prior; one Sector 98 transit at a scattered-light segment edge excluded) gives Rp/R* = 0.0203 +- 0.0010 (including the spread across detrending choices), b = 0.31 (+0.12/-0.17), T14 = 1.89 +- 0.06 h, and Rp = 1.38 +- 0.08 Re (1.37-1.39 Re with free limb darkening). With the density free, rho* <~ 3.6 rho_sun (97.7%). Per-sector depths agree within errors, odd and even transits differ by < 1 sigma, the flux at phase 0.5 is +20 +- 32 ppm (i.e. no dip), and a V-shaped transit is disfavored (delta-chi2 ~ 32)."

Abstract: "... corresponds to a ~1.4 Re planet ..."; "The 11.13 d signal (~2.2 Re) ...".

## Files

- Scripts: `v00_inspect.py`, `v00b_segments.py`, `v00c_epochs.py`, `v01_claret.py`, `vcommon.py`, `v02_prep.py`, `v03_fit.py`, `v04_tests.py`, `v05_vshape.py`, `v06_epoch_m1.py`, `v07_secscan.py`, `v07b_dip_check.py`, `v08_bic_order.py`, `v09_summary.py`, `run_fits*.sh`.
- Outputs (`out/`): `summary_verify_fits.csv`, `fit_*.json`, `chain_*.npz`, `tests_ms_template.json`, `vshape_*.json`, `s98_epoch_detrend_diag.json`, `bic_baseline_order.json`, `secscan_T14.{json,txt}`, `dip_check_phase0.368.json`, `lc_*.npz`, `prep_*.json`.
- Claret: `claret_tableab.ecsv`, `claret_table25.ecsv`, `claret_provenance_verify.json`.
- Logs: `logs/`.
