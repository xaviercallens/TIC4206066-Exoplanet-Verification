# METHOD_AUDIT — RNAAS manuscript, TIC 4206066

What was done for every number in `manuscript.tex`, by whom (this directory = independent re-analysis written from
scratch, `analysis/a1`–`a6`; "audit" = `../tic4206066_PUBLICATION_AUDIT/`, a separate earlier re-analysis; "deep" =
`../deep_followup_tic4206066/`), which primary source defines the method, and how the routes compare. Nothing in the
frozen search pipeline (`../../pipeline/`, `../../config/`, `../../candidates/`) was modified; it was only read.

Environment for `analysis/a1`–`a6`: `uv` venv, Python 3.12, numpy 2.x, astropy, batman-package 2.5.1, emcee, wotan,
tess-prf 0.1.3, lightkurve 2.5.1 (bitmask constants only); TRICERATOPS 1.1.0 in the audit's Python 3.11 venv
(`a4`). Full freeze: `tools/env_freeze.txt`.

## 1. Data (raw inputs only)
* Light curves: TESS-SPOC S6 PDCSAP (30 min), SPOC S32 PDCSAP (2 min), QLP S98 SYS_RM_FLUX (200 s) — the deep
  follow-up's MAST downloads, which the audit re-downloaded and found byte-identical (16/16 sha256). Quality: lightkurve
  DEFAULT_BITMASK (17087) for SPOC products; the frozen search mask (1610619071 = SPOC bits + QLP 2^29 and 2^30) for
  QLP S98, i.e. the Earth scattered-light intervals at the start of each S98 orbit (HEASARC S98 notes) are removed.
* Pixels: TESScut 21×21 cubes S6/S32/S98 (audit's downloads, identical to the deep follow-up's), the SPOC S32 2-min TPF.
* PRF: MAST `start_s0004` PRF FITS files (cam1_ccd1, cam1_ccd2), read with TESS_PRF 0.1.3, which un-interleaves the
  117×117 files into 9×9 sub-pixel PRFs. Same lesson as the TIC 158955499 audit (its superseded code is not used).
* Gaia DR3 census within 180″ (audit's TAP download); proper motions propagated to each sector epoch (a2 `scene`).
* Time system: all products declare TIMESYS = TDB, TIMEREF = SOLARSYSTEM, BJDREFI = 2457000 (`a1_results.json`
  `time_system_headers`). TESScut times carry the CCD-centre barycentric correction; the target-specific shifts
  (−11.8, −3.3, +0.5 s for S6, S32, S98) were computed by the audit (`s06`) and are applied in a1/a2; they are far below
  the timing errors (≥ 4 min).
* Own aperture photometry (a1 `own_phot`): r = 2.5 px (S6) / 2.0 px (S32, S98) around the brightest central pixel,
  background = median of the faint pixels > 5 px away. Used only as a second photometry (never as the primary).

## 2. Detrending and templates
* wotan biweight, window 0.75 d, both signals masked (±1.5 T14) in the trend fit, upward 4-MAD clipping (a1
  `detrend`). The audit showed depths change < 5 % for windows 0.5–1.5 d, cosine, or no detrending (`s17`).
* Transit template: batman, quadratic LD, circular orbit, supersampled to the cadence; shape fixed to the audit's
  Mann-prior joint fits (A: k 0.02116, b 0.35, ρ 2.507; B: k 0.03305, b 0.871, ρ 2.507; also B with the free-density
  shape k 0.0296, b 0.49, ρ 13.35). Per-event depths and times do not depend on the template at the quoted precision
  (B depths with the two templates differ by < 1 %).

## 3. Transit times, period, T0 (Step 1–2) — `a1`
* Per event: window ±0.3 d, other signal removed (±1 T14), linear baseline + template, mid-time on a 0.25-min grid.
* Errors: (i) white (Δχ² = 1 with the residual variance), (ii) cyclic residual permutation (200 shifts), (iii)
  injection–recovery: the template with the sector's median depth injected at 150 random times of the same light
  curve (away from both signals) and fitted identically → sd of the recovered time and depth. (iii) is adopted: it
  contains the real red noise and the finite-cadence sampling.
* Linear ephemeris: weighted least squares on the 20 events with coverage ≥ 0.7 and white S/N ≥ 2 (not at the scan
  edge). Result A: P = 3.1827856 ± 0.0000058 d, T0 = 3998.7375 ± 0.0023 BTJD, χ² = 12.6/18 (injection errors); with
  white errors χ² = 81/18, with permutation errors 20.8/18. So the earlier "χ² = 120/21" was a white-noise error
  artefact; no timing anomaly remains (Step 2, outcome A).
* Cross-checks: own photometry 3.1827828 ± 0.0000060 d (23 events); audit joint MCMC 3.1827853 ± 0.0000027 d (formal);
  audit per-event 3.1827793 ± 0.0000031 d; blind combined BLS 3.18279 d. All within 1σ of the adopted value.
* Cycle count: the S98-only linear ephemeris has σ_P = 0.00044 d; the alias spacing to S32 (569 cycles) is 0.0056 d =
  13σ, so the cycle count is unique (χ² for ±1 cycle offsets > 13,000 vs 12.6).
* Older value 3.182945 d = the frozen search's single-sector BLS grid value; it lies 0.4σ (S98-only precision) from
  the final period. The deep follow-up's 3.1827863 ± 0.0000028 d came from a joint fit in which B's transits were not
  masked (audit bug A1); the difference (1e-6 d) is far below the realistic error.

## 4. Signal B depths and sector comparison (Step 4) — `a1`, `a5`
* Same per-event method and injection errors in every sector and photometry, two templates.
* Pipeline set: S6 838 ± 114, S32 1003 ± 206, S98 1099 ± 93 ppm; S98 − archival = 1.6σ; Δχ² (common vs per-sector
  depth) = 3.1 for 2 extra parameters; ΔBIC favours the common depth (996 ± 68 ppm). Own photometry: 1.2σ (the own
  S32 event −164 is corrupted by a centroid excursion, fit at the scan edge, and is excluded). The audit's 2.2–2.7σ
  used white errors scaled by √χ²_red × β; with injection errors the tension is weaker. Dilution: CROWDSAP 0.993–0.997
  (< 1 %). The difference is reported as not significant.

## 5. Pixel localisation (Step 3) — `a2`, `a5`
* Cadences: DEFAULT_BITMASK removed; S98 straylight (2048) removed and, for the clean variant (adopted), the QLP-2^29
  scattered-light cadences removed by matching QLP cadence times; other signal removed (±1 T14).
* Scene: PRF model of all Gaia sources (T < 19) with a free global shift fitted to the median image → target stamp
  position. Difference image per transit (linear per-pixel baseline 0.75–2.5 T14 from mid-transit minus the in-transit
  mean, |dt| < 0.4 T14), averaged. Noise map: per-pixel sd of 300 difference images at random fake epochs.
* Fit: D = A·PRF(x, y) + c on pixels within 4 px; free position (0.1-px grid, parabolic refinement); fixed target and
  neighbour positions. Calibration: 300 model-free injections (real pixels within 3 px × (1 − depth) in fake transits)
  → mean and covariance of the recovered offset → z of the real offset; PRF-shaped injections on each neighbour able to
  host the signal (aperture-flux loss matched).
* Result A: S32 TPF z = (+0.77, −0.19); S32 FFI (+0.98, +0.07); S98 clean (+1.05, −1.41); combined S32 TPF + S98 clean
  (+7.9, −8.2)″ ± (5.3, 6.8)″, χ² = 5.3 for 2 dof (p ≈ 0.07) — consistent with the target at 2σ. S6: (+6.6, −3.4)″ at
  z = (+3.0, −0.6); every S6 transit is 0.4–3.5 h from a momentum dump (quality bit 32), so S6 is excluded from the
  headline (the audit reached the same conclusion). S98 with the scattered-light cadences kept: z = (+2.0, −1.8) — the
  earlier "19″" S98 offset was this effect. Neighbours: the nearest star able to produce the depth (43.5″) is 7.6σ from
  the combined position; all capable stars ≥ 7.1σ; PRF-shaped injections on them never reproduce the data
  (`a2_localisation_*.json`, `neighbours`).
* Result B: S6 z = (−0.55, +0.41), S32 TPF (−0.99, +1.22), S98 clean (+0.76, −0.12); combined S32 + S98 (+1.4, +0.1)″ ±
  (3.3, 4.4)″; all three (+0.2, +0.6)″ ± (2.3, 3.1)″.
* Comparison with the audit (`s12`, independent code): A S32 2-min z = (0.75, −0.10), S98 clean (0.07, −1.45), combined
  (+2.5, −7.5)″ ± 5.7″; B combined (+1.0, +1.4)″ ± (2.3, 3.0)″. Same conclusions; the S98 x-offset differs by 1σ
  (10 vs 9 clean epochs, different noise realisations).
* Superseded: the deep follow-up's shifted-OOT-image fits with bootstrap errors (S32 "10.5″ at 3.8σ", S98 "19″") —
  the errors were ~3× too small (audit A2) and the S98 set included scattered-light cadences (audit A5).

## 6. Signal B significance — `a3`
* Cross-sector link: statistic = maximum 1.6-h box S/N at the Sector 98 ephemeris over a ±0.004-d period scan
  (1601 steps), A masked in the detrending and removed. Null: 10,000 trials with the period drawn from 9–13.5 d
  (excluding ±0.1 d around P_B) and a random epoch, the scan width scaled to the same time window, B's own archival
  transits removed from the null data. Results: S6+S32 obs 10.03, FAP 2.0×10⁻³ (19/10,000; 9,528 distinct windows);
  S6 alone obs 7.08, 0/10,000; S32 alone 0.032; own photometry S6+S32 1.2×10⁻³, S6 0/10,000. The earlier "exceeds all
  300 epoch-shifted trials" (deep follow-up, previous draft) is withdrawn: those trials shared a handful of data windows
  (audit A14) and give no calibrated FAP.
* Rank test (distribution-free): in S6 the two predicted dips are the two strongest 1.6-h boxes outside A (5.8 and 4.1
  vs next 3.0 in TESS-SPOC; 5.2 and 6.3 vs 3.5 in own photometry); p = K·3/C(N, 2) = 1.6×10⁻³ (N = 325 independent
  boxes, K = 27 positions allowed by the scan window, 3 spacings).
* Sector 98 alone: astropy BLS (1–20 d, likelihood objective) on the A-removed light curve; null = B removed and the
  flux circularly shifted within each continuous segment (> 0.5 d gaps), 300 trials, maximum SDE recorded. QLP frozen
  mask: B SDE 6.07 (the global maximum), null median 5.47, FAP 0.35. In my own TESScut aperture photometry (straylight
  removed) B is not the strongest peak (SDE 4.25; a 14.3-d scattered-light peak reaches 8.87), FAP ≈ 1; the audit's own
  photometry gave 0.008. The S98-only significance is therefore photometry-dependent (0.008–1), and the Note quotes only
  the QLP result. The discovery MES 11.0 is LEO-vetter's pink-noise S/N, not a calibrated FAP.
* Audit values reproduced: link FAP 1.2×10⁻³ (S6+S32), S6 0/18,652, S32 0.033, rank p 1.6×10⁻³, S98 QLP FAP 0.33. Not
  reproduced: the audit's own-photometry S98 FAP 0.008 (mine ≈ 1, different aperture photometry).

## 7. Global shape parameters, EB tests, hidden companions (not recomputed here; sources)
* Rp/R★, b, T14, Rp, ρ★ limit: audit joint emcee fits (`../tic4206066_PUBLICATION_AUDIT/out/fits/`, Mann-prior and
  free-density chains; dynesty cross-check with free LD). Consistency check here: my weighted mean of the per-sector
  mid-transit depths (551 ± 34 ppm, pipeline) agrees with k² × limb-darkening factor (≈ 503 ppm) at 1.4σ.
* Odd/even (536 ± 40 vs 466 ± 43), secondary (−20 ± 32 ppm; empirical σ from random same-width boxes), V-shape Δχ² = 32:
  audit `s10`, `s10b`. Gaia RV limit (K < 3.2 km/s → M sin i < 0.016 M☉ at 3.18 d): audit `s14` from the Gaia DR3
  radial_velocity_error formula; recomputed here in `audit/number_audit.py` (0.016 M☉).
* Overluminosity: Way et al. (2026) catalogue row (O = +0.003 mag; local robust σ 0.068 for 5,205 comparable stars);
  2σ companion limit 0.38 M☉ (companion_check `overluminosity.py`, Pecaut & Mamajek dwarf table). The 0.54 M☉ host
  requirement: transit density < 3.5 ρ☉ (97.7 %) → host earlier than M1V (Pecaut & Mamajek: M0.5V 0.54 M☉, 3.35 ρ☉;
  M1V 0.50 M☉, 3.98 ρ☉), checked here in `audit/number_audit.py`.
* SED: astroARIADNE 1.5.0 (`../companion_check_tic4206066/sed_out/best_fit_average.dat`): R 0.6206 ± 0.015, Teff
  4012 ± 43 K; per-band residuals ≤ 0.05 mag, no IR excess. Not re-run here (grid download); values read from the files.

## 8. FPP — `a4` and audit `s13`
* TRICERATOPS 1.1.0, N = 10⁶ draws per run, `dilution_corrected=True`, exptime = 6 min (bin width), apertures = SPOC
  optimal apertures (S6, S32) and the QLP 3×3 (S98), stellar parameters Mann M/R, Teff 4030 K, ZP-corrected parallax.
  Default: companions unresolved within 2.2″ (no contrast curve).
* Folds: (1) audit fold (joint-fit ephemeris): A 0.0312 ± 0.0031, B 0.179 ± 0.011 (20 runs); (2) earlier fold (deep
  follow-up, S98 BLS epoch for B): A 0.0304, B 0.146 (reproduced by the audit to 7×10⁻⁵ per seed); (3) this work's
  fold (a1 linear ephemeris, own detrending): A 0.034 ± 0.005, B 0.195 ± 0.021 (10 runs, seeds 2000–2009).
  → quoted: A 0.031 ± 0.003 (0.031–0.034 across folds); B 0.13–0.20 (with the TIC stellar parameters 0.128).
* MOLUSC variant (audit): A 0.0047 ± 0.0027, B 0.064 ± 0.020 (age 5 Gyr; 1 Gyr: 0.0054, 0.059). Code reading
  (`marginal_likelihoods.py`): with `molusc_file` only companions with a(1−e) > 10 au are kept and the companion prior
  is set to zero in log (kept fraction 13 %). Stated in the manuscript.
* NFPP ≈ 10⁻¹² ignores the T = 13.06 star at 83.8″ (Gaussian PSF); that star is excluded by the pixel localisation.

## 9. Sector 110 — `a5`/`audit/number_audit.py`
HEASARC sector110_summary (2026-10-31 → 11-26, camera 1 boresight 75.88, −2.69); tess-point (audit `s16`) camera 1 /
CCD 2. Transit counts from the a1 ephemerides: A 8 (to Nov 26 00:00 UTC) or 9 (to Nov 27 00:00), B 3.

## 10. Novelty
Sealed priority re-check (`../priority_recheck_tic4206066_SEALED/`, cutoff 2026-09-25T04:46Z, 36 source groups incl.
T16, RAVEN, DTARPS), my own scripted search (`../rnaas_tic4206066_draft/novelty_check/`, 05:01–05:16Z, incl. ADS
full-text), and the refresh of the fast-changing sources in `novelty_refresh/` (ExoFOP TOI/CTOI, target page, NASA toi,
exo.MAST TCEs) shortly before this draft was finalised.
