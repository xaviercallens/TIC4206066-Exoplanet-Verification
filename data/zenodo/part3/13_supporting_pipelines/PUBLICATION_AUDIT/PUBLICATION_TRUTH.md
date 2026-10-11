# PUBLICATION TRUTH — TIC 4206066 (StKM 1-561)

**Authoritative source for manuscript writing.** It supersedes, for every number and statement, the earlier documents:
`deep_followup_tic4206066/deep_followup_report.md`, `…/fpp/FPP_REPORT.md`, `…/archival_checks/SUMMARY.md`,
`companion_check_tic4206066/SUMMARY.md`, `rnaas_tic4206066_draft/tex/manuscript.tex` and the frozen-search bundles. Those
files are kept unchanged as history. Every number here was recomputed in this audit from the lowest-level public data
(`scripts/`, `out/`). Machine-readable values: `FINAL_PARAMETERS.json`. Per-claim status and allowed wording:
`CLAIM_AUDIT.csv`. Methods and their primary sources: `METHOD_AUDIT.md`. Bugs: `BUG_REGRESSION_TESTS.md`. Open problems:
`UNRESOLVED_ISSUES.md`. Audit date: 2026-09-25 (UTC). Nothing was published or submitted; Jev and Opus were not run.

Conventions: times BTJD = BJD_TDB − 2457000. Depths in ppm of the undiluted target flux. "±" is the 68 % interval.
σ-values are stated with the error model used (white, red-noise-inflated or empirical).

---------------------------------------------------------------------------------------------------------------------
## 1. Status in one paragraph each

**Signal A — P = 3.1827853 d. Strong planet candidate; not validated.** Recovered independently in each of the three
sectors (2018, 2020, 2025), in four independently processed photometries of the same pixels, with 23 usable transits.
The transit is flat-bottomed (V-shape disfavoured, Δχ² = 32 with white errors rescaled to χ²_red = 1). Odd and even
depths agree (per-event fits ≤ 1.2σ; box comparisons +2.2σ and −2.1σ in two photometries, i.e. opposite signs). There is
no secondary eclipse (< 95 ppm, 3σ). The
transit-shape density is consistent with the star's (< 3.5 ρ☉ at 97.7 %). The source lies on the target within ~6″ (1σ)
in the data that can be calibrated (S32 2-min, S98 clean). Every Gaia star able to mimic it (43–84″) is excluded by the
pixels. An EB on the target is excluded by Gaia RVs. A circular-orbit planet on a bound companion is excluded by density
plus Gaia-XP overluminosity. TRICERATOPS 1.1.0 gives FPP = 0.031 ± 0.003 under its default assumptions: a "likely
planet" by the Giacalone et al. (2021) criteria, not "validated". Rp = 1.44 ± 0.06 R⊕.

**Signal B — P = 11.132745 d. Low-significance candidate (overall about 3σ); Sector 110 must confirm it.** Found in
S98 by the discovery pipeline (three clean transits, plus a fourth at the predicted time in a scattered-light-affected
interval). In S98 alone it is not secure. A calibrated BLS test, with the null built from segment-scrambled S98 data,
gives FAP 0.33 in the discovery (frozen-mask) data, 0.008 in the own photometry and ~1 with the QLP default mask, where
scattered light dominates. The S98 data contain systematics as strong as B. The independent support comes from
Sector 6: the two strongest dips of the whole sector (outside A's transits) fall at the two times predicted from S98, in
both photometries (distribution-free p = 1.6 × 10⁻³). S32 alone gives FAP 0.03. Combined as independent tests (Fisher),
the overall false-alarm probability is 5 × 10⁻⁵ to 1 × 10⁻² (3.9σ to 2.3σ), depending on which S98 product is used;
this is indicative. The earlier "FAP < 5 × 10⁻⁵" for the archival link rested on overlapping null trials and is
withdrawn (bug A14/C7). The S98 depth exceeds the archival depth at 2.2–2.7σ (unexplained). The transit shape does not
constrain the stellar density, so the radius depends on the assumed impact parameter: 2.25 ± 0.12 R⊕ if B orbits the
target, 2.0–2.1 R⊕ otherwise. The difference images are consistent with the target. TRICERATOPS FPP = 0.13–0.18 under
default assumptions (0.06 with MOLUSC); the dominant scenario, a planet around an unresolved bound companion, is not
excluded. The FPP assumes the signal is astrophysical and does not include the false-alarm probability above.

---------------------------------------------------------------------------------------------------------------------
## 2. The star (adopted values)

| quantity | adopted | route / source | alternatives (not adopted) |
|---|---|---|---|
| Gaia DR3 id | 3220388198192519424 | | TIC 4206066; TYC 4753-1591-1; 2MASS J05232636-0119380 |
| RA, Dec (ICRS, ep. 2016.0) | 80.85982666, −1.32696401 | Gaia DR3 | |
| μα*, μδ | −3.406, +70.806 mas/yr | Gaia DR3 | |
| ϖ | 28.108 ± 0.019 mas | DR3 28.0750 ± 0.0176 minus Lindegren+2021 ZP (−0.0332) | uncorrected 28.075 |
| distance | 35.58 ± 0.03 pc | 1/ϖ (corrected); Bailer-Jones+2021 35.58 | 35.62 (uncorrected) |
| G, BP−RP, T | 10.513, 1.739, 9.682 | Gaia DR3, TIC 8.2 (Stassun eq. gives 9.685) | |
| Ks, M_Ks | 7.694 ± 0.029, 4.938 ± 0.029 | 2MASS | |
| **R★** | **0.623 ± 0.019 R☉** | Mann+2015 R(M_Ks) (2016 erratum checked) | SED 0.621 ± 0.015; [Fe/H]-form 0.613–0.621; TIC 0.648 ± 0.065; Gaia DR2 0.618 |
| **M★** | **0.606 ± 0.015 M☉** | Mann+2019 n = 5, total relation error 2.4 % | official posterior code 0.611 ± 0.015; TIC 0.63 ± 0.08; Mann+2015 quartic 0.645 |
| **ρ★** | **2.51 +0.25/−0.22 ρ☉** (3.53 g cm⁻³) | M/R³ (MC) | TIC 2.31 |
| log g | 4.632 ± 0.028 | | |
| **Teff** | **4030 ± 90 K** | mean of TIC 4010, Mann V−J 4048, SED 4012, Gaia DR2 4069 | GSP-Spec 4182 (flagged, rejected) |
| L★ | 0.092 ± 0.010 L☉ | R²(Teff/5772)⁴ | SED 0.090, TIC 0.098, BC_K 0.107, BC_V 0.101 |
| spectral type | K7V–M0V | literature + M_Ks | |
| [Fe/H] | poorly constrained (−0.3 … 0) | GSP-Spec −0.27 (extrapolation flag), SED −0.08 ± 0.18 | |
| RV | 10.49 ± 0.22 km/s | Gaia DR3 (33 transits) | not RV-variable by the Gaia criterion (gof 3.0 < 4; p = 4e−5) |
| RUWE | 1.214 | Gaia DR3; 64th percentile of 828 comparable nearby stars | DR2 RUWE 1.04 |
| NSS / variability | none | Gaia DR3 | |
| rotation | not detected | TESS SAP, ASAS-SN | catalogue ~8 d not confirmed |
| X-ray | 1eRASS J052326.6-011938, 4.8″, DET_LIKE 8.0, L_X ≈ 5 × 10²⁷ erg s⁻¹, log L_X/L_bol ≈ −4.8 to −4.9 | eRASS1 | marginal single-survey detection |
| LD (TESS, quadratic) | prior centres u1 = 0.33, u2 = 0.30 (σ 0.1) | Claret 2017 table15 at 4030 K, 4.63 | |

---------------------------------------------------------------------------------------------------------------------
## 3. Data

* TESS sectors: **S6** (2018-12-15 → 2019-01-06; camera 1, CCD 1; 30-min FFI); **S32** (2020-11-20 → 12-16; camera 1,
  CCD 2; 2-min SPOC target + 10-min FFI); **S98** (2025-11-08 → 2026-01-05, four orbits; camera 1, CCD 2; 200-s FFI).
  No other sector has observed the star.
* Primary light curves: TESS-SPOC S6 PDCSAP; SPOC S32 PDCSAP; QLP S98 SYS_RM_FLUX with the frozen mask (removes QLP
  bits 2^29 and 2^30, i.e. the camera-1 Earth scattered light at the start of each S98 orbit documented by HEASARC).
  Independent photometry: own aperture photometry from TESScut in all three sectors, same method, plus the S32 2-min TPF.
* Detrending: biweight 0.75 d with both signals masked. The results do not depend on this choice (0.5–1.5 d biweight,
  cosine, none: depths within 5 %).
* Timing: all products are BJD_TDB. After the target-specific shift, TESScut, TESS-SPOC and QLP agree within 1–5 s.
  QLP S98 is not affected by the QLP barycentric errors reported for S74–79 and S99–104.

---------------------------------------------------------------------------------------------------------------------
## 4. Signal A

### 4.1 Ephemeris
| | value | route |
|---|---|---|
| **P** | **3.1827853 ± 0.0000027 d** | joint fit (batman + emcee), pipeline set |
| **T0** | **2460998.73679 ± 0.00112 BJD_TDB** (3998.73679 BTJD) | same; corr(T0, P) = 0.66 |
| cross-check | 3.1827853 +0.0000028/−0.0000027 d | dynesty, free LD, biweight 1.0 d |
| cross-check | 3.1827793 ± 0.0000031 d (χ² 47.5/19); 3.1827827 ± 0.0000032 d (own photometry) | per-event linear ephemeris |
| cross-check | 3.18279 d | blind BLS of all three sectors |

The cycle counts S98 → S32 (569) and S98 → S6 (792) are unique: the alias spacing is 13–18σ of the S98-only period
precision, and the blind combined BLS picks the same alias.

Blind recovery in single sectors (`s07`). The strongest BLS peak of S6 alone lies 0.02 period-resolution elements and
0.02 h from the ephemeris prediction. In S32 the peak (second, after one systematics peak) lies 0.04 elements and 0.03 h
away. For noise, the chance of landing that close in period (~2 of ~860 independent periods) and phase (~2.5 %) is
~10⁻⁴ per sector. A is therefore not a statistical fluke. Its astrophysical interpretation is treated in §4.4–4.6.

### 4.2 Shape and derived quantities (ρ★ prior = catalogue density)
| parameter | value |
|---|---|
| Rp/R★ | 0.0212 ± 0.0007 |
| depth | k² = 448 ± 30 ppm; mid-transit (limb-darkened) ≈ 503 ± 29 ppm; mean in 1.87 h box ≈ 470 ± 32 ppm |
| b | 0.35 +0.10/−0.14 (prior-driven) |
| a/R★ | 12.47 +0.37/−0.40 |
| i | 88.41 +0.66/−0.50° |
| T14 | 1.874 ± 0.048 h |
| **Rp** | **1.44 ± 0.06 R⊕** (free ρ★: 1.46 +0.09/−0.07; free LD + free ρ★: 1.50 +0.09/−0.08; own photometry: 1.41 ± 0.07) |
| a | 0.0358 ± 0.0003 AU |
| S | 72 +8/−7 S⊕ |
| Teq (A = 0, full redistribution) | 810 ± 22 K |
| K (per Earth mass) | 0.607 ± 0.010 m s⁻¹ M⊕⁻¹ |

Free-density fit: ρ★ = 2.75 +0.46/−0.97 ρ☉ (uniform prior) or 2.37 +0.72/−1.33 ρ☉ (ln-uniform prior). The upper side
is set by the data (< 3.5 ρ☉ at 97.7 % for both priors); the lower side by the prior. Consistent with 2.51 ± 0.23.

### 4.3 Events (23 usable: 6 in S6, 8 in S32, 9 in S98 clean)
Per-sector depths (template fits, red-noise errors): S6 521 ± 53, S32 476 ± 52, S98 511 ± 47 ppm. All six detrendings give
483–558 ppm per sector. Median single-event S/N 3.5. No event dominates: the largest single share of (S/N)² is 24–39 %
per sector, and leaving out any one event keeps S/N ≥ 7.2 per sector. Other photometries: TESS-SPOC S32 474 ± 50, SPOC SAP
S32 473 ± 53, own TESScut S32 480 ± 52, S98 468 ± 41 (15 events). S6 depends on the pipeline (TESS-SPOC 521 ± 53, own
653 ± 56, eleanor 857 ± 99, TGLC aperture 982 ± 101). Only the SPOC/TESS-SPOC value is used. Events inside S98
scattered-light intervals (QLP bit 2^29, SPOC straylight) give outlier depths in both directions and are excluded.

### 4.4 Eclipsing-binary tests (empirical red-noise errors)
* Box S/N (all sectors) 14.8 (empirical) vs 20.3 (white). The empirical σ comes from same-width boxes at random phases
  of the fold. These sample only ~37 independent phases, so σ itself is uncertain by ~12 %.
* Secondary at phase 0.5: −20 ± 32 ppm; 3σ upper limit 95 ppm.
* Phase scan: largest \|z\| = 2.2 among ~37 independent phases (random-box p = 0.01 per phase, not significant). In the
  own photometry the largest feature (\|z\| = 3.8 at phase 0.27) comes from one 4000-ppm, < 2-h dip at BTJD 2185.4 in
  the own S32 TESScut light curve. It is absent from SPOC PDCSAP and vanishes when that stretch is left out (`s10c`).
* P/3 phases (1/3, 2/3): −41 ± 32 and −7 ± 32 ppm → no dips. n mod 3 groups: consistent in the pipeline photometry
  (χ² 1.7/2). In the own photometry χ² = 10.6/2 (group 0: 269 ± 63 vs 519 ± 59 and 505 ± 57 ppm). All three groups show
  the dip at ≥ 4σ, so P = 3 P_A is excluded either way.
* Odd/even, per-event template fits: 536 ± 40 vs 466 ± 43 ppm (1.2σ, pipeline); 530 vs 530 ppm on matched events in the
  own photometry (0.0σ). Box comparison with empirical errors: odd deeper by 2.2σ (pipeline), even deeper by 2.1σ (own
  photometry). The sign flips between two photometries of the same pixels, so there is no consistent odd/even difference.
* Shape: a pure V (grazing EB) is disfavoured by Δχ² = 31.8 (pipeline) and 31.4 (own); ingress fraction < 0.18 (2σ).
  These Δχ² use white errors rescaled to χ²_red = 1. Red noise would lower them, but not below the ~4σ level for a
  red-noise factor up to ~1.4.

### 4.5 Where the dip comes from
* PRF localisation (official TESS PRF, injection-calibrated): S32 2-min offset z = (0.75, −0.10); S98 clean intervals
  z = (0.07, −1.45). Combined (S32 + S98 clean): **(+2.5, −7.5)″ ± 5.7″ from the target** (χ² = 2.1 for 2 dof).
* S6 is not usable for A: every S6 transit lies 0.4–3.5 h from a momentum dump. With S6 included, (+6.3, −3.9)″ ± (2.7,
  3.1)″ (χ² = 7.3, p = 0.03).
* Flux-weighted centroids (SPOC 2-min MOM_CENTR, empirical null): z = −1.7 and +0.1; pointing (POS_CORR) z = −1.1 and
  +0.3. No significant shift.
* Stars that could produce A: with full flux inclusion, the nearest are Gaia sources at 43.5″ (T 17.1, needs 45 %),
  50.3″, 50.9″, 54.0″ and 56.0″ (T 16.3–17.3). With the real PRF in the real apertures, the T = 13.06 star at 83.8″ could
  also produce it in the S32 SPOC aperture (28 % eclipse). All are ≥ 8σ from the localised position (the 83.8″ star ≥ 13σ).
  PRF-shaped injections on the six nearest capable stars never reproduce the real difference images (0 of 150 each, in
  S6, S32 and S98).

### 4.6 Hidden companions (A)
* EB on the target: excluded. Gaia RV scatter (robust amplitude 3.2 km/s) limits companions to M sin i < 0.016 M☉ at
  3.18 d and < 0.02 M☉ at 6.37 d.
* Planet on an unresolved bound companion (STP, TRICERATOPS's dominant scenario): on a circular orbit the transit
  density (< 3.5 ρ☉) requires a host of ≥ 0.54 M☉ (≥ M0.5V). Such a companion would raise the system's G-band
  luminosity by ≥ 0.46 mag; Gaia XP gives an overluminosity of +0.003 mag (Way et al. 2026; robust scatter 0.07 mag).
  → **circular-orbit STP excluded**. Eccentric configurations (e ≳ 0.6 near apoastron) remain.
* Eclipsing binaries diluted by the target (bound or background; all TRICERATOPS EB scenarios together): ≤ 2 × 10⁻⁶ in
  every variant.

---------------------------------------------------------------------------------------------------------------------
## 5. Signal B

### 5.1 Ephemeris
| | value | route |
|---|---|---|
| **P** | **11.132745 ± 0.000013 d** | joint fit, pipeline set |
| **T0** | **2461000.86587 ± 0.0016 BJD_TDB** (4000.86587 BTJD) | corr(T0, P) = 0.57 |
| cross-checks | 11.132743 ± 0.000013 (dynesty); 11.1327478 ± 0.0000133 (per-event, own); 11.13273 (blind combined BLS) | |
| S98 only | 11.13387 ± 0.0011 d (LEO-vetter trapezoid) | the cross-sector cycle count is unique (45–62σ) *if* the archival dips are B |

### 5.2 Does B exist outside S98? (null test corrected)
Statistic, as in the earlier test: the maximum box S/N (1.6 h) of the archival data over a ±0.004 d scan around the S98
period, at the S98 epoch. The S/N at the exact S98 period is ~0.7 because the S98 period is uncertain by ±0.0011 d (±6 h
at S6), so the scan is needed. Detrending masked only A, so the real B phase gets no advantage over null phases.

**The earlier null (300 trials) and this audit's first version (20,000 trials) are not valid.** They shifted only the
epoch at the fixed S98 period. The scan window spans ±0.9 d in S6, so trials whose epochs differ by less than that
re-use the same data. In S6 the 20,000 trials contain ~5–16 distinct windows, and their 95th, 99th percentile and
maximum are all the same feature (3.01). The correct null randomises the period as well (`s11b`).

| test (`s11b`, `s11c`) | observed | null | FAP |
|---|---|---|---|
| S6 + S32, pipeline; period (9–13.5 d) and epoch random, 20,000 trials | 10.38 | median 3.23, p99 7.32, p99.9 10.68, max 10.81 (18,411 distinct windows) | 1.2 × 10⁻³ (24 trials) |
| **S6 only, pipeline**, same null | **7.56** | median 2.52, p99 3.98, max 5.26 (2,623 distinct windows) | 0 of 18,652 |
| S6 only, own photometry, same null | 8.79 | median 2.60, max 4.81 | 0 of 18,741 |
| S32 only, pipeline, same null | 7.98 | p99 10.79 (an end-of-sector artefact of S/N 10.4) | 0.033 |
| S6 + S32, own photometry, same null | 10.08 | p99 19.65 (own S32 artefact at BTJD 2185.4) | 0.12 |
| **Rank test, S6** (distribution-free) | the two predicted B dips are the two strongest 1.6-h dips of the sector outside A: 6.2 and 4.3 vs next 3.2 (pipeline); 5.6 and 6.8 vs 3.7 (own) | 325 independent box positions; 27 allowed positions × 3 spacings | **1.6 × 10⁻³** |
| Gaussian noise model, S6 | 7.56 | per-position S/N sd 1.02; fraction > 3σ 0.22 % (Gaussian 0.13 %); 27 positions per window | ~2 × 10⁻¹² |

What can be claimed: in S6 the two dips at the S98-predicted times are the two strongest features of the whole sector
(outside A's transits), in both photometries. That is about 3σ without any noise model. Much smaller FAPs need the
Gaussian model, which the S6 per-position scatter supports but cannot prove in the tail. S32 does not add significance.
Its data contain artefacts as strong as B, and its two B dips rank 3rd and 4th among S32 features.

### 5.3 S98 transits
Clean: epochs 0 (1102 ± 132 ppm), 1 (1245 ± 139) and 4 (partial in QLP, 939 ± 242; complete in TESScut, 917 ± 196).
Epoch 2 (4023.131) lies in a QLP-flagged (2^29) but SPOC-unflagged interval: 1198 ± 198 (TESScut) / 1441 ± 204 (QLP)
ppm at the predicted time (O−C −3.0 ± 3.0 min). It is an out-of-sample confirmation of the S98 period, but it sits on a
background ramp. Epochs −1 and 3 are inside SPOC straylight intervals and are not used.

How strong is S98 alone? In the discovery (frozen-mask) data the two clean B events have box S/N 7.3 and 7.9, and one
non-B feature is stronger (BTJD 4017.68, S/N 8.7; `s11c`). An independent blind iterative BLS (`s07`) finds B with the
QLP default mask (SDE 9.8, but only after two scattered-light peaks have been found and masked) and in the own
photometry (SDE 7.4, after one other peak). With the frozen mask it is only the second peak (SDE 5.5, behind a 5.60-d
systematics peak). The S98 evidence therefore rests on the repetition of 3–4 events with consistent depth, duration and
shape.

Calibrated S98-only test (`s07b`, `s07c`). BLS as in `s07`, P = 1–20 d, A masked and removed. Null: B's transits removed
and every continuous segment circularly shifted (500 trials); the null records the maximum SDE over all periods.

| S98 product | B's SDE (period) | null max-SDE median / p95 / p99 | FAP |
|---|---|---|---|
| QLP, discovery (frozen) mask | 6.09 (11.1368 d; the strongest peak once A is removed) | 5.47 / 6.92 / 7.18 | 0.33 |
| QLP, default mask (keeps bit 2^29) | 1.48 (buried under scattered light; top peak 14.29 d) | 5.55 / 7.92 / 8.77 | 1.0 |
| own TESScut photometry, straylight removed | 10.36 (11.1329 d; strongest peak) | 5.75 / 9.51 / 10.16 | 0.008 |

Robustness for the discovery data (300 trials each): the circular-shift null gives FAP 0.22–0.30 for SDE, BLS power and
depth S/N. A 0.5-d block bootstrap gives 0.37 for SDE and 0.03–0.04 for power and depth S/N; this null breaks structure
longer than 0.5 d, so it is less conservative. Without any scrambling, the B-free S98 data already contain a peak with
SDE 6.48 and depth S/N 11.9, against B's 6.09 and 12.1. **S98 alone does not establish B.**

Indicative overall FAP (`s11d`, Fisher's method on independent data: S98 test × S6 rank test × S32 null): 4 × 10⁻³ /
1 × 10⁻³ (discovery mask; without / with S32), 1 × 10⁻² / 3 × 10⁻³ (default mask), 1.5 × 10⁻⁴ / 5 × 10⁻⁵ (own
photometry). That is 2.3σ to 3.9σ.

### 5.4 Shape, radius, depth
| parameter | with ρ★ = catalogue (if B orbits the target) | free ρ★ (ln prior) | free LD (dynesty) |
|---|---|---|---|
| Rp/R★ | 0.0331 ± 0.0015 | 0.0296 +0.0026/−0.0015 | 0.0302 +0.0037/−0.0022 |
| b | 0.87 ± 0.02 (prior-driven) | 0.49 ± 0.34 (≈ prior) | 0.53 +0.33/−0.39 |
| T14 | 1.66 ± 0.07 h | 1.55 +0.11/−0.09 h | 1.58 +0.13/−0.11 h |
| **Rp** | **2.25 ± 0.12 R⊕** | 2.01 +0.18/−0.12 R⊕ | 2.05 +0.26/−0.16 R⊕ |

ρ★ is **not constrained** by B: free-ρ median 13–17 ρ☉ depending on the prior, 2σ range ~1–27 ρ☉; Mann 2.5 ρ☉ is
inside. a = 0.0826 ± 0.0007 AU, S = 13.5 ± 1.5 S⊕, Teq = 533 ± 15 K, K = 0.400 m s⁻¹ M⊕⁻¹.
Per-sector depths (template, pipeline): S6 851 ± 101, S32 844 ± 124, S98 1139 ± 89 ppm. **S98 is deeper at 2.2–2.7σ**
in every detrending and both photometries. Odd/even: 1100 ± 103 vs 915 ± 72 (1.5σ). Phase 0.5: 61 ± 62 ppm, so
P ≠ 5.57 d. P/3 phases: no dips. V-shape disfavoured by Δχ² = 14.5 (pipeline) / 7.4 (own). The phase-scan maximum
(\|z\| = 3.6 at phase 0.24) is not significant after the look-elsewhere correction (~145 phases, global p ≈ 0.5).

### 5.5 Where the dip comes from
PRF localisation: S6 z = (−0.84, +0.34); S32 2-min (2 events) z = (−0.77, +1.5); S98 clean (3 events) z = (+2.14, −0.19).
Combined: **(+1.0, +1.4)″ ± (2.3, 3.0)″** (χ² 0.4), or (+6.1, +1.6)″ ± (3.6, 4.0)″ without S6. Consistent with the target.
In actual apertures only the 83.8″ T = 13.06 star could produce B (S32 aperture, 59 % eclipse). It is excluded by the
localisation.

### 5.6 Hidden companions (B)
The Gaia RVs exclude an EB on the target. The density argument does **not** work for B, because ρ★ is unconstrained.
A planet on a bound companion within ~1″ is the main open false-positive channel (TRICERATOPS STP ≈ 0.13–0.18). Before
that question, B's reality itself rests on about 3σ (§5.3).

---------------------------------------------------------------------------------------------------------------------
## 6. System-level facts
* P_B/P_A = 3.49774, 0.064 % inside the 7:2 commensurability. TTVs are not measurable with current timing (A rms O−C
  ~10 min).
* Companion space still open (both signals): bound stars within ~1″ (≲ 36 au) of mass ≲ 0.38–0.54 M☉ (q ≲ 0.6–0.9)
  that Gaia RUWE, Gaia imaging and XP overluminosity would not detect, and any companion on orbits wide enough to escape
  Gaia RVs (P ≳ 3 yr for M ≳ 0.1 M☉). For A such hosts would need e ≳ 0.6; for B no photometric restriction exists.
* No IR excess (Ks−W1 0.09, Ks−W3 0.13, Ks−W4 0.24 ± 0.13 mag, photospheric for K7–M0).

---------------------------------------------------------------------------------------------------------------------
## 7. False-positive probabilities (TRICERATOPS 1.1.0)

Every row: 20 runs × 10⁶ draws, seeds 1000–1019, TRICERATOPS 1.1.0, the same apertures as the earlier work. The table
is generated from `out/fpp/runs/*/summary.json` by `scripts/s13d_fpp_table.py` (copy in `out/fpp/fpp_table.md`).

| signal | variant | fold | exptime | FPP mean ± sd | range | NFPP | STP share of FPP | runs × draws |
|---|---|---|---|---|---|---|---|---|
| A | earlier configuration (Teff 4080 K, uncorrected ϖ), audit fold | S32+S98 | 3.33 min | **0.0304 ± 0.0031** | 0.0259–0.0356 | 1.5e-12 | 1.00 | 20 × 1e+06 |
| A | audit baseline: Mann M, R; Teff 4030 K; ZP-corrected ϖ | S32+S98 | 6 min | **0.0312 ± 0.0031** | 0.0265–0.0365 | 2.4e-12 | 1.00 | 20 × 1e+06 |
| A | baseline − neighbours excluded by the pixel localisation | S32+S98 | 6 min | **0.0312 ± 0.0031** | 0.0265–0.0365 | 0.0e+00 | 1.00 | 20 × 1e+06 |
| A | loc + MOLUSC companion limits (5 Gyr) | S32+S98 | 6 min | **0.0047 ± 0.0027** | 0.0009–0.0104 | 0.0e+00 | 1.00 | 20 × 1e+06 |
| A | loc + MOLUSC companion limits (1 Gyr) | S32+S98 | 6 min | **0.0054 ± 0.0032** | 0.0012–0.0108 | 0.0e+00 | 1.00 | 20 × 1e+06 |
| A | loc with TIC 8.2 stellar parameters | S32+S98 | 6 min | **0.0354 ± 0.0044** | 0.0274–0.0479 | 0.0e+00 | 1.00 | 20 × 1e+06 |
| A | audit baseline: Mann M, R; Teff 4030 K; ZP-corrected ϖ | S32only | 6 min | **0.0388 ± 0.0019** | 0.0360–0.0418 | 4.4e-07 | 1.00 | 20 × 1e+06 |
| A | audit baseline: Mann M, R; Teff 4030 K; ZP-corrected ϖ | S6only | 30 min | **0.0261 ± 0.0023** | 0.0230–0.0311 | 3.4e-07 | 0.99 | 20 × 1e+06 |
| A | audit baseline: Mann M, R; Teff 4030 K; ZP-corrected ϖ | S98only | 6 min | **0.0290 ± 0.0024** | 0.0250–0.0334 | 5.8e-09 | 1.00 | 20 × 1e+06 |
| B | earlier configuration (Teff 4080 K, uncorrected ϖ), audit fold | S32+S98 | 3.33 min | **0.1827 ± 0.0111** | 0.1640–0.2154 | 0.0e+00 | 1.00 | 20 × 1e+06 |
| B | audit baseline: Mann M, R; Teff 4030 K; ZP-corrected ϖ | S32+S98 | 6 min | **0.1791 ± 0.0111** | 0.1607–0.2117 | 0.0e+00 | 1.00 | 20 × 1e+06 |
| B | baseline − neighbours excluded by the pixel localisation | S32+S98 | 6 min | **0.1791 ± 0.0111** | 0.1607–0.2117 | 0.0e+00 | 1.00 | 20 × 1e+06 |
| B | loc + MOLUSC companion limits (5 Gyr) | S32+S98 | 6 min | **0.0638 ± 0.0200** | 0.0309–0.1012 | 0.0e+00 | 1.00 | 20 × 1e+06 |
| B | loc + MOLUSC companion limits (1 Gyr) | S32+S98 | 6 min | **0.0594 ± 0.0222** | 0.0173–0.1031 | 0.0e+00 | 1.00 | 20 × 1e+06 |
| B | loc with TIC 8.2 stellar parameters | S32+S98 | 6 min | **0.1282 ± 0.0155** | 0.0914–0.1556 | 0.0e+00 | 1.00 | 20 × 1e+06 |
| B | earlier configuration AND earlier fold (exact earlier run) | S32+S98 | 3.33 min | **0.1459 ± 0.0113** | 0.1275–0.1613 | 0.0e+00 | 1.00 | 20 × 1e+06 |
| B | audit baseline with the earlier fold | S32+S98 | 6 min | **0.1470 ± 0.0115** | 0.1285–0.1626 | 0.0e+00 | 1.00 | 20 × 1e+06 |
| B | audit baseline: Mann M, R; Teff 4030 K; ZP-corrected ϖ | S98only | 6 min | **0.1697 ± 0.0127** | 0.1507–0.2018 | 0.0e+00 | 1.00 | 20 × 1e+06 |

How to read it:
* **A**: 0.030–0.035 under the default assumptions (earlier configuration, audit baseline, localisation-cleaned, TIC
  parameters). The earlier 0.031 ± 0.003 is reproduced. Practically all of it is STP (a planet on an unresolved bound
  companion). All EB scenarios together ≤ 2 × 10⁻⁶; nearby-star scenarios ≤ 2.4 × 10⁻¹² (0 once the neighbour excluded
  by the localisation is removed). The "loc" rows equal the "base" rows because the removed scenarios carry ~10⁻¹² and
  the seeds are the same. With the MOLUSC companion limits: 0.005 ± 0.003 at either age.
* **B**: 0.18 ± 0.01 with the audit's fold (0.183 with the earlier stellar configuration); 0.13 ± 0.02 with the TIC
  stellar parameters (ρ★ 2.31 instead of 2.51 ρ☉), so B's FPP depends on the assumed host star at the ±0.03 level; 0.06 ±
  0.02 with MOLUSC. With the earlier fold and the earlier configuration the earlier value is reproduced exactly
  (0.1459 ± 0.011; largest per-seed difference 7 × 10⁻⁵). With the earlier fold and the audit configuration it is
  0.147 ± 0.012. So the 0.146 → 0.18 change comes entirely from the fold (epoch 3.3 min earlier, deeper in-transit
  mean), not from the stellar parameters or exposure time. EB scenarios 2–3 × 10⁻⁴.
* **Classification** (Giacalone et al. 2021 thresholds, with the caveats below): A is a "likely planet" in every variant
  and falls below the validation threshold (0.015) only in the MOLUSC variants. B is a "likely planet" in every variant
  and never below 0.015. Neither signal is validated.

Rules for the manuscript: (i) FPP is the probability that the signal is not a planet on the target *under TRICERATOPS's
scenario set and priors*. It is not the probability that the object is a planet. (ii) Quote run-to-run standard
deviations (Monte Carlo noise, 20 runs × 10⁶ draws); runs are reproducible given the seed. (iii) The multiplicity boost
(Guerrero et al. 2021: ≈ 20 for all SPOC multis, ≈ 54 below 6 R⊕) is not applied: those authors state it does not
apply to QLP-only candidates. (iv) The Giacalone et al. (2021) thresholds (validated: FPP < 0.015 and NFPP < 10⁻³;
likely planet: FPP < 0.5 and NFPP < 10⁻³) were calibrated with TRICERATOPS v1 and TRILEGAL backgrounds. 1.1.0 uses Gaia
DR3 backgrounds by default, which its own KNOWN_ISSUES list as not yet benchmarked. (v) The MOLUSC mode drops companions
with periastron < 10 au and sets the companion prior to the kept fraction (13 %). Its lower FPP partly reflects that
choice, not only data. (vi) NFPP ≈ 0 does not cover the 83.8″ star, because TRICERATOPS's Gaussian PSF gives it no
aperture flux; that star is excluded by the pixel analysis instead.

---------------------------------------------------------------------------------------------------------------------
## 8. Stellar activity
No rotation period is detected. TESS SAP periodograms are dominated by sector-length trends and, in S98, by the 7.3-d
orbit/scattered-light cadence. ASAS-SN g (1055 points, 2017–2024) peaks at 1.0 d (daily alias) and at 15.9 d and 8.2 d
(not significant). Neither transit period is a harmonic of these. A U-shaped, strictly periodic 1.9-h (A) or 1.7-h (B)
dip repeating over 7 years cannot be produced by spots or flares. The eRASS1 X-ray level (log L_X/L_bol ≈ −4.8 to −4.9)
is moderate and favours slower rotation than the ~8 d catalogue values; it is a single marginal detection.

---------------------------------------------------------------------------------------------------------------------
## 9. Two independent routes for the critical numbers
| quantity | route 1 | route 2 | agreement |
|---|---|---|---|
| P_A | joint transit fit 3.1827853 ± 0.0000027 | per-event times 3.1827793 ± 0.0000031 (pipeline) / 3.1827827 ± 0.0000032 (own); blind combined BLS 3.18279 | 1.5σ / 0.6σ |
| P_B | joint fit 11.132745 ± 0.000013 | per-event 11.1327478 ± 0.0000133 (own); blind BLS 11.13273 | 0.2σ |
| ρ★ | Mann M/R³ 2.51 +0.25/−0.22 | transit A (free) 2.75 +0.46/−0.97 (uniform prior) | consistent |
| R★ | Mann M_Ks 0.623 ± 0.019 | SED (ARIADNE) 0.621 ± 0.015; Gaia DR2 FLAME 0.618 | < 0.2σ |
| Rp_A | k R★ = 1.44 ± 0.06 | √(depth) R★ with depth = 448 ppm → 1.44; own photometry 1.41 ± 0.07 | yes |
| location A | PRF likelihood fit (S32+S98 clean) (+2.5, −7.5)″ ± 5.7″ | SPOC 2-min MOM_CENTR z = −1.7/+0.1 (empirical null); empirical-PRF shift (S32) z = +1.9 | consistent with target |
| contamination | catalogue census (full inclusion): nearest capable 43.5″ | PRF-integrated real apertures: 43.5–83.8″ capable in some apertures; all excluded by location | yes |
| FPP_A | TRICERATOPS 0.0304 (previous inputs) | 0.0312 (audit inputs, 6-min exptime) | yes |
| FPP_B | TRICERATOPS, earlier fold + earlier configuration: 0.1459 (earlier run: 0.1459) | audit fold: 0.179–0.183 | reproduced; fold-dependent |
| B beyond S98 | epoch-shifted scan statistic 10.38 (reproduced) | period+epoch-randomised null FAP 1.2 × 10⁻³; S6 rank test 1.6 × 10⁻³ | the earlier null method is invalid (§5.2) |

---------------------------------------------------------------------------------------------------------------------
## 10. What each statistic means
* **MES 10.6 / 11.0** (discovery): LEO-vetter's signal-to-pink-noise, dep/√(σ_w²/n_in + σ_r²/N_tr) (Pont+2006,
  Hartman & Bakos 2016), S98 only. Not a SPOC TPS MES.
* **Empirical S/N** (§4.4, §5.4): box depth / robust σ of 3000 same-width boxes at random phases. It includes red noise.
  The boxes sample only ~P/T14 independent phases (~37 for A, ~145 for B), so they give σ, not a p-value below ~1/37.
* **Per-event errors**: least-squares errors × √χ²_red × β (time-averaging red-noise factor).
* **FAP of the B link**: fraction of null scans with random period (9–13.5 d) and random epoch that reach the observed
  maximum, (n + 1)/(N + 1). The trials are not fully independent; the number of distinct windows is reported with it.
* **Rank-test p** (B in S6): the probability that the two strongest of N exchangeable, independent 1.6-h boxes form a
  pair at the predicted times: K × 3 / C(N, 2), with K = 27 positions inside the scan window and ±1 box in spacing.
* **Localisation z**: offset of the real fit from the mean of on-target injections, in units of their standard deviation
  (frequentist, calibrated on the star's own pixels).
* **ρ★, Rp, b intervals**: Bayesian credible intervals (16th–84th percentile) under the stated priors.
* **FPP / NFPP**: TRICERATOPS posterior scenario probabilities (Bayesian, model- and prior-dependent), with Monte Carlo
  standard deviations.
* **σ-equivalents** of FAPs are one-sided.

---------------------------------------------------------------------------------------------------------------------
## 11. Numbers that changed in this audit
| item | before | now |
|---|---|---|
| M★ error | ± 0.013 M☉ | ± 0.015 M☉ |
| Teff (adopted) | 4080 ± 100 K | 4030 ± 90 K |
| LD prior centres | 0.351, 0.289 | 0.33, 0.30 |
| P_A | 3.1827863 ± 0.0000028 d | 3.1827853 ± 0.0000027 d |
| Teq_A, S_A | 820 ± 13 K, 75 ± 5 S⊕ | 810 ± 22 K, 72 +8/−7 S⊕ |
| A odd/even | 1.5σ | 1.2σ (pipeline), 0.0σ (independent photometry) |
| A odd/even (box) | – | +2.2σ and −2.1σ in two photometries (opposite signs); per-event fits ≤ 1.2σ |
| A secondary error | ± 24 ppm (white) | ± 32 ppm (empirical) |
| A timing χ² | 120/21 (RNAAS) or 114/21 (report) | 47.5/19 (red-noise-inflated errors) |
| A median single-event S/N | 3.8 | 3.5 |
| A localisation | S6 1.6″ (0.2σ), S32 10.5″ (3.8σ), S98 19″ (1.6σ) | S32 0.75σ, S98 clean 1.45σ; combined (+2.5, −7.5)″ ± 5.7″; S6 unusable (momentum dumps) |
| S32 centroid / pointing shift | 4σ pointing shift, residual 3 ± 3″ | not significant (empirical null) |
| stars able to mimic A | five at 43–56″ | same, plus a T = 13.06 star at 83.8″ (via real PRF wings in SPOC apertures); all excluded |
| STP for A | not excluded | excluded on circular orbits (density + overluminosity) |
| B in S98 | 2 complete + 1 partial transit | 3 clean + 1 in a flagged interval (epoch 2) |
| B link null test | "exceeds all 300 trials" (null max 7.5) | earlier null invalid (overlapping trials). Period + epoch randomised: S6+S32 FAP 1.2 × 10⁻³; S6 alone 0 of 18,652; rank test in S6 p = 1.6 × 10⁻³; S32 alone 0.03 |
| B in S98 alone | "detected in S98" (MES 11.0) | not established by S98 alone: calibrated BLS FAP 0.33 (discovery mask), 0.008 (own photometry), ~1 (default mask) |
| B overall | implied FAP < 1/300 from the link alone | indicative 5 × 10⁻⁵ to 1 × 10⁻² (3.9σ to 2.3σ; S98 × S6 × S32, Fisher) |
| B depth difference S98 vs archival | 3.1σ (RNAAS) / 2.4σ (report) | 2.2–2.7σ |
| Rp_B | ≈ 2.0 (1.9–2.2) R⊕ | 2.25 ± 0.12 R⊕ if on the target; 2.0–2.1 without the density assumption |
| B localisation | 0.5/5.8/5.6″ within 1.3σ | combined (+1.0, +1.4)″ ± (2.3, 3.0)″ |
| FPP_B (default) | 0.146 ± 0.011 | 0.18 ± 0.01 with the audit's fold; 0.13 with TIC stellar parameters. The earlier fold reproduces 0.1459 exactly (per-seed difference < 10⁻⁴), so the change comes from the fold alone |
| TGLC S6 depth of A | 190 ± 171 ppm | 982 ± 101 ppm (aperture flux); S6 pipeline spread 521–982 ppm |
| S110 timing of B | ± 23 min | ± 2.7 min if the S6/S32 link holds (± 49 min on S98 alone) |
| rotation | "unknown" / catalogue ~8 d | not detected; ~8 d not confirmed |

---------------------------------------------------------------------------------------------------------------------
## 12. TESS Sector 110 (2026-10-31 → 11-26; official HEASARC/MIT pointing; camera 1, CCD 2, col 702.8, row 1577.0)
A (±2.0 min, 1σ): 2026-11-01 03:49, 11-04 08:13, 11-07 12:36, 11-10 16:59, 11-13 21:22, 11-17 01:45, 11-20 06:09,
11-23 10:32, 11-26 14:55 UTC (the last may fall after the end of the sector; downlink gaps unknown).
B (±2.7 min if the archival link is correct): 2026-11-01 11:31, 11-12 14:42, 11-23 17:53 UTC.
Uncertainties above are from the joint-fit chains. Per-event linear ephemerides with errors scaled by √(χ²/dof) (the
timing excess, §14 R3) give A ±2.1–2.2 min and B ±2.4 min (own photometry, 8 events) to ±5.1 min (pipeline, 4 events).
Use the larger values when planning windows.

---------------------------------------------------------------------------------------------------------------------
## 13. Priority wording (the only defensible form)
"As of 2026 September 25 (UTC), we found no prior public report of either signal in the searched sources (ExoFOP TOI and
CTOI lists, SPOC TCE records, the NASA Exoplanet Archive, T16/RAVEN/DTARPS candidate lists, SIMBAD, VizieR, ADS, arXiv)."
Do not write that the signals are new or unknown in an absolute sense.

---------------------------------------------------------------------------------------------------------------------
## 14. Adversarial review (a referee seeing this for the first time)

Each item: the objection → what the data say → what remains → the wording the manuscript needs.

**R1. "A is an instrumental artefact."** → Recovered alone in each of three sectors (2018, 2020, 2025) on two CCDs, with
23 events. No event dominates: the largest single share of (S/N)² is 24–39 % per sector, and every leave-one-out S/N is
≥ 7.2 per sector. Depths agree across sectors in the SPOC/TESS-SPOC/QLP products (521 ± 53, 476 ± 52, 511 ± 47 ppm). In
S6 every A transit lies 0.4–3.5 h from a momentum dump. The offsets drift from event to event, so a dump systematic
cannot fold coherently at P_A, and S32/S98 have different dump cadences. → What remains: why SPOC issued no TCE in S32
(S32-only S/N 6.7–9.1, near the 7.1 threshold; U13). → Wording: "four independent reductions of the same pixels", not
"four independent detections".

**R2. "B is an artefact of scattered light or noise."** → Largely justified for S98 alone. A calibrated BLS test gives
FAP 0.33 in the discovery (frozen-mask) data, ~1 with the QLP default mask and 0.008 in the own photometry. The B-free
S98 data contain a peak as strong as B (§5.3). The independent support is S6: its two strongest dips are the two
predicted events (p = 1.6 × 10⁻³, no noise model), and S32 adds FAP 0.03. The indicative overall FAP is 5 × 10⁻⁵ to
1 × 10⁻² (3.9σ to 2.3σ). → What remains: B's reality rests on about 3σ; S110 (three predicted transits) will settle it.
→ Wording: "a second, low-significance candidate signal"; never "confirmed", never "FAP < 5 × 10⁻⁵"; state the S98-only
result.

**R3. "The timing is inconsistent."** → All products are BJD_TDB. TESScut's CCD-centre barycentring is corrected
(−11.8, −3.3, +0.5 s). QLP S98 shows no barycentric drift (T5, T6). A's O−C scatter exceeds the errors (χ² 47.5/19;
42.6/22 in the independent photometry). → What remains: single-event timing errors at S/N ≈ 3.5 are probably ~1.5×
underestimated. The joint-fit P and T0 errors carry no red-noise term. The scaled per-event errors (P ± 3.1 × 10⁻⁶ d) are
~15 % larger than the joint fit's. → Wording: "no evidence for TTVs; the timing scatter is larger than the formal errors";
the S110 windows should use the scaled per-event ephemerides: A ±2.1–2.2 min instead of ±2.0 min; B ±2.4 min (own
photometry, 8 events) to ±5.1 min (pipeline, 4 events) instead of ±2.7 min.

**R4. "The dip is on another star."** → A: PRF localisation (S32 2-min + S98 clean, calibrated by injections into the
real pixels) puts the source at (+2.5, −7.5)″ ± 5.7″ from the target. Every Gaia DR3 star able to produce A (43–84″,
including the T = 13.06 star that only the real PRF puts into the SPOC apertures) is ≥ 8σ away. The flux-weighted
centroids show no shift. S6 cannot be used (momentum dumps). B: (+1.0, +1.4)″ ± (2.3, 3.0)″. → What remains:
everything within ~6″ (A) or ~3–4″ (B) that Gaia does not resolve, and any bound companion within ~1″. No high-resolution
imaging exists (U1). No known-planet control of the PRF method (U12). → Wording: "the transit source is consistent with
the target within ~6″ (1σ); all known neighbours able to produce the signals are excluded".

**R5. "It is an eclipsing binary."** → On the target: excluded for both signals by the Gaia RVs (M sin i < 0.016 M☉ at
3.18 d) together with the transit depth. A: flat bottom (V disfavoured at Δχ² ≈ 32 with white errors, ≳ 4σ with red
noise), no secondary (< 95 ppm), no P/3 or 2P signature (the box odd/even signs flip between photometries). All
TRICERATOPS EB scenarios together: ≤ 2 × 10⁻⁶ for A and 2–3 × 10⁻⁴ (at most 9 × 10⁻⁴ in a single run) for B. B:
odd/even within 1.5σ, no phase-0.5 dip (61 ± 62 ppm), V
disfavoured at Δχ² 14.5 (pipeline) but only 7.4 (own photometry). → What remains: for B, an EB on a bound companion is
restricted only by TRICERATOPS's priors. → Wording: "an eclipsing binary on the target is excluded".

**R6. "The uncertainties are underestimated."** → Found and fixed in this audit: bootstrap localisation errors (×3
too small), white-noise errors on secondary and odd/even, the B-link null, and the M★ error. Still open: joint-fit
ephemeris errors (see R3); FPP errors are Monte Carlo only, while the model and prior systematics are larger (B moves
0.146 → 0.18 with a 3-min fold shift; §7); Rp_A depends on the limb-darkening treatment at the 1σ level (1.44 vs 1.50
R⊕ with free LD); R★ carries the Mann-relation edge effect (U10). The free-density emcee chains stopped at N/τ = 36–50
(below the 60τ target), but dynesty with free LD gives the same density bound (ρ < 3.51 ρ☉ at 97.7 %). → Wording:
quote Rp_A = 1.44 ± 0.06 R⊕ and state that free limb darkening gives 1.50 +0.09/−0.08 R⊕.

**R7. "Results are prior-driven."** → A: b and a/R★ are set by the stellar-density prior. Rp_A is not (1.44–1.50 across
all priors). The density upper bound is data-driven and reproduced by two samplers. B: b, Rp and any density statement
are prior-driven (U2). FPP is prior-driven by construction. The MOLUSC-mode FPP partly reflects a code choice (the
companion prior is set to the kept fraction; §7). → Wording: give Rp_B as 2.0–2.4 R⊕ with the dependence stated, and
FPP values "under TRICERATOPS's default assumptions".

**R8. "Significance is inflated."** → Four inflations were found and corrected: the B-link FAP (A14/C7; < 5 × 10⁻⁵ →
~2 × 10⁻³); B's S98 "detection" (MES 11.0 is not a calibrated significance; the calibrated S98-only FAP is 0.008–0.33);
the reading of 3000 random boxes as 3000 independent trials (C9); and the S32 3.8σ localisation offset (A2). The discovery "MES 10.6/11.0" is LEO-vetter's pink-noise S/N on S98 only. It is not a SPOC MES and must not be
compared with the 7.1σ SPOC threshold. The "four photometries" share the same pixels. → Wording: state every
significance with its error model (§10).

**R9. "Citations are misused."** → Checked in the primary texts: Mann+2015 (with the 2016 erratum; the M_K-only
relation is unaffected), Mann+2019 Tables 6–7, Lindegren+2021 zero point, Claret 2017 table15, Stassun+2019 T-mag,
Giacalone+2021 thresholds, Guerrero+2021 boost, Way+2026 overluminosity, Gaia DR3 RV variability criteria, TESScut
barycentring, the HEASARC S98 scattered-light notes. Fragile points: the Giacalone+2021 thresholds were calibrated for
TRICERATOPS with TRILEGAL backgrounds, whereas 1.1.0 uses Gaia DR3 backgrounds that its own KNOWN_ISSUES call
unbenchmarked. The Guerrero+2021 boost does not apply to QLP-found candidates. QLP bit 2^29 has no official meaning;
cite the S98 release notes for the scattered light, not the QLP documentation. Way+2026 is most sensitive to
near-equal-mass companions, which is exactly the case used (q ≥ 0.9 for the circular-STP argument).

**R10. "The conclusions are stronger than the evidence."** → Not allowed: "validated" (FPP_A = 0.031 > 0.015; the
MOLUSC-mode 0.005 rests on a code-level choice); "B confirmed"; "new" or "first" without the searched-sources wording
(§13); "in 7:2 resonance" (only the period ratio is measured); a rotation period; "habitable" (Teq 810 and 530 K);
Rp_B = 2.25 R⊕ without its prior dependence; "EB excluded" for B's bound-companion channel. Allowed: A is a strong
planet candidate ("likely planet" by the Giacalone et al. 2021 criteria under stated assumptions); B is a
low-significance candidate signal (about 3σ overall) awaiting Sector 110.

---------------------------------------------------------------------------------------------------------------------
## 15. Readiness

**Verdict: technically ready for drafting a planet-candidate manuscript. NOT "ready for publication" by the criteria
set for this audit.**

| criterion | status |
|---|---|
| all critical numbers reproduced | Yes for the star, A (ephemeris, shape, radius, localisation, FPP) and B's ephemeris, shape and FPP (the earlier FPP_B is reproduced exactly with the earlier fold). Not reproduced: the earlier significance of B's archival link (withdrawn, §5.2) and B's detection significance in S98 alone (not established, §5.3). |
| contradictions resolved | **No.** B's S98 vs archival depth (2.2–2.7σ, U4); the S6 depth scale between photometries (U17); A's timing excess (U8). |
| methods documented | Yes (METHOD_AUDIT.md; every number traceable to `scripts/`, `out/`, `logs/`). |
| provenance complete | Yes (DATA_PROVENANCE.csv; outputs_sha256.json). |
| no known bug affects results | Yes. Every bug found (A1–A14, C1–C10) is fixed or shown to have no effect; 19/19 regression tests pass. |

What the package supports now:
1. **Signal A as a strong planet candidate:** P = 3.1827853 d, Rp = 1.44 ± 0.06 R⊕, on the target within ~6″, FPP ≈
   0.03 under TRICERATOPS defaults ("likely planet", not validated).
2. **Signal B only as a low-significance candidate signal** (about 3σ overall), with its Sector 110 prediction.

What would move it to "ready for publication" as a validation or two-planet paper:
1. High-resolution imaging (closes the bound-companion channel, U1).
2. TESS Sector 110 (2026-10-31 → 11-26): three predicted B transits (±2.7–5.1 min) confirm or reject B; eight to
   nine A transits refine A.
3. An explicit resolution or treatment of U4 and U17 in the text.

A candidate-level note (A as a candidate, B as a tentative signal) can be drafted now, provided it uses the allowed
wording in CLAIM_AUDIT.csv and states U1–U4, U8, U17 and U18.
