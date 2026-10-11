# Method audit — TIC 4206066 (StKM 1-561)

For every method behind a manuscript number: what was done earlier, which primary source defines the method, how the
audit re-did it, and the result. Scripts are in `scripts/`, outputs in `out/`, logs in `logs/`, figures in `plots/`.
Environment: separate `uv` venv (Python 3.12, `logs/audit_env_freeze.txt`); TRICERATOPS in its own Python 3.11 venv
(TRICERATOPS 1.1.0, numpy 1.23.5, setuptools < 81, the pins that TRICERATOPS 1.1.0 itself declares). The project `.venv`
and every frozen pipeline file were only read.

## 1. Data provenance (`s01`, `s02`, `s03`, `s21` → `DATA_PROVENANCE.csv`)

* All 16 time-series products on MAST for this star were re-downloaded on 2026-09-25 05:30 UTC and compared with the
  files used by the deep follow-up: **16/16 byte-identical** (sha256). No product exists beyond S6, S32, S98. There are no
  DV/TCE products and no SPOC or TESS-SPOC S98 product. TESScut cutouts (S6, S32, S98) were re-made: flux and time arrays
  are identical to the earlier cutouts.
* Versions: SPOC S32 2-min `spoc-5.0.20`, DATA_REL 48, TIC 8.1 in header; TESS-SPOC S6/S32 v1; QLP S6 (TICA), S32 (TICA
  v0.2.1), S98 (TICA v1.4.0, TGLC-era columns, DETMET = QSP); TESScut from SPOC FFIs (astrocut 0.11.1).
* Time system: every product declares TIMESYS = TDB, TIMEREF = SOLARSYSTEM, BJDREFI = 2457000 (`s06`).
* Catalogues re-queried fresh: Gaia DR3/DR2 (ESA TAP), TIC 8.2 (MAST), 2MASS, AllWISE, CatWISE, Tycho-2, APASS DR9,
  LAMOST DR5/DR7, Lépine & Gaidos 2011, GALEX, ROSAT 2RXS, eRASS1, VSX, WDS, Bailer-Jones 2021. The eRASS1 X-ray source
  (see §12) was not in the earlier reports.
* Re-used without re-download: the ASAS-SN Sky Patrol v2 file (pyasassn, 2026-09-25T03:20:56Z), the two Way et al. 2026
  catalogue rows, and the ARIADNE output files. Their sha256 are in DATA_PROVENANCE.csv.

## 2. Stellar parameters (`s04`)

| quantity | method and primary source (checked in `data/refs/`) | audit value | earlier value | status |
|---|---|---|---|---|
| M_Ks | 2MASS Ks 7.694 ± 0.029; Gaia DR3 ϖ = 28.0750 ± 0.0176 mas; Lindegren+2021 zero point −0.0332 mas (gaiadr3-zeropoint 0.1.0) | 4.938 ± 0.029 | 4.936 | VERIFIED (ZP shifts M_Ks by +0.003) |
| R★ | Mann+2015 Table 1: R = 1.9515 − 0.3520 M_K + 0.01680 M_K², 2.89 %. Checked in the arXiv v3 source and the **erratum ApJ 819, 87 (2016)**, pages 39–40 of the CHARA copy. The erratum fixes a press error in the *published* Table 1 (the [Fe/H] form printed f = −0.04458 instead of +0.04458). The M_K-only form is identical in all three versions. Validity 4.6 < M_K < 9.8, K7–M7, 2700–4100 K; this star sits at the hot/bright edge. | 0.623 ± 0.019 R☉ | 0.623 ± 0.019 | VERIFIED |
| R★ ([Fe/H] form) | same table, f = +0.04458; [Fe/H] = −0.27 (GSP-Spec, flagged) or 0 | 0.613 / 0.621 | – | alternative (±1.5 %) |
| M★ | Mann+2019 Table 6, n = 5, zp = 7.5 (arXiv 1811.06938v2 source). Relation error from their Table 7 (2.4 % at M_K = 5). Second route: official posterior code/data (github awmann/M_-M_K-, commit ffe5f19, `Mk-M_7_trim.fits`) gives 0.611 ± 0.015. | 0.606 ± 0.015 M☉ | 0.606 ± 0.013 | value VERIFIED, error CORRECTED |
| ρ★ | M/R³ by Monte Carlo (correlated M_K error included) | 2.51 +0.25/−0.22 ρ☉ (3.53 g cm⁻³) | 2.51 ± 0.23 | VERIFIED |
| log g | G M/R² | 4.632 ± 0.028 | 4.63 | VERIFIED |
| Teff | TIC 8.2 4010 ± 120; Mann+2015 erratum Table 2 V−J relation 4048 ± 86 (3984 with [Fe/H] = −0.27); ARIADNE 4012 +45/−40; Gaia DR2 Apsis 4069. GSP-Spec 4182 rejected (extrapolation flag, logg 5.01). | **4030 ± 90 K** | 4080 ± 100 | CORRECTED |
| L★ | R² (Teff/5772)⁴ | 0.092 ± 0.010 L☉ | – | new. Cross-checks: ARIADNE 0.090 ± 0.006, TIC 0.098, Mann BC_K 0.107, BC_V 0.101, FLAME 0.094 |
| Tmag | Stassun+2019 eq. (TIC 8 paper source, arXiv 1905.10694) | 9.685 (TIC 9.682) | – | VERIFIED |
| extinction | d = 35.6 pc | A_V ≈ 0 | – | negligible (Gaia MSC A_0 0.5 is the binary-model output, not used) |

Consistency: the SED radius (ARIADNE, re-read from `best_fit_average.dat`, not re-run) is 0.621 ± 0.015 R☉. Gaia DR2
FLAME gives 0.618 R☉. The TIC gives 0.648 ± 0.065 R☉ (ρ 2.31). Pecaut & Mamajek places M_K = 4.94 between K7V and K8V.

## 3. Light curves (`s05`, `lcio.py`, `s17`)

* Pipeline set (primary): TESS-SPOC S6 PDCSAP (30 min), SPOC S32 PDCSAP (2 min), QLP S98 SYS_RM_FLUX (200 s).
* Independent set: own aperture photometry from TESScut (S6 r = 2.5 px, S32 r = 2.0 px, S98 r = 2.0 px, same method
  in every sector; background = median of faint pixels > 5 px from the target) and from the S32 2-min TPF.
* Quality masks: lightkurve DEFAULT_BITMASK (17087) for SPOC/TESScut. Straylight (2048) is informational in the TESS
  SDPDD. For QLP S98 the official QLP pages (tess.mit.edu/qlp, 2023-11-21 and 2025-03-31 notes) document only bits 3, 6,
  12 (SPOC-derived) and 2^30 (QLP internal flag). They say "Bits 16–30 are unused", yet 2^29 is set in 44 % of S98
  cadences. It coincides with the Earth-scattered-light intervals at the start of each orbit that the official S98
  notes describe for camera 1. Baseline: the frozen mask (removes 2^29 and 2^30), as in discovery. Variants: keep
  2^29; drop straylight.
* Detrending: wotan biweight 0.75 d with **both** signals masked (±1.5 T14). Variants: 0.5, 1.0, 1.5 d biweight;
  robust cosine 1.0 d; none (local linear baseline only). Per-event depths change by < 5 % for A and < 5 % for B
  (`out/s17_detrend_robust.json`). **No result depends critically on the detrending setting.**
* Dilution: SPOC CROWDSAP 0.993–0.997 (applied in PDCSAP). QLP S98 is TGLC forward-modelled. TIC contamination ratio
  0.0027. Dilution effects are < 1 %, far below the depth errors.

## 4. Period / ephemeris (`s07`, `s08`, `s09`, `s19`)

* BLS (astropy BoxLeastSquares, not the frozen pipeline's implementation), log-uniform grid d ln P = dur_min/(5·baseline).
  All three sectors combined, blind (30-min bins): once two ~14-d peaks from S98 scattered light are removed, the
  strongest peaks are **B at 11.13273 d (SDE 19.5) and A at 3.18279 d (SDE 13.7)**. Per sector: A is recovered alone in
  S6 (3.1826 d, SDE 8.7), S32 (3.1825 d, SDE 7.6) and S98 (3.1816 d). B alone in S98 (11.1328 d, SDE 9.8, QLP default
  mask; 11.1326 d, SDE 7.4, own photometry). With the discovery (frozen) mask B is only the 2nd peak (SDE 5.5, behind a
  5.60-d systematics peak; `s07`). Its S98-only false-alarm probability (`s07b`, `s07c`; null from segment-scrambled
  S98 data) is 0.33 (discovery mask), 0.008 (own photometry) and ~1 (default mask). In S6 B appears at B/2
  (5.557 d, 789 ppm).
* Joint fit (batman + emcee, per-set offset and jitter, both signals masked in detrending, the other signal excluded
  from the windows): see PUBLICATION_TRUTH.md. Independent routes: per-event linear ephemeris (pipeline and own sets)
  and dynesty (free LD). They agree within 1.5σ (the per-event linear P for A is 1.5σ below the joint fit).
* Cycle count between sectors: the S98-only period (LEO trapezoid) has σ = 0.00031 d (A) and 0.0011 d (B). The
  alias spacing (P/N) is 13–18σ (A) and 45–62σ (B), so the S98→S32/S6 cycle count is unique. It is also chosen by the
  blind combined BLS.
* Time systems (`s06`): TESScut times are barycentred to the CCD centre (astrocut docs). Target-specific shift:
  −11.8 s (S6), −3.3 s (S32), +0.5 s (S98). After the shift, TESS-SPOC and QLP times agree with TESScut within 1–5 s.
  QLP S98 shows no drift across the sector (1.67/1.64/1.63/1.60 s), so it is **not affected** by the documented QLP
  barycentric errors (S74–79, S99–104). The SPOC 2-min and TESS-SPOC products carry target-specific BJD_TDB.

## 5. Events (`s09`, `s09b`, `s11`)

Per-event fits: fixed-shape template (joint-fit shape) + local linear baseline. Depth errors are scaled by
√(χ²_red) and by the time-averaging red-noise factor β (Pont+2006/Winn+2008). Timing by a Δχ² = 1 scan. Signal B's
link test is reproduced with the earlier statistic (maximum S/N over a ±0.004 d period scan at the S98 epoch).
Detrending masks only A, so the real phase gets no advantage over the null phases. The fast scan is validated against
brute force at every call.

**Null distribution (corrected).** The earlier null (and the audit's first version, `s11`) shifted only the epoch at
the fixed S98 period. The scan window spans ±0.9 d in S6 and ±0.66 d in S32, so such trials overlap: in S6 there are
~5–16 distinct windows among 20,000 trials. The valid null (`s11b`) draws the period uniformly in 9–13.5 d (|P − P_B| >
0.1 d) and the epoch uniformly, scales the scan grid so the time window is the same, and removes B's own archival
transits from the null data. Distinct argmax windows are counted to show how independent the trials are. A
distribution-free rank test (`s11c`) checks whether the two predicted S6 dips are the two strongest independent 1.6-h
boxes of the sector. The Gaussian-model FAP uses the measured per-position S/N scatter. The S98-only significance of B is
tested separately with BLS on segment-scrambled S98 data (`s07b`).

## 6. EB tests (`s10`, `s10b`)

Box depths with **empirical** errors from 3000 random-phase boxes of the same width (red noise included); phase scan
with look-elsewhere count; P/3 phases; n mod 3 groups (3P test); odd/even both as box depths and as per-event template
fits (matched events across photometries); trapezoid shape fit with ingress fraction profile (V vs U). The random boxes
sample only ~P/T14 independent phases (~37 for A), so they set σ, not p-values below ~1/37. They are drawn at phases
0.1–0.9, while the scan covers 0.05–0.95; re-deriving σ over 0.05–0.95 changes it by < 3 % (`s10c`). `s10c` also traces
each out-of-transit scan feature to the individual cycles that produce it.

## 7. Transit model (`s08`, `s08b`)

* batman 2.5.x (package `batman-package`), quadratic LD, circular orbit. Supersampling ceil(cadence/2 min), exposure
  time = cadence.
* LD: Claret (2017) J/A+A/600/A30 **table15** (quasi-spherical PHOENIX-COND, TESS, LSM), downloaded from CDS. The
  earlier hard-coded grid values match the table exactly. Bilinear interpolation at 4030 K, log g 4.63 gives
  u1 = 0.33, u2 = 0.30. Gaussian priors σ = 0.1.
* Priors: T0, P uniform (narrow); ln k uniform; b ~ U(0, 1+k); ρ★ either ln-uniform [0.1, 50] ("free") or Gaussian
  2.507 ± 0.23 ("Mann", with the Jacobian for sampling ln ρ). The free-ρ posterior is also importance-reweighted to a
  uniform-in-ρ prior to show prior sensitivity.
* Sampler: emcee, 48 walkers, convergence N > 60 τ_max and |Δτ|/τ < 2 %. Mann-prior fits converge in 28–36k steps
  (N/τ ≈ 90–110). Free-ρ fits hit the 80k-step cap at N/τ = 36–50 (A) and 29–43 (B). That is the b–ρ degeneracy, and
  for B it reflects an unconstrained density.
* Independent configuration: dynesty (600 live points, rwalk, dlogz 0.1), Kipping q1, q2 ~ U(0, 1), biweight 1.0 d.

## 8. Pixels / PRF (`s12`, `s12b`, `s12c`, `s12d`, tests)

Official PRF only: MAST `start_s0004` models (cam1_ccd1, cam1_ccd2; 50 files downloaded) read through TESS_PRF 0.1.3,
which un-interleaves `prf[8-i::9, 8-j::9]`. The interleaving and orientation are tested empirically (R1). The target's
position comes from a PRF fit of the median image with all Gaia DR3 sources (T < 19, proper motion propagated to each
sector). Difference images: per-transit linear per-pixel baseline, mean over transits. Noise map: per-pixel standard
deviation of 300–400 mean difference images at random fake epoch sets. Errors are **not** bootstrap: 150–300
model-free injections into the real pixels (background-subtracted pixels × (1 − δ) within the injection radius; radius
3, 6 and 30 px give identical results), fitted exactly like the data. The fake epochs of every noise map and injection set
are drawn one by one at random times that avoid both signals (±3 T14) and each other (> 6 T14), not at a common phase.
The trials therefore share windows but not a phase. The ~50–80 usable windows per sector limit the precision of the
injection σ to ~10 %, which does not change any conclusion (all excluded stars are ≥ 8σ away). Neighbour hypotheses: fixed-position fits and
PRF-shaped injections on each neighbour with the aperture flux loss matched. Three further routes: flux-weighted centroid
shift, empirical-PRF (shifted OOT image) fit, per-epoch fits. S32 SPOC MOM_CENTR/POS_CORR are tested with an empirical
null.

## 9. Sector-dependent depth (`s09`, `s17`)

The same per-event method in every sector and photometry, six detrendings. Dilution corrections are negligible (§3).
Cadence is handled by the supersampled template.

## 10. Hidden companions (`s14`)

Gaia DR3 census within 180″ with PRF-integrated aperture fluxes in the real apertures. RUWE is compared with 828 Gaia
DR3 stars of similar G, BP−RP and parallax (TAP query saved). Transit-density host constraint vs the Pecaut & Mamajek
dwarf sequence (v2022.04.16, downloaded). The Way+2026 overluminosity row is verified against the paper's definition
(arXiv 2602.20096 source: threshold O > 0.376 = half the equal-mass offset; best for near-equal-mass, ≲ 20 au). WISE
colours are compared with the dwarf sequence. The Gaia RV epoch scatter is derived from the documented
radial_velocity_error formula (Gaia DR3 datamodel: √(σ_med² + 0.11²), σ_med = √(π/2) σ/√N). Variability criteria
follow the Gaia DR3 documentation (p < 0.01 and gof > 4 and N ≥ 10).

## 11. FPP (`s13`)

TRICERATOPS 1.1.0 source read: FPP = 1 − (P_TP + P_PTP + P_DTP) and NFPP = Σ P(nearby scenarios) (triceratops.py
l. 1697–1701). Without contrast curves, companions are unresolved within 2.2″ and the Moe & Di Stefano 2017
bound-companion prior applies (marginal_likelihoods.py l. 483–493). With `molusc_file`, only companions with
a(1−e) > 10 au are kept and the companion prior is **set to 0 in log** (probability 1), with the array padded with
q = 0 entries to N. The implicit companion rate is therefore the kept fraction (13 % here). Version 1.1.0 changed two
defaults (CHANGES.rst): input assumed dilution-corrected, and the background population comes from Gaia DR3 instead
of TRILEGAL. KNOWN_ISSUES.rst says the Gaia-vs-TRILEGAL benchmark validation is **not done yet**. Validation thresholds
(Giacalone+2021 arXiv source): validated FPP < 0.015 and NFPP < 10⁻³; likely planet FPP < 0.5 and NFPP < 10⁻³. They
were calibrated with 2-min data; with 30-min data no validation region exists. The multiplicity boost (Guerrero+2021
text) is ≈ 20 for all SPOC multis and ≈ 54 below 6 R⊕. The authors state it does **not** apply to QLP-only candidates.
It is kept separate and not applied.

## 12. Activity (`s15`)

Lomb-Scargle of non-cotrended fluxes (SAP, QLP SAP, own raw) and ASAS-SN g/V, with spectral window and bootstrap levels
(shuffle bootstrap underestimates red-noise false alarms; stated). eRASS1 X-ray: DET_LIKE_0 = 8.04, 4.5 counts,
F(0.2–2.3 keV) = 3.4e−14 erg s⁻¹ cm⁻², 4.8″ from the target → L_X = 5.1e27 erg s⁻¹, log L_X/L_bol ≈ −4.8 to −4.9.

## 13. Novelty (`s18` + sealed package)

The sealed priority re-check (2026-09-25T04:46–05:25Z, 36 source groups, 552 queries, aliases, P/3–3P) was reviewed.
The fast-changing sources were re-queried at 06:30–06:34Z: ExoFOP TOI (8148) and CTOI (5137) lists, ExoFOP target
page (no TOI/CTOI/files/notes), NASA toi by TIC (0 rows), exo.MAST TCEs (`{"TCE":[]}`) — all negative. arXiv API and
SIMBAD TAP refused scripted access at 06:30Z. A web search for the identifiers found nothing.

## 14. Sector 110 (`s16`)

Official pointing from HEASARC TSSC `sector110_summary.html` and tess.mit.edu `sector-110`. Both give dates 2026-10-31
to 11-26, orbits 231–232, camera 1 boresight (75.883, −2.692), roll −38.493, spacecraft (106.067, −23.793, −121.301).
tess-point 0.9.5.post1 carries the same pointing (106.0674, −23.7934, 238.6986) and places the star on camera 1 /
CCD 2, column 702.8, row 1577.0.

## 15–17. Code, statistics and two-route checks

See BUG_REGRESSION_TESTS.md (code), CLAIM_AUDIT.csv (per-claim statistic definitions) and PUBLICATION_TRUTH.md §9
(two-route table).
