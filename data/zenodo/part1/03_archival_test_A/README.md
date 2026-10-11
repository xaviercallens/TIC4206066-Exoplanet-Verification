# 03_archival_test_A: archival false-alarm test for the 3.18 d signal

Printed in the Note (Sec. 2): *"An epoch search at the Sector 98-only period recovers the signal in Sectors 6 and 32 (false-alarm probability ≲1% per sector, ≲0.1% jointly, trials included)."* Fact `A.archival_FAP`.

## The test in plain steps

1. **Sector 98-only prediction** (`s1_s98_prediction.py` → `s1_s98_prediction.json`, `s1_boot_P_Tref.npy`, `s1_s98_prediction.log`). The QLP Sector 98 light curve alone (SYS_RM_FLUX, QUALITY == 0, biweight 0.75 d) is fitted with a symmetric trapezoid; then the period and a reference epoch are profiled with that shape fixed. The rule was written into the script before it was run, and nothing from Sectors 6 or 32 enters. Result: P = 3.18211 ± 0.00130 d (event bootstrap; ± 0.00057 d from the likelihood with red-noise factor 1.27). Propagated over ~800 cycles, the 1σ timing uncertainty in Sectors 6 and 32 is 25 h and 18 h, i.e. 4-13 transit durations, so the prediction fixes the period but not the phase.
2. **Statistic** (`s2_archival_fap.py`, helpers in `fap_common.py`). For a trial period and epoch, a box of width T14 = 1.87 h is placed at every predicted mid-time in the sector; S/N = mean box depth × √N_in / point-to-point scatter. The archival light curves are detrended once (biweight 0.75 d), with only the 11.13 d signal masked, so the 3.18 d phase is never favoured. Joint S6+S32 statistic: inverse-variance-weighted box depth over both sectors.
3. **Epoch search.** MAIN (fixed epoch at the predicted centre) finds no dip: S/N −0.5 (S6), +0.5 (S32); FAP 0.61 / 0.25 / 0.34. That test is uninformative because of step 1. OPT then maximises the S/N over epoch and period inside the 3σ prediction ellipse. Observed: S/N 9.84 (S6), 8.98 (S32), 13.22 (joint). The recovered period, 3.18276 d, is 0.02σ from the final multi-sector ephemeris.
4. **Null (identical machinery).** 20,000 trials with fake periods drawn from U(2.5, 4.5) d, random epochs, the true 3.18 d transits removed from the data, and exactly the same detrending, prediction uncertainty, lever arms, grids and windows. Null maxima: 4.59 (S6), 8.34 (S32), 7.37 (joint). The observed value is never exceeded (0 of 20,000 in each case).
5. **Counting independent trials (the N_eff caveat).** Neighbouring fake periods give correlated statistics, so 0/20,000 is not a 5 × 10⁻⁵ result. The FAP bound is (k + 1)/(N_eff + 1) with N_eff the number of independent null realisations.
   - The author's run counted N_eff heuristically: 164 (S6), 195 (S32), 3658 (joint) → bounds 0.61%, 0.51%, 0.027%; this gave the first printable sentence "<0.7% per sector, <0.03% jointly" (`a_fap_results.json` → `printable`), now superseded.
   - The pass-3 reviewer measured N_eff directly (`reviewer_neff/`): sample count divided by the integrated autocorrelation time of the OPT statistic on a log-spaced period grid (`x_neff2_main.json`: S6 77, S32 30, joint-with-fixed-anchors 43 out of 5000 periods), and by counting declustered exceedances above the 95/98/99% quantiles (`x_neff2_clusters_declustered.json`: pairs of [number of clusters, implied N_eff = clusters/(1 − q)] for minimum gaps of 1, 10, 30 and 60 grid samples, for the main and the 0.375 d edge-cut runs). With these, `PASS3_NOTE.json` records S6 1.0-1.3% and joint 0.006-0.06%. `FINAL_FACTS.csv` (`A.archival_FAP`) summarises the reviewer's N_eff as 77-140 (S6), 30-300 (S32) and 1.7 × 10³-1.6 × 10⁴ (joint, random relative anchors). The Note therefore prints the rounded, conservative "≲1% per sector, ≲0.1% jointly" (about 2.3-2.6σ per sector, 3.2-3.8σ jointly).
6. **Robustness** (`a_fap_edge0375.json`, `a_fap_w050.json`, `a_fap_w100.json`; `detrending_variants` in the results). Biweight windows of 0.5 / 0.75 / 1.0 d and a 0.375 d edge cut change the optimised S/N by −1.1 to +0.7 and never produce a null exceedance. An N2r null with periods in 3.10-3.26 d (bracketing the 3.125 d Sector 6 momentum-dump cadence) has lower maxima than the main null. The S32 null tail is set by segment-edge detrending artefacts. The cycle count S98→S6 (798) and S98→S32 (575) agrees with the final ephemeris.

All caveats are listed verbatim in `a_fap_results.json` → `caveats`; do not print "FAP X at the S98-predicted times" (`printable.do_not_print`).

## Exact statistical definition (added 2026-09-25 after the pass-3 review)

**What is computed.** A one-sided *permutation-type p-value* of the maximised box S/N, not a confidence bound.
For each data set the observed statistic is compared with the same statistic evaluated on null trials, and the printed
number is (k + 1)/(N_eff + 1), where k is the number of null trials at or above the observed value and N_eff the number
of statistically independent null realisations (Phipson & Smyth 2010 form; with k = 0 this is 1/(N_eff + 1)).
A 95 % binomial upper confidence bound with k = 0 would be about 3/N_eff (rule of three): ≈1.5-1.9 % per sector for
N_eff 160-200. The Note prints the p-value form, rounded up.

**Actual counts.**

| data set | observed S/N | k (of 20 000 trials) | author N_eff | author bound (k+1)/(N_eff+1) | reviewer N_eff (measured) | bound range |
|---|---|---|---|---|---|---|
| S6 | 9.837 | 0 | 164 | 0.61 % | 77 (IACT) - 100-140 (declustered tail) | 0.7-1.3 % |
| S32 | 8.982 | 0 | 195 | 0.51 % | 30 (IACT, inflated by segment-edge artefacts) - 140-300 (tail) | 0.3-3 % |
| joint | 13.218 | 0 | 3658 | 0.027 % | 1.7e3-1.6e4 | 0.006-0.06 % |

Printed in the Note: "≲1 % per sector, ≲0.1 % jointly" (the rounded-up upper ends of the reviewer's ranges).

**Where 160-200 comes from.** `s2_archival_fap.py`, `n_eff_periods()`: the number of distinct box trains inside the
trial-period range, (baseline / T14) × ln(P_max/P_min) with P in U(2.5, 4.5) d and T14 = 1.87 h: 164 for the S6
baseline, 195 for S32. This is a heuristic count. The pass-3 reviewer measured N_eff on a log-spaced period scan of
5 000 fake periods as (samples / integrated autocorrelation time) and as the number of declustered tail exceedances
(`reviewer_neff/x_neff2_main.json`, `reviewer_neff/x_neff2_clusters_declustered.json`); the S6 values are below 164,
which is why the printed bound was widened.

**Joint p-value.** The joint statistic is the inverse-variance-weighted box depth over S6 and S32 (S/N = Σ r/s² / √Σ 1/s²),
maximised over the epoch window. Its null is built from joint trials with an independent random anchor in each sector
(20 000 trials); N_eff(joint) = min(N_eff(S6), N_eff(S32)) × (number of independent relative S6-S32 phase positions inside
the joint window, 22 by the author's count, 33-65 by the reviewer's δ-scan). The joint bound is therefore a direct joint
p-value, not the product of the two per-sector bounds. A Fisher combination of the two per-sector p-values gives
0.035 % (3.4σ) as a cross-check.

**What this test is, and is not.** It is a *conditional* recovery test: given the 3.18 d period measured in Sector 98
alone, how often would noise in Sectors 6 and 32 produce an epoch-searched dip as strong as the observed one. It is not
the false-alarm probability of the original Sector 98 search over 126 246 stars (that search and its shortlist rule were
not changed and are not re-assessed here), and it does not by itself decide whether the periodic dips are planetary.
The fixed-epoch version of the test (no epoch search) is uninformative because the S98-only timing error propagated
to S6/S32 (18-25 h at 1σ) exceeds the transit duration: FAP 0.61 (S6), 0.25 (S32).

## Files

| file | content |
|---|---|
| `a_fap_results.json` | Assembled results: prediction, statistic, nulls, MAIN and OPT tests, detrending variants, recovered ephemeris, cycle count, momentum-dump check, caveats, first printable sentence. |
| `a_fap_main.json` | Full output of the main run (window 0.75 d, 20,000 trials, N2r 4,000). |
| `a_fap_edge0375.json`, `a_fap_w050.json`, `a_fap_w100.json` | Detrending variants (10,000 / 4,000 / 4,000 trials). |
| `s1_s98_prediction.json`, `s1_boot_P_Tref.npy`, `s1_s98_prediction.log` | Sector 98-only prediction, its bootstrap samples and log. |
| `s2_main.log`, `s2_edge0375.log`, `s2_w050.log`, `s2_w100.log` | Run logs of `s2_archival_fap.py`. |
| `a_fap_diagnostic.png` | Diagnostic figure (`s3_plot.py`). |
| `PASS3_NOTE.json` | The pass-3 reviewer finding and the printed wording. |
| `reviewer_neff/x_neff2_main.json`, `reviewer_neff/x_neff2_clusters_declustered.json` | Measured N_eff (see step 5). `reviewer_neff/x_neff2.py` is the script that wrote `x_neff2_main.json`. |
| `FACTS_IN_THIS_SECTION.csv` | The `A.archival_FAP` fact with its evidence mapped to package paths. |

## Scripts

| script | runs standalone from this folder? |
|---|---|
| `s4_results.py` | **Yes** (json, numpy): re-assembles `a_fap_results.json` from `s1_s98_prediction.json` and the four `a_fap_*.json` files. |
| `s1_s98_prediction.py`, `s2_archival_fap.py`, `s3_plot.py` | No: they read the three light curves through `fap_common.py`, which resolves the raw FITS path relative to the author's tree (`../../../../addenda/opus55_final_audit_tic4206066/00_data/mast/`). Download the three files listed in `../01_inputs/INPUT_PRODUCTS.csv` (TESS-SPOC S6 LC, SPOC S32 2-min LC, QLP S98 LC) and set `MAST` in `fap_common.py`. Command lines: `python s1_s98_prediction.py`; `python s2_archival_fap.py 0.75 10000 main 0.0` (arguments: window in days, number of fake periods, tag, edge cutoff). The main run used 20,000 trials. |
| `reviewer_neff/x_neff2.py` | No: it `exec`s the top of `s2_archival_fap.py` and needs the same raw data; arguments `<edge_cutoff> <tag>`. |

**Not included:** `_diag_main.npz` (0.6 MB diagnostic arrays), the reviewer's intermediate arrays and logs (`x_neff2_*_rows.npy`, `x_neff2_*_P0.npy`, `x_neff.py`, `x_neff_results.json`, duplicates of the files above). They are listed with sha256 in `../12_facts/EXCLUDED_FILES.csv`.
