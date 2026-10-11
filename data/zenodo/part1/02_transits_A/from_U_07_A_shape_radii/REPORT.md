# Audit 07_A_shape_radii — REPORT

> Rendered verbatim from the agent's structured result (`_structured/audit__07_A_shape_radii.json`); the sub-agent harness blocked agents from writing .md report files. Scripts, outputs and logs are in this folder.


**Verdict:** PASS_WITH_CHANGES

## Summary

I rebuilt the Sec. 2 numbers from the public data in 00_data/mast with my own code:
- quality masks: QUALITY==0; for QLP S98 this also drops the straylight 2048 bit, the 2^30 QLP flag and an undocumented 2^29 bit that brackets the straylight intervals;
- detrending: wotan biweight, 0.75 d window, both signals masked;
- limb darkening: Claret (2017) TESS coefficients from VizieR J/A+A/600/A30 (PHOENIX-q LSM u = 0.316, 0.313; ATLAS 0.375, 0.278);
- 16 joint batman+emcee fits, supersampled per cadence (S6 x15, S32 x2, S98 x4), all with N/tau >= 62.

The headline numbers are reproduced by one reasonable method, raw flux with per-transit linear baselines profiled in the likelihood:
- k = 0.02095 +- 0.0006, b = 0.33 (+0.09/-0.12), T14 = 1.88 +- 0.04 h, Rp = 1.42 +- 0.06 Re;
- with QLP DET_FLUX for S98: k = 0.0208;
- my own biweight light curves give only k = 0.0197 and Rp = 1.33 Re.

So the radius depends on how S98 is detrended at about the size of its quoted error. Most of the spread comes from S98 epoch -1, which sits at a segment edge. Red noise (beta = 1.32 in S32, 1.62 in S98) raises sigma_k to 0.0008.

The free-density limit rho* < 3.5 (97.7%) is neither prior- nor detrending-robust: I get 3.40 to 3.72.

Checks that hold:
- transit count: 23 transits have at least 50% of T14 covered; 24 have any in-transit data;
- per-sector depths agree (chi2 = 0.24/2 on raw flux, 2.9/2 on biweight);
- odd/even: 539 +- 47 vs 449 +- 51 ppm (1.3 sigma), close to the manuscript's 536/466, but the empirical errors are about 15% larger than quoted;
- no secondary eclipse: +12 +- 31 ppm flux at phase 0.5; the deepest dip in an all-phase scan is only 2.4-2.9 sigma after trials;
- V shape rejected: delta-chi2 = 34 with white errors, 21 with red-noise rescaling (the manuscript quotes 32);
- sqrt(depth) check: k = 0.0210 once the limb-darkening factor 1.14 is included;
- Signal B: 2.25 +- 0.12 Re and 2.01 Re are reproduced from k_B, but the abstract's "~2 Re" is the free-density value, not the adopted 2.25.

Problems with wording or labels:
- the phase-0.5 value has the wrong sign: its source is a depth of -20 ppm, i.e. flux +20;
- the "1.50 Re with free LD" value came from a fit with both LD and density free, and I could not reproduce it (my free-LD values are 1.35-1.45).

None of this is fatal. The fixes are wording plus wider error bars on Rp and on the rho limit.

The harness refused to write REPORT.md, so the full report content is in this structured result. Everything else is written in "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/07_A_shape_radii/":
- scripts: get_claret.py, common.py, s01_prep.py to s08_plot.py, run_batch*.sh;
- claret/ (VizieR tables, ReadMe, sha256 values);
- out/: fit_*.json and chain_*.npz, summary_fits.csv, tests_results.json, depths_localbaseline_mann_claret_jit_raw.json, detrend_variants.json, depth_check_local_baseline.json, beta_base.json, coverage_*.json, radii.json, secondary_scan_w*.txt, fold_and_secondary_scan.png;
- logs/.

No packages were installed. Documents cited:
- https://tess.mit.edu/qlp/ (QLP flag bits; SYS_RM_FLUX and DET_FLUX columns);
- the QLP data-product description in 00_data/docs;
- the VizieR J/A+A/600/A30 ReadMe;
- wotan 1.10 source, flatten.py (mask handling);
- arXiv:2603.22236 (QLP DRN 004), which has no 2^29 definition.

## Numeric rows

| quantity | manuscript | recomputed | status | method / note |
|---|---|---|---|---|
| Number of transits in joint fit | 23 | 23 with >=50% T14 coverage (S6 6, S32 8, S98 9); 24 with any in-transit point; 27 if the QLP 2^29 bit is kept | MATCH | s02_coverage.py on QUALITY==0 light curves  |
| Rp/R* (density prior, Claret LD) | 0.0212 +- 0.0007 | 0.02095 +- 0.0006 (raw + profiled local baselines); 0.0208 (DET_FLUX); 0.0207 (b29); 0.0197 +- 0.0006 (own biweight); 0.0200 +- 0.0008 (beta-rescaled) | MINOR_DIFF | emcee joint fits s03_jointfit.py detrending systematic ~0.001, comparable to the quoted error |
| impact parameter b | 0.35 +0.10/-0.14 | 0.33 +0.09/-0.12 (raw); 0.37 +0.08/-0.11 (biweight) | MATCH | joint fit  |
| T14 | 1.87 +- 0.05 h | 1.879 +- 0.041 h (raw); 1.854 +- 0.043 h (biweight); 1.890 +- 0.059 h (beta) | MATCH | joint fit derived  |
| Rp (density prior) | 1.44 +- 0.06 Re | 1.42 +- 0.06 (raw); 1.33 +- 0.06 (biweight); 1.36 +- 0.07 (beta); the manuscript's k propagates to 1.44 +- 0.065 | MINOR_DIFF | k x R*(N(0.623,0.019)) x 109.076  |
| Rp with free limb darkening | 1.50 Re | 1.45 (raw, density prior or uniform rho); 1.35 (biweight) | NOT_REPRODUCED | Kipping q1,q2 ~ U(0,1) the 1.50 source (earlier dynesty run) had density free as well |
| rho* 97.7% upper limit, density free | < 3.5 rho_sun | uniform rho: 3.47 raw / 3.65 biweight; log-uniform: 3.40 raw / 3.59 biweight; free LD + uniform: 3.59 raw / 3.72 biweight | MINOR_DIFF | emcee, rho prior variants  |
| per-sector depths agree | agree within errors | raw: 523+-69, 475+-77, 514+-45 ppm (chi2 0.24/2); biweight: 536+-58, 461+-74, 403+-53 ppm (chi2 2.9/2) | MATCH | template depth, empirical sigma from shifted ephemerides  |
| odd transit depth | 536 +- 40 ppm | 539 +- 47 ppm (raw local baseline); 453 +- 55 ppm (biweight) | MATCH | s06 / s04 template scaling with empirical sigma  |
| even transit depth | 466 +- 43 ppm | 449 +- 51 ppm (raw local baseline); 462 +- 52 ppm (biweight) | MATCH | as above odd-even differs by 1.3 sigma (raw) or 0.1 sigma (biweight) |
| flux at phase 0.5 | -20 +- 32 ppm | +12 +- 31 ppm (flux relative to local baseline; empirical sigma from random same-width boxes) | MINOR_DIFF | box of width T14, phase 0.5 masked in detrending source gives depth -20 ppm, i.e. flux +20: sign inverted in text |
| secondary scan over all phases (eccentric) | not quoted | deepest dip -66 to -99 ppm at 2.4-2.9 sigma_emp (widths 0.5-2 T14) | NOT_CHECKED | s04 scan no significant secondary eclipse |
| V-shape delta-chi2 | 32 | 34.2 (white errors); 21.2 (red-noise beta-rescaled) | MINOR_DIFF | supersampled trapezoid vs triangle  |
| transit depth (abstract) | ~500 ppm | mid-transit model depth 500 ppm (raw) / 440 ppm (biweight); k^2 = 439 ppm | MATCH | batman at mid-transit  |
| sqrt(depth) sanity check of k | (implicit) | 501 ppm / LD factor 1.141 gives k = 0.0210, Rp = 1.43 Re; naive sqrt(500 ppm) gives 1.52 Re | MATCH | LD-corrected sqrt(depth)  |
| Signal B Rp (density prior) | 2.25 +- 0.12 Re | 2.246 +- 0.121 Re (k 0.03305 +- 0.00147) | MATCH | k x R* Monte Carlo  |
| Signal B Rp (free density) | ~2.0 Re | 2.011 Re (k 0.0296) | MATCH | k x R*  |
| Signal B Rp in abstract | ~2 Re | the body's adopted value is 2.25 +- 0.12; ~2 is only the free-density value | MINOR_DIFF | consistency  |

## Findings

### F1 [MINOR] Rp/R* and Rp depend on S98 detrending at the level of the quoted error

The manuscript gives k = 0.0212 +- 0.0007 and Rp = 1.44 +- 0.06. These are reproduced only by some reasonable pipelines. Independent biweight detrending of the same data gives k = 0.0197 and Rp = 1.33. Red noise (beta 1.3-1.6) is not included in the quoted errors.

*Evidence:* Files: out/summary_fits.csv, out/detrend_variants.json, out/beta_base.json.
- k by method: mann_claret_jit_raw 0.02095 (Rp 1.42); det98 0.0208; b29 0.0207; base (biweight) 0.0197 (Rp 1.33); base without S98 epochs -1 and 8: 0.0203; beta-rescaled 0.0200 +- 0.0008.
- S98 epoch -1 (BTJD 3995.555, segment edge): biweight depth -234 ppm vs raw linear-baseline depth 898 ppm.

*Location:* Sec 2; abstract

*Proposed fix:* Quote Rp = 1.42 +- 0.08 Re, or state an additional detrending systematic of about 0.1 Re.

### F2 [MINOR] '1.50 Re with free limb darkening' is mislabelled and not reproduced

The 1.50 value comes from an earlier dynesty run in which both LD and density were free. My free-LD fits give 1.35-1.45 Re.

*Evidence:* - tic4206066_PUBLICATION_AUDIT/out/fits/dynesty_A_pipe.json: Rp median 1.503, rho median 1.59, b 0.62.
- My fits: mann_free_jit_raw 1.446; uni_free_jit_raw 1.447; mann_free_jit_base 1.349.

*Location:* Sec 2

*Proposed fix:* Write '(1.45-1.50 Re with free limb darkening and density)', or drop the parenthesis.

### F3 [MINOR] rho* < 3.5 rho_sun (97.7%) is not robust to prior or detrending

The 97.7% upper limit ranges from 3.40 to 3.72 across uniform and log-uniform density priors, two detrendings and free versus Claret-prior LD. The Sec. 4 host-mass argument (>0.54 Msun) uses this limit; with 3.7 the boundary moves to about 0.52 Msun, still earlier than M1V.

*Evidence:* From out/summary_fits.csv, rho97.7 by fit: uni_claret_jit_raw 3.47, logu_claret_jit_raw 3.40, uni_claret_jit_base 3.65, logu_claret_jit_base 3.59, uni_free_jit_raw 3.59, uni_free_jit_base 3.72. The earlier chain gave 3.47 (log-uniform) and 3.55 (reweighted to uniform).

*Location:* Sec 2 (and Sec 4)

*Proposed fix:* Write 'rho* <~ 3.6 rho_sun (97.7%; 3.4-3.7 depending on prior)' and check the 0.54 Msun statement.

### F4 [MINOR] Phase-0.5 value quoted with inverted sign

The text says 'flux at phase 0.5 is -20 +- 32 ppm'. The source value is a box depth of -20.4 ppm, i.e. flux +20 ppm. My measurement is a flux of +12 +- 31 ppm. Either way it is consistent with zero.

*Evidence:* - tic4206066_PUBLICATION_AUDIT/out/s10_eb_tests.json A_pipe.secondary_0p5.depth_ppm = -20.43.
- out/tests_results.json secondary.mask05.flux_ppm = 12.0, sig_emp_std_ppm = 31.3.

*Location:* Sec 2

*Proposed fix:* Write 'the depth at phase 0.5 is -20 +- 32 ppm (3 sigma < 95 ppm)', or use +12 +- 31 ppm flux.

### F5 [MINOR] V-shape delta-chi2 = 32 does not account for red noise

A supersampled trapezoid vs triangle gives delta-chi2 = 34.2 with white errors but 21.2 with beta-rescaled errors. The V shape is still strongly disfavored.

*Evidence:* out/tests_results.json vshape.white.dchi2 = 34.24; vshape.beta_scaled.dchi2 = 21.21; beta = 1.00, 1.32, 1.62 for S6, S32, S98.

*Location:* Sec 2

*Proposed fix:* Write 'disfavored (delta-chi2 ~ 21 after red-noise rescaling)'.

### F6 [MINOR] Abstract '~2 Re' for the 11.13 d signal is inconsistent with the adopted 2.25 +- 0.12

2.25 +- 0.12 (density prior) and ~2.0 (free density) are both reproduced. The abstract quotes only the free-density value.

*Evidence:* out/radii.json: 0.03305 x 0.623 x 109.076 = 2.246 +- 0.121; 0.0296 gives 2.011.

*Location:* Abstract

*Proposed fix:* Use '~2.2 Re' in the abstract.

### F7 [NOTE] S98 mask drops an undocumented QLP 2^29 quality bit

QUALITY==0 removes the 2^29 bit (10,676 of 24,162 cadences, bracketing the straylight intervals). I found no definition of this bit on the QLP pages or in DRN 004. Keeping these cadences adds 4 transits (27 total) and gives k = 0.0207, consistent.

*Evidence:* out/coverage_base.log vs out/coverage_b29.log; fit_mann_claret_jit_b29.json k median 0.02068.

*Location:* 

*Proposed fix:* Say 'cadences flagged by QLP (including straylight) removed'.

### F8 [NOTE] 'catalog density prior' wording; LD prior centre

The density prior comes from the Mann relations, not a catalog. The earlier fit centred LD at Teff 4030 K (0.325, 0.307) rather than 4012 K (0.316, 0.313); the effect is negligible.

*Evidence:* claret/tableab.ecsv interpolation in common.claret_ld(); tic4206066_PUBLICATION_AUDIT/scripts/s08_globalfit.py TEFF = 4030.

*Location:* 

*Proposed fix:* Write 'Mann-relation density prior'.

### F9 [NOTE] Odd/even errors slightly underestimated; result depends on method

On raw flux with local baselines, odd/even = 539 +- 47 vs 449 +- 51 ppm (1.3 sigma), matching the manuscript's values but with about 15% larger empirical errors. The biweight light curves give 453 vs 462 ppm. The earlier box method gave 534 vs 388 ppm (2.2 sigma). None of these indicates an EB.

*Evidence:* out/depths_localbaseline_mann_claret_jit_raw.json; out/tests_results.json depths.odd/even; s10_eb_tests.json odd_even_box.

*Location:* 

*Proposed fix:* Optionally quote 539 +- 47 and 449 +- 51 ppm.

## Recommended wording

Sec. 2: "A joint emcee fit of all 23 transits (circular orbit; quadratic limb darkening with Claret 2017 priors; Mann-relation density prior) gives Rp/R* = 0.0210 +- 0.0008, b = 0.33 (+0.09/-0.12), T14 = 1.88 +- 0.05 h, and Rp = 1.42 +- 0.08 Re (1.33-1.45 Re across detrending and limb-darkening choices). With the density free, rho* <~ 3.6 rho_sun (97.7%; 3.4-3.7 depending on prior). Per-sector depths agree within errors, odd and even transits give 539 +- 47 and 449 +- 51 ppm, the depth at phase 0.5 is -20 +- 32 ppm (consistent with zero), and a V-shaped transit is disfavored (delta-chi2 = 21 after red-noise rescaling)."

Abstract: "(23 transits, depth ~500 ppm) corresponds to a 1.4 Re planet" and "The 11.13 d signal (~2.2 Re)".

## Files

- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/07_A_shape_radii/common.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/07_A_shape_radii/get_claret.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/07_A_shape_radii/s01_prep.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/07_A_shape_radii/s02_coverage.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/07_A_shape_radii/s03_jointfit.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/07_A_shape_radii/s04_tests.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/07_A_shape_radii/s05_depth_check.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/07_A_shape_radii/s05b_detrend_variants.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/07_A_shape_radii/s06_depths_localbaseline.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/07_A_shape_radii/s07_radii_summary.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/07_A_shape_radii/s08_plot.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/07_A_shape_radii/claret/
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/07_A_shape_radii/out/summary_fits.csv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/07_A_shape_radii/out/tests_results.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/07_A_shape_radii/out/depths_localbaseline_mann_claret_jit_raw.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/07_A_shape_radii/out/detrend_variants.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/07_A_shape_radii/out/radii.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/07_A_shape_radii/out/beta_base.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/07_A_shape_radii/out/fold_and_secondary_scan.png
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/07_A_shape_radii/out/fit_*.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/07_A_shape_radii/out/chain_*.npz
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/07_A_shape_radii/logs/