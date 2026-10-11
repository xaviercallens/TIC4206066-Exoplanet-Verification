# Audit 06_stellar — REPORT

> Rendered verbatim from the agent's structured result (`_structured/audit__06_stellar.json`); the sub-agent harness blocked agents from writing .md report files. Scripts, outputs and logs are in this folder.


**Verdict:** PASS_WITH_CHANGES

## Summary

Verdict: PASS WITH CHANGES. I could reproduce every stellar number in the manuscript from fresh public data with my own code. Nothing is wrong in a way that affects a conclusion. Two wording issues need a fix before submission (both MINOR): the ±0.07 mag is credited to Way et al. (2026), which does not give it, and the ARIADNE error bars are narrower than the same run supports.

REPORT.md was not written: the harness blocked the Write tool for .md report files ("Subagents should return findings as text"). The report content is in this result. The numeric table is in out/AUDIT06_numeric_rows.csv.

DATA. All inputs come from fresh queries (UTC time and sha256 for each file in data/PROVENANCE.jsonl). The shared 00_data folder was still empty when I started, so I queried everything myself.
- Gaia DR3 via the Gaia TAP service: gaia_source, astrophysical_parameters, the XP sampled spectrum, and two 100-pc comparison samples (18,416 stars; one joined to 2MASS and AllWISE through the Gaia best-neighbour tables).
- Bailer-Jones+2021 distances.
- TIC 8.2 from MAST.
- VizieR: 2MASS, AllWISE, CatWISE2020, Tycho-2, APASS DR9.
- SIMBAD identifiers.
- Way+2026 catalogues from the Zenodo API. The md5 checksums match Zenodo; the files are kept in the session scratchpad, not in the project.
- Primary sources for every coefficient: arXiv sources of Mann+2015 v3 and Mann+2019, the official mk-mass code (commit ffe5f19), the colte code for Casagrande+2021, and the Pecaut & Mamajek table v2022.04.16.
- Added package: gaiadr3-zeropoint 0.1.0 (uv pip install into audit_env).

RESULTS.
- Distance: 1/parallax = 35.619 pc; with the Lindegren zero point (-0.0332 mas) 35.577 pc; Bailer-Jones geometric 35.579 pc. All give 35.6. MATCH.
- T: TIC value 9.6817, reproduced exactly by the Stassun+2019 relation. MATCH.
- M_Ks = 4.938 ± 0.029. Extinction is negligible: Way+2026 dereddened with a 3D dust map and list M_G 7.755, which equals my extinction-free 7.7546, and TIC E(B-V) = 0.
- Mann+2015 radius: 0.6229 ± 0.0188. The coefficients match the arXiv v3 source and the 2016 erratum; the error is the 2.89% scatter combined with the M_K error. MATCH.
- Mann+2019 mass: 0.6062 ± 0.0151 from Table 6 plus the 2.4% error in Table 7 plus the M_K error. MATCH. The official posterior code, which Mann et al. call more robust, gives 0.6109 ± 0.0153 (NOTE).
- Density: 2.507 +0.245/-0.218 ρ☉. MATCH (the +0.25 is a rounding of 0.245).
- Validity ranges: M_K 4.94 is inside Mann+2015 (4.6-9.8) and the safer Mann+2019 range (4.5-10.5). Teff (3990-4070 K by every route) is just below the Mann+2015 upper limit of 4100 K. A slightly subsolar [Fe/H] changes the radius by -1 to -1.6%.
- ARIADNE: the manuscript copies best_fit_average.dat exactly (0.6206 +0.0147/-0.0140; 4011.8 +45.3/-40.3). That output is a weighted sum of independent draws, which narrows the spread. The correct mixture of the per-grid posteriors from the same BMA.pkl gives R 0.620 +0.020/-0.018 and Teff 4017 +58/-64. Per-grid Teff runs from 3939 to 4038 K.
- Independent SED check: F_bol = 2.234e-9 erg/s/cm² from the Gaia XP spectrum plus 2MASS and WISE, so L = 0.0884 ± 0.0019 L☉. Mann+2015 and Pecaut & Mamajek bolometric corrections give 0.0874-0.0885.
- Independent Teff: Mann colour relations 3987-4047 ± 85-99 K; Casagrande colte 4068 ± 60 K (at the cool edge of its calibration); Pecaut & Mamajek 3981/4031 K.
- Stefan-Boltzmann radius: 0.60-0.62 R☉ (0.613 ± 0.023 for Teff 4020 K). The radius is confirmed independently.
- Infrared excess: compared with 100-pc stars of the same BP-RP, the offsets are Ks-W1 -0.012, Ks-W2 +0.012, Ks-W3 -0.027 ± 0.046 and Ks-W4 +0.029 ± 0.14. There is no excess to 22 µm; at 2σ, excesses above about 6% in W3 and 32% in W4 are ruled out.
- Way+2026 exists: AJ 171, 252, doi 10.3847/1538-3881/ae48ae (confirmed by Crossref); Zenodo record 18500082. The target row gives O = +0.003, not flagged, one component. The paper gives no per-star error. The 0.07 is the author's own scatter of O (1.4826×MAD) over 5205 comparable stars. I get 0.068. Other scatter definitions give 0.052-0.098.
- Companion limit: at 2σ the G flux ratio must be below 0.137, which means M < 0.380 M☉ from the Pecaut & Mamajek table (0.373 via the Mann+2019 relation). It is 0.35-0.42 M☉ depending on the scatter used. The paper's own flag threshold corresponds to 0.51 M☉.
- Independent CMD check: the star is +0.006 mag above a clipped local main-sequence fit (scatter 0.15 mag) and +0.03 vs Pecaut & Mamajek, i.e. on the single-star sequence.
- Circular-orbit host threshold: 0.53 M☉ (manuscript says 0.54).

WHICH COMPANIONS REMAIN ALLOWED. The data do not show that the star is single.
- Blended companions (Gaia-unresolved, ≲20 au ≈ 0.6″) up to about 0.35-0.42 M☉: all M3V and later dwarfs, brown dwarfs and white dwarfs.
- Faint companions between about 0.6″ and 2″, which the XP test does not cover. Gaia DR3 lists no other source within 11.1″, and ipd_frac_multi_peak = 0.
- The Gaia RV epochs are not strictly constant: p = 4.1e-5, renormalised goodness of fit 3.0, just under the Katz+2023 variability threshold of 4. This matters for the manuscript's RV-based exclusion and should be quantified by audit 09.
- An allowed M3V companion would bias the primary's M_Ks by 0.22 mag, which lowers R* by 6.5% and M* by 5.5%.

EARLIER-PACKAGE ERRORS (none reach the manuscript).
- tic4206066_PUBLICATION_AUDIT/scripts/s04_stellar.py has the wrong sign on the cubic bolometric-correction terms, giving wrong L (0.101/0.107) and R (0.666) in out/s04_stellar.json. They must not be reused.
- companion_check SUMMARY.md F4 says ARIADNE used DR2 photometry. It actually used DR3 magnitudes under DR2 passband names.

## Numeric rows

| quantity | manuscript | recomputed | status | method / note |
|---|---|---|---|---|
| Identity: Gaia DR3 3220388198192519424 = StKM 1-561 = TIC 4206066 | Gaia DR3 3220388198192519424 (StKM 1-561) | SIMBAD lists TIC 4206066, StKM 1-561, Gaia DR3 3220388198192519424, 2MASS J05232636-0119380 as the same object; SIMBAD SpT K6 | MATCH | SIMBAD query_objectids; Gaia TAP data/simbad_ids_StKM1-561.ecsv |
| Spectral class 'late-K dwarf' | late-K dwarf | BP-RP 1.739 and M_G 7.755 fall between K8V and K9V (Pecaut & Mamajek); Teff 3990-4070 K; SIMBAD K6 | MATCH | Pecaut & Mamajek table interpolation; colour-Teff relations out/s04_misc.json |
| Distance | 35.6 pc | 35.619 ± 0.022 (1/parallax, parallax 28.0750 ± 0.0176 mas); 35.577 (Lindegren+21 zero point -0.0332 mas); BJ21 geometric 35.579 [35.557, 35.596]; photogeometric 35.574 | MATCH | 1/parallax with and without gaiadr3-zeropoint 0.1.0; Bailer-Jones+2021 out/s01_relations.json |
| TESS magnitude T | 9.68 | 9.6817 ± 0.006 (TIC 8.2); reproduced to 1e-4 mag by the Stassun+2019 relation from the DR2 photometry stored in the TIC; 9.685 from DR3 photometry | MATCH | MAST TIC query; Stassun+2019 relation out/s04_misc.json |
| M_Ks | 4.94 | 4.938 ± 0.029 (zero-point-corrected parallax); 4.936 (uncorrected); 4.938 (Bailer-Jones distance) | MATCH | 2MASS Ks 7.694 ± 0.029 (quality AAA) plus distance modulus; A_Ks ≈ 0 Extinction ≈ 0: Way+2026 dereddened M_G 7.755 equals the extinction-free 7.7546 |
| R* (Mann+2015) | 0.623 ± 0.019 Rsun | 0.6229 ± 0.0188 (Monte Carlo, N=1.6e6; 2.89% scatter plus M_K error); with the [Fe/H] term: 0.618 ([Fe/H]=-0.08) / 0.613 (-0.27) | MATCH | Mann+2015 Table 1 coefficients from the arXiv v3 source (same as the 2016 erratum) M_K inside 4.6-9.8; Teff just below the 4100 K upper limit |
| M* (Mann+2019) | 0.606 ± 0.015 Msun | 0.6062 ± 0.0151 (Table 6 n=5 best fit, Table 7 2.4% error, plus M_K error); 0.6109 ± 0.0153 with the official posterior code | MATCH | Mann+2019 Table 6/7; posterior file Mk-M_7_trim.fits from github awmann/M_-M_K- (commit ffe5f19) The best fit sits 0.8% below the posterior median at this M_K |
| rho* | 2.51 +0.25/-0.22 rho_sun | 2.507 +0.245/-0.218 (Table 6 mass); 2.527 +0.247/-0.219 (official posterior mass); log g 4.632 ± 0.028 | MATCH | Monte Carlo of M/R^3 with shared M_K draws and independent relation scatters out/s01_relations.json |
| ARIADNE SED R* | 0.621 ± 0.015 Rsun | File value reproduced: 0.6206 +0.0147/-0.0140 (best_fit_average.dat). Correct model mixture from the same run: 0.620 +0.020/-0.018 (BMA.pkl, percentiles) or 0.620 +0.033/-0.030 (best_fit_sample.dat). Independent Stefan-Boltzmann radius: 0.613 ± 0.023 | MINOR_DIFF | Read the ARIADNE outputs and redraw the per-grid posteriors by weight; independent F_bol and Teff The central value is confirmed; the quoted error is too small |
| ARIADNE SED Teff | 4012 ± 43 K | File value reproduced: 4011.8 +45.3/-40.3. Mixture: 4017 +58/-64 (pickle) or 3989 +123/-78 (best_fit_sample.dat). Per-grid 3939-4038 K. Independent colour relations 3981-4068 K | MINOR_DIFF | As above; Mann+2015 Table 2 (with 60 K added); Casagrande+2021 colte; Pecaut & Mamajek Tayar+2022 floor of 2.4% = 96 K |
| No infrared excess | no infrared excess | Offsets from 100-pc stars of the same BP-RP (±0.03): Ks-W1 -0.012±0.04, Ks-W2 +0.012±0.03, Ks-W3 -0.027±0.046, Ks-W4 +0.029±0.14. 2σ excess limits: about 7% (W1/W2), 6% (W3), 32% (W4) | MATCH | Empirical locus from Gaia DR3 matched to 2MASS and AllWISE (best-neighbour tables) The ARIADNE fit only reached W2 (4.6 µm) and fitted free per-band jitter |
| RUWE; non-single-star solution | RUWE = 1.21; no NSS solution | ruwe 1.2143; non_single_star 0; ipd_frac_multi_peak 0; ipd_gof_harmonic_amplitude 0.022 | MATCH | Gaia DR3 gaia_source Gaia RV epochs: rv_chisq_pvalue 4.1e-5, rv_renormalised_gof 3.0 (see F8) |
| Way+2026 overluminosity | +0.003 mag | +0.003 (predicted M_G 7.758 minus measured 7.755); not flagged; systems catalogue: 1 component, no wide-binary catalogue membership | MATCH | Zenodo 10.5281/zenodo.18500082 v1 (md5 verified); AJ 171, 252 confirmed via Crossref out/way2026_target_row.csv, out/way2026_systems_row.csv |
| Overluminosity uncertainty ±0.07 (credited to Way+2026) | ±0.07 mag | Not in Way+2026: the paper and catalogue give no per-star error. It is the author's scatter (1.4826×MAD) of O among 5205 comparable catalogue stars; I get 0.068. Clipped core 0.052; clean-astrometry subset 0.053; 16-84% half-width 0.098 | DISCREPANT | Recomputed from the catalogue The attribution is wrong; the value itself is reproducible |
| Unresolved-companion mass limit | above 0.38 Msun excluded (2σ) | 0.380 Msun (σ=0.068, Pecaut & Mamajek masses) / 0.373 (via the Mann+2019 relation); G flux ratio < 0.137. 0.35 (σ≈0.052) to 0.42 (σ=0.10); 3σ 0.42-0.43; the paper's flag threshold (0.376 mag) corresponds to 0.51 | MATCH | Flux-sum conversion; Pecaut & Mamajek M_G-mass relation; Mann+2019 Applies only to companions blended in Gaia (≲20 au ≈ 0.6″) |
| Circular-orbit host mass threshold (rho ≤ 3.5 rho_sun) | above 0.54 Msun | 0.529 Msun (Mann+2015 R × Mann+2019 M, M_K = 5.45); about 0.53 from the Pecaut & Mamajek table (M0.5V 0.54 at 3.35; M1V 0.50 at 3.98) | MINOR_DIFF | Evaluate the relations at rho = 3.5 rho_sun Rounding only; no effect because 0.38 ≪ 0.53 |
| Bolometric luminosity | not stated (ARIADNE file: 0.090) | 0.0884 ± 0.0019 Lsun (integration of XP spectrum + 2MASS + WISE; F_bol 2.234e-9 erg/s/cm2); 0.0874-0.0885 (Mann+2015 bolometric corrections for V, J, K, r); 0.0881 (Pecaut & Mamajek) | MATCH | Independent F_bol Agrees with ARIADNE within 2% |
| Independent Gaia CMD position | not stated | +0.006 mag above a clipped local main-sequence fit (100 pc, clean astrometry, scatter 0.15 mag); +0.02 to +0.03 in narrow colour windows; +0.031 vs Pecaut & Mamajek | MATCH | Gaia DR3 100-pc sample queried via TAP The star lies on the single-star sequence |

## Findings

### F1 [MINOR] The ±0.07 mag overluminosity uncertainty is credited to Way et al. (2026) but comes from the author

The manuscript writes 'The Gaia XP overluminosity of Way2026, +0.003±0.07 mag'. Way et al. (2026) and their Zenodo catalogue give no per-star uncertainty and no numerical scatter for their predicted M_G. The 0.07 is the author's own robust scatter of O among catalogue stars with similar predicted M_G and colour. The value is reproducible (0.068) but depends on the definition: 0.052 (clipped core), 0.053 (clean-astrometry subset), 0.098 (16-84% half-width).

*Evidence:* out/s03_overluminosity.json, local_O_distribution: 1.4826×MAD = 0.0682 (n=5205), clipped σ 0.0518, near_clean_astrometry σ_MAD 0.0534 (n=2787). The origin is companion_check_tic4206066/overluminosity.py (sig_loc). Way+2026 (arXiv:2602.20096 HTML) states no per-star error or scatter value; the catalogue has 7 columns and no error column (Zenodo record JSON).

*Location:* Sect. 4 (Localization and False Positives), sentence on Way2026

*Proposed fix:* Write: 'The Gaia XP-based overluminosity from Way et al. (2026) is +0.003 mag; compared with the 0.07 mag scatter of stars with similar XP spectra, this excludes at 2σ an unresolved (≲20 au) companion contributing more than 14% of the G-band flux, i.e. above ≈0.38 M☉.'

### F2 [MINOR] The ARIADNE error bars (±0.015 Rsun, ±43 K) are too small

The quoted values come from ARIADNE's best_fit_average.dat. ARIADNE builds that 'weighted average' posterior as a weighted sum of independent draws from each grid's posterior, which narrows the spread by about sqrt(sum w²) = 0.52. The proper model-averaged posterior is the weighted mixture. From the same run it gives R = 0.620 +0.020/-0.018 Rsun and Teff = 4017 +58/-64 K; ARIADNE's own best_fit_sample.dat gives 0.620 +0.033/-0.030 and 3989 +123/-78 K. These errors are also formal only; the Tayar+2022 floors (4.2% in R, 2.4% in Teff) add about 0.026 Rsun and 96 K. The central values are robust: independent L = 0.0884 ± 0.0019 Lsun, Teff 3987-4068 K and Stefan-Boltzmann R = 0.613 ± 0.023 Rsun.

*Evidence:* companion_check_tic4206066/sed_out/best_fit_average.dat (rad 0.6206 +0.0147/-0.0140; teff 4011.8 +45.3/-40.3) vs best_fit_sample.dat (0.6203 +0.0329/-0.0296; 3988.7 +123.0/-78.0). BMA.pkl arrays: weighted_average R std 0.0087 vs mixture 0.019; Teff std 26 vs 58 K; weights 0.150/0.311/0.292/0.246; per-grid Teff 3939-4038 K, R 0.615-0.642. ARIADNE fitter.py bayesian_model_average: weighted_average += choice(o[k], nsamples)*weights[i]. Summary of my check in out/s02_sed.json, ariadne_bma_pickle.

*Location:* Sect. 1, 'a broadband SED fit with ARIADNE ... gives R*=0.621±0.015 Rsun, Teff=4012±43 K'

*Proposed fix:* Quote 'R* = 0.62 ± 0.02 Rsun and Teff = 4020 ± 60 K (model-averaged, formal errors)', or with systematic floors 'R* = 0.62 ± 0.03 Rsun, Teff = 4020 ± 110 K'.

### F3 [MINOR] The 0.38 Msun (2σ) limit depends on the scatter definition and holds only for companions blended in Gaia

The limit is reproduced: 0.380 Msun with Pecaut & Mamajek masses and 0.373 Msun via the Mann+2019 relation, for σ = 0.068. It moves to 0.35 Msun for the clipped or clean-single-star scatter (about 0.052) and to 0.42 Msun for σ = 0.10. Way et al. describe the method as sensitive to close (≲20 au) near-equal-mass pairs, which is ≲0.6″ at 35.6 pc. Their Table 2 predictions come from a final model whose sample overlaps its training stars. The manuscript's 'unresolved companions' is acceptable, but the separation range and the dependence on the adopted scatter are not stated.

*Evidence:* out/s03_overluminosity.json, companion_limits: sigma0.068 2σ: f_G < 0.137, M_G2 > 9.92, 0.380 (PM) / 0.373 (Mann); sigma_clipped 0.3487; sigma_0.10 0.4205; paper threshold 0.376 mag corresponds to 0.511 Msun. Way+2026 HTML: sensitive to 'close (≲20 AU), equal-mass systems'.

*Location:* Sect. 4, 'excludes unresolved companions above 0.38 Msun (2σ)'

*Proposed fix:* Add '(≲20 au)' after 'unresolved' and '≈' before 0.38. Optionally mention the range 0.35-0.42 Msun for different scatter estimates.

### F4 [NOTE] The official Mann+2019 posterior code gives 0.611 Msun, not 0.606

The manuscript uses the Table 6 best-fit coefficients with the Table 7 error. Mann+2019 (Sect. 5.1 and README) recommend the full posterior as 'a more robust method'. It gives 0.6109 ± 0.0153 Msun and rho = 2.53 +0.25/-0.22. The difference is 0.3σ and the manuscript's choice is legitimate.

*Evidence:* out/s01_relations.json: table6_plus_table7 0.6062 ± 0.0151; official_posterior_MC 0.6109 ± 0.0153; posterior relation at M_K0 without input errors 0.6110 ± 0.0081. Table 7 lists 0.6020 at M_K=5.0 where the Table 6 formula gives 0.5972.

*Location:* Sect. 1, M* = 0.606 ± 0.015

*Proposed fix:* Optional: quote 0.611 ± 0.015 Msun (posterior code), or keep 0.606 and state that the Table 6 coefficients were used.

### F5 [NOTE] 'No infrared excess' is correct, but the ARIADNE fit it rests on only reached 4.6 µm

ARIADNE fitted bands up to W2 only and included free per-band jitter terms of 0.09-0.13 mag in 2MASS and WISE, which can absorb an excess of about 10%. An independent comparison with 100-pc stars of the same BP-RP confirms there is no excess in W1-W4.

*Evidence:* out/s02_sed.json, ir_excess_empirical_locus: Ks-W3 target 0.134 ± 0.035 vs locus 0.161 (σ 0.046, n=1186); Ks-W4 0.243 ± 0.136 vs 0.214 (σ 0.080, n=47); W4 S/N 8.2. sed_out/best_fit_average.dat: Ks_noise 0.13, W1_noise 0.09, W2_noise 0.11.

*Location:* Sect. 1, 'and no infrared excess'

*Proposed fix:* Optional: 'and the 2MASS and WISE W1-W4 colours show no infrared excess'.

### F6 [NOTE] The star is inside the Mann relations' validity ranges, near the Teff edge

M_K = 4.94 is inside Mann+2015 (4.6-9.8) and the safer Mann+2019 range (4.5-10.5). Teff (3990-4070 K by every route) is just below the Mann+2015 upper limit of 4100 K. Mann+2015's late-K calibrators are mostly metal-rich; with [Fe/H] of -0.08 or -0.27 the [Fe/H]-dependent relation gives R = 0.618 / 0.613 (-0.8% / -1.6%, under 0.6σ).

*Evidence:* refs/m15/main.tex l.251 (validity ranges); refs/m19/Masses_1.tex l.1220; out/s01_relations.json with_FeH variants; out/s02_sed.json teff.

*Location:* Sect. 1

*Proposed fix:* None required.

### F7 [NOTE] The circular-orbit host threshold is 0.53 Msun, not 0.54

rho = 3.5 rho_sun is reached at M_K = 5.45, M = 0.529 Msun with the Mann relations, and at about 0.53 Msun by interpolating the Pecaut & Mamajek table. This is rounding and has no effect, since the 0.38 Msun limit is far below it.

*Evidence:* out/s04_misc.json host_density_threshold_mann and PM_density_by_SpT.

*Location:* Sect. 4, 'would need a host above 0.54 Msun'

*Proposed fix:* Optional: write '≳0.53 Msun'.

### F8 [NOTE] Gaia RV epochs are not strictly constant; allowed companions

Gaia DR3 gives rv_chisq_pvalue 4.1e-5 and rv_renormalised_gof 3.0 over 33 transits (rv_amplitude_robust 3.2 km/s). This is not formally RV-variable under Katz+2023 (which needs gof > 4), but it is not consistent with a constant RV either. The data still allow: blended companions ≲0.35-0.42 Msun (M3V and later, brown dwarfs, white dwarfs) within about 20 au; faint companions between about 0.6″ and 2″, where the XP test is uncalibrated and Gaia completeness is limited (Gaia lists nothing within 11.1″); and low-mass spectroscopic companions. An allowed M3V companion would bias M_Ks by 0.22 mag, lowering R* by 6.5% and M* by 5.5%. The manuscript does not claim the star is single, which is correct.

*Evidence:* data/gaia_dr3_source.ecsv; data/gaia_dr3_neighbours_60as.ecsv (nearest source 11.1″, G=18.7, parallax 0.22 mas); out/s03_overluminosity.json companion_grid_PM.

*Location:* Sect. 4, 'radial velocities that exclude an eclipsing stellar companion'

*Proposed fix:* Pass to audit 09 to quantify the RV-based exclusion given p = 4e-5.

### F9 [NOTE] Errors in the earlier package that do not reach the manuscript

(a) tic4206066_PUBLICATION_AUDIT/scripts/s04_stellar.py uses -4.084e-3 and -6.263e-3 for the cubic terms of the V- and K-band bolometric corrections; the correct sign is + (arXiv v3 and the 2016 erratum). This gives L_from_BC_V 0.101, L_from_BC_K 0.107 and radius_stefan_boltzmann 0.666 in out/s04_stellar.json; the correct values are 0.088 and about 0.60-0.62. (b) companion_check SUMMARY.md F4 says ARIADNE used DR2 photometry; star_inputs.json shows DR3 magnitudes (10.512974/11.353702/9.614222) under DR2 passband names (ARIADNE librarian _filtermap.py), with a global +0.037 mas parallax offset. (c) references.bib cites Mann+2015 without the 2016 erratum; the coefficients used are the corrected ones.

*Evidence:* tic4206066_PUBLICATION_AUDIT/out/s04_stellar.json; refs/m15/main.tex Table 3; companion_check_tic4206066/sed_out/star_inputs.json; ARIADNE librarian/_filtermap.py and _api.py (_PLX_ZEROPOINT = 0.037).

*Location:* Not in the manuscript (supporting package)

*Proposed fix:* Do not reuse those numbers; correct SUMMARY F4; optionally cite the 2016 erratum.

## Recommended wording

Sect. 1: "TIC 4206066 (Gaia DR3 3220388198192519424; T = 9.68) lies at 35.6 pc (Gaia DR3). From M_Ks = 4.94, the relations of Mann et al. (2015, 2019) give R* = 0.623 ± 0.019 R☉ and M* = 0.606 ± 0.015 M☉, hence ρ* = 2.51^{+0.25}_{-0.22} ρ☉; a broadband SED fit with ARIADNE (Vines & Jenkins 2022) gives R* = 0.62 ± 0.02 R☉ and T_eff = 4020 ± 60 K (model-averaged, formal errors), and the 2MASS and WISE W1–W4 colours show no infrared excess."

Optional variants: with the Tayar et al. (2022) floors, "R* = 0.62 ± 0.03 R☉, T_eff = 4020 ± 110 K". With the official Mann et al. (2019) posterior code, "M* = 0.611 ± 0.015 M☉, ρ* = 2.53^{+0.25}_{-0.22} ρ☉".

Sect. 4: "The Gaia XP-based overluminosity from Way et al. (2026) is +0.003 mag; compared with the 0.07 mag scatter of stars with similar XP spectra, this excludes at 2σ an unresolved (≲20 au) companion contributing more than 14% of the G-band flux, i.e. above ≈0.38 M☉. A planet on a bound companion in a circular orbit would need a host above ≈0.53 M☉ ..."

## Files

- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/q01_fetch.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/s01_relations.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/s02_sed.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/s03_overluminosity.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/s04_misc.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/s05_plots.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/data/PROVENANCE.jsonl
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/data/gaia_dr3_source.ecsv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/data/gaia_dr3_ap.ecsv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/data/gaia_dr3_neighbours_60as.ecsv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/data/gaia_dr3_xp_sampled.csv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/data/bailerjones2021.ecsv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/data/bj2021_vizier.ecsv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/data/tic82.ecsv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/data/2mass_psc.ecsv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/data/allwise.ecsv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/data/catwise2020.ecsv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/data/tycho2.ecsv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/data/apass9.ecsv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/data/gaia_100pc_bprp_1.54_1.94.ecsv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/data/gaia_100pc_bprp_1.64_1.84_2mass_allwise.ecsv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/data/zenodo_18500082_record.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/data/crossref_way2026.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/data/simbad_ids_StKM1-561.ecsv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/out/s01_relations.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/out/s02_sed.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/out/s02_colte_casagrande2021.csv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/out/s03_overluminosity.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/out/s04_misc.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/out/pecaut_mamajek_parsed.csv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/out/way2026_target_row.csv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/out/way2026_systems_row.csv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/out/AUDIT06_numeric_rows.csv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/out/s03_cmd_overluminosity.png
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/out/s05_sed_ir.png
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/logs/
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/refs/ (Mann+2015/2019 arXiv sources, mann2019_repo, colte_repo, EEM_dwarf_UBVIJHK_colors_Teff.txt)
- /private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/06_stellar_dl/ (Way+2026 Zenodo CSVs, ARIADNE source copies)