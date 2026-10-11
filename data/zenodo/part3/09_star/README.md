# 09_star: stellar parameters, Gaia RUWE and RV, rotation and activity

`FACTS_IN_THIS_SECTION.csv` lists every `STAR.*` fact: value, uncertainty, printed form, and evidence mapped to package paths.

| printed in the Note | value | where in the package |
|---|---|---|
| T = 9.68 | 9.6817 | `../05_companion_limit_038/from_U_06_stellar/data/tic82.ecsv` (Tmag); `out/s04_misc.json` |
| 35.6 pc | 35.577 pc (1/zero-point-corrected parallax 28.1083 mas); Bailer-Jones 35.579 | `from_T6_star/out/t6_05_parallax_zp.json` (script `t6_05_parallax_zp.py`); `../05_companion_limit_038/from_U_06_stellar/out/s01_relations.json` |
| M_Ks = 4.94 | 4.938 ± 0.029 | `../05_companion_limit_038/from_U_06_stellar/out/s01_relations.json` → `MK` |
| R = 0.623 ± 0.019, M = 0.606 ± 0.015, ρ = 2.51 +0.25/−0.22 | Mann et al. 2015 (R), Mann et al. 2019 (M), Monte Carlo | `from_T6_star/out/t6_04_density_mc.json` (script `t6_04_density_mc.py`); `../05_companion_limit_038/from_U_06_stellar/out/s01_relations.json` (script `s01_relations.py`) |
| RUWE = 1.21 | 1.2143 (DR3). Percentile among comparable stars: 61.5 (all), 74.4 (apparently single), 86.5-90 (RUWE < 1.4) | `../05_companion_limit_038/from_U_06_stellar/data/gaia_dr3_source.ecsv` (`ruwe`); `from_T6_star/out/t6_02_ruwe_pct.json` with sample `data/ruwe_sample_fresh.csv` (script `t6_02_ruwe_pct.py`) |
| RV scatter p = 4 × 10^-5 | `rv_chisq_pvalue` 4.06e-5, `rv_renormalised_gof` 3.0, 33 transits; below the DR3 RV-variability rule (p < 0.01 **and** gof > 4 **and** N ≥ 10) | `gaia_dr3_source.ecsv` (`rv_*`); context in `from_T6_star/out/t6_02_ruwe_pct.json` → `rv_context`, `from_T6_star/verify/v01_ruwe_rv_check.json`; deblending check in `from_extra_science_rv_check/rv01_deblend.py` / `.json` (Gaia TAP sample `rv01_sample.csv`) |
| RVs exclude an eclipsing stellar companion of the target | 95% K < 1.16 km/s vs 14.9 km/s for 0.08 Msun at 3.18 d | `../05_companion_limit_038/from_U_09_companion/out/rv_limits.json`; `../04_fpp/from_T4_fpp_companion/t4_02_checks.json` → `rv` |
| no non-single-star solution | NSS tables empty | `../01_inputs/catalog_queries/gaiadr3_nss_*_target.csv` |

**Documented, not printed.** The ARIADNE SED clause was removed from the Note in the pass-3 revision (2026-09-25 evening); the facts `STAR.R_ARIADNE` and `STAR.Teff` stay in the package as supporting values.

| quantity | value | where in the package |
|---|---|---|
| ARIADNE SED: R = 0.6205 +0.0200/−0.0181 Rsun, Teff = 4017 +58/−64 K (were printed as 0.62 ± 0.02 and 4020 ± 60 until the pass-3 revision) | mixture of the Bayesian-model-average posterior: R 0.6205 +0.0200/−0.0181, Teff 4017 +58/−64 | `from_T6_star/out/t6_01_ariadne_bma.json` → `par.rad.mixture_pct`, `par.teff.mixture_pct`. Script `t6_01_ariadne_bma.py` reads the ARIADNE `BMA.pkl`, which is not included (35.6 MB pickle; listed in `../12_facts/EXCLUDED_FILES.csv`). The fit set-up and the `.dat` summaries are in `../05_companion_limit_038/from_companion_check/sed_out/`. |

**Rotation and activity (not printed in the Note).** No rotation period is established.
- A 7.97 d, ~550 ppm modulation appears only in the Sector 32 PDCSAP light curve; it is absent in Sector 6: `from_T6_star/verify/v02_s32_rotation_ls.json`, `v03_other_sectors_ls.json`.
- The star is chromospherically quiet.
- It has a marginal eRASS1 X-ray detection with log L_X/L_bol ≈ −4.8: `from_T6_star/out/t6_03_xray_rossby.json`.

Facts: `STAR.rotation`, `STAR.activity_status`, `STAR.xray`, `STAR.chromosphere`.

**Source folders.**
- `from_T6_star/`: reconciliation reruns `t6_01`–`t6_05` with logs, the fresh Gaia comparison sample, and verifier checks.
- `from_extra_science_rv_check/`: Gaia RV-scatter context.
- The opus55 stellar step (`s01_relations.py`, `s02_sed.py`, catalogue data) is in `../05_companion_limit_038/from_U_06_stellar/`, because it also contains the overluminosity derivation.

**Not included.** The Wright et al. (2011) paper text used by the verifier (third-party).
