# 05_companion_limit_038: limits on an unresolved bound companion

Start with **`DERIVATION.md`**. It walks through, step by step and with equations, where these numbers come from:

- XP overluminosity → ≈0.38 Msun;
- the 0.35-0.42 Msun spread for different scatter definitions;
- the BP/RP flux excess (C*) argument for 0.2-2″;
- the circular-orbit host mass ≈0.52 Msun;
- the host range 0.15-0.38 Msun with e ≳ 0.2-0.6.

| file / folder | content |
|---|---|
| `DERIVATION.md` | Step-by-step derivation, with the file and key for every number. |
| `recompute_mass_limit.py` | Standalone (numpy, pandas). Reproduces the mass limits from the Way et al. catalogue row and the Pecaut & Mamajek table. Tested: 0.381 Msun (2σ, σ = 0.068), 0.349 (σ = 0.052), 0.420 (σ = 0.10), 0.423 (3σ), 0.511 (Way flag threshold). |
| `FACTS_IN_THIS_SECTION.csv` | `COMP.*` and `STAR.overluminosity` facts, with evidence mapped to package paths. |
| `from_companion_check/` | First overluminosity check (`overluminosity.py`, `out/overluminosity.json`, `plots/`) and ARIADNE SED check (`sed_fit.py`, `sed_residuals.py`, `sed_out/*.dat`, `out/sed_residuals.json`). `data/way2026_*_row.csv` hold the target's rows from the Way et al. (2026) Zenodo catalogue (doi:10.5281/zenodo.18500082, CC-BY-4.0). |
| `from_U_06_stellar/` | opus55 audit stellar step. `s01_relations.py` (Mann relations; also the source of R, M, ρ, M_Ks in `../09_star`), `s02_sed.py`, `s03_overluminosity.py` (independent overluminosity limit; `out/s03_overluminosity.json`, `out/pecaut_mamajek_parsed.csv`), `s04_misc.py`, catalogue data (`data/`: Gaia DR3 source/AP/XP, 100 pc comparison samples, TIC, 2MASS, AllWISE, CatWISE, APASS, Tycho-2, SIMBAD ids, Bailer-Jones distances, `PROVENANCE.jsonl`). |
| `from_U_09_companion/` | Companion space. `c04_profile.py` fits the transit shape against host mass with dilution (`out/profile_A.json`, `out/companion_space.json`). `c05_rv.py` gives Gaia RV limits (`out/rv_limits.json`: 95% K < 1.16 km/s). `c07_trice_prior.py` checks the TRICERATOPS prior. `verify/` has the independent fits (`v3_fit_A.json`) and the RUWE percentile. |
| `from_U_referee_resolve_R4/` | `r4_01_way_oos.py` (scatter definitions, out-of-sample RUWE), `r4_02_speckle_calib.py` (C* calibration with SOAR speckle binaries × Gaia DR3; `data/soar_gaia_*.csv`), `r4_03_molusc_regimes.py`, `r4_04_host_ecc.py` (circular-host threshold, e_min). `out/r4_summary.json` collects the key numbers. |

Also related:
- `../04_fpp/from_T4_fpp_companion/t4_02_checks.json` → `overlum_fluxratio_2sig` (flux-ratio arithmetic) and `rv` (Gaia RV mass function);
- `../09_star/` (RUWE and RV context).

**Not included.**
- ARIADNE posterior pickles (`sed_out/*.pkl`, 5-36 MB each) and a partial NetCDF file.
- The Way et al. paper text and the full Zenodo catalogue.
- The Mann et al. source tarballs and repositories (`06_stellar/refs/`).

They are listed in `../12_facts/EXCLUDED_FILES.csv`; the public sources are cited in `DERIVATION.md`.

## Gaia radial-velocity limit K < 1.2 km/s (printed in Sec. 4) — provenance

Source: `from_U_09_companion/c05_rv.py` → `out/rv_limits.json` (fact `COMP.Gaia_RV`).

- **Input (Gaia DR3 gaia_source, this star):** radial_velocity_error 0.221 km/s, rv_nb_transits 33,
  rv_visibility_periods_used 15, rv_time_duration 922.9 d, rv_amplitude_robust 3.22 km/s,
  rv_chisq_pvalue 4.06e-05, rv_renormalised_gof 3.00. Gaia DR3 publishes no epoch radial velocities for this star,
  so the epoch scatter is inferred from these summary fields: σ_epoch(observed) = 0.88 km/s from the p-value (χ² = 73.6, ν = 32),
  typical epoch scatter of single stars in the same (G_RVS, template T_eff) bin = 0.75 km/s from rv_renormalised_gof, excess = 0.46 km/s.
- **Observing schedule: simulated, not recovered.** The Gaia epochs are not published; `c05_rv.py` simulates 4 000 Gaia-like schedules with
  33 transits spread over 15 visibility periods within 923 d, transit pairs 106.5 min apart (the two fields of view) and 6-h spin steps.
- **Noise model:** white noise per epoch equal to the typical single-star epoch scatter (0.75 km/s).
- **Orbit model:** circular, random phase, sin i = 1 for the quoted mass function (edge-on is the conservative choice for an eclipsing companion).
- **Definition of the limit:** K_95(P) is the semi-amplitude for which 95 % of the simulations give an epoch standard deviation larger than the observed one.
  At P = 3.1828 d: **K_95 = 1.16 km/s** (printed "K < 1.2"). Using rv_amplitude_robust instead of the standard deviation gives 1.05 km/s.
- **The two numbers 1.16 and 1.24 km/s:** 1.16 km/s is the 95 % limit from the simulation above; 1.24 km/s = √2 × σ_epoch(observed) is the value if *all* observed scatter
  were orbital (a model-free upper bound, `K_conservative_all_scatter_kms`); 0.64 km/s would be the amplitude if only the excess scatter were orbital.
- **Comparison:** an eclipsing stellar companion of the target at 3.18 d needs K = 14.9 km/s for 0.08 Msun (sin i = 1); every version of the limit is ≥ 12× smaller,
  which is why the Note says "exclude". The limit weakens for P ≳ 3 yr (beyond the Gaia baseline) and says nothing about wide bound companions.
- The star fails only the gof term of the Gaia DR3 RV-variability rule (p < 0.01 and gof > 4 and N ≥ 10): p = 4.1e-5, gof = 3.0, N = 33.
