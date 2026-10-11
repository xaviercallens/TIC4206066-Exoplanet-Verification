# Derivation of the companion limits used in the Note (Sec. 4)

This file explains, step by step, where four statements in the Note come from:

- Gaia-blended companions above about 0.38 Msun are disfavored (2 sigma);
- at 0.2-2 arcsec such companions would raise the BP/RP flux excess, which is normal;
- a circular orbit needs a companion host above about 0.52 Msun;
- hosts of 0.15-0.38 Msun need e ≳ 0.2-0.6.

All numbers below are copied from files in this package. `recompute_mass_limit.py` (in this folder) reproduces steps 1-4 using only numpy and pandas.

## A. XP overluminosity → mass limit (~0.38 Msun)

**Step 1. Measured overluminosity.**
Way et al. (2026, AJ, arXiv:2602.20096) trained a random-forest model on single K/M dwarfs. It predicts the absolute magnitude M_G of a *single* star from its Gaia XP (RP) spectrum. The catalogue row for the target (Gaia DR3 3220388198192519424; `from_companion_check/data/way2026_overluminous_row.csv`) gives:

- measured M_G = 7.755 (`absG_dr`);
- predicted single-star M_G = 7.758 (`pabsG_dr`).

So

    O = M_G,predicted − M_G,measured = +0.003 mag

The catalogue flag is `ol = False`. The systems catalogue lists 1 component. Sources: `from_companion_check/out/overluminosity.json` (`O_obs_mag`) and `from_U_06_stellar/out/way2026_target_row.csv`.

**Step 2. Scatter of O among comparable stars.**
Way et al. give no per-star error. The scatter is the author's own estimate, taken from catalogue stars with the same predicted M_G (±0.25 mag) and the same predicted colour (±0.05):

- σ = 1.4826 × MAD = **0.068 mag** (n = 5205; median O = +0.016). Source: `from_companion_check/out/overluminosity.json` → `local_bin`.
- Alternative definitions (`from_U_referee_resolve_R4/out/r4_01_way_oos.json` → `local_O_by_ruwe`):
  - 0.052 mag: clipped core;
  - 0.053 mag: RUWE 1.0-1.2, the regression-eligible stars;
  - 0.054 mag: RUWE 1.2-1.3, never in the training set. The target's RUWE is 1.214, so it lies outside Way's regression sample (`ruwe < 1.2`).
  - 0.10 mag: a conservative value.
- The Note prints this as "0.05-0.07 mag".

**Step 3. Upper limit on the companion's share of the G-band flux.**
A companion adds light, so the measured star looks brighter than the single-star prediction. At k sigma:

    O_lim = O + k·σ
    f_G   = 10^(0.4·O_lim) − 1          (largest allowed companion/primary flux ratio in G)

For k = 2 and σ = 0.068: O_lim = 0.139 mag and f_G = 0.137.

This is the same flux-sum logic Way et al. use for the equal-mass offset, 2.5 log10 2. A companion also reddens the XP spectrum, which makes the predicted M_G fainter and O larger. The simple conversion is therefore on the conservative side (`from_companion_check/SUMMARY.md` §1).

**Step 4. Flux ratio → companion absolute magnitude → mass.**

    M_G,2 = M_G,predicted − 2.5·log10(f_G) = 7.758 − 2.5·log10(0.137) = 9.92

M_G,2 is converted to a mass by linear interpolation in the Pecaut & Mamajek (2013) dwarf table, version 2022.04.16 (`../01_inputs/reference_tables/mamajek_table.txt`, columns `M_G` and `Msun`). This gives **M_2 < 0.38 Msun** (q < 0.63).

Results (`recompute_mass_limit.py`, `from_companion_check/out/overluminosity.json` → `limits`, `from_U_referee_resolve_R4/out/r4_summary.json` → `a_attribution_and_mass_limit`):

| scatter σ (mag) | 2σ limit (Msun) | 3σ limit (Msun) |
|---|---|---|
| 0.052 (clipped core) | 0.349 | 0.393 |
| 0.053-0.054 (RUWE-matched, in and out of the training set) | 0.353-0.355 | 0.395-0.396 |
| 0.068 (1.4826 MAD, adopted) | **0.381** | 0.423 |
| 0.10 (conservative) | 0.420 | 0.469 |
| Way et al. flag threshold O = 0.376 mag | 0.511 | – |

So the 2σ limit is **0.35-0.42 Msun, depending on how the scatter is defined**. The Note uses ≈0.38 Msun, the conservative end of the 0.05-0.07 mag range it prints. In `FINAL_FACTS.csv` the fact is `COMP.XP_overluminosity`; its uncertainty column also lists 0.373-0.383 for σ = 0.068 when Mann instead of Pecaut & Mamajek masses are used.

**Step 5. Which separations this covers.**
The test only sees companions whose light is blended into the Gaia measurement. That means ≲0.2″, the Gaia resolution limit that Way et al. discuss (≈7 au at 35.6 pc). Way et al.'s "≲20 au" refers to 100 pc and is not used here. The limit says nothing about fainter companions at any separation.

## B. BP/RP flux excess (C*) for 0.2-2″

- Target C* = −0.0025. Using Riello et al. (2021), eq. 18, σ_C*(G = 10.51) = 0.0065, so C* is −0.39σ: normal. Source: `from_U_referee_resolve_R4/out/r4_02_speckle_calib.json` → `target`.
- **Calibration.** SOAR speckle binaries (Ziegler et al. 2020, 2021) were matched to Gaia DR3 (`from_U_referee_resolve_R4/data/soar_gaia_*.csv`, script `r4_02_speckle_calib.py`). Among hosts with one companion at 0.2-2″:
  - with ΔI ≤ 2: 0 of 38 have a C* as normal as the target (95% upper limit 7.6%);
  - with ΔI ≤ 1.76 (≥0.38 Msun): 0 of 32 (upper limit 8.9%).
  - For comparison, 38% of target-like records without companions look like the target (n = 647).
- **Conclusion.** A ≥0.38 Msun companion at 0.2-2″ would very likely raise C*, and C* is normal. The calibration only covers ΔI ≤ 2, so it says nothing about 0.15-0.38 Msun companions. Fact: `COMP.BPRP_excess`.
- Beyond about 1.2″, companions of this brightness would be resolved Gaia sources, and there is no other Gaia source within 10″ (`A.nearest_gaia_neighbour`).

## C. Circular-orbit host mass (~0.52 Msun)

The idea: if the 3.18 d transit were on a companion star, the transit duration fixes the host's density, which gives a minimum host mass.

1. Fit the transit with the stellar density free. The data allow at most ρ_max ≈ 3.5-3.6 ρsun: 3.53 at Δχ² = 4 (`from_U_09_companion/REPORT.md`), and 97.7% limits of 3.4-3.7 (fact `A.rho_free`).
2. For a circular orbit, the host must have ρ_host ≤ ρ_max. The Mann et al. (2015) R(M_K) relation and the Mann et al. (2019) M(M_K) relation give ρ(M_K) = M/R³. Solving ρ(M_K) = ρ_max gives:

   | ρ_max | host mass |
   |---|---|
   | 3.5 | 0.529 Msun |
   | 3.6 | 0.523 Msun |
   | 3.64 | 0.520 Msun |

   Source: `from_U_referee_resolve_R4/r4_04_host_ecc.py`, output `out/r4_04_host_ecc.json` → `undiluted_threshold_Mann`, summarised in `out/r4_summary.json` → `d_threshold`.
3. Fits that include dilution by the target's light give the same answer:
   - 0.523 Msun at Δχ² = 4 and 0.506 at Δχ² = 9 (`from_U_09_companion/out/profile_A.json`, `out/companion_space.json`, script `c04_profile.py`);
   - 0.517 Msun from the verifier (`from_U_09_companion/verify/v3_fit_A.json`).

   Result: **≈0.52 Msun** (0.517-0.529). Fact: `COMP.host_mass_circular`. This holds for the 3.18 d signal only; the 11.13 d shape allows any host ≥0.14 Msun.

## D. Lower host-mass bound and eccentricity (0.15-0.38 Msun, e ≳ 0.2-0.6)

**1. Lower bound on the host mass.** A very faint host needs a large, heavily diluted planet, and that makes the transit V-shaped. The diluted-host fits disfavour hosts below:

- 0.168 Msun at Δχ² > 4 (`from_U_09_companion/out/profile_A.json`);
- 0.13-0.14 Msun in the verifier.

Result: **≈0.15 Msun** (0.13-0.17, depending on method).

**2. Eccentricity needed by hosts denser than ρ_max** (photo-eccentric effect; Kipping 2010, Dawson & Johnson 2012). An eccentric orbit changes the transit speed relative to a circular one by

    g = (1 + e sin ω) / sqrt(1 − e²)

A circular fit then returns ρ_circ = ρ_host · g³. A host with ρ_host > ρ_max therefore needs

    g ≤ (ρ_max / ρ_host)^(1/3)

The most favourable case is a transit at apoapsis (ω = 270°), where g = sqrt((1−e)/(1+e)). This gives

    e_min = (1 − g²)/(1 + g²)

**3. Result** for ρ_max = 3.5 (`from_U_referee_resolve_R4/out/r4_04_host_ecc.json` → `eccentricity_requirement`, summary `r4_summary.json` → `e_ecc`):

| host mass (Msun) | e_min |
|---|---|
| 0.14 | 0.61 |
| 0.17 | 0.53 |
| 0.20 | 0.47 |
| 0.30 | 0.31 |
| 0.38 | 0.20 |

Result: **e ≳ 0.2-0.6**. Under the TRICERATOPS eccentricity prior, Beta(0.867, 3.03) with uniform ω, only 0.5-11% of draws satisfy the requirement (`r4_summary.json` → `e_ecc.*.prior_frac`). Fact: `COMP.host_lower_bound_and_ecc`.

## E. How this enters the FPP (see `../04_fpp`)

TRICERATOPS 1.1.0 models unresolved bound companions at 0.1-2.2″ (P > 2500 d), with masses up to the target's mass. About 90% (84-96%) of the default FPP comes from companions above 0.38 Msun (`../04_fpp/FPP_DECOMPOSITION_038.csv`, `../04_fpp/recompute_share_038.py`). Removing those companions from the prior at all separations is a post hoc cut and gives FPP ≈ 0.002-0.007. That cut is justified by:

- A (≲0.2″);
- B (0.2-2″);
- resolved Gaia sources (≳1.2″).

## Scripts and outputs in this folder

- `from_companion_check/overluminosity.py` computes the local scatter and the 2σ/3σ limits (`out/overluminosity.json`, `plots/overluminosity.png`). It needs the Way et al. Zenodo catalogue files (doi:10.5281/zenodo.18500082, CC-BY-4.0: `overluminous_catalog.csv`, `catalog_of_systems_100pc_lite.csv`), which are not included here. Download them and pass their folder as the first argument. The script also reads the Pecaut & Mamajek table from the author's tree; a copy is in `../01_inputs/reference_tables/`.
- `from_U_06_stellar/s03_overluminosity.py` is an independent re-derivation (`out/s03_overluminosity.json`, `out/pecaut_mamajek_parsed.csv`). It also includes the Gaia 100 pc comparison samples (`data/gaia_100pc_*.ecsv`).
- `from_U_referee_resolve_R4/`:
  - `r4_01_way_oos.py`: scatter definitions and the out-of-sample check;
  - `r4_02_speckle_calib.py`: C* calibration;
  - `r4_03_molusc_regimes.py`;
  - `r4_04_host_ecc.py`: host-mass threshold and e_min.
- `from_U_09_companion/`: transit-shape profile fits vs host mass (`c04_profile.py`), Gaia RV limits (`c05_rv.py`, `out/rv_limits.json`) and the TRICERATOPS prior check (`c07_trice_prior.py`).
- `recompute_mass_limit.py` is standalone; it reproduces the table in A (tested: 0.381 / 0.349 / 0.420 / 0.423 / 0.511 Msun).
- The flux-ratio arithmetic is also re-derived in `../04_fpp/from_T4_fpp_companion/t4_02_checks.json` → `overlum_fluxratio_2sig`.
