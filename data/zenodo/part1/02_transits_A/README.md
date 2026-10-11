# 02_transits_A: the 3.18 d signal (timing, ephemeris, depths, shape)

## Summary tables (start here)

| file | content |
|---|---|
| `per_transit_A.csv` | All 23 transits in the shape fit (6 in S6, 8 in S32, 9 in S98). For each: predicted time, fitted offset, mid-time (BTJD and BJD_TDB), per-sector injection-calibrated error, white-noise error, coverage, pre/post baseline, Δχ² vs no transit (≈ S/N²), a per-event depth, and O−C. `used_in_19_timing_fit` flags the 19 transits in the ephemeris; `selection_note` explains the other 4. Times: U 02 verifier, ±60 min search (`reproduce_ephemeris/v2_real_events.csv`). |
| `ephemeris_A.json` | Adopted linear ephemeris: P = 3.18278495 ± 0.00000566 d, T0 = 2460998.735992 ± 0.002883 BJD_TDB, χ² = 9.51 for 17 dof, n = 19. Printed in the Note as P = 3.182785 ± 0.000006 d, T0 = 2460998.736 ± 0.003, χ²/dof ≈ 9.5/17. |
| `depths_A_per_sector.csv` | Template depths per sector, with β-covariance errors: S6 521 ± 53, S32 476 ± 52, S98 496 ± 48 ppm (S98 segment-edge event E-1 excluded). The adopted S32 error is the empirical 75 ppm (`A.depth_S32` in `FACTS_IN_THIS_SECTION.csv`). Combined 498 ± 29 ppm (503 ± 29 with all 23 events). χ² of a constant depth 0.37 for 2 dof (0.08-3.06 across six method sets). Source: `from_T2_A_shape/out/r01_depths.json`, key `audit_s09_beta_excl98e-1`. |
| `FACTS_IN_THIS_SECTION.csv` | Every signal-A fact from `../12_facts/FINAL_FACTS.csv`: value, uncertainty, printed form, and the evidence files mapped to paths in this package. |

The archival false-alarm test of this signal (Sector 98-only period, epoch search in Sectors 6 and 32; printed "≲1% per sector, ≲0.1% jointly") is in `../03_archival_test_A/`.

## Reproduce the ephemeris (standalone)

```
cd reproduce_ephemeris
python v4_ephem.py        # needs numpy, pandas, scipy, astropy, batman-package
```

The script reads:
- `v2_real_events.csv`: the 23 measured transit times;
- `v3_inj_S{6,32,98}.csv`: injection-recovery tests in the real light curves, which give the per-sector timing error = SD of the recovered offsets: S6 12.9, S32 16.4, S98 12.0 min (±60 min search).

It writes `v4_results.json`. Key `ephem.drop_onesided_573_and_m1_R60` is the adopted 19-transit fit. Tested 2026-09-25 with `.venv` (Python 3.12.13): the output equals the stored `v4_results.json`, since the Monte Carlo uses fixed seeds.

`v2_real.py` and `v3_inj.py` measured the times and ran the injections. They need the TESS FITS files listed in `../01_inputs` and the path `D` in `v1_common.py`.

**Selection rule for the 19 timing transits.** A transit is used if it has ≥70% coverage, S/N ≥ 2, and ≥15 min of data after egress. Excluded:
- S32 E-573: baseline truncated by the 11.13 d mask after about 4 min;
- S98 E-1: about 13 min before a data gap;
- S32 E-567: S/N < 2;
- S32 E-566: fails the S/N/coverage cut; it abuts an 11.13 d transit.

The rule is coded as the explicit event lists at lines 33-36 of `v4_ephem.py` (set `two`).

## Source folders (copied for transparency; paths inside point to the author's tree)

| folder | what it holds | used for |
|---|---|---|
| `from_U_02_timing_A/` | Timing pipeline of the opus55 final audit: auditor `t3`–`t9` plus an independent verifier (`verify/`). | P, T0, χ², TTV limits, injection errors |
| `from_U_01_ephemeris_A/` | Epoch table (`out_s03_epochs.csv`: 23 transits, coverage and quality-flag fractions), S98 quality intervals, cycle count, BLS aliases, S110 pointing. | N_all = 23; sector list |
| `from_U_07_A_shape_radii/` | Joint shape fits: emcee, exposure-integrated, Claret LD, Mann density prior or free density (`out/summary_fits.csv`, `verify/out/summary_verify_fits.csv`), odd/even, V-shape, phase 0.5. | Rp/R*, b, T14, ρ limit |
| `from_U_referee_resolve_R-S6/` | Refits and timing without Sector 6. | "Without Sector 6 … ≲1σ" |
| `from_T1_A_timing/` | Reconciliation rerun of all stored per-transit tables (`r1_ephemerides.json`: 26 pipeline/selection combinations and the drop-S6 test) and its verifier. | cross-pipeline range of P, T0, χ² |
| `from_T2_A_shape/` | Reconciliation reruns: depths from stored rows (`r01`), V-shape with red noise (`r02`), derived Rp and ρ (`r03`), and verifier checks. | depth ≈500 ppm, Rp = 1.4 ± 0.1 R⊕, Δχ² = 32 / ≥21 |

Other evidence for the shape parameters is in the independent audit pipeline, `../13_supporting_pipelines/PUBLICATION_AUDIT/`:
- `FINAL_PARAMETERS.json` → `signal_A.shape_mann_prior` (b, T14) and `signal_A.free_density` (ρ limit);
- `out/s10_eb_tests.json` (phase-0.5 depth −20.4 ± 31.7 ppm; V-shape Δχ² 31.84).

**Not included.** MCMC chains (`chain_*.npz`) and derived light-curve arrays (`lc_*.npz`, `01_ephemeris_A/lc/`) are left out; they are listed with sha256 in `../12_facts/EXCLUDED_FILES.csv` and can be recreated with the scripts from the FITS inputs.
