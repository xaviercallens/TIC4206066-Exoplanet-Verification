# 06_pixels_localisation: where the signals come from on the sky

## Summary tables (start here)

| file | content |
|---|---|
| `OFFSETS.csv` | Difference-image PRF offsets from the target (arcsec, RA/Dec) for each data set, and the combined positions. **3.18 d:** S32 SPOC 2-min, POS_CORR-decorrelated, 8 transits: (+5.21 ± 3.08, +0.39 ± 2.78); S98 TESScut, 9 transits: (+5.65 ± 5.27, −0.32 ± 4.52); combined (+5.32 ± 2.66, +0.20 ± 2.37), printed as "(+5 ± 3, 0 ± 3)″". **11.13 d:** S98 only, 3 transits: (+1.19 ± 2.29, −1.02 ± 3.12), printed as "(+1 ± 3, −1 ± 3)″". The rows for Sector 6 and for S32 11.13 d say why those sets are not used. Source: `from_U_referee_resolve_R3/out/r4_panel_d_spec.json` (the same values are plotted in Fig. 1b). |
| `HOST_PROBABILITY_e3.csv` | Per-star relative host probability for A and B. It combines the PRF Δχ² of forcing the source onto each star, the injection-calibrated error scale, and depth feasibility. Script and JSON: `host_probability_e3/`. Caveat: a companion unresolved from the target (≲1-2″) cannot be told apart from the target here. |
| `FACTS_IN_THIS_SECTION.csv` | Localisation and pixel re-check facts (`A.offset*`, `A.capable_*`, `A.S6_diffimg`, `PIX.*`, `PRF.method`, `B.offset`), with evidence mapped to package paths. |

## Method (fact `PRF.method`)

- **PRF model:** the TESS pixel response function, fitted with `TESS_PRF` 0.1.3 (Bell & Higgins 2022), using the MAST `start_s0004` PRF files.
- **Cutouts:** TESScut/Astrocut (Brasseur et al. 2019).
- **Pointing correction:** POS_CORR decorrelation is applied to S32 SPOC 2-min data only. The TESScut POS_CORR columns are zero.
- **Errors:** per axis, the largest of three:
  - the SD over 300 model-free injections into the real pixels (`from_U_03_pixels/inj/*.json`, `inject.py`);
  - the jackknife SE;
  - the robust bootstrap SD (`from_U_referee_resolve_R3/r1_boot.py`, `out/boot_*.npy`).
- **Combination:** `from_U_referee_resolve_R3/r2_combine.py`, log `logs/r2_combine.log`. Across Sector 98 reductions the 3.18 d offset is at most 2.6σ from the target (`from_T3_pixels_localisation/t3_combine.json`, `t3_s98_driver.json`).
- **Sector 6:** all 6 transits lie 0.41-3.54 h from momentum dumps, so its difference images are not used (`from_U_03_pixels/s6_dumps.json`).

## Neighbour census

**Which stars can produce the depth.** `from_U_03_pixels/neighbours.py` lists every Gaia DR3 star within the apertures (`neighbour_census*.csv`, `neighbours_vs_fitted_position.csv`, `neighbour_census_summary.json`). For each star it gives the eclipse depth it would need in the real S6/S32/S98 apertures.

- 5 stars could produce the 3.18 d depth. The nearest is TIC 4206077 at 42.8-43.5″. The others are at 50.5″, 53.9″, 56.3″ and 84.1″.
- With rounded ±3″ errors, TIC 4206077 lies 14.3σ from the fitted position; with the exact covariance, 17.6-18.0σ. The Note prints "≳43″ (>10σ)".

**Δχ² for each star:** `from_U_03_pixels/neighbours_dchi2.py` → `neighbours_dchi2.json`.

**Nearest Gaia source of any brightness:** 10.35-10.86″ at the TESS epochs; it can produce at most 406-418 ppm. Source: `from_U_referee_resolve_R3/logs/r4_panel_d_spec.log`.

## Verified pixel re-check (aperture test)

`from_tic4206066_pixel_recheck/` is the author's exploratory re-check, byte-identically reproducible. `from_T8_pixel_recheck/` is its verification.

**Verified outcome** (facts `PIX.*`):
- **Signal B, Sector 32 "off-target" pattern: a methodological false positive.** The first pass found a 4.5-5.0σ pattern, caused by global detrending during a pointing drift. Per-event fits give flat depths across apertures: r1 − r3 = +8 ± 120 ppm, and r2 − r3 = −81 ± 59 ppm. Sources: `from_T8_pixel_recheck/out/t8_b_and_paired.json`, `t8_b_drift.json`, `from_tic4206066_pixel_recheck/out/after_decorrelated.json`.
- **Signal A, aperture (central-pixel) effect: not confirmed.** The only >3σ value (S98, 226.6 ppm) comes entirely from one event, epoch 8, which sits next to a gap with a momentum dump. S32+S98 give −15.9 ± 56.7 ppm. Sources: `from_T8_pixel_recheck/out/t8_s6_and_clipping.json`, `t8_t9_decomposition.json`, `t8_a_events_and_nulls.json`.
- **A source at TIC 4206077 would make the core shallower, not deeper** (`t8_prf_offtarget_signature.json`). The aperture test therefore has no bearing on the ≳43″ exclusion.
- **The Note makes no aperture claim.**
- **Superseded text in the re-check.** `from_tic4206066_pixel_recheck/SUMMARY.md` records the first-pass "~4σ" B pattern and B centroid numbers that predate its v3 fixes. The JSON outputs and the T8 verification replace them; see `../11_superseded/SUPERSEDED.md`.

## Key figures

- `from_U_03_pixels/fig_diffimages.png`: difference images.
- `from_U_03_pixels/fig_localisation.png`: fitted positions.
- `from_tic4206066_pixel_recheck/plots/*.png`: aperture tests.
- Fig. 1b of the Note: `../08_figure/`.

## Source folders

| folder | content |
|---|---|
| `from_U_03_pixels/` | PRF localisation: `pixloc.py`, `run_real.py`, `poscorr_decorr.py`, `inject.py`, `summarize.py`, `g1_prf_check.py`, neighbours, Sector 6 dumps; results `real_results.json`, `summary.json`; `verify/`. |
| `from_U_referee_resolve_R3/` | Bootstrap errors, combination, Gaia contrast, Fig. 1b specification (`out/r4_panel_d_spec.json`). |
| `from_T3_pixels_localisation/` | Reconciliation reruns: combination variants, S98 drivers for A and B, cube comparison, the 3995 gap event. |
| `from_T8_pixel_recheck/` | Verification of the pixel re-check (see above). |
| `from_tic4206066_pixel_recheck/` | The re-check package itself: scripts `s01`–`s06`, `run_all.sh`, outputs, plots, reviews, its own `SHA256SUMS` and `env/pip_freeze.txt`. |
| `host_probability_e3/` | `e3_host_probability.py` / `.json`. Extra analysis, not printed in the Note. |

**Not included.** The TESScut cubes and TPFs (see `../01_inputs/INPUT_PRODUCTS.csv`). Two folders in T8, `pkg_copy/` and `orig_out_snapshot/`, are duplicates of `from_tic4206066_pixel_recheck/` and are not copied.
