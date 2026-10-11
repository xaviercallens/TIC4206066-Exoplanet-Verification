# 08_figure: Figure 1 of the Note

| file | content |
|---|---|
| `figure.pdf`, `figure.png` | The figure exactly as it goes with the submission-candidate manuscript. The PDF is byte-identical to `submission_candidate/figure.pdf` (sha256 in `../SHA256SUMS`). |
| `figure_values_final.json` | Every plotted value. **Panel (a):** 3.18 d fold ephemeris, per-sector bin centres, means and errors, template depths S6 521 / S32 476 / S98 496 ppm. **Panel (b):** S32 and S98 offsets, errors, the 1σ rings, and the capable Gaia DR3 stars. **Panel (c):** TRICERATOPS FPP mean, SD and range for the three 6-min-bin folds, TIC parameters, and companions capped at 0.38 Msun (5-min fold). It also holds `manuscript_consistency` (plotted vs printed numbers) and `sources`: every input file with its sha256. All 13 sources were re-hashed when this package was built and still match. |
| `make_figure_v2.py` | The script that produced `figure.pdf`/`.png` and `figure_values_final.json` (2026-09-25 17:15 UTC; Python 3.12.13, numpy 2.5.3, matplotlib 3.11.2). |
| `make_figure_v2.patch.diff` | Difference from the U proposed-revision version. The per-sector template depths now use FINAL_FACTS `A.depth_S6/S32/S98`; the S98 segment-edge epoch −1 is omitted from the fold; the capable-star label is "≳43″ (>10σ)"; output paths changed. |
| `original_U_proposed_revision_figure_src/` | The U version of the script and its value file, for comparison. |

**Reproducing.** `make_figure_v2.py` does not run standalone from this package:
- it reads its inputs from the author's tree (absolute paths at the top of the script);
- it imports `rnaas_tic4206066_pubprep/analysis/a1_timing_depths.py` (in this package: `../13_supporting_pipelines/pubprep/analysis/`) for light-curve loading, detrending and the template;
- it needs the three TESS light curves listed in `../01_inputs`.

To run it, set `REV`/`SC_OUT` and the `a1` paths to local copies. The expected output is `figure_values_final.json`.

**Sources of the panels:**
- (a) `../02_transits_A` (ephemeris, depths);
- (b) `../06_pixels_localisation/OFFSETS.csv` and `from_U_referee_resolve_R3/out/r4_panel_d_spec.json`;
- (c) `../04_fpp/FPP_RUN_CENSUS.csv`.
