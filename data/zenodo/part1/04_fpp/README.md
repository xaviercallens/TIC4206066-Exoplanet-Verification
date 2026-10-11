# 04_fpp: TRICERATOPS false-positive probabilities (3.18 d signal)

TRICERATOPS 1.1.0 (Giacalone et al. 2021) was run in `trice_env` (Python 3.11.15; see `../10_environment/pip_freeze_trice_env.txt`). Every run uses N = 10^6 draws unless the variant name says otherwise.

## Summary tables (start here)

| file | content |
|---|---|
| `FPP_RUN_CENSUS.csv` | Every stored run set, for signals A and B: source, variant, number of runs, FPP mean, SD, min, max, NFPP max, and `package_path` to the per-run file. Copy of `from_T4_fpp_companion/t4_01_fpp_run_census.csv`, which recounts the runs from the run files, plus the package path column. |
| `FPP_DECOMPOSITION_038.csv` | FPP split by companion mass (signal A). Two kinds of rows. (1) U 08 re-weighting of the same draws, per variant: FPP, FPP with companions ≤ 0.38 Msun only, share from companions > 0.38 Msun, and STP evidence fraction > 0.38 Msun. (2) Direct runs with a modified companion prior: truncated at 0.38 or 0.45 Msun, and the Gaia/MOLUSC-constrained runs. |
| `recompute_share_038.py` | Standalone and stdlib-only. Recomputes the share and the post-cut FPP from `from_U_08_fpp/out_fpp_summary.csv`. Output (tested): share by FPP drop 0.837-0.951; STP weight > 0.38 Msun 0.854-0.956; FPP after the cut 0.0020-0.0067. |
| `FACTS_IN_THIS_SECTION.csv` | The FPP facts from `FINAL_FACTS.csv`, with evidence mapped to package paths. |

## How the printed numbers follow

- **"FPP = 0.03-0.04 (4-20 runs of 10^6 draws)"** (`FPP.A_default`, `FPP.A_runs`).
  - These are default priors with S32+S98 folds in 5-10 min bins. The variant means are 0.030-0.044 over 15 configurations, and single runs range 0.026-0.048.
  - The configurations come from the audit (20 runs each), U 08 (5 or 20), the U 08 verifier (4 or 8), resolve R5 (5) and pubprep (10). See `from_T4_fpp_companion/t4_02_checks.json` → `A_default_S32S98` and `A_default_S32S98_summary`.
  - The census also lists a 2-seed reproduction (`default_A_adoptedfold`, mean 0.028). It is not one of the 15 configurations.
- **"nearly all from bound companions"** (`FPP.A_scenario_composition`). The STP scenario (planet on an unresolved bound companion) is ≥ 99.9% of FPP_A for S32+S98 folds. See `from_T4_fpp_companion/t4_06_eb_totals.json`.
- **"about 90% from those above 0.38 Msun"** (`FPP.A_share_above_0.38`; 84-96% across variants).
  - `from_U_08_fpp/f02_trice_run.py` wraps the TRICERATOPS likelihood functions with recording hooks: `q_hook`, `e_hook`, `lme_hook` and `summarize()`. For each run these store the companion mass and eccentricity of every STP/PTP draw with its evidence weight, the weighted mass fractions (`wfrac_mass_gt_0.38`), and the evidence restricted to companions ≤ 0.38 Msun (`lnZ_mass_le_0.38`).
  - From these, FPP is recomputed with the restricted evidence (`FPP_comp_mass_le_0.38`).
  - `f04_summarize.py` averages over runs into `out_fpp_summary.csv`.
  - `from_T4_fpp_companion/t4_02_checks.py` forms the share = 1 − FPP(≤ 0.38)/FPP (keys `share_gt038_A_default`, `share_gt038_A_range`). `recompute_share_038.py` repeats this step standalone.
- **"FPP ≈ 0.002-0.007 after a post hoc prior cut"** (`FPP.A_constrained`).
  - Re-weighting per variant gives 0.0020-0.0067; the 0.0067 is the single-sector S32 fold, and S32+S98 folds give 0.002-0.005.
  - Direct runs agree: companion prior truncated at 0.38 Msun, 0.0043 ± 0.0006 (8 runs, 5-min fold; `from_U_08_fpp/verify/runs_v/trunc0.38_A_v_A_S32+S98_5min.json`); companions > 0.38 Msun removed at all separations, 0.0024 ± 0.0006 (5 runs; `from_U_referee_resolve_R5/runs/ol0.38_all_A_adoptedfold.json`).
  - Fig. 1c plots the 5-min truncated runs.
- **"NFPP < 10^-6"** (`FPP.A_NFPP`). The maximum over all stored signal-A runs is 4.6 × 10^-7 (`t4_02_checks.json` → `NFPP_max_all_default_A`).
- **Validation threshold** "FPP < 0.015, NFPP < 10^-3" is a literature value (Giacalone et al. 2021). The paper is not redistributed. No validation is claimed.
- **Companions at 0.1-2.2″** (`FP.remaining_scenario_A`). The TRICERATOPS bound-companion prior counts only P > 2500 d, i.e. 3.6-78 au = 0.1-2.2″ at 35.6 pc.

## Documented facts only (not printed in the Note)

- **MOLUSC mode** (`FPP.A_molusc_mode`). Running TRICERATOPS with a MOLUSC companion file gives ~0.005 at N = 10^6, but the value scales with N (for example, `A_moluscFresh5_N5e5` and `A_moluscFresh5_N2e6` in the census). It is **not** a physical constraint and is not quoted. Runs: `from_U_08_fpp/runs/A_moluscFresh5*`, `../13_supporting_pipelines/PUBLICATION_AUDIT/out/fpp/runs/A_molusc{1,5}_S32+S98`. The MOLUSC draw tables (`molusc_*_kept.csv.gz`, 36 MB each) are listed in `../12_facts/EXCLUDED_FILES.csv`.
- **Gaia RUWE + contrast constrained, properly normalised** (`FPP.A_gaia_constrained`): 0.0030 ± 0.0008 (20 runs, range 0.0019-0.0046). Source: `from_U_referee_resolve_R5/runs/gaia_only_A_adoptedfold.json`, scripts `r01_molusc_all.sh`, `r05_direct_runs.py`, `r07_ruwe_robustness.py`.
- **Blended-only cut** (< 0.2″) alone: 0.024 ± 0.004 (`ol0.38_bl0.2_A_adoptedfold`). This is not the constrained FPP.
- **Signal B**: FPP 0.13-0.21 with default priors. Not printed.

## Source folders

| folder | content |
|---|---|
| `from_U_08_fpp/` | `f01_build_folds.py` builds the folds (`folds/*.csv`: 2, 6, 10 min bins, S32/S98/S32+S98, ±T0 shifts). `f02_trice_run.py` runs TRICERATOPS with the recording hooks, apertures and stellar parameters (arguments: M 0.606, R 0.623, Teff, parallax 28.075, `--ap spoc` or 5×5). `f03_priors.py`, `f04_summarize.py`, `f05_plot.py`; `run_queue*.sh` hold the exact command lines. `runs/*/runs.json` has per-run FPP, NFPP, scenario probabilities and STP draw summaries. The verifier (`verify/`) has its own folds, direct runs and summaries. `4206066_gaia_background.csv` is the Gaia background star list; `runs/stars_32_98_TIC.csv` is the TIC star list used by TRICERATOPS. |
| `from_U_referee_resolve_R5/` | Constrained-prior runs: Gaia-only, overluminosity ≤ 0.38 at all separations or blended-only, Teff 4012. Summary `out_r06_summary.json`. |
| `from_T4_fpp_companion/` | Reconciliation census (`t4_01`), checks (`t4_02`: RV mass function, shares, ranges, bin width, NFPP, overluminosity arithmetic), SDSS separation (`t4_03`), Gaia 60″ census (`t4_04`), SED colour sign (`t4_05`), EB totals (`t4_06`), and verifier. |
| `../13_supporting_pipelines/PUBLICATION_AUDIT/` | Independent audit FPP runs: `out/fpp/runs/*/probs_runNN.csv` (20 runs per variant), `out/fpp/inputs/` (binned folds, `apertures.json`, `inputs_meta.json`), `scripts/s13_fpp_inputs.py`, `s13_fpp_run.py`, `run_fpp_all.sh`, logs `logs/fpp_*.log`. |
| `../13_supporting_pipelines/pubprep/analysis/` | `a4_fpp_ownfold.py` and `out/fpp/fpp_ownfold_{A,B}.json` (10 runs each). |

**Not included.** The TRICERATOPS target pickles (`target_32_98.pkl`, `target_v.pkl`) and the 70 MB queue log. They are listed with sha256 in `../12_facts/EXCLUDED_FILES.csv`; the pickles can be rebuilt with `f02_trice_run.py`. The run scripts import TRICERATOPS and read inputs from the author's tree, so they do not run standalone from this package.
