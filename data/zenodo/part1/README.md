# Reproducibility package: *Two Transit-like Signals in TESS Photometry of the Nearby K Dwarf TIC 4206066*

> **Deposit status.** Prepared for Zenodo; the DOI below was reserved on 2026-09-25 for this record (it resolves once the record is published).
> **DOI: `10.5281/zenodo.22967099`** (https://doi.org/10.5281/zenodo.22967099)

**Author:** Pavel Rabtsevich (author-confirmed). The affiliation printed in the Note is "Independent Researcher".

**What is NOT deposited** (details and sha256 in `12_facts/EXCLUDED_FILES.csv`):
- raw TESS light-curve/pixel FITS and TESScut cubes: listed instead with MAST URLs, retrieval times and sha256 in `01_inputs/INPUT_PRODUCTS.csv`;
- MCMC chains and derived light-curve arrays (re-creatable with the copied scripts);
- the TRICERATOPS target pickles and the ARIADNE posterior pickles;
- the frozen Sector 98 search pipeline over 126,246 stars and its LEO-Vetter runs: only the recount (`tic_recount.py`) and the LEO-Vetter re-run for this star (`leo_rerun.py`) are included, in `13_supporting_pipelines/U_10_novelty_protocol/`;
- the TRICERATOPS run scripts do not run standalone from the package (they read the author's tree and the raw inputs); the run outputs are included;
- the validation thresholds FPP < 0.015 and NFPP < 10⁻³ are literature values (Giacalone et al. 2021), not results of this package.

**Record metadata (Zenodo record 10.5281/zenodo.22967099):** title "Supporting data and code for *Two Transit-like Signals in TESS Photometry of the Nearby K Dwarf TIC 4206066*"; resource type Dataset; version 1.0.0; language English; keywords: TESS, TIC 4206066, exoplanet candidates, transit photometry, false-positive analysis. ORCID: not provided. Related identifier: the Note, to be added once published.

## Licensing

- **Data products, results, figures and documentation created by the author** (all `*.csv`, `*.json`, `*.ecsv`, `*.md`, `*.txt`, `*.png`, `*.pdf` files that are outputs of this study, including `FINAL_FACTS.csv`, the per-transit tables, the FPP run summaries, the localisation outputs and the figure data): **Creative Commons Attribution 4.0 International (CC BY 4.0)**, `LICENSE-CC-BY-4.0.txt`, to the extent of the author's rights in them.
- **Original code written for this study** (every `*.py` and `*.sh` file in the package, including the copies in `from_*` folders and `13_supporting_pipelines/`): **MIT License**, `LICENSE-MIT.txt`. Copyright (c) 2026 Pavel Rabtsevich.
- **Third-party materials** are included only as small excerpts or query results and keep their original terms and attribution; nothing in this package re-licenses them:
  - *Gaia DR3* catalogue rows (`01_inputs/catalog_queries/gaiadr3_*`, the Gaia rows in `09_star/`, `05_companion_limit_038/`): ESA/Gaia/DPAC. Gaia data are "open and free to use, provided credit is given to ESA/Gaia/DPAC"; the required acknowledgement is: "This work has made use of data from the European Space Agency (ESA) mission Gaia (https://www.cosmos.esa.int/gaia), processed by the Gaia Data Processing and Analysis Consortium (DPAC, https://www.cosmos.esa.int/web/gaia/dpac/consortium). Funding for the DPAC has been provided by national institutions, in particular the institutions participating in the Gaia Multilateral Agreement." (https://gea.esac.esa.int/archive/documentation/GDR3/Miscellaneous/sec_credit_and_citation_instructions/).
  - *TESS Input Catalog (TIC 8.2)* rows and the *TESS/MAST* product metadata in `01_inputs/` (light curves, TPFs, TESScut cutouts and PRF files are listed, not deposited): NASA/MAST public data; cite Ricker et al. 2015 (TESS), Stassun et al. 2019 (TIC), Jenkins et al. 2016 (SPOC), Caldwell et al. 2020 (TESS-SPOC), Huang et al. 2020 and Petitpas et al. 2026 (QLP), Brasseur et al. 2019 (TESScut), Bell & Higgins 2022 (TESS_PRF); MAST data DOIs are in `01_inputs/INPUT_PRODUCTS.csv`.
  - *Claret (2017) TESS limb-darkening tables* (`02_transits_A/from_U_07_A_shape_radii/claret/`, VizieR J/A+A/600/A30): copyright of the author/CDS; cite Claret 2017, A&A 600, A30.
  - *Mamajek's "Modern Mean Dwarf Stellar Color and Effective Temperature Sequence"* (`01_inputs/reference_tables/mamajek_table.txt`, version 2022.04.16, http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt) and the parsed rows in `05_companion_limit_038/`: Eric Mamajek; cite Pecaut & Mamajek 2013, ApJS 208, 9.
  - *Way et al. 2026* catalogue rows (`way2026_*` files): cite Way et al. 2026, AJ 171, 252 (their Zenodo catalogue record is referenced in `05_companion_limit_038/`).
  - *Mann et al. 2015, 2019* relations are used through their published coefficients (cite the papers); no third-party code is included.
  - *TRICERATOPS* (Giacalone et al. 2021; MIT-licensed software, not included) and *MOLUSC* (Wood et al. 2021, not included): only run outputs and the author's driver scripts are deposited; `04_fpp/from_U_08_fpp/molusc/audit08_molusc.diff` is a one-line command-line patch to MOLUSC and stays under MOLUSC's terms.
  - Python dependencies (numpy, scipy, astropy, batman, emcee, wotan, matplotlib, pandas, lightkurve, astroquery) are listed in `10_environment/` and are not redistributed.

## Running the scripts: which paths to replace

All scripts were written and run inside the author's project tree. Every hard-coded path is listed, with file, line number and variable name, in **`12_facts/PATH_CONSTANTS.csv`** (167 lines in 109 scripts). To rerun a script:

1. Recreate the tree root: choose a local folder `<ROOT>` and replace every occurrence of `/Users/paulsmith/space <3/tess-jev` with `<ROOT>` (the constants are usually named `ROOT`, `U`, `D`, `A`, `R`, `F`, `T`, `W`, `ALL`, or appear inline in `open()` / `Path()` calls; `PY` points to the interpreter, replace it with your own Python).
2. Put the raw inputs where the scripts expect them: `<ROOT>/addenda/opus55_final_audit_tic4206066/00_data/{mast,tesscut,prf,catalogs}/` for the fresh downloads (file names and URLs in `01_inputs/INPUT_PRODUCTS.csv`), and `<ROOT>/addenda/deep_followup_tic4206066/{lc,pixels}/` for the earlier copies (byte-identical; same table).
3. Keep the folder layout of the `from_*` copies: each `from_<step>` folder mirrors `<ROOT>/addenda/<package>/<step>/`, and scripts import their siblings (`common.py`, `fap_common.py`, `t2_common.py`, `b_lib.py`) from the same folder.
4. Paths under `/private/tmp/claude-501/...` are the two auxiliary Python environments used by the audits (`10_environment/pip_freeze_trice_env.txt`, `pip_freeze_audit_env.txt`); rebuild them from the freeze files and point `PY` at them.
5. Scripts that already run standalone from the package (no path edits) are listed in section "Scripts that run standalone from this package" below.

Historical logs (`*.log`, `provenance.jsonl`) keep the original paths on purpose; they are records, not inputs.

## What this is

This package supports the RNAAS Research Note on the two transit-like signals in TESS photometry of TIC 4206066 (StKM 1-561, Gaia DR3 3220388198192519424):

- a 3.18 d signal (P = 3.182785 ± 0.000006 d, Rp = 1.4 ± 0.1 R⊕ if on the target);
- a lower-significance 11.13 d signal.

It contains:
- the input product list (TESS light curves, pixel files, TESScut cutouts, PRF models, catalogue queries), with URLs, retrieval times and sha256;
- the per-transit measurements and fits;
- the TRICERATOPS false-positive runs and their split by companion mass;
- the derivation of the companion limits;
- pixel-level localisation;
- the signal-B tests and chronology;
- the exact data behind Figure 1;
- stellar parameters;
- software environments;
- the list of superseded values;
- the reconciled fact table (`12_facts/FINAL_FACTS.csv`), with a map from every number printed in the Note to a file and key in this package.

**The Note version this package matches** (`submission_candidate/` of the author's reconciliation folder, rebuilt 2026-09-25 evening after the pass-3 review):
- `manuscript.tex`, sha256 `1d0c2964996cbc65f0bf954decf71de7c8758390635619ced3908fce52cf2a70`;
- `manuscript.pdf`, sha256 `9f826c9fbf5c0869db38527d0b54c1def6266dc9047af7e16883c130e313fa0f`.

The manuscript itself is not included.

**Tools.** As disclosed in the Note, Claude Code (Anthropic) assisted with code, computations, figures, references and drafting under the author's direction.

## Folders

| folder | content (each folder has its own README.md) |
|---|---|
| `01_inputs/` | Every public input product: MAST product id or URL, retrieval UTC (2026-09-25 11:29-11:37), bytes, sha256 (`INPUT_PRODUCTS.csv`). Includes catalogue query results (Gaia DR3, TIC 8.2, 2MASS, AllWISE, …), retrieval scripts and logs. Raw FITS are listed, not copied. |
| `02_transits_A/` | 3.18 d signal: 23 per-transit times with errors and selection flags (19 used in the timing fit), per-sector depths, the ephemeris fit (P, T0, χ²), shape fits, and a standalone ephemeris reproduction. |
| `03_archival_test_A/` | Archival false-alarm test of the 3.18 d signal: Sector 98-only period, epoch search in Sectors 6 and 32, identical null machinery (20,000 fake-period trials), printed as "≲1% per sector, ≲0.1% jointly"; includes the pass-3 reviewer's measured N_eff. |
| `04_fpp/` | TRICERATOPS 1.1.0: folds, set-up, every run set (census with means, SD, ranges, run counts, NFPP), and the split of the FPP by companion mass (≈90% from > 0.38 Msun; 0.002-0.007 after the post hoc cut). |
| `05_companion_limit_038/` | `DERIVATION.md`: XP overluminosity → ≈0.38 Msun (0.35-0.42 depending on the scatter definition), BP/RP excess for 0.2-2″, circular-orbit host ≈0.52 Msun, e ≳ 0.2-0.6. Scripts and JSON, plus a standalone recomputation. |
| `06_pixels_localisation/` | Difference-image offsets (A and B, per sector and combined), injection error calibration, neighbour census and per-star Δχ², host probabilities, and the verified aperture re-check. |
| `07_signal_B/` | 11.13 d signal: ephemeris, Sector 98 events, archival recovery test (null, FAP 2×10⁻³), post hoc chronology with file times, depths, Sector 110 predictions. |
| `08_figure/` | Figure 1 (`figure.pdf`/`.png`), the plotted values (`figure_values_final.json`) and the plotting script. |
| `09_star/` | Stellar parameters (Mann relations; the ARIADNE SED values are documented, not printed), Gaia RUWE and RV context, rotation/activity notes. |
| `10_environment/` | pip freeze of the three Python environments; Python, TeX and OS versions. |
| `11_superseded/` | `SUPERSEDED.md`: old values, where they appeared, what replaced them, why. |
| `12_facts/` | `FINAL_FACTS.csv` (unchanged), `MANUSCRIPT_NUMBER_AUDIT.csv`, `NOTE_NUMBER_MAP.csv`, `PATH_MAP.csv`, `EXCLUDED_FILES.csv`. |
| `13_supporting_pipelines/` | Earlier independent pipelines cited in several sections (publication audit, publication-prep analysis, novelty/search recount). |
| `SHA256SUMS` | sha256 of every file in the package (check with `shasum -a 256 -c SHA256SUMS`). |

Folders named `from_<source>/` are copies of the author's analysis folders. File contents are unchanged; the excluded files are listed in `12_facts/EXCLUDED_FILES.csv`. Source names:
- `U` = the opus55 final audit;
- `T1`-`T9` = topic reruns of the final reconciliation;
- `companion_check`, `tic4206066_pixel_recheck`, `deep_followup` = earlier stages.

`12_facts/PATH_MAP.csv` translates the author's paths (as written in `FINAL_FACTS.csv`) into package paths.

## Every number printed in the Note → package file and key

Source: `12_facts/NOTE_NUMBER_MAP.csv`, built from `MANUSCRIPT_NUMBER_AUDIT.csv` and the evidence paths of `FINAL_FACTS.csv`. The *file* column gives the first evidence file present in the package; more are listed in the CSV. The *key* column gives the JSON key, column or line to look at. The value quoted in FINAL_FACTS is in the CSV column `final_value`.

The *check* column is the result of an automatic check (`12_facts/COVERAGE_CHECK.json`, see below):
- "found in mapped file": the listed evidence contains a number that rounds to each printed number;
- "found in section folder": a file in the same folder does;
- "(upper bound)": the printed value is an upper limit ("≲1%", "K<1.2"); the file holds the bounded value itself (e.g. 0.0061 for "≲1%"); these rows have `check_type = BOUND` in `NOTE_NUMBER_MAP.csv`;
- "not a fact": citations and labels;
- "meta": acknowledgment text;
- "literature": thresholds quoted from Giacalone et al. 2021;
- "threshold": the RUWE 1.4 (Lindegren et al. 2018) and gof 4 (Gaia DR3 datamodel) binarity thresholds; the target's own values are in the mapped file;
- "placeholder": the Zenodo DOI (now 10.5281/zenodo.22967099, reserved).

| # | where | printed | fact id | file | key | check |
|---|---|---|---|---|---|---|
| 1 | Title | TIC 4206066 | STAR.identifiers | `05_companion_limit_038/from_U_06_stellar/data/simbad_ids_StKM1-561.ecsv` | (8 ids) | found in mapped file |
| 2 | Title | Two Transit-like Signals | DATA.search.MES_A | `13_supporting_pipelines/U_10_novelty_protocol/leo_rerun.json` |  | no number |
| 3 | Abstract | TIC 4206066 (StKM 1-561) | STAR.identifiers | `05_companion_limit_038/from_U_06_stellar/data/simbad_ids_StKM1-561.ecsv` | (8 ids) | found in mapped file |
| 4 | Abstract | two transit-like signals | DATA.search.MES_A | `13_supporting_pipelines/U_10_novelty_protocol/leo_rerun.json` |  | no number |
| 5 | Abstract | Sector 98 | DATA.sectors | `02_transits_A/from_U_01_ephemeris_A/out_s13_s110.json` | tesspoint_internal | found in mapped file |
| 6 | Abstract | Sectors 6 and 32 | DATA.sectors | `02_transits_A/from_U_01_ephemeris_A/out_s13_s110.json` | tesspoint_internal | found in mapped file |
| 7 | Abstract | 3.18 d | A.P | `02_transits_A/from_U_02_timing_A/verify/v4_results.json` | ephem.drop_onesided_573_and_m1_R60 (P 3.18278495, eP 5.66e-6 | found in mapped file |
| 8 | Abstract | 23 transits | A.N_all | `02_transits_A/from_U_01_ephemeris_A/out_s03_epochs.csv` | rows 2-34 | found in mapped file |
| 9 | Abstract | ~500 ppm | A.depth_combined | `02_transits_A/from_T2_A_shape/out/r01_depths.json` | audit_s09_beta_excl98e-1.combined_all_events = 497.6 +- 29.4 | found in mapped file |
| 10 | Abstract | 1.4 +- 0.1 R_E | A.Rp | `02_transits_A/from_T2_A_shape/out/r03_derived.json` | Rp_k0.0204_0.0010 = [1.385, 0.080], Rp_k0.0195 = 1.324, Rp_k0.02116 = 1.437 | found in mapped file |
| 11 | Abstract | 0.15-0.38 Msun | FP.remaining_scenario_A | `04_fpp/from_T4_fpp_companion/t4_06_eb_totals.json` |  | found in mapped file |
| 12 | Abstract | FPP = 0.03-0.04 | FPP.A_default | `04_fpp/from_T4_fpp_companion/t4_01_fpp_run_census.csv` |  | found in mapped file |
| 13 | Abstract | 11.13 d | B.P | `07_signal_B/from_U_04_B_existence/out/b08_ephem.json` | ephemerides['S6x2 + S98(0,1) [manuscript-like, 4 events]'].P=11.1327392+-1.49e-5 and ['... | found in mapped file |
| 14 | Abstract | ~2.2 R_E | B.Rp_on_target | `13_supporting_pipelines/PUBLICATION_AUDIT/FINAL_PARAMETERS.json` | signal_B.shape_mann_prior.Rp_Rearth | found in mapped file |
| 15 | Abstract | ~3 sigma | B.significance | `13_supporting_pipelines/pubprep/analysis/out/a3_signalB.json` | T1_link[0] (obs 10.03, n_null_ge_obs 19, fap 0.0019998, null_max 10.77) | found in mapped file |
| 16 | Sec1 | TIC 4206066 (Gaia DR3 3220388198192519424 | STAR.identifiers | `05_companion_limit_038/from_U_06_stellar/data/simbad_ids_StKM1-561.ecsv` | (8 ids) | found in mapped file |
| 17 | Sec1 | T = 9.68 | STAR.Tmag | `05_companion_limit_038/from_U_06_stellar/data/tic82.ecsv` | Tmag, e_Tmag | found in mapped file |
| 18 | Sec1 | 35.6 pc | STAR.distance | `09_star/from_T6_star/out/t6_05_parallax_zp.json` | d_corrected_pc, d_err_stat | found in mapped file |
| 19 | Sec1 | (Gaia Collaboration et al. 2016; 2023) |  | – |  | not a fact |
| 20 | Sec1 | M_Ks = 4.94 | STAR.M_Ks | `05_companion_limit_038/from_U_06_stellar/out/s01_relations.json` | MK.MK_zpcorr_plx=4.93817, MK.sigma_total=0.02906 | found in mapped file |
| 21 | Sec1 | Mann et al. (2015, 2019) |  | – |  | not a fact |
| 22 | Sec1 | R* = 0.623 +- 0.019 Rsun | STAR.R | `05_companion_limit_038/from_U_06_stellar/out/s01_relations.json` | radius_mann2015.MC_scatter_plus_inputs (0.62292 +0.01881/-0.01878) | found in mapped file |
| 23 | Sec1 | M* = 0.606 +- 0.015 Msun | STAR.M | `05_companion_limit_038/from_U_06_stellar/out/s01_relations.json` | mass_mann2019.table6_plus_table7_2.4pct_MC (0.60621 +- 0.01514) | found in mapped file |
| 24 | Sec1 | rho* = 2.51 +0.25/-0.22 rho_sun | STAR.rho | `09_star/from_T6_star/out/t6_04_density_mc.json` |  | found in mapped file |
| 25 | Sec1 | (Vines & Jenkins 2022) |  | – |  | not a fact |
| 26 | Sec1 | Sector 98 | DATA.sectors | `02_transits_A/from_U_01_ephemeris_A/out_s13_s110.json` | tesspoint_internal | found in mapped file |
| 27 | Sec1 | (Huang et al. 2020) |  | – |  | not a fact |
| 28 | Sec1 | 126,246 stars | DATA.search.n_stars | `13_supporting_pipelines/U_10_novelty_protocol/tic_recount.json` | n_tmag_le_11p5 | found in mapped file |
| 29 | Sec1 | (T <= 11.5) | DATA.search.n_stars | `13_supporting_pipelines/U_10_novelty_protocol/tic_recount.json` | n_tmag_le_11p5 | found in mapped file (upper bound) |
| 30 | Sec1 | (Kovacs et al. 2002) |  | – |  | not a fact |
| 31 | Sec1 | (Kunimoto et al. 2025) |  | – |  | not a fact |
| 32 | Sec1 | 10.6, 11.0 | DATA.search.MES_A | `13_supporting_pipelines/U_10_novelty_protocol/leo_rerun.json` |  | found in mapped file |
| 33 | Sec1 | 15 | DATA.search.shortlist | `13_supporting_pipelines/U_10_novelty_protocol/search_protocol_extract.json` | (shortlist_v1.1, amendment_01) | found in mapped file |
| 34 | Sec1 | (Ricker et al. 2015) |  | – |  | not a fact |
| 35 | Sec1 | Sector 6, 30 min PDCSAP | DATA.products.S6 | `01_inputs/from_U_00_data/README.md` | section 2 table | found in mapped file |
| 36 | Sec1 | Caldwell et al. 2020 |  | – |  | not a fact |
| 37 | Sec1 | Sector 32, 2 min | DATA.products.S32 | `01_inputs/from_U_00_data/README.md` | section 2 table | found in mapped file |
| 38 | Sec1 | Jenkins et al. 2016 |  | – |  | not a fact |
| 39 | Sec1 | Sector 98, 200 s | DATA.products.S98_cadence | `07_signal_B/from_T7_data_future/s98_quality_bits.json` | TIMEDEL_s, median_dt_s | found in mapped file |
| 40 | Sec1 | Han & Brandt 2023; Petitpas et al. 2026 |  | – |  | not a fact |
| 41 | Sec1 | Sector 98 | DATA.products.S98_quality | `07_signal_B/from_T7_data_future/s98_quality_bits.json` |  | found in mapped file |
| 42 | Sec2 | 3.18 Day | A.P | `02_transits_A/from_U_02_timing_A/verify/v4_results.json` | ephem.drop_onesided_573_and_m1_R60 (P 3.18278495, eP 5.66e-6 | found in mapped file |
| 43 | Sec2 | (Kreidberg 2015) |  | – |  | not a fact |
| 44 | Sec2 | P = 3.182785 +- 0.000006 d | A.P | `02_transits_A/from_U_02_timing_A/verify/v4_results.json` | ephem.drop_onesided_573_and_m1_R60 (P 3.18278495, eP 5.66e-6 | found in mapped file |
| 45 | Sec2 | T0 = 2460998.736 +- 0.003 BJD_TDB | A.T0 | `02_transits_A/from_U_02_timing_A/verify/v4_results.json` | ephem.drop_onesided_573_and_m1_R60 (T0 2460998.735992, eT0 0.002883) | found in mapped file |
| 46 | Sec2 | 19 transits | A.N_timing | `02_transits_A/from_U_02_timing_A/verify/v4_ephem.py` | L33-36 (set 'two') | found in mapped file |
| 47 | Sec2 | >=70% coverage | A.N_timing | `02_transits_A/from_U_02_timing_A/verify/v4_ephem.py` | L33-36 (set 'two') | found in mapped file |
| 48 | Sec2 | S/N >= 2 | A.N_timing | `02_transits_A/from_U_02_timing_A/verify/v4_ephem.py` | L33-36 (set 'two') | found in mapped file |
| 49 | Sec2 | >=15 min of post-egress data | A.N_timing | `02_transits_A/from_U_02_timing_A/verify/v4_ephem.py` | L33-36 (set 'two') | found in mapped file |
| 50 | Sec2 | chi2/dof ~ 9.5/17 | A.timing_chi2 | `02_transits_A/from_U_02_timing_A/verify/v4_results.json` | ephem.drop_onesided_573_and_m1_R60.chi2=9.51, _R144.chi2=8.51 | found in mapped file |
| 51 | Sec2 | (Foreman-Mackey et al. 2013) |  | – |  | not a fact |
| 52 | Sec2 | 23 transits | A.N_all | `02_transits_A/from_U_01_ephemeris_A/out_s03_epochs.csv` | rows 2-34 | found in mapped file |
| 53 | Sec2 | Claret 2017 |  | – |  | not a fact |
| 54 | Sec2 | b = 0.35 +0.10/-0.14 | A.b | `13_supporting_pipelines/PUBLICATION_AUDIT/FINAL_PARAMETERS.json` | signal_A.shape_mann_prior.b (0.346 +0.096/-0.139) | found in mapped file |
| 55 | Sec2 | T14 = 1.87 +- 0.05 h | A.T14 | `13_supporting_pipelines/PUBLICATION_AUDIT/FINAL_PARAMETERS.json` | signal_A.shape_mann_prior.T14_h 1.874 +0.048/-0.047 | found in mapped file |
| 56 | Sec2 | Rp/R* = 0.0204 +- 0.0010 | A.k | `02_transits_A/from_U_07_A_shape_radii/out/summary_fits.csv` | rows mann_claret_jit_raw_excl-1_8 (0.02062), mann_claret_jit_base_excl-1_8 (0.02034) | found in mapped file |
| 57 | Sec2 | Rp = 1.4 +- 0.1 R_E | A.Rp | `02_transits_A/from_T2_A_shape/out/r03_derived.json` | Rp_k0.0204_0.0010 = [1.385, 0.080], Rp_k0.0195 = 1.324, Rp_k0.02116 = 1.437 | found in mapped file |
| 58 | Sec2 | rho* <~ 3.6 rho_sun (97.7%) | A.rho_free | `13_supporting_pipelines/PUBLICATION_AUDIT/FINAL_PARAMETERS.json` | signal_A.free_density.rho_logprior.p97_7=3.472, rho_uniformprior.p97_7=3.553 | found in mapped file (upper bound) |
| 59 | Sec2 | phase 0.5 | A.phase05 | `13_supporting_pipelines/PUBLICATION_AUDIT/scripts/s10_eb_tests.py` | L37-42 (box returns 1 - sum(w f)/sum(w)), L95-97 (d5 = box(...,0.5,T14) | found in mapped file |
| 60 | Sec2 | -20 +- 32 ppm | A.phase05 | `13_supporting_pipelines/PUBLICATION_AUDIT/scripts/s10_eb_tests.py` | L37-42 (box returns 1 - sum(w f)/sum(w)), L95-97 (d5 = box(...,0.5,T14) | found in mapped file |
| 61 | Sec2 | Delta chi2 = 32 | A.Vshape | `13_supporting_pipelines/PUBLICATION_AUDIT/scripts/s10_eb_tests.py` | L141-161 | found in mapped file |
| 62 | Sec2 | >=21 with red noise | A.Vshape | `13_supporting_pipelines/PUBLICATION_AUDIT/scripts/s10_eb_tests.py` | L141-161 | found in mapped file |
| 63 | Sec2 | Without Sector 6 ... <~1 sigma | A.S6_effect | `02_transits_A/from_T1_A_timing/r1_ephemerides.json` | extra.dropS6:* | found in mapped file (upper bound) |
| 64 | Sec3 | 11.13 Day | B.P | `07_signal_B/from_U_04_B_existence/out/b08_ephem.json` | ephemerides['S6x2 + S98(0,1) [manuscript-like, 4 events]'].P=11.1327392+-1.49e-5 and ['... | found in mapped file |
| 65 | Sec3 | 11.13 d | B.P | `07_signal_B/from_U_04_B_existence/out/b08_ephem.json` | ephemerides['S6x2 + S98(0,1) [manuscript-like, 4 events]'].P=11.1327392+-1.49e-5 and ['... | found in mapped file |
| 66 | Sec3 | ~1000 ppm | B.depth_common | `07_signal_B/from_U_05_B_depth/out/s04_grid.csv` | (BASE row | found in mapped file |
| 67 | Sec3 | Rp ~ 2.2 R_E | B.Rp_on_target | `13_supporting_pipelines/PUBLICATION_AUDIT/FINAL_PARAMETERS.json` | signal_B.shape_mann_prior.Rp_Rearth | found in mapped file |
| 68 | Sec3 | offset (+1 +- 3, -1 +- 3)'' | B.offset | `06_pixels_localisation/from_U_referee_resolve_R3/logs/r2_combine.log` | 'B_S32PC+S98_auditor', 'B_S98only_auditor' | found in mapped file |
| 69 | Sec3 | two clean Sector 98 transits | B.S98_transits | `07_signal_B/from_U_04_B_existence/logs/b01.log` | L46 (epoch 2: b29 1.0, bkg 40166) and epoch-4 line (frac_b30 0.55) | found in mapped file |
| 70 | Sec3 | Sector 6 | B.significance_S6only | `13_supporting_pipelines/pubprep/analysis/out/a3_signalB.json` | T2_rank[0].p_rank=0.0015521 | found in mapped file |
| 71 | Sec3 | (marginally in Sector 32) | B.significance_S32only | `13_supporting_pipelines/pubprep/analysis/out/a3_signalB.json` | T1_link[2].fap=0.0319 | found in mapped file |
| 72 | Sec3 | ~3 sigma | B.significance | `13_supporting_pipelines/pubprep/analysis/out/a3_signalB.json` | T1_link[0] (obs 10.03, n_null_ge_obs 19, fap 0.0019998, null_max 10.77) | found in mapped file |
| 73 | Sec3 | ~2e-3 | B.significance | `13_supporting_pipelines/pubprep/analysis/out/a3_signalB.json` | T1_link[0] (obs 10.03, n_null_ge_obs 19, fap 0.0019998, null_max 10.77) | found in mapped file |
| 74 | Sec3 | Sector 32 ... 3.18 d | B.archival_recovery | `07_signal_B/from_U_05_B_depth/verify/v01_info.json` | (bit-512 cadence 1484.870 | found in mapped file |
| 75 | Sec3 | P = 11.13274 +- 0.00003 d | B.P | `07_signal_B/from_U_04_B_existence/out/b08_ephem.json` | ephemerides['S6x2 + S98(0,1) [manuscript-like, 4 events]'].P=11.1327392+-1.49e-5 and ['... | found in mapped file |
| 76 | Sec3 | T0 = 2461000.865 +- 0.003 BJD_TDB | B.T0 | `07_signal_B/from_U_04_B_existence/out/b08_ephem.json` | ephemerides[...4 events].T0=4000.86552+-0.00158 | found in mapped file |
| 77 | Sec4 | Sector 32 | PRF.method | `06_pixels_localisation/from_U_03_pixels/g1_prf_check.json` |  | found in section folder |
| 78 | Sec4 | (Sector 98; Brasseur et al. 2019) | PRF.method | `06_pixels_localisation/from_U_03_pixels/g1_prf_check.json` |  | found in section folder |
| 79 | Sec4 | (Bell & Higgins 2022) |  | – |  | not a fact |
| 80 | Sec4 | 3.18 d | A.P | `02_transits_A/from_U_02_timing_A/verify/v4_results.json` | ephem.drop_onesided_573_and_m1_R60 (P 3.18278495, eP 5.66e-6 | found in mapped file |
| 81 | Sec4 | (+5 +- 3, 0 +- 3)'' | A.offset | `06_pixels_localisation/from_U_referee_resolve_R3/logs/r2_combine.log` | line 'A_eventscatter_auditor_pos' | found in mapped file |
| 82 | Sec4 | consistent within ~2 sigma | A.offset_significance | `08_figure/original_U_proposed_revision_figure_src/figure_values_v2.json` | key panel_b.combined.A.exact (z_E 1.999663721290341, p 0.1349557197582979) | found in mapped file |
| 83 | Sec4 | <=2.6 sigma across Sector 98 reductions | A.offset_significance | `08_figure/original_U_proposed_revision_figure_src/figure_values_v2.json` | key panel_b.combined.A.exact (z_E 1.999663721290341, p 0.1349557197582979) | found in mapped file (upper bound) |
| 84 | Sec4 | Sector 6 | A.S6_diffimg | `06_pixels_localisation/from_U_03_pixels/s6_dumps.json` | S6_ffi.transits | found in mapped file |
| 85 | Sec4 | 0.4-3.5 h | A.S6_diffimg | `06_pixels_localisation/from_U_03_pixels/s6_dumps.json` | S6_ffi.transits | found in mapped file |
| 86 | Sec4 | Gaia DR3 | A.capable_stars | `06_pixels_localisation/from_U_03_pixels/neighbours_vs_fitted_position.csv` | rows with dreqA_min<=1 | found in mapped file |
| 87 | Sec4 | >~43'' | A.capable_stars | `06_pixels_localisation/from_U_03_pixels/neighbours_vs_fitted_position.csv` | rows with dreqA_min<=1 | found in mapped file |
| 88 | Sec4 | (>10 sigma) | A.capable_sigma | `06_pixels_localisation/from_U_referee_resolve_R3/logs/r2_combine.log` | 'capable-star distance corrected_rounded(+5+-3,0+-3)' | found in mapped file |
| 89 | Sec4 | Gaia DR3 | STAR.RUWE | `05_companion_limit_038/from_U_06_stellar/data/gaia_dr3_source.ecsv` | ruwe | found in mapped file |
| 90 | Sec4 | RUWE = 1.21 | STAR.RUWE | `05_companion_limit_038/from_U_06_stellar/data/gaia_dr3_source.ecsv` | ruwe | found in mapped file |
| 91 | Sec4 | radial-velocity scatter (p = 4e-5) | STAR.Gaia_RV | `05_companion_limit_038/from_U_06_stellar/data/gaia_dr3_source.ecsv` | (radial_velocity, rv_*) | found in mapped file |
| 92 | Sec4 | within 10'' | A.nearest_gaia_neighbour | `06_pixels_localisation/from_U_referee_resolve_R3/logs/r4_panel_d_spec.log` | nearest_other_source | found in mapped file |
| 93 | Sec4 | Way et al. (2026) |  | – |  | not a fact |
| 94 | Sec4 | 0.05-0.07 mag | COMP.XP_overluminosity | `05_companion_limit_038/from_companion_check/out/overluminosity.json` | keys O_obs_mag, local_bin.robust_sigma, limits | found in mapped file |
| 95 | Sec4 | +0.003 mag | STAR.overluminosity | `05_companion_limit_038/from_U_06_stellar/out/way2026_target_row.csv` | absG_dr, pabsG_dr, ol | found in mapped file |
| 96 | Sec4 | (<~0.2'') | COMP.XP_overluminosity | `05_companion_limit_038/from_companion_check/out/overluminosity.json` | keys O_obs_mag, local_bin.robust_sigma, limits | found in mapped file (upper bound) |
| 97 | Sec4 | ~0.38 Msun (2 sigma | COMP.XP_overluminosity | `05_companion_limit_038/from_companion_check/out/overluminosity.json` | keys O_obs_mag, local_bin.robust_sigma, limits | found in mapped file |
| 98 | Sec4 | Pecaut & Mamajek 2013 |  | – |  | not a fact |
| 99 | Sec4 | at 0.2-2'' | COMP.BPRP_excess | `05_companion_limit_038/from_U_referee_resolve_R4/out/r4_02_speckle_calib.json` | keys target, single_companion_sep0.2-2.0_dI_le1.76, single_companion_sep0.15-2.5_dI_le1.76 | found in mapped file |
| 100 | Sec4 | (Riello et al. 2021) |  | – |  | not a fact |
| 101 | Sec4 | 3.18 d | A.P | `02_transits_A/from_U_02_timing_A/verify/v4_results.json` | ephem.drop_onesided_573_and_m1_R60 (P 3.18278495, eP 5.66e-6 | found in mapped file |
| 102 | Sec4 | ~0.52 Msun | COMP.host_mass_circular | `05_companion_limit_038/from_U_09_companion/out/companion_space.json` |  | found in mapped file |
| 103 | Sec4 | 0.15-0.38 Msun | COMP.host_lower_bound_and_ecc | `05_companion_limit_038/from_U_referee_resolve_R4/out/r4_summary.json` | key e_ecc | found in mapped file |
| 104 | Sec4 | e >~ 0.2-0.6 | COMP.host_lower_bound_and_ecc | `05_companion_limit_038/from_U_referee_resolve_R4/out/r4_summary.json` | key e_ecc | found in mapped file |
| 105 | Sec4 | TRICERATOPS 1.1.0 | FPP.A_default | `04_fpp/from_T4_fpp_companion/t4_01_fpp_run_census.csv` |  | found in mapped file |
| 106 | Sec4 | (Giacalone et al. 2021) |  | – |  | not a fact |
| 107 | Sec4 | Sector 32+98 folds | FPP.A_default | `04_fpp/from_T4_fpp_companion/t4_01_fpp_run_census.csv` |  | found in mapped file |
| 108 | Sec4 | (5-10 min bins) | FPP.A_default | `04_fpp/from_T4_fpp_companion/t4_01_fpp_run_census.csv` |  | found in mapped file |
| 109 | Sec4 | at 0.1-2.2'' | FP.remaining_scenario_A | `04_fpp/from_T4_fpp_companion/t4_06_eb_totals.json` |  | found in mapped file |
| 110 | Sec4 | FPP = 0.03-0.04 | FPP.A_default | `04_fpp/from_T4_fpp_companion/t4_01_fpp_run_census.csv` |  | found in mapped file |
| 111 | Sec4 | (4-20 runs of 10^6 draws) | FPP.A_runs | `04_fpp/from_T4_fpp_companion/t4_01_fpp_run_census.csv` | (n_runs column counted from list lengths / probs_run files) | found in mapped file |
| 112 | Sec4 | about 90% | FPP.A_share_above_0.38 | `04_fpp/from_T4_fpp_companion/t4_02_checks.json` | keys share_gt038_A_default, share_gt038_A_range | found in mapped file |
| 113 | Sec4 | 0.38 Msun | FPP.A_share_above_0.38 | `04_fpp/from_T4_fpp_companion/t4_02_checks.json` | keys share_gt038_A_default, share_gt038_A_range | found in mapped file |
| 114 | Sec4 | FPP ~ 0.002-0.007 | FPP.A_constrained | `04_fpp/from_U_08_fpp/out_fpp_summary.csv` | column FPP_comp_mass_le_0.38_mean | found in mapped file |
| 115 | Sec4 | (FPP < 0.015, NFPP < 1e-3 | FPP.validation_rule | – |  | literature |
| 116 | Sec4 | NFPP < 1e-6 | FPP.A_NFPP | `04_fpp/from_T4_fpp_companion/t4_02_checks.json` | keys NFPP_max_all_default_A, A_default_single_sector, A_default_S32S98_summary.NFPP_max | found in mapped file (upper bound) |
| 117 | Sec5 | 3.18 d | A.P | `02_transits_A/from_U_02_timing_A/verify/v4_results.json` | ephem.drop_onesided_573_and_m1_R60 (P 3.18278495, eP 5.66e-6 | found in mapped file |
| 118 | Sec5 | 11.13 d ... ~3 sigma | B.status_wording | – |  | found in section folder |
| 119 | Sec5 | Sector 110 (2026 October 31-November 26) | FUTURE.S110_dates | `07_signal_B/from_T7_data_future/web/TESS_orbit_times.csv` | L334-L336 | found in mapped file |
| 120 | Sec5 | eight or nine 3.18 d transit times (one in a scheduled data gap) and three 11.13 d ones | FUTURE.S110_counts | `07_signal_B/from_T7_data_future/s110_predict.json` | ephemerides.*.counts and events[].status | found in mapped file |
| 121 | Sec5 | As of 2026 September 25 (UTC) | NOVELTY.sources_queried | `13_supporting_pipelines/U_10_novelty_protocol/queries.log` |  | found in mapped file |
| 122 | Caption | (a) 3.18 d | A.P | `02_transits_A/from_U_02_timing_A/verify/v4_results.json` | ephem.drop_onesided_573_and_m1_R60 (P 3.18278495, eP 5.66e-6 | found in mapped file |
| 123 | Caption | 11.13 d and one Sector 98 edge transit removed | A.depth_S98 | `02_transits_A/from_T2_A_shape/out/r01_depths.json` | audit_s09_beta_excl98e-1.per_sector.S98 = 495.9 +- 47.7 (8 events) | found in section folder |
| 124 | Caption | inset: central 22'' |  | – |  | not a fact |
| 125 | Caption | 1 sigma rings (3.18 d: Sectors 32+98; 11.13 d: Sector 98) | B.offset | `06_pixels_localisation/from_U_referee_resolve_R3/logs/r2_combine.log` | 'B_S32PC+S98_auditor', 'B_S98only_auditor' | found in mapped file |
| 126 | Caption | 3.18 d depth (crosses; a fifth at 84'') | A.capable_stars | `06_pixels_localisation/from_U_03_pixels/neighbours_vs_fitted_position.csv` | rows with dreqA_min<=1 | found in mapped file |
| 127 | Caption | (c) |  | – |  | not a fact |
| 128 | Caption | three 6-min-bin folds | FPP.A_default | `04_fpp/from_T4_fpp_companion/t4_01_fpp_run_census.csv` |  | found in mapped file |
| 129 | Caption | capped at 0.38 Msun (post hoc; 5-min fold) | FPP.A_constrained | `04_fpp/from_U_08_fpp/out_fpp_summary.csv` | column FPP_comp_mass_le_0.38_mean | found in mapped file |
| 130 | Ack | GPT-5.6 Sol and GPT-6 Pro | META.chatgpt_models | – |  | meta |
| 131 | Sec2 | false-alarm probability <~1% per sector, <~0.1% jointly | A.archival_FAP | `03_archival_test_A/a_fap_results.json` | keys OPTIMISED, nulls, printable | found in mapped file (upper bound) |
| 132 | Sec4 | RUWE (1.4) and gof (4) thresholds | STAR.Gaia_RV | `05_companion_limit_038/from_U_06_stellar/data/gaia_dr3_source.ecsv` | (radial_velocity, rv_*) | threshold |
| 133 | Sec4 | K<1.2 km/s | COMP.Gaia_RV | `05_companion_limit_038/from_U_09_companion/out/rv_limits.json` |  | found in mapped file (upper bound) |
| 134 | Ack | doi:10.5281/zenodo.22967099 | META.zenodo | – |  | reserved DOI |
| 135 | Ack | GPT-5.6 Sol, GPT-6 Pro | META.chatgpt_models | – |  | meta |

**Coverage check result** (`12_facts/COVERAGE_CHECK.json`; method in its `summary.method` field):
- 135 rows in total. 108 rows have every printed number found in package files: 104 of them in the mapped evidence file itself, the rest in another file of the same section folder.
- 8 of the OK rows are printed upper limits, matched by the BOUND rule.
- The remaining rows: 20 citations/labels, 2 wording rows without digits ("two … signals"), 2 acknowledgment rows, 1 literature threshold (TRICERATOPS validation), 1 binarity-threshold row (RUWE 1.4, gof 4) and 1 DOI placeholder.
- **Gaps: 0.**

**Limitations of the check.**
- It confirms that a matching number is present in the file. It does not re-derive the number.
- Small integers (sector numbers, counts) match easily.
- The ARIADNE clause (R = 0.62 ± 0.02 Rsun, Teff = 4020 ± 60 K) was removed from the Note in the pass-3 revision; its facts stay in the package as documented, not printed (`09_star/README.md`).

## How to reproduce

**1. Get the inputs.** Download the products in `01_inputs/INPUT_PRODUCTS.csv` (column `url`) and check them against the `sha256` column. The retrieval scripts in `01_inputs/retrieval_scripts/` (astroquery 0.4.11) do this automatically.

**2. Scripts that run standalone from this package** (tested 2026-09-25 with the main `.venv`, Python 3.12.13):

```
cd 02_transits_A/reproduce_ephemeris && python v4_ephem.py
    # 19-transit ephemeris; output equals the stored v4_results.json
cd 03_archival_test_A && python s4_results.py
    # re-assembles a_fap_results.json from the stored s1/s2 outputs (no raw data needed)
python 04_fpp/recompute_share_038.py
    # share of FPP from companions > 0.38 Msun: 0.837-0.951; FPP after cut 0.0020-0.0067
python 05_companion_limit_038/recompute_mass_limit.py
    # 2σ mass limit 0.381 Msun (σ = 0.068), 0.349-0.420 for σ = 0.052-0.10
```

**3. All other scripts** are copied unchanged for transparency. They contain absolute paths to the author's project tree (`/Users/paulsmith/space <3/tess-jev/...`) and read the raw FITS inputs. To rerun them, download the inputs and edit the path constants at the top of each script (`ROOT`, `U`, `D`, `REV`, …). Each `from_*` folder has the original `REPORT.md` or `result.json`, which records how its step was run.

**Environments:**
- main analysis: Python 3.12.13 (numpy 2.5.3, scipy 1.18.1, astropy 8.0.1, batman 2.5.3, emcee 3.1.6, wotan 1.10, TESS_PRF 0.1.3);
- TRICERATOPS: Python 3.11.15 (triceratops 1.1.0, numpy 1.26.4, lightkurve 2.5.1);
- audit env: Python 3.12.13 (+ dynesty 3.1.0, tess-point 0.9.5.post1).

Full lists are in `10_environment/`. OS: macOS 26.3.1 (arm64). Manuscript build: Tectonic 0.17.0, AASTeX 7.0.2; word count with TeXcount 3.1.1.

## Not included (and why)

Every excluded file is listed with size and sha256 in `12_facts/EXCLUDED_FILES.csv`:
- **Raw TESS FITS and TESScut cubes** (up to 466 MB each). Download them from MAST via `01_inputs/INPUT_PRODUCTS.csv`.
- **MCMC chains and derived light-curve arrays.** They can be recreated with the scripts.
- **Python pickles:** the TRICERATOPS targets and the ARIADNE posteriors (up to 36 MB).
- **Large logs and MOLUSC draw tables** (36-70 MB).
- **Third-party papers, documentation snapshots and code repositories.** They are cited instead.
- **Bulk catalogue tables of the novelty check** (TOI/CTOI lists, SPOC TCE tables, QLP target list).
- **Internal search-cohort bookkeeping** of the author's wider TESS search. A minimal extract is kept in `13_supporting_pipelines/U_10_novelty_protocol/search_protocol_extract.json`.

The package is 41.7 MB (before SHA256SUMS) in 1888 files (including SHA256SUMS; 1887 entries in SHA256SUMS).
