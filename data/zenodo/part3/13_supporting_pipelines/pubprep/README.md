# RNAAS manuscript, TIC 4206066 — publication-prep package (NOT SUBMITTED)

Nothing was submitted, uploaded or published. The frozen search pipeline, configs, cohorts and earlier evidence were only read.
This directory supersedes `../rnaas_tic4206066_draft/` (kept as history).

## Deliverables
| file | what |
|---|---|
| `manuscript.pdf`, `manuscript.tex`, `references.bib` | the Note (AASTeX v7.0.2, `aastex702.cls`, `aasjournalv7.bst`); 3 pages, 1488 words (texcount), abstract 132 |
| `figure.pdf`, `figure.png` | the single 5-panel figure (`analysis/a6_figure.py`) |
| `MANUSCRIPT_NUMBER_AUDIT.csv` | every number: source artifact → recalculated value → manuscript value (67 rows; `audit/number_audit.py`) |
| `REFERENCE_AUDIT.md` | every citation resolved (CrossRef/ASCL records in `citations/`) and what it establishes |
| `METHOD_AUDIT.md` | methods, primary sources, and the comparison between this re-analysis and the earlier audit |
| `CHANGES_FROM_CURRENT_DRAFT.md` | every number/statement changed since the previous draft and why |
| `UNRESOLVED_LIMITATIONS.md` | L1–L14, ordered by impact |
| `REFEREE_REPORT_1.md`, `REFEREE_REPORT_2.md` | two mock hostile-referee passes; pass 1 changed the text |
| `WORDCOUNT.txt` | official texcount command and breakdown |
| `outputs_sha256.json` | sha256 of every file in this directory (`tools/hashes.py`) |

## Independent re-analysis (written from scratch here, `analysis/`)
`a1_timing_depths.py` (per-transit times/depths, injection-calibrated errors, ephemerides, B sector depths) →
`out/a1_results.json`, `out/a1_per_transit.csv`, `out/a1_oc.png`; `a2_prf_localisation.py` (official TESS PRF,
model-free injections) → `out/a2_localisation_*.json`, `out/a2_diffimages_*.png`; `a3_signalB_tests.py` (B link null
test, S6 rank test, S98-only test) → `out/a3_signalB.json`; `a4_fpp_ownfold.py` (TRICERATOPS with a third fold) →
`out/fpp/`; `a5_summaries.py` → `out/a5_summary.json`; `a6_figure.py`. Logs `out_a*.log`. Environment:
`tools/env_freeze.txt`.

## Before any submission
1. Fill the placeholders in `manuscript.tex`: affiliation city/country, e-mail (mandatory in AASTeX v7), ORCID if any.
   Word margin is 12 (texcount).
2. Re-run `python3 novelty_refresh/refresh.py` (ExoFOP TOI/CTOI, target page, NASA toi, exo.MAST) and check `refresh_result.json`.
3. Rebuild: `sh tools/build.sh` (Tectonic; the figure needs the analysis env: `analysis/a6_figure.py`).
4. Optional if the editorial word count allows: add back Kunimoto et al. (2022, RNAAS 6, 236) for the QLP S98 data
   and Stassun et al. (2019) for the TIC magnitude (see REFERENCE_AUDIT.md).
