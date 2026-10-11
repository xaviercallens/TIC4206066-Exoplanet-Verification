# 13_supporting_pipelines: earlier independent pipelines cited as evidence in several sections

| folder | what it is | used in the Note for |
|---|---|---|
| `PUBLICATION_AUDIT/` | Independent publication audit (2026-09-25, before the opus55 final audit). Scripts `s01`–`s23`, outputs `out/`, logs, plots, regression tests, `FINAL_PARAMETERS.json`, `DATA_PROVENANCE.csv`, `PUBLICATION_TRUTH.md`. | **Shape of the 3.18 d signal:** b, T14 and ρ limit (`FINAL_PARAMETERS.json` → `signal_A.shape_mann_prior`, `signal_A.free_density`), with the fitting script `scripts/s08_globalfit.py`. **Eclipsing-binary tests:** phase-0.5 depth and V-shape Δχ² (`out/s10_eb_tests.json`, `scripts/s10_eb_tests.py`). **Per-event depths:** `out/s09_events.json`. **B radius:** `signal_B.shape_mann_prior`. **FPP runs:** `out/fpp/runs/*/probs_runNN.csv`, `out/fpp/inputs/`, `scripts/s13_*`. |
| `pubprep/analysis/` | The author's publication-prep analysis. `a1_timing_depths.py` (timing, depths; imported by the figure script), `a2_prf_localisation.py`, `a3_signalB_tests.py` (B link test, S6 rank, S98-only), `a4_fpp_ownfold.py`, `a5_summaries.py`, `a6_figure.py`; outputs in `out/`. | **B significance:** FAP 2×10⁻³ (`out/a3_signalB.json` → `T1_link[0]`); Sector 6 / Sector 32 alone. **B depth:** `out/a5_summary.json`. **FPP ownfold runs:** `out/fpp/`. Per-transit rows (`out/a1_per_transit.csv`) re-used by the reconciliation reruns. |
| `U_10_novelty_protocol/` | Discovery-search recount and novelty queries of the opus55 audit. `tic_recount.py` (126,246 stars with T ≤ 11.5), `leo_rerun.py` (LEO-Vetter MES 10.6 / 11.0), `fresh_queries.py` with `queries.log`, `ads_queries.json`, `novelty_summary.json` and `raw/` query results (ExoFOP target page, NASA Exoplanet Archive, SIMBAD, VizieR, arXiv, T16/RAVEN records). `search_protocol_extract.json` is a minimal extract of the shortlist rule (see below). | "126,246 stars", "MES 10.6, 11.0", "shortlist threshold (15)", and "As of 2026 September 25 … no prior report". |
| `extra_science_e4_injection_recovery/` | Injection-recovery completeness of a simple BLS search in the Sector 98 QLP light curve (`e4_injection_recovery.py`, `e4_completeness_table.json`, `e4_injections.json`). | Not cited in the Note. Included for context. |

**Not included:**
- `PUBLICATION_AUDIT/data/` (raw MAST/TESScut/PRF mirror; byte-identical to `../01_inputs` products);
- `out/molusc/` and `out/phot/`, plus the MCMC chains;
- bulk catalogue tables of the novelty check (ExoFOP TOI/CTOI lists, SPOC TCE tables, RAVEN/T16 full tables, QLP S98 target list; each 2-37 MB);
- `protocol_recount.{json,py,stdout}`: search-cohort bookkeeping and a remote-host path, replaced by `search_protocol_extract.json`.

All are listed with sha256 in `../12_facts/EXCLUDED_FILES.csv`.

**Paths.** The scripts in these folders use absolute paths to the author's project tree. They are provided for transparency and do not run standalone from this package.
