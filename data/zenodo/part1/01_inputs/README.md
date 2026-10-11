# 01_inputs: public input data

This folder lists every public product the analysis used: TESS light curves, target pixel files, TESScut cutouts, PRF models and catalogue queries. Large binary products are **not** included. Each one is listed with its archive URL, retrieval time (UTC) and sha256, so it can be downloaded again and checked byte for byte.

| file | content |
|---|---|
| `INPUT_PRODUCTS.csv` | One row per file retrieved on 2026-09-25 between 11:29 and 11:37 UTC. Columns: category, file name, MAST product id / VizieR table, sector, camera/CCD, cadence, download URL (the requested URL too, if MAST redirected), retrieval UTC, bytes, sha256 (plus the zip sha256 for TESScut), and whether the file is in this package. Built from `from_U_00_data/DATA_MANIFEST.json`. |
| `from_U_00_data/` | `DATA_MANIFEST.json` (173 entries, including HTTP headers and key FITS header values), and `README.md` / `REPORT.md` of the retrieval step. |
| `retrieval_scripts/` | The scripts that downloaded everything: `a_mast_query.py`, `a2_mast_download.py`, `b_tesscut.py`, `c_prf.py`, `d_catalogs.py`, `g_docs.py`, plus the header/compare/manifest helpers. They use astroquery 0.4.11. |
| `retrieval_logs/` | Per-download provenance (`provenance.jsonl`), FITS header dump (`fits_headers.json`), momentum-dump times from QUALITY bit 6 (`momentum_dumps_from_quality.csv/json`), and the comparison with earlier copies. |
| `catalog_queries/` | Catalogue query results as CSV/ECSV plus the exact query text (`*.query.txt`). Included: Gaia DR3 `gaia_source` (240″ cone), NSS tables, `vari_summary` and `astrophysical_parameters` for the target; Bailer-Jones et al. 2021 distances; TIC 8.2 (300″ cone); 2MASS, AllWISE, CatWISE2020, Tycho-2, APASS DR9. |
| `reference_tables/mamajek_table.txt` | Pecaut & Mamajek (2013) dwarf colour/Teff/mass table, version 2022.04.16, as used for the companion mass limit. The table header asks users to cite Pecaut & Mamajek (2013, ApJS 208, 9). |

## Light curves and pixels used in the Note

| Sector | product used | file name (MAST) |
|---|---|---|
| 6 | TESS-SPOC FFI light curve, 1800 s, PDCSAP (cam 1 / CCD 1) | `hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits` (TP: `..._tp.fits`) |
| 32 | SPOC 2-min target light curve, PDCSAP (cam 1 / CCD 2) | `tess2020324010417-s0032-0000000004206066-0200-s_lc.fits` (TPF: `..._tp.fits`) |
| 98 | QLP FFI light curve, 200 s, all QLP-flagged cadences removed (cam 1 / CCD 2) | `hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits` |
| 6, 32, 98 | TESScut SPOC-FFI cutouts, 21×21 and 31×31 pixels | `tess-s00NN-1-C_80.859830_-1.326960_{21x21,31x31}_astrocut.fits` |
| – | PRF models, `start_s0004`, cam1_ccd1 and cam1_ccd2 (25 files each) | `tess2019107181900-prf-1-1-*.fits`, `tess2019107181901-prf-1-2-*.fits` |

The exact URLs and sha256 values are in `INPUT_PRODUCTS.csv`. The TESScut URL pattern is
`https://mast.stsci.edu/tesscut/api/v0.1/astrocut?x=21&y=21&units=px&sector=98&ra=80.85983&dec=-1.32696`.

Some analysis scripts read earlier local copies of the same files, from the folders `deep_followup_tic4206066/lc` and `tic4206066_PUBLICATION_AUDIT/data`. The manifest's `matches_existing` field records that these copies are byte-identical to the 2026-09-25 downloads.

**Discovery search input (not included).** The Sector 98 QLP light curves of the 126,246 stars with T ≤ 11.5 are available from the QLP high-level science product at MAST. The target list (`qlp_s0098_target_list.csv`, 28.7 MB) is listed with its sha256 in `../12_facts/EXCLUDED_FILES.csv`. The search configuration and recount are in `../13_supporting_pipelines/U_10_novelty_protocol/`.

**Not redistributed.** Documentation snapshots (MAST, HEASARC, QLP and SPOC pages) and third-party papers are not included. Their URLs are in `from_U_00_data/DATA_MANIFEST.json` (entries `docs/...`).
