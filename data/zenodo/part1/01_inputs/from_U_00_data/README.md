# 00_data: fresh public data for the TIC 4206066 audit

Retrieved 2026-09-25, 11:29–11:37 UTC. All files were downloaded fresh from the public archives by the scripts in `scripts/`. No earlier analysis code or derived tables were used.
`DATA_MANIFEST.json` has one entry per file (173 entries). Each entry gives the relpath, URL, UTC time, sha256, bytes, key FITS header values and a `matches_existing` field. The raw per-download log is `logs/provenance.jsonl`.

**Summary.** Every MAST light curve and pixel file is byte-identical to the earlier copies in `deep_followup_tic4206066/lc` and `tic4206066_PUBLICATION_AUDIT/data/mast`. The TESScut 21x21 arrays are identical to both earlier copies. All 50 PRF files are byte-identical. The Gaia DR3 target rows match the earlier copies value for value. Nothing in the public data changed between the earlier work and this audit. One thing downstream agents need to know: in the QLP Sector 32 light curve, QUALITY bit 6 (value 32) does not mark the real momentum dumps (see section 6). Use the SPOC/TESS-SPOC/TESScut flags for dump times.

Python: `audit_env` (astroquery 0.4.11, astropy 8.0.1, numpy 2.5.3, requests 2.34.2). No packages were installed.
Target: TIC 4206066 = Gaia DR3 3220388198192519424. Query centre (80.85983, -1.32696) is the Gaia DR3 position at epoch 2016.0. The TIC/QLP J2000 position is (80.8598410, -1.3272783), about 1.1" south, because pmdec is +70.8 mas/yr.

## 1. Layout

| dir | content | script |
|---|---|---|
| `mast/` | time-series FITS (11 files) + full MAST observation/product listings | `a_mast_query.py`, `a2_mast_download.py` |
| `tesscut/21x21/`, `tesscut/31x31/` | TESScut FFI cutouts for S6, S32, S98 | `b_tesscut.py` |
| `prf/start_s0004/cam1_ccd1/`, `cam1_ccd2/` | full 25-file PRF sets | `c_prf.py` |
| `catalogs/` | Gaia DR3, Bailer-Jones 2021, TIC 8.2, VizieR photometry (ECSV + CSV + `.query.txt`) | `d_catalogs.py` |
| `docs/` | official documentation snapshots (raw + `.txt` extraction) | `g_docs.py` (PRF docs in `c_prf.py`) |
| `logs/` | `provenance.jsonl`, `fits_headers.json`, `momentum_dumps_from_quality.{json,csv}`, `comparison_with_existing.json`, job logs | `f_headers.py`, `h_compare.py` |

`i_manifest.py` builds `DATA_MANIFEST.json`.

## 2. MAST inventory (a)

Queries: a 10" cone search plus a target_name search for the TIC and Gaia ids. This returned 26 observations, 12 of them time series, and 19 products. Full lists: `mast/mast_observations_all.csv`, `mast/mast_products_all.csv`, `mast/mast_products_listing_urls.json`.

Downloaded (all byte-identical to both earlier copies):

| file | sector | cadence | CAM/CCD | pipeline / version |
|---|---|---|---|---|
| hlsp_tess-spoc_..._s0006_tess_v1_lc.fits / _tp.fits | 6 | 1800 s | 1/1 | spoc-5.0.20-20201120, DATA_REL 42 |
| tess2020324010417-s0032-...-0200-s_lc.fits / _tp.fits | 32 | 120 s | 1/2 | spoc-5.0.20-20201120, DATA_REL 48 |
| hlsp_tess-spoc_..._s0032_tess_v1_lc.fits / _tp.fits | 32 | 600 s | 1/2 | spoc-5.0.21-20210107, DATA_REL 48 |
| hlsp_qlp_tess_ffi_s0006-..._llc.fits | 6 | 1800 s | 1/1 | QLP, CALIB='TICA', TICVER 8 |
| hlsp_qlp_tess_ffi_s0032-..._llc.fits | 32 | 600 s | 1/2 | QLP, CALIB='TICA v0.2.1', TICVER 8.1 |
| hlsp_qlp_tess_ffi_s0098-..._llc.fits | 98 | 200 s | 1/2 | QLP, CALIB='TICA v1.4.0', TICVER 8.2 |
| hlsp_gsfc-eleanor-lite_..._s0006_..._v1.0_lc.fits | 6 | 1800 s | 1/1 | eleanor VERSION 2.0.3 |
| hlsp_tglc_..._gaiaid-3220388198192519424-s0006-..._llc.fits | 6 | 1800 s | 1/1 | TGLC |

The MAST Download API redirected the SPOC 2-min files to `stpubdata.s3.us-east-1.amazonaws.com`. The manifest records both the requested URL and the final URL.

Listed but not downloaded (URLs in `mast_products_listing_urls.json`; earlier copies are in `deep_followup_tic4206066/lc`): T16 S6, TASOC S6 (cbv + ens), TARS S6 and S32 (plus preview JPEGs), and the QLP S6 `.txt` light curve. Non-time-series observations are ELEANOR postcards (S6), TICA FFIs (S32, S98), SPOC FFIs (S6, S32, S98), GALEX AIS, PS1, and SDSS images.

Not available at MAST as of 11:30 UTC:
- no SPOC 2-min data for S6 or S98;
- no TESS-SPOC light curve for S98;
- no SPOC DV/TCE products for S32; only the LC and TP were listed;
- eleanor-lite, TGLC, T16 and TASOC exist for S6 only.

TESScut `get_sectors` lists only S6 (1-1), S32 (1-2) and S98 (1-2). Sector 110 is not listed yet.

## 3. TESScut (b)

Downloaded with `astroquery.mast.Tesscut.download_cutouts` (SPOC FFIs, `inflate=False`). The zip sha256 was recorded and each FITS was then extracted. The exact URLs are in the manifest, e.g. `https://mast.stsci.edu/tesscut/api/v0.1/astrocut?x=21&y=21&units=px&sector=98&ra=80.85983&dec=-1.32696`.

| sector | CAM/CCD | astrocut (PROCVER) | cube DATE | DATA_REL | n_cad | TIMEDEL | 21x21 bytes | 31x31 bytes |
|---|---|---|---|---|---|---|---|---|
| 6 | 1/1 | 0.11.1 | 2021-08-12 | 42 | 993 | 1800 s | 8.9 MB | 19.2 MB |
| 32 | 1/2 | 0.11.1 | 2021-02-03 | 48 | 3600 | 600 s | 32.0 MB | 69.4 MB |
| 98 | 1/2 | 0.11.1 | 2026-01-13 | 134 | 24162 | 200 s | 214.7 MB | 466.0 MB |

Comparison with earlier copies:
- **21x21 cutouts:** TIME, TIMECORR, CADENCENO, RAW_CNTS, FLUX, FLUX_ERR, FLUX_BKG, FLUX_BKG_ERR, QUALITY and POS_CORR1/2 are identical to both earlier copies (`deep_followup_tic4206066/pixels`, and `tic4206066_PUBLICATION_AUDIT/data/tesscut`, whose centre differed by about 0.01"). The WCS keywords are also the same.
- **31x31 cutouts:** no earlier copy exists. Their central 21x21 pixels are identical to the fresh 21x21 cutouts.
- **Keyword caveat:** CROWDSAP and FLFRCSAP are undefined in the cutouts, which is expected for FFI cutouts.

## 4. PRF (c)

Camera/CCD from the headers: S6 = cam 1 / CCD 1; S32 and S98 = cam 1 / CCD 2.

Which set applies: the MAST TESS Archive Manual states that `start_s0001/` applies to Sectors 1–3 and `start_s0004/` to Sector 4 and beyond (https://outerspace.stsci.edu/spaces/TESS/pages/14563420/2.0+-+Data+Product+Overview; snapshot `docs/spoc/mast_tess_archive_manual_2.0_data_product_overview.html`). The archive's own `start_s0004/00README.txt` says the same (snapshot `docs/prf/prf_index_start_s0004_00README.txt`). The FITS keyword is VERSION='UPDATED_2.0'.

Downloaded: all 25 `*-row????-col????.fits` files for `cam1_ccd1` and `cam1_ccd2` from https://archive.stsci.edu/missions/tess/models/prf_fitsfiles/start_s0004/. All 50 are byte-identical to `tic4206066_PUBLICATION_AUDIT/data/prf/`.

Note: the archive has only these two PRF sets. They are built from 2018 July commissioning data, so S98 (2025) necessarily uses the same `start_s0004` models. The README says interpolation is accurate to 1–5% and that PRF column coordinates are offset by 44 columns from FFI/TPF coordinates.

## 5. Catalogues (d, e)

| file (catalogs/) | source | rows |
|---|---|---|
| gaiadr3_gaia_source_cone240as | ESA TAP `SELECT *` from gaiadr3.gaia_source, 240" cone, plus `dist_arcsec_j2016pos` | 167 (153 columns) |
| gaiadr3_nss_two_body_orbit / nss_acceleration_astro / nss_non_linear_spectro / nss_vim_fl _target | ESA TAP, source_id = target | 0 each |
| gaiadr3_vari_summary_target | ESA TAP | 0 |
| gaiadr3_astrophysical_parameters_target | ESA TAP `SELECT *` | 1 (226 columns) |
| bailerjones2021_I352 | VizieR I/352/gedr3dis, Source = target | 1 |
| tic82_cone300as | MAST Catalogs TIC, 300" | 254 (125 columns) |
| 2mass_psc_II246 / allwise_II328 / catwise2020_II365 / tycho2_I259 / apass_dr9_II336 | VizieR, 30" cone, all columns + `_r` | 1 / 1 / 2 / 1 / 1 |
| galex (II/335/galex_ais) | VizieR, 30" | **no rows**; no file written |

These are raw target values for convenience. None of them was evaluated here:
- **Gaia DR3 astrometry and photometry:**
  - parallax 28.0750 ± 0.0176 mas;
  - pmra/pmdec −3.406 / +70.806 mas/yr;
  - G 10.513, BP−RP 1.739;
  - RUWE 1.2143, astrometric_excess_noise 0.0935 mas (sig 15.76);
  - ipd_gof_harmonic_amplitude 0.0222, ipd_frac_multi_peak 0, ipd_frac_odd_win 0;
  - duplicated_source False, non_single_star 0, phot_variable_flag NOT_AVAILABLE, visibility_periods_used 21.
- **Gaia DR3 radial velocity:**
  - RV 10.49 ± 0.22 km/s, rv_nb_transits 33, rv_amplitude_robust 3.22 km/s;
  - **rv_chisq_pvalue 4.06e-5, rv_renormalised_gof 3.00**, rv_template_teff 4250 K.
- **Bailer-Jones 2021:** rgeo 35.579 pc (35.557–35.596).
- **2MASS:** J/H/Ks 8.527/7.932/7.694 (±0.029/0.034/0.029), Qflg AAA, _r = 1.0".
- **AllWISE:** W1–W4 7.605/7.626/7.56/7.451, ccf 0.
- **Tycho-2:** TYC 4753-1591-1, BT/VT 13.424/11.241.
- **APASS:** V 11.126, B 12.481.

Comparison with earlier copies:
- **Gaia target rows:** all 152 common columns of the gaia_source row and all 226 astrophysical_parameters columns agree with `tic4206066_PUBLICATION_AUDIT/data/gaia/*.csv`. Some agree only to float32 print precision because the older CSV was written with fewer digits.
- **Gaia 180" neighbours:** the same 96 source_ids as the earlier list. The NSS and vari_summary tables are empty in both.
- **TIC:** the target row is equal apart from `dstArcSec`, which depends on the query centre. The earlier 0.0834-degree cone also contained TIC 4206079, which lies at 300.09" and so falls outside the 300" radius used here.

## 6. Key FITS header facts (f)

Full per-file values are in `logs/fits_headers.json` and in the manifest `header` blocks. All products have TIMESYS=TDB, TIMEREF=SOLARSYSTEM, TIMEUNIT=d and BJDREFI=2457000 (so BTJD = BJD_TDB − 2457000).

| product | TSTART–TSTOP (BTJD) | INT_TIME / TIMEDEL | CROWDSAP / FLFRCSAP | flux columns |
|---|---|---|---|---|
| TESS-SPOC S6 LC | 1468.2766–1490.0469 | 1.98 s / 1800 s | 0.99504 / 0.92497 | SAP_FLUX, PDCSAP_FLUX (+ errors, SAP_BKG, centroids) |
| SPOC S32 2-min LC | 2174.2329–2200.2304 | 1.98 s / 120 s | 0.99695 / 0.92180 | SAP_FLUX, PDCSAP_FLUX |
| TESS-SPOC S32 LC | 2174.2273–2200.2345 | 1.98 s / 600 s | 0.99291 / 0.92180 | SAP_FLUX, PDCSAP_FLUX |
| QLP S6 | 1468.2870–1490.0365 | 1800 s | – | SAP_FLUX, KSPSAP_FLUX(+_ERR,_SML,_LAG) |
| QLP S32 | 2174.2308–2200.2311 | 600 s | – | SAP_FLUX, KSPSAP_FLUX(+_ERR,_SML,_LAG) |
| QLP S98 | 3988.4989–4045.9020 | 200 s | – | SAP_FLUX, **DET_FLUX**(+_ERR,_SML,_LAG), **SYS_RM_FLUX** |
| TESScut S6/S32/S98 | see table above | 1.98 s | undefined | RAW_CNTS, FLUX, FLUX_ERR, FLUX_BKG(+_ERR) |

Header conventions that differ between products:
- QLP and TGLC write `TELAPS` and `BJDREFR` instead of `TELAPSE` and `BJDREFF`.
- eleanor-lite's `TELAPSE` (0.0208 d) is the cadence length, not the elapsed time.
- `EXPOSURE` is the total livetime in days in SPOC LC/TP files but the per-cadence value in days in TESScut files.
- The QLP S98 format differs from S6 and S32: it has DET_FLUX instead of KSPSAP_FLUX and an extra SYS_RM_FLUX column.

### Momentum dumps from QUALITY bit 6 (value 32)

The full table is in `logs/momentum_dumps_from_quality.csv`. Times are BTJD; add 2457000 for BJD_TDB.

- **S6** (TESS-SPOC LC/TP, TESScut, QLP and eleanor-lite agree): 1468.3912, 1471.5161, 1474.6411, 1481.2868, 1484.4117, 1487.5366.
  - One 30-min FFI is flagged per dump, so each dump time is known only to that 30-min window. Timestamps differ between products by up to 2e-4 d.
  - This matches the DRN: dumps every 3.125 d.
- **S32, SPOC 2-min:** 2179.4865–2179.4948 (7 cadences), 2186.9449–2186.9491 (4 cadences; the first cadences of orbit 72, which starts at 2186.93981 per the DRN), and 2193.4867–2193.4950 (7 cadences).
- **S32, 10-min FFI** (TESScut and TESS-SPOC): 2179.4740–2179.4949 and 2193.4742–2193.4950. The DRN says there was one dump halfway through each orbit.
- **S98** (TESScut 200 s; QLP flags the same cadences): 3988.4987, 3995.8208, 4003.0109, 4010.0366, 4017.0992, 4024.5251, 4031.8029 and 4038.8444, with 2–6 cadences each.
  - Each dump is at the start of a data segment, right after a 5.1-h downlink gap.
  - This matches the DRN: two dumps per orbit (orbits 205–208), performed at the end of each downlink.
- **QLP S32 problem:** bit 6 is set on 8 isolated cadences (2175.50, 2179.65, 2181.97, 2184.12, 2187.74, 2192.86, 2196.84, 2197.93). None of them is a SPOC-flagged dump. At the real dumps (2179.474–2179.495 and 2193.474–2193.495) QLP sets 2048 or 4096, or nothing.
  - **Do not use QLP S32 QUALITY&32 to locate momentum dumps.**
  - QLP S32 is built from TICA v0.2.1 FFIs. Its CADENCENO does not match the SPOC cadence numbers, and its TIME differs from TESScut by up to 4.3 s.
- **Bit 13 means different things in QLP and SPOC.** In QLP, bit 13 (value 4096) means "low precision points". The QLP data-product description misprints the value as "4196". In SPOC, bit 13 (value 4096) means "Scattered Light Exclude" (SDPDD Table 32).
- **QLP S98 bit 12 (2048, straylight):** set on 5757 of 24162 cadences. The DRN and HEASARC notes say Earth scattered light affects Camera 1 near the start of each orbit, and the Moon affects it in orbits 205 and 207.

## 7. Documentation snapshots (g)

Each file under `docs/` has a `.txt` extraction. URL, UTC time and sha256 are in the manifest and in `docs/docs_download_log.json`.

- **SPOC quality bits:** SDPDD Rev F, §9 Table 32. Sources: https://archive.stsci.edu/files/live/sites/mast/files/home/missions-and-data/active-missions/tess/_documents/EXP-TESS-ARC-ICD-TM-0014-Rev-F.pdf, plus the archive copy `.../missions/tess/doc/EXP-TESS-ARC-ICD-TM-0014-Rev-F-FINAL.pdf`. Bit 6 = 32 = reaction-wheel desaturation.
- **DRNs:**
  - S6: https://archive.stsci.edu/missions/tess/doc/tess_drn/tess_sector_06_drn08_v02.pdf
  - S32: .../tess_sector_32_drn48_v02.pdf
  - S98: .../tess_sector_98_drn134_v01.pdf
  - the DRN index;
  - HEASARC condensed DRNs: https://heasarc.gsfc.nasa.gov/docs/tess/data_release_notes.html
- **HEASARC S98 summary (scattered light):** https://heasarc.gsfc.nasa.gov/docs/tess/sector98_summary.html. The sector table is also saved.
- **QLP:**
  - https://archive.stsci.edu/hlsp/qlp
  - the readme: https://archive.stsci.edu/hlsps/qlp/hlsp_qlp_tess_ffi_all_tess_v1_readme.txt
  - quality flags, §2.3.3: https://archive.stsci.edu/hlsps/qlp/hlsp_qlp_tess_ffi_all_tess_v1_data-prod-desc.pdf
  - https://tess.mit.edu/qlp/
- **PRF:** the archive manual page, the MAST data-products page, both PRF directory indexes and both `00README.txt` files.
- **TESScut:** https://mast.stsci.edu/tesscut/docs/getting_started.html
- **Failed:** https://heasarc.gsfc.nasa.gov/docs/tess/observing-technical.html returned HTTP 404 and was not needed.

Six doc snapshots are byte-identical to the earlier copies in `tic4206066_PUBLICATION_AUDIT/data/refs`: the QLP readme, the QLP data-prod-desc, the MAST QLP page, the MIT QLP page, the HEASARC S98 summary and the HEASARC sector table.

## 8. Problems and notes for downstream agents

1. **(NOTE, data quality)** Momentum-dump flags in QLP S32 are misplaced (section 6). Anyone who took S32 dump times from QLP would get wrong times. The SPOC, TESS-SPOC and TESScut flags agree with each other and with the DRN.
2. **(NOTE)** S6 dump times are known only to the 30-min FFI window (±15 min). If the manuscript's statement "every transit falls 0.4–3.5 h from a reaction-wheel desaturation" needs better precision, a S6 2-min light curve of any other target would help, because dumps are spacecraft-wide. That file was not downloaded because it was not requested.
3. **(NOTE)** No PRF model exists for the extended mission. S98 uses the 2018 `start_s0004` model, which is the only official option.
4. **(NOTE)** The Gaia DR3 RV diagnostics above (p-value 4.1e-5, renormalised gof 3.0, amplitude 3.2 km/s over 33 transits) bear on the manuscript's claim that the RVs "exclude an eclipsing stellar companion". They should be evaluated by the companion/stellar agents against the Gaia DR3 RV-variability criteria (Katz et al. 2023). They were not evaluated here.
5. The SPOC and TESS-SPOC products (S6, S32) carry TICVER 8.1. QLP S98 and the TIC query are 8.2. The target's Tmag is 9.6817 in both.
6. Nothing changed in the archive relative to the earlier copies (section 2–5). Downstream differences from earlier analyses therefore cannot be blamed on different input data.
