# 00_data: provisioning report

Scope: data provisioning only; no signal analysis. Full details, tables and URLs are in `README.md`. Per-file provenance (URL, UTC, sha256, bytes, headers, match with earlier copies) is in `DATA_MANIFEST.json`.

What was done:
- Listed all MAST observations and products for TIC 4206066 (26 observations, 19 products).
- Downloaded 11 time-series FITS: the TESS-SPOC S6 LC and TP; the SPOC S32 2-min LC and TPF; the TESS-SPOC S32 LC and TP; the QLP S6, S32 and S98 LCs; eleanor-lite S6; TGLC S6.
- Downloaded 6 TESScut cutouts: S6, S32 and S98, each at 21x21 and 31x31 (astrocut 0.11.1).
- Downloaded 50 PRF files: start_s0004 for cam1/ccd1 and cam1/ccd2.
- Queried Gaia DR3 (240" cone plus NSS, vari_summary and AP tables for the target), Bailer-Jones 2021, TIC 8.2 (300" cone), 2MASS, AllWISE, CatWISE2020, Tycho-2, APASS9 and GALEX (no match).
- Saved 17 official documentation snapshots.
- Extracted the header keywords and momentum-dump lists.

Comparison with earlier copies: every MAST product and every PRF file is byte-identical. The TESScut 21x21 arrays are numerically identical. The Gaia and TIC target rows match. The input data did not change.

Findings (data level):

| # | severity | finding |
|---|---|---|
| 1 | NOTE | In the QLP S32 light curve, QUALITY bit 6 (value 32) is set on 8 cadences that are not momentum dumps and is absent at the two real dumps (BTJD 2179.47–2179.50 and 2193.47–2193.50). Take dump times from SPOC, TESS-SPOC or TESScut flags, which agree with each other and with the DRN. |
| 2 | NOTE | Bit 13 (value 4096) means "low-precision point" in QLP but "Scattered Light Exclude" in SPOC. The QLP document misprints the value as 4196. |
| 3 | NOTE | S6 has only 30-min FFI data, so dump times are known only to ±15 min. This limits how precisely the manuscript's "0.4–3.5 h from a desaturation" can be checked. |
| 4 | NOTE | The only official PRF models are from 2018 (start_s0004). S98 has to use them. |
| 5 | NOTE | Gaia DR3 gives rv_chisq_pvalue 4.1e-5, rv_renormalised_gof 3.0, rv_amplitude_robust 3.2 km/s (33 transits), RUWE 1.214 and astrometric_excess_noise_sig 15.8. These are passed to the companion/stellar agents for evaluation against the manuscript's RV claim; they were not evaluated here. |
| 6 | NOTE | Header conventions differ: QLP and TGLC use TELAPS and BJDREFR; the eleanor-lite TELAPSE holds the cadence length; the QLP S98 flux columns are DET_FLUX and SYS_RM_FLUX, not KSPSAP_FLUX. |
