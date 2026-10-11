# Audit 03_pixels — REPORT

> Rendered verbatim from the agent's structured result (`_structured/audit__03_pixels.json`); the sub-agent harness blocked agents from writing .md report files. Scripts, outputs and logs are in this folder.


**Verdict:** PASS_WITH_CHANGES

## Summary

Verdict: PASS WITH CHANGES. The qualitative localisation claims hold. The 3.18 d position number does not reproduce, the 11.13 d Dec value is transcribed wrongly, and Figure 1d contradicts the text on Sector 6.

REPORT.md was not written: this sub-agent is not allowed to write report .md files, so the full report is in this structured result. The scripts, JSON/CSV outputs, logs and two figures are in /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels.

What I did:
- All code is new. The only inputs were the 00_data files; all 64 match the DATA_MANIFEST sha256 values (input_provenance.json).
- G1, PRF files:
  - Each MAST PRF file is 117x117 = 9x9 interleaved 13x13 sub-PRFs, as the MAST TESS Archive Manual and start_s0004/00README say. Numerically the file sums to 81 and each sub-PRF to 1.
  - TESS_PRF.locate moves the centroid by the requested (dx, dy) to rms 0.01 px (max 0.036 px).
  - Its interpolation across the 5x5 grid of files equals my own bilinear interpolation to 1e-16.
  - start_s0004 applies to S6, S32 and S98.
- G2, coordinates:
  - Gaia DR3 positions were propagated to each sector's epoch.
  - A global WCS shift was fitted to the out-of-transit scene and is at most 0.054 px (about 1").
  - The target is never counted as its own neighbour.
- Difference images: one per transit, with a per-pixel linear baseline plus a box term, then a PRF fit with free x, y, amplitude and background.
- G3, injections:
  - 300 model-free injections into the real pixels per sector, product and signal, using the real events' cadence patterns.
  - For off-target hypotheses, the real target image shifted by 1–2 whole pixels; the shifted source is recovered at 0.99–1.02 of the true offset.
- G4, neighbour flux: each star's flux fraction inside the actual apertures (SPOC/TESS-SPOC 18-pixel apertures; the QLP S98 3x3 BESTAP, which is TGLC-based from Sector 94), TIC 8.2 magnitudes or Stassun+2019 eqs. 1–2, and the CROWDSAP correction. My aperture model reproduces SPOC's own FLFRCSAP/CROWDSAP.

Main result:
- **S32 pointing artefact.** During the eight 3.18 d transits of S32, the SPOC pointing correction POS_CORR2 is displaced by +0.37 ± 0.16 millipixels. That is a 2.7 sigma outlier against random epochs, and it predicts a spurious offset of (+8.9, −12.2)".
  - After decorrelating the pixels against POS_CORR, S32 moves from (+9.1, −8.3)" to (+5.2, +0.4)".
  - The signal-to-noise ratio at the target doubles, from 3.5 to 7.7.
- **S98 edge event.** The package's S98 value of −9" in Dec is driven by one transit at BTJD 3995.555 that ends at a downlink gap. Dropping it changes Dec from −8.1" to −1.6".
- **3.18 d position.** The headline S32 (decorrelated) + S98 (clean, flags 2048 + 2^29) position is (+5.5 ± 2.2, +0.2 ± 2.4)", chi2 6.0 for 2 dof (p = 0.049). With the injection bias removed it is (+4.5 ± 2.2, +0.4 ± 2.4)", p = 0.13.
  - Without decorrelation the position is (+5.7 ± 3.2, −0.9 ± 4.3)", p = 0.16.
  - "Consistent at 2 sigma" is honest but borderline.
  - All three sectors show a consistent eastward offset of about 5". The 11.13 d signal shows none.
- **11.13 d position.** (+1.1 ± 2.2, −1.2 ± 2.9)". The package's own value is (+1.4 ± 3.3, +0.1 ± 4.4)", so the manuscript's "+1 ± 4" in Dec is a rounding or transcription error.
- **Sector 6.** "0.4–3.5 h" is confirmed: 0.41–3.54 h, with ±0.25 h resolution from the flags.
  - The dump pointing steps are 1–2 millipixels, and the raw S6 difference image is corrupted.
  - After decorrelation S6 agrees with the target, but the result is unstable when single transits are dropped. "Unreliable" fits better than "unusable".
  - Figure 1d nevertheless plots an S6 point for the 3.18 d signal, 3.0 sigma east of the target. Adding that point would give p = 7e-4, so excluding S6 decides whether "consistent at 2 sigma" holds.
- **Neighbours.**
  - "≥43"" is confirmed: the nearest capable star is TIC 4206077 at 43.1".
  - "≥7 sigma" is confirmed: 7.15 from the manuscript's own position, and at least 18 from mine.
  - The nearest Gaia star (11", Tmag 18.13) cannot produce the depth: it would need a 120–123% eclipse.
  - The Tmag 13.06 star at 84" can produce the depth only through its PRF wings, and it lies outside the plotted range of Fig. 1d.

## Numeric rows

| quantity | manuscript | recomputed | status | method / note |
|---|---|---|---|---|
| 3.18 d source offset, East (S32+S98) | +8 +- 5 arcsec | +5.5 +- 2.2 (S32 POS_CORR-decorrelated + S98 clean29); +4.5 +- 2.2 bias-corrected; +5.7 +- 3.2 without decorrelation | MINOR_DIFF | per-transit per-pixel linear-baseline difference images, TESS_PRF fit, 300 model-free on-target injections per set, inverse-covariance combination The package's own a5 value +7.9 +- 5.3 rounds correctly to +8 +- 5. My errors are smaller because the S98 injection SD is 3.3" and S32 is 3.1" after decorrelation. |
| 3.18 d source offset, North (S32+S98) | -8 +- 7 arcsec | +0.2 +- 2.4 (decorrelated); +0.4 +- 2.4 bias-corrected; -0.9 +- 4.3 without decorrelation | DISCREPANT | as above; jackknife and pointing diagnostics The -8" comes from (i) POS_CORR2 pointing jitter correlated with the S32 transit windows (+0.37 +- 0.16 millipixels; 2.7 sigma vs random epochs) and (ii) the S98 transit at BTJD 3995.555, which ends at a data gap. It lies within 1.2 sigma of the manuscript's own error. |
| 3.18 d position consistent with target | consistent at 2 sigma | chi2 = 6.0/2 dof, p = 0.049 (decorrelated); p = 0.13 bias-corrected; p = 0.024 with robust covariance; p = 0.16 without decorrelation; p = 0.10 with decorrelated S6 added | MATCH | Mahalanobis chi2 with the injection covariance Borderline. A consistent eastward offset of about 5" appears in S6, S32 and S98. |
| Sectors used for the 3.18 d localisation | Sectors 32 and 98 | S32 (SPOC 2-min TPF) + S98 (TESScut 200 s, 2048 + QLP 2^29 excluded) | MATCH | reproduced the package's set choice  |
| S6 transit distance to nearest reaction-wheel desaturation | 0.4-3.5 h (every transit) | 0.76, 2.15, 3.54, 3.19, 1.80, 0.41 h (6 transits in the data; the 7th is in the mid-sector gap) | MATCH | QUALITY bit 6 (value 32) in the TESScut and TESS-SPOC S6 files; S6 DRN gives dumps every 3.125 d The 30-min FFI flags give only +-0.25 h resolution. |
| Sector 6 usability for the 3.18 d localisation | unusable | raw: (-8.7, -34.7)", p < 0.003 vs injections (corrupted); +-0.5 h dump excision: (+6.3 +- 4.0, -1.3 +- 4.9)"; POS_CORR-decorrelated: (+1.0 +- 2.7, -1.9 +- 4.0)", but a jackknife swings it by up to 40" | MINOR_DIFF | difference imaging with dump excision, a step model, or POS_CORR decorrelation 'Unreliable' is more accurate than 'unusable'. |
| Figure 1d S6 3.18 d point | (+8.2 +- 2.7, -2.4 +- 3.6)" (a5_summary.json, plotted as a filled circle) | not reproducible; depends on dump treatment: (-8.7, -34.7), (+6.3, -1.3), (+1.0, -1.9) | DISCREPANT | as above Plotted although the text says S6 is unusable. It is 3.0 sigma east in the package's own numbers. |
| Figure 1d S98 3.18 d point | (+7.0 +- 5.6, -9.3 +- 7.3)" (a5) | (+5.6 +- 3.3, -0.3 +- 4.5)" (9 events, clean29) | DISCREPANT | injection-calibrated PRF difference image My rebuild of the package's baseline recipe gives (+1.7, -8.1); dropping the event at 3995.555 gives (+5.2, -1.6). |
| Figure 1d S32 2-min 3.18 d point | (+13.5 +- 17.2, -2.3 +- 21.6)" (a5) | raw (+9.1 +- 13.8, -8.3 +- 15.3); POS_CORR-decorrelated (+5.2 +- 3.1, +0.4 +- 2.8) | MINOR_DIFF | as above  |
| Nearest Gaia DR3 star able to produce the 3.18 d depth | >= 43 arcsec | 43.1" (TIC 4206077, Tmag 17.02; required eclipse depth 0.67 / 0.90 / 1.11 in the S6 / S32 / QLP-S98 apertures) | MATCH | PRF flux fraction in the actual apertures, TIC 8.2 Tmag, CROWDSAP-corrected depth definition, delta = 500 ppm, capable if required depth <= 1 Five capable stars (43, 50, 54, 56 and 84"). If the star must work in all three apertures, or with a 0.5 threshold, only the Tmag 13.06 star at 84" qualifies. |
| Sigma-distance of capable stars from the fitted position | >= 7 sigma | min 7.15 (manuscript position and errors; 6.8 sigma 2-D Gaussian equivalent); min 18.3 (audit position) | MATCH | Mahalanobis distance with full covariance  |
| 11.13 d source offset, East | +1 +- 3 arcsec | +1.1 +- 2.2 (S32 decorrelated + S98 clean29); +0.5 +- 2.0 with S6 added | MATCH | as for the 3.18 d signal  |
| 11.13 d source offset, North | +1 +- 4 arcsec | -1.2 +- 2.9 (audit); the package's own S32+S98 value is +0.1 +- 4.4 | MINOR_DIFF | as above The package's own value rounds to 0 +- 4, not +1. The '+1' matches neither the package nor this audit. |
| Gaia DR3 stars able to produce the 11.13 d depth | (not stated) | only TIC 4206057 (Tmag 13.06, 84"; required depth 0.61 in the S32 SPOC aperture), >= 37 sigma from the 11.13 d position | NOT_CHECKED | G4 census, delta = 1000 ppm Informational; the manuscript makes no claim here. |
| Nearest Gaia star (11", Tmag 18.13) maximum depth | (not stated; package: at most 409 ppm) | required eclipse depth 1.20-1.23 at 500 ppm, i.e. at most about 410 ppm; 5.6 sigma from the audit position | MATCH | G4 census  |

## Findings

### F1 [MINOR] 3.18 d position (+8 +- 5, -8 +- 7)" does not reproduce; the -8" Dec part is a pointing and single-event artefact

An independent analysis gives (+5.5 +- 2.2, +0.2 +- 2.4)" from S32 (POS_CORR-decorrelated) plus S98 (clean), and (+5.7 +- 3.2, -0.9 +- 4.3)" without decorrelation. Two causes produce the -8" Dec: SPOC pointing jitter in S32 that happens to line up with the transit windows, and one S98 transit that ends at a downlink gap. The conclusion 'consistent with the target at about 2 sigma' survives (p = 0.02-0.16 across variants), and the capable-star exclusion becomes stronger. It is rated MINOR rather than MAJOR because the quoted errors cover the corrected value within about 1.2 sigma. The text and Fig. 1d should still be updated.

*Evidence:* S32 2-min POS_CORR2 in-transit shift +0.37 +- 0.16 millipixels (random-epoch SD 0.135), predicting (+8.9, -12.2)"; raw S32 offset (+9.1, -8.3)" (pointing.json). After decorrelation S32 = (+5.2, +0.4)" and S/N 3.5 -> 7.7 (poscorr_decorr_S32_tpf.json). S98 with the package-like baseline = (+1.7, -8.1)"; dropping BTJD 3995.555 gives (+5.2, -1.6)" (diagnose_s98.json). Combined values in summary.json, keys A:S32tpfPC+S98clean29:*.

*Location:* Sec. 4, sentence 3; Fig. 1d

*Proposed fix:* Decorrelate against POS_CORR (or drop the compromised events) and quote (+5 +- 2, 0 +- 2)" (chi2 = 6.0/2 dof). Otherwise state that the Dec component is sensitive to pointing jitter and baseline choice at the 5-8" level.

### F2 [MINOR] 11.13 d Dec offset '+1 +- 4' does not match the package's own number

The package's S32+S98 result for the 11.13 d signal is (+1.4 +- 3.3, +0.1 +- 4.4)", which rounds to (+1 +- 3, 0 +- 4)". The manuscript's '+1' in Dec seems to come from a different analysis (the internal audit's +1.4). This audit finds (+1.1 +- 2.2, -1.2 +- 2.9)". All versions are consistent with the target.

*Evidence:* rnaas_tic4206066_pubprep/analysis/out/a5_summary.json localisation.B; MANUSCRIPT_NUMBER_AUDIT.csv row 52; this audit summary.json B:S32tpfPC+S98clean29:bias0:rob0 = (+1.12 +- 2.21, -1.23 +- 2.94), chi2 0.52.

*Location:* Sec. 4, 'The 11.13 d signal lies at (+1 +- 3, +1 +- 4)"'

*Proposed fix:* Write (+1 +- 3, 0 +- 4)" (package) or (+1 +- 2, -1 +- 3)" (this audit).

### F3 [MINOR] Fig. 1d plots an S6 3.18 d offset that the text calls unusable, and excluding S6 decides the headline

Panel (d) draws the S6 3.18 d point, (+8.2 +- 2.7, -2.4 +- 3.6)", with the smallest error bars in the panel, although the text says S6 is unusable. In the package's numbers that point is 3.0 sigma east of the target, and adding it gives chi2 = 14.5 (p = 7e-4). So the S6 exclusion decides whether 'consistent at 2 sigma' holds. The exclusion itself is justified: each momentum dump shifts the pointing by 1-2 millipixels inside every transit window, and the raw S6 difference image is corrupted. Once pointing is decorrelated, S6 agrees with the target (+1.0 +- 2.7, -1.9 +- 4.0)", though a jackknife over single transits is unstable.

*Evidence:* a6_figure.py panel (d) loops over per_set including S6_ffi; a5_summary.json S6_ffi z = (+3.02, -0.58), combined_all_independent chi2 14.55. This audit: raw S6 (-8.7, -34.7)", p < 0.003 vs injections; dump centroid steps 0.4-1.8 millipixels (pointing.json); decorrelated S6 in poscorr_decorr_S6_tspoc_S6_tspoc_dumpx0.5h_S32_tspoc.json.

*Location:* Fig. 1d and caption; Sec. 4

*Proposed fix:* Drop the S6 3.18 d point from panel (d), or plot the decorrelated value and mark it as excluded. Word the text as 'S6 is not used because a dump-induced pointing step falls inside every transit window' ('unreliable' rather than 'unusable').

### F4 [NOTE] '0.4-3.5 h from a desaturation' confirmed

All six S6 transits in the data fall 0.41, 0.76, 1.80, 2.15, 3.19 and 3.54 h from the nearest quality-bit-6 dump. The dumps come every 3.125 d (S6 DRN Sec. 1.2), and P_A drifts only 1.39 h per 3.125 d relative to them.

*Evidence:* s6_dumps.json (TESScut and TESS-SPOC LC flags agree); https://archive.stsci.edu/missions/tess/doc/tess_drn/tess_sector_06_drn08_v02.pdf

*Location:* Sec. 4

*Proposed fix:* None.

### F5 [NOTE] 'All Gaia DR3 stars able to produce the depth lie >= 43" away, >= 7 sigma' confirmed

This follows G4: flux fractions in the real S6/S32 SPOC apertures and the QLP S98 3x3 aperture, TIC 8.2 magnitudes, CROWDSAP-corrected depths, and delta = 500 ppm. Five stars can reach a required eclipse depth <= 1 in at least one aperture; the nearest is TIC 4206077 at 43.1". The capable stars lie at Mahalanobis distance >= 7.15 from the manuscript's own position (6.8 sigma 2-D Gaussian equivalent) and >= 18.3 from this audit's. Only the Tmag 13.06 star at 84" is capable in all apertures or with a 0.5 threshold. The nearest Gaia star (11", Tmag 18.13) would need a 120-123% eclipse.

*Evidence:* neighbour_census_per_product.csv, neighbours_vs_fitted_position.csv, summary.json 'neighbours'; the model reproduces SPOC FLFRCSAP/CROWDSAP (0.918/0.9953 vs 0.925/0.9950 in S6; 0.916/0.9936 vs 0.922/0.9970 in S32).

*Location:* Sec. 4

*Proposed fix:* Optional: say that '7 sigma' is a Mahalanobis distance.

### F6 [NOTE] Consistent ~5" eastward offset of the 3.18 d signal in all three sectors

S6 (decorrelated or dump-excised), S32 (decorrelated) and S98 each place the 3.18 d source 1-6" east, about 2 sigma combined. The 11.13 d signal shows no offset. Injections show a method bias of about +2" E in S98 only. The pixel data cannot exclude an unresolved source within about 5-10" of the target; that question needs Gaia contrast limits or high-resolution imaging.

*Evidence:* summary.json sets S98_clean29_A, S32_tpf_poscorr_A, S6_tspoc_poscorr_A and S6_dumpx0.5h_A; combined S6+S32+S98 decorrelated = (+3.6 +- 1.7, -0.6 +- 2.0)", p = 0.10.

*Location:* Sec. 4

*Proposed fix:* Optionally add 'a ~5" eastward offset (<~2 sigma) is not significant'.

### F7 [NOTE] S98 scattered-light definition is undocumented

Only QLP bit 2048 plus bit 2^29 reproduces the manuscript's 3 clean + 1 scattered-light 11.13 d transits in S98. Bit 2^29 is not described in the QLP data-product documentation or on the MIT QLP page, which document 2048 and 2^30; in flagged cadences the sky is 672 e/s against 135 e/s. The localisation is insensitive to the choice: all data, 2048-only and 2048+2^29 agree within 1 sigma.

*Evidence:* 00_data/docs/qlp/*.txt; https://tess.mit.edu/qlp/ ; real_results.json S98_clean / S98_clean29 / S98_all.

*Location:* Sec. 1 ('Sector 98 scattered-light intervals removed')

*Proposed fix:* State the flag definition, e.g. in the data section.

### F8 [NOTE] Fig. 1d caption: capable star at 84" not shown; Bell2022 bib entry has an empty title

TIC 4206057 (Tmag 13.06, 84"; required depth 0.31 in the S32 SPOC aperture) can produce the depth through its PRF wings but lies outside the plotted x-range (x_lim 62 to -76; the star is at dE = +66"). The bib entry Bell2022 has title = {\relax}, which is for the reference audit.

*Evidence:* a6_figure.py line 69 (nb filtered to sep < 62"); references.bib lines 189-195.

*Location:* Fig. 1 caption; references.bib

*Proposed fix:* Extend the axes or add '(one more at 84" off-panel)'.

### F9 [NOTE] Weak 11.13 d signal in the S32 pixels (for audit 05)

The S32 11.13 d box depth at the target is 330 ppm raw (S/N 1.9) and 561 ppm after decorrelation (S/N 4.8), against the manuscript's S32 depth of 1003 +- 206 ppm. There are only two events, and the one at 2175.096 lies 2.3 h from a 3.18 d transit (2175.001), which masks its pre-transit baseline.

*Evidence:* real_results.json S32_tpf B; poscorr_decorr_S32_tpf.json

*Location:* Sec. 3 per-sector depths

*Proposed fix:* Cross-check in the depth audit.

## Recommended wording

I fitted difference images with the TESS pixel response function (Bell & Higgins 2022) after removing spacecraft pointing jitter with the SPOC POS_CORR vectors; uncertainties come from transits injected into the real pixels. For the 3.18 d signal, Sectors 32 and 98 place the source at (+5 +- 2, 0 +- 2)" from the target (chi2 = 6.0 for 2 dof, consistent at 2 sigma). Sector 6 is not used, because every transit falls 0.4-3.5 h from a reaction-wheel desaturation whose pointing step dominates the difference image; after decorrelation it agrees with the target. All Gaia DR3 stars able to produce the depth lie >= 43" away, > 7 sigma from the fitted position. The 11.13 d signal lies at (+1 +- 2, -1 +- 3)". [Minimal-change alternative: keep the 3.18 d numbers, change the 11.13 d position to (+1 +- 3, 0 +- 4)", and remove the Sector 6 3.18 d point from Fig. 1d.]

## Files

- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/common.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/pixloc.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/g1_prf_check.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/g1_prf_check.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/run_real.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/real_results.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/real_results_S6_tspoc_S6_tspoc_dumpx0.5h_S32_tspoc.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/real_images.npz
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/inject.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/inj/
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/summarize.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/summary.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/neighbours.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/neighbour_census.csv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/neighbour_census_per_product.csv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/neighbour_census_summary.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/neighbours_vs_fitted_position.csv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/neighbours_dchi2.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/neighbours_dchi2.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/pointing.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/pointing.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/poscorr_decorr.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/poscorr_decorr_S32_tpf.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/poscorr_decorr_S6_tspoc_S6_tspoc_dumpx0.5h_S32_tspoc.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/diagnose_s98.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/diagnose_s98.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/s6_dumps.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/s6_dumps.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/make_figures.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/fig_localisation.png
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/fig_diffimages.png
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/input_provenance.json