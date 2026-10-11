# Audit 04_B_existence — REPORT

> Rendered verbatim from the agent's structured result (`_structured/audit__04_B_existence.json`); the sub-agent harness blocked agents from writing .md report files. Scripts, outputs and logs are in this folder.


**Verdict:** PASS_WITH_CHANGES

## Summary

Signal B (P ≈ 11.13 d) holds up when tested from scratch on the fresh public data. The "≈3σ, low-significance candidate" framing is supported and a little conservative, but several numbers and phrases need fixing.

Sector 98 (QLP, quality=0, A masked):
- Two clean transits: E0 1085±138 ppm (S/N 7.9) and E1 1189±154 ppm (S/N 7.7).
- A 'partial' E4: 951±243 ppm, S/N 3.9. It is covered across the whole transit, but 55% of its in-transit cadences carry QLP's internal poor-quality flag (bit value 2^30), interleaved with good cadences.
- A scattered-light E2: 1212±200 ppm.
- E−1 and E3 fall in 2048-straylight intervals and cannot be used.

S98-only search:
- BLS and TLS find B as the global maximum in QLP DET_FLUX + biweight: BLS SDE 12.3, TLS SDE 13.6.
- It is rank 16 (BLS) and rank 4 (TLS) in own TESScut photometry.
- The 3.71 d and 5.57 d peaks in candidate.json are B's P/3 and P/2 aliases.
- Odd/even differ by 1.6σ; there is no secondary eclipse.

S98-only null (segment-cyclic shifts, B removed): the FAP depends on the method. It is 1.6e-3 (SDE) or 2e-4 (S/N) for QLP DET; 0.37 (SDE) or 0.010 (S/N) for SYS_RM, which reproduces the manuscript's 0.35; and about 1 for own photometry.

Archival link test (S6+S32, pipeline):
- My statistic is 9.75 at P = 11.132767 (manuscript: 10.0).
- 0 of 50,000 null trials reach it, with B's archival windows kept or removed. There are 36,471 distinct configurations; an exponential tail fit gives about 5e-6.
- The result holds for 6 photometries and for windows ±0.004, ±0.0045 and ±0.0113 d.
- S6 alone: 0/50,000, but only 1,418 distinct configurations, so FAP ≲ 7e-4 (resolution-limited). S32 alone: 0.02.
- The other independent S98 peaks (9.97 d, 1.06 d) do not link (FAP 0.14–0.49).
- The manuscript's 2e-3 is therefore conservative; its heavier null tail comes from box artifacts at Sector 32 data edges.
- The "kept" null is the correct permutation null; the manuscript's "removed" null is formally anti-conservative, although here both give 0/50,000.
- The ±0.004 d scan window is not derived from the S98-only period error (σ_P 0.0015–0.0038 d; S98-only estimates span 11.1316–11.1416 d).

Sector 6 dips:
- They are the two strongest 1.6-h boxes in PDCSAP, SAP, QLP and eleanor-lite, under all 4 detrendings.
- They are not the top two in TGLC aperture or own photometry, where artifacts at momentum dumps and gap edges rank higher.
- The TGLC PSF light curve is corrupted near 1471–1474 BTJD.
- No background or centroid anomalies; momentum dumps are 21.8 h and 10.9 h away.
- Rank-test p ≈ 2e-3 (window ±0.004 d) to 5.7e-3 (±0.011 d).

Ephemeris: the manuscript's four events (S6×2 + S98 E0/E1) give P = 11.132739±0.000015 d and T0 = BTJD 4000.8655±0.0016, matching the manuscript's values; its error bars are about 2× larger than mine (conservative).

Combined (Fisher, per-candidate): 2.5–3.9σ depending on method. Recommendation: keep B as a low-significance candidate, with the wording fixes F1–F6.

REPORT.md was not written: the harness blocked writing report .md files, so the full report content is in this result. Scripts, outputs, logs and the diagnostic figure are in the audit subdirectory. No packages were installed.

## Numeric rows

| quantity | manuscript | recomputed | status | method / note |
|---|---|---|---|---|
| S98 clean transit depths (E0, E1, E4 partial) | 1103±161, 1267±161, 927±161 ppm | 1085±138, 1189±154, 951±243 ppm (QLP DET, trapezoid, red-noise errors) | MINOR_DIFF | b01_s98_events.py fixed-centre trapezoid + local linear baseline Values agree within 1σ; the error on the 'partial' E4 is understated (55% of cadences flagged 2^30) |
| S98 scattered-light transit | one in a scattered-light interval | E2 4023.134: 1212±200 ppm (2^29 half-orbit, background 1.9×); E−1 and E3 also in S98 but in 2048 straylight, unusable | MATCH | b01 Two further predicted S98 epochs are not mentioned |
| S98-only FAP (segment-scrambled) | 0.35 (300 trials) | 0.37 (SYS_RM, SDE, 3000); 1.6e-3 (QLP DET, SDE, 5000); 2e-4 (QLP DET S/N); ~1.0 (own TESScut) | MINOR_DIFF | b03_s98_null.py cyclic shifts within segments, B removed, BLS 1–20 d Reproduced only for SYS_RM+SDE; strongly method-dependent |
| Link-test peak box S/N (S6+S32) | 10.0 | 9.75 at P=11.132767 (pipeline); QLP 8.48, SAP 8.95, own 7.33 | MINOR_DIFF | b04_link.py flank-baseline box, empirical-β errors Estimator differences |
| Scan window | ±0.004 d around S98 value | S98-only σ_P = 0.0015 d (3 events) / 0.0038 d (2 full); ±3σ = ±0.0045–0.011 d; S98 BLS P varies 11.1316–11.1416 | DISCREPANT | b01/b08 timing fits with injection errors Window not derived from S98 error; results robust in my null |
| Link FAP (S6+S32) | 2e-3 (19/10,000) | <2e-5 formal (0/50,000 kept and removed; 36,471 distinct configs; tail fit ~5e-6); own phot 4e-5 | MINOR_DIFF | b04/b04b random P 9–13.5 d + epoch, scaled window Manuscript value is conservative (its null has S32 edge artifacts) |
| S6-only link FAP | 0 of 10,000 trials | 0/50,000 but only 1,418 distinct configurations → FAP ≲ 7e-4 | MINOR_DIFF | b04b distinct-configuration count Precision overstated |
| S6 predicted dips are two strongest 1.6 h dips | yes | Yes for PDCSAP/SAP/QLP/eleanor (all 4 detrendings); no for TGLC-aper (ranks 3,1) and own TESScut (ranks 2–4), where boxes at dumps and gap edges rank higher; TGLC-PSF corrupted | MINOR_DIFF | b05_s6_rank.py Rank p = 2.0e-3 (W1) to 5.7e-3 (W3) |
| S32 predicted events | (included in link) | 2175.10: 732±211 (S/N 3.5); 2186.23 in gap; 2197.36: 861±232 (S/N 3.7); S32-only FAP 0.02 | MATCH | b04/b05  |
| P_B | 11.13274±0.00003 d | 11.132739±0.000015 (same 4 events); 11.132751±0.000017 (7 events) | MATCH | b08_ephem.py injection-calibrated timing  |
| T0_B | 2461000.865±0.003 BJD_TDB | 2461000.8655±0.0016 (4 events); 2461000.8660±0.0018 (7) | MATCH | b08  |
| Overall significance | ~3σ | 2.5σ (manuscript p's, Fisher) to 3.9σ; link alone 2.9σ (ms) / ≳3.7σ (mine) | MATCH | Fisher combination of independent S98 and archival p Per-candidate; excludes search-wide trials |
| Recovery FAP attribution | two predicted transits in Sector 6 (recovery FAP 2e-3) | 2e-3 is the S6+S32 test; S6 alone ≲1e-3 | DISCREPANT | b04  |

## Findings

### F1 [MINOR] Abstract attributes the 2e-3 FAP to Sector 6 alone

The 2e-3 recovery FAP is for the combined S6+S32 test. S6 alone gives 0/10,000 in the manuscript and 0/50,000 here, which is resolution-limited to ≲1e-3.

*Evidence:* a3_signalB.json T1 (S6+S32 19/10000; S6 0/10000); out/b04_link_main.json; out/b04b_link_extra.json

*Location:* Abstract

*Proposed fix:* 'two predicted transits in Sector 6 (archival recovery in Sectors 6+32, FAP ≲2e-3)'

### F2 [MINOR] '0 of 10,000 trials' overstates precision

The S6-only null samples only about 1,400–1,750 distinct box configurations, so the FAP cannot be resolved below ~1e-3.

*Evidence:* out/b04b_link_extra.json: 1,418 distinct configs in 50,000 trials; earlier run: 1,745 distinct windows

*Location:* Sec. 3

*Proposed fix:* Report FAP ≲1e-3

### F3 [MINOR] S98-only FAP 0.35 is method-dependent

0.35 is reproduced only with SYS_RM_FLUX and the SDE statistic (0.37 in 3,000 trials). With QLP DET_FLUX the FAP is 1.6e-3 (SDE) or 2e-4 (S/N); own TESScut photometry gives ~1.

*Evidence:* out/b03_null_QLP_DET_bw075_shift_s1.json, out/b03_null_QLP_SYSRM_bw075_shift_s5.json, out/b03_null_OWN_bhalves_bw075_shift_s3.json

*Location:* Sec. 3

*Proposed fix:* Quote the range 2e-4 to ~1 and the method

### F4 [MINOR] Scan window ±0.004 d not justified by S98-only σ_P

S98-only σ_P is 0.0015 d (3 events) or 0.0038 d (2 full events). Method-dependent S98 period estimates span 11.1316–11.1416 d; my BLS peak (11.1416) lies outside the manuscript window. A ±0.011 d window leaves my FAP unchanged, but would raise the manuscript-null FAP to up to ~6e-3.

*Evidence:* out/b01_s98_events.json eph98; out/b08_ephem.json; out/b02_s98_bls.json; out/b04_link_windows.json

*Location:* Sec. 3

*Proposed fix:* State that the window is ≈3σ of the S98 period error, or widen it

### F5 [MINOR] 'Partial' S98 transit mislabelled; error understated

E4 is covered in time, but 55% of its in-transit cadences carry the QLP 2^30 flag, interleaved. Its depth error is ~240 ppm, not 161 ppm.

*Evidence:* logs/b01.log epoch 4: frac_b30 0.55, 951±243 ppm

*Location:* Sec. 3

*Proposed fix:* 'the last with half its cadences flagged', 927±240 ppm

### F6 [MINOR] 'Two strongest S6 dips' depends on the photometry

True for TESS-SPOC PDCSAP/SAP, QLP and eleanor-lite. In TGLC aperture and own TESScut photometry, boxes at dumps and gap edges (1481.27, 1478.15/1478.56) outrank one B dip. TGLC PSF photometry is corrupted near event 1.

*Evidence:* logs/b05.log; out/b05_s6_rank.json

*Location:* Sec. 3 / Fig. 1c

*Proposed fix:* 'in the TESS-SPOC light curve'

### F7 [NOTE] Unmentioned S98 epochs; the 'removed' null is formally anti-conservative

Epochs −1 and 3 fall in 2048-straylight intervals and are uninformative. Removing B's archival windows from the null is anti-conservative under H0, but the correct 'kept' null gives the same 0/50,000.

*Evidence:* out/b01_s98_events.json; out/b04_link_main.json

*Location:* Sec. 3

*Proposed fix:* 

### F8 [NOTE] Existence robust; ≈3σ supported

The link passes in 6 photometries and 3 windows; the other S98 peaks fail it (FAP 0.14–0.49); no pixel anomalies; the ephemeris is reproduced; the combined per-candidate significance is 2.5–3.9σ.

*Evidence:* out/b04b_link_extra.json; out/b08_ephem.json; b09_diagnostics.png

*Location:* Abstract / Sec. 5

*Proposed fix:* 

## Recommended wording

Abstract: "The 11.13 d signal (≈2 R⊕) rests on two clean Sector 98 transits (plus one flagged and one in scattered light) and two predicted transits recovered in Sector 6. The archival recovery (Sectors 6 and 32) has FAP ≲2×10⁻³, and the signal is significant at about 3σ overall."

Sec. 3: "Alone they do not establish it: depending on photometry and statistic, the S98-only FAP ranges from 2×10⁻⁴ to ≈1 (0.35 for SYS_RM_FLUX with SDE). Scanning periods within ±0.004 d (≈3σ of the S98-only period)… FAP = 2×10⁻³ (conservative). In Sector 6 the two predicted dips are the two strongest 1.6 h dips of the TESS-SPOC light curve outside the 3.18 d transits (FAP ≲10⁻³)." Keep B as a low-significance candidate.

## Files

- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/04_B_existence/common_b.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/04_B_existence/s98data.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/04_B_existence/archdata.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/04_B_existence/b00_ownphot.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/04_B_existence/b01_s98_events.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/04_B_existence/b02_s98_search.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/04_B_existence/b03_s98_null.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/04_B_existence/b04_link.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/04_B_existence/b04b_link_extra.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/04_B_existence/b05_s6_rank.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/04_B_existence/b08_ephem.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/04_B_existence/b09_plots.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/04_B_existence/b09_diagnostics.png
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/04_B_existence/out/
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/04_B_existence/logs/