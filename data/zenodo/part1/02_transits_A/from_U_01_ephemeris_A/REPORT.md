# Audit 01_ephemeris_A — REPORT

> Rendered verbatim from the agent's structured result (`_structured/audit__01_ephemeris_A.json`); the sub-agent harness blocked agents from writing .md report files. Scripts, outputs and logs are in this folder.


**Verdict:** PASS_WITH_CHANGES

## Summary

I reproduced the Signal A ephemeris independently, with my own code on the fresh 00_data files. It holds within the quoted errors. The FAIL criterion is not triggered: the adopted P reproduces the transit times of all three sectors within errors.

NOTE ON DELIVERABLES: the harness blocked writing REPORT.md ("Subagents should return findings as text, not write report files"). The full report is in this summary, the findings and numeric_rows. All scripts, JSON/CSV outputs, figures, logs and docs are in /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/01_ephemeris_A. The compiled numbers are in out_summary.json.

DATA. I used TESS-SPOC S6 30-min PDCSAP, SPOC S32 2-min PDCSAP (plus TESS-SPOC 10-min as a check), QLP S98 200-s SYS_RM_FLUX (SAP_FLUX and DET_FLUX as checks) and TESScut time columns, all with QUALITY==0.
- Per the MIT QLP news page, DET_FLUX replaced KSPSAP_FLUX as the detrended default from S56, and SYS_RM_FLUX is the systematics-removed, not-detrended flux.
- HEASARC S110 pages were downloaded to docs/, with sha256 and UTC time in docs/download_log.txt.
- No packages were installed.

1. TIME SYSTEMS (s01, s01b)
- All products have TIMESYS=TDB, TIMEREF=SOLARSYSTEM and BJDREFI=2457000. SPOC and TESScut have TIMEPIXR=0.5, i.e. mid-exposure (SDPDD Rev F).
- My own astropy Roemer delay toward the target matches SPOC and TESS-SPOC TIMECORR to -0.2..+1.3 s.
- QLP S98 TIME equals SPOC spacecraft time plus the target delay to +1.1..+2.5 s.
- TESScut uses the CCD reference-pixel correction and is off by at most 20 s (S6 +2.7..+19.3 s; S98 -9.9..+9.1 s). This is far below the timing errors.
- Known timing issues are irrelevant here:
  - S6 DRN clock offset of about 2 s.
  - The 2-s SPOC off-by-one error, fixed from Sector 20 (heasarc.gsfc.nasa.gov/docs/tess/sector20_summary.html).
  - QLP barycentric bugs in S74-79 and orbits 210-220 (archive.stsci.edu/hlsp/qlp). S98 is orbits 205-208, so it is not affected.

2. MASKS AND DETRENDING (s02, s04)
- The S98 QLP scattered-light intervals are flagged by an undocumented bit 2^29: 10,676 cadences covering the first ~6 d of each orbit. Documented bit 2048 covers 5,757 cadences; the 2^30 QLP flag covers 5,474.
- A biweight detrend (0.5/0.75/1.0 d, A and B masked) fails in the short S98 segments left after masking; it collapses to a constant.
- So the primary fits use raw flux with per-transit linear baselines fitted jointly with the transit. Biweight detrending is used only as a sensitivity check.

3. EPOCH CENSUS (s03)
- S6: 6 epochs with data. E-795 has coverage 0.80; all six lie 0.4-3.5 h from a momentum dump.
- S32: 8 epochs. E-566 abuts a predicted B transit (B ingress about 10 min after A egress).
- S98: 18 predicted epochs, of which 9 are clean, 4 are fully covered but only in scattered-light intervals, 2 are partial in scattered light, E8 has coverage 0.09 and 2 have no usable data.
- '23' = 6+8+9 transits with coverage >= 0.5 (in fact >= 0.7) after masking. Reproduced exactly.
- '20' is reproduced only with the package's unstated rule (coverage >= 0.7, S/N >= 2, not at the scan edge). My implementation of that rule gives 21 transits (6/6/9); 22 without the S/N cut.

4. ROUTE 1: PER-TRANSIT TIMING (s06-s08)
- Injection-recovery used 400 injections per sector. Errors in minutes, as 16-84% half-width / MAD / std: S6 12.2 / 11.5 / 19.4; S32 14.5 / 9.1 / 29.8; S98 9.8 / 7.4 / 16.5.
- The error distributions are heavy-tailed: 6-15% of injections are off by more than 30 min. Formal errors underestimate the true scatter by up to 6x.
- Linear ephemerides, with T0 in BTJD at the manuscript reference epoch:
  - T1 (21 transits): P = 3.1827826 +- 0.0000048, T0 = 3998.7338 +- 0.0022, chi2 41.7/19. The chi2 is driven by S32 E-573 at -79 min, whose chi2 surface is bimodal.
  - T1 without E-573 (20 transits): P = 3.1827789 +- 0.0000049, T0 = 3998.7344, chi2 11.9/18, rms 8.2 min.
  - Robust (Huber) fit: P = 3.1827813, T0 = 3998.7352.
  - Manuscript-like 20 with std errors: P = 3.1827825 +- 0.0000085, T0 = 3998.7339 +- 0.0038, chi2 11.1/18.

5. ROUTE 2: BLS AND GLOBAL FIT (s05, s09, s10)
- Combined 3-sector BLS: best P = 3.182784, dchi2 = 257. Every other peak within +-0.02 d is at least 157 lower; the +-1 and +-2 S98<->S32 aliases are at least 180 lower.
- Each sector alone detects the signal near 3.18 d: S6 dchi2 104, S32 98, S98 79.
- Global batman least-squares fit (23 windows): P = 3.1827858 +- 0.0000020, T0 = 3998.73666 +- 0.00102, corr(T0,P) = 0.79.
- emcee (32 walkers x 20,000 steps, red-noise beta 1.01/1.42/1.20; autocorrelation time 171-190 steps for P and T0): P = 3.1827834 +0.0000036/-0.0000040, T0 = 3998.73559 +0.00169/-0.00209, T0 at the decorrelated epoch = 2789.9112 +- 0.0013.
- Out of scope: T14 = 1.90 +0.10/-0.08 h, consistent with the manuscript's 1.87.

6. CYCLE COUNT (s11)
- S98-only: P = 3.18276 +- 0.00024. Extrapolated 576 cycles to S32, the timing uncertainty is 0.135 d, so P/sigma = 23.5. The prediction lands 0.10 sigma from the adopted ephemeris.
- S32 and S98 together predict the S6 transits to +-2.2 min; the S6 data lie at -0.96 sigma.
- Alias spacing: 0.00557 d (S98-S32), 0.0145 d (S32-S6), 0.00402 d (S98-S6).
- Alternative S98-S32 cycle counts (k = +-1 to +-3) fit worse by dchi2 265-349 in the global fit.
- Per-sector shift at the manuscript P and T0: S6 +1.0 (+3.3/-2.5) min, S32 -1.5 (+1.0/-2.8) min, S98 -0.75 +- 1.25 min.

7. HISTORICAL 3.182945 d (s12)
- This is the frozen S98-only BLS peak on a coarse grid (step ~2e-4 d at 3.18 d). My S98-only BLS peaks at 3.18271-3.18273.
- The difference from the adopted P is 1.6e-4 d, which is 0.37-0.87 sigma of the S98-only period uncertainty (2.35e-4 template fit, 1.84e-4 bootstrap, 4.3e-4 per-transit).

8. SECTOR 110 (s13)
- HEASARC gives 2026 Oct 31 - Nov 26.
- tess-point 0.9.5.post1 (pointing matches HEASARC) puts the target on Camera 1, CCD 2, column 703, row 1577, i.e. on silicon.
- A transits (UTC, sigma about 3.4 min): Nov 1 03:45 ... Nov 23 10:26 gives 8 transits, plus a 9th on Nov 26 14:49. The end time on Nov 26 is still TBD, so 'eight or nine' is correct.
- B transits (sigma about 4.5 min): Nov 1 11:24, Nov 12 14:35, Nov 23 17:45. 'Three' is correct.
- These counts ignore data gaps and scattered light.

FINDINGS (no FATAL or MAJOR)
- Minor wording and reproducibility items: the unstated 20-transit selection rule, the undefined S98 scattered-light mask, and the S110 wording.
- Notes: T0 depends on the method at the 1-5 min level; the per-transit errors are heavy-tailed; the cycle count is decisive; the historical period is explained; time systems are verified.

## Numeric rows

| quantity | manuscript | recomputed | status | method / note |
|---|---|---|---|---|
| P_A (d) | 3.1827856 +- 0.0000058 | Global LSQ 3.1827858 +- 0.0000020; emcee 3.1827834 +0.0000036/-0.0000040; per-transit T1 (21) 3.1827826 +- 0.0000048; T1 minus the S32 E-573 outlier (20) 3.1827789 +- 0.0000049; Huber 3.1827813; manuscript-like 20 with std errors 3.1827825 +- 0.0000085; BLS 3.182784 | MATCH | Own batman template timing with injection errors; global LSQ and emcee with profiled per-window baselines; astropy BLS Every variant is within 0.9 sigma (out_s05, out_s08, out_s09_emcee_long, out_s10) |
| sigma_P (d) | 0.0000058 | 0.0000048 (injection 16-84% errors) / 0.0000081-0.0000085 (injection std) / 0.0000039 (emcee, beta-inflated) / 0.0000020 (LSQ formal) | MATCH | Weighted LSQ and MCMC The manuscript value lies inside the range of my error models |
| T0_A (BJD_TDB) | 2460998.7375 +- 0.0023 | LSQ 2460998.73666 +- 0.00102; emcee 2460998.73559 +0.00169/-0.00209; per-transit 2460998.7338-.7355 (+- 0.0022-0.0038) | MINOR_DIFF | Same as P Recomputed T0 is 1.2-5.3 min earlier (at most 1.5 sigma). Driven by ambiguous S98 events E-1 (chi2 minima at -3 and +19 min) and E10 |
| sigma_T0 (d) | 0.0023 | 0.0022 (16-84%) / 0.0037-0.0038 (std) / 0.0019 (emcee) / 0.0010 (LSQ) | MATCH | As above corr(T0,P) = 0.74-0.79 at the quoted epoch; decorrelated epoch about -290 to -400 |
| Number of transits in the timing fit | 20 | 21 (6/6/9) under the analogous rule: coverage >= 0.7, no B overlap, S/N >= 2, not at grid edge; 22 without the S/N cut; 20 = 5/6/9 only with the package rule | MINOR_DIFF | Epoch census (out_s03) and per-transit fits (out_s06) The rule is not stated in the manuscript. Differences: S6 E-795 is clean in my fit (S/N 4.7); S32 E-566 abuts a B transit |
| Timing-fit chi2 / dof | 12.6 / 18 | 11.9/18 (20 transits without S32 E-573, 16-84% errors); 11.1/18 (manuscript-like 20, std errors); 41.7/19 (21 transits, 16-84% errors, one bimodal outlier) | MINOR_DIFF | Weighted LSQ Depends on the error statistic and on how one bimodal S32 event is handled |
| rms O-C (min) | not quoted (package: 9.0) | 8.2 (20 transits without the outlier); 19.0 (21 transits) | MATCH | Weighted LSQ residuals  |
| Transits in the joint fit | 23 | 23 = 6 (S6) + 8 (S32) + 9 (S98) with coverage >= 0.5 after the quality mask | MATCH | Epoch census 4 more fully covered S98 transits sit in scattered-light intervals (E2, 3, 6, 11); 27 or more were observed in total |
| Cycle counts S6 -> S32 -> S98 | implied: E-791 -> E-571 is 220 cycles, E-571 -> E0 is 571 cycles | Unique. Alternatives k = +-1..3 are worse by dchi2 265-349 in the global fit; BLS aliases are at least 157-180 lower; the S98-only extrapolation to S32 has P/sigma = 23.5 | MATCH | Global template fit, BLS, single-sector ephemerides S32+S98 predict S6 to +-2.2 min; S6 lies at -0.96 sigma |
| Per-sector mean O-C at the manuscript ephemeris (min) | implied ~0 | S6 +1.0 (+3.3/-2.5); S32 -1.5 (+1.0/-2.8); S98 -0.75 +- 1.25 (sector-shift fit); per-transit weighted means +1.8 +- 5.0, +3.8 +- 6.5, -4.7 +- 3.3 | MATCH | Sector-level time-shift fit with beta-inflated sigma FAIL criterion not met |
| Frozen S98 BLS period vs adopted P | 3.182945 (candidate.json) vs 3.1827856 | Difference 1.6e-4 d = 0.37-0.87 sigma of the S98-only sigma_P (2.35e-4 / 1.84e-4 / 4.3e-4); my S98-only BLS peaks at 3.18271-3.18273 | MATCH | S98-only BLS following config/prejev_search.json; template fit; bootstrap Explained by single-sector precision and the coarse grid (~2e-4 d) |
| Time system | BJD_TDB | All products TDB/SOLARSYSTEM, BJDREFI 2457000, mid-exposure; target barycentric corrections validated to <= 2.5 s | MATCH | FITS headers and astropy Roemer delay TESScut CCD-reference offsets <= 20 s; S98 not affected by the QLP timing bugs |
| Sector 110 dates | 2026 October 31 to November 26 | Oct 31 to Nov 26 2026 (HEASARC); target on Camera 1, CCD 2 (col 703, row 1577) | MATCH | HEASARC pages; tess-point 0.9.5.post1  |
| Signal A transits in S110 | eight or nine | 8 up to Nov 26 00:00 UTC (Nov 1 03:45 ... Nov 23 10:26); 9th on Nov 26 14:49 (+-3.4 min) | MATCH | Manuscript ephemeris propagated to UTC at the geocentre Data gaps and scattered light not considered |
| Signal B transits in S110 | three | 3: Nov 1 11:24, Nov 12 14:35, Nov 23 17:45 UTC (+-4.5 min) | MATCH | Manuscript B ephemeris  |
| T14 (h) | 1.87 +- 0.05 | 1.90 +0.10/-0.08 (emcee, free LD and density) | MATCH | emcee Cross-check only |

## Findings

### F01 [MINOR] The rule behind the 20-transit timing set is not stated and depends on the implementation

The manuscript gives '20 transits' for the timing fit and '23' for the joint fit without saying how 3 transits were excluded. The package rule (coverage >= 0.7, S/N >= 2, not at the scan edge) gives 5/6/9. My implementation of the same rule gives 21 (6/6/9). The ephemeris is insensitive to the choice: all selections give P within 0.9 sigma.

*Evidence:* out_s03_epochs.csv; out_s06_per_transit_all_variants.csv: S6 E-795 has S/N 4.71 and O-C +6.9 min; S32 E-567 has S/N 1.48; S32 E-566 abuts B (dt_B = +1.93 h). out_s08_ephemeris_fits.json. Package: analysis/out/a1_per_transit.csv (diagnosis only).

*Location:* Section 2, first sentence

*Proposed fix:* State the rule, e.g. 'from the 20 transits with >=70% coverage and S/N >= 2'.

### F02 [MINOR] The S98 scattered-light mask is not defined and relies on an undocumented QLP bit

'Sector 98 scattered-light intervals removed' is in practice QUALITY != 0. The scattered-light part is bit 2^29, which marks the first ~6 d of each orbit and is not defined in any QLP document (the QLP-update page lists bits 3, 6, 12 and 31=2^30 only). This mask removes 4 fully covered A transits. Including them changes P by less than 0.3 sigma.

*Evidence:* out_s02_s98_quality_intervals.json: 2^29 flags 10,676 cadences in 3988.50-3995.16, 4003.23-4009.39, 4017.71-4023.49 and 4032.03-4038.13; bit 2048 flags 5,757. Global LSQ with scattered-light transits: P = 3.1827862 vs 3.1827858. Per-transit O-C of E2, E6, E11: +1.9, +2.3, +0.3 min.

*Location:* Section 1, last sentence

*Proposed fix:* '...with all QLP-flagged Sector 98 cadences (including the scattered-light first half of each orbit) removed'.

### F03 [MINOR] The S110 transit counts are calendar counts, and 'eight or nine and three transits' is ambiguous

The numbers are correct for times within the HEASARC dates. However, the sentence does not say which count belongs to which signal, and it ignores downlink gaps and Camera 1 scattered light. On the same camera in S98, only 9 of 18 A epochs survived the quality mask.

*Evidence:* out_s13_s110.json: A has 8 transits up to Nov 26 00:00 UTC plus Nov 26 14:49; B at Nov 1 11:24, Nov 12 14:35, Nov 23 17:45. Target on Camera 1, CCD 2 (tess-point 0.9.5.post1). HEASARC sector110_summary.html gives 'Days of Science Data Collection: TBD'.

*Location:* Section 5 (Status)

*Proposed fix:* 'TESS Sector 110 (2026 October 31 to November 26) spans eight or nine predicted 3.18 d and three 11.13 d transits.'

### F04 [NOTE] T0 depends on the method at the 1-5 min level, and T0 and P are correlated

All my estimates put T0 earlier than the manuscript, by 1.2-5.3 min, i.e. at most 1.5 sigma. The spread comes from ambiguous per-transit times, especially S98 E-1, which is at the end of orbit 205a with 0.4 h of post-egress baseline. At the quoted epoch corr(T0,P) is about 0.75-0.79. The effect on predictions is negligible (S110 shift about 1-3 min).

*Evidence:* out_s05_globalfit_lsq.json: T0 = 3998.73666. out_s09_emcee_long.json: 3998.73559. out_s08: 3998.7338-3998.7355. out_s06c_chi2_surfaces.json: S98 E-1 chi2 minima at -3 min (0) and +19 min (0.35).

*Location:* Section 2

*Proposed fix:* No change required. Optionally quote T0 near the weighted-centre epoch, or give the covariance.

### F05 [NOTE] Per-transit timing errors are heavy-tailed; the chi2 value depends on the error statistic

Injection-recovery shows 6-15% of events off by more than 30 min; formal errors underestimate the scatter by up to 6x. With core errors, one bimodal S32 2-min event (E-573: -79 min; the 10-min FFI gives -9 min) raises chi2 from ~12 to ~42. Injection calibration is essential; the manuscript's approach is appropriate.

*Evidence:* out_s07_injection_summary_base.json; out_s06c_chi2_surfaces.json (S32 E-573 minima at -79, -90, -5 min with dchi2 0, 0.9, 1.2)

*Location:* Section 2

*Proposed fix:* None required.

### F06 [NOTE] The cycle count is unique and all sectors agree with the adopted ephemeris

Each sector alone detects P ~3.18 d. The S98-only ephemeris fixes the S32 cycle count at P/sigma = 23.5. S32+S98 predict S6 to +-2.2 min, and S6 agrees at -0.96 sigma. Alternative cycle counts are strongly excluded. Per-sector shifts are at most 0.8 sigma.

*Evidence:* out_s10_bls_aliases.json (combined dchi2 257; next peak 157 lower; S6/S32/S98 alone 104/98/79); out_s11_cycle_count.json (alternatives dchi2 >= 265; shifts +1.0/-1.5/-0.75 min)

*Location:* Section 2

*Proposed fix:* None.

### F07 [NOTE] The historical 3.182945 d is the coarse-grid S98-only BLS peak

The difference from the adopted P, 1.6e-4 d, is 0.37-0.87 sigma of the S98-only period uncertainty. The BLS grid step at 3.18 d is ~2e-4 d.

*Evidence:* out_s12_s98_frozen.json; candidates/tic4206066_p3.182945/candidate.json; config/prejev_search.json

*Location:* not in manuscript

*Proposed fix:* None.

### F08 [NOTE] Time systems are verified; no timing-product issue

All products are BJD_TDB at mid-exposure. Target barycentric corrections agree to 2.5 s or better. TESScut differs by up to 20 s (CCD reference pixel). S98 is outside the QLP barycentric-bug sectors.

*Evidence:* out_s01_time_systems.json; out_s01b_bary_validation.json; https://archive.stsci.edu/hlsp/qlp; https://heasarc.gsfc.nasa.gov/docs/tess/sector20_summary.html; S6 DRN section 3.4

*Location:* Section 2

*Proposed fix:* None.

### F09 [NOTE] Two edge cases among the 23 transits; biweight detrending fails in short S98 segments

S32 E-566 abuts a predicted B transit (B ingress ~10 min after A egress). My timing of this event fails, falling into the B mask. S98 E-1 is at an orbit end. After masking, S98 segments shorter than the biweight window leave ramps untreated. This matters for any detrended-flux depth or shape analysis (other audits).

*Evidence:* out_s03_epochs.csv; fig_s06_events_baseline.png; out_s06 variants detr_w0.5 and detr_w0.75 (S98 E-1 at +131 min)

*Location:* Sections 1-2

*Proposed fix:* None for the ephemeris; flagged for the shape and depth audits.

## Recommended wording

Section 2: '...give P = 3.1827856 +- 0.0000058 d and T0 = 2460998.7375 +- 0.0023 BJD_TDB from the 20 transits with >=70% coverage and S/N >= 2 (chi2 = 12.6 for 18 degrees of freedom).' Section 1: '...with all QLP-flagged Sector 98 cadences (including the scattered-light first half of each orbit) removed.' Section 5: 'TESS Sector 110 (2026 October 31 to November 26) spans eight or nine predicted 3.18 d and three 11.13 d transits.' The ephemeris numbers themselves need no change.

## Files

- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/01_ephemeris_A/out_summary.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/01_ephemeris_A/out_s08_ephemeris_fits.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/01_ephemeris_A/out_s09_emcee_long.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/01_ephemeris_A/out_s11_cycle_count.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/01_ephemeris_A/out_s13_s110.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/01_ephemeris_A/out_s06_per_transit_all_variants.csv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/01_ephemeris_A/out_s07_injection_summary_base.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/01_ephemeris_A/fig_s14_oc.png
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/01_ephemeris_A/fig_s10_bls_aliases.png
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/01_ephemeris_A/docs/download_log.txt