# Verification of audit 01_ephemeris_A (Signal A ephemeris)

Verifier scripts: v01_qlp_bits.py, v02_s110.py, v03_grid.py (own code; inputs only 00_data raw FITS and HEASARC snapshots in ../docs). Outputs: v01_qlp_bits.json, v02_s110.json, v03_grid.json, v03_sector_shifts.json, v03_grid.log.

## Independent re-derivations
- v03 (trapezoid template, T14 = 1.87 h, per-window linear baseline profiled, B masked +-1 h, QUALITY == 0, coverage >= 0.5 with baseline on both sides, 30-min cadence supersampled): 23 windows = 6 (S6) + 8 (S32) + 9 (S98). Grid best fit P = 3.1827856 d (the manuscript value to 1e-6 d grid step), T0 = BTJD 3998.7368 (-1.0 min from manuscript), depth 502 ppm, reduced chi2 0.98. Per-sector shifts at the manuscript ephemeris: S6 +1.5 (-2..+4) min, S32 -3.5 (-4.5..0) min, S98 -1.0 (-1.5..0) min (white-noise dchi2 = 1 intervals, so they are underestimates).
- v01: QLP S98 QUALITY bits: 2^11 (2048) 5757 cadences, 2^29 10676, 2^30 5474, 2^5 24, 2^2 4; QUALITY == 0 keeps 12067 of 24162. 2^29 is not defined in the QLP data-product description (sec. 2.3.3 lists bits 1,3,4,5,6,8,11,12,13) nor in the MIT QLP update notes (bits 3, 6, 12, and 2^30; "Bits 16-30 are unused").
- v02: S110 (HEASARC: Oct 31 - Nov 26 2026, Camera 1 boresight 75.88, -2.69): A transits (UTC, geocentric) Nov 1 03:45 ... Nov 23 10:26 (8), 9th Nov 26 14:49 (+-3.5 min); B Nov 1 11:24, Nov 12 14:34, Nov 23 17:45 (+-4.5 min). Identical to the auditor.
- Package reconstruction (diagnosis only): the manuscript P, T0, chi2 (3.1827856, 3998.73754, 12.63/18) are exactly the package a1 'e_dt_inj_min' fit, which uses ONE injection error per sector (S6 13.8, S32 24.6, S98 10.1 min) for every transit; its 20 = S6 5 (E-795 dropped because its window touches a flagged cadence) + S32 6 (E-567 S/N 1.64 < 2, E-566 coverage 0.64 < 0.7) + S98 9. The '23' comes from a different window definition (6+8+9, coverage >= 0.5), as the package's own number audit admits.

## Verdicts
| item | auditor | verifier |
|---|---|---|
| F01 20-transit rule unstated | MINOR | CONFIRMED, MINOR |
| F02 S98 mask relies on undocumented bit 2^29 | MINOR | CONFIRMED, MINOR |
| F03 S110 wording | MINOR | CONFIRMED, MINOR (counts correct) |
| T0 MINOR_DIFF | MINOR_DIFF | PARTIAL: real but method-dependent spread (-1.0 min in my fit; package white/perm variants +1.9 min); within 1 sigma of the quoted 3.3 min; NOTE |
| chi2 MINOR_DIFF | MINOR_DIFF | PARTIAL: 12.6/18 exactly reproduced from the package; the auditor's 11-12/18 agree; difference is error-model choice; NOTE |
| n=20 MINOR_DIFF | MINOR_DIFF | CONFIRMED as MINOR (same issue as F01) |
| P_A MATCH | MATCH | CONFIRMED (independent grid gives the same P) |
| Cycle count MATCH | MATCH | CONFIRMED (no circularity; auditor code reads only 00_data and its own template_params.npy) |
| 23 transits MATCH | MATCH | CONFIRMED (independent window rule gives 6+8+9) |

## Missed issues
1. NOTE: the adopted timing errors are a single per-sector injection scatter applied to every transit, not per-event errors; the method sentence ('errors calibrated by injecting the template at random times') is accurate but the chi2 = 12.6 depends on this choice.
2. NOTE: S32 E-566 (in the 23-transit joint fit) lies 1.93 h before a predicted B transit; the propagated B timing error back 162 cycles is ~8 min (sqrt(0.003^2 + (162*3e-5)^2) d), comparable to the ~10 min gap between A egress and B ingress, so B contamination of that event's post-egress baseline cannot be excluded.
3. NOTE: the abstract's '23 transits' counts usable transits; at least 27 A transits fall on S98/S6/S32 data including 4 fully covered in QLP-flagged S98 intervals.
