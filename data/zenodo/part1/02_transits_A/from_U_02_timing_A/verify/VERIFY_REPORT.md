# VERIFY_REPORT: adversarial check of audit 02_timing_A (Signal A timing / TTV)

## What I did
I wrote independent code (v1_common.py, v2_real.py, v3_inj.py, v4_ephem.py) that uses a deliberately different method from the auditor's:
- **Inputs:** raw MAST files from 00_data/mast, sha256-checked against DATA_MANIFEST.json. S6 TESS-SPOC PDCSAP and S32 SPOC 2-min PDCSAP use lightkurve DEFAULT_BITMASK = 17087, printed from lightkurve 2.6.0. For S98 QLP SYS_RM_FLUX I keep only QUALITY == 0, which removes the 2^29/2^30 scattered-light bits.
- **Detrending:** none; instead a quadratic polynomial is fitted jointly with a fixed batman template (k = 0.0212, b = 0.35, rho = 2.51, u = 0.37/0.28, supersampled) in a +-0.3 d window.
- **B masking:** B cadences within T14_B/2 + 24 min are removed.
- **Scan:** mid-time scan over +-144 min, with the +-60-min sub-range, and **depth constrained >= 0**.
- **Injections:** 1500 per sector (500-ppm template at random times; windows containing a real A transit excluded), with the same selection as the real events (coverage >= 0.7, Delta chi2 >= 4, not at the scan edge).
- **Ephemeris:** WLS linear ephemeris with per-sector injection SDs, plus a Monte Carlo that draws per-event errors from the injection pools.

No auditor or package code was imported, and no derived product was used as input. Total compute was about 3 minutes.

## Results (manuscript vs recomputed)
| set / range | P (d) | sigma_P | T0 (BJD_TDB) | sigma_T0 | chi2/dof | MC sd_P |
|---|---|---|---|---|---|---|
| manuscript | 3.1827856 | 5.8e-6 | 2460998.7375 | 0.0023 | 12.6/18 | - |
| package's 20 events, R60 | 3.1827883 | 5.7e-6 | 2460998.73753 | 0.0027 | 12.1/18 | 5.8e-6 |
| package's 20 events, R144 | 3.1827906 | 6.2e-6 | 2460998.73714 | 0.0028 | 21.9/18 (p 0.24) | 6.2e-6 |
| auditor's 20 events, R60 | 3.1827852 | 5.6e-6 | 2460998.73592 | 0.0029 | 9.7/18 | 5.7e-6 |
| auditor's 20 events, R144 | 3.1827871 | 6.0e-6 | 2460998.73549 | 0.0029 | 19.6/18 (p 0.36) | 6.0e-6 |
| 19 events (drop S32 E-573 and S98 E-1), R60 / R144 | 3.1827850 / 3.1827852 | 5.7e-6 / 6.0e-6 | 2460998.7360 | 0.0029 | 9.5/17 / 8.5/17 | 5.8-6.0e-6 |

Injection SDs (R60, selected) are 12.9 / 16.4 / 12.0 min for S6 / S32 / S98. The robust sigmas are 10.4 / 10.4 / 7.8 min. For R144, the S32 SD rises to 23.5 min (kurtosis 12, with 2.9% of recoveries beyond 60 min). White (Delta chi2 = 1) errors are too small by a factor of 1.3 / 4.3 / 2.6.

## Verdicts on the auditor's claims
- **F1 (chi2 depends on the search range via S32 E-573): CONFIRMED, MINOR.** My independent fit moves E-573 from -8.75 min (R60) to -81.25 min (R144), as the auditor found. The size of the effect is smaller in my code: chi2 = 19.6-21.9/18, p = 0.24-0.36, versus the auditor's 26.6/18.
- **F2 (sigma_T0 = 0.0023 is about 20% too small): CONFIRMED, MINOR.** My S98 SD is 12.0 min. The package's 10.1 min sits at about the 3rd percentile of 150-injection subsets, and my sigma_T0 is 0.0027-0.0029 d. Caveat: the package's own pipeline could legitimately have a smaller scatter, so this is a rounding-level issue. Quoting 0.003 d is safe.
- **F3 (S98 E-1 has no post-egress baseline): PARTIAL, MINOR.**
  - The coverage problem is real. Clean data stop at +1.02 h (egress is at +0.94 h), with two isolated points at +1.13 and +1.24 h and then a gap to +6.4 h. That gives about 0.12 T14 of post-egress data.
  - The auditor's "my fit prefers an anti-transit" (S/N -4.6 at -75.7 min) is an artefact of its method: wotan chunk detrending plus a depth allowed to go negative. My fit gives +18.75 min, depth 700 ppm, Delta chi2 = 19 over the no-transit model (the package got +19.5 min). Even my sign-free fit lands at +18.75 min.
  - **The proposed fix is incomplete.** S32 E-573, which the auditor keeps, has an even shorter post-egress baseline: 0.035 T14, about 4 min. The B transit mask starts about 6 min after A's egress. A "baseline on both sides" rule therefore removes E-573 as well, and the count becomes 19, not "stays 20".
  - Removing both one-sided events gives a clean result: P = 3.1827850 +- 0.0000057 d, chi2 = 9.5/17 at R60 and 8.5/17 at R144. The search-range dependence of F1 then disappears.
- **Numeric row T0 (MINOR_DIFF): CONFIRMED.** The value itself is reproduced: with the package's event set I get 2460998.73753, matching the manuscript's 2460998.7375. The auditor's 0.7361 comes from its different event set. The error is about 20% low.
- **Numeric row "number of transits = 20" (MINOR_DIFF): PARTIAL.** The difference between the two event sets is real. My fit supports the package's inclusion of E-1 on the grounds of signal, but both E-1 and E-573 should be treated alike.
- **Numeric row chi2 (MINOR_DIFF): CONFIRMED.** With the package's event set I reproduce 12.1/18 against the manuscript's 12.6.
- **Numeric row per-sector injection sigma (MINOR_DIFF): CONFIRMED.** The per-sector SDs differ between pipelines, but sigma_P does not.

## Spot-checks of the auditor's MATCH conclusions
1. **sigma_P = 5.8e-6 is not inflated: CONFIRMED.** My independent Monte Carlo gives sd(P) = 5.7-6.2e-6. The observed chi2 falls at the 17th-74th percentile of the MC distribution, depending on the set and search range.
2. **chi2 = 12.6 is not suspiciously low: CONFIRMED.** For the package's set, 29% of MC draws fall below it.
3. **No TTV evidence: CONFIRMED (plausible).** The linear-ephemeris chi2 is acceptable in all variants (p >= 0.24). A back-of-envelope sensitivity estimate, Delta chi2 of about N A^2 / (2 sigma^2) = 6 with N = 20 and sigma of about 13 min, gives A of about 10 min, consistent with the auditor's limit of about 10-13 min. I did not reproduce the periodogram.

There is no circularity: the auditor's scripts do not import package code. The only place package numbers are read is the explicit what-if in t8, which uses the package's SDs.

## Missed issues
- **MINOR:** S32 E-573 is one-sided in the same way as S98 E-1, because the B transit mask truncates its post-egress baseline. It is the event that drives F1. The recommended wording or selection rule should cover both events, and the count may become 19.
- **NOTE:** the auditor's scan leaves the sign of the depth free. On low-S/N or one-sided windows this can lock onto anti-transits, as it did for E-1, and it drove the auditor's exclusion of E-1. A positive-depth constraint is the more defensible choice for timing a known transit.
- **NOTE:** the size of the white-error chi2 depends on the method (my white chi2 is 33-56/18 against the auditor's 120). F4's qualitative point stands (white errors are too small; the old chi2 of 120 is not caused by grid quantisation), but the auditor's 120.4 figure is specific to its own pipeline.

## Files
v1_common.py, v2_real.py (-> v2_real_events.csv, v2_real.log), v3_inj.py (-> v3_inj_S{6,32,98}.csv), v4_ephem.py (-> v4_results.json, v4.log), v0_inspect.py (raw S98 E-1 quality runs).
