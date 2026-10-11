# Verification of audit 03_pixels — VERIFY_REPORT

> Rendered verbatim from the agent's structured result (`_structured/verify__03_pixels.json`); the sub-agent harness blocked agents from writing .md report files. Scripts, outputs and logs are in this folder.


## Overall assessment

The auditor's main results survive an independent re-derivation, and I found no circularity: the auditor's code reads only the 00_data raw files. My own difference-image and PRF code used a different pointing treatment (per-event POS_CORR regressors) and bootstrap errors. It reproduces S32 raw (+9.3, -8.4)\", S32 decorrelated (+5.1, +0.2)\", S98 clean (+5.4, -0.2)\", S6 raw (-8.6, -34.5)\", S6 dump-excised (+6.3, -1.2)\", and the S6 dump distances of 0.41-3.54 h.

All three MINOR findings are confirmed:
- F1: the manuscript's -8\" Dec is not supported. One S98 event (3995.555) has no post-transit baseline because a data gap starts at 3995.608, and S32 jitter pulls the Dec through event 2175.001.
- F2: the 11.13 d Dec should be 0 +- 4 by the package's own numbers.
- F3: Fig. 1d plots an unstable S6 3.18 d point that the text excludes.

Two weaknesses are in the auditor's supporting evidence, not in its conclusions:
- The '2.7 sigma' POS_CORR2 correlation is really 2.0 sigma by the auditor's own numbers and about 1.1 sigma unweighted. The '(+8.9, -12.2)\" predicted' offset comes from the centroid, not from POS_CORR.
- The injection-based errors behind the recommended '(+5 +- 2, 0 +- 2)\"' are smaller than the event-to-event scatter; about ±3\" is safer.

The manuscript's qualitative claims hold: the 3.18 d position is consistent with the target at about 2 sigma, capable stars lie at 43\" or more and 7 sigma or more away, and the 11.13 d signal is on target. The numbers in Sec. 4 and Fig. 1d need the changes the auditor proposed, with slightly more conservative error bars.

## Checks

### F1: 3.18 d position (+8 +- 5, -8 +- 7)" does not reproduce; the Dec -8" comes from pointing and one bad S98 event

**CONFIRMED**, corrected severity **MINOR**

I rewrote the method independently, with different choices from the auditor: POS_CORR1/2 enter as extra regressors in each per-event, per-pixel fit instead of a global pre-decorrelation; the target reference comes from a single-PRF fit to the median image; errors come from an event bootstrap and jackknife. The results match the auditor's. S32 2-min raw gives (+9.3, -8.4)" (auditor +9.1, -8.3). S32 with POS_CORR regressors gives (+5.1, +0.2)" (auditor +5.2, +0.4); the bootstrap SD falls from (11.9, 11.2)" to (3.1, 4.4)". In the raw S32 jackknife, dropping event 2175.001 moves Dec from -8.4" to -2.1"; that event has the largest POS_CORR2 excursion (+0.82 millipixels). S98 clean (2048 + QLP 2^29) gives (+5.4, -0.2)" from 9 events. The package's 10th event, at BTJD 3995.555, ends at 3995.594. The data gap starts at 3995.608 and runs to 3995.821, so that event has no post-transit baseline, and standard acceptance rules reject it. So the manuscript's -8" Dec is not supported. The conclusion stays qualitatively the same, and the manuscript's own error bars cover the corrected value, so MINOR stands. One part of the auditor's evidence is overstated (see missed issues): the '2.7 sigma' POS_CORR2 correlation and the 'predicted (+8.9, -12.2)"' offset.

*Evidence:* verify/v1_loc.py, v1_loc.json, v1_loc.log; S98 gap list from the TESScut S98 TIME column: (3995.608, 3995.821).

### F2: 11.13 d Dec '+1 +- 4' does not match the package's own value

**CONFIRMED**, corrected severity **MINOR**

a5_summary.json localisation.B.combined_S32tpf_S98clean gives offset (+1.41, +0.09)" and sigma (3.27, 4.43)", which rounds to (+1 +- 3, 0 +- 4)". MANUSCRIPT_NUMBER_AUDIT.csv row 52 shows the source of the '+1' in Dec: the internal audit's (+1.0, +1.4) value, combined with the package's error bars. Two sources were mixed. My independent S98 11.13 d value is (+1.2, -0.8)" with bootstrap SD (3.1, 2.9)". Every version is consistent with the target, so the error is cosmetic.

*Evidence:* rnaas_tic4206066_pubprep/analysis/out/a5_summary.json; MANUSCRIPT_NUMBER_AUDIT.csv line 52; verify/v1_loc.json (S98_ffi B)

### F3: Fig. 1d plots the S6 3.18 d point although the text calls S6 unusable

**CONFIRMED**, corrected severity **MINOR**

In a6_figure.py panel (d), mk includes 'S6_ffi', and the loop over loc[sig]['per_set'] draws it as a filled circle for the 3.18 d signal at (+8.19 +- 2.70, -2.39 +- 3.57)". The caption also says circles are Sector 6. The package's combined_all_independent chi2 is 14.55 for 2 dof. My own code gives raw S6 = (-8.6, -34.5)" (auditor -8.7, -34.7). Leave-one-out is extremely unstable: dropping any one of the last three transits gives about (+7..+10, -1..+1)", close to the package value, while dropping any of the first three keeps about (-8, -31)". The ±0.5 h dump-excised value is (+6.3, -1.2)" (auditor +6.3, -1.3). So the plotted S6 point is one arbitrary state of an unstable estimator, and the figure contradicts the text.

*Evidence:* a6_figure.py lines 59-64; verify/v3_s6.py, v3_s6.json

### Row: Figure 1d S98 3.18 d point (+7.0 +- 5.6, -9.3 +- 7.3)" (DISCREPANT)

**CONFIRMED**, corrected severity **MINOR**

My independent S98 value from 9 clean events is (+5.4, -0.2)". In the jackknife, E ranges +1.7 to +8.3" and N ranges -1.4 to +2.0". The package's 10-event set includes the 3995.555 event, which has no post-transit baseline. The auditor's rebuild of the package's recipe gives (+1.7, -8.1)"; the package itself gives (+7.0, -9.3)".

*Evidence:* verify/v1_loc.json; 03_pixels/diagnose_s98.log

### Row: Figure 1d S6 3.18 d point not reproducible (DISCREPANT)

**CONFIRMED**, corrected severity **MINOR**

The value depends on which transits are included and how the dumps are treated: (-8.6, -34.5)" raw, (+6.3, -1.2)" dump-excised, about (+8, 0)" when one late transit is dropped. The package's (+8.2, -2.4)" sits within this range but is not a robust result.

*Evidence:* verify/v3_s6.json

### Row: S6 usability 'unusable' vs 'unreliable' (MINOR_DIFF)

**CONFIRMED**, corrected severity **MINOR**

I independently confirmed that every transit in the data lies 0.41-3.54 h from a quality-bit-32 flag (0.76, 2.15, 3.54, 3.19, 1.80, 0.41 h; the 7th predicted transit falls in the gap). Because the raw result is unstable, 'not used/unreliable' is a fair description. 'Unusable' is also defensible, so this is a wording choice only.

*Evidence:* verify/v3_s6.json dist_h

### Row: S32 2-min 3.18 d point (MINOR_DIFF)

**CONFIRMED**, corrected severity **NOTE**

My raw value (+9.3, -8.4)" and POS_CORR-regressed value (+5.1, +0.2)" reproduce the auditor's. The package's plotted (+13.5 +- 17.2, -2.3 +- 21.6)" has large errors and is statistically compatible with both.

*Evidence:* verify/v1_loc.json

### Row: 3.18 d East offset (+8 +- 5 vs +5.5 +- 2.2) (MINOR_DIFF)

**PARTIAL**, corrected severity **MINOR**

The central value is confirmed: independently, S32(PC) is +5.1" and S98 is +5.4". The auditor's injection-based errors look optimistic. The auditor's own S98 jackknife sigma in East is 5.5" (diagnose_s98.log), against an injection SD of 3.3"; my S32(PC) bootstrap Dec SD is 4.4", against an injection value of 2.8". A realistic combined error is about (+-3, +-3)" rather than (+-2.2, +-2.4)". So the recommended '(+5 +- 2, 0 +- 2)"' may overstate precision; '(+5 +- 3, 0 +- 3)"' is safer. Using event-scatter errors also moves 'consistent at 2 sigma' further away from the borderline.

*Evidence:* diagnose_s98.log 'jackknife sigma=[5.5 3.1]'; verify/v1_loc.json boot_sd

### Row: 11.13 d North offset (MINOR_DIFF)

**CONFIRMED**, corrected severity **MINOR**

Same as F2.

*Evidence:* a5_summary.json

### MATCH spot-check: 'consistent at 2 sigma' (3.18 d)

**CONFIRMED**, corrected severity **NONE**

My independent positions, S32(PC) (+5.1, +0.2)" and S98 (+5.4, -0.2)", give an eastward offset of about 5" and no Dec offset. With realistic errors of 3-4" per sector, the combination lies at or below about 2 sigma. I found no circularity: common.py reads only 00_data files, and no earlier derived products are imported.

*Evidence:* verify/v1_loc.json; 03_pixels/common.py FILES/DATA

### MATCH spot-check: capable stars >= 43" and >= 7 sigma

**CONFIRMED**, corrected severity **NONE**

TIC 4206077 lies at (dE, dN) = (+4.5, +42.8)" (sep 43.06"). From the manuscript's own position (+8 +- 5, -8 +- 7)", the diagonal distance is sqrt(0.7^2 + 7.26^2) = 7.3 sigma; from the corrected position it is more than 15 sigma. For a Tmag 17.02 star against Tmag 9.68, the flux ratio is 1.2e-3, so a partial-aperture required depth of about 0.7-1.1 is plausible.

*Evidence:* 03_pixels/neighbours_vs_fitted_position.csv

### MATCH spot-check: S6 '0.4-3.5 h' from desaturation

**CONFIRMED**, corrected severity **NONE**

I recomputed the distances from the TESScut S6 QUALITY bit 32 with the manuscript ephemeris: 0.41-3.54 h. They also match the package's own transit_to_dump_h values.

*Evidence:* verify/v3_s6.json

## Issues the auditor missed

### [MINOR] The auditor's pointing evidence is overstated: POS_CORR2 is not a 2.7 sigma transit correlation, and the (+8.9, -12.2)" 'prediction' is not independent

The auditor's summary says POS_CORR2 is displaced by +0.37 +- 0.16 millipixels, '2.7 sigma vs random epochs'. The auditor's own pointing.json gives z_vs_random = 2.0 for POS_CORR2, because 0.368/0.135 ignores the random-epoch mean. My unweighted in-transit POS_CORR2 shift is +0.16 millipixels against a random mean of 0.055 and SD 0.097, only about 1.1 sigma, and it is dominated by event 2175.001 (+0.82 millipixels). In addition, the 'implied (+8.9, -12.2)"' offset in pointing.py is computed from the flux-weighted centroid of the target pixels, not from POS_CORR. That makes it essentially another centroid-based measurement of the same in-transit offset, not an independent prediction from pointing. The decorrelated result itself is sound and reproduces independently. The better way to describe the cause is: S32 pointing jitter inflates the noise and pulls Dec through event 2175.001. It is not a significant transit-correlated pointing shift.

*Evidence:* 03_pixels/pointing.json S32_tpf.A.z_vs_random = [-1.38, 1.51, -0.49, 2.00]; pointing.py implied_px = -mean[:2]/depth (centroid columns); verify/v2_poscorr.json

### [MINOR] The recommended wording '(+5 +- 2, 0 +- 2)"' understates the uncertainty

The event-to-event scatter (jackknife or bootstrap) is larger than the injection SD: S98 East 5.5" against 3.3", and S32(PC) Dec 4.4" against 2.8". Quoting ±2" invites a claim of a 2.5 sigma eastward offset. About ±3" per axis is more defensible, and 'consistent with the target within ~2 sigma' still holds.

*Evidence:* 03_pixels/diagnose_s98.log; verify/v1_loc.json boot_sd

### [NOTE] S32 11.13 d rests on at most one clean event

Under standard acceptance rules (enough pre- and post-transit baseline, with a ±0.5 h mask around 3.18 d transits), only the 11.13 d event at 2197.361 survives in S32. The event at 2175.096 lies about 2.3 h from the 3.18 d transit at 2175.001, so its baseline is compromised. The S32 11.13 d point in Fig. 1d, (-7.9, +13.0)", is effectively a single-event or compromised measurement and should be de-emphasised or dropped.

*Evidence:* verify/v1_loc.json S32_tpf B tcs = [2197.361]

### [NOTE] VERIFY_REPORT.md not written

The operator-level instruction for this sub-agent forbids writing report .md files, so the verification report is this structured result. Scripts and outputs are in 03_pixels/verify/: v1_loc.py/.json/.log, v2_poscorr.py/.json and v3_s6.py/.json.

*Evidence:* 
