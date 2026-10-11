# Audit 5: signal B (11.13 d) depth, TIC 4206066

Audit ID: 05_B_depth. Date: 2026-09-25. Everything was written from scratch in this directory. The inputs are only the fresh public files in `../00_data` (sha256 in `../00_data/DATA_MANIFEST.json`), and I did not import any code or tables from earlier analyses. I read the earlier outputs only to diagnose discrepancies (`rnaas_tic4206066_pubprep/analysis/a1_timing_depths.py`, `out/a5_summary.json`).

## Bottom line

- **Radius and depth scale: reproduced.**
  - With the density prior I get Rp = 2.23 ± 0.12 R⊕ (manuscript 2.25 ± 0.12).
  - With the density free I get Rp ≈ 2.0 R⊕.
  - The light curve does not constrain b or ρ*. This is confirmed.
  - The overall depth is ≈950–1000 ppm.
- **Sector 32 depth (1003 ± 206 ppm): not robust.** In Sector 32 both B transits start only 0.4–0.8 h after the end of a transit of signal A, whose mid-times are 2.26 h and 1.93 h earlier. The manuscript's analysis drops every cadence within ±T14_A = ±1.87 h of an A mid-time, so it discards 26% and 46% of the B in-transit cadences. It then fits a free time shift to each truncated event. If A is subtracted or modelled and the whole B transit is used, S32 comes out at ≈710–770 ppm.
- **Sector 6 depth (838 ppm): depends on one cadence.** The mid-transit cadence of the second S6 event is the *only* cadence in the whole TESS-SPOC S6 light curve that carries SPOC bit 10 (value 512, "impulsive outlier removed before cotrending"). The lightkurve default mask removes it. Keeping it gives S6 ≈ 1130 ppm. Every FFI product that does not drop that cadence gives 1060–1200 ppm (QLP, TGLC, eleanor, own photometry).
- **The "1.6σ" figure and the "common depth not rejected (Δχ² = 3.1)" wording are not a robust summary.**
  - My consistent baseline (BASE) gives S98 − archival = **2.2σ** with injection-calibrated errors, or 1.9σ with β-scaled errors.
  - The Δχ² for a common depth is 5.1 (p = 0.08), or 3.9 (p = 0.14) with β-scaled errors.
  - Across 38 defensible configurations the excess ranges from 0.9σ to 3.3σ (median 2.2σ); 68% of them exceed 2σ.
  - The manuscript's 1.6σ is reproduced only with the manuscript's own combination of choices (my MSLIKE: 1.5σ).
  - A defensible statement is: "Sector 98 appears deeper than the archival sectors by ≈2σ (1–3σ depending on reduction choices); the evidence for a sector-dependent depth is marginal."
- **Dilution and flux scale do not explain the difference.** Contamination is below 1%. Signal A, used as a calibrator, gives S98/archival = 1.04 ± 0.12 for the manuscript's products.

Verdict: **PASS_WITH_CHANGES**. The Rp numbers are fine. The per-sector depths, the 1.6σ and the "common depth not rejected" wording must be revised.

## Data and method

- **Products (as in the manuscript):**
  - S6: TESS-SPOC PDCSAP, 30 min (sha256 ac7736ea…).
  - S32: SPOC PDCSAP, 2 min (764dad34…).
  - S98: QLP SYS_RM_FLUX, 200 s (a0341366…).
  - Alternatives: TESS-SPOC 10-min S32, SPOC SAP with CROWDSAP/FLFRCSAP correction, QLP KSPSAP/SAP for S6 and S32, QLP DET/SAP for S98, TGLC S6, eleanor S6, and my own TESScut aperture photometry with r = 1.5/2.0/2.5/3.0 px (background from star-free pixels at 6–10 px, because TESScut FLUX_BKG = 0).
- **Masks:**
  - SPOC products: lightkurve 2.6.0 DEFAULT_BITMASK 17087 (the manuscript's choice). Variant: 16575, which keeps bit 512.
  - S98 QLP: spacecraft bits plus 2048 (predicted straylight), 2^29 (set on the high-background intervals: median background 51k vs 21k when clean) and 2^30 (QLP internal bad-quality flag, https://tess.mit.edu/qlp/). This is the manuscript's "scattered-light intervals removed".
- **Signal A:** removed by subtracting a fixed template (k = 0.0212, b = 0.35, mid-depth 522 ppm) at P_A = 3.1827856, T0_A = BTJD 3998.7375. Variants: A modelled with a free depth, and A masked at ±0.6 or ±1.0 T14_A.
- **Detrending:** wotan biweight 0.75 d with A and B cadences excluded from the trend (wotan recommends masking known transits). Variants: biweight 0.5/1.0/1.5 d, cosine (robust), rspline, hspline, and none (local 1st/2nd-order polynomial per event).
- **Depth:** amplitude of a fixed-shape batman template, quadratic limb darkening u = 0.32/0.31 (Claret 2017 TESS: ATLAS 0.375/0.278, PHOENIX 0.316/0.313 and 0.299/0.333; `out/ldc_claret2017.json`). Supersampled to ≤1-min sub-exposures. It is a linear least-squares fit in ±0.3 d windows with per-event linear baselines, at the fixed manuscript ephemeris (P = 11.13274, T0 = BTJD 4000.865). "Depth" means the model mid-transit depth, the same quantity the manuscript uses.
- **Errors:**
  - White: robust out-of-transit scatter.
  - Injection-calibrated κ: 300–400 injections per sector per configuration (`s03_injections.py`), either as the whole transit train at a random phase or as single events. The pipeline is re-run for each injection, including detrending.
  - β factor: time-averaging (Winn et al. 2008).
- **Joint shape fits** (`s02_joint_fit.py`, emcee, 40 walkers): k, b, ρ*, u1, u2, T0, P over all 7 events.
  - (i) ρ* ~ N(2.51, 0.24) ρ☉.
  - (ii) ρ* free: log-uniform, with a uniform-prior variant.
  - Rp uses R* = 0.623 ± 0.019 R☉ and R☉/R⊕ = 109.076.

## Events used

| Sector | Epochs | Coverage | Note |
|---|---|---|---|
| S6 | −227, −226 | 1.0, 0.9 | −226 loses its mid cadence under the default mask (coverage 0.6) |
| S32 | −164, −162 | 1.0, 1.0 | A mid-times 2.26 h and 1.93 h before the B mid-time. With A masked at ±1.87 h, coverage drops to 0.74 and 0.54 |
| S98 | 0, 1, 4 | 1.0, 1.0, 0.5 | Epoch 4 loses half its cadences to QLP bit 2^30 (normal background) |
| S98 flagged | −1, 3 (bit 2048 + 2^29); 2 (2^29 only) | – | Three epochs fall in flagged scattered-light intervals, not one as the manuscript says |

Plot of every event: `plots/s01_events.png`.

## Headline results: manuscript vs recomputed

| Quantity | Manuscript | Recomputed | Status |
|---|---|---|---|
| S6 depth | 838 ± 114 | 866 ± 143 (BASE); 812–1204 across configs; 1129 ± 122 if the bit-512 cadence is kept | MINOR_DIFF (method-fragile) |
| S32 depth | 1003 ± 206 | 766 ± 164 (BASE, full transits); 966 ± 268 reproduced only with the A mask and per-event timing | DISCREPANT |
| S98 depth | 1099 ± 93 | 1136 ± 91 (BASE); 1104–1310 | MATCH (within 1σ) |
| S98 per-event | 1103, 1267, 927 ± 161 | 1112 ± 153, 1275 ± 153, 894 ± 224 (partial) | MATCH (the partial event's error is understated) |
| S98 − archival | 1.6σ | 2.2σ (inj) / 1.9σ (β); range 0.9–3.3σ; MSLIKE 1.5σ | DISCREPANT (fragile) |
| Common depth | 996 ± 68 | 1005 ± 70 (BASE, κ); 955 ± 77 (β) | MATCH |
| Δχ² common vs per-sector | 3.1 | 5.1 (p = 0.077) κ; 3.9 (p = 0.14) β; range 2.3–11.2 | MINOR_DIFF |
| Previous draft: 3.1σ (white errors) | – | white-error BASE 2.7σ; MSLIKE 1.75σ | superseded |
| Rp, density prior | 2.25 ± 0.12 | 2.23 +0.12/−0.12 (k = 0.0329 ± 0.0015, b = 0.880 ± 0.017, T14 = 1.62 h); 2.37 ± 0.12 with the bit-512 cadence kept | MATCH |
| Rp, density free | ≈2.0 | 2.02 +0.32/−0.14 (log-uniform ρ; 98th percentile 3.9, grazing tail); 1.95 +0.13/−0.10 (uniform ρ) | MATCH |
| b, ρ* unconstrained | yes | b = 0.61 +0.31/−0.41, ρ* = 10.7 +9.6/−9.2 ρ☉ (free) | MATCH |
| Abstract "≈2 R⊕" | – | on-target value 2.2–2.4 | MINOR (wording) |

**Sanity check with limb darkening.** The mid depth is ≈950 ppm. At b = 0.88, I(μ)/⟨I⟩ ≈ 0.88, so k ≈ √(950e-6/0.88) = 0.0329 and Rp = 2.23 R⊕. At b = 0.5, the ratio is ≈1.13, so k ≈ 0.029 and Rp ≈ 1.97 R⊕. The prior-vs-free difference is therefore purely a limb-darkening effect, and a plain √depth gives 2.1 R⊕.

## Sensitivity grid (`out/s04_grid.csv`, `plots/s07_depths_by_config.png`)

Errors are injection-scaled. κ comes from the configuration's own injections for BASE, MSLIKE and KEEP512, and from BASE otherwise.

| Config | S6 | S32 | S98 | Archival | nσ white | nσ κ | Common | Δχ² (p) |
|---|---|---|---|---|---|---|---|---|
| BASE | 866 ± 143 | 766 ± 164 | 1136 ± 91 | 823 ± 108 | 2.68 | 2.22 | 1005 ± 70 | 5.1 (0.077) |
| MSLIKE | 880 ± 137 | 966 ± 268 | 1141 ± 106 | 897 ± 122 | 1.75 | 1.51 | 1036 ± 80 | 2.3 (0.309) |
| MSLIKE_fixeph | 866 ± 143 | 830 ± 196 | 1136 ± 91 | 854 ± 115 | 2.27 | 1.92 | 1028 ± 72 | 3.7 (0.155) |
| KEEP512 | 1129 ± 122 | 766 ± 178 | 1136 ± 98 | 1013 ± 101 | 1.56 | 0.88 | 1076 ± 70 | 3.6 (0.165) |
| A_MODEL | 866 ± 143 | 724 ± 165 | 1136 ± 91 | 806 ± 108 | 2.86 | 2.34 | 999 ± 70 | 5.9 (0.052) |
| A_MASK06 | 866 ± 143 | 710 ± 166 | 1136 ± 91 | 800 ± 108 | 2.92 | 2.38 | 996 ± 70 | 6.2 (0.046) |
| A_MASK10 | 866 ± 143 | 830 ± 196 | 1136 ± 91 | 854 ± 115 | 2.27 | 1.92 | 1028 ± 72 | 3.7 (0.155) |
| BW050 | 830 ± 140 | 785 ± 163 | 1153 ± 90 | 811 ± 106 | 2.90 | 2.46 | 1009 ± 69 | 6.1 (0.048) |
| BW100 | 875 ± 145 | 781 ± 165 | 1162 ± 92 | 834 ± 109 | 2.77 | 2.30 | 1025 ± 70 | 5.5 (0.065) |
| BW150 | 873 ± 151 | 824 ± 167 | 1183 ± 94 | 851 ± 112 | 2.69 | 2.27 | 1045 ± 72 | 5.2 (0.075) |
| COSINE | 873 ± 146 | 780 ± 164 | 1170 ± 90 | 832 ± 109 | 2.89 | 2.40 | 1033 ± 69 | 5.9 (0.052) |
| RSPLINE | 875 ± 144 | 787 ± 164 | 1171 ± 94 | 837 ± 108 | 2.79 | 2.34 | 1027 ± 71 | 5.6 (0.060) |
| HSPLINE | 875 ± 141 | 805 ± 164 | 1161 ± 96 | 845 ± 107 | 2.60 | 2.19 | 1020 ± 72 | 4.9 (0.085) |
| LOCALPOLY2 | 918 ± 151 | 934 ± 176 | 1182 ± 100 | 925 ± 115 | 1.95 | 1.69 | 1071 ± 75 | 2.9 (0.238) |
| LOCALPOLY1 | 883 ± 141 | 853 ± 165 | 1171 ± 94 | 870 ± 107 | 2.47 | 2.11 | 1041 ± 71 | 4.5 (0.107) |
| POLY0 | 866 ± 143 | 766 ± 164 | 1137 ± 91 | 823 ± 108 | 2.69 | 2.22 | 1006 ± 70 | 5.2 (0.076) |
| DT_SECTOR | 867 ± 142 | 769 ± 164 | 1137 ± 91 | 825 ± 107 | 2.67 | 2.21 | 1006 ± 69 | 5.1 (0.078) |
| DT_EVENT | 880 ± 142 | 773 ± 164 | 1141 ± 91 | 834 ± 108 | 2.64 | 2.18 | 1013 ± 70 | 5.0 (0.083) |
| FREE_TPL | 849 ± 141 | 747 ± 161 | 1120 ± 90 | 805 ± 106 | 2.74 | 2.27 | 989 ± 68 | 5.4 (0.068) |
| NO_SUPERSAMPLE | 812 ± 136 | 766 ± 164 | 1135 ± 91 | 793 ± 105 | 2.90 | 2.47 | 988 ± 69 | 6.1 (0.047) |
| NO_PARTIAL | 909 ± 180 | 766 ± 164 | 1193 ± 101 | 831 ± 121 | 2.83 | 2.29 | 1044 ± 78 | 5.6 (0.061) |
| S98_KEEP_B30 | 866 ± 143 | 766 ± 164 | 1104 ± 87 | 823 ± 108 | 2.48 | 2.03 | 993 ± 68 | 4.3 (0.114) |
| S98_WITH_EP2 | 866 ± 143 | 766 ± 164 | 1210 ± 82 | 823 ± 108 | 3.47 | 2.85 | 1067 ± 65 | 8.3 (0.015) |
| S32_TESSSPOC10 | 866 ± 143 | 768 ± 178 | 1136 ± 91 | 828 ± 111 | 2.59 | 2.14 | 1012 ± 71 | 4.8 (0.091) |
| S32_SAPCROWD | 866 ± 143 | 822 ± 177 | 1136 ± 91 | 849 ± 111 | 2.37 | 2.00 | 1021 ± 70 | 4.0 (0.132) |
| S6_SAPCROWD | 905 ± 143 | 766 ± 164 | 1136 ± 91 | 845 ± 108 | 2.54 | 2.06 | 1015 ± 70 | 4.7 (0.097) |
| S6_QLP | 1059 ± 123 | 766 ± 164 | 1136 ± 91 | 954 ± 98 | 1.84 | 1.36 | 1052 ± 67 | 3.9 (0.143) |
| S6_TGLC | 1204 ± 170 | 766 ± 164 | 1136 ± 91 | 978 ± 118 | 1.68 | 1.06 | 1077 ± 72 | 4.6 (0.102) |
| S6_ELEANOR | 1077 ± 168 | 766 ± 164 | 1136 ± 91 | 918 ± 117 | 2.04 | 1.47 | 1054 ± 72 | 3.9 (0.142) |
| QLP_ALL | 1059 ± 123 | 602 ± 170 | 1136 ± 91 | 902 ± 100 | 2.42 | 1.73 | 1030 ± 67 | 7.7 (0.021) |
| QLP_ALL_SAP | 1053 ± 124 | 633 ± 170 | 1134 ± 92 | 907 ± 100 | 2.31 | 1.67 | 1031 ± 68 | 6.7 (0.034) |
| S98_DET | 866 ± 143 | 766 ± 164 | 1153 ± 90 | 823 ± 108 | 2.84 | 2.35 | 1018 ± 69 | 5.8 (0.056) |
| S98_SAP | 866 ± 143 | 766 ± 164 | 1134 ± 92 | 823 ± 108 | 2.65 | 2.20 | 1003 ± 70 | 5.0 (0.081) |
| S98_OWN_r2.0 | 866 ± 143 | 766 ± 164 | 1309 ± 99 | 823 ± 108 | 3.90 | 3.31 | 1085 ± 73 | 11.2 (0.004) |
| OWN_r1.5 | 1133 ± 168 | 558 ± 256 | 1308 ± 106 | 960 ± 140 | 2.72 | 1.98 | 1182 ± 84 | 7.4 (0.024) |
| OWN_r2.0 | 1091 ± 141 | 692 ± 189 | 1309 ± 99 | 948 ± 113 | 3.10 | 2.39 | 1152 ± 75 | 8.6 (0.014) |
| OWN_r2.5 | 1137 ± 134 | 705 ± 196 | 1273 ± 98 | 999 ± 111 | 2.50 | 1.85 | 1153 ± 73 | 6.7 (0.034) |
| OWN_r3.0 | 1130 ± 149 | 704 ± 199 | 1310 ± 101 | 977 ± 119 | 2.83 | 2.13 | 1172 ± 77 | 7.5 (0.024) |

What the grid shows:

- **Detrending, timing and template:** biweight windows, cosine, splines, per-sector or per-event timing, the density-prior vs free-density template, and supersampling all leave the excess at **2.2–2.5σ**. The local polynomial (no global detrending) gives 1.7–2.1σ.
- **A treatment:** dropping ±T14_A cadences lowers the excess to 1.9σ and raises S32 to 830 ppm, or 966 ppm with per-event timing.
- **Partial S98 event:** dropping it gives 2.3σ. Keeping its QLP-bad cadences gives 2.0σ.
- **Adding S98 epoch 2** (the 2^29-flagged event): 2.85σ.
- **Products:**
  - S6 FFI products (QLP, TGLC, eleanor) put S6 at 1060–1200 ppm, so the excess drops to 1.1–1.5σ.
  - QLP for S32 gives 600–630 ppm (QLP-classic S32 also gives about 20% shallower A transits).
  - My own photometry puts S98 at 1270–1310 ppm, with the excess at 1.9–2.4σ.

## Red noise (`out/s03_inj_*.json`, `out/s06_beta.json`)

| Sector | κ, sector (inj) | κ, event | β at T14 | sd of a full-event depth (inj) |
|---|---|---|---|---|
| S6 | 1.20 | 1.09–1.15 | 1.01 | 160–171 ppm |
| S32 (2 min) | 1.52–1.91 (seed-dependent) | 1.64 | 1.62 | 247–250 ppm |
| S98 | 0.97–1.17 | 1.07 | 1.34 | 154 ppm |

The injection bias is ≤ 25 ppm. The manuscript's per-event errors (S6 ~162, S98 161) agree with my full-event scatter. Its S98 partial event, however, is given the same error as the full events, although its error is about 40% larger (≈220 ppm).

## Dilution and flux scale (`out/s05_flux_scale.json`)

- **SPOC:** PDCSAP is crowding-corrected (CROWDSAP 0.9950 for S6 and 0.9970 for S32 2-min; FLFRCSAP 0.925/0.922). Applying the correction to SAP myself changes B by ≤ 60 ppm.
- **QLP S98:** built with TGLC photometry since Sector 94. The QLP data release notes (Petitpas et al. 2026, arXiv:2603.22236; copy in `docs/`) say it "forward models the FFIs with the effective point spread function (ePSF) to remove contamination from nearby stars". TGLC also rescales the aperture median to the Gaia-predicted flux (Han & Brandt 2023, AJ 165, 71). So QLP S98 is nominally undiluted, like PDCSAP.
- **Empirical check with signal A (fixed template):**
  - S6 PDC 538 ± 57, S32 2-min PDC 454 ± 49, S98 QLP 511 ± 45 ppm, so S98/archival = 1.04 ± 0.12.
  - My own S98 photometry gives 514–565 ppm.
  - No flux-scale offset is detectable at the 10–12% level. Contamination (<1%) cannot produce a 250–300 ppm (~30%) difference.
  - QLP-classic S32 (FIPHOT/TICA) gives 20–25% shallower A depths, 350–387 ppm. This is one reason not to mix QLP S32 with SPOC S32.
  - S6 FFI products give inconsistent A depths (630–970 ppm). The A transits in S6 fall near momentum dumps, as the manuscript notes.

## Findings

1. **MAJOR: the S32 depth is an artefact of the A mask.**
   - Manuscript: 1003 ± 206 ppm.
   - Both S32 B transits start 0.4–0.8 h after an A transit ends.
   - The A mask of ±1.87 h (`a1_timing_depths.py`: `~in_tr(t, other, 1.0)`) removes 26% and 46% of the in-transit B cadences, and the per-event free timing is then fitted to truncated events.
   - With A subtracted or modelled I get 710–770 ppm (BASE 766 ± 164).
2. **MAJOR: the "1.6σ" and "common depth not rejected (Δχ² = 3.1)" statements are not robust.**
   - The significance depends on three undisclosed choices: the A mask, removing one flagged S6 cadence, and the product used for S6.
   - BASE gives 2.2σ (inj) or 1.9σ (β), with a Δχ² of 5.1 or 3.9. The range is 0.9–3.3σ, and 68% of configurations exceed 2σ.
   - The conclusion that the depths are not significantly different at 3σ survives in nearly all configurations; the maximum of 3.3σ comes only from my own S98 photometry.
   - The single number 1.6σ is at the low end of the range, and "not rejected" hides the fact that the three per-sector depths scatter more than their errors in most reductions (p ≈ 0.02–0.15).
3. **MINOR: the S6 depth hinges on one cadence.**
   - The mid-transit cadence of S6 epoch −226 (BTJD 1484.870, PDCSAP −1960 ppm; QLP −1456, TGLC −1664) is the only bit-512 cadence in the whole S6 light curve.
   - Dropping it gives 866 ppm; keeping it gives 1129 ppm, and Rp_prior becomes 2.37 R⊕.
   - This should be disclosed, or the per-sector depths given as a range.
4. **MINOR: S98 event count and per-event errors.**
   - Three S98 B epochs, not one, fall in flagged scattered-light intervals (−1 and 3 in bit 2048 plus 2^29; 2 in 2^29).
   - The partial epoch 4 should carry ≈220 ppm, not 161 ppm.
5. **MINOR: the abstract "≈2 R⊕".** It is the free-density value. The on-target (density-prior) value is 2.23–2.25 R⊕, so "≈2.2 R⊕" is the fairer summary.
6. **NOTE: radius.** Rp = 2.23 ± 0.12 R⊕ (prior) and 2.0 R⊕ (free) are reproduced. The free-density radius depends on the prior (2.02 log-uniform, 1.95 uniform) and has a grazing tail (98th percentile 3.9 R⊕).
7. **NOTE: flux scale.** The products are on a common flux scale within ±12% (A calibrator), and dilution is <1%.
8. **NOTE: MCMC convergence.** The free-density chains have long autocorrelation times (τ ≈ 600 over 20,000 steps, about 22τ after burn-in). The medians are stable, but the tails are approximate.

## Recommended wording (Sec. 3)

"Per-sector depths depend on reduction choices: with signal A subtracted and injection-calibrated errors I find 870 ± 140, 770 ± 160 and 1140 ± 90 ppm for Sectors 6, 32 and 98 (Sector 6 rises to ≈1100 ppm if a mid-transit cadence flagged by PDC as an outlier is retained). Sector 98 is deeper than the archival mean by ≈2σ (1–3σ across detrending, product and masking choices), so a sector-dependent depth is possible but not established; the common depth is 1000 ± 70 ppm. The light curve does not constrain b or ρ*: Rp = 2.2 ± 0.1 R⊕ if the signal orbits the target and ≈2.0 R⊕ with the density free." Abstract: "(≈2.2 R⊕)".

## Files

- **Scripts:** `b_lib.py`, `b_run.py`, `s_configs.py`, `s00_ldc_claret2017.py`, `s01_diag_events.py`, `s02_joint_fit.py`, `s03_injections.py`, `s04_grid.py`, `s05_flux_scale.py`, `s06_beta.py`, `s07_summary.py`.
- **Outputs:** `out/`.
- **Logs:** `logs/`.
- **Plots:** `plots/`.
- **Docs:** `docs/arxiv_2603.22236_QLP_DRN004.pdf` (sha256 b0e41fab…, fetched 2026-09-25).

No packages were installed. Environment: audit_env (Python 3.12, batman, wotan, emcee, lightkurve 2.6.0).
