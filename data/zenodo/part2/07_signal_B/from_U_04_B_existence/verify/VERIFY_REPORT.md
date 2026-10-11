# Verification of audit 04_B_existence — VERIFY_REPORT

> Rendered verbatim from the agent's structured result (`_structured/verify__04_B_existence.json`); the sub-agent harness blocked agents from writing .md report files. Scripts, outputs and logs are in this folder.


## Overall assessment

The auditor's overall verdict stands (PASS_WITH_CHANGES, keep B as a low-significance ≈3σ candidate), but two of its supporting claims are wrong and need correcting.

My independent code on the raw public files (sha256 match the manifest) reproduces:
- the ephemeris: P = 11.132747 ± 0.000012 d, T0 = BTJD 4000.8658 ± 0.0015;
- the S98 depths;
- the understated error on the half-masked E4 transit (≈250 ppm, not 161);
- the method dependence of the S98-only FAP: SYS_RM gives 0.25 and DET_FLUX <5e-3;
- that the manuscript's S6+S32 FAP of 2e-3 is conservative, because its null is inflated by an edge artefact at the end of S32.

Corrections:
- Both the manuscript's null and the auditor's 'kept' null exclude the ephemerides that can hit B's own dips. For Sector 6 alone the correct chance-alignment FAP is ≈1–2e-3. The manuscript's '0 of 10,000 trials' is therefore overstated about 15x, and the auditor's replacement ('≲7e-4 / ≲1e-3') is also too strong. It should read '≈2×10⁻³ (rank test)'.
- F4's claim that the BLS peak lies outside the window rests on the auditor's coarse BLS grid; the missing justification for ±0.004 d is still a valid MINOR point.
- F6 is a NOTE: the only counter-examples are momentum-dump artefacts in TGLC and own photometry.
- The high S98-only FAP is caused by segment-edge detrending artefacts in SYS_RM_FLUX.
- The S6+S32 FAP is better stated as ~1e-5 to 1e-4 than as the tail-fit 5e-6.

No FATAL or MAJOR issues. Scripts, JSON outputs and logs are in 04_B_existence/verify/ (v01–v07). No packages were installed. VERIFY_REPORT.md was not written because the session's system rules forbid report .md files; this structured result carries the report content.

## Checks

### F1: Abstract attributes the 2e-3 recovery FAP to Sector 6 alone

**PARTIAL**, corrected severity **MINOR**

The wording point holds. In Sec. 3, 2e-3 is the combined S6+S32 test (a3_signalB.json T1: 19/10,000). The auditor's support for it does not: 'S6 alone ≲1e-3' is wrong. When the null is allowed to hit B's own dips, S6 alone gives about 1.1e-3 (verify full-prior null), 1.9e-3 (my rank test), 1.55e-3 (manuscript's own T2 rank test) and 2.0e-3 (auditor's rank test). So '2e-3' is numerically about right even for S6 alone, and the abstract is loose wording rather than a wrong number. The proposed fix ('archival recovery in Sectors 6+32, FAP ≲2e-3') is fine.

*Evidence:* verify/out/v02_link.json (S6 NEAR strip 1280/50,000 → FAP 1.14e-3); verify/out/v03_s6_rank.json p_rank 1.9e-3; rnaas_tic4206066_pubprep/analysis/out/a3_signalB.json T2 p_rank 1.55e-3

### F2: '0 of 10,000 trials' (S6) overstates precision

**CONFIRMED**, corrected severity **MINOR**

Confirmed, and the overstatement is larger than the auditor says. The manuscript's null removes B's archival windows, and the auditor's 'kept' null excludes |P0-Pc|<0.1 d. Neither can reproduce the observed chance alignment. I sampled that excluded strip directly (P0 ~ U(11.033, 11.233)); 2.56% of it reaches the observed S6 statistic, which gives an S6-alone FAP of about 1.1e-3 under the 9–13.5 d null prior. This agrees with the rank tests (1.5–2e-3). '0 of 10,000' implies <1e-4, so it overstates by about 15x. The auditor's 'FAP ≲7e-4 / ≲1e-3' has the inequality the wrong way round: 1/N_distinct is a floor, not a ceiling. The fix should read 'chance probability ≈2×10⁻³ (rank test)'.

*Evidence:* verify/out/v02_link.json runs[1] (EXCL 0/50,000, NEAR 1280/50,000, REM 0/50,000); verify/logs/v02_link.log; verify/out/v03_s6_rank.json

### F3: S98-only FAP 0.35 is method-dependent

**CONFIRMED**, corrected severity **MINOR**

Reproduced with my own BLS grid and null (200 segment-cyclic shifts, B removed). QLP DET_FLUX gives B as the global maximum, SDE 8.58, 0/200 null ≥ (FAP <5e-3). SYS_RM_FLUX gives SDE 6.04, 50/200 (FAP 0.25), the manuscript's column. I also found why the SYS_RM null is inflated, which the auditor did not identify. After biweight detrending, SYS_RM keeps segment-edge dips: BTJD 4017.67 (S/N 8.2, within ~40 min of the segment end), 3995.18 and 4002.08. The cyclic-shift null moves these into the segment interior, so the null S/N maximum reaches about 11 (DET null: about 5.9). The manuscript's 0.35 is therefore an artefact-inflated, conservative value. Its conclusion ('alone they do not establish it') still holds after the 126k-star search.

*Evidence:* verify/out/v06_null_DET_FLUX.json, verify/out/v06_null_SYS_RM_FLUX.json, verify/logs/v06_*.log, verify/out/v07_s98_features.json

### F4: Scan window ±0.004 d not justified by S98-only σ_P (numeric row status DISCREPANT)

**PARTIAL**, corrected severity **MINOR**

The core point is valid. From my own E0/E1 times, the S98-only σ_P is 0.0030 d (0.0015 d if E4 is added), so ±0.004 d is only 1.3–2.7σ, and the manuscript does not say how the window was chosen. One piece of evidence is an artefact, though. The claim that the auditor's BLS peak (11.1416) lies outside the window comes from the auditor's coarse period grid (step dP = dmin·P²/(T·5) ≈ 0.021 d at 11.13 d). My finer grid puts the S98 BLS peak at 11.1310 (DET) and 11.1356 (SYS_RM), both inside the window. 'DISCREPANT' is too strong; 'missing justification' is the right status.

*Evidence:* 04_B_existence/s98data.py period_grid; logs/b02.log; verify/out/v06_null_*.json obs.P_B; verify/logs/v04.log S98 E0/E1 times

### F5: 'Partial' S98 transit mislabelled; error understated

**CONFIRMED**, corrected severity **MINOR**

The manuscript's own per-transit table gives E4 coverage 0.50, e_D_white 245 ppm and e_D_perm 216 ppm. The quoted ±161 ppm is the set-wide full-coverage injection error, applied to all three events. With the manuscript ephemeris, 48% of E4's in-transit cadences carry QLP bit 2^30, interleaved in time; QLP docs define this as the internal poor-quality flag. My empirical random-epoch error for E4 on the q==0 data is 247–260 ppm, against about 125 ppm for full events. 'Partial' means half the cadences are masked, not partial coverage in time.

*Evidence:* rnaas_tic4206066_pubprep/analysis/out/a1_per_transit.csv row B,B_mann,pipe,S98,epoch 4; verify/out/v01_s98_events.json; 00_data/docs/qlp/mit_qlp.html.txt (2^30 = QLP quality flag)

### F6: 'Two strongest S6 dips' depends on photometry

**PARTIAL**, corrected severity **NOTE**

With my own estimator, the two predicted dips are the top two boxes in TESS-SPOC PDCSAP, SAP and QLP KSPSAP under all 4 detrendings. In TGLC aperture photometry they are top two in 2 of 4 detrendings. In the other 2, the only box that outranks event 1 is at 1481.24–1481.27, 0.5–1.2 h from the 1481.29 momentum dump: an instrumental artefact, not a competing dip. The manuscript states the claim for its own S6 light curve (TESS-SPOC, Fig. 1c), where it is true. Adding 'TESS-SPOC' costs nothing but is not required.

*Evidence:* verify/out/v03_s6_rank.json; verify/logs/v03.log

### F7: removed-null anti-conservative; 'kept' null is the correct permutation null

**PARTIAL**, corrected severity **NOTE**

Two parts are right. Epochs −1 and 3 fall in straylight intervals (bits 2048, 2^29 and 2^30 set; v01). Removing B's windows makes the null lighter: my REM null maximum is 7.9, against 11.42 for the kept/EXCL null. The claim that the auditor's 'kept' null is the correct permutation null is wrong. It excludes |P0-Pc|<0.1 d, which removes exactly the ephemerides that can reproduce the alignment. For S6 alone that share dominates (FAP ≈1.1e-3, not 0/50,000); for S6+S32 it is small (≈7e-6).

*Evidence:* verify/out/v02_link.json; verify/out/v01_s98_events.json

### Row: S98 clean transit depths (MINOR_DIFF)

**CONFIRMED**, corrected severity **MINOR**

My box depths use a local linear baseline, the manuscript ephemeris and q==0. E0 = 902–919 ± 128, E1 = 1094–1112 ± 123, E4 = 866–890 ± 250 ppm. Dividing by about 0.85 to convert box depth to trapezoid depth gives about 1060, 1300 and 1040 ppm. Both the manuscript and auditor values agree with these within 1σ. The only real issue is E4's understated error (F5).

*Evidence:* verify/out/v01_s98_events.json

### Row: link-test peak S/N 9.75 vs 10.0 (MINOR_DIFF)

**CONFIRMED**, corrected severity **NOTE**

The value depends on the estimator. My in-box, beta-calibrated estimator gives 11.50 at P = 11.132731; the auditor's flank-baseline estimator gives 9.75 and the manuscript's gives 10.03. Neither is an error; what matters is the significance of each value against its own null.

*Evidence:* verify/out/v02_link.json runs[0].obs

### Row: link FAP S6+S32 (manuscript 2e-3 conservative; auditor <2e-5, tail fit ~5e-6)

**CONFIRMED**, corrected severity **NOTE**

With my implementation: EXCL null 0/50,000, REM 0/50,000, and a chance-alignment share of 8/50,000 of the strip (about 7e-6 of the full prior). The manuscript's 2e-3 is therefore conservative. The auditor's explanation also holds: the manuscript's S32 null tail is driven by a box at BTJD 2200.00 (S/N 10.4) in the last ~1 h of S32 data (data end 2200.045, QUALITY = 0), a detrending edge effect. However, my EXCL null maximum (11.42) almost reaches the observed 11.50, and the auditor's (9.37 vs 9.75) is similarly close. The auditor's tail-fit ~5e-6 is optimistic; a safer statement is FAP ~1e-5 to 1e-4 with clean data handling.

*Evidence:* verify/out/v02_link.json runs[0]; verify/out/v05_s32_edge.json; pubprep a3_signalB.json T2 'pipe S32' top_other 2200.00 S/N 10.38

### Row: S6-only link FAP (auditor 'FAP ≲7e-4, resolution-limited')

**REFUTED**, corrected severity **MINOR**

1/N_distinct is a lower bound on the resolvable p-value, not an upper bound on the FAP. Measured directly by including the chance-alignment strip, the S6-alone FAP is about 1.1e-3 (rank tests: 1.5–2e-3). The manuscript's '0 of 10,000' is still overstated, as in F2, but the auditor's replacement value is too strong by about 2x.

*Evidence:* verify/out/v02_link.json runs[1]; verify/out/v03_s6_rank.json

### Spot-check MATCH: P_B and T0_B

**CONFIRMED**, corrected severity **NONE**

Independent fit of the same 4 events: a cadence-integrated trapezoid on S6 PDCSAP and on S98 SYS_RM q==0 with my own detrending, with injection-calibrated errors. Result: P = 11.132747 ± 0.000012 d, T0 = BTJD 4000.8658 ± 0.0015 (χ² 1.9/2). This agrees with the manuscript (11.13274 ± 0.00003, 4000.865 ± 0.003) and with the auditor. The manuscript's errors are about 2x larger, which is conservative. No circularity: b08 imports only the auditor's own modules and raw FITS.

*Evidence:* verify/out/v04_ephem.json; verify/logs/v04.log; 04_B_existence/b08_ephem.py imports

### Spot-check MATCH: overall significance ≈3σ

**PARTIAL**, corrected severity **NOTE**

The arithmetic checks out. Fisher with (0.35, 2e-3) gives p = 5.8e-3, or 2.52σ one-sided. (1.6e-3, 2e-3) gives 3.92σ. The link alone at 2e-3 gives 2.88σ. The method is questionable, though: the S98 per-star p is not informative after the target was picked from a 126,246-star search, so Fisher-combining it dilutes the result, and the overall evidence should rest on the archival link test. That gives ≈2.9σ with the manuscript's null and ≈3.8–4.1σ with cleaner nulls. The '≈3σ, low-significance candidate' wording is supported and conservative.

*Evidence:* scipy chi2/norm computation in this verification; verify/out/v02_link.json

### Spot-check MATCH: S32 predicted events

**CONFIRMED**, corrected severity **NONE**

My box depths at the S32-only best period are 734 ± 203 (2175.09) and 916 ± 238 ppm (2197.36); the middle event (2186.22) falls in the data gap. The S32-only FAP is about 0.05 (manuscript 0.032, auditor 0.02). The events are consistent but weak on their own.

*Evidence:* verify/out/v02_link.json runs[2]

## Issues the auditor missed

### [MINOR] Exclusion-zone nulls omit the chance-alignment share (affects the S6 claim and the auditor's recommended wording)

Two nulls cannot reproduce the observed configuration: the manuscript's (B windows removed) and the auditor's 'kept' null (|P0-Pc|>0.1 d excluded). The honest permutation FAP adds the fraction of ephemeris space whose scan hits B's own archival dips. For S6 alone this share dominates: FAP ≈1.1e-3, matching the rank tests' 1.5–2e-3. So the auditor's proposed Sec. 3 text '(FAP ≲10⁻³)' for the S6 dips is too strong, and the correct wording is '≈2×10⁻³ (rank test)'. For S6+S32 the share is small (≈7e-6) and the conclusion is unchanged.

*Evidence:* verify/v02_link.py NEAR mode; verify/out/v02_link.json

### [NOTE] Cause of the high S98-only FAP (0.35): detrending edge artefacts in SYS_RM_FLUX

With SYS_RM_FLUX plus a 0.75-d biweight, the strongest non-B box in S98 is an edge dip at 4017.67 (S/N 8.2) that outranks B's E0; further edge dips sit at 3995.18 and 4002.08. QLP DET_FLUX has none of these. The cyclic-shift null moves them into the segment interiors and inflates the null. With DET_FLUX the S98-only FAP is <5e-3 (verifier) and 1.6e-3 (auditor). The manuscript's 0.35 is conservative, and a single sentence naming the flux column would pre-empt a referee query.

*Evidence:* verify/out/v07_s98_features.json; verify/out/v06_null_*.json

### [NOTE] Auditor's coarse BLS period grid

The auditor's BLS grid steps ≈0.021 d at 11.13 d. That makes its S98 BLS periods (11.1416 etc.) and the stated '11.1316–11.1416 spread' partly discretisation artefacts. They should not be used to argue the ±0.004 d window is too narrow.

*Evidence:* 04_B_existence/s98data.py period_grid (os_=5, dmin=1/24, T≈50 d)

### [NOTE] Selection history and trials for the archival test

An earlier project-level archival check classified B as 'signal absent/not locked in other sectors': an unconstrained narrow BLS on S6+S32 drifted 61.7 h in epoch. The anchored ±0.004 d link test in the manuscript was designed afterwards. A strict reader could count two archival tests, a trials factor of about 2, making the manuscript FAP about 4e-3. The per-star S98 p also carries no weight after the 126k-star search. Neither changes the '≈3σ low-significance candidate' framing, but the per-candidate nature of the FAP should be explicit.

*Evidence:* addenda/pool_multisector/SUMMARY.md and multisector_results.json (tic4206066_p11.132976: dT0_vs_s98_h 61.7)

### [NOTE] Tail of the S6+S32 'kept' null almost reaches the observed value

In both implementations the kept/EXCL null maximum almost reaches the observed statistic: 11.42 vs 11.50 (verifier) and 9.37 vs 9.75 (auditor). This comes from partial hits on one B dip. A FAP of ~1e-5 to 1e-4 is defensible; the auditor's exponential-tail extrapolation of 5e-6 is not.

*Evidence:* verify/out/v02_link.json runs[0].EXCL.q; 04_B_existence/logs/b04_main.log
