# Verification of audit 08_fpp — VERIFY_REPORT

> Rendered verbatim from the agent's structured result (`_structured/verify__08_fpp.json`); the sub-agent harness blocked agents from writing .md report files. Scripts, outputs and logs are in this folder.


## Overall assessment

The FPP numbers are sound. I reproduced them independently with my own fold, a fresh TRICERATOPS target and different methods, and they agree with both the manuscript and the auditor within systematics:
- 3.18 d signal: default FPP 0.036 +- 0.004 with 5-min bins and 0.044 with 10-min bins. MOLUSC mode 0.006-0.009.
- 11.13 d signal: default 0.213 +- 0.019; MOLUSC mode 0.05-0.09.
- NFPP is at most 1.4e-11 and STP makes up about 99% of the FPP in every run.

The auditor's two MAJOR findings should both be downgraded.

**F1 is only partly right.**
- Correct: 'drops the companion-occurrence prior' does not lower the FPP (n_kept/N = 0.13 against f_comp = 0.14).
- Wrong: the claim that the Gaia limits, not the 10 AU cut, cause the drop. Gaia RUWE and contrast survivors without the 10 AU cut give a higher FPP than the default: 0.045-0.050 for A and 0.24-0.27 for B. The periastron > 10 AU cut is what removes the close massive companions, and the manuscript already says so.
- Consequence: the manuscript needs only a MINOR wording fix, dropping or rephrasing the occurrence-prior clause. The auditor's proposed replacement would make the text wrong and must not be used.

**F2's facts are confirmed.**
- About 90% of the STP evidence for A comes from companions above 0.38 Msun.
- Truncating the companion prior at the manuscript's own 0.38 Msun limit gives FPP_A = 0.004 (0.007-0.009 at 0.45 Msun, 0.017 at 0.54 Msun). B stays at 0.16.
- I rate it an upper-end MINOR. 'Neither signal meets FPP < 0.015' is literally true for the tool's calibrated configuration, the error is conservative, and not claiming validation remains correct.
- A caveat sentence is still strongly recommended.

**F3 and F4 are confirmed as MINOR.** The fold and binning range is too narrow, and the criterion should be quoted as FPP < 0.015 and NFPP < 1e-3, with the note that it was calibrated for circular orbits.

**F5 is confirmed but is only a NOTE.**

**Spot checks of MATCH claims.** I found no circularity or code bugs in the auditor's hooks or indices. The MOLUSC survivor file is the auditor's own fresh run; I reused it (sha256 recorded) and did not re-run MOLUSC.

Net verdict for audit 08: PASS_WITH_MINOR_CHANGES. The recommended MOLUSC and overluminosity wording is in missed_issues[0].

## Checks

### F1 (MAJOR): the stated reason MOLUSC lowers the FPP is incorrect; the real cause is the Gaia RUWE and contrast limits

**PARTIAL**, corrected severity **MINOR**

The code facts check out. The a(1-e) > 10 AU filter is at marginal_likelihoods.py l.948, the clamp and pad at l.950-951, and lnprior_companion = zeros at l.1038. The auditor is also right on one point: 'drops the companion-occurrence prior' is not a lowering factor. n_kept/N = 0.1295, roughly the default f_comp(<2.2") = 0.1413, and assuming a companion for every star can only push the FPP up.

The auditor's central claim is refuted, though. The auditor says the 10 AU cut is not the cause and that the real cause is the Gaia RUWE and contrast limits. I checked this with the Gaia survivors and no periastron cut:
- Using the same MOLUSC survivors (all 429797, Gaia limits applied) without the 10 AU cut, the prior mass in companions above 0.38 Msun is 0.085. That is higher than the default 0.046.
- 87% of the surviving >0.38 Msun companions have a < 10 AU, and 60% have a < 1 AU. Gaia astrometry and imaging cannot exclude them. The 10 AU cut removes them and brings the prior mass to 0.0083.
- Direct TRICERATOPS runs with all Gaia survivors (no 10 AU cut) give FPP_A = 0.045 +- 0.006 (3 runs) and FPP_B = 0.237 (2 runs). The fixed-mass grid gives 0.049-0.050 and 0.267. Both are higher than the default FPPs on the same fold (0.036 and 0.213).
- Relative to the default, the Gaia-only STP evidence is 1.73x for A and 1.87x for B. With the 10 AU cut as coded it is 0.163-0.167x for A and 0.31x for B.

So the 10 AU clause in the manuscript is correct and is the dominant cause. Only the second clause ('drops the companion-occurrence prior') is wrong as a reason for the drop. The auditor's proposed fix ('mainly because the Gaia limits remove companions above ~0.4 Msun') would introduce a factual error and must not be adopted. The manuscript fix is a wording change (MINOR): replace 'and drops the companion-occurrence prior' with 'which removes close, massive companions that Gaia astrometry cannot exclude; it also assumes every star has a companion'.

*Evidence:* Files in /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp/verify/:
- v05_molusc_dist.json (MOLUSC file sha256 63cdae36..., reused from the auditor's fresh run):
  - Effective prior above 0.38 Msun: default 0.0456; MOLUSC kept 0.0083; MOLUSC all survivors 0.0846.
  - Massive survivors: frac_a_lt_10AU 0.871, frac_a_lt_1AU 0.595.
- runs_v/molusc_gaiaonly_A_*.json: FPP 0.0406 / 0.0432 / 0.0516.
- runs_v/molusc_gaiaonly_B_*.json: FPP 0.232 / 0.242.
- runs_v/default_A (8 runs, FPP 0.0364) and default_B (6 runs, FPP 0.213).
- v06_combine.json and v06_combine_N2e6.json: Z_STP_ratio_molusc_gaia_only_over_default 1.72-1.73 (A) and 1.87 (B); coded-mode ratio 0.163-0.167 (A) and 0.31 (B).
- Method: v07_molusc_direct.py reproduces MOLUSC mode by patching the q sampler and setting the prior to zero. It does not use the molusc_file code path or a pandas shim.

### F2 (MAJOR): the unconstrained FPP is dominated by companions the manuscript itself excludes, and 'neither signal meets FPP<0.015' is misleading for A

**CONFIRMED**, corrected severity **MINOR**

I reproduced this with two methods independent of the auditor's per-draw hooks, on my own fold and target (Teff 4012, 5-min bins):
- (a) Direct TRICERATOPS runs with TRICERATOPS' own companion-mass prior truncated by rejection at M <= 0.38 Msun, then adding ln P(M <= 0.38) to the PTP/PEB/STP/SEB evidence. This gives FPP_A = 0.0043 +- 0.0006 (8 runs). At M <= 0.45 it gives 0.0085 +- 0.0015. For B it gives 0.158 +- 0.006.
- (b) Fixed-companion-mass likelihood grids. These reproduce the default lnZ_STP to 0.01-0.06 dex. They give an STP evidence share above 0.38 Msun of 90-91% for A (auditor: 94%) and 30% for B (auditor: 29%). Truncated FPP_A is 0.0036-0.0039 at 0.38, 0.0071-0.0074 at 0.45, and 0.017 at 0.54.

The auditor's hook code has no bug: sample_q_companion, sample_ecc and _log_mean_exp are each called once per scenario on full-length aligned arrays; the prior is constant; the FPP indices 0/3/9 are correct.

The facts stand, but I downgrade the severity to MINOR (an upper-end MINOR):
- 'Without it neither signal meets FPP < 0.015' is literally true for the tool's calibrated default configuration.
- The error is in the conservative direction, and the candidate status is unaffected.
- The result depends on a hard 2-sigma mass cut. At 0.54 Msun FPP_A is already 0.017, above 0.015.

One imprecision in F2: 'the manuscript itself excludes both ... circular orbits'. The auditor's own data (weight of M <= 0.38 with e < 0.1 is about 1e-9) confirm the manuscript's statement that circular orbits need hosts above 0.54 Msun. Only the 0.38 Msun limit does work here.

A caveat sentence is still strongly recommended.

*Evidence:* Files in /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp/verify/:
- runs_v/trunc0.38_A_v_A_S32+S98_5min.json: 8 runs, 0.0034-0.0051.
- runs_v/trunc0.45_A_*.json: 0.0063-0.0095.
- runs_v/trunc0.38_B_*.json: 0.151-0.167.
- v06_combine.json and v06_combine_N2e6.json: STP_evidence_share_gt0.38_default A 0.902/0.912, B 0.295; FPP_default_trunc0.38 A 0.0036/0.0039; FPP_default_trunc0.54 A 0.0168/0.0173.
- Code: v03_runs.py (truncated mode), v04_grid.py.

### F3 (MINOR): fold and method sensitivity of FPP_A understated (manuscript 0.031 +- 0.003, '0.031-0.034 across three folds')

**CONFIRMED**, corrected severity **MINOR**

On my independent fold (1.0 d biweight window, mask of +-0.75 T14, Teff 4012, SPOC aperture from get_spoc_apertures, QLP 3x3 aperture):
- 5-min bins: FPP_A = 0.0364 +- 0.0036 (8 runs, range 0.032-0.043).
- 10-min bins: 0.0437 +- 0.0037 (4 runs).
This matches the auditor (0.034 with 6-min bins, 0.043 with 10-min bins).

The bin width, which also sets the exptime passed to TRICERATOPS, alone moves the FPP by about 25%. The manuscript does not state it. The headline 0.031 sits at the low end, and the quoted range is too narrow. 'FPP = 0.03-0.04' is the fair wording.

*Evidence:* - verify/runs_v/default_A_v_A_S32+S98_5min.json
- verify/runs_v/default_A_v_A_S32+S98_10min.json
- verify/folds/folds_meta.json: sigma_oot 136/95 ppm, depth 472/468 ppm.

### F4 (MINOR): validation criterion quoted incompletely; calibration caveat missing

**CONFIRMED**, corrected severity **MINOR**

The two sources say:
- The Giacalone et al. (2021) abstract: 'validated planets (FPP < 0.015 and NFPP < 10^-3)'.
- The paper's Section 2 assumes 'a fixed orbital period and zero eccentricity'.

TRICERATOPS 1.1.0 instead samples planet eccentricities from beta(0.867, 3.030) (priors.py sample_ecc). The eccentric draws are exactly what let low-mass hosts survive here, so the calibration caveat matters.

The NFPP omission is harmless for 'neither meets'. NFPP is at most 1.4e-11 in all my runs. But NFPP must be cited if the text ever says a signal meets the criterion.

*Evidence:* - 08_fpp/docs/G21.txt: abstract (l.~81) and the Section 2 text on zero eccentricity.
- trice_env triceratops/priors.py l.134-155.
- verify/v_summary.json: NFPP_max values.

### F5 (MINOR): the MOLUSC-mode FPP depends on N and has high scatter

**CONFIRMED**, corrected severity **NOTE**

The code confirms it. Kept q values are zero-padded to N (l.951), and _log_mean_exp divides by N_total, so Z_STP scales as n_kept/N.

The scatter is large. My direct coded-mode runs give A 0.0094 +- 0.0015 (3 runs) and B 0.069 +- 0.027 (2 runs). The grids give 0.0063-0.0065 and 0.076. The auditor got 0.0055 +- 0.0036 and 0.069 +- 0.016. All are consistent with the manuscript's 0.005 +- 0.003 and 0.06 +- 0.02.

N = 1e6 equals MOLUSC's n_generated, which is the natural choice, and the manuscript's '(drops the companion-occurrence prior)' already gestures at this. Stating N and the every-star-has-a-companion assumption is a reproducibility note, not a needed correction.

*Evidence:* - marginal_likelihoods.py l.949-951 and l.1095.
- _numerics.py _log_mean_exp.
- verify/runs_v/molusc_coded_A_*.json and molusc_coded_B_*.json.
- verify/v06_combine*.json: FPP_molusc_as_coded.

### Row MINOR_DIFF: FPP of the 3.18 d signal, headline (manuscript 0.031 +- 0.003; auditor 0.0342 +- 0.0035)

**CONFIRMED**, corrected severity **MINOR**

My independent value is 0.0364 +- 0.0036 with 5-min bins and 0.0437 with 10-min bins. The manuscript value lies within the fold and binning systematics but at their low end. Same issue as F3.

*Evidence:* verify/runs_v/default_A_*.json

### Row MINOR_DIFF: FPP range across folds, 3.18 d (manuscript 0.031-0.034)

**CONFIRMED**, corrected severity **MINOR**

My S32+S98 folds alone span 0.032-0.048 in single runs and 0.036-0.044 in means. The quoted range is too narrow. Merge with F3.

*Evidence:* verify/v_summary.json direct_runs

### Row DISCREPANT: MOLUSC-mode effective companion prior vs default ('drops the occurrence prior' as the cause)

**PARTIAL**, corrected severity **MINOR**

The numbers are right: f_comp is 0.1413 against n_kept/N of 0.1295, and prior mass above 0.38 Msun is 0.0456 against 0.0083. The auditor's causal reading is wrong: the 0.0083 comes from the 10 AU cut, not from the Gaia limits. Gaia survivors without the cut carry 0.0846. See F1.

*Evidence:* verify/v05_molusc_dist.json

### Row MINOR_DIFF: validation criterion (manuscript quotes FPP < 0.015 only)

**CONFIRMED**, corrected severity **MINOR**

G21 defines validated planets as FPP < 0.015 and NFPP < 1e-3. The omission is harmless here. See F4.

*Evidence:* 08_fpp/docs/G21.txt abstract

### Row DISCREPANT: FPP_A with the manuscript's own companion constraint M <= 0.38 Msun (auditor 0.002-0.007)

**CONFIRMED**, corrected severity **MINOR**

Independent values:
- Direct runs with a truncated prior: 0.0043 +- 0.0006.
- Grid: 0.0036-0.0039.
- M <= 0.45: 0.0071-0.0085.
- M <= 0.54: 0.017, above 0.015.

A stays below 0.015 only for mass limits at or below about 0.5 Msun. The auditor's quoted range should read about 0.002-0.009 for limits of 0.38-0.45 Msun. B stays at 0.16 (confirmed).

*Evidence:* - verify/runs_v/trunc0.38_A_*.json and trunc0.45_A_*.json.
- verify/v06_combine*.json.

### Spot-check MATCH: MOLUSC-mode FPP 0.005 +- 0.003 (A) and 0.06 +- 0.02 (B)

**CONFIRMED**, corrected severity **NONE**

I reproduced these by a route independent of the auditor's pandas-3 shim: a q-sampler patch with the prior set to zero inside the bound-companion scenarios, plus the fixed-mass grid. Results are A 0.006-0.009 and B 0.05-0.09, consistent within the large scatter.

Circularity: the MOLUSC survivor file is the auditor's own fresh MOLUSC b26efbf output, not an earlier-package product. Its survivor count (429797) differs from the earlier run (429613), so it is a new realization. I did not re-run MOLUSC myself, so MOLUSC reproduction itself rests on the auditor.

*Evidence:* - verify/runs_v/molusc_coded_*.json
- verify/v06_combine*.json
- sha256 63cdae3625723da694d6ecb1d654e1bdaf9ab7d56008a4660bb2b85db2990998

### Spot-check MATCH: FPP of the 11.13 d signal, 0.13-0.20

**CONFIRMED**, corrected severity **NOTE**

My own fold and target give 0.213 +- 0.019 (6 runs, range 0.190-0.243). The auditor's variants also reach 0.21. The manuscript's upper end of 0.20 is slightly low; 0.13-0.22 would be safer. This is rounding level.

*Evidence:* verify/runs_v/default_B_v_B_S32+S98_5min.json

### Spot-check MATCH: 2.2 arcsec default, STP dominance, NFPP negligible, TRICERATOPS network setup

**CONFIRMED**, corrected severity **NONE**

- lnprior_bound_TP(M_s, plx, 2.2) is hard-coded when there is no contrast curve (l.1006-1007). G21 attributes 2.2" to Brown et al. 2018.
- My own target build (TIC query, TESScut and Gaia DR3 background, done fresh over the network in 46 s; 121 stars; SPOC S32 aperture of 18 px via get_spoc_apertures) gives STP >= 99% of the FPP and NFPP <= 1.4e-11.
- No reuse of earlier derived products was found in the auditor's f01/f02. Inputs are the fresh 00_data files only.

*Evidence:* - verify/v02_target.log
- verify/v_summary.json
- marginal_likelihoods.py l.1002-1010

## Issues the auditor missed

### [MAJOR] The auditor's recommended wording for the MOLUSC sentence would introduce a factual error

Both the auditor's minimum edit ('mainly because the Gaia limits remove companions above ~0.4 Msun, which dominate the default FPP') and the full recommended paragraph ('removing them, or adding the Gaia RUWE and contrast limits through MOLUSC, gives FPP ~0.002-0.007') attribute the MOLUSC drop to the Gaia limits. With the Gaia RUWE and contrast limits alone (all MOLUSC survivors, no 10 AU cut), TRICERATOPS gives FPP_A of about 0.045-0.050 and FPP_B of about 0.24-0.27, both higher than the defaults.

The drop comes from the periastron > 10 AU cut. The manuscript already says this. The cut removes close (mostly a < 1 AU) massive companions that Gaia astrometry cannot exclude but that the XP overluminosity limit does. Suggested wording:

'MOLUSC (Wood et al. 2021) with the Gaia RUWE and imaging limits lowers these to 0.005 +- 0.003 and 0.06 +- 0.02, but only because this mode keeps companions with periastron beyond 10 AU, which removes close massive companions that Gaia astrometry cannot exclude (it also assumes every star has a companion). For the 3.18 d signal ~90% of the default FPP comes from companions above 0.38 Msun, which the overluminosity limit excludes; truncating the companion prior there gives FPP ~0.004 (0.16 for the 11.13 d signal). These constrained values lie outside the configuration for which Giacalone et al. (2021) calibrated FPP < 0.015 and NFPP < 10^-3, so I do not claim validation.'

*Evidence:* - verify/runs_v/molusc_gaiaonly_{A,B}_*.json
- verify/v05_molusc_dist.json: survivors_M_gt_0.38 a-quantiles, median 0.70 AU.
- verify/v06_combine*.json

### [MINOR] The MOLUSC-mode FPP of A is itself dominated by companions above 0.38 Msun

Even in MOLUSC mode as coded, 82-84% of the STP evidence for the 3.18 d signal comes from companions above 0.38 Msun. The auditor's own A_moluscFresh5 run shows 83%. Combining MOLUSC mode with the overluminosity limit gives FPP_A of about 0.001.

So the manuscript's constrained values are all driven by the same channel, close or massive companions, and the 0.38 Msun limit is the decisive constraint. This is worth one clause if the text is revised.

*Evidence:* - verify/v06_combine*.json: STP_evidence_share_gt0.38_molusc_kept 0.82-0.84; FPP_molusc_kept_and_M_le_0.38 0.0010-0.0012.
- 08_fpp/out_fpp_summary.json: A_moluscFresh5 STP_mass_bins.

### [MINOR] The constrained FPP_A depends on where the companion-mass limit is placed

FPP_A rises from 0.004 at M <= 0.38 to 0.007-0.009 at 0.45 and 0.017 at 0.54 Msun. The 0.38 Msun limit is a 2-sigma bound applied as a hard cut, which is optimistic.

Any sentence implying that A would meet FPP < 0.015 should give a range (about 0.004-0.009 for limits of 0.38-0.45 Msun) and not a single number. The auditor's 0.002-0.007 is slightly low at the top end.

*Evidence:* - verify/v06_combine.json: FPP_default_trunc0.30/0.38/0.45/0.54 for A = 0.0011/0.0036/0.0071/0.0168.
- verify/runs_v/trunc0.45_A_*.json: 0.0085 +- 0.0015.

### [MINOR] Bin width and exptime used for the TRICERATOPS input are not stated

FPP_A changes by about 25% between 5-6 min bins (0.034-0.036) and 10-min bins (0.043-0.044). The exptime passed to TRICERATOPS equals the bin width. For reproducibility the manuscript or package should state the bin width, the exptime and the sigma definition.

*Evidence:* - verify/runs_v/default_A_v_A_S32+S98_{5min,10min}.json
- 08_fpp/out_fpp_summary.txt: A_base vs A_10min.

### [NOTE] F2's 'excludes circular orbits' phrasing is imprecise

The near-circular STP weight, 51% for A, sits on hosts above 0.54 Msun. This agrees with the manuscript's statement that circular orbits need hosts above 0.54 Msun; the weight of M <= 0.38 with e < 0.1 is about 1e-9. Those draws are excluded only through the 0.38 Msun overluminosity limit. The manuscript's circular-orbit statement is thus independently supported, not contradicted.

*Evidence:* 08_fpp/out_fpp_summary.json: A_base STP_ecc_bins; the auditor's wfrac_mass_le_0.38_ecc_lt_0.1

### [NOTE] VERIFY_REPORT.md not written

The subagent rules forbid writing report .md files (the auditor's REPORT.md was blocked the same way). The full verification is in this structured output.

Scripts and outputs are in /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp/verify/:
- v01_fold.py to v08_summary.py
- folds/, runs_v/, v05_molusc_dist.json, v06_combine.json, v06_combine_N2e6.json, v_summary.json, v02_target.log, q*.log

The B high-N grid (runs_v/grid_B_*_N2e6.json) is partial and was stopped because of machine load. It was not used.

*Evidence:* verify/ directory listing
