# Verification of audit 09_companion — VERIFY_REPORT

> Rendered verbatim from the agent's structured result (`_structured/verify__09_companion.json`); the sub-agent harness blocked agents from writing .md report files. Scripts, outputs and logs are in this folder.


## Overall assessment

The auditor's numbers hold up under independent re-derivation. The main finding holds: the 'circular orbit needs > 0.54 Msun, survives only eccentric' sentence applies to the 3.18 d signal only, and for the 11.13 d signal a circular orbit on any companion of at least 0.14 Msun fits. That remains a MINOR wording fix in Sec. 4.

What I recomputed independently:
- Own detrending, a different S98 flux column, and TRICERATOPS-based dilution give rho < 3.60 rho_sun, a circular minimum host mass of about 0.517 Msun, e_min of about 0.20-0.53 for 0.17-0.37 Msun hosts, and a free circular fit for signal B.
- The Gaia RV error formula agrees with the official datamodel.
- TRICERATOPS 1.1.0 STP does draw eccentric orbits: Beta(0.867, 3.03), with 32% of companions above 0.38 Msun.
- The overluminosity limit re-derives to 0.37-0.38 Msun.
- The RUWE percentile the auditor could not finish is now done: 61st percentile, typical.

Severity changes:
- F2 downgraded to NOTE. The manuscript never claims TRICERATOPS is circular-only, and 'untested' refers to the observational tests.
- F3 and the 3.2 km/s and 0.016 Msun rows downgraded to NOTE. The critique is right (K95 is about 1.2 km/s, and 3.2 km/s is peak-to-peak), but those numbers appear only in package notes, not in the manuscript.

Nothing FATAL or MAJOR. The manuscript's Sec. 4 numbers (RUWE 1.21, no NSS solution, RV exclusion of an eclipsing companion, 0.38 Msun, 0.54 Msun as approximately 0.52) are accurate. The Way et al. overluminosity value itself (+0.003 ± 0.07) remains unverified.

VERIFY_REPORT.md was not written because this session's instructions prohibit writing report .md files. Scripts, JSON outputs, logs and PROVENANCE.json are in /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/09_companion/verify/ (v1_trice_prior.py/.json, v2_hosts.py/.json, v3_fit.py, v3_fit_A.json, v3_fit_B.json, v3_A.log, v3_B.log, v4_ruwe_sample.csv, v4_ruwe_pct.json, gaia_source_dm.html).

## Checks

### F1 (MINOR): the circular-orbit host-mass argument holds only for the 3.18 d signal; for 11.13 d a circular STP is allowed

**CONFIRMED**, corrected severity **MINOR**

I refit independently with my own per-transit quadratic detrending instead of wotan, used QLP DET_FLUX instead of SYS_RM_FLUX, took the dilution from TRICERATOPS flux_relation, and used Mann15 radii. For signal B, a circular orbit on a diluted companion host fits for every host of at least 0.14 Msun (dchi2 = +1.5 at 0.14, and -0.2 to +0.2 for 0.17-0.606). Only 0.10 Msun is disfavoured (dchi2 = 22.7). Sec. 3 of the manuscript already says the 11.13 d light curve does not constrain b or rho. However, the Sec. 4 sentence ('survives only for eccentric orbits') comes straight after the 11.13 d localisation and reads as covering both signals. For signal A my circular-orbit dchi2 is 8.3 at 0.50, 3.2 at 0.52, 1.3 at 0.53 and 0.02 at 0.54 Msun, so the dchi2 = 4 crossing is about 0.517 Msun, against the auditor's 0.523. The auditor's scope point stands.

*Evidence:* verify/v3_fit_B.json, verify/v3_fit_A.json, verify/v3_B.log; manuscript.tex l.35 and l.39

### F2 (MINOR): TRICERATOPS STP already includes eccentric orbits, so 'untested' framing is imprecise

**PARTIAL**, corrected severity **NOTE**

The code facts check out. lnZ_STP draws eccs = sample_ecc(planet=True), which is Beta(0.867, 3.03), with uniform argp. These are passed to lnL_TP and on to simulate_TP_transit, which calls evaluate_ps(e=ecc, w=(90-argp)). Companion masses span 0.100-0.606 Msun, and 32.3% of draws are above 0.38 Msun (my own run of sample_q_companion with M_s = 0.606). Under the Beta prior, only 1.4-11.6% of draws satisfy the photo-eccentric slowdown needed for e_min = 0.5-0.19. The auditor's premise is wrong, though. The manuscript never says TRICERATOPS assumes circular orbits. 'Main untested false-positive scenario' in the abstract means not excluded by the observational tests listed, and the FPP is explicitly qualified as 'under its default assumptions', with the MOLUSC variant discussed. The only substantive gap is that the default FPP does not apply the Gaia overluminosity cut. Mentioning that would be a helpful clarification, but the existing wording is not inaccurate.

*Evidence:* trice_env triceratops/marginal_likelihoods.py lnZ_STP (sample_q_companion call, lnprior_bound_TP with seps = 2.2, eccs = sample_ecc(...planet=True...), lnL_TP(..., eccs[i], argps[i], ...)); likelihoods.py simulate_TP_transit evaluate_ps(e=ecc, w=(90-argp)); verify/v1_trice_prior.json

### F3 (MINOR): 3.2 km/s is peak-to-peak noise, not a K limit; calibrated K95 about 1.16 km/s; context for rv_chisq_pvalue

**CONFIRMED**, corrected severity **NOTE**

I checked the official Gaia DR3 datamodel HTML (MathML parsed). For rv_method_used = 1, sigma_med = sqrt(pi/2)·sigma(V_j)/sqrt(N) with a 0.11 km/s floor. That gives sigma(V_j) = 0.877 km/s from the fresh 00_data values (radial_velocity_error 0.2207, N = 33, rv_method_used = 1). Even the most conservative bound, which attributes all the scatter to a sinusoid (K <= sqrt(2)·sigma = 1.24 km/s), is about 2.6x tighter than 3.2 km/s. rv_amplitude_robust is documented as max minus min after 3xIQR outlier removal, so it is a peak-to-peak range, not K. The Katz variability selection needs rv_renormalised_gof > 4, and the observed value is 3.0. The finding is right, but the manuscript text quotes neither 3.2 km/s nor 0.016 Msun. They appear only in MANUSCRIPT_NUMBER_AUDIT.csv, METHOD_AUDIT.md and CHANGES_FROM_CURRENT_DRAFT.md, and the manuscript sentence ('exclude an eclipsing stellar companion') is correct under either K limit. So this is a package-notes correction (NOTE), not a manuscript fix.

*Evidence:* verify/gaia_source_dm.html (sha256 800bf953...), 00_data/catalogs/gaiadr3_gaia_source_cone240as.ecsv; grep of manuscript.tex finds no '3.2' or '0.016'

### Numeric row: K limit 3.2 km/s (MINOR_DIFF)

**CONFIRMED**, corrected severity **NOTE**

As for F3. Calibrated K95 is about 1.16 km/s (auditor's Monte Carlo) and the conservative sqrt(2)·sigma(V_j) is 1.24 km/s (mine). Package-only number.

*Evidence:* auditor out/rv_limits.json; own sigma(V_j) recomputation

### Numeric row: M sin i < 0.016 Msun at 3.18 d (MINOR_DIFF)

**CONFIRMED**, corrected severity **NOTE**

0.0161 follows arithmetically from K = 3.2 km/s with M1 = 0.606. K <= 1.24 km/s gives about 0.006 Msun. The stellar-companion exclusion holds either way (a 0.08 Msun companion at 3.18 d gives K = 14.9 km/s). Package-only number.

*Evidence:* auditor c05_rv.py mass function checked by reading; the same formula at K = 1.24 gives about 0.0062 Msun

### Numeric row: minimum circular-orbit host mass 0.54 vs 0.523 Msun (MINOR_DIFF)

**CONFIRMED**, corrected severity **NOTE**

My independent fit gives about 0.517 Msun (dchi2 = 4). The manuscript's 0.54 is roughly where rho_host = M/R^3 (Mann) equals 3.35-3.5 rho_sun, i.e. the undiluted density argument. Dilution lengthens T14 slightly and lowers the threshold by about 0.02 Msun. That is a method difference with no consequence, because everything above 0.38 Msun is already excluded by the overluminosity. Writing 'about 0.52-0.54' is optional.

*Evidence:* verify/v3_fit_A.json hosts (0.50: 8.26, 0.52: 3.18, 0.53: 1.31, 0.54: 0.02)

### Numeric row: 'survives only for eccentric orbits' (DISCREPANT for signal B)

**CONFIRMED**, corrected severity **MINOR**

Same as F1. It is true for signal A: for hosts of 0.17-0.37 Msun my circular dchi2 is 53-142 and the rho_app limit is about 3.6, giving e_min of about 0.20-0.53 at omega = 270 degrees. It is false for signal B.

*Evidence:* verify/v3_fit_A.json, verify/v3_fit_B.json

### Spot-check MATCH: RUWE = 1.21, no NSS solution, no neighbour within 10 arcsec

**CONFIRMED**, corrected severity **NONE**

Independently read from 00_data/catalogs/gaiadr3_gaia_source_cone240as.ecsv, not from the auditor's files: ruwe 1.2143, non_single_star 0, ipd_frac_multi_peak 0, AEN 0.0935 mas; nearest Gaia source at 11.08 arcsec with G = 18.75. I also completed the RUWE-percentile step the auditor could not finish (Gaia TAP answered). Among 383 stars with G 10.3-10.7, BP-RP 1.6-1.9 and plx > 15 mas, RUWE 1.214 is at the 61st percentile (median 1.14), and AEN 0.093 mas is at the 42nd percentile (median AEN_sig 15.9, so AEN_sig 15.8 is typical for this brightness). The astrometry is unremarkable.

*Evidence:* verify/v4_ruwe_sample.csv (sha256 e41f3cdc...), verify/v4_ruwe_pct.json

### Spot-check MATCH: rho* < 3.5 rho_sun (97.7%) and T14 = 1.87 h

**CONFIRMED**, corrected severity **NONE**

My independent pipeline (different detrending, different S98 flux column, 23-transit windows) gives best rho = 3.18 with b = 0 and a dchi2 = 4 upper limit of rho = 3.60 rho_sun. The auditor has 3.13 and 3.53, the manuscript 3.5. The agreement is within method scatter and there is no circularity: the auditor's c02-c04 read only raw MAST FITS from 00_data. Minor caveat: Claret LD is read from the TRICERATOPS package data file, which is a copy of the published table and acceptable.

*Evidence:* verify/v3_fit_A.json target_free, target_rho_max_dchi2_4 = 3.604

### Spot-check MATCH: overluminosity -> 0.38 Msun (2 sigma)

**CONFIRMED**, corrected severity **NOTE**

I re-derived this with my own parse of the Mamajek EEM table. With the primary at the Way-predicted M_G (7.758, the auditor's convention), the limit is 0.380-0.383 Msun. With the primary at observed M_G + O_lim (7.755 + 0.14), it is 0.369-0.372. So 0.37-0.38, and the manuscript's 0.38 is the conservative end. Partial circularity: the auditor took MG1 = 7.758 from the Way catalogue row via earlier products without re-downloading it. It is consistent with the fresh observed M_G = 7.755 plus O = +0.003, but the Way value +0.003 ± 0.07 itself remains unverified by both of us.

*Evidence:* inline computation in this session (Mamajek EEM_dwarf.txt M_G-Msun: M_G 9.67 -> 0.40, 10.05 -> 0.37)

## Issues the auditor missed

### [NOTE] RUWE context step was left undone; it is now done and supports the text

The RUWE and excess-noise percentile step (task step 1) was not completed. I ran it: RUWE 1.214 is at the 61st percentile of 383 stars of similar G, colour and distance (median 1.14), and AEN is at the 42nd percentile. The manuscript could add 'typical for stars of similar magnitude and colour' to pre-empt a referee reading RUWE 1.21 as elevated.

*Evidence:* verify/v4_ruwe_pct.json

### [NOTE] Tidal-circularisation plausibility (step 7) not quantified or cited

The auditor's proposed wording calls e >= 0.2-0.5 at 3.18 d 'disfavoured' without a calculation or citation. My order-of-magnitude estimate uses the Goldreich & Soter tau_e = (4/63) Q_p (M_p/M_*) (a/R_p)^5 / n for a roughly 2 R_earth, roughly 5 M_earth planet on a 0.3 Msun host at P = 3.18 d (a about 0.028 au, a/R_p about 330). It gives tau_e of about 2 Myr for Q_p = 100 and about 200 Myr for Q_p = 1e4, much shorter than a field-star age. Such an orbit therefore needs ongoing excitation. The claim is reasonable, but a citation should be added if it enters the text (e.g., Van Eylen et al. 2019, AJ 157, 61, for low eccentricities of small planets).

*Evidence:* hand calculation in this verification

### [NOTE] Lower host-mass bound from transit shape is method dependent (0.14-0.17 Msun)

For signal A the auditor finds that hosts below 0.168 Msun are disfavoured at dchi2 > 4. My independent fit gives dchi2 = 3.8 at 0.14 and 8.6 at 0.10 Msun, so the bound is about 0.13-0.14 Msun. The difference comes from the dilution prescription and the detrending. If the optional sentence 'the transit shape disfavours hosts below 0.17 Msun' is added, it should read '~0.15 Msun'.

*Evidence:* verify/v3_fit_A.json hosts M = 0.10, 0.14, 0.17 dchi2_free = 8.56, 3.77, 1.85

### [NOTE] Companion-space inner edges (RV about 2-3 au, astrometry about 20 au) are crude and must not be quoted

The auditor's astrometric exclusion uses a heuristic photocentre rms (0.58 alpha, times 1.47 (T/P)^2 for P > T) compared with the total observed AEN, with no statistical margin, even though AEN for G ~ 10.5 stars is dominated by calibration noise (median AEN_sig about 16). The RV edge is where P reaches the 923 d baseline, and there the sinusoid-based K limit extrapolates poorly. Neither number is in the manuscript and neither should be added without a scanning-law or MOLUSC-style simulation.

*Evidence:* auditor c06_space.py lines 60-70; verify/v4_ruwe_pct.json median_aen_sig = 15.9
