# Audit 09_companion — REPORT

> Rendered verbatim from the agent's structured result (`_structured/audit__09_companion.json`); the sub-agent harness blocked agents from writing .md report files. Scripts, outputs and logs are in this folder.


**Verdict:** PASS_WITH_CHANGES

## Summary

Verdict: the Sec. 4 companion sentences hold up numerically. Three small wording changes are needed. They are mainly about scope and about TRICERATOPS. None of the problems is FATAL or MAJOR.

Figures checked with fresh Gaia DR3 data (TAP, 2026-09-25):
- RUWE is 1.214. There is no NSS solution (non_single_star = 0, all NSS tables empty).
- ipd_frac_multi_peak = 0, ipd_gof_harmonic_amplitude = 0.022, duplicated_source = False.
- astrometric_excess_noise = 0.093 mas (sig 15.8).
- No Gaia source lies within 10″; the nearest is at 11.06″ with G = 18.75. Within 10″, 2MASS and AllWISE each have one source, and it is the target.

Radial velocities (from the rv_method_used=1 error formula):
- The epoch scatter σ(V_j) is 0.877 km/s. rv_chisq_pvalue is 4e-5, but rv_renormalised_gof = 3.0 is below the Katz >4 variability cut. This means the formal epoch errors are too small for this magnitude/temperature bin.
- In a Monte Carlo of Gaia-like sampling, 7.8% of constant stars scatter more than this star. rv_amplitude_robust = 3.2 km/s is just the noise peak-to-peak.
- The 95% limit is K < 1.16 km/s, which gives M sin i < 0.0058 Msun at 3.18 d and < 0.0088 Msun at 11.13 d.
- The package's figures (K < 3.2 km/s → 0.016 and 0.024 Msun) reproduce arithmetically but are about 3x too loose.
- An eclipsing stellar companion of the target is excluded: a 0.08 Msun companion would give K = 14.9 km/s at 3.18 d and 11.8 km/s at 6.37 d.

Overluminosity: re-deriving the conversion gives 0.374–0.378 Msun at 2σ (manuscript 0.38). I did not re-download the Way et al. catalogue itself.

Transit shape vs host (own batman profile fits of S6+S32+S98, with dilution included):
- With the target as host, the upper density limit (Δχ²=4) is ρ < 3.53 ρsun. This matches the manuscript's 3.5 (97.7%).
- Signal A, circular orbit on a companion host: needs M_c > 0.523 Msun at Δχ²=4, or > 0.506 at Δχ²=9. The manuscript's 0.54 is about right.
- Signal A, any orbit: hosts below 0.168 Msun are disfavoured (Δχ² > 4) because the strongly diluted transit becomes V-shaped.
- For hosts allowed by the overluminosity (0.17–0.37 Msun), the planet needs e ≥ 0.19–0.50 even in the most favourable orientation (ω = 270°).
- Signal B: the transit shape is unconstrained, and a circular orbit on any companion of at least 0.14 Msun fits (Δχ² < 1).

TRICERATOPS 1.1.0 source:
- The STP scenario is not circular-only. `sample_ecc` (priors.py l.134–155) draws e from Beta(0.867, 3.03), and ω is uniform. The transit model (likelihoods.py l.37–91) uses both.
- Companion masses run 0.10–0.606 Msun (`sample_q_companion`); 32% of draws are above 0.38 Msun.
- The bound-companion prior (`lnprior_bound_TP`) counts only P > 2500 d out to 2.2″; its value here is 0.141.
- So the default FPP already includes eccentric orbits. What it does not include is the overluminosity cut, which removes the best-fitting draws (0.52–0.61 Msun, near-circular).
- My order-of-magnitude estimate is that roughly 10% of the STP evidence survives the Gaia constraints. This is an estimate, not a recomputed FPP.

Remaining companion space:
- Signal A: 0.17–0.37 Msun, e ≳ 0.2–0.5. RVs rule out companions inside about 2–3 au. A crude astrometric model suggests nothing inside about 20 au, but this is not rigorous. The upper bound is the 2.2″ (78 au) limit.
- Signal B: 0.14–0.37 Msun, any eccentricity.

Not completed: the RUWE percentile among similar stars. Both the Gaia archive and TAPVizieR timed out. The scripts `c01c_ruwe_sample.py` and `c08_ruwe_pct.py` are ready to run.

REPORT.md could not be written: the harness blocks subagents from writing report files. Its full content is in this structured result.

Files are in /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/09_companion/:
- scripts: c00_relations.py, c01_gaia_fetch.py, c01b_vizier.py, c01c_ruwe_sample.py, c02_lc_prep.py, c03_transit_model.py, c04_profile.py, c05_rv.py, c06_space.py, c07_trice_prior.py, c08_ruwe_pct.py
- data/ (with PROVENANCE.jsonl)
- out/: profile_A.json, profile_B.json, rv_limits.json, companion_space.json, companion_space_table.csv, triceratops_priors.json, lc_windows_meta.json
- plots/companion_space.png
- logs/

## Numeric rows

| quantity | manuscript | recomputed | status | method / note |
|---|---|---|---|---|
| Gaia DR3 RUWE | 1.21 | 1.2143 | MATCH | Gaia TAP gaia_source  |
| Gaia non-single-star solution | none | non_single_star=0; nss_two_body_orbit, nss_acceleration_astro, nss_non_linear_spectro, nss_vim_fl: 0 rows | MATCH | Gaia TAP  |
| RV excludes eclipsing stellar companion of target | excluded | 95% K<1.16 km/s (P<=1000 d); a 0.08 Msun companion gives K=14.9 km/s (3.18 d) and 11.8 km/s (6.37 d) | MATCH | sigma(V_j) from radial_velocity_error formula + Monte Carlo of Gaia sampling (c05_rv.py)  |
| K limit (package, not in manuscript text) | 3.2 km/s | rv_amplitude_robust=3.225 (peak-to-peak noise); calibrated 95% K<1.16 km/s | MINOR_DIFF | c05_rv.py 3.2 is conservative, not a K measurement |
| M sin i limit at 3.18 d (package) | <0.016 Msun | 0.0161 for K=3.2; 0.0058 for calibrated 95% K | MINOR_DIFF | mass function, M1=0.606  |
| Gaia XP overluminosity | +0.003 +- 0.07 mag | not re-downloaded; conversion re-derived with O=+0.003, sigma=0.068 | NOT_CHECKED |   |
| Overluminosity companion mass limit (2 sigma) | 0.38 Msun | 0.374 (sigma 0.068) to 0.378 (sigma 0.07) Msun; 3 sigma 0.43 | MATCH | flux ratio 10^(0.4 O_lim)-1 -> M_G -> Mann19/Mamajek chain (c06_space.py)  |
| Minimum circular-orbit host mass (signal A) | 0.54 Msun | 0.523 Msun (dchi2=4), 0.506 (dchi2=9), with dilution | MINOR_DIFF | batman profile fit with dilution, Mann radii, Claret LD (c04_profile.py)  |
| rho* upper limit (density free) | <3.5 rho_sun (97.7%) | 3.53 rho_sun (profile dchi2=4); best fit 3.13 | MATCH | c04_profile.py  |
| T14 of signal A | 1.87 +- 0.05 h | 1.866 h (free rho), 1.868 h (Mann rho) | MATCH | c04_profile.py  |
| M_Ks | 4.94 | 4.935 (2MASS Ks 7.694, plx 28.075) | MATCH |   |
| STP scenario survives only for eccentric orbits | eccentric only | true for signal A (e>=0.19-0.50 for 0.17-0.37 Msun hosts); false for signal B (circular allowed for all hosts >=0.14 Msun) | DISCREPANT | c04_profile.py A and B scope: holds for 3.18 d only |
| Almost all FPP is STP (TRICERATOPS draws) | almost all | not re-run; source shows STP draws eccentric orbits and companions 0.10-0.606 Msun, 32% of draws >0.38 Msun | NOT_CHECKED |   |

## Findings

### F1 [MINOR] Circular-orbit host-mass argument holds only for the 3.18 d signal

In Sec. 4 the sentence 'A planet on a bound companion in a circular orbit would need a host above 0.54 Msun ... survives only for eccentric orbits' comes right after the 11.13 d localisation, so it reads as covering both signals. For signal B the transit shape is unconstrained, and a circular orbit around any companion of at least 0.14 Msun fits the data. The recomputed minimum host mass for signal A is 0.52 Msun at 2σ, so 0.54 is about right.

*Evidence:* out/profile_A.json: circular dchi2 = 58.8 at 0.38 Msun, 10.9 at 0.50, 4.7 at 0.52, 0.5 at 0.54; 2σ crossing at 0.523 Msun (out/companion_space.json). out/profile_B.json: circular dchi2 0.4-0.8 for 0.16-0.606 Msun.

*Location:* Sec. 4, first paragraph, last sentence

*Proposed fix:* Write 'For the 3.18 d signal, a planet on such a companion would need a host above ~0.52-0.54 Msun in a circular orbit ...' and add 'the 11.13 d transit shape does not constrain the host.'

### F2 [MINOR] TRICERATOPS STP already includes eccentric orbits; 'untested' framing is imprecise

The abstract calls the eccentric bound-companion case 'the main untested false-positive scenario'. In TRICERATOPS 1.1.0 the STP scenario draws planet eccentricities from Beta(0.867, 3.03) with uniform ω and models them in the light curve. So FPP = 0.031 already contains eccentric STP. What the default FPP leaves out is the overluminosity cut on companions above 0.38 Msun, and those companions supply the best-fitting, near-circular STP draws. The surviving channel is 0.17-0.37 Msun hosts with e ≥ 0.19-0.50. Under TRICERATOPS' own eccentricity prior only 2-12% of draws satisfy this.

*Evidence:* trice_env site-packages/triceratops/marginal_likelihoods.py lnZ_STP l.892-1113 (eccs = sample_ecc(..., planet=True) l.1043); priors.py sample_ecc l.134-155; likelihoods.py simulate_TP_transit l.37-91 (evaluate_ps e=ecc, w=90-argp); out/triceratops_priors.json: 32% of companion draws >0.38 Msun; out/companion_space.json ecc_prior_fraction.

*Location:* Abstract; Sec. 4

*Proposed fix:* Say the default FPP includes companions above 0.38 Msun that the Gaia overluminosity excludes. Say the surviving STP needs a 0.17-0.37 Msun companion and e ≳ 0.2-0.5, which is disfavoured for a small planet at 3.18 d but only high-resolution imaging can close it.

### F3 [MINOR] Gaia RV: 3.2 km/s is peak-to-peak noise, not a K limit; low p-value needs context if quoted

The RV sentence in the manuscript is correct. The supporting package, however, treats rv_amplitude_robust (3.2 km/s, the max-min of 33 noisy epochs) as a K limit. The calibrated 95% limit is K < 1.16 km/s. Separately, rv_chisq_pvalue = 4e-5 could look like RV variability to a referee. It comes from formal epoch errors that are about 1.3x too small for this bin: rv_renormalised_gof = 3.0 is below the Katz >4 cut, and 7.8% of simulated constant stars scatter more than this star.

*Evidence:* out/rv_limits.json: sigma_Vj 0.877 km/s, typical 0.749, RUWE_RV 1.17, K0_sim frac_sd_ge_obs 0.078, frac_amp_ge_obs 0.39, K95 1.16 km/s -> Msini 0.0058 (3.18 d), 0.0088 (11.13 d). Gaia DR3 datamodel definitions (radial_velocity_error, rv_amplitude_robust, rv_renormalised_gof).

*Location:* Sec. 4; METHOD_AUDIT sec. 7

*Proposed fix:* If a number is quoted, use 'K ≲ 1.2 km/s (95%, 33 epochs)'. Update the package notes.

### F4 [NOTE] Overluminosity mass limit reproduced

Converting O = +0.003 with a 2σ width of 0.068-0.07 mag gives a G-band flux ratio of 0.137-0.141, a companion M_G ≥ 9.89-9.92, and M ≤ 0.374-0.378 Msun. The Way catalogue row itself was not re-downloaded.

*Evidence:* out/companion_space.json overlum_2sig_*

*Location:* Sec. 4

*Proposed fix:* None.

### F5 [NOTE] rho* < 3.5 rho_sun confirmed

With the target as host, the profile-likelihood upper limit (Δχ² = 4) is 3.53 ρsun; the best fit is 3.13 with T14 = 1.866 h.

*Evidence:* out/profile_A.json target_free_rho; companions M=0.606 rho_app_max_2sig 3.53

*Location:* Sec. 2

*Proposed fix:* None.

### F6 [NOTE] Transit shape also sets a lower host-mass bound (signal A)

A very faint host needs a large, heavily diluted planet, which produces a V-shaped transit. For signal A this disfavours hosts below 0.168 Msun at Δχ² > 4 and below 0.124 Msun at Δχ² > 9. Combined with the overluminosity cut, the remaining host range is 0.17-0.37 Msun. RVs exclude such companions inside about 2-3 au.

*Evidence:* out/profile_A.json dchi2_ecc: 13.7 (0.10), 9.6 (0.12), 4.6 (0.16), 3.1 (0.18); out/companion_space_table.csv; plots/companion_space.png

*Location:* Sec. 4

*Proposed fix:* Optionally add 'the transit shape disfavours hosts below 0.17 Msun'.

### F7 [NOTE] Audit limitations

Several parts of this audit are incomplete or approximate. (1) The RUWE and excess-noise percentiles among similar stars were not computed: the Gaia TAP and TAPVizieR comparison queries timed out. (2) The astrometric limit comes from a crude photocentre model, not a scanning-law simulation or the MOLUSC/Kraus RUWE calibration. (3) The fits assume white noise (χ²/dof = 1.04). (4) Claret limb darkening is clamped at 3500 K for cooler hosts. (5) TESS-band fluxes of late-M companions use Stassun Eq. 1 extrapolated beyond its validity range. (6) The harness blocked writing REPORT.md, so the report content is in this structured result.

*Evidence:* logs/c01_ruwe.log, the bmy1sce92 task output (TAPVizieR 504 Gateway Time-out), c06_space.py docstring

*Location:* n/a

*Proposed fix:* Run c01c_ruwe_sample.py and c08_ruwe_pct.py when the archive responds.

## Recommended wording

Gaia DR3 gives RUWE = 1.21, no non-single-star solution, and radial velocities (33 epochs; K ≲ 1.2 km/s) that exclude an eclipsing stellar companion of the target. The Gaia XP overluminosity of Way et al. (2026), +0.003 ± 0.07 mag, excludes unresolved companions above 0.38 M☉ (2σ). For the 3.18 d signal, a planet on a bound companion in a circular orbit would need a host above ≈0.52 M☉, and the transit shape disfavours hosts below 0.17 M☉. The scenario therefore survives only for 0.17–0.38 M☉ companions with e ≳ 0.2–0.5. The 11.13 d transit shape does not constrain the host. [FPP paragraph: add that the default TRICERATOPS STP scenario includes eccentric orbits and companions up to 0.61 M☉, and does not apply the overluminosity limit.]

## Files

- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/09_companion/c00_relations.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/09_companion/c01_gaia_fetch.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/09_companion/c01b_vizier.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/09_companion/c01c_ruwe_sample.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/09_companion/c02_lc_prep.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/09_companion/c03_transit_model.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/09_companion/c04_profile.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/09_companion/c05_rv.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/09_companion/c06_space.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/09_companion/c07_trice_prior.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/09_companion/c08_ruwe_pct.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/09_companion/data/PROVENANCE.jsonl
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/09_companion/out/profile_A.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/09_companion/out/profile_B.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/09_companion/out/rv_limits.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/09_companion/out/companion_space.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/09_companion/out/companion_space_table.csv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/09_companion/out/triceratops_priors.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/09_companion/out/lc_windows_meta.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/09_companion/plots/companion_space.png