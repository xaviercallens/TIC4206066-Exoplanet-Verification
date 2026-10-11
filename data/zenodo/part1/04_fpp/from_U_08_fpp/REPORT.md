# Audit 08_fpp — REPORT

> Rendered verbatim from the agent's structured result (`_structured/audit__08_fpp.json`); the sub-agent harness blocked agents from writing .md report files. Scripts, outputs and logs are in this folder.


**Verdict:** PASS_WITH_CHANGES

## Summary

I reproduced the FPP numbers from scratch, and they agree with the manuscript. Two interpretive statements are wrong or misleading and must be reworded before submission.

**What I did.** I used only the fresh 00_data files:
- S32 SPOC 2-min PDCSAP (sha256 764dad34...) and S98 QLP SAP_FLUX (sha256 a0341366...), keeping QUALITY==0 cadences. This drops the S98 straylight bit 2048 and the QLP bits 2^29 and 2^30.
- Detrended with wotan biweight (0.75 d window, both signals masked), removed the other signal, folded within +-0.25 d and binned to 6 min.
- sigma = empirical out-of-transit scatter: A 113 ppm, B 232 ppm. Depths: A 487 ppm, B 1037 ppm.

TRICERATOPS 1.1.0 is byte-identical to the PyPI wheel (sha256 94c090b1..., released 2026-09-02). Setup:
- Sectors [32,98]; the SPOC S32 optimal aperture taken from the TPF (18 px); a 3x3 aperture for S98 (QLP header BESTAP='3x3 square'; QLP uses TGLC PSF photometry from S94 on).
- Star: M 0.606, R 0.623, Teff 4030, parallax 28.075 mas.
- Network steps (TIC query, TESScut, Gaia DR3 background) all worked. The target object was pickled so every run saw the same field.
- My own driver records the evidence of each STP/PTP draw.

**Results (N=1e6).**
- Baseline, 20 runs each: A FPP 0.0342+-0.0035 (range 0.029-0.040); B FPP 0.200+-0.012 (0.172-0.223).
- Sensitivity runs (5 each: TIC stellar parameters, Teff 4012, 5x5 apertures, T0 shifts, 10-min bins, formal sigma, single-sector folds): A 0.032-0.043, B 0.14-0.21.
- NFPP <= 1e-7 in every variant (1.5e-16 for A in the baseline). STP is >=99.8% of the FPP.
- Fresh MOLUSC run (github woodml/MOLUSC b26efbf; needed a one-line CLI patch and scipy<1.14): 429797/1e6 companions survive and 129548 have periastron > 10 AU. MOLUSC-mode FPP: A 0.0055+-0.0036, B 0.069+-0.016 (10 runs each). All of these match the manuscript within scatter.

**Finding 1 (MAJOR): the explanation of the MOLUSC drop is wrong.** Both code facts the manuscript cites are true: rows with a(1-e) > 10 AU are kept at marginal_likelihoods.py l.466/684/948/1176, and lnprior_companion is set to 0 at l.521/750/1038/1284. But they are not why the FPP falls.
- The effective companion prior in MOLUSC mode is n_kept/N = 0.130, against the default Moe & Di Stefano value within 2.2" (78 AU) of 0.141.
- Dropping the occurrence prior amounts to assuming every star has a companion, which pushes the FPP up, not down.
- The real cause: the Gaia RUWE and contrast limits cut the prior mass in companions above 0.38 Msun from 0.046 to 0.0083 (5.5x). Only 6.4% of kept MOLUSC companions are above 0.38 Msun, against 32% in the default draws. The FPP ratio (about 6) matches.

**Finding 2 (MAJOR): the headline FPP contradicts the manuscript's own constraints.** For signal A:
- 94% of the STP evidence comes from companion hosts above 0.38 Msun, and 62% from near-twins of 0.54-0.61 Msun.
- 51% comes from orbits with e < 0.1.
- The manuscript itself excludes both: companions above 0.38 Msun by overluminosity, and circular orbits unless the host is above 0.54 Msun.
- Keeping only companions of 0.38 Msun or less in the same draws gives FPP_A = 0.002-0.007; a 0.45 Msun cut gives 0.005-0.012.
- So "without it neither signal meets FPP<0.015" is misleading for A. For B it still holds (0.10-0.17 after the cut). Declining to claim validation remains correct.

**Minor points.**
- The A fold range "0.031-0.034" understates the sensitivity; my variants span 0.032-0.043.
- The validation criterion should also cite NFPP < 1e-3, and note that the thresholds were calibrated with the 2021 code (circular orbits, TRILEGAL).
- The MOLUSC-mode value scales as 1/N: 0.0125, 0.0055 and 0.0022 for N = 5e5, 1e6 and 2e6.

**Checked and OK.**
- "FPP presumes a real signal".
- The FPP is never turned into a planet probability.
- 2.2" is the TRICERATOPS default (hard-coded, from Brown et al. 2018).
- The Guerrero paraphrase is accurate (their words: estimates "do not apply to planet candidates detected only by the QLP"). The stronger reason not to apply a boost is that B is not secure.

**Reproducibility note.** TRICERATOPS 1.1.0 MOLUSC mode crashes under pandas 3 ("assignment destination is read-only"). I used a shim that passes the same numbers as writable arrays.

**REPORT.md was not written.** The harness blocked writing it with "Subagents should return findings as text, not write report files". The full report text, tables and wording are in this output instead. Scripts, runs and outputs are in /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp/:
- f01_build_folds.py, f02_trice_run.py, f03_priors.py, f04_summarize.py, f05_plot.py
- out_fpp_summary.{json,csv,txt}, out_f03_priors.json, fig_fpp_audit08.png
- runs/*/runs.json, folds/, molusc/, docs/, logs/

## Numeric rows

| quantity | manuscript | recomputed | status | method / note |
|---|---|---|---|---|
| TRICERATOPS version | 1.1.0 | 1.1.0 (installed = PyPI wheel sha256 94c090b1...; diff -r empty) | MATCH | pip metadata + wheel diff  |
| Runs x draws (baseline) | 20 runs of 1e6 | 20 x 1e6 per signal reproduced | MATCH |   |
| Unresolved-companion radius | 2.2 arcsec | 2.2 arcsec is hard-coded default without contrast curve (marginal_likelihoods.py l.490,719,1007,1247); = 78 AU at 35.6 pc; f_comp=0.141 | MATCH | source reading  |
| FPP, 3.18 d signal (headline) | 0.031 +- 0.003 | 0.0342 +- 0.0035 (20x1e6, own S32+S98 fold; range 0.029-0.040) | MINOR_DIFF | TRICERATOPS 1.1.0 on own fold within fold systematics |
| FPP range across folds, 3.18 d | 0.031-0.034 across three folds | 0.032-0.043 across my variants (T0 +-0.0005 d 0.033/0.037; S98-only 0.032; S32-only 0.041; 10-min bins 0.043; TIC params 0.040; 5x5 ap 0.034; formal sigma 0.037) | MINOR_DIFF |  quoted range too narrow |
| FPP, 11.13 d signal across folds and stellar params | 0.13-0.20 | 0.14-0.21 (baseline 0.200+-0.012; TIC params 0.141; T0+0.003 d 0.149; T0-0.003 d 0.211; 5x5 0.213) | MATCH |   |
| Share of FPP from bound-companion (STP) scenario | almost all | STP/FPP >= 0.998 in all runs (A: STP 0.034; B: 0.199) | MATCH |   |
| NFPP | not quoted | A 1.5e-16, B 0 (baseline); <= 1e-7 in all variants | NOT_CHECKED |  needed for the full Giacalone criterion; it passes |
| FPP with MOLUSC, 3.18 d | 0.005 +- 0.003 | 0.0055 +- 0.0036 (10x1e6, fresh MOLUSC 5 Gyr, RUWE+Gaia contrast; range 0.0002-0.013) | MATCH |  scales as 1/N: 0.0125 (N=5e5), 0.0022 (N=2e6) |
| FPP with MOLUSC, 11.13 d | 0.06 +- 0.02 | 0.069 +- 0.016 (10x1e6) | MATCH |   |
| MOLUSC-mode periastron cut | periastron beyond 10 AU | a(1-e) > 10 AU (l.466/684/948/1176); 129548 of 429797 fresh survivors kept | MATCH |   |
| MOLUSC-mode effective companion prior vs default | 'drops the companion-occurrence prior' (implied cause of lower FPP) | default f_comp=0.141; MOLUSC n_kept/N=0.130; prior mass >0.38 Msun 0.046 vs 0.0083 | DISCREPANT |  the causal claim is wrong; the drop comes from Gaia removing massive companions |
| Validation criterion | FPP < 0.015 | Giacalone+2021: FPP < 0.015 AND NFPP < 1e-3 (likely planet: FPP < 0.5, NFPP < 1e-3) | MINOR_DIFF |  incomplete quote |
| FPP_A with manuscript's own companion constraint (M <= 0.38 Msun) | not given ('neither signal meets FPP<0.015') | 0.0021 baseline; 0.002-0.007 across variants (M <= 0.45: 0.005-0.012) | DISCREPANT |  contradicts 'neither signal meets' for A |
| FPP_B with M <= 0.38 Msun | not given | 0.157 baseline; 0.10-0.17 | NOT_CHECKED |  B stays above threshold |
| Guerrero multiplicity-boost statement | does not apply to candidates found only by QLP | Guerrero+2021: estimates 'do not apply to planet candidates detected only by the QLP'; boost derived from Prime Mission SPOC 2-min TOIs | MATCH |  accurate paraphrase; a better reason is that B is not secure |

## Findings

### F1 [MAJOR] The stated reason MOLUSC lowers the FPP is incorrect

The manuscript says MOLUSC mode lowers the FPP 'only because this mode keeps companions with periastron beyond 10 AU alone and drops the companion-occurrence prior'. Both code facts are true. They do not explain the drop.
- The effective companion prior in MOLUSC mode is n_kept/N = 0.130, almost the same as the default Moe & Di Stefano prior within 2.2 arcsec (0.141).
- Setting lnprior_companion = 0 amounts to assuming every star has a companion, which pushes the FPP up.
- The actual cause: the Gaia RUWE and contrast limits remove massive companions, which dominate the STP evidence. The prior mass above 0.38 Msun falls from 0.046 to 0.0083 (5.5x), and the FPP falls about 6x.

*Evidence:* - TRICERATOPS source: trice_env/.../triceratops/marginal_likelihoods.py l.466/684/948/1176 (a(1-e)>10 filter), l.468 (clamp to 0.1 Msun), l.469 (zero-pad to N), l.521/750/1038/1284 (lnprior_companion = zeros); priors.py l.580-782 (lnprior_bound_TP).
- Prior comparison: 08_fpp/out_f03_priors.json. default_f_comp_within_2.2arcsec = 0.141; fresh MOLUSC effective_prior_at_N1e6 = 0.1295; frac_kept_mass_gt_0.38 = 0.064 vs default 0.323; effective prior above 0.38 Msun 0.0083 vs 0.0456.
- FPP runs: A default 0.0342 vs MOLUSC 0.0055 (08_fpp/out_fpp_summary.txt). In MOLUSC mode the FPP scales as 1/N: 0.0125, 0.0055 and 0.0022 for N = 5e5, 1e6 and 2e6.

*Location:* Sec. 4, paragraph 2, sentence 2

*Proposed fix:* Replace the clause with: 'mainly because the Gaia limits remove companions above ~0.4 Msun, which dominate the default FPP; this mode also assumes a companion for every star and keeps only periastra > 10 AU'. Alternatively use the full wording in recommended_wording.

### F2 [MAJOR] Unconstrained FPP is dominated by companions the manuscript itself excludes; 'neither signal meets FPP<0.015' is misleading for the 3.18 d signal

The manuscript uses Gaia XP overluminosity to exclude unresolved companions above 0.38 Msun (2 sigma). It also says circular orbits need hosts above 0.54 Msun, so the scenario survives only for eccentric orbits. The TRICERATOPS default FPP_A = 0.03 comes almost entirely from exactly those excluded configurations. Applying the manuscript's own 0.38 Msun cut to the same TRICERATOPS draws gives FPP_A = 0.002-0.007, below 0.015. So quoting 0.031 as the key number, and saying neither signal meets the criterion without MOLUSC, is internally inconsistent for A. For B, 71% of the STP evidence comes from hosts of 0.38 Msun or less, and the FPP stays at 0.10-0.17 after the cut. Not claiming validation is still defensible: the cut is post hoc, the 0.38 Msun limit is a 2 sigma bound, and TRICERATOPS thresholds were calibrated on the default method. But the stated reasoning has to change.

*Evidence:* Evidence-weighted STP draws recorded by hooks in 08_fpp/f02_trice_run.py; runs/A_base/runs.json and runs/B_base/runs.json (20 runs each); out_fpp_summary.json.
- A: STP share by host mass (Msun) = 0.10-0.20: 0.003, 0.20-0.30: 0.017, 0.30-0.38: 0.037, 0.38-0.45: 0.077, 0.45-0.54: 0.251, 0.54-0.61: 0.616.
- A: share with e < 0.1 = 0.51; share with M <= 0.38 and e < 0.1 = 1e-9.
- FPP_comp_mass_le_0.38: A_base 0.0021, variants 0.002-0.007. FPP_comp_mass_le_0.45: 0.005-0.012.
- B: STP share above 0.38 Msun = 0.29; FPP(M <= 0.38) = 0.157.
- Plot: 08_fpp/fig_fpp_audit08.png.

*Location:* Abstract ('TRICERATOPS gives FPP = 0.031 +- 0.003 under its default assumptions'); Sec. 4 paragraph 2 ('without it neither signal meets ...')

*Proposed fix:* State that about 95% of the default FPP of the 3.18 d signal comes from companions above 0.38 Msun that the overluminosity limit disfavors. Give the constrained value (~0.002-0.007) alongside the MOLUSC value. Explain that validation is not claimed because these are post hoc constraints outside the calibrated TRICERATOPS configuration, and because high-resolution imaging is needed. Exact wording is in recommended_wording.

### F3 [MINOR] Fold/method sensitivity of FPP_A understated

The manuscript gives 0.031 +- 0.003 and '0.031-0.034 across three folds'. My independent fold gives 0.034 +- 0.004. Reasonable method variants (binning, single sector, stellar parameters) span 0.032-0.043. The headline is fine within systematics, but the quoted range is too narrow.

*Evidence:* 08_fpp/out_fpp_summary.txt:
- A_base 0.0342+-0.0035 (20 runs); A_10min 0.0431+-0.0024; A_S32only 0.0409+-0.0028; A_tic 0.0398+-0.0053; A_T0p 0.0366; A_T0m 0.0327; A_S98only 0.0322; A_teff4012 0.0335; A_ap5x5 0.0344; A_sigformal 0.0367.
- Earlier package values for comparison: 0.0312 and 0.0343 (tic4206066_PUBLICATION_AUDIT/out/fpp/fpp_table.md; rnaas_tic4206066_pubprep/analysis/out/fpp/fpp_ownfold_A.json).

*Location:* Sec. 4 paragraph 2; abstract

*Proposed fix:* Quote 'FPP = 0.03-0.04' (or 0.031 +- 0.003 for the adopted fold, 0.03-0.04 across folds and binning).

### F4 [MINOR] Validation criterion quoted incompletely; calibration caveat missing

Giacalone et al. (2021) define validated planets as FPP < 0.015 AND NFPP < 1e-3. The manuscript cites only FPP < 0.015; NFPP is about 1e-16 here, so this is harmless. The criterion was also calibrated with the 2021 code. That code assumed zero eccentricity and a TRILEGAL background. TRICERATOPS 1.1.0 samples eccentric orbits and uses a Gaia DR3 background by default. Its own KNOWN_ISSUES.rst says this default is not yet benchmarked against the paper sample.

*Evidence:* 08_fpp/docs/G21.txt (pdftotext of arXiv:2002.00691):
- l.81: 'validated planets (FPP < 0.015 and NFPP < 10^-3)'.
- l.657-676: 'assume a fixed orbital period and zero eccentricity'.
TRICERATOPS 1.1.0 source and docs:
- priors.py l.148: beta.rvs(0.867, 3.030) for planet eccentricities.
- docs/triceratops_CHANGES.rst: Gaia DR3 background is the default in 1.1.0.
- docs/triceratops_KNOWN_ISSUES.rst: 'Validate the Gaia DR3 background-population default ... Still to do'.

*Location:* Sec. 4 paragraph 2

*Proposed fix:* Write 'FPP < 0.015 and NFPP < 10^-3' and, if space allows, note that the criterion was calibrated on the original code configuration.

### F5 [MINOR] MOLUSC-mode FPP is an N-dependent, high-scatter number

In MOLUSC mode TRICERATOPS pads the kept companions with q=0 up to N and averages over N, so the effective prior is n_kept/N. The FPP therefore scales as 1/N. The quoted value is meaningful only because N equals MOLUSC's n_generated (1e6), which implicitly assumes a 100% companion fraction. Run-to-run scatter is large. 63% of kept companions are brown dwarfs that are clamped to 0.1 Msun.

*Evidence:* - runs/A_moluscFresh5_N5e5: 0.0125 +- 0.0135.
- runs/A_moluscFresh5 (N=1e6): 0.0055 +- 0.0036, range 0.0002-0.0129.
- runs/A_moluscFresh5_N2e6: 0.0022 +- 0.0021.
- out_f03_priors.json: frac_kept_q_below_0.1Msun_clamped = 0.626.
- Code: marginal_likelihoods.py l.468-469.

*Location:* Sec. 4 paragraph 2

*Proposed fix:* If the MOLUSC value is kept, state that N matched the MOLUSC sample and that the mode assumes a companion for every star.

### F6 [NOTE] Reproducibility: MOLUSC mode crashes under pandas 3; MOLUSC needs patch and scipy<1.14

- TRICERATOPS 1.1.0 with molusc_file raises 'ValueError: assignment destination is read-only' at marginal_likelihoods.py:468 under pandas 3.0.6 (copy-on-write). I worked around it with a read_csv shim that returns writable copies of the same numbers (f02_trice_run.py).
- MOLUSC b26efbf's CLI raises AttributeError 'star_jitter' without a one-line patch (molusc/audit08_molusc.diff).
- MOLUSC calls scipy.interpolate.interp2d, which was removed in scipy 1.14.
- The package should record the pandas and scipy versions and the MOLUSC patch.

*Evidence:* 08_fpp/logs/queue2.log (traceback); /private/tmp/.../scratchpad/a08_cache/molusc_run and 08_fpp/molusc/audit08_molusc.diff; trice_env pandas 3.0.6, scipy 1.13.1.

*Location:* reproducibility package

*Proposed fix:* Record the software versions and patches in the data/software availability notes.

### F7 [NOTE] Multiplicity-boost sentence accurate but the best justification is omitted

The paraphrase of Guerrero et al. (2021) is accurate. Their boost (about 20) was derived from Prime Mission SPOC 2-min TOIs. These signals are not TOIs of any pipeline; they were found by the author in QLP light curves, and SPOC 2-min S32 data exist without a TCE. The strongest reason not to apply a boost is that the 11.13 d signal (about 3 sigma) does not establish the system as multi-transiting. For scale: applying a boost of about 20 would put FPP_A near 0.002, so leaving it out is the conservative choice.

*Evidence:* 08_fpp/docs/Gu21.txt: l.1575-1577 'We consider only TOIs found by the SPOC pipeline'; l.1634 '(75/930)/(930/232705) ~ 20'; l.1639-1641 'these estimates for the multiplicity boost do not apply to planet candidates detected only by the QLP'.

*Location:* Sec. 4 paragraph 2, last sentence

*Proposed fix:* 'No multiplicity boost is applied, because the 11.13 d signal is not secure (cf. Guerrero et al. 2021).'

### F8 [NOTE] Checked OK: 'FPP presumes a real signal', no planet-probability conversion, 2.2 arcsec default, STP dominance

- TRICERATOPS tests only astrophysical scenarios. Giacalone et al. discard instrumental false alarms, 'which our tool does not test for'.
- The text never turns FPP into a planet probability.
- 2.2 arcsec is the default when there is no contrast curve (Brown et al. 2018).
- STP accounts for at least 99.8% of the FPP.
- DTP counts as a planet and the Gaia background contributes about 1e-5, so the choice of background source is irrelevant here.

*Evidence:* docs/G21.txt l.1242 and l.690-691; triceratops.py l.1697; runs/A_base/runs.json probs: TP 0.836, PTP 0.126, STP 0.034, DTP 0.0036, BTP 1e-5.

*Location:* Sec. 4 paragraph 2

*Proposed fix:* 

### F9 [NOTE] For audit 09: Gaia DR3 RV statistics are borderline

Outside this audit's scope. The Gaia DR3 row for the target has rv_chisq_pvalue = 4.1e-5, rv_renormalised_gof = 3.0, rv_amplitude_robust = 3.2 km/s and rv_nb_transits = 33. The manuscript's 'radial velocities that exclude an eclipsing stellar companion' should be checked against the published Gaia RV-variability criteria by the companion audit.

*Evidence:* 00_data/catalogs/gaiadr3_gaia_source_cone240as.csv, source_id 3220388198192519424.

*Location:* Sec. 4 paragraph 1

*Proposed fix:* 

## Recommended wording

Sec. 4, paragraph 2 (about the same length):
"TRICERATOPS 1.1.0 \citep{Giacalone2021}, with Sector 32+98 folds, 20 runs of $10^6$ draws, and its default of bound companions unresolved within $2.2''$, gives FPP $=0.03$--$0.04$ for the 3.18\,d signal and 0.13--0.21 for the 11.13\,d signal across folds and stellar parameters (FPP presumes a real signal); almost all of it is the bound-companion scenario. For the 3.18\,d signal, $\approx$95\% of this comes from companions above $0.38\,M_\odot$, which the overluminosity limit disfavors; removing them, or adding the Gaia RUWE and contrast limits through MOLUSC \citep{Wood2021}, gives FPP $\approx0.002$--$0.007$. For the 11.13\,d signal the companion hosts are mostly below $0.38\,M_\odot$, and FPP remains 0.06--0.17. These constrained values rely on a post hoc mass cut or on MOLUSC-mode settings (a companion assumed for every star, periastron $>10$\,AU) outside the configuration for which \citet{Giacalone2021} calibrated FPP $<0.015$ and NFPP $<10^{-3}$, so I do not claim validation. No multiplicity boost is applied, because the 11.13\,d signal is not secure \citep[cf.][]{Guerrero2021}."

Abstract: replace "TRICERATOPS gives FPP $=0.031\pm0.003$ under its default assumptions" with "TRICERATOPS gives FPP $\approx0.03$ under its default priors, almost all of it from bound companions that Gaia disfavors."

Minimum edit if the word budget is binding:
- Replace "but only because this mode keeps companions with periastron beyond 10\,AU alone and drops the companion-occurrence prior" with "mainly because the Gaia limits remove companions above $\sim$0.4\,$M_\odot$, which dominate the default FPP; this mode also assumes a companion for every star".
- Change "FPP $<0.015$" to "FPP $<0.015$ and NFPP $<10^{-3}$".

## Files

- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp/f01_build_folds.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp/f02_trice_run.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp/f03_priors.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp/f04_summarize.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp/f05_plot.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp/out_fpp_summary.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp/out_fpp_summary.csv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp/out_fpp_summary.txt
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp/out_f03_priors.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp/fig_fpp_audit08.png
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp/folds/folds_meta.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp/runs/
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp/molusc/molusc_fresh_age5_kept.csv.gz
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp/molusc/audit08_molusc.diff
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp/docs/
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp/logs/
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp/run_queue1.sh
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp/4206066_gaia_background.csv