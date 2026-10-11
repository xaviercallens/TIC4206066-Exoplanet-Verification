# Verification of Audit 06 (stellar parameters)

Scripts: v01_ariadne.py (BMA.pkl re-read), v02_overlum.py (Way+2026 catalogue, own scatter and mass limits from the raw Pecaut & Mamajek EEM table), v03_relations.py (own Monte Carlo: Lindegren ZP, M_K, Mann+2015 R, Mann+2019 M, rho, rho=3.5 threshold). Coefficients checked in the arXiv sources (m15/main.tex l.1684; m19/Masses_1.tex l.753-857; zp=7.5).
Way+2026 text: arXiv HTML 2602.20096 fetched 2026-09-25 12:17 UTC, sha256 6a16843a...f1f6 (way2026_arxiv.html).

## Results
- Distance 35.619 (raw) / 35.577 pc (ZP -0.0332 mas); M_K 4.938 +-0.029; R 0.6229 +-0.019; M 0.6062 +-0.015; rho 2.507 +0.243/-0.217; M_G 7.757. All reproduce the manuscript. Not circular: inputs are fresh Gaia/2MASS, coefficients from the primary sources.
- F2 (ARIADNE errors) CONFIRMED, MINOR: ARIADNE 'weighted_average' = sum_i w_i * choice(post_i) (fitter.py l.1638) -> narrowed. Mixture from BMA.pkl: R 0.620 +0.020/-0.018, Teff 4017 +58/-64 (pickle 'weighted_samples' KDE agrees). Caveat: best_fit_sample.dat errors (+123/-78 K) do not match simple percentiles of the pickle samples; do not quote them as 'the mixture'.
- F1 (+-0.07 attributed to Way+2026) CONFIRMED, MINOR: no 0.07 or per-star error in the paper text; robust local scatter of O = 0.065-0.074 (1.4826 MAD) for several windows.
- F3 (0.38 Msun) CONFIRMED number: 0.380 (sigma 0.068), 0.349 (0.052), 0.420 (0.10); empirical lower-side 2.3% tail of O (-0.127 mag) gives ~0.37. Local median O is +0.014, so the target is slightly below the local median (conservative).
- Host threshold rho=3.5: M_K 5.446, M 0.529 Msun (manuscript 0.54): rounding, NOTE.

## Issue the auditor missed
- The auditor's proposed wording '(<~20 au)' and 'blended <~20 au ~ 0.6"' misreads Way+2026: their '<~20 AU' is the Gaia 0.2" resolution limit at 100 pc. At 35.6 pc, 0.2" ~ 7 au. The fix should read 'unresolved by Gaia (<~0.2", ~7 au at 35.6 pc)' or be stated in arcsec. MINOR (the fix itself would introduce an error).
- W1 = 7.61 is brighter than the WISE W1 saturation onset (~8.1 mag); W1 (and marginally W2) are PSF-wing photometry of a partially saturated star. No effect on the no-excess conclusion (W3/W4 decisive). NOTE.
