# Hidden-companion check — TIC 4206066 (StKM 1-561) (EXPLORATORY ADDENDUM, 2026-09-25)

Question: is there an unresolved bound companion star bright enough to host the transits? This is the scenario (STP) that
carries almost all of the TRICERATOPS false-positive probability (`../deep_followup_tic4206066/fpp/FPP_REPORT.md`).
No frozen pipeline output, config, cohort label or selection rule was opened or changed. Separate env in the session
scratchpad (astroARIADNE 1.5.0, numpy 2.4.6, pandas 3.0.6); project `.venv` not used.

## 1. Overluminosity in Gaia DR3 (`overluminosity.py`, `out/overluminosity.json`, `plots/overluminosity.png`)

Catalogue: Way, Lépine, Gagné & Medan 2026, AJ (doi 10.3847/1538-3881/ae48ae, arXiv:2602.20096), Zenodo 18500082,
CC-BY-4.0. `overluminous_catalog.csv` md5 d27fb6e7…0731, `catalog_of_systems_100pc_lite.csv` md5 3f1eb18b…c6ccef666
(both match Zenodo; raw files kept out of the project, the two rows for this star are in `data/`). Method: a random-forest
model predicts the absolute magnitude M_G of a *single* K/M dwarf from its Gaia XP (RP) spectrum; overluminosity
O = predicted − measured M_G. Paper threshold O > 0.376 mag (half the equal-mass offset 2.5 log10 2).

* **FACT** — TIC 4206066 is in the catalogue: measured M_G 7.755, predicted 7.758 → **O = +0.003 mag**; measured
  G−RP 0.899 vs predicted 0.896. Flag `ol = False`. Systems catalogue: 1 component, not in BANYAN Σ young groups.
* **FACT** — 5 205 stars with the same predicted M_G (±0.25) and colour (±0.05): median O = +0.016, robust σ = 0.068
  mag; 58 % of them are more overluminous than our star.
* **INFERENCE** — limit on companion light in G: O < 0.139 (2σ) → flux ratio < 0.137 → companion M_G > 9.92 →
  **mass < 0.38 M☉, q < 0.63** (Pecaut & Mamajek dwarf table v2022.04.16). At 3σ: < 0.42 M☉, q < 0.70.
  Conversion uses the same flux-sum logic as the paper (companion light adds to the predicted single-star magnitude);
  a companion also reddens the XP spectrum, which makes the predicted single-star M_G fainter and O larger, so the
  simple conversion is on the conservative side.
* Scope (from the paper): most sensitive to near-equal-mass companions within ≲ 20 AU; companions with q ≲ 0.6 are
  not detectable this way. The method validates at 98 % recovery for SB2 (q ≈ 1) and 65 % for RV-variable binaries.

**Relation to the FPP result.** In the earlier TRICERATOPS runs the best-fit STP companion (column `M_s` of the STP row
in every `probs_runNN.csv`) is heavier than 0.38 M☉ in 80 % (3.18 d, variant A), 100 % (B), 70 % (C) and 57 % / 80 % /
35 % (11.13 d, A/B/C) of runs. Those companions are now excluded at 2σ. This is a statement about best-fit draws only;
the FPP itself integrates over all companion masses and was **not recomputed**.

## 2. SED fit (`sed_fit.py`, `sed_residuals.py`, `sed_out/`, `out/sed_residuals.json`, `plots/sed_residuals.png`)

astroARIADNE 1.5.0 (Vines & Jenkins 2022), Bayesian model averaging over the grids the package allowed for this star:
btsettl (weight 0.31), btnextgen (0.29), btcond (0.25), phoenix (0.15); kurucz / ck04 / bosz were skipped by the
package's own gate (its catalogue lookup returned no Gaia Teff, so the star is treated as < 4000 K). dynesty, 500 live
points, Av fixed to 0 (35.6 pc), distance prior Bailer-Jones 35.58 ± 0.02 pc, radius uniform. Photometry (10 bands,
`sed_out/mags.dat`): Tycho-2 B_T V_T, Gaia DR2 BP G RP, 2MASS J H Ks, WISE W1 W2. Removed before the fit: SDSS and
Pan-STARRS (saturated at G = 10.5: SDSS saturates at r ≈ 14, PS1 at g,r,i ≈ 13.5, z ≈ 13, y ≈ 12) and TESS T (derived
from Gaia, not independent).

| quantity | SED (BMA) | earlier value | source of earlier value |
|---|---|---|---|
| Teff | 4012 ± 43 K | 4010 / 4035 / 4182; adopted 4080 ± 100 | TIC 8.2 / LAMOST LRS / Gaia GSP-Spec |
| R★ | **0.621 ± 0.015 R☉** | **0.623 ± 0.019 R☉** | Mann et al. 2015 from M_Ks |
| L★ | 0.090 ± 0.006 L☉ | – | |
| [Fe/H] | −0.08 ± 0.18 | −0.27 | Gaia GSP-Spec |
| log g | 5.17 ± 0.35 | 4.63 | derived from Mann M, R (SED does not constrain log g) |

* **FACT** — the single-star model reproduces every band from 0.5 to 4.6 µm to within 0.05 mag (BMA medians, catalogue
  errors): BP −0.009, V_T −0.043 (±0.095), G +0.003, RP −0.018, J −0.005, H +0.030, Ks −0.038 (±0.029), W1 +0.001,
  W2 +0.050 (±0.020, 2.5σ; the largest deviation). Tycho B_T is +0.62 ± 0.40 mag (1.5σ; B_T = 13.4 is near the
  Tycho-2 faint limit). No trend with wavelength; the J, H, Ks mean residual is −0.004 mag. **No infrared excess.**
* **FACT** — the SED radius agrees with the Mann relation radius used everywhere else (0.621 vs 0.623 R☉); Teff agrees
  with TIC and LAMOST and is 0.7σ below the adopted 4080 ± 100 K.
* **INFERENCE** — a fully blended companion would make G − Ks redder than a single star: by +0.115 mag for M2V
  (0.44 M☉), +0.067 for M4V (0.23 M☉), +0.038 for M5V (0.16 M☉) (Pecaut & Mamajek M_G, M_Ks; primary M_G 7.76,
  M_Ks 4.94). The observed G − Ks residual is −0.041 ± 0.03 mag, i.e. *bluer* than the model. Companions ≳ M4V
  (≳ 0.23 M☉) are therefore disfavoured at ~2σ. This is not a rigorous limit: Teff and R are free in the fit and absorb
  part of any excess, and the fitted per-band noise terms are non-zero. §1 (q < 0.63, published method) is the limit to
  quote; §2 is consistent with it and extends the argument to somewhat lower masses.

Failure modes: (F1) the final NetCDF export raised `ValueError: unable to infer dtype on variable 'posterior'`
(arviz/xarray incompatibility) after all results had been written — `sed_out/ariadne_result.nc` is a partial file and
must be ignored; `BMA.pkl` and the `.dat` files are complete. (F2) three of the seven recommended grids not fitted (see
above); the four fitted grids agree on R to ±0.014 R☉ and on Teff to ±50 K. (F3) radius and Teff errors are the fit's
formal errors only; no systematic floor added. (F4) Gaia DR2 photometry was used by the package's Librarian, not DR3.

## 3. Verdict

* An unresolved companion with q ≳ 0.63 (M ≳ 0.38 M☉) is excluded at 2σ by the Gaia XP overluminosity test (§1);
  q ≳ 0.4 (M ≳ 0.23 M☉) is disfavoured by the SED colour (§2). Neither test reaches the M5–M8 companions
  (0.08–0.16 M☉) that could still host the transits; speckle/AO imaging remains the only way to close that channel.
* The star's radius is confirmed independently (0.621 ± 0.015 R☉), so the 1.44 R⊕ planet radius does not rest on the
  Mann relation alone.
* The FPP numbers in `../deep_followup_tic4206066/fpp/FPP_REPORT.md` were **not** recomputed with these limits.
