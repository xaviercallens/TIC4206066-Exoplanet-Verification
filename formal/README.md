# Formal Lean 4 Verification: Exoplanet Transit & Keplerian Invariants

This directory contains the machine-checked mathematical specification and formal proofs for the **TIC 4206066 Exoplanet Verification** pipeline, written in **Lean 4** and checked against **Mathlib**.

## 1. Formalized Theorems in `ANSE/ExoplanetTransit.lean`

All theorems are formally proved with **zero axioms** and **zero `sorry` placeholders** (`zero_sorry_soundness: true`):

| Theorem / Definition | Mathematical Formulation | Epistemic Significance |
| :--- | :--- | :--- |
| `kepler_third_law` | $a^3 = \frac{G M_* P^2}{4 \pi^2}$ | Kepler's Third Law governing 2-body orbital mechanics. |
| `stellar_density` | $\rho_* = \frac{M_*}{\frac{4}{3} \pi R_*^3}$ | Physical mean stellar density of a spherical star. |
| `transit_density` | $\rho_{\text{transit}} = \frac{3 \pi}{G P^2} \left(\frac{a}{R_*}\right)^3$ | Mean density derived purely from transit observables ($P$, $a/R_*$). |
| `stellar_density_pos` | $M_* > 0 \land R_* > 0 \implies \rho_* > 0$ | Density positivity invariant. |
| `kepler_transit_stellar_density_invariant` | $\text{kepler\_third\_law} \implies \rho_{\text{transit}} = \rho_*$ | **Master Equivalence Theorem**: Proves algebraic identity between transit-derived density and host stellar density. |
| `transit_depth` | $\delta = \left(\frac{R_p}{R_*}\right)^2$ | Geometric transit depth relation. |
| `transit_depth_nonneg` | $\forall R_p, R_*, \delta \ge 0$ | Transit depth is strictly non-negative. |
| `transit_depth_pos` | $R_p > 0 \land R_* > 0 \implies \delta > 0$ | Transit depth positivity. |
| `odd_even_transit_depth_symmetric` | $\delta(2k+1) - \delta(2k) = 0$ | **Socratic Elenchus Null Hypothesis**: Constant depth across odd and even epochs rules out simple alternating EBs. |
| `eclipsing_binary_odd_even_disparity` | $\delta_{\text{prim}} \ne \delta_{\text{sec}} \implies \delta(2k) \ne \delta(2k+1)$ | **Falsification Criterion**: Disparity in odd-even depths refutes the single-planet hypothesis. |

---

## 2. Compilation and Verification

To compile and verify the Lean 4 specifications locally:

```bash
# Enter formal directory
cd formal

# Fetch mathlib cache (if mathlib precompiled cache is desired)
lake exe cache get

# Compile the specification
lake build ANSE.ExoplanetTransit

# Verify clean elaboration
lake env lean ANSE/ExoplanetTransit.lean
```

## 3. Toolchain Specification
- **Lean Version:** `leanprover/lean4:v4.34.0-rc2` (specified in `lean-toolchain`)
- **Package Manager:** Lake (`lakefile.toml`)
- **Dependencies:** Mathlib v4.34.0-rc2

## 4. Certification gate (fail-closed)
Workflow step 8 passes only if Lean actually elaborates `ExoplanetTransit.lean` cleanly (exit 0, empty stderr) and the source has no `sorry`. If `lake` is missing, times out (900 s) or Mathlib is unavailable, the step reports `lean_elaboration_exit_code = -1` and **fails**. To reuse a prebuilt Mathlib cache set `LAKE_CACHE_DIR=/path/to/.lake`. Covered by `tests/test_certification_gate.py`.

Scope: Lean checks the stated algebraic invariants only; it does not validate the observational data or the planet.
