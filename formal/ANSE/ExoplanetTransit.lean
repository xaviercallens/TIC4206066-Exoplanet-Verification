import Mathlib.Data.Real.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real

namespace ANSE.ExoplanetTransit

/-!
# Exoplanet Transit & Keplerian Density Invariants in Lean 4
Formal mathematical specification of the Socratic Elenchus transit vetting framework
for TIC 4206066:
1. Kepler's Third Law orbital dynamics
2. Transit-derived mean stellar density invariant (rho_* = (3*pi / G*P^2) * (a/R_*)^3)
3. Equivalence theorem between Keplerian transit density and physical stellar density
4. Odd-even transit depth symmetry for planetary transits (Elenchus null hypothesis)
-/

/-- Astrophysical Stellar-Planetary System parameters -/
structure TransitSystem where
  M_star : ℝ      -- Stellar mass (kg or M_sun)
  R_star : ℝ      -- Stellar radius (m or R_sun)
  P : ℝ           -- Orbital period (s or days)
  a : ℝ           -- Semi-major axis (m or AU)
  G : ℝ           -- Gravitational constant
  hM : 0 < M_star
  hR : 0 < R_star
  hP : 0 < P
  ha : 0 < a
  hG : 0 < G

/-- Kepler's Third Law semi-major axis cube relation: a^3 = (G * M * P^2) / (4 * pi^2) -/
def kepler_third_law (sys : TransitSystem) : Prop :=
  sys.a ^ 3 = (sys.G * sys.M_star * sys.P ^ 2) / (4 * Real.pi ^ 2)

/-- Mean stellar density rho = M / ((4/3) * pi * R^3) -/
noncomputable def stellar_density (M R : ℝ) : ℝ :=
  M / ((4 / 3) * Real.pi * R ^ 3)

/-- Transit-derived mean stellar density rho_transit = (3 * pi / (G * P^2)) * (a / R)^3 -/
noncomputable def transit_density (G P a R : ℝ) : ℝ :=
  ((3 * Real.pi) / (G * P ^ 2)) * ((a / R) ^ 3)

/-- Theorem: Mean stellar density is strictly positive for positive mass and radius -/
theorem stellar_density_pos (sys : TransitSystem) :
    0 < stellar_density sys.M_star sys.R_star := by
  dsimp [stellar_density]
  have hpi : 0 < Real.pi := Real.pi_pos
  have hconst : 0 < (4 : ℝ) / 3 * Real.pi := mul_pos (by norm_num) hpi
  have hR3 : 0 < sys.R_star ^ 3 := pow_pos sys.hR 3
  have hdenom : 0 < (4 : ℝ) / 3 * Real.pi * sys.R_star ^ 3 := mul_pos hconst hR3
  exact div_pos sys.hM hdenom

/--
Master Invariant Theorem:
Under Kepler's Third Law, the transit-derived density calculated solely from
the orbital period P and normalized separation (a/R_*) identically equals
the true physical spherical stellar density M / ((4/3)*pi*R^3).
-/
theorem kepler_transit_stellar_density_invariant (sys : TransitSystem)
    (h_kepler : kepler_third_law sys) :
    transit_density sys.G sys.P sys.a sys.R_star =
      stellar_density sys.M_star sys.R_star := by
  dsimp [transit_density, stellar_density, kepler_third_law] at *
  have hpi_pos : Real.pi ≠ 0 := Real.pi_ne_zero
  have hG_pos : sys.G ≠ 0 := ne_of_gt sys.hG
  have hP_pos : sys.P ≠ 0 := ne_of_gt sys.hP
  have hR_pos : sys.R_star ≠ 0 := ne_of_gt sys.hR
  have hP2_pos : sys.P ^ 2 ≠ 0 := by positivity
  have hR3_pos : sys.R_star ^ 3 ≠ 0 := by positivity

  -- (a / R)^3 = a^3 / R^3
  rw [div_pow]
  rw [h_kepler]
  -- Algebraic simplification of field division
  field_simp

/-- Geometric transit depth delta = (R_p / R_*)^2 -/
noncomputable def transit_depth (R_p R_star : ℝ) : ℝ :=
  (R_p / R_star) ^ 2

/-- Planetary transit depth is non-negative -/
theorem transit_depth_nonneg (R_p R_star : ℝ) :
    0 ≤ transit_depth R_p R_star := by
  dsimp [transit_depth]
  exact sq_nonneg (R_p / R_star)

/-- Strict positivity of transit depth for non-zero planetary radius -/
theorem transit_depth_pos (R_p R_star : ℝ) (hp : 0 < R_p) (hs : 0 < R_star) :
    0 < transit_depth R_p R_star := by
  dsimp [transit_depth]
  have hdiv : 0 < R_p / R_star := div_pos hp hs
  exact sq_pos_of_pos hdiv

/-- Planetary transit periodic signal model parameterized by constant depth delta_0 -/
def planet_transit_depth_at_epoch (delta_0 : ℝ) (_epoch : ℤ) : ℝ :=
  delta_0

/--
Elenchus Odd-Even Null Hypothesis Theorem:
For any genuine planetary transit around a single non-variable star,
the transit depth at any odd epoch (2*k + 1) identically equals the transit depth
at any even epoch (2*k), ruling out alternating temperature eclipsing binary models.
-/
theorem odd_even_transit_depth_symmetric (delta_0 : ℝ) (k : ℤ) :
    planet_transit_depth_at_epoch delta_0 (2 * k + 1) -
    planet_transit_depth_at_epoch delta_0 (2 * k) = 0 := by
  dsimp [planet_transit_depth_at_epoch]
  ring

/--
Eclipsing Binary Asymmetry Definition:
An eclipsing binary exhibits depth disparity delta_primary != delta_secondary.
-/
def eclipsing_binary_depth_at_epoch (delta_prim delta_sec : ℝ) (epoch : ℤ) : ℝ :=
  if epoch % 2 = 0 then delta_prim else delta_sec

/--
Theorem: If an eclipsing binary has primary depth delta_prim and secondary depth delta_sec
with delta_prim ≠ delta_sec, the odd-even difference is non-zero.
-/
theorem eclipsing_binary_odd_even_disparity
    (delta_prim delta_sec : ℝ) (h_diff : delta_prim ≠ delta_sec) (k : ℤ) :
    eclipsing_binary_depth_at_epoch delta_prim delta_sec (2 * k) ≠
    eclipsing_binary_depth_at_epoch delta_prim delta_sec (2 * k + 1) := by
  dsimp [eclipsing_binary_depth_at_epoch]
  have h_even : (2 * k) % 2 = 0 := by
    omega
  have h_not_zero : ¬((2 * k + 1) % 2 = 0) := by
    omega
  simp only [h_even, ite_true]
  simp only [h_not_zero, ite_false]
  exact h_diff

end ANSE.ExoplanetTransit
