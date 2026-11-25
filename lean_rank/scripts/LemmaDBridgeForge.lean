import LeanRank.LemmaDBridge

/-!
Negative tests for `LeanRank/LemmaDBridge.lean` (written by an AI model; not read by a person).
Each `example` marked REFUSED must produce exactly one error; the ones marked ACCEPTED must compile.
-/

namespace MagicSquareOfSquares.Family

-- ACCEPTED: dropping the abscissa 12 + 8√2 from the Mathlib-form hypothesis makes it false
-- (the image (12 + 8√2, 80 + 56√2) of (−1, −4) is a point of `E1K`).
example (h : ∀ x y : Q2, E1K.Equation x y → x = 0 ∨ x = 4 ∨ x = -4 ∨ x = -12 - 8 * s2) : False := by
  have hN := E1K.equation_iff_nonsingular.mpr (E1K_nonsingular_of_E1Q2 E1_points_on.2.2.2.2.2.1)
  have h2 := Real.sqrt_nonneg 2
  rcases h _ _ hN with hx | hx | hx | hx <;> have hx' := congrArg (algebraMap Q2 ℝ) hx <;>
    simp only [map_add, map_mul, map_ratCast, algebraMap_s2, map_zero, map_ofNat, map_neg,
      map_sub] at hx' <;> push_cast at hx' <;> linarith

-- REFUSED (forged point of `E1K`): (4, 4√2).
example : E1K.Equation 4 (4 * s2) := by
  rw [WeierstrassCurve.Affine.equation_iff]
  apply (algebraMap Q2 ℝ).injective
  have hr : (√2 : ℝ) ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  simp only [E1K, map_add, map_mul, map_pow, map_neg, map_zero, map_ofNat, algebraMap_s2]
  linear_combination 96 * hr

-- REFUSED (forged discriminant).
example : E1K.Δ = 98304 * (2 + s2) ^ 5 := by
  simp only [E1K, WeierstrassCurve.Δ, WeierstrassCurve.b₂, WeierstrassCurve.b₄, WeierstrassCurve.b₆,
    WeierstrassCurve.b₈]
  ring

-- ACCEPTED (control): the honest point (4, 8√2) and the honest statement.
example : E1K.Equation 4 (8 * s2) := by
  rw [WeierstrassCurve.Affine.equation_iff]
  apply (algebraMap Q2 ℝ).injective
  have hr : (√2 : ℝ) ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  simp only [E1K, map_add, map_mul, map_pow, map_neg, map_zero, map_ofNat, algebraMap_s2]
  linear_combination 96 * hr
example (h : E1PointsK) : LemmaD := lemmaD_of_E1PointsK h

end MagicSquareOfSquares.Family
