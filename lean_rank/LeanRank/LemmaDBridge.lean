import LeanRank.LemmaDReduction

/-!
# `E1` over the number field `ℚ(√2)` as a Mathlib Weierstrass curve

Written by an AI model (Claude, made by Anthropic) on 25 November 2026 (UTC); not read by a person.
Experiment, not counted in the census.

`Q2 = ℚ⟮√2⟯ ⊆ ℝ` (a Mathlib `IntermediateField`, with a `NumberField` instance), `s2 = √2 ∈ Q2`, and
`E1K : y² = x³ + 8(2 + √2)x² − 8(2 + √2)²x` as a `WeierstrassCurve.Affine Q2`. `E1PointsK` says that every
affine point of `E1K` (Mathlib's `Equation`) has `x ∈ {0, ±4, ±(12 + 8√2)}`; `lemmaD_of_E1PointsK` derives
`LemmaD` from it (through `E1Points_of_E1PointsK` and `lemmaD_of_E1Points`). This puts the one remaining
hypothesis of Lemma D in the language of Stoll's EllipticCurves library (Mordell–Weil over number fields,
2-descent, reduction at primes), where it could be proved.

Also: `E1K.IsElliptic` (`Δ_E1K : Δ = 98304(2 + √2)⁶`), `E1K_fg` (`E1(ℚ(√2))` is finitely generated, by Stoll's
`fg_point_of_numberField`), and `E1K_nonsingular_of_E1Q2` (a point `(a + b√2, c + e√2)` of `E1Q2` is a point of
Mathlib's `E1K`; with `E1_points_on` this gives the nine affine points of the group of order 10 in `E1K`).
-/

namespace MagicSquareOfSquares.Family

open Polynomial in
theorem isIntegral_sqrt_two : IsIntegral ℚ (√2 : ℝ) :=
  ⟨X ^ 2 - C 2, monic_X_pow_sub_C 2 two_ne_zero, by simp⟩

/-- `ℚ(√2)` inside `ℝ`. -/
noncomputable abbrev Q2 : IntermediateField ℚ ℝ := IntermediateField.adjoin ℚ {(√2 : ℝ)}

instance : FiniteDimensional ℚ Q2 := IntermediateField.adjoin.finiteDimensional isIntegral_sqrt_two

instance : NumberField Q2 := ⟨⟩

/-- `√2` as an element of `Q2`. -/
noncomputable def s2 : Q2 := ⟨√2, IntermediateField.mem_adjoin_simple_self ℚ _⟩

theorem algebraMap_s2 : algebraMap Q2 ℝ s2 = √2 := rfl

/-- `E1 : y² = x³ + 8dx² − 8d²x`, `d = 2 + √2`, over `Q2`. -/
noncomputable def E1K : WeierstrassCurve.Affine Q2 := ⟨0, 8 * (2 + s2), 0, -8 * (2 + s2) ^ 2, 0⟩

/-- **Hypothesis (Mathlib form)**: the affine points of `E1K` have `x ∈ {0, ±4, ±(12 + 8√2)}`. -/
def E1PointsK : Prop :=
  ∀ x y : Q2, E1K.Equation x y → x = 0 ∨ x = 4 ∨ x = -4 ∨ x = 12 + 8 * s2 ∨ x = -12 - 8 * s2

theorem E1Points_of_E1PointsK (h : E1PointsK) : E1Points := by
  intro a b c e hE
  have hxR : algebraMap Q2 ℝ ((a : Q2) + (b : Q2) * s2) = (a : ℝ) + (b : ℝ) * √2 := by
    simp only [map_add, map_mul, map_ratCast, algebraMap_s2]
  have hyR : algebraMap Q2 ℝ ((c : Q2) + (e : Q2) * s2) = (c : ℝ) + (e : ℝ) * √2 := by
    simp only [map_add, map_mul, map_ratCast, algebraMap_s2]
  have heq : E1K.Equation ((a : Q2) + (b : Q2) * s2) ((c : Q2) + (e : Q2) * s2) := by
    rw [WeierstrassCurve.Affine.equation_iff]
    apply (algebraMap Q2 ℝ).injective
    simp only [E1K, map_add, map_mul, map_pow, map_neg, map_zero, map_ofNat, map_ratCast, algebraMap_s2]
    unfold E1Q2 E1eq at hE
    linear_combination hE
  rcases h _ _ heq with hx | hx | hx | hx | hx
  · have hx' := congrArg (algebraMap Q2 ℝ) hx
    simp only [map_add, map_mul, map_ratCast, algebraMap_s2, map_zero] at hx'
    exact Or.inl (q2_eq_zero hx')
  · have hx' := congrArg (algebraMap Q2 ℝ) hx
    simp only [map_add, map_mul, map_ratCast, algebraMap_s2, map_ofNat] at hx'
    obtain ⟨h1, h2⟩ := q2_eq_zero (a := a - 4) (b := b) (by push_cast; linear_combination hx')
    exact Or.inr (Or.inl ⟨by linarith, h2⟩)
  · have hx' := congrArg (algebraMap Q2 ℝ) hx
    simp only [map_add, map_mul, map_ratCast, algebraMap_s2, map_neg, map_ofNat] at hx'
    obtain ⟨h1, h2⟩ := q2_eq_zero (a := a + 4) (b := b) (by push_cast; linear_combination hx')
    exact Or.inr (Or.inr (Or.inl ⟨by linarith, h2⟩))
  · have hx' := congrArg (algebraMap Q2 ℝ) hx
    simp only [map_add, map_mul, map_ratCast, algebraMap_s2, map_ofNat] at hx'
    obtain ⟨h1, h2⟩ := q2_eq_zero (a := a - 12) (b := b - 8) (by push_cast; linear_combination hx')
    exact Or.inr (Or.inr (Or.inr (Or.inl ⟨by linarith, by linarith⟩)))
  · have hx' := congrArg (algebraMap Q2 ℝ) hx
    simp only [map_add, map_mul, map_ratCast, algebraMap_s2, map_sub, map_neg, map_ofNat] at hx'
    obtain ⟨h1, h2⟩ := q2_eq_zero (a := a + 12) (b := b + 8) (by push_cast; linear_combination hx')
    exact Or.inr (Or.inr (Or.inr (Or.inr ⟨by linarith, by linarith⟩)))

/-- **Lemma D from the Mathlib-form hypothesis on `E1K` over the number field `Q2`.** -/
theorem lemmaD_of_E1PointsK (h : E1PointsK) : LemmaD :=
  lemmaD_of_E1Points (E1Points_of_E1PointsK h)

/-! ### `E1K` is an elliptic curve, so Stoll's Mordell–Weil theorem applies -/

instance : E1K.IsCharNeTwoNF := ⟨rfl, rfl⟩

theorem Δ_E1K : E1K.Δ = 98304 * (2 + s2) ^ 6 := by
  simp only [E1K, WeierstrassCurve.Δ, WeierstrassCurve.b₂, WeierstrassCurve.b₄, WeierstrassCurve.b₆,
    WeierstrassCurve.b₈]
  ring

theorem two_add_s2_ne_zero : (2 + s2 : Q2) ≠ 0 := by
  intro h
  have h' := congrArg (algebraMap Q2 ℝ) h
  rw [map_add, map_ofNat, algebraMap_s2, map_zero] at h'
  have : (0 : ℝ) < 2 + √2 := by positivity
  linarith

instance : E1K.IsElliptic := by
  rw [WeierstrassCurve.isElliptic_iff, isUnit_iff_ne_zero, Δ_E1K]
  exact mul_ne_zero (by norm_num) (pow_ne_zero _ two_add_s2_ne_zero)

/-- `E1(ℚ(√2))` is finitely generated (Stoll's `fg_point_of_numberField`). -/
theorem E1K_fg : IsAddFG E1K.Point := WeierstrassCurve.Affine.fg_point_of_numberField

/-- A point `(a + b√2, c + e√2)` of `E1Q2` is a (nonsingular) point of Mathlib's `E1K`. -/
theorem E1K_nonsingular_of_E1Q2 {a b c e : ℚ} (hE : E1Q2 a b c e) :
    E1K.Nonsingular ((a : Q2) + (b : Q2) * s2) ((c : Q2) + (e : Q2) * s2) := by
  have hxR : algebraMap Q2 ℝ ((a : Q2) + (b : Q2) * s2) = (a : ℝ) + (b : ℝ) * √2 := by
    simp only [map_add, map_mul, map_ratCast, algebraMap_s2]
  have hyR : algebraMap Q2 ℝ ((c : Q2) + (e : Q2) * s2) = (c : ℝ) + (e : ℝ) * √2 := by
    simp only [map_add, map_mul, map_ratCast, algebraMap_s2]
  apply E1K.equation_iff_nonsingular.mp
  rw [WeierstrassCurve.Affine.equation_iff]
  apply (algebraMap Q2 ℝ).injective
  simp only [E1K, map_add, map_mul, map_pow, map_neg, map_zero, map_ofNat, map_ratCast, algebraMap_s2]
  unfold E1Q2 E1eq at hE
  linear_combination hE

end MagicSquareOfSquares.Family

#print axioms MagicSquareOfSquares.Family.E1Points_of_E1PointsK
#print axioms MagicSquareOfSquares.Family.lemmaD_of_E1PointsK
#print axioms MagicSquareOfSquares.Family.E1K_fg
#print axioms MagicSquareOfSquares.Family.E1K_nonsingular_of_E1Q2
