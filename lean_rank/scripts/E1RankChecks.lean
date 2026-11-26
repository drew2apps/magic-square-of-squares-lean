import LeanRank.E1RankZero

/-!
# Statement audit for `E1K_finrank_eq_zero` (no errors expected)

Written by an AI model (Claude, made by Anthropic) on 26 November 2026 (UTC); not read by a person.

The curve is the `E1K` of `LemmaDBridge.lean` (`y² = x³ + 8(2 + √2)x² − 8(2 + √2)²x` over `ℚ⟮√2⟯ ⊆ ℝ`), the
statements have no hypotheses, and the group is not trivial (it contains `T = (0, 0)` and `(4, 8√2)`), so
`E1K_finite` is about a nontrivial finite group.
-/

open MagicSquareOfSquares.Family MagicSquareOfSquares.E1Descent

example : E1K.a₁ = 0 ∧ E1K.a₂ = 8 * (2 + s2) ∧ E1K.a₃ = 0 ∧ E1K.a₄ = -8 * (2 + s2) ^ 2 ∧ E1K.a₆ = 0 :=
  ⟨rfl, rfl, rfl, rfl, rfl⟩

example : ((s2 : Q2) : ℝ) = √2 := rfl

example : Q2 = IntermediateField.adjoin ℚ {(√2 : ℝ)} := rfl

/-- `(4, 8√2)` is a point of `E1K` (from `E1_points_on` of `LemmaDReduction.lean`). -/
noncomputable def P4 : E1K.Point :=
  .some _ _ (E1K_nonsingular_of_E1Q2 E1_points_on.2.1)

example : P4 ≠ 0 := WeierstrassCurve.Affine.Point.some_ne_zero _

example : Nontrivial E1K.Point := ⟨⟨P4, 0, WeierstrassCurve.Affine.Point.some_ne_zero _⟩⟩

#check (E1K_finrank_eq_zero : Module.finrank ℤ E1K.Point = 0)
#check (E1K_isOfFinAddOrder : ∀ P : E1K.Point, IsOfFinAddOrder P)
#check (E1K_finite : Finite E1K.Point)
#check (MagicSquareOfSquares.E1Rank.E1K_finrank_eq_zero_of_descent)

#print axioms E1K_finrank_eq_zero
#print axioms E1K_isOfFinAddOrder
#print axioms E1K_finite
