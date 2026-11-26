import LeanRank.E1Descent

/-!
# `E1(ℚ(√2))` is a finite group: every point has finite order

Written by an AI model (Claude, made by Anthropic) on 26 November 2026 (UTC); not read by a person.
Experiment, not counted in the census.

From `E1K_finrank_eq_zero` (rank 0, `E1Descent.lean`) and Stoll's Mordell–Weil `E1K_fg`:
`E1K_isOfFinAddOrder : ∀ P : E1K.Point, IsOfFinAddOrder P` (the form the torsion computation needs to
conclude `E1PointsK`, hence Lemma D), and `E1K_finite : Finite E1K.Point`. `lemmaD_of_torsion`: Lemma D
follows from the statement that the affine points of finite order have `x ∈ {0, ±4, ±(12 + 8√2)}`.
-/

namespace MagicSquareOfSquares.E1Descent

open MagicSquareOfSquares.Family

instance : IsAddFG E1K.Point := E1K_fg

/-- **Every point of `E1(ℚ(√2))` has finite order.** -/
theorem E1K_isOfFinAddOrder (P : E1K.Point) : IsOfFinAddOrder P := by
  have h : Module.IsTorsion ℤ E1K.Point :=
    Module.finrank_eq_zero_iff_isTorsion.mp E1K_finrank_eq_zero
  obtain ⟨a, ha⟩ := @h P
  exact isOfFinAddOrder_iff_zsmul_eq_zero.mpr
    ⟨a, nonZeroDivisors.coe_ne_zero a, by rwa [Submonoid.smul_def] at ha⟩

/-- **`E1(ℚ(√2))` is finite.** -/
theorem E1K_finite : Finite E1K.Point :=
  AddCommGroup.finite_of_fg_isAddTorsion _ E1K_isOfFinAddOrder

/-- **`E1PointsK` from the torsion alone**: if every affine point of finite order of `E1K` has
`x ∈ {0, ±4, ±(12 + 8√2)}`, then every affine point has (rank 0 is `E1K_isOfFinAddOrder`). -/
theorem E1PointsK_of_torsion
    (h : ∀ (x y : Q2) (hxy : E1K.Nonsingular x y),
      IsOfFinAddOrder (WeierstrassCurve.Affine.Point.some x y hxy) →
      x = 0 ∨ x = 4 ∨ x = -4 ∨ x = 12 + 8 * s2 ∨ x = -12 - 8 * s2) : E1PointsK := by
  intro x y hxy
  exact h x y (E1K.equation_iff_nonsingular.mp hxy) (E1K_isOfFinAddOrder _)

/-- **Lemma D from the torsion of `E1(ℚ(√2))` alone** (the rank part is now proved). -/
theorem lemmaD_of_torsion
    (h : ∀ (x y : Q2) (hxy : E1K.Nonsingular x y),
      IsOfFinAddOrder (WeierstrassCurve.Affine.Point.some x y hxy) →
      x = 0 ∨ x = 4 ∨ x = -4 ∨ x = 12 + 8 * s2 ∨ x = -12 - 8 * s2) : LemmaD :=
  lemmaD_of_E1PointsK (E1PointsK_of_torsion h)

end MagicSquareOfSquares.E1Descent

#print axioms MagicSquareOfSquares.E1Descent.E1K_isOfFinAddOrder
#print axioms MagicSquareOfSquares.E1Descent.E1K_finite
#print axioms MagicSquareOfSquares.E1Descent.E1PointsK_of_torsion
#print axioms MagicSquareOfSquares.E1Descent.lemmaD_of_torsion
