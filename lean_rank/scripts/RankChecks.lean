import LeanRank.Index

/-!
# Checks and negative controls for `LeanRank/Rank.lean` (not part of the library)

Written by an AI model; not read by a person.  Run from `lean_rank/`:
`lake env lean scripts/RankChecks.lean`.
-/

namespace MagicSquareOfSquares.Rank

open Module WeierstrassCurve WeierstrassCurve.Affine MagicSquareOfSquares.TwoDescent

/-- Negative control: a dependent triple is refused by the independence hypothesis.  With
`Q := R`, `1·R + 1·R + 0·T = R + R` has odd coefficients. -/
theorem forged_dependent {G : Type*} [AddCommGroup G] (R T : G) :
    ¬ ∀ c₁ c₂ c₃ : ℤ, ∀ W : G, c₁ • R + c₂ • R + c₃ • T = W + W →
      Even c₁ ∧ Even c₂ ∧ Even c₃ := by
  intro h
  exact Int.not_even_one (h 1 1 0 R (by simp)).1

/-- The same on the curve `E⁺(7, 1)`: replacing `Q` by `R` in `Theorem2At 7 1 5` gives a false
statement, so a certificate with two dependent points cannot be accepted. -/
theorem forged_7_1 (hR : (Eplus 7 1).toAffine.Nonsingular (xR 7 1) (yR 7 1))
    (hT : (Eplus 7 1).toAffine.Nonsingular (eP 7 1) 0) :
    ¬ ∀ c₁ c₂ c₃ : ℤ, ∀ W : (Eplus 7 1).toAffine.Point,
      c₁ • Point.some _ _ hR + c₂ • Point.some _ _ hR + c₃ • Point.some _ _ hT = W + W →
        Even c₁ ∧ Even c₂ ∧ Even c₃ :=
  forged_dependent _ _

/-- `2·T = O + O`: `T` has order 2 (the factor `#E[2] = 2` in `[E : 2E] = 2 ^ (rank + 1)`). -/
example (hT : (Eplus 7 1).toAffine.Nonsingular (eP 7 1) 0) :
    (2 : ℤ) • Point.some _ _ hT = (0 : (Eplus 7 1).toAffine.Point) + 0 := by
  rw [two_zsmul, T_add_T, add_zero]

/-- Consistency of the group-theory lemma with a group of rank 1: in `ℤ` no three elements are
independent modulo `2ℤ` (since `finrank ℤ ℤ = 1`). -/
theorem no_indep_in_int (R Q T : ℤ) :
    ¬ ∀ c₁ c₂ c₃ : ℤ, ∀ W : ℤ, c₁ • R + c₂ • Q + c₃ • T = W + W →
      Even c₁ ∧ Even c₂ ∧ Even c₃ := by
  intro h
  have htors : ∀ P : ℤ, P + P = 0 → P = 0 ∨ P = 0 := fun P hP => Or.inl (by omega)
  have := two_le_finrank_of_indep (T₀ := 0) htors h
  rw [Module.finrank_self] at this
  omega

/-- `[E⁺(7, 1) : 2E⁺(7, 1)] ≥ 8`. -/
example : 8 ≤ (nsmulAddMonoidHom (α := (Eplus 7 1).toAffine.Point) 2).range.index :=
  eight_le_index_of_theorem2At (by norm_num) (by norm_num) (by norm_num) (by norm_num) t2_7_1_p

end MagicSquareOfSquares.Rank

#print axioms MagicSquareOfSquares.Rank.forged_dependent
#print axioms MagicSquareOfSquares.Rank.forged_7_1
#print axioms MagicSquareOfSquares.Rank.no_indep_in_int
