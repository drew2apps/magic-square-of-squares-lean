import LeanRank.E1RankZero
import LeanRank.TorsionE1
import LeanRank.LemmaDEquiv
import LeanRank.Theorem1Family
import LeanRank.Theorem1Order

/-!
# Theorems 1 and 2 of `UNIFORM_RANK.md`, with no hypothesis left

Written by an AI model; not read by a person.

Lemma D (`LemmaD`, `Theorem2Family.lean`) is proved here from the two halves of `E1(ℚ(√2))`:
rank 0 (`E1Descent.E1K_isOfFinAddOrder`, `E1Descent.lean`/`E1RankZero.lean`: a 2-isogeny descent over
`ℤ[√2]`) and torsion (`TorsionE1.torsion_abscissa`, `TorsionE1.lean`: reduction at a prime above 7), through
the reduction `lemmaD_of_E1PointsK` (`LemmaDReduction.lean`, `LemmaDBridge.lean`). The family theorems of
`Theorem1Family.lean`, `Theorem1Order.lean` and `LemmaDEquiv.lean`, which take `LemmaD` as their only
hypothesis, are restated below without it.

The models: `TwoDescent.Eplus m n` is `E⁺` at the ratio `m : n` and `TwoDescent.Eplus m (−n)` is `E⁻`
(`TwoDescentMap.lean`). Pairs are coprime integers.
-/

namespace MagicSquareOfSquares.UniformRank

open MagicSquareOfSquares.Family MagicSquareOfSquares.G1 MagicSquareOfSquares.TwoDescent

/-- **Lemma D** of `UNIFORM_RANK.md`, kernel-checked. -/
theorem lemmaD : LemmaD :=
  MagicSquareOfSquares.E1Descent.lemmaD_of_torsion TorsionE1.torsion_abscissa

/-- `E1(ℚ(√2))` has exactly 10 points. -/
theorem card_E1K : Nat.card E1K.Point = 10 :=
  TorsionE1.natCard_E1K_of_finrank_eq_zero MagicSquareOfSquares.E1Descent.E1K_finrank_eq_zero

/-- **Theorem 1** (rank `≥ 1`): for coprime `m, n` with `mn(m² − n²) ≠ 0`, both `E⁺` and `E⁻` have rank at
least 1. -/
theorem theorem1 {m n : ℤ} (hmn : IsCoprime m n) (hm : m ≠ 0) (hn : n ≠ 0) (hne : m ^ 2 ≠ n ^ 2) :
    1 ≤ Module.finrank ℤ (TwoDescent.Eplus (m : ℚ) (n : ℚ)).toAffine.Point ∧
      1 ≤ Module.finrank ℤ (TwoDescent.Eplus (m : ℚ) (-(n : ℚ))).toAffine.Point :=
  rank_one_family lemmaD hmn hm hn hne

/-- Theorem 1, the point: the point `P` of `E⁺` is nonsingular and of infinite order. -/
theorem theorem1_P {m n : ℤ} (hmn : IsCoprime m n) (hm : m ≠ 0) (hn : n ≠ 0) (hne : m ^ 2 ≠ n ^ 2) :
    ∃ hP : (TwoDescent.Eplus (m : ℚ) (n : ℚ)).toAffine.Nonsingular (xP m n) (yP m n),
      ∀ N : ℕ, 0 < N → N • WeierstrassCurve.Affine.Point.some _ _ hP ≠ 0 :=
  P_infinite_order lemmaD hmn hm hn hne

/-- **Theorem 2** (rank `≥ 2`): for coprime `m, n` with `m² + n² = 2k²` and `m² ≠ n²`, the points `R, Q, T`
are independent in `E(ℚ)/2E(ℚ)` on both curves, and both have rank at least 2. -/
theorem theorem2 {m n k : ℤ} (hmn : IsCoprime m n) (hk : m ^ 2 + n ^ 2 = 2 * k ^ 2) (hne : m ^ 2 ≠ n ^ 2) :
    TwoDescent.Theorem2At (m : ℚ) (n : ℚ) (k : ℚ) ∧ TwoDescent.Theorem2At (m : ℚ) (-(n : ℚ)) (k : ℚ) ∧
      2 ≤ Module.finrank ℤ (TwoDescent.Eplus (m : ℚ) (n : ℚ)).toAffine.Point ∧
      2 ≤ Module.finrank ℤ (TwoDescent.Eplus (m : ℚ) (-(n : ℚ))).toAffine.Point :=
  rank_family_all lemmaD hmn hk hne

end MagicSquareOfSquares.UniformRank

#print axioms MagicSquareOfSquares.UniformRank.lemmaD
#print axioms MagicSquareOfSquares.UniformRank.card_E1K
#print axioms MagicSquareOfSquares.UniformRank.theorem1
#print axioms MagicSquareOfSquares.UniformRank.theorem1_P
#print axioms MagicSquareOfSquares.UniformRank.theorem2
