import LeanRank.Rank

/-!
# Forged claims: examples 1–4 must be REFUSED by Lean (expected: exactly 4 errors, at lines 19, 23, 27, 33)

Written by an AI model; not read by a person.  Run from `lean_rank/`:
`lake env lean scripts/RankForge.lean` must exit non-zero with one error in each of examples 1–4;
example 5 (the honest version of 4) must be accepted.
The proved negative controls (a dependent triple is refused, rank-1 consistency) are in
`scripts/RankChecks.lean`.
-/

namespace MagicSquareOfSquares.Rank

open Module WeierstrassCurve WeierstrassCurve.Affine MagicSquareOfSquares.TwoDescent
  MagicSquareOfSquares.G1

-- 1. rank ≥ 3 does not follow from the rank ≥ 2 proof
example : 3 ≤ finrank ℤ (Eplus 7 1).toAffine.Point := rank_7_1.1

-- 2. the degenerate ratio n = 0 (singular curve, Δ = 0) is refused: `(0 : ℚ) ≠ 0` cannot be proved
example (h : Theorem2At 7 0 5) : 2 ≤ finrank ℤ (Eplus 7 0).toAffine.Point :=
  two_le_rank_of_theorem2At (by norm_num) (by norm_num) (by norm_num) (by norm_num) h

-- 3. a ratio that is not twice-square ((12, 5): 144 + 25 = 169 = 13²) is not in the list
example : 2 ≤ finrank ℤ (Eplus 12 5).toAffine.Point :=
  (rank_twiceSquare (12, 5, 13) (by decide)).1

-- 4. a forged certificate: Theorem 2 at (7, 1) with a wrong square-root bracket for H(7, 1) = 1828
example : 2 ≤ finrank ℤ (Eplus 7 1).toAffine.Point :=
  two_le_rank_of_theorem2At (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (theorem2_of 7 1 5 (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (nsq 1828 43 (by unfold genH; norm_num) (by norm_num) (by norm_num))
      (nsq 3172 56 (by unfold genH; norm_num) (by norm_num) (by norm_num))
      (nsq 5798416 2407 (by unfold genH; norm_num) (by norm_num) (by norm_num))
      (nsq 14 3 (by norm_num) (by norm_num) (by norm_num))
      (TwoDescent.not_isSquare_of_neg (by norm_num)))

-- 5. control for 4: the same term with the correct bracket 42 (42² < 1828 < 43²) is accepted
example : 2 ≤ finrank ℤ (Eplus 7 1).toAffine.Point :=
  two_le_rank_of_theorem2At (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (theorem2_of 7 1 5 (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (nsq 1828 42 (by unfold genH; norm_num) (by norm_num) (by norm_num))
      (nsq 3172 56 (by unfold genH; norm_num) (by norm_num) (by norm_num))
      (nsq 5798416 2407 (by unfold genH; norm_num) (by norm_num) (by norm_num))
      (nsq 14 3 (by norm_num) (by norm_num) (by norm_num))
      (TwoDescent.not_isSquare_of_neg (by norm_num)))

end MagicSquareOfSquares.Rank
