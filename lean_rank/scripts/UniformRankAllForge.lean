import LeanRank.UniformRankAll
import LeanRank.CensusClosureMSS

/-!
Forged claims against `EminusModel.lean`, `UniformRankAll.lean` and `CensusClosureMSS.lean` (run with
`lake env lean scripts/UniformRankAllForge.lean` in `lean_rank/`). Written by an AI model; not read by a person.
Each claim uses a theorem outside its hypotheses, or states a false identity; Lean must refuse each.
Expected: exactly 9 errors, one per claim, at lines 26 30 34 38 41 45 56 60 64 (and nothing else fails).

1 Theorem 2 at the square-type ratio `(12, 5)` (`k = 13/2` does not give `m² + n² = 2k²`)
2 Theorem 1 at `m = n`
3 Lemma B at `m = n = 1/2`
4 the parity fact at `(1, 2)`, which is not a twice-square pair (`2·1·2 = 4` is a square)
5 `E⁻ = E⁺` at `(7, 1)` (the coefficient `a₂`)
6 the recipe at the wrong root `z0 = 2` for `q⁺` gives `E⁺` (the coefficient `a₂`)
7 the map of the recipe at `z = z0`
8 `(5, 1)` claimed to be one of the six census-closure ratios
9 the census closure at `4 : 3` without distinct entries (the square of `±1`'s)
-/

open Module MagicSquareOfSquares.TwoDescent MagicSquareOfSquares.EminusModel
open MagicSquareOfSquares.UniformRankAll LeanRank.CensusClosureMSS

/-! 1 -/
example : 2 ≤ finrank ℤ (Eplus 12 5).toAffine.Point :=
  (theorem2_all (m := 12) (n := 5) (k := 13 / 2) (by norm_num) (by norm_num)).2.2.1

/-! 2 -/
example : 1 ≤ finrank ℤ (Eplus 3 3).toAffine.Point :=
  (theorem1_all (m := 3) (n := 3) (by norm_num) (by norm_num) (by norm_num)).1

/-! 3 -/
example : ¬ IsSquare ((1 / 2 : ℚ) ^ 8 + 14 * (1 / 2) ^ 4 * (1 / 2) ^ 4 + (1 / 2) ^ 8) :=
  lemmaB_rat_all (by norm_num) (by norm_num) (by norm_num)

/-! 4 -/
example : ¬ IsSquare (2 * (1 : ℚ) * 2) :=
  not_isSquare_two_mul_all (k := 1) (by norm_num) (by norm_num)

/-! 5 -/
example : (Eminus 7 1).a₂ = (Eplus 7 1).a₂ := by
  rw [Eminus_eq]; simp only [Eplus, curve, aP, eP]; norm_num

/-! 6 -/
example : (weierstrassOf (qPoly 2 7 1) 2).a₂ = (Eplus 7 1).a₂ := by
  simp only [weierstrassOf, Eplus, curve, aP, eP, qPoly, gPoly, Polynomial.derivative_mul,
    Polynomial.derivative_add, Polynomial.derivative_sub, Polynomial.derivative_X, Polynomial.derivative_C,
    Polynomial.derivative_X_sq, Polynomial.derivative_one, Polynomial.derivative_zero, Polynomial.eval_mul,
    Polynomial.eval_add, Polynomial.eval_sub, Polynomial.eval_C, Polynomial.eval_X, Polynomial.eval_pow,
    Polynomial.eval_one, Polynomial.eval_zero]
  norm_num

/-! 7 -/
example : (Eplus 7 1).toAffine.Equation (mapA (qPoly 2 7 1) (-2) / (-2 - (-2)))
    (mapA (qPoly 2 7 1) (-2) * 0 / (-2 - (-2)) ^ 2) :=
  EplusModel_map (m := 7) (n := 1) (z := -2) (Y := 0) (by norm_num) (by simp [qPoly])

/-! 8 -/
example : ((5 : ℤ), (1 : ℤ)) ∈ censusClosureRatios := by
  decide

/-! 9 -/
example : ((1 : ℤ) + -1) * 3 ≠ (1 + -1) * 4 :=
  no_mss_4_3 control_needs_distinct.1 (by decide)

/-! Positive control (must compile): Theorem 2 at the non-coprime pair `(14, 2)`, `k = 10`. -/
example : 2 ≤ finrank ℤ (Eplus 14 (-2)).toAffine.Point :=
  (theorem2_all (m := 14) (n := 2) (k := 10) (by norm_num) (by norm_num)).2.2.2
