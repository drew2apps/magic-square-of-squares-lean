import LeanRank.Theorem1Family
import LeanRank.LemmaDEquiv
import LeanRank.Theorem1Order

/-!
# Forged claims for the family theorems: examples 1–6, 9, 10 must be REFUSED by Lean; 7, 8, 11 accepted

Written by an AI model; not read by a person.  Run from `lean_rank/`:
`lake env lean scripts/FamilyForge.lean` must exit non-zero with at least one error in each of examples 1–6, 9
and 10, and none in examples 7, 8 and 11.
-/

namespace MagicSquareOfSquares.Family

open Module WeierstrassCurve WeierstrassCurve.Affine MagicSquareOfSquares.TwoDescent
  MagicSquareOfSquares.G1 MagicSquareOfSquares.LemmaB

-- 1. Lemma B at the trivial solution m = n = 1 (1 + 14 + 1 = 4²): `1 ^ 2 ≠ 1 ^ 2` cannot be proved
example : ¬ IsSquare ((1 : ℤ) ^ 8 + 14 * 1 ^ 4 * 1 ^ 4 + 1 ^ 8) :=
  lemmaB isCoprime_one_left one_ne_zero one_ne_zero (by norm_num)

-- 2. Pocklington at a = b = 1 (1 − 1 + 1 = 1²): `Odd (1 + 1)` cannot be proved
example : False :=
  pocklington (a := 1) (b := 1) (c := 1) isCoprime_one_left one_ne_zero one_ne_zero (by decide) (by norm_num)

-- 3. Theorem 2 at the square-type ratio (12, 5) (144 + 25 = 13², not twice a square): `169 = 2·13²` is false
example (hD : LemmaD) : Theorem2At 12 5 13 ∧ Theorem2At 12 (-5) 13 := by
  have := theorem2_family hD (m := 12) (n := 5) (k := 13) (Int.isCoprime_iff_gcd_eq_one.mpr rfl)
    (by norm_num) (by norm_num) (by norm_num)
  exact_mod_cast this

-- 4. Theorem 2 at the non-coprime ratio (14, 2) (196 + 4 = 2·10²): `IsCoprime 14 2` cannot be proved
example (hD : LemmaD) : Theorem2At 14 2 10 ∧ Theorem2At 14 (-2) 10 := by
  have := theorem2_family hD (m := 14) (n := 2) (k := 10) (Int.isCoprime_iff_gcd_eq_one.mpr rfl)
    (by norm_num) (by norm_num) (by norm_num)
  exact_mod_cast this

-- 5. rank ≥ 3 does not follow from the family theorem
example (hD : LemmaD) : 3 ≤ finrank ℤ (Eplus 7 1).toAffine.Point := by
  have := (rank_family hD (m := 7) (n := 1) (k := 5) (Int.isCoprime_iff_gcd_eq_one.mpr rfl)
    (by norm_num) (by norm_num) (by norm_num)).1
  exact_mod_cast this

-- 6. Theorem 2 without Lemma D: the hypothesis `¬ IsSquare H(4, 1)` is false (H(4, 1) = 13²), so it cannot be
--    supplied from `control_H_4_1`'s negation
example : ¬ IsSquare (genH 4 1) := control_H_4_1

-- 7. control for 3–5: the twice-square ratio (7, 1) is accepted (hypothesis: Lemma D)
example (hD : LemmaD) : 2 ≤ finrank ℤ (Eplus 7 1).toAffine.Point ∧
    2 ≤ finrank ℤ (Eplus 7 (-1)).toAffine.Point := by
  have := rank_family hD (m := 7) (n := 1) (k := 5) (Int.isCoprime_iff_gcd_eq_one.mpr rfl)
    (by norm_num) (by norm_num) (by norm_num)
  exact_mod_cast this

-- 8. control: rank ≥ 1 at the square-type ratio (12, 5) is accepted (Theorem 1 needs no twice-square)
example (hD : LemmaD) : 1 ≤ finrank ℤ (Eplus 12 5).toAffine.Point := by
  have := (rank_one_family hD (m := 12) (n := 5) (Int.isCoprime_iff_gcd_eq_one.mpr rfl)
    (by norm_num) (by norm_num) (by norm_num)).1
  exact_mod_cast this

-- 9. (refused) Theorem 2 at the degenerate twice-square ratio (1, 1): `1 ^ 2 ≠ 1 ^ 2` cannot be proved
example (hD : LemmaD) : 2 ≤ finrank ℤ (Eplus 1 1).toAffine.Point := by
  have := (rank_family_all hD (m := 1) (n := 1) (k := 1) isCoprime_one_left (by norm_num)
    (by norm_num)).2.2.1
  exact_mod_cast this

-- 10. (refused) infinite order of `R` at the degenerate ratio (1, −1)
example (hD : LemmaD) (hR : (Eplus 1 (-1)).toAffine.Nonsingular (xR 1 (-1)) (yR 1 (-1))) :
    (1 : ℕ) • Point.some _ _ hR ≠ 0 := by
  have := R_infinite_order hD (m := 1) (n := -1) isCoprime_one_left one_ne_zero (by norm_num)
    (by norm_num) (by exact_mod_cast hR) 1 one_pos
  exact_mod_cast this

-- 11. (accepted) control for 9–10: `R` has infinite order on `E⁺` at (7, 1), from Lemma D
example (hD : LemmaD) (hR : (Eplus 7 1).toAffine.Nonsingular (xR 7 1) (yR 7 1)) :
    (5 : ℕ) • Point.some _ _ hR ≠ 0 := by
  have := R_infinite_order hD (m := 7) (n := 1) (Int.isCoprime_iff_gcd_eq_one.mpr rfl) (by norm_num)
    (by norm_num) (by norm_num) (by exact_mod_cast hR) 5 (by norm_num)
  exact_mod_cast this

end MagicSquareOfSquares.Family
