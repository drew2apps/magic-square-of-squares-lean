import LeanRank.CensusClosure

/-!
# Sanity checks for `LeanRank.CensusClosure` (each must SUCCEED)

Written by an AI model; not read by a person.  Not part of the library.  Run with
`lake env lean scripts/CensusClosureChecks.lean` from `lean_rank/`.  Expected: no error.

The hypotheses of the headline theorems are not contradictory, and the excluded cases are needed:
* the degenerate square at the ratio `4 : 3` (`τ = 1`) satisfies the eight line equations and
  `(b + d)·3 = (h + f)·4` and `b² ≠ d²`; only `b² ≠ f²` fails, so `ratio_ne_4_3` is false without
  it;
* `tau_m_n` is attained at `τ = 1` and `τ = −1` (both forms equal `m² + n²`, a square for these
  three pairs) and, for `(7, 1)`, at `τ = 0`;
* `mss_point` applies to that degenerate square (`τ = 1`);
* its reflection `1 7 5 / 7 5 1 / 5 1 7` (roots `c b a / f e d / i (−h) g`), used by
  `ratio_ne_partner`, is again a magic square of squares, at the partner ratio `14 : 2 = 7 : 1`.
-/

open WeierstrassCurve WeierstrassCurve.Affine
open LeanRank.RankZeroGen LeanRank.CensusClosure

/-- The degenerate square `25 49 1 / 1 25 49 / 49 1 25` (roots `5 7 1 / 1 5 7 / 7 −1 5`, magic sum
`75`) at the ratio `4 : 3`: the eight line equations, the ratio, `b² ≠ d²`, and `b² = f²`. -/
example :
    (5 : ℤ) ^ 2 + 7 ^ 2 + 1 ^ 2 = 75 ∧ (1 : ℤ) ^ 2 + 5 ^ 2 + 7 ^ 2 = 75 ∧
    (7 : ℤ) ^ 2 + (-1) ^ 2 + 5 ^ 2 = 75 ∧ (5 : ℤ) ^ 2 + 1 ^ 2 + 7 ^ 2 = 75 ∧
    (7 : ℤ) ^ 2 + 5 ^ 2 + (-1) ^ 2 = 75 ∧ (1 : ℤ) ^ 2 + 7 ^ 2 + 5 ^ 2 = 75 ∧
    (5 : ℤ) ^ 2 + 5 ^ 2 + 5 ^ 2 = 75 ∧ (1 : ℤ) ^ 2 + 5 ^ 2 + 7 ^ 2 = 75 ∧
    ((7 : ℤ) + 1) * 3 = ((-1 : ℤ) + 7) * 4 ∧ (7 : ℤ) ^ 2 ≠ 1 ^ 2 ∧ (7 : ℤ) ^ 2 = 7 ^ 2 := by
  norm_num

/-- `tau_4_3` is attained at `τ = 1` and `τ = −1`. -/
example : ∃ r0 r8 : ℚ, r0 ^ 2 = 16 + 9 * (1 : ℚ) ^ 2 ∧ r8 ^ 2 = 9 + 16 * (1 : ℚ) ^ 2 :=
  ⟨5, 5, by norm_num, by norm_num⟩
example : ∃ r0 r8 : ℚ, r0 ^ 2 = 16 + 9 * (-1 : ℚ) ^ 2 ∧ r8 ^ 2 = 9 + 16 * (-1 : ℚ) ^ 2 :=
  ⟨5, 5, by norm_num, by norm_num⟩

/-- `tau_7_1` is attained at `τ = 0`. -/
example : ∃ r0 r8 : ℚ, r0 ^ 2 = 49 + 1 * (0 : ℚ) ^ 2 ∧ r8 ^ 2 = 1 + 49 * (0 : ℚ) ^ 2 :=
  ⟨7, 1, by norm_num, by norm_num⟩

/-- The theorems apply: from the hypotheses of `tau_12_5` at `τ = 1`, its conclusion. -/
example : (1 : ℚ) = 0 ∨ (1 : ℚ) = 1 ∨ (1 : ℚ) = -1 :=
  tau_12_5 (τ := 1) (r0 := 13) (r8 := 13) (by norm_num) (by norm_num)

/-- `E(ℚ)` for `(24, 7)` has 16 points. -/
example : Nat.card (curve (-331776) (-2401)).Point = 16 := card_24_7

/-- `mss_point` applies to the degenerate square at `4 : 3` (`τ = (h + f)/(b − d) = 6/6`). -/
example : ∃ τ r0 r8 : ℚ, τ = ((-1 : ℤ) + (7 : ℤ) : ℚ) / ((7 : ℤ) - (1 : ℤ) : ℚ) ∧
    r0 ^ 2 = ((4 : ℕ) : ℚ) ^ 2 + ((3 : ℕ) : ℚ) ^ 2 * τ ^ 2 ∧
    r8 ^ 2 = ((3 : ℕ) : ℚ) ^ 2 + ((4 : ℕ) : ℚ) ^ 2 * τ ^ 2 ∧
    (curve (-((4 : ℕ) : ℤ) ^ 4) (-((3 : ℕ) : ℤ) ^ 4)).Equation
      (((4 : ℕ) : ℚ) ^ 2 * ((3 : ℕ) : ℚ) ^ 2 * τ ^ 2)
      (((4 : ℕ) : ℚ) ^ 2 * ((3 : ℕ) : ℚ) ^ 2 * τ * r0 * r8) :=
  mss_point (m := 4) (n := 3) (by norm_num) (a := 5) (b := 7) (c := 1) (d := 1) (e := 5) (f := 7)
    (g := 7) (h := -1) (i := 5) (T := 75) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)

/-- The reflected degenerate square `1 49 25 / 49 25 1 / 25 1 49` (roots `1 7 5 / 7 5 1 / 5 1 7`,
i.e. `c b a / f e d / i (−h) g` of the square above) has `(b + d) : (h + f) = 14 : 2 = 7 : 1`. -/
example :
    (1 : ℤ) ^ 2 + 7 ^ 2 + 5 ^ 2 = 75 ∧ (7 : ℤ) ^ 2 + 5 ^ 2 + 1 ^ 2 = 75 ∧
    (5 : ℤ) ^ 2 + 1 ^ 2 + 7 ^ 2 = 75 ∧ (1 : ℤ) ^ 2 + 7 ^ 2 + 5 ^ 2 = 75 ∧
    (7 : ℤ) ^ 2 + 5 ^ 2 + 1 ^ 2 = 75 ∧ (5 : ℤ) ^ 2 + 1 ^ 2 + 7 ^ 2 = 75 ∧
    (1 : ℤ) ^ 2 + 5 ^ 2 + 7 ^ 2 = 75 ∧ (5 : ℤ) ^ 2 + 5 ^ 2 + 5 ^ 2 = 75 ∧
    ((7 : ℤ) + 7) * 1 = ((1 : ℤ) + 1) * 7 := by
  norm_num
