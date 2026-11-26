import LeanRank.E1Descent

/-!
# Forged claims against `E1Descent.lean` (each must FAIL)

Written by an AI model (Claude, made by Anthropic) on 26 November 2026 (UTC); not read by a person.

Expected: the 4 forged claims (marked FORGED) are all refused, with exactly **5 errors**: lines 20, 25, 35
once each, line 30 twice (forged claim 3 fails at two of its hypotheses: `1² = ε²` and `Odd (1 : ℤ√2).im`
are both false). The last theorem is honest and must be accepted: it shows that the second alternative of
`DescentW` (`b/x` a square) cannot be dropped.
-/

namespace MagicSquareOfSquares.E1Descent

open Zsqrtd MagicSquareOfSquares.Zsqrt2 MagicSquareOfSquares.Family MagicSquareOfSquares.E1Rank

/-- FORGED 1: the kept class `3` of `E1'` "excluded" at the first real place (`σ₁(3) = 3 > 0`). -/
example {n r s t : ℤ√2} (hr : r ≠ 0) (ht : t ≠ 0) : n ^ 2 ≠ r ^ 2 * Fq A2 B2 three s t :=
  real_excl sig1 sig1_inj hr ht (by rw [sig1_three]; norm_num) sig1_q

/-- FORGED 2: the class `2 = √2·√2 ≡ 1` of `E1` "excluded" by the odd-valuation test. -/
example {n r s t : ℤ√2} (hr : r ≠ 0) (hst : ¬(pi2 ∣ s ∧ pi2 ∣ t)) :
    n ^ 2 ≠ r ^ 2 * Fq A1 B1 (pi2 * pi2) s t :=
  V1 (by rw [pi2_dvd_iff]; decide) hr hst

/-- FORGED 3: the kept class `1` of `E1` "excluded" by the mod-2 test (`u = 1`). -/
example {n r s t : ℤ√2} (hr : r ≠ 0) (hst : ¬(pi2 ∣ s ∧ pi2 ∣ t)) :
    n ^ 2 ≠ r ^ 2 * Fq A1 B1 (pi2 ^ 4 * 1) s t :=
  M1 (by decide) ⟨0, rfl⟩ (by rw [pi2_dvd_iff]; decide) hr hst

/-- FORGED 4: the descent statement for `E1` without its second alternative. -/
example : ∀ x y : Q2, y ^ 2 = x * (x ^ 2 + 8 * (2 + s2) * x + -8 * (2 + s2) ^ 2) → x ≠ 0 →
    ∃ u, x = u ^ 2 :=
  fun x y h hx => (descentE1 x y h hx).resolve_right (by simp)

/-- Honest (accepted): `(−4, 16 + 8√2)` is a point of `E1` whose abscissa is not a square, so the second
alternative of `DescentW` is needed. -/
theorem descentE1_sharp : ∃ x y : Q2, y ^ 2 = x * (x ^ 2 + 8 * (2 + s2) * x + -8 * (2 + s2) ^ 2) ∧
    x ≠ 0 ∧ ¬∃ u : Q2, x = u ^ 2 := by
  refine ⟨-4, 16 + 8 * s2, ?_, by norm_num, ?_⟩
  · linear_combination (32 : Q2) * s2_mul_s2
  · rintro ⟨u, hu⟩
    have h := congrArg (algebraMap Q2 ℝ) hu
    rw [map_neg, map_ofNat, map_pow] at h
    nlinarith [sq_nonneg (algebraMap Q2 ℝ u)]

end MagicSquareOfSquares.E1Descent

#print axioms MagicSquareOfSquares.E1Descent.descentE1_sharp
