import MagicSquareOfSquares.TwoDescentCertificates

/-!
Checks for `TwoDescentMap.lean` and `TwoDescentCertificates.lean`, not part of the library (run with
`lake env lean scripts/TwoDescentChecks.lean`). Each `example` states, as a theorem, that a forged certificate
cannot be built. Written by an AI model; not read by a person.

The forged certificates themselves (a wrong value of `H(7, 1)`, a bracket for the square `169 = H(4, 1)`, a wrong
root `43`, `H(7, 1)` claimed negative, Theorem 2 at the square-type ratio `(12, 5)`, the `E⁻` certificate of
`(7, 1)` used for `E⁺`, `control_bracket` flipped, the square `16` bracketed) were each written once and refused by
Lean (25 November 2026, slice a1); they are not kept here because a file that must fail cannot be compiled in CI.
-/

namespace MagicSquareOfSquares.TwoDescent

open MagicSquareOfSquares.G1 WeierstrassCurve WeierstrassCurve.Affine

/-- `H(4, 1) = 169 = 13²` (`UNIFORM_RANK.md` section 6): no bracket exists, so `nsq` cannot certify it. -/
example : genH 4 1 = 169 ∧ ¬ ∃ s : ℕ, s * s < 169 ∧ 169 < (s + 1) * (s + 1) := by
  refine ⟨by unfold genH; norm_num, ?_⟩
  rintro ⟨s, h1, h2⟩
  rcases le_or_gt s 12 with h | h
  · have := Nat.mul_self_le_mul_self (show s + 1 ≤ 13 by omega)
    omega
  · have := Nat.mul_self_le_mul_self (show 13 ≤ s from h)
    omega

/-- At the square-type ratio `(12, 5)` there is no `k` with `m² + n² = 2k²`, so `theorem2_of` does not apply
(`UNIFORM_RANK.md`: there is no point `Q` on the square type). -/
example : ¬ ∃ k : ℚ, (12 : ℚ) ^ 2 + 5 ^ 2 = 2 * k ^ 2 := by
  rintro ⟨k, hk⟩
  exact not_isSquare_of_bracket 338 18 (by norm_num) (by norm_num)
    ⟨2 * k, by push_cast; linear_combination 2 * hk⟩

/-- On the double `P = 2R` of `(7, 1)`, both tests return squares (`control_7_1`). -/
example (hR : (Eplus 7 1).toAffine.Nonsingular (xR 7 1) (yR 7 1)) :
    IsSquare (alpha (eP 7 1) (aP 7 1) (bP 7 1) (Point.some _ _ hR + Point.some _ _ hR)) := by
  obtain ⟨hP, h1, h2, -, -⟩ := control_7_1 hR
  rw [h1, h2]
  exact ⟨300, rfl⟩

/-- The same from the general theorem: `α(W + W)` is always a square. -/
example (W : (Eplus 7 1).toAffine.Point) : IsSquare (alpha (eP 7 1) (aP 7 1) (bP 7 1) (W + W)) :=
  alpha_double (by rw [qe_P]; unfold genH; norm_num) W

end MagicSquareOfSquares.TwoDescent
