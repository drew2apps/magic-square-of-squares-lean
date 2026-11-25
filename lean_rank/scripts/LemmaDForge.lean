import LeanRank.LemmaDReduction

/-!
Negative tests for `LeanRank/LemmaDReduction.lean` (written by an AI model; not read by a person).
Each `example` marked REFUSED must produce exactly one error; the ones marked ACCEPTED must compile.
-/

namespace MagicSquareOfSquares.Family

-- ACCEPTED: dropping the abscissa 12 + 8√2 from `E1Points` gives a false hypothesis (the image of (−1, 4)).
example (h : ∀ a b c e : ℚ, E1Q2 a b c e →
    (a = 0 ∧ b = 0) ∨ (a = 4 ∧ b = 0) ∨ (a = -4 ∧ b = 0) ∨ (a = -12 ∧ b = -8)) : False := by
  rcases h 12 8 80 56 E1_points_on.2.2.2.2.2.1 with ⟨h1, -⟩ | ⟨h1, -⟩ | ⟨h1, -⟩ | ⟨h1, -⟩ <;> norm_num at h1

-- ACCEPTED: dropping the abscissa 4 is false too (the image of (1, 4)).
example (h : ∀ a b c e : ℚ, E1Q2 a b c e →
    (a = 0 ∧ b = 0) ∨ (a = -4 ∧ b = 0) ∨ (a = 12 ∧ b = 8) ∨ (a = -12 ∧ b = -8)) : False := by
  rcases h 4 0 0 8 E1_points_on.2.1 with ⟨h1, -⟩ | ⟨h1, -⟩ | ⟨h1, -⟩ | ⟨h1, -⟩ <;> norm_num at h1

-- REFUSED (forged point): (4, 4√2) is not on E1.
example : E1Q2 4 0 0 4 := by
  have hr : (√2 : ℝ) ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  unfold E1Q2 E1eq
  push_cast
  linear_combination 96 * hr

-- REFUSED (forged map, 4w in place of 8w, the factor-2 slip of the LMFDB remark).
example {r x w : ℝ} (hr : r ^ 2 = 2) (hz : x - 1 + r ≠ 0)
    (hw : w ^ 2 = 2 * (x ^ 2 + 1) * (x ^ 4 - 2 * x ^ 3 + 2 * x ^ 2 + 2 * x + 1)) :
    E1eq (2 + r) (4 * (x ^ 2 + 1) / (x - 1 + r) ^ 2) (4 * w / (x - 1 + r) ^ 3) :=
  d2_to_E1 hr hz hw

-- REFUSED (forged conclusion x = 1 only; (−1, ±4) is a point): the proof of `lemmaD_of_E1Points` with
-- `x = 1` in place of `x = 1 ∨ x = −1` must fail in the case X = 12 + 8√2.
example (hE : E1Points) : ∀ x w : ℚ,
    w ^ 2 = 2 * (x ^ 2 + 1) * (x ^ 4 - 2 * x ^ 3 + 2 * x ^ 2 + 2 * x + 1) → x = 1 := by
  intro x w hw
  have hD : (x - 1) ^ 2 - 2 ≠ 0 := sub_ne_zero.mpr (rat_sq_ne_two (x - 1))
  have hD2 : ((x - 1) ^ 2 - 2) ^ 2 ≠ 0 := pow_ne_zero 2 hD
  rcases hE _ _ _ _ (pullback hw) with ⟨h1, h2⟩ | ⟨h1, h2⟩ | ⟨h1, h2⟩ | ⟨h1, h2⟩ | ⟨h1, h2⟩ <;>
    rw [div_eq_iff hD2] at h1 h2
  · rw [zero_mul] at h1
    exact absurd h1 (by positivity)
  · have h : (x - 1) * (x ^ 2 + 1) = 0 := by linear_combination (-1/8 : ℚ) * h2
    rcases mul_eq_zero.1 h with h | h
    · linarith
    · exact absurd h (by positivity)
  · have h : (x - 1) * (x ^ 2 + 1) = 0 := by linear_combination (-1/8 : ℚ) * h2
    rcases mul_eq_zero.1 h with h | h
    · linarith
    · exact absurd h (by positivity)
  · nlinarith [sq_nonneg (x - 2), sq_nonneg x, sq_nonneg (x + 1)]
  · exact absurd (by linear_combination (-x^3/8 + 5*x^2/8 - 11*x/32 - 9/16) * h1
      + (x^3/4 - 7*x^2/8 + 9*x/16 + 29/32) * h2 : (1 : ℚ) = 0) one_ne_zero

-- REFUSED (forged pull-back: (−1, 4) sent to abscissa 4).
example : 4 * ((-1:ℚ) ^ 2 + 1) * ((-1 - 1) ^ 2 + 2) / ((-1 - 1) ^ 2 - 2) ^ 2 = 4 := by norm_num

-- REFUSED (forged 2-torsion: a second point of order 2 at X = 4).
example : E1Q2 4 0 0 0 := by
  have hr : (√2 : ℝ) ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  unfold E1Q2 E1eq
  push_cast
  linear_combination 96 * hr

-- ACCEPTED (control): the honest statements.
example (hE : E1Points) : LemmaD := lemmaD_of_E1Points hE
example (a b : ℚ) (h : E1Q2 a b 0 0) : a = 0 := (E1_two_torsion a b h).1

end MagicSquareOfSquares.Family
