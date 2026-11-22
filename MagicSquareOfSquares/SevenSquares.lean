import MagicSquareOfSquares.Defs
import Mathlib.Tactic

/-!
# Squares with seven square entries

* Bremner's square (1999), the only known 3×3 magic square with seven square entries, is a magic
  square (`bremner_isMagic`), and its two other entries are not squares (`not_isSquare_360721`,
  `not_isSquare_222121`).
* Two facts behind the search program `scripts/search/seven.c`. That program is not verified in
  Lean, and its result (no other such square with centre root up to `10⁹`) is not a Lean theorem.
  - `full_line_triangle`: if `x² − d` and `x² + d` are both squares, `p²` and `q²`, then
    `(q + p, q − p, 2x)` is a Pythagorean triple and `2d = (q + p)(q − p)`.
  - `rest_of_two_lines`: the four lines through the centre have deviations `s`, `t`, `s + t`,
    `s − t` with `0 < t < s`. If two of them are `a`, `b`, the other two `r1`, `r2` are given by one
    of six formulas. These are the six cases of `seven.c`; its first version had only four.
-/

namespace MagicSquareOfSquares

/-- Bremner's square: rows `373² 289² 565²`, `360721 425² 23²`, `205² 527² 222121`. -/
theorem bremner_isMagic :
    IsMagic (373 ^ 2) (289 ^ 2) (565 ^ 2) 360721 (425 ^ 2) (23 ^ 2) (205 ^ 2) (527 ^ 2) 222121
      541875 := by
  refine IsMagic.mk ?_ ?_ ?_ ?_ ?_ ?_ ?_ ?_ <;> norm_num

theorem not_isSquare_360721 : ¬ IsSquare (360721 : ℕ) := by
  rintro ⟨r, hr⟩
  have h1 : 600 < r := by
    by_contra h
    have h' : r ≤ 600 := Nat.le_of_not_lt h
    nlinarith
  have h2 : r < 601 := by
    by_contra h
    have h' : 601 ≤ r := Nat.le_of_not_lt h
    nlinarith
  omega

theorem not_isSquare_222121 : ¬ IsSquare (222121 : ℕ) := by
  rintro ⟨r, hr⟩
  have h1 : 471 < r := by
    by_contra h
    have h' : r ≤ 471 := Nat.le_of_not_lt h
    nlinarith
  have h2 : r < 472 := by
    by_contra h
    have h' : 472 ≤ r := Nat.le_of_not_lt h
    nlinarith
  omega

/-- A line through the centre `x²` with both ends square is a Pythagorean triangle with
hypotenuse `2x` (halve it: the legs are `(q ± p)/2` and the hypotenuse is `x`). -/
theorem full_line_triangle (x d p q : ℤ) (h1 : x ^ 2 - d = p ^ 2) (h2 : x ^ 2 + d = q ^ 2) :
    (q + p) ^ 2 + (q - p) ^ 2 = (2 * x) ^ 2 ∧ 2 * d = (q + p) * (q - p) := by
  constructor
  · linear_combination (-2) * h1 + (-2) * h2
  · linear_combination h2 - h1

/-- The six cases for two lines out of four. The hypothesis lists the ways of choosing two of the
deviations `s`, `t`, `s + t`, `s − t` as `a`, `b` (for `0 < t < s` and `a < b` these are all the
ways), with `r1`, `r2` the other two. The identity itself needs no inequality. -/
theorem rest_of_two_lines (s t a b r1 r2 : ℤ)
    (h : (a = t ∧ b = s ∧ r1 = s + t ∧ r2 = s - t) ∨
         (a = s ∧ b = s + t ∧ r1 = t ∧ r2 = s - t) ∨
         (a = t ∧ b = s + t ∧ r1 = s ∧ r2 = s - t) ∨
         (a = s - t ∧ b = s ∧ r1 = t ∧ r2 = s + t) ∨
         (a = s - t ∧ b = s + t ∧ r1 = s ∧ r2 = t) ∨
         (a = t ∧ b = s - t ∧ r1 = s ∧ r2 = s + t) ∨
         (a = s - t ∧ b = t ∧ r1 = s ∧ r2 = s + t)) :
    (r1 = a + b ∧ r2 = b - a) ∨
    (r1 = b - a ∧ (r2 = 2 * a - b ∨ r2 = b - 2 * a)) ∨
    (r1 = b - a ∧ r2 = 2 * b - a) ∨
    (2 * r1 = a + b ∧ 2 * r2 = b - a) ∨
    (r1 = a + b ∧ r2 = 2 * a + b) ∨
    (r1 = a + b ∧ r2 = a + 2 * b) := by
  rcases h with ⟨rfl, rfl, rfl, rfl⟩ | ⟨rfl, rfl, rfl, rfl⟩ | ⟨rfl, rfl, rfl, rfl⟩ |
    ⟨rfl, rfl, rfl, rfl⟩ | ⟨rfl, rfl, rfl, rfl⟩ | ⟨rfl, rfl, rfl, rfl⟩ | ⟨rfl, rfl, rfl, rfl⟩
  · exact Or.inl ⟨by ring, by ring⟩
  · exact Or.inr (Or.inl ⟨by ring, Or.inl (by ring)⟩)
  · exact Or.inr (Or.inl ⟨by ring, Or.inr (by ring)⟩)
  · exact Or.inr (Or.inr (Or.inl ⟨by ring, by ring⟩))
  · exact Or.inr (Or.inr (Or.inr (Or.inl ⟨by ring, by ring⟩)))
  · exact Or.inr (Or.inr (Or.inr (Or.inr (Or.inl ⟨by ring, by ring⟩))))
  · exact Or.inr (Or.inr (Or.inr (Or.inr (Or.inr ⟨by ring, by ring⟩))))

end MagicSquareOfSquares
