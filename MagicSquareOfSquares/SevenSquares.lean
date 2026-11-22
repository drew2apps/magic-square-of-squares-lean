import Mathlib.Tactic

/-!
# Two facts behind the search for squares with seven square entries

The program `scripts/search/seven.c` is not verified in Lean, and its result is not a Lean theorem.
This file checks the two pieces of algebra it rests on. (Bremner's square itself is in
`Examples.lean`.)

* `full_line_legs`: if `x² − d = p²` and `x² + d = q²`, there are integers `u`, `v` with
  `p = u − v`, `q = u + v`, `u² + v² = x²` and `d = 2uv`. So a line through the centre `x²` with
  both ends square is a Pythagorean triangle with hypotenuse `x`, and its deviation is twice the
  product of the legs. (`full_line_triangle` is the same without the halving.)
* The four lines through the centre have deviations `s`, `t`, `s + t`, `s − t`, `0 < t < s`.
  `TwoLines s t a b r1 r2` lists the ways in which `a < b` can be two of them and `r1`, `r2` the
  other two. `two_lines_exhaust` proves that the list is complete: if `a < b`, `r1`, `r2` are the
  four deviations in some order, then `TwoLines` holds for `(r1, r2)` or for `(r2, r1)`.
  `rest_of_two_lines` gives the six formulas of `seven.c` for `r1`, `r2` in terms of `a`, `b`
  (`SixCases`), and `six_cases` combines the two statements. The first version of `seven.c` had
  four of the six formulas.
-/

namespace MagicSquareOfSquares

theorem full_line_triangle (x d p q : ℤ) (h1 : x ^ 2 - d = p ^ 2) (h2 : x ^ 2 + d = q ^ 2) :
    (q + p) ^ 2 + (q - p) ^ 2 = (2 * x) ^ 2 ∧ 2 * d = (q + p) * (q - p) := by
  constructor
  · linear_combination (-2) * h1 + (-2) * h2
  · linear_combination h2 - h1

theorem full_line_legs (x d p q : ℤ) (h1 : x ^ 2 - d = p ^ 2) (h2 : x ^ 2 + d = q ^ 2) :
    ∃ u v : ℤ, p = u - v ∧ q = u + v ∧ u ^ 2 + v ^ 2 = x ^ 2 ∧ d = 2 * u * v := by
  have h4 : Even ((q + p) ^ 2) := ⟨x ^ 2 + p * q, by linear_combination (-1) * h1 + (-1) * h2⟩
  have h5 : Even (q + p) := (Int.even_pow.mp h4).1
  obtain ⟨u, hu⟩ := h5
  have hq : q = 2 * u - p := by linarith
  subst hq
  refine ⟨u, u - p, by ring, by ring, ?_, ?_⟩
  · nlinarith [h1, h2]
  · nlinarith [h1, h2]

/-- The ways in which `a`, `b` can be two of the deviations `s`, `t`, `s + t`, `s − t` (with
`a < b` when `0 < t < s`) and `r1`, `r2` the other two, in a fixed order. -/
def TwoLines (s t a b r1 r2 : ℤ) : Prop :=
  (a = t ∧ b = s ∧ r1 = s + t ∧ r2 = s - t) ∨
  (a = s ∧ b = s + t ∧ r1 = t ∧ r2 = s - t) ∨
  (a = t ∧ b = s + t ∧ r1 = s ∧ r2 = s - t) ∨
  (a = s - t ∧ b = s ∧ r1 = t ∧ r2 = s + t) ∨
  (a = s - t ∧ b = s + t ∧ r1 = s ∧ r2 = t) ∨
  (a = t ∧ b = s - t ∧ r1 = s ∧ r2 = s + t) ∨
  (a = s - t ∧ b = t ∧ r1 = s ∧ r2 = s + t)

/-- The six formulas of `seven.c` for the other two deviations. -/
def SixCases (a b r1 r2 : ℤ) : Prop :=
  (r1 = a + b ∧ r2 = b - a) ∨
  (r1 = b - a ∧ (r2 = 2 * a - b ∨ r2 = b - 2 * a)) ∨
  (r1 = b - a ∧ r2 = 2 * b - a) ∨
  (2 * r1 = a + b ∧ 2 * r2 = b - a) ∨
  (r1 = a + b ∧ r2 = 2 * a + b) ∨
  (r1 = a + b ∧ r2 = a + 2 * b)

theorem rest_of_two_lines {s t a b r1 r2 : ℤ} (h : TwoLines s t a b r1 r2) :
    SixCases a b r1 r2 := by
  unfold TwoLines at h
  unfold SixCases
  rcases h with ⟨rfl, rfl, rfl, rfl⟩ | ⟨rfl, rfl, rfl, rfl⟩ | ⟨rfl, rfl, rfl, rfl⟩ |
    ⟨rfl, rfl, rfl, rfl⟩ | ⟨rfl, rfl, rfl, rfl⟩ | ⟨rfl, rfl, rfl, rfl⟩ | ⟨rfl, rfl, rfl, rfl⟩
  · exact Or.inl ⟨by ring, by ring⟩
  · exact Or.inr (Or.inl ⟨by ring, Or.inl (by ring)⟩)
  · exact Or.inr (Or.inl ⟨by ring, Or.inr (by ring)⟩)
  · exact Or.inr (Or.inr (Or.inl ⟨by ring, by ring⟩))
  · exact Or.inr (Or.inr (Or.inr (Or.inl ⟨by ring, by ring⟩)))
  · exact Or.inr (Or.inr (Or.inr (Or.inr (Or.inl ⟨by ring, by ring⟩))))
  · exact Or.inr (Or.inr (Or.inr (Or.inr (Or.inr ⟨by ring, by ring⟩))))

set_option maxHeartbeats 4000000 in
/-- The list `TwoLines` is complete: if `0 < t < s`, `a < b`, and `a`, `b`, `r1`, `r2` are
pairwise distinct members of `{s, t, s + t, s − t}`, then `TwoLines` holds for `(r1, r2)` or for
`(r2, r1)`. -/
theorem two_lines_exhaust (s t a b r1 r2 : ℤ) (ht : 0 < t) (hst : t < s) (hab : a < b)
    (ha : a = s ∨ a = t ∨ a = s + t ∨ a = s - t) (hb : b = s ∨ b = t ∨ b = s + t ∨ b = s - t)
    (h1 : r1 = s ∨ r1 = t ∨ r1 = s + t ∨ r1 = s - t)
    (h2 : r2 = s ∨ r2 = t ∨ r2 = s + t ∨ r2 = s - t)
    (d1 : a ≠ r1) (d2 : a ≠ r2) (d3 : b ≠ r1) (d4 : b ≠ r2) (d5 : r1 ≠ r2) :
    TwoLines s t a b r1 r2 ∨ TwoLines s t a b r2 r1 := by
  unfold TwoLines
  rcases ha with ha | ha | ha | ha <;> rcases hb with hb | hb | hb | hb <;>
    rcases h1 with h1 | h1 | h1 | h1 <;> rcases h2 with h2 | h2 | h2 | h2 <;>
    subst ha hb h1 h2 <;> first | (exfalso; omega) | simp

/-- For two lines `a < b` out of the four, the other two deviations are given by one of the six
formulas, in one of the two orders. -/
theorem six_cases (s t a b r1 r2 : ℤ) (ht : 0 < t) (hst : t < s) (hab : a < b)
    (ha : a = s ∨ a = t ∨ a = s + t ∨ a = s - t) (hb : b = s ∨ b = t ∨ b = s + t ∨ b = s - t)
    (h1 : r1 = s ∨ r1 = t ∨ r1 = s + t ∨ r1 = s - t)
    (h2 : r2 = s ∨ r2 = t ∨ r2 = s + t ∨ r2 = s - t)
    (d1 : a ≠ r1) (d2 : a ≠ r2) (d3 : b ≠ r1) (d4 : b ≠ r2) (d5 : r1 ≠ r2) :
    SixCases a b r1 r2 ∨ SixCases a b r2 r1 :=
  (two_lines_exhaust s t a b r1 r2 ht hst hab ha hb h1 h2 d1 d2 d3 d4 d5).elim
    (fun h => Or.inl (rest_of_two_lines h)) (fun h => Or.inr (rest_of_two_lines h))

end MagicSquareOfSquares
