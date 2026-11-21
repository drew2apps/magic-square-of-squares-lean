import MagicSquareOfSquares.Defs

/-!
# Magic squares of squares: identities, zero entries and parity

## Main results

* `IsMSS.sum_eq`, `IsMSS.opposite`, `IsMSS.corner`: the linear identities, for squares.
* `IsMSS.eq_zero_of_exists_eq_zero`: if one entry is `0` then all nine are `0`.
* `IsMSS.parity`: the nine roots are all odd or all even.
* `IsMSS.odd`: in a primitive square all nine roots are odd.
-/

namespace MagicSquareOfSquares

variable {a b c d e f g h i T : ℤ}

namespace IsMSS

/-- Reflecting in the main diagonal. -/
theorem transpose (H : IsMSS a b c d e f g h i T) : IsMSS a d g b e h c f i T :=
  IsMagic.transpose H

/-- Rotating by a quarter turn. -/
theorem rotate (H : IsMSS a b c d e f g h i T) : IsMSS g d a h e b i f c T :=
  IsMagic.rotate H

/-- The magic sum is three times the centre entry. -/
theorem sum_eq (H : IsMSS a b c d e f g h i T) : T = 3 * e ^ 2 :=
  IsMagic.sum_eq H

/-- Opposite entries sum to twice the centre entry. -/
theorem opposite (H : IsMSS a b c d e f g h i T) :
    a ^ 2 + i ^ 2 = 2 * e ^ 2 ∧ c ^ 2 + g ^ 2 = 2 * e ^ 2 ∧
    b ^ 2 + h ^ 2 = 2 * e ^ 2 ∧ d ^ 2 + f ^ 2 = 2 * e ^ 2 :=
  IsMagic.opposite H

/-- Twice a corner entry is the sum of the two edge entries that do not touch it. -/
theorem corner (H : IsMSS a b c d e f g h i T) :
    2 * a ^ 2 = f ^ 2 + h ^ 2 ∧ 2 * c ^ 2 = d ^ 2 + h ^ 2 ∧
    2 * g ^ 2 = b ^ 2 + f ^ 2 ∧ 2 * i ^ 2 = b ^ 2 + d ^ 2 :=
  IsMagic.corner H

end IsMSS

/-- `x² = 2y²` has only the zero solution in integers. -/
theorem eq_zero_of_sq_eq_two_mul_sq {x y : ℤ} (hxy : x ^ 2 = 2 * y ^ 2) : x = 0 ∧ y = 0 := by
  by_cases hy : y = 0
  · subst hy
    refine ⟨?_, rfl⟩
    simpa using hxy
  · exfalso
    have hy' : (y : ℝ) ≠ 0 := by exact_mod_cast hy
    have h2 : ((x : ℝ) / y) ^ 2 = 2 := by
      rw [div_pow, div_eq_iff (pow_ne_zero 2 hy')]
      exact_mod_cast hxy
    have hs : Real.sqrt 2 = |(x : ℝ) / y| := by
      rw [← h2, Real.sqrt_sq_eq_abs]
    have hirr : Irrational |(x : ℝ) / y| := hs ▸ irrational_sqrt_two
    exact hirr ⟨|(x : ℚ) / y|, by push_cast; rfl⟩

/-- A sum of two integer squares is zero only if both are zero. -/
theorem eq_zero_of_sq_add_sq_eq_zero {x y : ℤ} (hxy : x ^ 2 + y ^ 2 = 0) : x = 0 ∧ y = 0 := by
  have hx : x ^ 2 = 0 := by nlinarith [sq_nonneg x, sq_nonneg y]
  have hy : y ^ 2 = 0 := by nlinarith [sq_nonneg x, sq_nonneg y]
  exact ⟨by simpa using hx, by simpa using hy⟩

namespace IsMSS

/-- If the centre is `0` then all nine roots are `0`. -/
theorem eq_zero_of_center_eq_zero (H : IsMSS a b c d e f g h i T) (he : e = 0) :
    a = 0 ∧ b = 0 ∧ c = 0 ∧ d = 0 ∧ e = 0 ∧ f = 0 ∧ g = 0 ∧ h = 0 ∧ i = 0 := by
  obtain ⟨h1, h2, h3, h4⟩ := H.opposite
  subst he
  simp only [ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, mul_zero] at h1 h2 h3 h4
  obtain ⟨ha, hi⟩ := eq_zero_of_sq_add_sq_eq_zero h1
  obtain ⟨hc, hg⟩ := eq_zero_of_sq_add_sq_eq_zero h2
  obtain ⟨hb, hh⟩ := eq_zero_of_sq_add_sq_eq_zero h3
  obtain ⟨hd, hf⟩ := eq_zero_of_sq_add_sq_eq_zero h4
  exact ⟨ha, hb, hc, hd, rfl, hf, hg, hh, hi⟩

/-- If one of the nine entries is `0` then all nine are `0`. So a magic square of squares
with distinct entries has no zero entry. -/
theorem eq_zero_of_exists_eq_zero (H : IsMSS a b c d e f g h i T)
    (h0 : a = 0 ∨ b = 0 ∨ c = 0 ∨ d = 0 ∨ e = 0 ∨ f = 0 ∨ g = 0 ∨ h = 0 ∨ i = 0) :
    a = 0 ∧ b = 0 ∧ c = 0 ∧ d = 0 ∧ e = 0 ∧ f = 0 ∧ g = 0 ∧ h = 0 ∧ i = 0 := by
  apply H.eq_zero_of_center_eq_zero
  obtain ⟨h1, h2, h3, h4⟩ := H.opposite
  have key : ∀ {x y : ℤ}, x ^ 2 + y ^ 2 = 2 * e ^ 2 → x = 0 → e = 0 := by
    intro x y hxy hx
    subst hx
    exact (eq_zero_of_sq_eq_two_mul_sq (x := y) (y := e) (by linarith)).2
  have key' : ∀ {x y : ℤ}, x ^ 2 + y ^ 2 = 2 * e ^ 2 → y = 0 → e = 0 := by
    intro x y hxy hy
    exact key (x := y) (y := x) (by linarith) hy
  rcases h0 with h0 | h0 | h0 | h0 | h0 | h0 | h0 | h0 | h0
  · exact key h1 h0
  · exact key h3 h0
  · exact key h2 h0
  · exact key h4 h0
  · exact h0
  · exact key' h4 h0
  · exact key' h2 h0
  · exact key' h3 h0
  · exact key' h1 h0

end IsMSS

/-- If `x² + y² = 2e²` then `x`, `y`, `e` are all odd or all even. -/
theorem parity_of_sq_add_sq {x y e : ℤ} (hxy : x ^ 2 + y ^ 2 = 2 * e ^ 2) :
    (Odd x ∧ Odd y ∧ Odd e) ∨ (Even x ∧ Even y ∧ Even e) := by
  rcases Int.even_or_odd' x with ⟨k, rfl | rfl⟩ <;>
  rcases Int.even_or_odd' y with ⟨m, rfl | rfl⟩ <;>
  rcases Int.even_or_odd' e with ⟨n, rfl | rfl⟩
  · right; simp
  · exfalso; ring_nf at hxy; omega
  · exfalso; ring_nf at hxy; omega
  · exfalso; ring_nf at hxy; omega
  · exfalso; ring_nf at hxy; omega
  · exfalso; ring_nf at hxy; omega
  · exfalso; ring_nf at hxy; omega
  · left; simp

namespace IsMSS

/-- The nine roots are all odd or all even. -/
theorem parity (H : IsMSS a b c d e f g h i T) :
    (Odd a ∧ Odd b ∧ Odd c ∧ Odd d ∧ Odd e ∧ Odd f ∧ Odd g ∧ Odd h ∧ Odd i) ∨
    (Even a ∧ Even b ∧ Even c ∧ Even d ∧ Even e ∧ Even f ∧ Even g ∧ Even h ∧ Even i) := by
  obtain ⟨h1, h2, h3, h4⟩ := H.opposite
  rcases Int.even_or_odd e with he | he
  · right
    have key : ∀ {y z : ℤ}, y ^ 2 + z ^ 2 = 2 * e ^ 2 → Even y ∧ Even z := by
      intro y z hyz
      rcases parity_of_sq_add_sq hyz with ⟨_, _, ho⟩ | ⟨hy, hz, _⟩
      · exact absurd he (Int.not_even_iff_odd.mpr ho)
      · exact ⟨hy, hz⟩
    exact ⟨(key h1).1, (key h3).1, (key h2).1, (key h4).1, he, (key h4).2, (key h2).2,
      (key h3).2, (key h1).2⟩
  · left
    have key : ∀ {y z : ℤ}, y ^ 2 + z ^ 2 = 2 * e ^ 2 → Odd y ∧ Odd z := by
      intro y z hyz
      rcases parity_of_sq_add_sq hyz with ⟨hy, hz, _⟩ | ⟨_, _, hev⟩
      · exact ⟨hy, hz⟩
      · exact absurd hev (Int.not_even_iff_odd.mpr he)
    exact ⟨(key h1).1, (key h3).1, (key h2).1, (key h4).1, he, (key h4).2, (key h2).2,
      (key h3).2, (key h1).2⟩

/-- In a primitive magic square of squares all nine roots are odd. -/
theorem odd (H : IsMSS a b c d e f g h i T) (hprim : Primitive a b c d e f g h i) :
    Odd a ∧ Odd b ∧ Odd c ∧ Odd d ∧ Odd e ∧ Odd f ∧ Odd g ∧ Odd h ∧ Odd i := by
  rcases H.parity with hodd | ⟨ha, hb, hc, hd, he, hf, hg, hh, hi⟩
  · exact hodd
  · exfalso
    have two : ∀ {x : ℤ}, Even x → ((2 : ℕ) : ℤ) ∣ x := fun hx => by
      exact_mod_cast even_iff_two_dvd.mp hx
    exact hprim 2 Nat.prime_two
      ⟨two ha, two hb, two hc, two hd, two he, two hf, two hg, two hh, two hi⟩

end IsMSS

end MagicSquareOfSquares
