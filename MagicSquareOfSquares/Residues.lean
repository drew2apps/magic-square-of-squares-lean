import MagicSquareOfSquares.Primes

/-!
# Magic squares of squares: residues modulo 24 and 72

## Main results

* `IsMSS.not_three_dvd`: in a primitive square no root is divisible by `3`.
* `IsMSS.sq_emod_24`: in a primitive square every entry is `≡ 1 (mod 24)`.
* `IsMSS.sum_emod_72`: in a primitive square the magic sum is `≡ 3 (mod 72)`.
-/

namespace MagicSquareOfSquares

variable {a b c d e f g h i T : ℤ}

/-- The square of an integer coprime to `6` is `≡ 1 (mod 24)`. -/
theorem sq_emod_24_of_odd_of_not_three_dvd {x : ℤ} (h2 : Odd x) (h3 : ¬(3 : ℤ) ∣ x) :
    x ^ 2 % 24 = 1 := by
  obtain ⟨k, rfl⟩ := h2
  have h : (2 * k + 1) % 24 = 1 ∨ (2 * k + 1) % 24 = 5 ∨ (2 * k + 1) % 24 = 7 ∨
      (2 * k + 1) % 24 = 11 ∨ (2 * k + 1) % 24 = 13 ∨ (2 * k + 1) % 24 = 17 ∨
      (2 * k + 1) % 24 = 19 ∨ (2 * k + 1) % 24 = 23 := by omega
  rw [sq, Int.mul_emod]
  rcases h with h | h | h | h | h | h | h | h <;> rw [h] <;> norm_num

namespace IsMSS

/-- In a primitive magic square of squares no root is divisible by `3`. -/
theorem not_three_dvd (H : IsMSS a b c d e f g h i T) (hprim : Primitive a b c d e f g h i) :
    ¬(3 : ℤ) ∣ a ∧ ¬(3 : ℤ) ∣ b ∧ ¬(3 : ℤ) ∣ c ∧ ¬(3 : ℤ) ∣ d ∧ ¬(3 : ℤ) ∣ e ∧
    ¬(3 : ℤ) ∣ f ∧ ¬(3 : ℤ) ∣ g ∧ ¬(3 : ℤ) ∣ h ∧ ¬(3 : ℤ) ∣ i := by
  have : Fact (Nat.Prime 3) := ⟨Nat.prime_three⟩
  have key := H.not_dvd_of_three_mod_eight hprim (p := 3) (by norm_num)
  exact_mod_cast key

/-- In a primitive magic square of squares every entry is `≡ 1 (mod 24)`. -/
theorem sq_emod_24 (H : IsMSS a b c d e f g h i T) (hprim : Primitive a b c d e f g h i) :
    a ^ 2 % 24 = 1 ∧ b ^ 2 % 24 = 1 ∧ c ^ 2 % 24 = 1 ∧ d ^ 2 % 24 = 1 ∧ e ^ 2 % 24 = 1 ∧
    f ^ 2 % 24 = 1 ∧ g ^ 2 % 24 = 1 ∧ h ^ 2 % 24 = 1 ∧ i ^ 2 % 24 = 1 := by
  obtain ⟨oa, ob, oc, od, oe, of, og, oh, oi⟩ := H.odd hprim
  obtain ⟨ta, tb, tc, td, te, tf, tg, th, ti⟩ := H.not_three_dvd hprim
  exact ⟨sq_emod_24_of_odd_of_not_three_dvd oa ta, sq_emod_24_of_odd_of_not_three_dvd ob tb,
    sq_emod_24_of_odd_of_not_three_dvd oc tc, sq_emod_24_of_odd_of_not_three_dvd od td,
    sq_emod_24_of_odd_of_not_three_dvd oe te, sq_emod_24_of_odd_of_not_three_dvd of tf,
    sq_emod_24_of_odd_of_not_three_dvd og tg, sq_emod_24_of_odd_of_not_three_dvd oh th,
    sq_emod_24_of_odd_of_not_three_dvd oi ti⟩

/-- In a primitive magic square of squares the magic sum is `≡ 3 (mod 72)`. -/
theorem sum_emod_72 (H : IsMSS a b c d e f g h i T) (hprim : Primitive a b c d e f g h i) :
    T % 72 = 3 := by
  have he := (H.sq_emod_24 hprim).2.2.2.2.1
  have hT := H.sum_eq
  omega

end IsMSS

end MagicSquareOfSquares
