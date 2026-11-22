import MagicSquareOfSquares.G1Norm
import MagicSquareOfSquares.Primes

/-!
# The centre and the number `m² + n²`

Let `m : n` be the ratio `(b + d) : (h + f)` of a magic square of squares in lowest terms and `e`
the centre root.

* `IsMSS.norm_dvd_center`: `m² + n²` divides `8e²`. Indeed `b + d = m·u`, `h + f = n·u`,
  `b − d = n·v`, `f − h = m·v`, and `8e² = (m² + n²)(u² + v²)`.
* `IsMSS.norm_isSquare_or_prime_dvd_center`: either `m² + n²` is a square or twice a square (the
  ratios for which the G1 system has a degenerate rational point), or some prime `p ≡ 1 (mod 8)`
  divides `m² + n²` to an odd power, and this prime divides the centre root.

So a magic square of squares whose centre root has no prime factor `≡ 1 (mod 8)` has
`m² + n²` a square or twice a square. The proof combines `norm_dvd_center` with the statement on
primes `≡ 5 (mod 8)` (`G1Primes.lean`) and the fact that no prime `≡ 3 (mod 4)` divides a sum of
two coprime squares.
-/

namespace MagicSquareOfSquares

open G1

/-- A positive natural number in which every odd prime has even valuation is a square or twice
a square. -/
theorem isSquare_or_two_mul_of_even_val {N : ℕ} (hN : N ≠ 0)
    (h : ∀ p : ℕ, p.Prime → p ≠ 2 → Even (padicValNat p N)) : IsSquare N ∨ IsSquare (2 * N) := by
  obtain ⟨a, b, hab, hsq⟩ := Nat.sq_mul_squarefree N
  have ha0 : a ≠ 0 := by
    rintro rfl
    exact hN (by rw [← hab, mul_zero])
  have hb0 : b ≠ 0 := by
    rintro rfl
    exact hN (by rw [← hab]; simp)
  have hprime : ∀ p : ℕ, p.Prime → p ∣ a → p = 2 := by
    intro p hp hpa
    by_contra hp2
    have : Fact p.Prime := ⟨hp⟩
    have hev := h p hp hp2
    rw [← hab, sq, padicValNat.mul (mul_ne_zero hb0 hb0) ha0, padicValNat.mul hb0 hb0] at hev
    have hle : padicValNat p a ≤ 1 := by
      by_contra hgt
      have h2 : 2 ≤ padicValNat p a := by omega
      have hd : p ^ 2 ∣ a := (padicValNat_dvd_iff_le ha0).mpr h2
      rw [sq] at hd
      exact hp.one_lt.ne' (Nat.isUnit_iff.mp (hsq p hd))
    have hge : 1 ≤ padicValNat p a := one_le_padicValNat_of_dvd ha0 hpa
    obtain ⟨k, hk⟩ := hev
    omega
  have ha : a = 1 ∨ a = 2 := by
    by_cases ha1 : a = 1
    · exact Or.inl ha1
    · right
      have hmf := hprime a.minFac (Nat.minFac_prime ha1) (Nat.minFac_dvd a)
      have h2 : 2 ∣ a := hmf ▸ Nat.minFac_dvd a
      obtain ⟨a', rfl⟩ := h2
      by_cases ha'1 : a' = 1
      · omega
      · exfalso
        have hq := hprime a'.minFac (Nat.minFac_prime ha'1)
          (Dvd.dvd.mul_left (Nat.minFac_dvd a') 2)
        have h2' : 2 ∣ a' := hq ▸ Nat.minFac_dvd a'
        obtain ⟨a'', rfl⟩ := h2'
        exact absurd (Nat.isUnit_iff.mp (hsq 2 ⟨a'', by ring⟩)) (by norm_num)
  rcases ha with rfl | rfl
  · exact Or.inl ⟨b, by rw [← hab]; ring⟩
  · exact Or.inr ⟨2 * b, by rw [← hab]; ring⟩

/-- In a magic square of squares with `(b + d)·n = (h + f)·m` and `m`, `n` coprime, the number
`m² + n²` divides `8e²`, where `e` is the centre root. -/
theorem IsMSS.norm_dvd_center {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T)
    {m n : ℕ} (hcop : Nat.Coprime m n) (hmn : (b + d) * n = (h + f) * m) (hne : h + f ≠ 0) :
    ((m : ℤ) ^ 2 + (n : ℤ) ^ 2) ∣ 8 * e ^ 2 := by
  obtain ⟨-, -, hbh, hdf⟩ := IsMagic.opposite
    (show IsMagic (a ^ 2) (b ^ 2) (c ^ 2) (d ^ 2) (e ^ 2) (f ^ 2) (g ^ 2) (h ^ 2) (i ^ 2) T from hs)
  rcases Nat.eq_zero_or_pos m with hm0 | hm0
  · subst hm0
    have hn1 : n = 1 := by simpa using hcop
    subst hn1
    simp
  · have hmz : (m : ℤ) ≠ 0 := by exact_mod_cast hm0.ne'
    have hnz : (n : ℤ) ≠ 0 := by
      intro h0
      rw [h0, mul_zero] at hmn
      exact hne ((mul_eq_zero.mp hmn.symm).resolve_right hmz)
    have hc : IsCoprime (m : ℤ) (n : ℤ) := Nat.isCoprime_iff_coprime.mpr hcop
    obtain ⟨u, hu⟩ : (m : ℤ) ∣ b + d := hc.dvd_of_dvd_mul_right ⟨h + f, by rw [hmn]; ring⟩
    have hY : h + f = n * u := by
      have h1 : (m : ℤ) * (h + f) = m * (n * u) := by
        linear_combination (-1 : ℤ) * hmn + (n : ℤ) * hu
      exact mul_left_cancel₀ hmz h1
    have hu0 : u ≠ 0 := by
      rintro rfl
      exact hne (by rw [hY, mul_zero])
    have hZW : (b - d) * m = (f - h) * n := by
      have h1 : u * ((b - d) * m) = u * ((f - h) * n) := by
        linear_combination hbh - hdf + (-(b - d)) * hu + (f - h) * hY
      exact mul_left_cancel₀ hu0 h1
    obtain ⟨v, hv⟩ : (n : ℤ) ∣ b - d := hc.symm.dvd_of_dvd_mul_right ⟨f - h, by rw [hZW]; ring⟩
    have hW : f - h = m * v := by
      have h1 : (n : ℤ) * (f - h) = n * (m * v) := by
        linear_combination (-1 : ℤ) * hZW + (m : ℤ) * hv
      exact mul_left_cancel₀ hnz h1
    exact ⟨u ^ 2 + v ^ 2, by
      linear_combination (-2 : ℤ) * hbh + (-2 : ℤ) * hdf + (b + d + m * u) * hu +
        (b - d + n * v) * hv + (h + f + n * u) * hY + (f - h + m * v) * hW⟩

/-- An odd prime dividing `m² + n²` divides the centre root. -/
theorem IsMSS.prime_dvd_center_of_dvd_norm {a b c d e f g h i T : ℤ}
    (hs : IsMSS a b c d e f g h i T) {m n : ℕ} (hcop : Nat.Coprime m n)
    (hmn : (b + d) * n = (h + f) * m) (hne : h + f ≠ 0) {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2)
    (hd : p ∣ m * m + n * n) : (p : ℤ) ∣ e := by
  have h1 : (p : ℤ) ∣ 8 * e ^ 2 := by
    refine dvd_trans ?_ (hs.norm_dvd_center hcop hmn hne)
    have hc : ((m * m + n * n : ℕ) : ℤ) = (m : ℤ) ^ 2 + (n : ℤ) ^ 2 := by
      push_cast
      ring
    rw [← hc]
    exact_mod_cast hd
  have hpz : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  rcases hpz.dvd_or_dvd h1 with h8 | he
  · exfalso
    have h8' : p ∣ 8 := by exact_mod_cast h8
    have h8'' : p ∣ 2 ^ 3 := by
      norm_num
      exact h8'
    exact hp2 ((Nat.prime_dvd_prime_iff_eq hp Nat.prime_two).mp (hp.dvd_of_dvd_pow h8''))
  · exact hpz.dvd_of_dvd_pow he

/-- **In a magic square of squares with `(b + d) : (h + f) = m : n` in lowest terms, either
`m² + n²` is a square or twice a square, or some prime `p ≡ 1 (mod 8)` divides the centre root
and divides `m² + n²` to an odd power.** -/
theorem IsMSS.norm_isSquare_or_prime_dvd_center {a b c d e f g h i T : ℤ}
    (hs : IsMSS a b c d e f g h i T) {m n : ℕ} (hn : 0 < n) (hcop : Nat.Coprime m n)
    (hmn : (b + d) * n = (h + f) * m) (hne : h + f ≠ 0) (hbd : b ≠ d) :
    IsSquare (m * m + n * n) ∨ IsSquare (2 * (m * m + n * n)) ∨
      ∃ p : ℕ, p.Prime ∧ p % 8 = 1 ∧ (p : ℤ) ∣ e ∧ Odd (padicValNat p (m * m + n * n)) := by
  by_cases hex : ∃ p : ℕ, p.Prime ∧ p ≠ 2 ∧ Odd (padicValNat p (m * m + n * n))
  · obtain ⟨p, hp, hp2, hodd⟩ := hex
    right
    right
    have : Fact p.Prime := ⟨hp⟩
    have hpd : p ∣ m * m + n * n :=
      dvd_of_one_le_padicValNat (by obtain ⟨k, hk⟩ := hodd; omega)
    refine ⟨p, hp, ?_, hs.prime_dvd_center_of_dvd_norm hcop hmn hne hp hp2 hpd, hodd⟩
    have h5 : p % 8 ≠ 5 := fun h5 => (Nat.not_even_iff_odd.mpr hodd)
      (hs.even_padicValNat_ratio_norm hn hcop hmn hne hbd hp h5)
    have h3 : p % 4 ≠ 3 := by
      intro h3
      have hd : (p : ℤ) ∣ (m : ℤ) ^ 2 + (n : ℤ) ^ 2 := by
        have hc : ((m * m + n * n : ℕ) : ℤ) = (m : ℤ) ^ 2 + (n : ℤ) ^ 2 := by
          push_cast
          ring
        rw [← hc]
        exact_mod_cast hpd
      obtain ⟨h1, h2⟩ := dvd_of_dvd_sq_add_sq h3 hd
      have h1' : p ∣ m := by exact_mod_cast h1
      have h2' : p ∣ n := by exact_mod_cast h2
      have hg := Nat.dvd_gcd h1' h2'
      rw [hcop.gcd_eq_one] at hg
      exact hp.one_lt.ne' (Nat.dvd_one.mp hg)
    have hodd2 : p % 2 = 1 := (hp.eq_two_or_odd).resolve_left hp2
    omega
  · push Not at hex
    have hN : m * m + n * n ≠ 0 := by
      have := Nat.mul_pos hn hn
      omega
    rcases isSquare_or_two_mul_of_even_val hN
      (fun p hp hp2 => Nat.not_odd_iff_even.mp (hex p hp hp2)) with h | h
    · exact Or.inl h
    · exact Or.inr (Or.inl h)

end MagicSquareOfSquares
