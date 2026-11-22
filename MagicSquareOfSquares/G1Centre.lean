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

* `IsMSS.centre_two_ratios`: the same for the second ratio `m' : n' = (b + d) : (f − h)`, which is
  the first ratio of the square with `h` replaced by `−h`, and the two together. Since
  `8e² = (m² + n²)·k²·(m'² + n'²)` with `k ≠ 0` (`IsMSS.center_factorization`), the number
  `2(m² + n²)(m'² + n'²)` is a square (`IsMSS.two_mul_norms_isSquare`). Hence exactly one of the
  two numbers is even (`IsMSS.norms_parity`), and every odd prime divides them to exponents of
  the same parity (`IsMSS.norms_padicVal_parity`). So either a prime `≡ 1 (mod 8)` divides the
  centre root and divides both numbers to an odd power, or one of `m² + n²`, `m'² + n'²` is a
  square and the other is twice a square.
* `IsMSS.edge_sums_identity`: the same factorisation in terms of the entries.
* `IsMSS.prime_dvd_edge_line`: an odd prime dividing `m² + n²` to an odd power divides the centre
  root and one of the two middle lines entirely (`b` and `h`, or `d` and `f`). With Rabern's
  theorem that no prime `≡ 5 (mod 8)` divides a mid-edge entry of a primitive square
  (`IsMSS.not_dvd_edge`), this gives the statement on primes `≡ 1 (mod 8)` again.
* `IsMSS.primitive_factorization`: for a primitive square, `(m² + n²)(m'² + n'²) = 2e²` and the
  edge roots are `±(mm' ± nn')`, `±(nm' ± mn')`.
* The statements on two ratios use only `b² + h² = d² + f² = 2e²` (`TwoLinesSq`), so they hold
  for any two of the four lines through the centre: `TwoLinesSq.edge_params_two`,
  `TwoLinesSq.two_mul_norms_isSquare`, `TwoLinesSq.norms_parity`,
  `TwoLinesSq.norms_padicVal_parity`, `TwoLinesSq.prime_dvd_line`, with `IsMSS.twoLines` and
  `IsMSS.twoLines_diagonals`, `IsMSS.twoLines_column_diagonal`, … for the six pairs. The statement
  on primes `≡ 1 (mod 8)` is for the pair of middle lines only.

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

/-- Two lines through the centre of a magic square of squares, in terms of the roots: `(b, h)` and
`(d, f)` with `b² + h² = d² + f² = 2e²`. Everything in this file about two ratios uses only this,
so it applies to any two of the four lines through the centre (`IsMSS.twoLines` and the five
theorems after it). -/
def TwoLinesSq (b d e f h : ℤ) : Prop :=
  b ^ 2 + h ^ 2 = 2 * e ^ 2 ∧ d ^ 2 + f ^ 2 = 2 * e ^ 2

/-- The middle column `(b, h)` and the middle row `(d, f)`. -/
theorem IsMSS.twoLines {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T) :
    TwoLinesSq b d e f h := by
  obtain ⟨-, -, hbh, hdf⟩ := IsMagic.opposite
    (show IsMagic (a ^ 2) (b ^ 2) (c ^ 2) (d ^ 2) (e ^ 2) (f ^ 2) (g ^ 2) (h ^ 2) (i ^ 2) T from hs)
  exact ⟨hbh, hdf⟩

/-- The two diagonals `(a, i)` and `(c, g)`. -/
theorem IsMSS.twoLines_diagonals {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T) :
    TwoLinesSq a c e g i := by
  obtain ⟨hai, hcg, -, -⟩ := IsMagic.opposite
    (show IsMagic (a ^ 2) (b ^ 2) (c ^ 2) (d ^ 2) (e ^ 2) (f ^ 2) (g ^ 2) (h ^ 2) (i ^ 2) T from hs)
  exact ⟨hai, hcg⟩

/-- The middle column `(b, h)` and the diagonal `(a, i)`. -/
theorem IsMSS.twoLines_column_diagonal {a b c d e f g h i T : ℤ}
    (hs : IsMSS a b c d e f g h i T) : TwoLinesSq b a e i h := by
  obtain ⟨hai, -, hbh, -⟩ := IsMagic.opposite
    (show IsMagic (a ^ 2) (b ^ 2) (c ^ 2) (d ^ 2) (e ^ 2) (f ^ 2) (g ^ 2) (h ^ 2) (i ^ 2) T from hs)
  exact ⟨hbh, hai⟩

/-- The middle column `(b, h)` and the diagonal `(c, g)`. -/
theorem IsMSS.twoLines_column_antidiagonal {a b c d e f g h i T : ℤ}
    (hs : IsMSS a b c d e f g h i T) : TwoLinesSq b c e g h := by
  obtain ⟨-, hcg, hbh, -⟩ := IsMagic.opposite
    (show IsMagic (a ^ 2) (b ^ 2) (c ^ 2) (d ^ 2) (e ^ 2) (f ^ 2) (g ^ 2) (h ^ 2) (i ^ 2) T from hs)
  exact ⟨hbh, hcg⟩

/-- The middle row `(d, f)` and the diagonal `(a, i)`. -/
theorem IsMSS.twoLines_row_diagonal {a b c d e f g h i T : ℤ}
    (hs : IsMSS a b c d e f g h i T) : TwoLinesSq d a e i f := by
  obtain ⟨hai, -, -, hdf⟩ := IsMagic.opposite
    (show IsMagic (a ^ 2) (b ^ 2) (c ^ 2) (d ^ 2) (e ^ 2) (f ^ 2) (g ^ 2) (h ^ 2) (i ^ 2) T from hs)
  exact ⟨hdf, hai⟩

/-- The middle row `(d, f)` and the diagonal `(c, g)`. -/
theorem IsMSS.twoLines_row_antidiagonal {a b c d e f g h i T : ℤ}
    (hs : IsMSS a b c d e f g h i T) : TwoLinesSq d c e g f := by
  obtain ⟨-, hcg, -, hdf⟩ := IsMagic.opposite
    (show IsMagic (a ^ 2) (b ^ 2) (c ^ 2) (d ^ 2) (e ^ 2) (f ^ 2) (g ^ 2) (h ^ 2) (i ^ 2) T from hs)
  exact ⟨hdf, hcg⟩

/-- The parametrisation behind `norm_dvd_center`, for two lines: with `m > 0`, there are integers `u`, `v` with
`b + d = m·u`, `h + f = n·u`, `b − d = n·v`, `f − h = m·v`, and then
`8e² = (m² + n²)(u² + v²)`. -/
theorem TwoLinesSq.edge_params {b d e f h : ℤ} (hl : TwoLinesSq b d e f h)
    {m n : ℕ} (hm0 : 0 < m) (hcop : Nat.Coprime m n) (hmn : (b + d) * n = (h + f) * m)
    (hne : h + f ≠ 0) :
    ∃ u v : ℤ, b + d = m * u ∧ h + f = n * u ∧ b - d = n * v ∧ f - h = m * v ∧
      8 * e ^ 2 = ((m : ℤ) ^ 2 + (n : ℤ) ^ 2) * (u ^ 2 + v ^ 2) := by
  have hbh := hl.1
  have hdf := hl.2
  have hmz : (m : ℤ) ≠ 0 := by exact_mod_cast hm0.ne'
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
  exact ⟨u, v, hu, hY, hv, hW, by
    linear_combination (-2 : ℤ) * hbh + (-2 : ℤ) * hdf + (b + d + m * u) * hu +
      (b - d + n * v) * hv + (h + f + n * u) * hY + (f - h + m * v) * hW⟩

/-- `TwoLinesSq.edge_params` for the middle column and the middle row. -/
theorem IsMSS.edge_params {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T)
    {m n : ℕ} (hm0 : 0 < m) (hcop : Nat.Coprime m n) (hmn : (b + d) * n = (h + f) * m)
    (hne : h + f ≠ 0) :
    ∃ u v : ℤ, b + d = m * u ∧ h + f = n * u ∧ b - d = n * v ∧ f - h = m * v ∧
      8 * e ^ 2 = ((m : ℤ) ^ 2 + (n : ℤ) ^ 2) * (u ^ 2 + v ^ 2) :=
  hs.twoLines.edge_params hm0 hcop hmn hne

/-- In a magic square of squares with `(b + d)·n = (h + f)·m` and `m`, `n` coprime, the number
`m² + n²` divides `8e²`, where `e` is the centre root. -/
theorem IsMSS.norm_dvd_center {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T)
    {m n : ℕ} (hcop : Nat.Coprime m n) (hmn : (b + d) * n = (h + f) * m) (hne : h + f ≠ 0) :
    ((m : ℤ) ^ 2 + (n : ℤ) ^ 2) ∣ 8 * e ^ 2 := by
  rcases Nat.eq_zero_or_pos m with hm0 | hm0
  · subst hm0
    have hn1 : n = 1 := by simpa using hcop
    subst hn1
    simp
  · obtain ⟨u, v, -, -, -, -, h8⟩ := hs.edge_params hm0 hcop hmn hne
    exact ⟨u ^ 2 + v ^ 2, h8⟩

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

/-- The two ratios of two lines through the centre. If `b² + h² = d² + f² = 2e²`, `m : n` is the
ratio `(b + d) : (h + f)` and `m' : n'` the ratio `(b + d) : (f − h)`, both in lowest terms, then
for some `k ≠ 0`: `b + d = m·k·m'`, `h + f = n·k·m'`, `b − d = n·k·n'`, `f − h = m·k·n'` and
`8e² = (m² + n²)·k²·(m'² + n'²)`. -/
theorem TwoLinesSq.edge_params_two {b d e f h : ℤ} (hl : TwoLinesSq b d e f h)
    {m n m' n' : ℕ} (hcop : Nat.Coprime m n) (hcop' : Nat.Coprime m' n')
    (hmn : (b + d) * n = (h + f) * m) (hmn' : (b + d) * n' = (f - h) * m')
    (hne : h + f ≠ 0) (hne' : f - h ≠ 0) :
    ∃ k : ℤ, k ≠ 0 ∧ b + d = m * (k * m') ∧ h + f = n * (k * m') ∧ b - d = n * (k * n') ∧
      f - h = m * (k * n') ∧
      8 * e ^ 2 = ((m : ℤ) ^ 2 + (n : ℤ) ^ 2) * k ^ 2 * ((m' : ℤ) ^ 2 + (n' : ℤ) ^ 2) := by
  have hbh := hl.1
  have hdf := hl.2
  have hkey : (b - d) * (b + d) = (f - h) * (h + f) := by linear_combination hbh - hdf
  have hbd0 : b + d ≠ 0 := by
    intro h0
    rw [h0, mul_zero] at hkey
    exact mul_ne_zero hne' hne hkey.symm
  have hm : 0 < m := by
    rcases Nat.eq_zero_or_pos m with hm0 | hm0
    · exfalso
      subst hm0
      have hn1 : n = 1 := by simpa using hcop
      subst hn1
      exact hbd0 (by simpa using hmn)
    · exact hm0
  obtain ⟨u, v, hu, hY, hv, hW, h8⟩ := hl.edge_params hm hcop hmn hne
  have hmz : (m : ℤ) ≠ 0 := by exact_mod_cast hm.ne'
  have huv : u * n' = v * m' := by
    have h1 : (m : ℤ) * (u * n') = m * (v * m') := by
      linear_combination hmn' + (-(n' : ℤ)) * hu + (m' : ℤ) * hW
    exact mul_left_cancel₀ hmz h1
  have hc' : IsCoprime (m' : ℤ) (n' : ℤ) := Nat.isCoprime_iff_coprime.mpr hcop'
  obtain ⟨x, y, hxy⟩ := hc'
  have hku : (x * u + y * v) * m' = u := by linear_combination u * hxy - y * huv
  have hkv : (x * u + y * v) * n' = v := by linear_combination v * hxy + x * huv
  refine ⟨x * u + y * v, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro hk0
    rw [hk0, zero_mul] at hku
    exact hne (by rw [hY, ← hku, mul_zero])
  · rw [hku]
    exact hu
  · rw [hku]
    exact hY
  · rw [hkv]
    exact hv
  · rw [hkv]
    exact hW
  · linear_combination h8 - ((m : ℤ) ^ 2 + (n : ℤ) ^ 2) * (u + (x * u + y * v) * m') * hku -
      ((m : ℤ) ^ 2 + (n : ℤ) ^ 2) * (v + (x * u + y * v) * n') * hkv

/-- `8e² = (m² + n²)·k²·(m'² + n'²)` for some `k ≠ 0`. -/
theorem TwoLinesSq.center_factorization {b d e f h : ℤ} (hl : TwoLinesSq b d e f h)
    {m n m' n' : ℕ} (hcop : Nat.Coprime m n) (hcop' : Nat.Coprime m' n')
    (hmn : (b + d) * n = (h + f) * m) (hmn' : (b + d) * n' = (f - h) * m')
    (hne : h + f ≠ 0) (hne' : f - h ≠ 0) :
    ∃ k : ℤ, k ≠ 0 ∧
      8 * e ^ 2 = ((m : ℤ) ^ 2 + (n : ℤ) ^ 2) * k ^ 2 * ((m' : ℤ) ^ 2 + (n' : ℤ) ^ 2) := by
  obtain ⟨k, hk, -, -, -, -, h8⟩ := hl.edge_params_two hcop hcop' hmn hmn' hne hne'
  exact ⟨k, hk, h8⟩

/-- The factorisation of `center_factorization` in terms of the entries: a polynomial identity,
given that opposite edge entries sum to twice the centre. -/
theorem TwoLinesSq.edge_sums_identity {b d e f h : ℤ} (hl : TwoLinesSq b d e f h) :
    ((b + d) ^ 2 + (h + f) ^ 2) * ((b + d) ^ 2 + (f - h) ^ 2) = 8 * e ^ 2 * (b + d) ^ 2 := by
  have hbh := hl.1
  have hdf := hl.2
  linear_combination (2 * e ^ 2 + b ^ 2 + 4 * b * d + 2 * d ^ 2 - 2 * f ^ 2 + h ^ 2) * hbh +
    (-2 * e ^ 2 + 4 * b ^ 2 + 4 * b * d + d ^ 2 + f ^ 2) * hdf

/-- If `8e² = N·k²·N'` with `N`, `N'` positive and `k ≠ 0`, then `N` and `N'` are not both
squares and not both twice a square: otherwise `8e²` or `32e²` would be a nonzero square. -/
theorem not_both_same_type {N N' k e : ℤ} (hN : 0 < N) (hN' : 0 < N') (hk : k ≠ 0)
    (h8 : 8 * e ^ 2 = N * k ^ 2 * N') :
    ¬((∃ r : ℤ, N = r * r) ∧ ∃ r' : ℤ, N' = r' * r') ∧
      ¬((∃ r : ℤ, 2 * N = r * r) ∧ ∃ r' : ℤ, 2 * N' = r' * r') := by
  constructor
  · rintro ⟨⟨r, hr⟩, ⟨r', hr'⟩⟩
    have h := (eq_zero_of_sq_eq_two_mul_sq (x := r * k * r') (y := 2 * e) (by
      rw [hr, hr'] at h8
      linear_combination (-1 : ℤ) * h8)).1
    rcases mul_eq_zero.mp h with h | h
    · rcases mul_eq_zero.mp h with h | h
      · rw [h] at hr
        linarith
      · exact hk h
    · rw [h] at hr'
      linarith
  · rintro ⟨⟨r, hr⟩, ⟨r', hr'⟩⟩
    have h := (eq_zero_of_sq_eq_two_mul_sq (x := r * k * r') (y := 4 * e) (by
      linear_combination (-4 : ℤ) * h8 - 2 * k ^ 2 * N' * hr - k ^ 2 * (r * r) * hr')).1
    rcases mul_eq_zero.mp h with h | h
    · rcases mul_eq_zero.mp h with h | h
      · rw [h] at hr
        linarith
      · exact hk h
    · rw [h] at hr'
      linarith

/-- `2(m² + n²)(m'² + n'²)` is a square: it is `(4e/k)²` in the notation of
`center_factorization`. This is the condition that the centre is a square, written
multiplicatively. -/
theorem TwoLinesSq.two_mul_norms_isSquare {b d e f h : ℤ} (hl : TwoLinesSq b d e f h) {m n m' n' : ℕ} (hcop : Nat.Coprime m n)
    (hcop' : Nat.Coprime m' n') (hmn : (b + d) * n = (h + f) * m)
    (hmn' : (b + d) * n' = (f - h) * m') (hne : h + f ≠ 0) (hne' : f - h ≠ 0) :
    IsSquare (2 * ((m : ℤ) ^ 2 + (n : ℤ) ^ 2) * ((m' : ℤ) ^ 2 + (n' : ℤ) ^ 2)) := by
  obtain ⟨k, hk, h8⟩ := hl.center_factorization hcop hcop' hmn hmn' hne hne'
  have hdvd : k ^ 2 ∣ (4 * e) ^ 2 :=
    ⟨2 * ((m : ℤ) ^ 2 + (n : ℤ) ^ 2) * ((m' : ℤ) ^ 2 + (n' : ℤ) ^ 2), by linear_combination 2 * h8⟩
  have hk4 : k ∣ 4 * e := (Int.pow_dvd_pow_iff (by norm_num)).mp hdvd
  obtain ⟨q, hq⟩ := hk4
  refine ⟨q, ?_⟩
  have hk2 : k ^ 2 ≠ 0 := pow_ne_zero 2 hk
  apply mul_left_cancel₀ hk2
  linear_combination (-2 : ℤ) * h8 + (4 * e + k * q) * hq

/-- Twice a product of two odd numbers is not a square. -/
theorem two_mul_odd_mul_odd_ne_sq {S S' q : ℤ} (hS : Odd S) (hS' : Odd S')
    (h : 2 * S * S' = q * q) : False := by
  have hqq : Even (q * q) := ⟨S * S', by rw [← h]; ring⟩
  have hq : Even q := by
    rcases Int.even_mul.mp hqq with h | h <;> exact h
  obtain ⟨c, rfl⟩ := hq
  have h2 : S * S' = 2 * (c * c) := by
    have h4 : (2 : ℤ) * (S * S') = 2 * (2 * (c * c)) := by linear_combination h
    exact mul_left_cancel₀ (by norm_num) h4
  exact (Int.not_even_iff_odd.mpr (hS.mul hS')) ⟨c * c, by rw [h2]; ring⟩

/-- Eight times a product of two odd numbers is not a square. -/
theorem two_mul_two_odd_ne_sq {s s' q : ℤ} (hs : Odd s) (hs' : Odd s')
    (h : 2 * (2 * s) * (2 * s') = q * q) : False := by
  have hqq : Even (q * q) := ⟨2 * s * (2 * s'), by rw [← h]; ring⟩
  have hq : Even q := by
    rcases Int.even_mul.mp hqq with h | h <;> exact h
  obtain ⟨c, rfl⟩ := hq
  have h2 : 2 * s * s' = c * c := by
    have h4 : (4 : ℤ) * (2 * s * s') = 4 * (c * c) := by linear_combination h
    exact mul_left_cancel₀ (by norm_num) h4
  exact two_mul_odd_mul_odd_ne_sq hs hs' h2

/-- For coprime `m`, `n`, if `m² + n²` is even it is twice an odd number. -/
theorem norm_eq_two_mul_odd {m n : ℕ} (hcop : Nat.Coprime m n)
    (hev : Even ((m : ℤ) ^ 2 + (n : ℤ) ^ 2)) :
    ∃ s : ℤ, Odd s ∧ (m : ℤ) ^ 2 + (n : ℤ) ^ 2 = 2 * s := by
  have hiff := Int.even_add.mp hev
  have hm : Odd (m : ℤ) := by
    by_contra hm
    rw [Int.not_odd_iff_even] at hm
    have hn : Even (n : ℤ) := (Int.even_pow.mp (hiff.mp (Int.even_pow.mpr ⟨hm, by norm_num⟩))).1
    have h2m : (2 : ℤ) ∣ (m : ℤ) := hm.two_dvd
    have h2n : (2 : ℤ) ∣ (n : ℤ) := hn.two_dvd
    have h2m' : 2 ∣ m := by exact_mod_cast h2m
    have h2n' : 2 ∣ n := by exact_mod_cast h2n
    have hg := Nat.dvd_gcd h2m' h2n'
    rw [hcop.gcd_eq_one] at hg
    omega
  have hn : Odd (n : ℤ) := by
    by_contra hn
    rw [Int.not_odd_iff_even] at hn
    have h1 : Even ((m : ℤ) ^ 2) := hiff.mpr (Int.even_pow.mpr ⟨hn, by norm_num⟩)
    exact (Int.not_even_iff_odd.mpr hm) (Int.even_pow.mp h1).1
  obtain ⟨x, hx⟩ := hm
  obtain ⟨y, hy⟩ := hn
  exact ⟨2 * (x ^ 2 + x + y ^ 2 + y) + 1, ⟨x ^ 2 + x + y ^ 2 + y, rfl⟩, by rw [hx, hy]; ring⟩

/-- Exactly one of the numbers `m² + n²`, `m'² + n'²` is even. -/
theorem TwoLinesSq.norms_parity {b d e f h : ℤ} (hl : TwoLinesSq b d e f h)
    {m n m' n' : ℕ} (hcop : Nat.Coprime m n) (hcop' : Nat.Coprime m' n')
    (hmn : (b + d) * n = (h + f) * m) (hmn' : (b + d) * n' = (f - h) * m')
    (hne : h + f ≠ 0) (hne' : f - h ≠ 0) :
    (Even ((m : ℤ) ^ 2 + (n : ℤ) ^ 2) ∧ Odd ((m' : ℤ) ^ 2 + (n' : ℤ) ^ 2)) ∨
      (Odd ((m : ℤ) ^ 2 + (n : ℤ) ^ 2) ∧ Even ((m' : ℤ) ^ 2 + (n' : ℤ) ^ 2)) := by
  obtain ⟨q, hq⟩ := hl.two_mul_norms_isSquare hcop hcop' hmn hmn' hne hne'
  rcases Int.even_or_odd ((m : ℤ) ^ 2 + (n : ℤ) ^ 2) with h1 | h1 <;>
    rcases Int.even_or_odd ((m' : ℤ) ^ 2 + (n' : ℤ) ^ 2) with h2 | h2
  · exfalso
    obtain ⟨s, hs1, e1⟩ := norm_eq_two_mul_odd hcop h1
    obtain ⟨s', hs2, e2⟩ := norm_eq_two_mul_odd hcop' h2
    rw [e1, e2] at hq
    exact two_mul_two_odd_ne_sq hs1 hs2 hq
  · exact Or.inl ⟨h1, h2⟩
  · exact Or.inr ⟨h1, h2⟩
  · exact (two_mul_odd_mul_odd_ne_sq h1 h2 hq).elim

/-- For an odd prime `p`, the exponents of `p` in `m² + n²` and in `m'² + n'²` have the same
parity. So the two numbers have the same odd squarefree part. -/
theorem TwoLinesSq.norms_padicVal_parity {b d e f h : ℤ} (hl : TwoLinesSq b d e f h) {m n m' n' : ℕ} (hn : 0 < n) (hn' : 0 < n')
    (hcop : Nat.Coprime m n) (hcop' : Nat.Coprime m' n') (hmn : (b + d) * n = (h + f) * m)
    (hmn' : (b + d) * n' = (f - h) * m') (hne : h + f ≠ 0) (hne' : f - h ≠ 0)
    {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) :
    Even (padicValNat p (m * m + n * n) + padicValNat p (m' * m' + n' * n')) := by
  obtain ⟨q, hq⟩ := hl.two_mul_norms_isSquare hcop hcop' hmn hmn' hne hne'
  have hN : 2 * (m * m + n * n) * (m' * m' + n' * n') = q.natAbs * q.natAbs := by
    apply Nat.cast_injective (R := ℤ)
    push_cast
    first
      | rw [Int.natAbs_mul_self']
      | rw [abs_mul_abs_self]
    linear_combination hq
  have : Fact p.Prime := ⟨hp⟩
  have hS : m * m + n * n ≠ 0 := by
    have := Nat.mul_pos hn hn
    omega
  have hS' : m' * m' + n' * n' ≠ 0 := by
    have := Nat.mul_pos hn' hn'
    omega
  have hQ : q.natAbs ≠ 0 := by
    intro h0
    rw [h0] at hN
    exact mul_ne_zero (mul_ne_zero two_ne_zero hS) hS' (hN.trans (by simp))
  have hv := congrArg (padicValNat p) hN
  rw [padicValNat.mul (mul_ne_zero two_ne_zero hS) hS', padicValNat.mul two_ne_zero hS,
    padicValNat.mul hQ hQ, padicValNat_primes hp2] at hv
  exact ⟨padicValNat p q.natAbs, by omega⟩

/-- If an odd prime `p` divides `m² + n²` to an odd power, it divides the centre root and one of
the two lines entirely: `b` and `h`, or `d` and `f`. (Indeed `4bd = m²(u² + v²) − (m² + n²)v²` in
the notation of `edge_params`, and `p` divides `u² + v²`.) -/
theorem TwoLinesSq.prime_dvd_line {b d e f h : ℤ} (hl : TwoLinesSq b d e f h) {m n : ℕ}
    (hn : 0 < n) (hcop : Nat.Coprime m n) (hmn : (b + d) * n = (h + f) * m) (hne : h + f ≠ 0)
    {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) (hodd : Odd (padicValNat p (m * m + n * n))) :
    (p : ℤ) ∣ e ∧ (((p : ℤ) ∣ b ∧ (p : ℤ) ∣ h) ∨ ((p : ℤ) ∣ d ∧ (p : ℤ) ∣ f)) := by
  have : Fact p.Prime := ⟨hp⟩
  have hpz : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  have hpS : p ∣ m * m + n * n :=
    dvd_of_one_le_padicValNat (by obtain ⟨k, hk⟩ := hodd; omega)
  have hm0 : 0 < m := by
    rcases Nat.eq_zero_or_pos m with h0 | h0
    · exfalso
      subst h0
      have hn1 : n = 1 := by simpa using hcop
      subst hn1
      simp at hodd
    · exact h0
  obtain ⟨u, v, hu, hY, hv, hW, h8⟩ := hl.edge_params hm0 hcop hmn hne
  have hbh := hl.1
  have hdf := hl.2
  have hS : ((m * m + n * n : ℕ) : ℤ) = (m : ℤ) ^ 2 + (n : ℤ) ^ 2 := by
    push_cast
    ring
  have hpSz : (p : ℤ) ∣ (m : ℤ) ^ 2 + (n : ℤ) ^ 2 := by
    rw [← hS]
    exact_mod_cast hpS
  have not2 : ∀ {j : ℕ}, (p : ℤ) ∣ (2 : ℤ) ^ j → False := by
    intro j hj
    have hj' : p ∣ 2 ^ j := by exact_mod_cast hj
    exact hp2 ((Nat.prime_dvd_prime_iff_eq hp Nat.prime_two).mp (hp.dvd_of_dvd_pow hj'))
  have hpe : (p : ℤ) ∣ e := by
    have h1 : (p : ℤ) ∣ 2 ^ 3 * e ^ 2 := by
      have e8 : (2 : ℤ) ^ 3 * e ^ 2 = ((m : ℤ) ^ 2 + (n : ℤ) ^ 2) * (u ^ 2 + v ^ 2) := by
        linear_combination h8
      rw [e8]
      exact Dvd.dvd.mul_right hpSz _
    rcases hpz.dvd_or_dvd h1 with h2 | he
    · exact (not2 h2).elim
    · exact hpz.dvd_of_dvd_pow he
  have hnz : (0 : ℤ) < n := by exact_mod_cast hn
  have hSpos : (0 : ℤ) < (m : ℤ) ^ 2 + (n : ℤ) ^ 2 := by
    nlinarith [sq_nonneg (m : ℤ), mul_pos hnz hnz]
  have he0 : e ≠ 0 := by
    rintro rfl
    have h0 : u ^ 2 + v ^ 2 = 0 := by
      have h8' : ((m : ℤ) ^ 2 + (n : ℤ) ^ 2) * (u ^ 2 + v ^ 2) = 0 := by linear_combination -h8
      exact (mul_eq_zero.mp h8').resolve_left hSpos.ne'
    have hu2 : u ^ 2 = 0 := by nlinarith [sq_nonneg u, sq_nonneg v]
    have hu0 : u = 0 := by simpa using hu2
    exact hne (by rw [hY, hu0, mul_zero])
  -- `p` divides `u² + v²`, by counting the exponent of `p` in `8e² = (m² + n²)(u² + v²)`
  have hnat : 2 ^ 3 * e.natAbs ^ 2 = (m * m + n * n) * (u.natAbs ^ 2 + v.natAbs ^ 2) := by
    zify
    simp only [sq_abs]
    linear_combination h8
  have hen : e.natAbs ≠ 0 := Int.natAbs_ne_zero.mpr he0
  have hSn : m * m + n * n ≠ 0 := by
    have := Nat.mul_pos hn hn
    omega
  have hNn : u.natAbs ^ 2 + v.natAbs ^ 2 ≠ 0 := by
    intro h0
    rw [h0, mul_zero] at hnat
    exact mul_ne_zero (by norm_num) (pow_ne_zero 2 hen) hnat
  have hp8 : padicValNat p (2 ^ 3) = 0 :=
    padicValNat.eq_zero_of_not_dvd (fun hd => hp2
      ((Nat.prime_dvd_prime_iff_eq hp Nat.prime_two).mp (hp.dvd_of_dvd_pow hd)))
  have hval := congrArg (padicValNat p) hnat
  rw [padicValNat.mul (by norm_num) (pow_ne_zero 2 hen), sq e.natAbs, padicValNat.mul hen hen,
    padicValNat.mul hSn hNn, hp8] at hval
  have hpN' : p ∣ u.natAbs ^ 2 + v.natAbs ^ 2 :=
    dvd_of_one_le_padicValNat (by obtain ⟨k, hk⟩ := hodd; omega)
  have hpN : (p : ℤ) ∣ u ^ 2 + v ^ 2 := by
    have hc : ((u.natAbs ^ 2 + v.natAbs ^ 2 : ℕ) : ℤ) = u ^ 2 + v ^ 2 := by
      push_cast
      simp only [sq_abs]
    rw [← hc]
    exact_mod_cast hpN'
  have hbd : (p : ℤ) ∣ 2 ^ 2 * (b * d) := by
    have e4 : (2 : ℤ) ^ 2 * (b * d) =
        (m : ℤ) ^ 2 * (u ^ 2 + v ^ 2) - ((m : ℤ) ^ 2 + (n : ℤ) ^ 2) * v ^ 2 := by
      linear_combination (b + d + m * u) * hu - (b - d + n * v) * hv
    rw [e4]
    exact dvd_sub (Dvd.dvd.mul_left hpN _) (Dvd.dvd.mul_right hpSz _)
  have hbd' : (p : ℤ) ∣ b * d := by
    rcases hpz.dvd_or_dvd hbd with h2 | h2
    · exact (not2 h2).elim
    · exact h2
  refine ⟨hpe, ?_⟩
  rcases hpz.dvd_or_dvd hbd' with hb | hd
  · left
    refine ⟨hb, hpz.dvd_of_dvd_pow (n := 2) ?_⟩
    have e1 : h ^ 2 = 2 * e ^ 2 - b ^ 2 := by linear_combination hbh
    rw [e1]
    exact dvd_sub (Dvd.dvd.mul_left (dvd_pow hpe two_ne_zero) _) (dvd_pow hb two_ne_zero)
  · right
    refine ⟨hd, hpz.dvd_of_dvd_pow (n := 2) ?_⟩
    have e1 : f ^ 2 = 2 * e ^ 2 - d ^ 2 := by linear_combination hdf
    rw [e1]
    exact dvd_sub (Dvd.dvd.mul_left (dvd_pow hpe two_ne_zero) _) (dvd_pow hd two_ne_zero)

/-- `TwoLinesSq.center_factorization` for the middle column and the middle row. -/
theorem IsMSS.center_factorization {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T)
    {m n m' n' : ℕ} (hcop : Nat.Coprime m n) (hcop' : Nat.Coprime m' n')
    (hmn : (b + d) * n = (h + f) * m) (hmn' : (b + d) * n' = (f - h) * m')
    (hne : h + f ≠ 0) (hne' : f - h ≠ 0) :
    ∃ k : ℤ, k ≠ 0 ∧
      8 * e ^ 2 = ((m : ℤ) ^ 2 + (n : ℤ) ^ 2) * k ^ 2 * ((m' : ℤ) ^ 2 + (n' : ℤ) ^ 2) :=
  hs.twoLines.center_factorization hcop hcop' hmn hmn' hne hne'

/-- The factorisation in terms of the entries of the middle column and the middle row. -/
theorem IsMSS.edge_sums_identity {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T) :
    ((b + d) ^ 2 + (h + f) ^ 2) * ((b + d) ^ 2 + (f - h) ^ 2) = 8 * e ^ 2 * (b + d) ^ 2 :=
  hs.twoLines.edge_sums_identity

/-- `2(m² + n²)(m'² + n'²)` is a square, for the middle column and the middle row. -/
theorem IsMSS.two_mul_norms_isSquare {a b c d e f g h i T : ℤ}
    (hs : IsMSS a b c d e f g h i T) {m n m' n' : ℕ} (hcop : Nat.Coprime m n)
    (hcop' : Nat.Coprime m' n') (hmn : (b + d) * n = (h + f) * m)
    (hmn' : (b + d) * n' = (f - h) * m') (hne : h + f ≠ 0) (hne' : f - h ≠ 0) :
    IsSquare (2 * ((m : ℤ) ^ 2 + (n : ℤ) ^ 2) * ((m' : ℤ) ^ 2 + (n' : ℤ) ^ 2)) :=
  hs.twoLines.two_mul_norms_isSquare hcop hcop' hmn hmn' hne hne'

/-- Exactly one of the numbers `m² + n²`, `m'² + n'²` is even. -/
theorem IsMSS.norms_parity {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T)
    {m n m' n' : ℕ} (hcop : Nat.Coprime m n) (hcop' : Nat.Coprime m' n')
    (hmn : (b + d) * n = (h + f) * m) (hmn' : (b + d) * n' = (f - h) * m')
    (hne : h + f ≠ 0) (hne' : f - h ≠ 0) :
    (Even ((m : ℤ) ^ 2 + (n : ℤ) ^ 2) ∧ Odd ((m' : ℤ) ^ 2 + (n' : ℤ) ^ 2)) ∨
      (Odd ((m : ℤ) ^ 2 + (n : ℤ) ^ 2) ∧ Even ((m' : ℤ) ^ 2 + (n' : ℤ) ^ 2)) :=
  hs.twoLines.norms_parity hcop hcop' hmn hmn' hne hne'

/-- For an odd prime `p`, the exponents of `p` in `m² + n²` and in `m'² + n'²` have the same
parity. -/
theorem IsMSS.norms_padicVal_parity {a b c d e f g h i T : ℤ}
    (hs : IsMSS a b c d e f g h i T) {m n m' n' : ℕ} (hn : 0 < n) (hn' : 0 < n')
    (hcop : Nat.Coprime m n) (hcop' : Nat.Coprime m' n') (hmn : (b + d) * n = (h + f) * m)
    (hmn' : (b + d) * n' = (f - h) * m') (hne : h + f ≠ 0) (hne' : f - h ≠ 0)
    {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) :
    Even (padicValNat p (m * m + n * n) + padicValNat p (m' * m' + n' * n')) :=
  hs.twoLines.norms_padicVal_parity hn hn' hcop hcop' hmn hmn' hne hne' hp hp2

/-- An odd prime dividing `m² + n²` to an odd power divides the centre root and one of the two
middle lines entirely: `b` and `h`, or `d` and `f`. -/
theorem IsMSS.prime_dvd_edge_line {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T)
    {m n : ℕ} (hn : 0 < n) (hcop : Nat.Coprime m n) (hmn : (b + d) * n = (h + f) * m)
    (hne : h + f ≠ 0) {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2)
    (hodd : Odd (padicValNat p (m * m + n * n))) :
    (p : ℤ) ∣ e ∧ (((p : ℤ) ∣ b ∧ (p : ℤ) ∣ h) ∨ ((p : ℤ) ∣ d ∧ (p : ℤ) ∣ f)) :=
  hs.twoLines.prime_dvd_line hn hcop hmn hne hp hp2 hodd

/-- A prime whose square divides `2x²` divides `x`. -/
theorem prime_dvd_of_sq_dvd_two_mul_sq {q : ℕ} (hq : q.Prime) {x : ℤ}
    (h : (q : ℤ) ^ 2 ∣ 2 * x ^ 2) : (q : ℤ) ∣ x := by
  have hqz : Prime (q : ℤ) := Nat.prime_iff_prime_int.mp hq
  by_cases h2 : q = 2
  · subst h2
    obtain ⟨c, hc⟩ := h
    push_cast at hc
    have hx2 : Even (x ^ 2) := ⟨c, by linarith⟩
    have hx : (2 : ℤ) ∣ x := (Int.even_pow.mp hx2).1.two_dvd
    exact_mod_cast hx
  · have h1 : (q : ℤ) ∣ 2 * x ^ 2 := dvd_trans (dvd_pow_self _ two_ne_zero) h
    rcases hqz.dvd_or_dvd h1 with h' | h'
    · exfalso
      have h'' : q ∣ 2 := by exact_mod_cast h'
      exact h2 ((Nat.prime_dvd_prime_iff_eq hq Nat.prime_two).mp h'')
    · exact hqz.dvd_of_dvd_pow h'

/-- **Primitive squares.** If no prime divides all nine roots, then `k = ±2` in
`edge_params_two`: the four edge roots are `±(mm' + nn')`, `±(mm' − nn')`, `±(nm' + mn')`,
`±(nm' − mn')` with one common sign, and `(m² + n²)(m'² + n'²) = 2e²`. -/
theorem IsMSS.primitive_factorization {a b c d e f g h i T : ℤ}
    (hs : IsMSS a b c d e f g h i T) (hprim : Primitive a b c d e f g h i)
    {m n m' n' : ℕ} (hcop : Nat.Coprime m n) (hcop' : Nat.Coprime m' n')
    (hmn : (b + d) * n = (h + f) * m) (hmn' : (b + d) * n' = (f - h) * m')
    (hne : h + f ≠ 0) (hne' : f - h ≠ 0) :
    ∃ s : ℤ, (s = 1 ∨ s = -1) ∧ b = s * (m * m' + n * n') ∧ d = s * (m * m' - n * n') ∧
      f = s * (n * m' + m * n') ∧ h = s * (n * m' - m * n') ∧
      ((m : ℤ) ^ 2 + (n : ℤ) ^ 2) * ((m' : ℤ) ^ 2 + (n' : ℤ) ^ 2) = 2 * e ^ 2 := by
  have hl := hs.twoLines
  obtain ⟨k, hk, h1, h2, h3, h4, h8⟩ := hl.edge_params_two hcop hcop' hmn hmn' hne hne'
  -- `k` is even: otherwise `8e²` would be twice an odd number
  have hkev : Even k := by
    by_contra hko
    rw [Int.not_even_iff_odd] at hko
    rcases hl.norms_parity hcop hcop' hmn hmn' hne hne' with ⟨hE, hO⟩ | ⟨hO, hE⟩
    · obtain ⟨s, hsodd, es⟩ := norm_eq_two_mul_odd hcop hE
      have hodd : Odd (s * k ^ 2 * ((m' : ℤ) ^ 2 + (n' : ℤ) ^ 2)) := (hsodd.mul hko.pow).mul hO
      have h4e : s * k ^ 2 * ((m' : ℤ) ^ 2 + (n' : ℤ) ^ 2) = 4 * e ^ 2 := by
        rw [es] at h8
        have h2' : (2 : ℤ) * (s * k ^ 2 * ((m' : ℤ) ^ 2 + (n' : ℤ) ^ 2)) = 2 * (4 * e ^ 2) := by
          linear_combination (-1 : ℤ) * h8
        exact mul_left_cancel₀ (by norm_num) h2'
      rw [h4e] at hodd
      exact (Int.not_even_iff_odd.mpr hodd) ⟨2 * e ^ 2, by ring⟩
    · obtain ⟨s, hsodd, es⟩ := norm_eq_two_mul_odd hcop' hE
      have hodd : Odd (((m : ℤ) ^ 2 + (n : ℤ) ^ 2) * k ^ 2 * s) := (hO.mul hko.pow).mul hsodd
      have h4e : ((m : ℤ) ^ 2 + (n : ℤ) ^ 2) * k ^ 2 * s = 4 * e ^ 2 := by
        rw [es] at h8
        have h2' : (2 : ℤ) * (((m : ℤ) ^ 2 + (n : ℤ) ^ 2) * k ^ 2 * s) = 2 * (4 * e ^ 2) := by
          linear_combination (-1 : ℤ) * h8
        exact mul_left_cancel₀ (by norm_num) h2'
      rw [h4e] at hodd
      exact (Int.not_even_iff_odd.mpr hodd) ⟨2 * e ^ 2, by ring⟩
  obtain ⟨k1, rfl⟩ := hkev
  have eb : b = k1 * (m * m' + n * n') := by
    have h' : (2 : ℤ) * b = 2 * (k1 * (m * m' + n * n')) := by linear_combination h1 + h3
    exact mul_left_cancel₀ (by norm_num) h'
  have ed : d = k1 * (m * m' - n * n') := by
    have h' : (2 : ℤ) * d = 2 * (k1 * (m * m' - n * n')) := by linear_combination h1 - h3
    exact mul_left_cancel₀ (by norm_num) h'
  have ef : f = k1 * (n * m' + m * n') := by
    have h' : (2 : ℤ) * f = 2 * (k1 * (n * m' + m * n')) := by linear_combination h2 + h4
    exact mul_left_cancel₀ (by norm_num) h'
  have eh : h = k1 * (n * m' - m * n') := by
    have h' : (2 : ℤ) * h = 2 * (k1 * (n * m' - m * n')) := by linear_combination h2 - h4
    exact mul_left_cancel₀ (by norm_num) h'
  have h2e : 2 * e ^ 2 =
      ((m : ℤ) ^ 2 + (n : ℤ) ^ 2) * k1 ^ 2 * ((m' : ℤ) ^ 2 + (n' : ℤ) ^ 2) := by
    have h' : (4 : ℤ) * (2 * e ^ 2) =
        4 * (((m : ℤ) ^ 2 + (n : ℤ) ^ 2) * k1 ^ 2 * ((m' : ℤ) ^ 2 + (n' : ℤ) ^ 2)) := by
      linear_combination h8
    exact mul_left_cancel₀ (by norm_num) h'
  -- a prime dividing `k1` would divide all nine roots
  have hna : k1.natAbs = 1 := by
    by_contra hne1
    obtain ⟨q, hq, hqd⟩ := Nat.exists_prime_and_dvd hne1
    have hqk : (q : ℤ) ∣ k1 := Int.natCast_dvd.mpr hqd
    have qb : (q : ℤ) ∣ b := by rw [eb]; exact Dvd.dvd.mul_right hqk _
    have qd : (q : ℤ) ∣ d := by rw [ed]; exact Dvd.dvd.mul_right hqk _
    have qf : (q : ℤ) ∣ f := by rw [ef]; exact Dvd.dvd.mul_right hqk _
    have qh : (q : ℤ) ∣ h := by rw [eh]; exact Dvd.dvd.mul_right hqk _
    obtain ⟨ca, cc, cg, ci⟩ := IsMagic.corner
      (show IsMagic (a ^ 2) (b ^ 2) (c ^ 2) (d ^ 2) (e ^ 2) (f ^ 2) (g ^ 2) (h ^ 2) (i ^ 2) T from hs)
    have sq : ∀ {x y : ℤ}, (q : ℤ) ∣ x → (q : ℤ) ∣ y → (q : ℤ) ^ 2 ∣ x ^ 2 + y ^ 2 :=
      fun hx hy => dvd_add (pow_dvd_pow_of_dvd hx 2) (pow_dvd_pow_of_dvd hy 2)
    have qe : (q : ℤ) ∣ e := prime_dvd_of_sq_dvd_two_mul_sq hq (by rw [← hl.1]; exact sq qb qh)
    have qa : (q : ℤ) ∣ a := prime_dvd_of_sq_dvd_two_mul_sq hq (by rw [ca]; exact sq qf qh)
    have qc : (q : ℤ) ∣ c := prime_dvd_of_sq_dvd_two_mul_sq hq (by rw [cc]; exact sq qd qh)
    have qg : (q : ℤ) ∣ g := prime_dvd_of_sq_dvd_two_mul_sq hq (by rw [cg]; exact sq qb qf)
    have qi : (q : ℤ) ∣ i := prime_dvd_of_sq_dvd_two_mul_sq hq (by rw [ci]; exact sq qb qd)
    exact hprim q hq ⟨qa, qb, qc, qd, qe, qf, qg, qh, qi⟩
  have hpm : k1 = 1 ∨ k1 = -1 := by
    rcases Int.natAbs_eq k1 with hk1 | hk1
    · left
      rw [hna] at hk1
      exact_mod_cast hk1
    · right
      rw [hna] at hk1
      exact_mod_cast hk1
  refine ⟨k1, hpm, eb, ed, ef, eh, ?_⟩
  rcases hpm with rfl | rfl <;> linear_combination (-1 : ℤ) * h2e

/-- **The two edge ratios together.** `m : n` is the ratio `(b + d) : (h + f)` and `m' : n'` the
ratio `(b + d) : (f − h)`, both in lowest terms. Either a prime `≡ 1 (mod 8)` divides the centre
root `e` and divides both `m² + n²` and `m'² + n'²` to an odd power, or one of the two numbers is
a square and the other is twice a square. -/
theorem IsMSS.centre_two_ratios {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T)
    {m n m' n' : ℕ} (hn : 0 < n) (hn' : 0 < n') (hcop : Nat.Coprime m n)
    (hcop' : Nat.Coprime m' n') (hmn : (b + d) * n = (h + f) * m)
    (hmn' : (b + d) * n' = (f - h) * m') (hne : h + f ≠ 0) (hne' : f - h ≠ 0) (hbd : b ≠ d) :
    (∃ p : ℕ, p.Prime ∧ p % 8 = 1 ∧ (p : ℤ) ∣ e ∧ Odd (padicValNat p (m * m + n * n)) ∧
        Odd (padicValNat p (m' * m' + n' * n'))) ∨
      (IsSquare (m * m + n * n) ∧ IsSquare (2 * (m' * m' + n' * n'))) ∨
      (IsSquare (2 * (m * m + n * n)) ∧ IsSquare (m' * m' + n' * n')) := by
  have hs' : IsMSS a b c d e f g (-h) i T := by
    unfold IsMSS
    rw [neg_sq]
    exact hs
  have e1 : -h + f = f - h := by ring
  have h1 := hs.norm_isSquare_or_prime_dvd_center hn hcop hmn hne hbd
  have h2 := hs'.norm_isSquare_or_prime_dvd_center hn' hcop' (by rw [e1]; exact hmn')
    (by rw [e1]; exact hne') hbd
  obtain ⟨k, hk, h8⟩ := hs.center_factorization hcop hcop' hmn hmn' hne hne'
  have hnz : (0 : ℤ) < n := by exact_mod_cast hn
  have hnz' : (0 : ℤ) < n' := by exact_mod_cast hn'
  have hN : (0 : ℤ) < (m : ℤ) ^ 2 + (n : ℤ) ^ 2 := by
    nlinarith [sq_nonneg (m : ℤ), mul_pos hnz hnz]
  have hN' : (0 : ℤ) < (m' : ℤ) ^ 2 + (n' : ℤ) ^ 2 := by
    nlinarith [sq_nonneg (m' : ℤ), mul_pos hnz' hnz']
  obtain ⟨hA, hB⟩ := not_both_same_type hN hN' hk h8
  have cast1 : ∀ {x y r : ℕ}, x * x + y * y = r * r →
      (x : ℤ) ^ 2 + (y : ℤ) ^ 2 = (r : ℤ) * r := by
    intro x y r hr
    have h := congrArg (Nat.cast : ℕ → ℤ) hr
    push_cast at h
    linear_combination h
  have cast2 : ∀ {x y r : ℕ}, 2 * (x * x + y * y) = r * r →
      2 * ((x : ℤ) ^ 2 + (y : ℤ) ^ 2) = (r : ℤ) * r := by
    intro x y r hr
    have h := congrArg (Nat.cast : ℕ → ℤ) hr
    push_cast at h
    linear_combination h
  -- a prime with odd exponent in one of the two numbers has odd exponent in the other
  have both : ∀ {p : ℕ}, p.Prime → p % 8 = 1 →
      (Odd (padicValNat p (m * m + n * n)) ∨ Odd (padicValNat p (m' * m' + n' * n'))) →
      Odd (padicValNat p (m * m + n * n)) ∧ Odd (padicValNat p (m' * m' + n' * n')) := by
    intro p hp hp8 hodd
    have hp2 : p ≠ 2 := by omega
    obtain ⟨r, hr⟩ := hs.norms_padicVal_parity hn hn' hcop hcop' hmn hmn' hne hne' hp hp2
    rcases hodd with ⟨j, hj⟩ | ⟨j, hj⟩
    · exact ⟨⟨j, hj⟩, ⟨r - j - 1, by omega⟩⟩
    · exact ⟨⟨r - j - 1, by omega⟩, ⟨j, hj⟩⟩
  rcases h1 with ⟨r, hr⟩ | ⟨r, hr⟩ | ⟨p, hp, hp8, hpe, hodd⟩
  · rcases h2 with ⟨r', hr'⟩ | ⟨r', hr'⟩ | ⟨p, hp, hp8, hpe, hodd⟩
    · exact absurd ⟨⟨r, cast1 hr⟩, ⟨r', cast1 hr'⟩⟩ hA
    · exact Or.inr (Or.inl ⟨⟨r, hr⟩, ⟨r', hr'⟩⟩)
    · exact Or.inl ⟨p, hp, hp8, hpe, both hp hp8 (Or.inr hodd)⟩
  · rcases h2 with ⟨r', hr'⟩ | ⟨r', hr'⟩ | ⟨p, hp, hp8, hpe, hodd⟩
    · exact Or.inr (Or.inr ⟨⟨r, hr⟩, ⟨r', hr'⟩⟩)
    · exact absurd ⟨⟨r, cast2 hr⟩, ⟨r', cast2 hr'⟩⟩ hB
    · exact Or.inl ⟨p, hp, hp8, hpe, both hp hp8 (Or.inr hodd)⟩
  · exact Or.inl ⟨p, hp, hp8, hpe, both hp hp8 (Or.inl hodd)⟩

end MagicSquareOfSquares
