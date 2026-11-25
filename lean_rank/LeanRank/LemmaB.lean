import Mathlib.NumberTheory.PythagoreanTriples
import Mathlib.Data.Rat.Lemmas
import Mathlib.Tactic

/-!
# Lemma B of `UNIFORM_RANK.md` (Euler, `d = 14`) for every ratio

Written by an AI model (Claude, made by Anthropic) on 25 November 2026 (UTC); not read by a person.
Experiment, not counted in the census.

* `pocklington`: `a⁴ − a²b² + b⁴ = c²` has no solution in integers with `a, b` coprime, `ab ≠ 0` and
  `a + b` odd (H. C. Pocklington's equation; the descent below was written for this file).
* `euler14`: `x⁴ + 14x²y² + y⁴ = z²` has no solution in integers with `x, y` coprime, `xy ≠ 0`, `x² ≠ y²`
  (Euler's case `d = 14`; Dickson, *History* II, p. 635).
* `lemmaB`: for coprime integers `m, n` with `mn ≠ 0`, `m² ≠ n²`, `m⁸ + 14m⁴n⁴ + n⁸` is not a square in `ℤ`,
  nor (`lemmaB_rat`) in `ℚ`.

## The descent (`pock_aux`)

Let `a` be even, `b` odd, coprime, `a ≠ 0`, `a⁴ − a²b² + b⁴ = c²`. The primitive triple `(a² − b², ab, c)` gives
`a² − b² = M² − N²`, `ab = 2MN`, `c = M² + N²` with `M, N` coprime. The four-number lemma on `(a/2)·b = M·N`
(`Nat.gcd_mul_gcd_of_coprime_of_mul_eq_mul`) gives `|a/2| = rs`, `|b| = tu`, `|M| = rt`, `|N| = su`; then
`s²(4r² + u²) = t²(r² + u²)`, so `r² + u² = s²κ`, `4r² + u² = t²κ` with `κ ∣ 3`, and `κ = 3` is impossible
modulo 3. So `u² + r² = s²` and `u² + (2r)² = t²`; the second triple gives `u = P² − Q²`, `r = PQ`, hence
`P⁴ − P²Q² + Q⁴ = s²` with `P, Q` coprime of opposite parity and `s < c`. Strong induction on `c`.

`euler14` reduces to `pocklington`: if `x, y` are both odd, `a = (x + y)/2`, `b = (x − y)/2` give
`x⁴ + 14x²y² + y⁴ = 16(a⁴ − a²b² + b⁴)`; if they have opposite parity, the primitive triple
`(x² − y², 4xy, z)` gives `x² − y² = M² − N²`, `2xy = MN`, so `M⁴ − M²N² + N⁴ = (x² + y²)²`.
-/

namespace MagicSquareOfSquares.LemmaB

/-- `a² − b²` and `ab` are coprime when `a`, `b` are. -/
theorem isCoprime_sq_sub_mul {a b : ℤ} (h : IsCoprime a b) : IsCoprime (a ^ 2 - b ^ 2) (a * b) := by
  have h1 : IsCoprime (a ^ 2 - b ^ 2) a := by
    have h2 := ((h.symm.pow_left (m := 2)).neg_left).add_mul_left_left a
    convert h2 using 1; ring
  have h2 : IsCoprime (a ^ 2 - b ^ 2) b := by
    have h3 := (h.pow_left (m := 2)).add_mul_left_left (-b)
    convert h3 using 1; ring
  exact h1.mul_right h2

/-- An odd integer is coprime to `2`. -/
theorem isCoprime_two_of_odd {n : ℤ} (h : Odd n) : IsCoprime n 2 := by
  obtain ⟨k, rfl⟩ := h
  exact ⟨1, -k, by ring⟩

/-- **The four-number lemma** in `ℕ`. -/
theorem four_numbers {A B M N : ℕ} (hAB : A.Coprime B) (hMN : M.Coprime N) (h : A * B = M * N) :
    A = A.gcd M * A.gcd N ∧ B = B.gcd M * B.gcd N ∧ M = A.gcd M * B.gcd M ∧ N = A.gcd N * B.gcd N := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · have := Nat.gcd_mul_gcd_of_coprime_of_mul_eq_mul hAB h.symm
    rw [Nat.gcd_comm M A, Nat.gcd_comm N A] at this
    exact this.symm
  · have := Nat.gcd_mul_gcd_of_coprime_of_mul_eq_mul hAB.symm (h.symm.trans (mul_comm A B))
    rw [Nat.gcd_comm M B, Nat.gcd_comm N B] at this
    exact this.symm
  · exact (Nat.gcd_mul_gcd_of_coprime_of_mul_eq_mul hMN h).symm
  · exact (Nat.gcd_mul_gcd_of_coprime_of_mul_eq_mul hMN.symm (h.trans (mul_comm M N))).symm

/-- If `3 ∣ r² + u²` then `3 ∣ r` and `3 ∣ u`. -/
theorem three_dvd_of_sq_add_sq {r u : ℕ} (h : 3 ∣ r ^ 2 + u ^ 2) : 3 ∣ r ∧ 3 ∣ u := by
  have key : ∀ x y : ZMod 3, x ^ 2 + y ^ 2 = 0 → x = 0 ∧ y = 0 := by decide
  have h0 : ((r ^ 2 + u ^ 2 : ℕ) : ZMod 3) = 0 := (CharP.cast_eq_zero_iff (ZMod 3) 3 _).mpr h
  push_cast at h0
  obtain ⟨hr, hu⟩ := key _ _ h0
  exact ⟨(CharP.cast_eq_zero_iff (ZMod 3) 3 r).mp hr, (CharP.cast_eq_zero_iff (ZMod 3) 3 u).mp hu⟩

/-- The size comparison of the descent: `s < (rt)² + (su)²` for `r, s, t, u ≥ 1`. -/
theorem lt_aux {r s t u : ℤ} (hr : 1 ≤ r) (hs : 1 ≤ s) (ht : 1 ≤ t) (hu : 1 ≤ u) :
    s < (r * t) ^ 2 + (s * u) ^ 2 := by
  have h1 : 1 ≤ r * t := by nlinarith
  have h2 : s ≤ s * u := by nlinarith
  have h3 : s * u ≤ (s * u) ^ 2 := by nlinarith
  nlinarith

/-- **The descent.** No solution of `a⁴ − a²b² + b⁴ = c²` with `a` even, `b` odd, `a, b` coprime, `a ≠ 0`. -/
theorem pock_aux : ∀ c : ℕ, ∀ a b : ℤ, IsCoprime a b → Even a → Odd b → a ≠ 0 →
    a ^ 4 - a ^ 2 * b ^ 2 + b ^ 4 = (c : ℤ) ^ 2 → False := by
  intro c
  induction c using Nat.strong_induction_on with
  | _ c ih =>
  intro a b hab ha hb ha0 h
  have hb0 : b ≠ 0 := by
    rintro rfl
    exact (Int.not_even_iff_odd.mpr hb) ⟨0, by ring⟩
  have hc : (0 : ℤ) < c := by
    rcases Nat.eq_zero_or_pos c with hc | hc
    · exfalso
      subst hc
      rw [Nat.cast_zero] at h
      have h0 : (a ^ 2 - b ^ 2) ^ 2 + (a * b) ^ 2 = 0 := by linear_combination h
      have : 0 < (a * b) ^ 2 := by positivity
      nlinarith [sq_nonneg (a ^ 2 - b ^ 2)]
    · exact_mod_cast hc
  have hpt : PythagoreanTriple (a ^ 2 - b ^ 2) (a * b) c := by
    unfold PythagoreanTriple; linear_combination h
  have hcop : Int.gcd (a ^ 2 - b ^ 2) (a * b) = 1 :=
    Int.isCoprime_iff_gcd_eq_one.mp (isCoprime_sq_sub_mul hab)
  obtain ⟨a', rfl⟩ := ha
  have hodd : (( a' + a') ^ 2 - b ^ 2) % 2 = 1 :=
    Int.odd_iff.mp (Even.sub_odd ⟨2 * a' ^ 2, by ring⟩ hb.pow)
  obtain ⟨M, N, hx, hy, hz, hMN, -, -⟩ :=
    PythagoreanTriple.coprime_classification' hpt hcop hodd hc
  have hy' : a' * b = M * N :=
    mul_left_cancel₀ (two_ne_zero : (2 : ℤ) ≠ 0) (by linear_combination hy)
  have ha'0 : a' ≠ 0 := by rintro rfl; exact ha0 (by ring)
  -- pass to absolute values
  have hAB : a'.natAbs * b.natAbs = M.natAbs * N.natAbs := by
    rw [← Int.natAbs_mul, ← Int.natAbs_mul, hy']
  have copAB : Nat.Coprime a'.natAbs b.natAbs := by
    have h1 : IsCoprime a' b := by
      have : IsCoprime (2 * a') b := by convert hab using 1; ring
      exact this.of_mul_left_right
    have := Int.isCoprime_iff_gcd_eq_one.mp h1
    rwa [Int.gcd_eq_natAbs] at this
  have copMN : Nat.Coprime M.natAbs N.natAbs := by
    have := hMN
    rwa [Int.gcd_eq_natAbs] at this
  obtain ⟨hA, hB, hM, hN⟩ := four_numbers copAB copMN hAB
  set r := a'.natAbs.gcd M.natAbs with hr
  set s := a'.natAbs.gcd N.natAbs with hs
  set t := b.natAbs.gcd M.natAbs with ht
  set u := b.natAbs.gcd N.natAbs with hu
  have hA0 : a'.natAbs ≠ 0 := Int.natAbs_ne_zero.mpr ha'0
  have hB0 : b.natAbs ≠ 0 := Int.natAbs_ne_zero.mpr hb0
  have hr0 : r ≠ 0 := by intro h0; apply hA0; rw [hA, h0, zero_mul]
  have hs0 : s ≠ 0 := by intro h0; apply hA0; rw [hA, h0, mul_zero]
  have ht0 : t ≠ 0 := by intro h0; apply hB0; rw [hB, h0, zero_mul]
  have hu0 : u ≠ 0 := by intro h0; apply hB0; rw [hB, h0, mul_zero]
  have copst : s.Coprime t :=
    Nat.Coprime.coprime_dvd_right ⟨u, hB⟩ (Nat.Coprime.coprime_dvd_left ⟨r, by rw [hA, mul_comm]⟩ copAB)
  have copru : r.Coprime u :=
    Nat.Coprime.coprime_dvd_right ⟨t, by rw [hB, mul_comm]⟩ (Nat.Coprime.coprime_dvd_left ⟨s, hA⟩ copAB)
  -- the equation `4A² − B² = M² − N²` in terms of `r, s, t, u`
  have hE : (4 * ((r : ℤ) * s) ^ 2 - ((t : ℤ) * u) ^ 2) = ((r : ℤ) * t) ^ 2 - ((s : ℤ) * u) ^ 2 := by
    have e1 := Int.natAbs_sq a'
    have e2 := Int.natAbs_sq b
    have e3 := Int.natAbs_sq M
    have e4 := Int.natAbs_sq N
    rw [hA] at e1; rw [hB] at e2; rw [hM] at e3; rw [hN] at e4
    push_cast at e1 e2 e3 e4
    linear_combination 4 * e1 - e2 - e3 + e4 + hx
  have hkey : s ^ 2 * (4 * r ^ 2 + u ^ 2) = t ^ 2 * (r ^ 2 + u ^ 2) := by
    have : ((s ^ 2 * (4 * r ^ 2 + u ^ 2) : ℕ) : ℤ) = ((t ^ 2 * (r ^ 2 + u ^ 2) : ℕ) : ℤ) := by
      push_cast; linear_combination hE
    exact_mod_cast this
  obtain ⟨κ, hκ⟩ : s ^ 2 ∣ r ^ 2 + u ^ 2 :=
    Nat.Coprime.dvd_of_dvd_mul_left (Nat.Coprime.pow 2 2 copst) ⟨4 * r ^ 2 + u ^ 2, hkey.symm⟩
  have hκ' : 4 * r ^ 2 + u ^ 2 = t ^ 2 * κ := by
    apply Nat.eq_of_mul_eq_mul_left (pow_pos (Nat.pos_of_ne_zero hs0) 2)
    rw [hkey, hκ]; ring
  have hκ3 : κ ∣ 3 := by
    have d1 : (κ : ℤ) ∣ ((3 * r ^ 2 : ℕ) : ℤ) := ⟨(t : ℤ) ^ 2 - (s : ℤ) ^ 2, by
      have e1 : ((r ^ 2 + u ^ 2 : ℕ) : ℤ) = ((s ^ 2 * κ : ℕ) : ℤ) := by rw [hκ]
      have e2 : ((4 * r ^ 2 + u ^ 2 : ℕ) : ℤ) = ((t ^ 2 * κ : ℕ) : ℤ) := by rw [hκ']
      push_cast at e1 e2 ⊢
      linear_combination e2 - e1⟩
    have d2 : (κ : ℤ) ∣ ((3 * u ^ 2 : ℕ) : ℤ) := ⟨4 * (s : ℤ) ^ 2 - (t : ℤ) ^ 2, by
      have e1 : ((r ^ 2 + u ^ 2 : ℕ) : ℤ) = ((s ^ 2 * κ : ℕ) : ℤ) := by rw [hκ]
      have e2 : ((4 * r ^ 2 + u ^ 2 : ℕ) : ℤ) = ((t ^ 2 * κ : ℕ) : ℤ) := by rw [hκ']
      push_cast at e1 e2 ⊢
      linear_combination 4 * e1 - e2⟩
    have := Nat.dvd_gcd (Int.natCast_dvd_natCast.mp d1) (Int.natCast_dvd_natCast.mp d2)
    rwa [Nat.gcd_mul_left, Nat.Coprime.gcd_eq_one (Nat.Coprime.pow 2 2 copru), mul_one] at this
  have hκ1 : κ = 1 := by
    rcases (Nat.dvd_prime Nat.prime_three).mp hκ3 with h1 | h3
    · exact h1
    · exfalso
      subst h3
      obtain ⟨h3r, h3u⟩ := three_dvd_of_sq_add_sq (r := r) (u := u) ⟨s ^ 2, by rw [hκ]; ring⟩
      have := Nat.dvd_gcd h3r h3u
      rw [Nat.Coprime.gcd_eq_one copru] at this
      exact absurd (Nat.le_of_dvd one_pos this) (by norm_num)
  subst hκ1
  rw [mul_one] at hκ hκ'
  -- `u` is odd
  have hu_odd : Odd u := by
    have : Odd b.natAbs := Int.natAbs_odd.mpr hb
    rw [hB] at this
    exact (Nat.odd_mul.mp this).2
  -- the triple `(u, 2r, t)`
  have hpt2 : PythagoreanTriple (u : ℤ) (2 * (r : ℤ)) (t : ℤ) := by
    unfold PythagoreanTriple
    have : ((4 * r ^ 2 + u ^ 2 : ℕ) : ℤ) = ((t ^ 2 : ℕ) : ℤ) := by rw [hκ']
    push_cast at this
    linear_combination this
  have hcop2 : Int.gcd (u : ℤ) (2 * (r : ℤ)) = 1 := by
    have : ((2 * r : ℕ) : ℤ) = 2 * (r : ℤ) := by push_cast; ring
    rw [← this, Int.gcd_natCast_natCast]
    exact Nat.Coprime.mul_right (Nat.coprime_two_right.mpr hu_odd) copru.symm
  have hodd2 : (u : ℤ) % 2 = 1 := by
    obtain ⟨j, hj⟩ := hu_odd
    omega
  have ht' : (0 : ℤ) < t := by exact_mod_cast Nat.pos_of_ne_zero ht0
  obtain ⟨P, Q, hu', hr', ht'', hPQ, hpar, -⟩ :=
    PythagoreanTriple.coprime_classification' hpt2 hcop2 hodd2 ht'
  have hrPQ : (r : ℤ) = P * Q :=
    mul_left_cancel₀ (two_ne_zero : (2 : ℤ) ≠ 0) (by linear_combination hr')
  have hnew : P ^ 4 - P ^ 2 * Q ^ 2 + Q ^ 4 = (s : ℤ) ^ 2 := by
    have : ((r ^ 2 + u ^ 2 : ℕ) : ℤ) = ((s ^ 2 : ℕ) : ℤ) := by rw [hκ]
    push_cast at this
    rw [hrPQ, hu'] at this
    linear_combination this
  -- `s < c`
  have hlt : s < c := by
    have e3 := Int.natAbs_sq M
    have e4 := Int.natAbs_sq N
    rw [hM] at e3; rw [hN] at e4
    push_cast at e3 e4
    have hc' : (c : ℤ) = ((r : ℤ) * t) ^ 2 + ((s : ℤ) * u) ^ 2 := by
      rw [hz]; linear_combination -e3 - e4
    have h1 : (1 : ℤ) ≤ r := by exact_mod_cast Nat.pos_of_ne_zero hr0
    have h2 : (1 : ℤ) ≤ t := by exact_mod_cast Nat.pos_of_ne_zero ht0
    have h3 : (1 : ℤ) ≤ u := by exact_mod_cast Nat.pos_of_ne_zero hu0
    have h4 : (1 : ℤ) ≤ s := by exact_mod_cast Nat.pos_of_ne_zero hs0
    have : (s : ℤ) < c := by rw [hc']; exact lt_aux h1 h4 h2 h3
    exact_mod_cast this
  have hPQ' : IsCoprime P Q := Int.isCoprime_iff_gcd_eq_one.mpr hPQ
  have hr0' : (r : ℤ) ≠ 0 := by exact_mod_cast hr0
  rcases hpar with ⟨hP, hQ⟩ | ⟨hP, hQ⟩
  · exact ih s hlt P Q hPQ' (Int.even_iff.mpr hP) (Int.odd_iff.mpr hQ)
      (by rintro rfl; apply hr0'; rw [hrPQ]; ring) hnew
  · exact ih s hlt Q P hPQ'.symm (Int.even_iff.mpr hQ) (Int.odd_iff.mpr hP)
      (by rintro rfl; apply hr0'; rw [hrPQ]; ring) (by linear_combination hnew)

/-- **Pocklington's equation.** `a⁴ − a²b² + b⁴ = c²` has no integer solution with `a, b` coprime, `ab ≠ 0`
and `a + b` odd. -/
theorem pocklington {a b c : ℤ} (hab : IsCoprime a b) (ha : a ≠ 0) (hb : b ≠ 0) (hpar : Odd (a + b))
    (h : a ^ 4 - a ^ 2 * b ^ 2 + b ^ 4 = c ^ 2) : False := by
  have h' : a ^ 4 - a ^ 2 * b ^ 2 + b ^ 4 = ((c.natAbs : ℕ) : ℤ) ^ 2 := by rw [Int.natAbs_sq]; exact h
  rcases Int.even_or_odd a with hae | hao
  · have hbo : Odd b := by
      have := Odd.sub_even hpar hae
      rwa [add_sub_cancel_left] at this
    exact pock_aux _ a b hab hae hbo ha h'
  · have hbe : Even b := by
      have := Odd.sub_odd hpar hao
      rwa [add_sub_cancel_left] at this
    exact pock_aux _ b a hab.symm hbe hao hb (by linear_combination h')

/-- **Euler's case `d = 14`.** `x⁴ + 14x²y² + y⁴ = z²` has no integer solution with `x, y` coprime, `xy ≠ 0`
and `x² ≠ y²`. -/
theorem euler14 {x y z : ℤ} (hxy : IsCoprime x y) (hx : x ≠ 0) (hy : y ≠ 0) (hne : x ^ 2 ≠ y ^ 2)
    (h : x ^ 4 + 14 * x ^ 2 * y ^ 2 + y ^ 4 = z ^ 2) : False := by
  rcases Int.even_or_odd x with hxe | hxo <;> rcases Int.even_or_odd y with hye | hyo
  · -- both even: impossible
    obtain ⟨i, rfl⟩ := hxe
    obtain ⟨j, rfl⟩ := hye
    obtain ⟨u, v, huv⟩ := hxy
    have : (2 : ℤ) * (u * i + v * j) = 1 := by linear_combination huv
    omega
  · -- `x` even, `y` odd
    exact euler14_opp hxy hx hy (Or.inl ⟨hxe, hyo⟩) h
  · exact euler14_opp hxy hx hy (Or.inr ⟨hxo, hye⟩) h
  · -- both odd
    obtain ⟨i, rfl⟩ := hxo
    obtain ⟨j, rfl⟩ := hyo
    have h16 : (4 : ℤ) ^ 2 ∣ z ^ 2 :=
      ⟨(i + j + 1) ^ 4 - (i + j + 1) ^ 2 * (i - j) ^ 2 + (i - j) ^ 4, by rw [← h]; ring⟩
    obtain ⟨w, rfl⟩ := (Int.pow_dvd_pow_iff two_ne_zero).mp h16
    have hab : IsCoprime (i + j + 1) (i - j) := by
      obtain ⟨u, v, huv⟩ := hxy
      exact ⟨u + v, u - v, by linear_combination huv⟩
    refine pocklington (c := w) hab ?_ ?_ ⟨i, by ring⟩
      (mul_left_cancel₀ (by norm_num : (16 : ℤ) ≠ 0) (by linear_combination h))
    · intro h0; apply hne; linear_combination 4 * (i - j) * h0
    · intro h0; apply hne; linear_combination 4 * (i + j + 1) * h0
where
  euler14_opp {x y z : ℤ} (hxy : IsCoprime x y) (hx : x ≠ 0) (hy : y ≠ 0)
      (hp : (Even x ∧ Odd y) ∨ (Odd x ∧ Even y)) (h : x ^ 4 + 14 * x ^ 2 * y ^ 2 + y ^ 4 = z ^ 2) : False := by
    have hodd : Odd (x ^ 2 - y ^ 2) := by
      rcases hp with ⟨h1, h2⟩ | ⟨h1, h2⟩
      · obtain ⟨i, rfl⟩ := h1; exact Even.sub_odd ⟨2 * i ^ 2, by ring⟩ h2.pow
      · obtain ⟨j, rfl⟩ := h2; exact Odd.sub_even h1.pow ⟨2 * j ^ 2, by ring⟩
    have hz0 : z ≠ 0 := by
      rintro rfl
      have : 0 < (x * y) ^ 2 := by positivity
      nlinarith [sq_nonneg (x ^ 2 - y ^ 2), sq_nonneg (x ^ 2 + y ^ 2)]
    have hpt : PythagoreanTriple (x ^ 2 - y ^ 2) (2 * 2 * x * y) ((z.natAbs : ℕ) : ℤ) := by
      unfold PythagoreanTriple; linear_combination h - Int.natAbs_sq z
    have hcop : Int.gcd (x ^ 2 - y ^ 2) (2 * 2 * x * y) = 1 := by
      apply Int.isCoprime_iff_gcd_eq_one.mp
      have h2 := isCoprime_two_of_odd hodd
      have := (h2.mul_right h2).mul_right (isCoprime_sq_sub_mul hxy)
      convert this using 1; ring
    have hzpos : (0 : ℤ) < ((z.natAbs : ℕ) : ℤ) := by
      have := Int.natAbs_pos.mpr hz0
      exact_mod_cast this
    obtain ⟨M, N, hx', hy', -, hMN, hpar, -⟩ :=
      PythagoreanTriple.coprime_classification' hpt hcop (Int.odd_iff.mp hodd) hzpos
    have hMN' : M * N = 2 * x * y :=
      mul_left_cancel₀ (two_ne_zero : (2 : ℤ) ≠ 0) (by linear_combination -hy')
    have hM0 : M ≠ 0 := by rintro rfl; apply mul_ne_zero (mul_ne_zero two_ne_zero hx) hy; linear_combination -hMN'
    have hN0 : N ≠ 0 := by rintro rfl; apply mul_ne_zero (mul_ne_zero two_ne_zero hx) hy; linear_combination -hMN'
    have hparMN : Odd (M + N) := by
      rcases hpar with ⟨h1, h2⟩ | ⟨h1, h2⟩
      · exact Int.odd_iff.mpr (by omega)
      · exact Int.odd_iff.mpr (by omega)
    exact pocklington (c := x ^ 2 + y ^ 2) (Int.isCoprime_iff_gcd_eq_one.mpr hMN) hM0 hN0 hparMN
      (by linear_combination (-(x ^ 2 - y ^ 2) - (M ^ 2 - N ^ 2)) * hx' + (M * N + 2 * x * y) * hMN')

/-- **Lemma B of `UNIFORM_RANK.md`.** For coprime integers `m, n` with `mn ≠ 0` and `m² ≠ n²`,
`m⁸ + 14m⁴n⁴ + n⁸` (`= H(m, n)·H(m, −n)`) is not a square in `ℤ`. -/
theorem lemmaB {m n : ℤ} (hmn : IsCoprime m n) (hm : m ≠ 0) (hn : n ≠ 0) (hne : m ^ 2 ≠ n ^ 2) :
    ¬ IsSquare (m ^ 8 + 14 * m ^ 4 * n ^ 4 + n ^ 8) := by
  rintro ⟨z, hz⟩
  refine euler14 (x := m ^ 2) (y := n ^ 2) (z := z) (hmn.pow_left.pow_right) (pow_ne_zero 2 hm)
    (pow_ne_zero 2 hn) ?_ (by linear_combination hz)
  intro h4
  have h0 : (m ^ 2 - n ^ 2) * (m ^ 2 + n ^ 2) = 0 := by linear_combination h4
  rcases mul_eq_zero.mp h0 with h1 | h1
  · exact hne (by linear_combination h1)
  · have : 0 < m ^ 2 := by positivity
    nlinarith [sq_nonneg n]

/-- **Lemma B over `ℚ`.** -/
theorem lemmaB_rat {m n : ℤ} (hmn : IsCoprime m n) (hm : m ≠ 0) (hn : n ≠ 0) (hne : m ^ 2 ≠ n ^ 2) :
    ¬ IsSquare ((m : ℚ) ^ 8 + 14 * (m : ℚ) ^ 4 * (n : ℚ) ^ 4 + (n : ℚ) ^ 8) := by
  have : ((m : ℚ) ^ 8 + 14 * (m : ℚ) ^ 4 * (n : ℚ) ^ 4 + (n : ℚ) ^ 8) =
      ((m ^ 8 + 14 * m ^ 4 * n ^ 4 + n ^ 8 : ℤ) : ℚ) := by push_cast; ring
  rw [this, Rat.isSquare_intCast_iff]
  exact lemmaB hmn hm hn hne

end MagicSquareOfSquares.LemmaB

#print axioms MagicSquareOfSquares.LemmaB.pocklington
#print axioms MagicSquareOfSquares.LemmaB.euler14
#print axioms MagicSquareOfSquares.LemmaB.lemmaB
#print axioms MagicSquareOfSquares.LemmaB.lemmaB_rat
