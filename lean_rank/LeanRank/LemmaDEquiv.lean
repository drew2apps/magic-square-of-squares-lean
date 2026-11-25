import LeanRank.Theorem2Family

/-!
# Lemma D is exactly the missing input: `LemmaD ↔ TwiceSquareH`

Written by an AI model (Claude, made by Anthropic) on 25 November 2026 (UTC); not read by a person.
Experiment, not counted in the census.

`TwiceSquareH`: for every coprime pair of integers `m, n` with `m² + n² = 2k²` and `m² ≠ n²`, `H(m, n)` is not a
square in `ℚ`. This is the input that `theorem2_family_of_H` needs (at `(m, n)` and `(m, −n)`), and
`lemmaD_iff` shows it is equivalent to `LemmaD` (the rational points of `D2` have `x = ±1`): the
reduction of Theorem 2 to Lemma D loses nothing.

Proof of `TwiceSquareH → LemmaD` (`d2_integral`): write `x = p/q` in lowest terms; then
`2(p² + q²)H(p, q)` is the square of the integer `wq³`. If `p + q` is odd, `p² + q²` and `H(p, q)` are odd, and
`2·odd` is not a square. If `p, q` are odd, `p² + q² = 2S₁`, `H(p, q) = 4H₁` with
`H₁ = 2S₁(i − j)² + pq³` (`p = 2i + 1`, `q = 2j + 1`; from `H = (p² + q²)(p − q)² + 4pq³`), `S₁` and `H₁` are
coprime and positive, and `S₁H₁` is a square, so `p² + q² = 2a²` and `H(p, q) = (2b)²`.
-/

namespace MagicSquareOfSquares.Family

open MagicSquareOfSquares.G1

/-- `H(m, n)` over `ℤ`. -/
def Hz (m n : ℤ) : ℤ := m ^ 4 - 2 * m ^ 3 * n + 2 * m ^ 2 * n ^ 2 + 2 * m * n ^ 3 + n ^ 4

theorem genH_cast (m n : ℤ) : genH (m : ℚ) (n : ℚ) = ((Hz m n : ℤ) : ℚ) := by
  unfold genH Hz; push_cast; ring

/-- The twice-square form of Lemma D: what `theorem2_family_of_H` needs. -/
def TwiceSquareH : Prop :=
  ∀ m n k : ℤ, IsCoprime m n → m ^ 2 + n ^ 2 = 2 * k ^ 2 → m ^ 2 ≠ n ^ 2 → ¬ IsSquare (genH (m : ℚ) (n : ℚ))

theorem twiceSquareH_of_lemmaD (hD : LemmaD) : TwiceSquareH := by
  intro m n k hmn hk hne
  obtain ⟨-, hno⟩ := odd_of_twice_square hmn hk
  have hn : (n : ℚ) ≠ 0 := by
    intro h
    have : n = 0 := by exact_mod_cast h
    rw [this] at hno
    exact (Int.not_even_iff_odd.mpr hno) ⟨0, rfl⟩
  refine not_isSquare_genH hD (k := (k : ℚ)) hn (by exact_mod_cast hk) ?_ ?_
  · intro h; apply hne; have : m = n := by exact_mod_cast h
    rw [this]
  · intro h; apply hne; have : m = -n := by exact_mod_cast h
    rw [this]; ring

/-- **The integral analysis of `D2`.** If `p, q` are coprime and `2(p² + q²)H(p, q) = z²`, then
`p² + q² = 2a²` and `H(p, q) = (2b)²` for some integers `a, b`. -/
theorem d2_integral {p q z : ℤ} (hpq : IsCoprime p q) (hz : 2 * (p ^ 2 + q ^ 2) * Hz p q = z * z) :
    ∃ a b : ℤ, p ^ 2 + q ^ 2 = 2 * a ^ 2 ∧ Hz p q = (2 * b) ^ 2 := by
  have hid : Hz p q = (p ^ 2 + q ^ 2) * (p - q) ^ 2 + 4 * p * q ^ 3 := by unfold Hz; ring
  rcases Int.even_or_odd (p + q) with hev | hodd
  · -- `p + q` even: `p, q` both odd
    have hpo : Odd p := by
      rcases Int.even_or_odd p with ⟨i, hi⟩ | hp
      · exfalso
        obtain ⟨s, hs⟩ := hev
        obtain ⟨u, v, huv⟩ := hpq
        have : (2 : ℤ) * (u * i + v * (s - i)) = 1 := by
          have hq' : q = 2 * (s - i) := by linear_combination hs - hi
          rw [hi, hq'] at huv; linear_combination huv
        generalize u * i + v * (s - i) = X at this
        omega
      · exact hp
    have hqo : Odd q := by
      have := Even.sub_odd hev hpo
      rwa [add_sub_cancel_left] at this
    obtain ⟨i, rfl⟩ := hpo
    obtain ⟨j, rfl⟩ := hqo
    set S₁ := 2 * i ^ 2 + 2 * i + 2 * j ^ 2 + 2 * j + 1 with hS₁
    set H₁ := 2 * S₁ * (i - j) ^ 2 + (2 * i + 1) * (2 * j + 1) ^ 3 with hH₁
    have eS : (2 * i + 1) ^ 2 + (2 * j + 1) ^ 2 = 2 * S₁ := by rw [hS₁]; ring
    have eH : Hz (2 * i + 1) (2 * j + 1) = 4 * H₁ := by rw [hid, hH₁, hS₁]; ring
    have h16 : (4 : ℤ) ^ 2 ∣ z ^ 2 := ⟨S₁ * H₁, by rw [sq z, ← hz, eS, eH]; ring⟩
    obtain ⟨z', rfl⟩ := (Int.pow_dvd_pow_iff two_ne_zero).mp h16
    have hz' : S₁ * H₁ = z' ^ 2 := by
      have : (16 : ℤ) * (S₁ * H₁) = 16 * z' ^ 2 := by
        rw [eS, eH] at hz; linear_combination hz
      exact mul_left_cancel₀ (by norm_num) this
    -- coprimality of `S₁` and `H₁`
    have cS : IsCoprime ((2 * i + 1) ^ 2 + (2 * j + 1) ^ 2) (2 * i + 1) := by
      have := (hpq.symm.pow_left (m := 2)).add_mul_left_left (2 * i + 1)
      convert this using 1; ring
    have cS' : IsCoprime ((2 * i + 1) ^ 2 + (2 * j + 1) ^ 2) (2 * j + 1) := by
      have := (hpq.pow_left (m := 2)).add_mul_left_left (2 * j + 1)
      convert this using 1; ring
    rw [eS] at cS cS'
    have c1 : IsCoprime S₁ (2 * i + 1) := cS.of_mul_left_right
    have c2 : IsCoprime S₁ (2 * j + 1) := cS'.of_mul_left_right
    have cSH : IsCoprime S₁ H₁ := by
      have := (c1.mul_right (c2.pow_right (n := 3))).add_mul_left_right (2 * (i - j) ^ 2)
      rw [hH₁]; convert this using 1; ring
    have hSpos : 0 < S₁ := by
      have hi : 0 ≤ i * (i + 1) := by rcases le_or_gt 0 i with h | h <;> nlinarith
      have hj : 0 ≤ j * (j + 1) := by rcases le_or_gt 0 j with h | h <;> nlinarith
      rw [hS₁]; nlinarith
    have hHpos : 0 < H₁ := by
      have h1 : 0 < Hz (2 * i + 1) (2 * j + 1) := by
        have e : Hz (2 * i + 1) (2 * j + 1) =
            ((2 * i + 1) ^ 2 - (2 * j + 1) ^ 2 - (2 * i + 1) * (2 * j + 1)) ^ 2 +
              3 * ((2 * i + 1) * (2 * j + 1)) ^ 2 := by unfold Hz; ring
        have h2 : 0 < ((2 * i + 1) * (2 * j + 1)) ^ 2 := by
          have : (2 * i + 1) * (2 * j + 1) ≠ 0 := by
            apply mul_ne_zero <;> omega
          positivity
        rw [e]; nlinarith [sq_nonneg ((2 * i + 1) ^ 2 - (2 * j + 1) ^ 2 - (2 * i + 1) * (2 * j + 1))]
      rw [eH] at h1; linarith
    obtain ⟨a, ha⟩ := Int.sq_of_gcd_eq_one (Int.isCoprime_iff_gcd_eq_one.mp cSH) hz'
    obtain ⟨b, hb⟩ := Int.sq_of_gcd_eq_one (Int.isCoprime_iff_gcd_eq_one.mp cSH.symm)
      (by rw [mul_comm]; exact hz')
    have ha' : S₁ = a ^ 2 := by
      rcases ha with h | h
      · exact h
      · exfalso; have := sq_nonneg a; linarith
    have hb' : H₁ = b ^ 2 := by
      rcases hb with h | h
      · exact h
      · exfalso; have := sq_nonneg b; linarith
    exact ⟨a, b, by rw [eS, ha'], by rw [eH, hb']; ring⟩
  · -- `p + q` odd: `p² + q²` and `H(p, q)` odd, and `2·odd` is not a square
    exfalso
    have hS : Odd (p ^ 2 + q ^ 2) := by
      have := Odd.sub_even (hodd.pow (n := 2)) (⟨p * q, by ring⟩ : Even (2 * p * q))
      convert this using 1; ring
    have hd : Odd (p - q) := by
      have := Odd.sub_even hodd (⟨q, by ring⟩ : Even (2 * q))
      convert this using 1; ring
    have hH : Odd (Hz p q) := by
      rw [hid]
      exact Odd.add_even (hS.mul (hd.pow (n := 2))) ⟨2 * p * q ^ 3, by ring⟩
    have hO := hS.mul hH
    have hze : Even (z * z) := ⟨(p ^ 2 + q ^ 2) * Hz p q, by rw [← hz]; ring⟩
    obtain ⟨w, rfl⟩ : Even z := (Int.even_mul.mp hze).elim id id
    have : (p ^ 2 + q ^ 2) * Hz p q = 2 * (w * w) :=
      mul_left_cancel₀ (two_ne_zero : (2 : ℤ) ≠ 0) (by linear_combination hz)
    exact (Int.not_even_iff_odd.mpr hO) ⟨w * w, by rw [this]; ring⟩

theorem lemmaD_of_twiceSquareH (h : TwiceSquareH) : LemmaD := by
  intro x w hw
  have hq0 : (x.den : ℚ) ≠ 0 := Nat.cast_ne_zero.mpr x.den_pos.ne'
  have hx : x = (x.num : ℚ) / ((x.den : ℤ) : ℚ) := by rw [Int.cast_natCast]; exact (Rat.num_div_den x).symm
  have hpq : IsCoprime x.num (x.den : ℤ) := by
    rw [Int.isCoprime_iff_gcd_eq_one, Int.gcd_eq_natAbs, Int.natAbs_natCast]
    exact x.reduced
  set p := x.num with hp
  set q : ℤ := (x.den : ℤ) with hq
  have hq0' : (q : ℚ) ≠ 0 := by rw [hq, Int.cast_natCast]; exact hq0
  have hN : (w * (q : ℚ) ^ 3) ^ 2 = ((2 * (p ^ 2 + q ^ 2) * Hz p q : ℤ) : ℚ) := by
    rw [mul_pow, hw, hx]
    unfold Hz; push_cast
    field_simp
  have hsq : IsSquare (2 * (p ^ 2 + q ^ 2) * Hz p q) :=
    Rat.isSquare_intCast_iff.mp ⟨w * (q : ℚ) ^ 3, by rw [← hN, sq]⟩
  obtain ⟨z, hz⟩ := hsq
  obtain ⟨a, b, ha, hb⟩ := d2_integral hpq hz
  by_cases hpq2 : p ^ 2 = q ^ 2
  · have h0 : (p - q) * (p + q) = 0 := by linear_combination hpq2
    rcases mul_eq_zero.mp h0 with h1 | h1
    · left
      have : p = q := by linear_combination h1
      rw [hx, this, div_self hq0']
    · right
      have : p = -q := by linear_combination h1
      rw [hx, this, Int.cast_neg, neg_div, div_self hq0']
  · exfalso
    exact h p q a hpq ha hpq2 ⟨2 * (b : ℚ), by rw [genH_cast, hb]; push_cast; ring⟩

/-- **Lemma D is equivalent to its twice-square form.** -/
theorem lemmaD_iff : LemmaD ↔ TwiceSquareH := ⟨twiceSquareH_of_lemmaD, lemmaD_of_twiceSquareH⟩

/-- Theorem 2 for every twice-square ratio from `TwiceSquareH` (equivalently `LemmaD`). -/
theorem rank_family' (h : TwiceSquareH) {m n k : ℤ} (hmn : IsCoprime m n) (hn : 1 ≤ n) (hnm : n < m)
    (hk : m ^ 2 + n ^ 2 = 2 * k ^ 2) :
    2 ≤ Module.finrank ℤ (TwoDescent.Eplus (m : ℚ) (n : ℚ)).toAffine.Point ∧
      2 ≤ Module.finrank ℤ (TwoDescent.Eplus (m : ℚ) (-(n : ℚ))).toAffine.Point :=
  rank_family (lemmaD_of_twiceSquareH h) hmn hn hnm hk

/-- **Theorem 2 at every twice-square ratio, any signs and order** (conditional on `LemmaD`): for coprime
integers `m, n` with `m² + n² = 2k²` and `m² ≠ n²` (no condition `0 < n < m`), `R, Q, T` are independent in
`E(ℚ)/2E(ℚ)` on `Eplus m n` and on `Eplus m (−n)`, and both have rank `≥ 2`. -/
theorem rank_family_all (hD : LemmaD) {m n k : ℤ} (hmn : IsCoprime m n) (hk : m ^ 2 + n ^ 2 = 2 * k ^ 2)
    (hne : m ^ 2 ≠ n ^ 2) :
    TwoDescent.Theorem2At (m : ℚ) (n : ℚ) (k : ℚ) ∧ TwoDescent.Theorem2At (m : ℚ) (-(n : ℚ)) (k : ℚ) ∧
      2 ≤ Module.finrank ℤ (TwoDescent.Eplus (m : ℚ) (n : ℚ)).toAffine.Point ∧
      2 ≤ Module.finrank ℤ (TwoDescent.Eplus (m : ℚ) (-(n : ℚ))).toAffine.Point := by
  obtain ⟨hmo, hno⟩ := odd_of_twice_square hmn hk
  have hm0 : (m : ℚ) ≠ 0 := by
    intro h; have : m = 0 := by exact_mod_cast h
    rw [this] at hmo; exact (Int.not_even_iff_odd.mpr hmo) ⟨0, rfl⟩
  have hn0 : (n : ℚ) ≠ 0 := by
    intro h; have : n = 0 := by exact_mod_cast h
    rw [this] at hno; exact (Int.not_even_iff_odd.mpr hno) ⟨0, rfl⟩
  have h1 : (m : ℚ) ≠ n := by
    intro h; apply hne; have : m = n := by exact_mod_cast h
    rw [this]
  have h2 : (m : ℚ) ≠ -n := by
    intro h; apply hne; have : m = -n := by exact_mod_cast h
    rw [this]; ring
  have hp : (m : ℚ) + n ≠ 0 := by intro h; apply h2; linear_combination h
  have hd : (m : ℚ) - n ≠ 0 := by intro h; apply h1; linear_combination h
  have hp' : (m : ℚ) + -n ≠ 0 := by intro h; apply h1; linear_combination h
  have hd' : (m : ℚ) - -n ≠ 0 := by intro h; apply h2; linear_combination h
  have hkQ : (m : ℚ) ^ 2 + (n : ℚ) ^ 2 = 2 * (k : ℚ) ^ 2 := by exact_mod_cast hk
  have hkQ' : (m : ℚ) ^ 2 + (-(n : ℚ)) ^ 2 = 2 * (k : ℚ) ^ 2 := by rw [neg_sq]; exact hkQ
  have hH : ¬ IsSquare (genH (m : ℚ) (n : ℚ)) := not_isSquare_genH hD hn0 hkQ h1 h2
  have hH' : ¬ IsSquare (genH (m : ℚ) (-(n : ℚ))) :=
    not_isSquare_genH hD (neg_ne_zero.mpr hn0) hkQ' h2 (by rw [neg_neg]; exact h1)
  have hHH : ¬ IsSquare (genH (m : ℚ) n * genH m (-n)) := by
    rw [gen_H_mul]
    exact LemmaB.lemmaB_rat hmn (by exact_mod_cast hm0) (by exact_mod_cast hn0) hne
  have s2 : ¬ IsSquare (2 * (m : ℚ) * n) := not_isSquare_two_mul hmo hno
  have s2' : ¬ IsSquare (2 * (m : ℚ) * -(n : ℚ)) := by
    have := not_isSquare_two_mul hmo hno.neg
    rwa [Int.cast_neg] at this
  have t1 : TwoDescent.Theorem2At (m : ℚ) (n : ℚ) (k : ℚ) :=
    TwoDescent.theorem2_of _ _ _ hkQ hm0 hn0 hp hd hH hH' hHH s2
      (by rw [show (-2 : ℚ) * m * n = 2 * m * -n by ring]; exact s2')
  have t2 : TwoDescent.Theorem2At (m : ℚ) (-(n : ℚ)) (k : ℚ) :=
    TwoDescent.theorem2_of _ _ _ hkQ' hm0 (neg_ne_zero.mpr hn0) hp' hd' hH' (by rwa [neg_neg])
      (by rwa [neg_neg, mul_comm]) s2' (by rw [show (-2 : ℚ) * m * -n = 2 * m * n by ring]; exact s2)
  exact ⟨t1, t2, Rank.two_le_rank_of_theorem2At hm0 hn0 hp hd t1,
    Rank.two_le_rank_of_theorem2At hm0 (neg_ne_zero.mpr hn0) hp' hd' t2⟩

/-- `LemmaD` cannot be sharpened to `x = 1`: `(−1, 4)` is a point of `D2`. -/
theorem control_lemmaD_sharp :
    ¬ ∀ x w : ℚ, w ^ 2 = 2 * (x ^ 2 + 1) * (x ^ 4 - 2 * x ^ 3 + 2 * x ^ 2 + 2 * x + 1) → x = 1 := by
  intro h
  have := h (-1) 4 (by norm_num)
  norm_num at this

end MagicSquareOfSquares.Family

#print axioms MagicSquareOfSquares.Family.rank_family_all
#print axioms MagicSquareOfSquares.Family.control_lemmaD_sharp
#print axioms MagicSquareOfSquares.Family.d2_integral
#print axioms MagicSquareOfSquares.Family.lemmaD_iff
#print axioms MagicSquareOfSquares.Family.rank_family'
