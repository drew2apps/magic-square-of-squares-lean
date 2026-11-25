import LeanRank.RankZero

/-!
# A rank-0 certificate checker for `y² = x(x − e₂)(x − e₃)` over `ℚ`

Written by an AI model; not read by a person.  Experiment, not counted.

The argument of `LeanRank.RankZero`, for any curve `E : y² = x(x − e₂)(x − e₃)` with integers
`e₂ ≠ e₃`, both nonzero.  A point with `x(x − e₂)(x − e₃) ≠ 0` has `x = n/e²`, and
`n = d₁a²`, `n − e₂e² = d₂b²`, `d₁d₂(n − e₃e²) = t²`, where `d₁` is a signed product of primes of
`S₁ ⊇ primes(e₂e₃)` and `d₂` one of primes of `S₂ ⊇ primes(e₂(e₂ − e₃))` (`descent_gen`).  A
certificate lists the classes `(d₁, d₂)` that are allowed (`good`, at most four) and tables of
squares modulo prime powers; `Cert.check` (evaluated by `decide +kernel` for each curve) verifies
that every other class is excluded at the real place (`realExcluded`) or modulo one of the prime
powers (`excludedAt`, with the condition `¬ (p ∣ n ∧ p ∣ e)`).  With the images of the three
points of order 2 (witnesses), `#(im μ) ≤ 4 = #E(ℚ)[2]` and `rank E(ℚ) = 0` (`finrank_eq_zero_of`).
-/

open WeierstrassCurve WeierstrassCurve.Affine Polynomial
open LeanRank.RankZero (int_sq_of_rat_sq cls cls_eq cls_one)

namespace LeanRank.RankZeroGen

/-! ### 1. Local tests -/

theorem sq_mem_of_table {M : ℕ} (hM : 0 < M) {L : List ℤ}
    (hL : ∀ k : ℕ, k < M → ((k : ℤ) * k) % (M : ℤ) ∈ L) (z : ℤ) : z ^ 2 % (M : ℤ) ∈ L := by
  have hM' : (0 : ℤ) < M := by exact_mod_cast hM
  have h0 : 0 ≤ z % (M : ℤ) := Int.emod_nonneg z hM'.ne'
  have h1 : z % (M : ℤ) < M := Int.emod_lt_of_pos z hM'
  obtain ⟨k, hk⟩ := Int.eq_ofNat_of_zero_le h0
  have hkM : k < M := by omega
  have := hL k hkM
  rw [← hk] at this
  rw [sq, Int.mul_emod]
  exact this

def pass (e₂ e₃ d₁ d₂ p M s q : ℤ) (L : List ℤ) : Bool :=
  !((d₁ * s) % p == 0 && q % p == 0) &&
    L.any (fun r => (d₁ * s - e₂ * q - d₂ * r) % M == 0) &&
    L.any (fun r => (d₁ * d₂ * (d₁ * s - e₃ * q) - r) % M == 0)

/-- The class `(d₁, d₂)` is excluded modulo `M` (a power of the prime `p`, with `L` the squares
modulo `M`) when no pair of residues `(a², e²)` passes. -/
def excludedAt (e₂ e₃ d₁ d₂ p M : ℤ) (L : List ℤ) : Bool :=
  !(L.any fun s => L.any fun q => pass e₂ e₃ d₁ d₂ p M s q L)

theorem excludedAt_eq_false {e₂ e₃ d₁ d₂ p : ℤ} {M : ℕ} {L : List ℤ} (hM : 0 < M)
    (hp : Prime p) (hpM : p ∣ (M : ℤ)) (hL : ∀ k : ℕ, k < M → ((k : ℤ) * k) % (M : ℤ) ∈ L)
    {a e b t : ℤ} (h1 : d₁ * a ^ 2 - e₂ * e ^ 2 = d₂ * b ^ 2)
    (h2 : d₁ * d₂ * (d₁ * a ^ 2 - e₃ * e ^ 2) = t ^ 2) (hprim : ¬ (p ∣ d₁ * a ^ 2 ∧ p ∣ e)) :
    excludedAt e₂ e₃ d₁ d₂ p M L = false := by
  have hA := Int.mod_modEq (a ^ 2) (M : ℤ)
  have hE := Int.mod_modEq (e ^ 2) (M : ℤ)
  have hB := Int.mod_modEq (b ^ 2) (M : ℤ)
  have hT := Int.mod_modEq (t ^ 2) (M : ℤ)
  have k1 : (d₁ * (a ^ 2 % M) - e₂ * (e ^ 2 % M) - d₂ * (b ^ 2 % M)) % (M : ℤ) = 0 := by
    have := ((hA.mul_left d₁).sub (hE.mul_left e₂)).sub (hB.mul_left d₂)
    have h0 : d₁ * a ^ 2 - e₂ * e ^ 2 - d₂ * b ^ 2 = 0 := by rw [h1]; ring
    rw [h0] at this
    exact this
  have k2 : (d₁ * d₂ * (d₁ * (a ^ 2 % M) - e₃ * (e ^ 2 % M)) - t ^ 2 % M) % (M : ℤ) = 0 := by
    have := (((hA.mul_left d₁).sub (hE.mul_left e₃)).mul_left (d₁ * d₂)).sub hT
    have h0 : d₁ * d₂ * (d₁ * a ^ 2 - e₃ * e ^ 2) - t ^ 2 = 0 := by rw [h2]; ring
    rw [h0] at this
    exact this
  have k3 : ¬ ((d₁ * (a ^ 2 % M)) % p = 0 ∧ (e ^ 2 % M) % p = 0) := by
    rintro ⟨c1, c2⟩
    have m1 : (d₁ * (a ^ 2 % M)) % p = (d₁ * a ^ 2) % p := (hA.mul_left d₁).of_dvd hpM
    have m2 : (e ^ 2 % M) % p = (e ^ 2) % p := hE.of_dvd hpM
    apply hprim
    refine ⟨Int.dvd_of_emod_eq_zero (m1 ▸ c1), ?_⟩
    exact hp.dvd_of_dvd_pow (Int.dvd_of_emod_eq_zero (m2 ▸ c2))
  have hk3 : (!((d₁ * (a ^ 2 % M)) % p == 0 && (e ^ 2 % M) % p == 0)) = true := by
    rw [Bool.not_eq_true', Bool.and_eq_false_iff]
    by_contra hc
    simp only [not_or, Bool.not_eq_false, beq_iff_eq] at hc
    exact k3 hc
  have hpass : pass e₂ e₃ d₁ d₂ p M (a ^ 2 % M) (e ^ 2 % M) L = true := by
    simp only [pass, Bool.and_eq_true, beq_iff_eq, List.any_eq_true]
    exact ⟨⟨hk3, b ^ 2 % M, sq_mem_of_table hM hL b, k1⟩, t ^ 2 % M, sq_mem_of_table hM hL t, k2⟩
  have hany : (L.any fun s => L.any fun q => pass e₂ e₃ d₁ d₂ p M s q L) = true :=
    List.any_eq_true.mpr ⟨_, sq_mem_of_table hM hL a,
      List.any_eq_true.mpr ⟨_, sq_mem_of_table hM hL e, hpass⟩⟩
  simp [excludedAt, hany]

/-! ### 2. The real place -/

def realList (e₂ e₃ d₁ d₂ : ℤ) : List (ℤ × Bool) :=
  [(0, decide (0 < d₁)), (e₂, decide (0 < d₂)), (e₃, decide (0 < d₁ * d₂))]

/-- Excluded at the real place: the three sign conditions on `x` are incompatible. -/
def realExcluded (e₂ e₃ d₁ d₂ : ℤ) : Bool :=
  (realList e₂ e₃ d₁ d₂).any fun u => (realList e₂ e₃ d₁ d₂).any fun v =>
    u.2 && !v.2 && decide (v.1 ≤ u.1)

theorem real_sound {c : List (ℤ × Bool)} {n E : ℤ} (hE : 0 < E)
    (hc : ∀ w ∈ c, (w.2 = true → w.1 * E < n) ∧ (w.2 = false → n < w.1 * E)) :
    (c.any fun u => c.any fun v => u.2 && !v.2 && decide (v.1 ≤ u.1)) = false := by
  by_contra hne
  rw [Bool.not_eq_false, List.any_eq_true] at hne
  obtain ⟨u, hu, huv⟩ := hne
  rw [List.any_eq_true] at huv
  obtain ⟨v, hv, h⟩ := huv
  simp only [Bool.and_eq_true, Bool.not_eq_true', decide_eq_true_eq] at h
  obtain ⟨⟨h1, h2⟩, h3⟩ := h
  have := (hc u hu).1 h1
  have := (hc v hv).2 h2
  nlinarith

theorem sq_pos' {a : ℤ} (ha : a ≠ 0) : 0 < a ^ 2 :=
  lt_of_le_of_ne (sq_nonneg a) (Ne.symm (pow_ne_zero 2 ha))

theorem realExcluded_eq_false {e₂ e₃ d₁ d₂ n e a b t : ℤ} (he : e ≠ 0) (ha : a ≠ 0)
    (hb : b ≠ 0) (ht : t ≠ 0) (hna : n = d₁ * a ^ 2) (hgb : n - e₂ * e ^ 2 = d₂ * b ^ 2)
    (hgt : d₁ * d₂ * (n - e₃ * e ^ 2) = t ^ 2) : realExcluded e₂ e₃ d₁ d₂ = false := by
  have ha2 := sq_pos' ha
  have hb2 := sq_pos' hb
  have ht2 := sq_pos' ht
  have hd12 : d₁ * d₂ ≠ 0 := by
    intro h0; rw [h0, zero_mul] at hgt; linarith
  have hd1 : d₁ ≠ 0 := left_ne_zero_of_mul hd12
  have hd2 : d₂ ≠ 0 := right_ne_zero_of_mul hd12
  apply real_sound (n := n) (E := e ^ 2) (sq_pos' he)
  intro w hw
  simp only [realList, List.mem_cons, List.not_mem_nil, or_false] at hw
  rcases hw with rfl | rfl | rfl
  · simp only [decide_eq_true_eq, decide_eq_false_iff_not, zero_mul]
    constructor
    · intro hd; rw [hna]; exact mul_pos hd ha2
    · intro hd
      have : d₁ < 0 := lt_of_le_of_ne (not_lt.mp hd) hd1
      rw [hna]; exact mul_neg_of_neg_of_pos this ha2
  · simp only [decide_eq_true_eq, decide_eq_false_iff_not]
    constructor
    · intro hd; have := mul_pos hd hb2; linarith
    · intro hd
      have : d₂ < 0 := lt_of_le_of_ne (not_lt.mp hd) hd2
      have := mul_neg_of_neg_of_pos this hb2; linarith
  · simp only [decide_eq_true_eq, decide_eq_false_iff_not]
    constructor
    · intro hd
      by_contra hc
      have := mul_nonpos_of_nonneg_of_nonpos hd.le (by linarith : n - e₃ * e ^ 2 ≤ 0)
      linarith
    · intro hd
      have : d₁ * d₂ < 0 := lt_of_le_of_ne (not_lt.mp hd) hd12
      by_contra hc
      have := mul_nonpos_of_nonpos_of_nonneg this.le (by linarith : 0 ≤ n - e₃ * e ^ 2)
      linarith

/-! ### 3. Square classes supported on a list of primes -/

/-- Products of the sublists of `S` (structural recursion, so that the kernel evaluates it fast;
`List.sublists` goes through `Array` and does not). -/
def prods : List ℕ → List ℕ
  | [] => [1]
  | p :: S => prods S ++ (prods S).map (p * ·)

/-- Signed products of sublists of `S`. -/
def classesOf (S : List ℕ) : List ℤ :=
  (prods S).map (fun v : ℕ => (v : ℤ)) ++ (prods S).map (fun v : ℕ => -(v : ℤ))

/-- Strip the primes of `S` from `N`: `N = u · w² · N'` with `u` a product of a sublist of `S`,
no prime of `S` dividing `N'`, and `N'` coprime to `u · w²`. -/
theorem strip (S : List ℕ) (hS : ∀ p ∈ S, p.Prime) (N : ℕ) (hN : N ≠ 0) :
    ∃ u ∈ prods S, ∃ w N' : ℕ, N = u * w ^ 2 * N' ∧ (∀ p ∈ S, ¬ p ∣ N') ∧
      Nat.Coprime N' (u * w ^ 2) := by
  induction S generalizing N with
  | nil => exact ⟨1, by simp [prods], 1, N, by simp, by simp, by simp⟩
  | cons p S ih =>
    have hp : p.Prime := hS p (by simp)
    obtain ⟨α, N₁, hpN₁, hN₁⟩ := Nat.exists_eq_pow_mul_and_not_dvd hN p hp.one_lt.ne'
    have hN₁0 : N₁ ≠ 0 := by rintro rfl; exact hN (by rw [hN₁, mul_zero])
    obtain ⟨u, hu, w, N', hN', hS', hcop⟩ := ih (fun q hq => hS q (by simp [hq])) N₁ hN₁0
    have hpN' : ¬ p ∣ N' := fun h => hpN₁ (by rw [hN']; exact Dvd.dvd.mul_left h _)
    have hcopp : Nat.Coprime N' (p ^ α) :=
      Nat.Coprime.pow_right _ (Nat.Coprime.symm ((Nat.Prime.coprime_iff_not_dvd hp).mpr hpN'))
    have e2 : p ^ α = p ^ (α % 2) * (p ^ (α / 2)) ^ 2 := by
      rw [← pow_mul, ← pow_add]; congr 1; omega
    have hS'' : ∀ q ∈ p :: S, ¬ q ∣ N' := by
      intro q hq
      rcases List.mem_cons.mp hq with rfl | hq
      · exact hpN'
      · exact hS' q hq
    rcases Nat.mod_two_eq_zero_or_one α with h0 | h1
    · refine ⟨u, List.mem_append_left _ hu, p ^ (α / 2) * w, N', ?_, hS'', ?_⟩
      · rw [hN₁, hN', e2, h0]; ring
      · have : u * (p ^ (α / 2) * w) ^ 2 = p ^ α * (u * w ^ 2) := by
          rw [e2, h0]; ring
        rw [this]; exact Nat.Coprime.mul_right hcopp hcop
    · refine ⟨p * u, List.mem_append_right _ (List.mem_map.mpr ⟨u, hu, rfl⟩),
        p ^ (α / 2) * w, N', ?_, hS'', ?_⟩
      · rw [hN₁, hN', e2, h1]; ring
      · have : p * u * (p ^ (α / 2) * w) ^ 2 = p ^ α * (u * w ^ 2) := by
          rw [e2, h1]; ring
        rw [this]; exact Nat.Coprime.mul_right hcopp hcop

/-- If `g · G` is a square and every prime dividing both is in `S`, then `g = d · a²` with
`d ∈ classesOf S`. -/
theorem split_gen (S : List ℕ) (hS : ∀ p ∈ S, p.Prime) {g G Z : ℤ} (hg : g ≠ 0)
    (hZ : g * G = Z ^ 2)
    (hcop : ∀ p : ℕ, p.Prime → (p : ℤ) ∣ g → (p : ℤ) ∣ G → p ∈ S) :
    ∃ d a : ℤ, d ∈ classesOf S ∧ g = d * a ^ 2 := by
  obtain ⟨u, hu, w, N', hN, hS', hcopN⟩ := strip S hS g.natAbs (Int.natAbs_ne_zero.mpr hg)
  obtain ⟨σ, hσ, hgσ⟩ : ∃ σ : ℤ, (σ = 1 ∨ σ = -1) ∧
      g = σ * (u : ℤ) * (w : ℤ) ^ 2 * (N' : ℤ) := by
    have hc : ((g.natAbs : ℕ) : ℤ) = (u : ℤ) * (w : ℤ) ^ 2 * (N' : ℤ) := by
      rw [hN]; push_cast; ring
    rcases Int.natAbs_eq g with h | h
    · exact ⟨1, Or.inl rfl, by rw [h, hc]; ring⟩
    · exact ⟨-1, Or.inr rfl, by rw [h, hc]; ring⟩
  have hgcd : Int.gcd (N' : ℤ) (σ * (u : ℤ) * (w : ℤ) ^ 2 * G) = 1 := by
    rw [Int.gcd_eq_natAbs, Int.natAbs_natCast]
    apply Nat.coprime_of_dvd
    intro k hk hkN hkG
    rw [Int.natAbs_mul, Int.natAbs_mul, Int.natAbs_mul, Int.natAbs_pow, Int.natAbs_natCast,
      Int.natAbs_natCast] at hkG
    have hσ1 : σ.natAbs = 1 := by rcases hσ with rfl | rfl <;> rfl
    rw [hσ1, one_mul] at hkG
    rcases (Nat.Prime.dvd_mul hk).mp hkG with h | h
    · -- `k` divides `u · w²`, which is coprime to `N'`
      have h1 : k ∣ Nat.gcd N' (u * w ^ 2) := Nat.dvd_gcd hkN (by
        rcases (Nat.Prime.dvd_mul hk).mp h with h' | h'
        · exact Dvd.dvd.mul_right h' _
        · exact Dvd.dvd.mul_left h' _)
      rw [hcopN] at h1
      exact hk.one_lt.ne' (Nat.dvd_one.mp h1)
    · have hkg : (k : ℤ) ∣ g := by
        rw [Int.natCast_dvd, hN]
        exact Dvd.dvd.mul_left hkN _
      have hkG' : (k : ℤ) ∣ G := Int.natCast_dvd.mpr h
      exact hS' k (hcop k hk hkg hkG') hkN
  have hprod : (N' : ℤ) * (σ * (u : ℤ) * (w : ℤ) ^ 2 * G) = Z ^ 2 := by
    rw [← hZ, hgσ]; ring
  obtain ⟨c, hc⟩ := Int.sq_of_gcd_eq_one hgcd hprod
  have hl' : (u : ℤ) ∈ (prods S).map (fun v : ℕ => (v : ℤ)) :=
    List.mem_map_of_mem (f := fun v : ℕ => (v : ℤ)) hu
  have hl'' : -(u : ℤ) ∈ (prods S).map (fun v : ℕ => -(v : ℤ)) :=
    List.mem_map_of_mem (f := fun v : ℕ => -(v : ℤ)) hu
  rcases hc with h | h <;> rcases hσ with rfl | rfl
  · exact ⟨_, (w : ℤ) * c, List.mem_append_left _ hl', by rw [hgσ, h]; ring⟩
  · exact ⟨_, (w : ℤ) * c, List.mem_append_right _ hl'', by rw [hgσ, h]; ring⟩
  · exact ⟨_, (w : ℤ) * c, List.mem_append_right _ hl'', by rw [hgσ, h]; ring⟩
  · exact ⟨_, (w : ℤ) * c, List.mem_append_left _ hl', by rw [hgσ, h]; ring⟩

/-! ### 4. The descent -/

/-- **Descent** on `y² = x(x − e₂)(x − e₃)`: the integer data of a point. -/
theorem descent_gen {e₂ e₃ : ℤ} (S₁ S₂ : List ℕ) (hS₁ : ∀ p ∈ S₁, p.Prime)
    (hS₂ : ∀ p ∈ S₂, p.Prime)
    (hd₁ : ∀ p : ℕ, p.Prime → (p : ℤ) ∣ e₂ * e₃ → p ∈ S₁)
    (hd₂ : ∀ p : ℕ, p.Prime → (p : ℤ) ∣ e₂ * (e₂ - e₃) → p ∈ S₂)
    {x y : ℚ} (h : y ^ 2 = x * (x - e₂) * (x - e₃)) (hx0 : x ≠ 0) (hx2 : x - e₂ ≠ 0)
    (hx3 : x - e₃ ≠ 0) :
    ∃ n e a b t d₁ d₂ : ℤ, d₁ ∈ classesOf S₁ ∧ d₂ ∈ classesOf S₂ ∧ e ≠ 0 ∧ a ≠ 0 ∧ b ≠ 0 ∧
      t ≠ 0 ∧ x * (e : ℚ) ^ 2 = n ∧ n = d₁ * a ^ 2 ∧ n - e₂ * e ^ 2 = d₂ * b ^ 2 ∧
      d₁ * d₂ * (n - e₃ * e ^ 2) = t ^ 2 ∧ IsCoprime n e := by
  obtain ⟨n, q, hq0, hxq, hcop⟩ : ∃ n q : ℤ, 0 < q ∧ x * q = n ∧ IsCoprime n q :=
    ⟨x.num, x.den, by exact_mod_cast x.den_pos, by exact_mod_cast Rat.mul_den_eq_num x,
      Int.isCoprime_iff_gcd_eq_one.mpr x.reduced⟩
  obtain ⟨Y, hY⟩ : ∃ Y : ℤ, Y ^ 2 = q * (n * (n - e₂ * q) * (n - e₃ * q)) := by
    apply int_sq_of_rat_sq (r := y * q ^ 2)
    push_cast
    linear_combination (q : ℚ) ^ 4 * h + (q : ℚ) * ((x * q) ^ 2 + x * q * n + n ^ 2 -
      (e₂ + e₃) * q * (x * q + n) + e₂ * e₃ * q ^ 2) * hxq
  obtain ⟨u, v, huv⟩ := hcop
  have hcop2 : IsCoprime q (n * (n - e₂ * q) * (n - e₃ * q)) := by
    have c1 : IsCoprime q n := ⟨v, u, by linear_combination huv⟩
    have c2 : IsCoprime q (n - e₂ * q) := ⟨v + e₂ * u, u, by linear_combination huv⟩
    have c3 : IsCoprime q (n - e₃ * q) := ⟨v + e₃ * u, u, by linear_combination huv⟩
    exact (c1.mul_right c2).mul_right c3
  obtain ⟨e, he⟩ := Int.sq_of_gcd_eq_one (Int.isCoprime_iff_gcd_eq_one.mp hcop2) hY.symm
  have hqe : q = e ^ 2 := by
    rcases he with he | he
    · exact he
    · nlinarith [sq_nonneg e]
  have he0 : e ≠ 0 := by rintro rfl; rw [hqe] at hq0; norm_num at hq0
  have hne : IsCoprime n e := ⟨u, v * e, by rw [hqe] at huv; linear_combination huv⟩
  have hxe : x * (e : ℚ) ^ 2 = n := by rw [← hxq, hqe]; push_cast; ring
  obtain ⟨Z, hZ⟩ : ∃ Z : ℤ, Z ^ 2 = n * (n - e₂ * e ^ 2) * (n - e₃ * e ^ 2) := by
    apply int_sq_of_rat_sq (r := y * e ^ 3)
    push_cast
    linear_combination (e : ℚ) ^ 6 * h +
      ((x * e ^ 2) ^ 2 + x * e ^ 2 * n + n ^ 2 - (e₂ + e₃) * e ^ 2 * (x * e ^ 2 + n) +
        e₂ * e₃ * (e : ℚ) ^ 4) * hxe
  have heQ : (e : ℚ) ≠ 0 := by exact_mod_cast he0
  have hn0 : n ≠ 0 := by
    intro h0
    have : x * (e : ℚ) ^ 2 = 0 := by rw [hxe, h0]; simp
    rcases mul_eq_zero.mp this with h1 | h1
    · exact hx0 h1
    · exact heQ (pow_eq_zero_iff (n := 2) (by norm_num) |>.mp h1)
  have hg2 : n - e₂ * e ^ 2 ≠ 0 := by
    intro h0
    have : (x - e₂) * (e : ℚ) ^ 2 = 0 := by
      have : ((n - e₂ * e ^ 2 : ℤ) : ℚ) = 0 := by rw [h0]; simp
      push_cast at this
      linear_combination this + hxe
    rcases mul_eq_zero.mp this with h1 | h1
    · exact hx2 h1
    · exact heQ (pow_eq_zero_iff (n := 2) (by norm_num) |>.mp h1)
  have hg3 : n - e₃ * e ^ 2 ≠ 0 := by
    intro h0
    have : (x - e₃) * (e : ℚ) ^ 2 = 0 := by
      have : ((n - e₃ * e ^ 2 : ℤ) : ℚ) = 0 := by rw [h0]; simp
      push_cast at this
      linear_combination this + hxe
    rcases mul_eq_zero.mp this with h1 | h1
    · exact hx3 h1
    · exact heQ (pow_eq_zero_iff (n := 2) (by norm_num) |>.mp h1)
  have hpe : ∀ p : ℕ, p.Prime → (p : ℤ) ∣ n → ¬ (p : ℤ) ∣ e := by
    intro p hp hpn hpe'
    have h1 : (p : ℤ) ∣ 1 := (hne.isUnit_of_dvd' hpn hpe').dvd
    have h1' : p ∣ 1 := by exact_mod_cast h1
    exact hp.one_lt.ne' (Nat.dvd_one.mp h1')
  have hpe2 : ∀ p : ℕ, p.Prime → (p : ℤ) ∣ n → ∀ c : ℤ, (p : ℤ) ∣ c * e ^ 2 → (p : ℤ) ∣ c := by
    intro p hp hpn c hc
    have hpZ : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
    rcases hpZ.dvd_or_dvd hc with h | h
    · exact h
    · exact absurd (hpZ.dvd_of_dvd_pow h) (hpe p hp hpn)
  obtain ⟨d₁, a, hd₁m, hna⟩ := split_gen S₁ hS₁ (G := (n - e₂ * e ^ 2) * (n - e₃ * e ^ 2)) hn0
    (by rw [hZ]; ring) (by
      intro p hp hpn hpG
      have hpZ : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
      apply hd₁ p hp
      rcases hpZ.dvd_or_dvd hpG with h | h
      · have : (p : ℤ) ∣ e₂ * e ^ 2 := by
          have := dvd_sub hpn h; convert this using 1; ring
        exact Dvd.dvd.mul_right (hpe2 p hp hpn e₂ this) _
      · have : (p : ℤ) ∣ e₃ * e ^ 2 := by
          have := dvd_sub hpn h; convert this using 1; ring
        exact Dvd.dvd.mul_left (hpe2 p hp hpn e₃ this) _)
  obtain ⟨d₂, b, hd₂m, hgb⟩ := split_gen S₂ hS₂ (G := n * (n - e₃ * e ^ 2)) hg2
    (by rw [hZ]; ring) (by
      intro p hp hpg hpG
      have hpZ : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
      apply hd₂ p hp
      rcases hpZ.dvd_or_dvd hpG with h | h
      · have : (p : ℤ) ∣ e₂ * e ^ 2 := by
          have := dvd_sub h hpg; convert this using 1; ring
        exact Dvd.dvd.mul_right (hpe2 p hp h e₂ this) _
      · have h4 : (p : ℤ) ∣ (e₃ - e₂) * e ^ 2 := by
          have := dvd_sub hpg h; convert this using 1; ring
        rcases hpZ.dvd_or_dvd h4 with h5 | h5
        · have : (p : ℤ) ∣ e₂ - e₃ := by
            have := dvd_neg.mpr h5; convert this using 1; ring
          exact Dvd.dvd.mul_left this _
        · have hpe' := hpZ.dvd_of_dvd_pow h5
          have hpn : (p : ℤ) ∣ n := by
            have := dvd_add hpg (Dvd.dvd.mul_left (dvd_pow hpe' two_ne_zero) e₂)
            convert this using 1; ring
          exact absurd hpe' (hpe p hp hpn))
  have ha0 : a ≠ 0 := by rintro rfl; apply hn0; rw [hna]; ring
  have hb0 : b ≠ 0 := by rintro rfl; apply hg2; rw [hgb]; ring
  obtain ⟨t, ht⟩ : ∃ t : ℤ, t ^ 2 = d₁ * d₂ * (n - e₃ * e ^ 2) := by
    apply int_sq_of_rat_sq (r := (Z : ℚ) / (a * b))
    have hint : Z ^ 2 = d₁ * d₂ * (n - e₃ * e ^ 2) * (a * b) ^ 2 := by
      linear_combination hZ + (n - e₃ * e ^ 2) * (n - e₂ * e ^ 2) * hna +
        (n - e₃ * e ^ 2) * d₁ * a ^ 2 * hgb
    have hab : ((a : ℚ) * b) ≠ 0 := by
      have : (a : ℚ) ≠ 0 := by exact_mod_cast ha0
      have : (b : ℚ) ≠ 0 := by exact_mod_cast hb0
      positivity
    rw [div_pow, div_eq_iff (pow_ne_zero 2 hab)]
    exact_mod_cast hint
  have ht0 : t ≠ 0 := by
    rintro rfl
    have h0 : d₁ * d₂ * (n - e₃ * e ^ 2) = 0 := by rw [← ht]; ring
    have hd₁0 : d₁ ≠ 0 := by rintro rfl; apply hn0; rw [hna]; ring
    have hd₂0 : d₂ ≠ 0 := by rintro rfl; apply hg2; rw [hgb]; ring
    rcases mul_eq_zero.mp h0 with h1 | h1
    · rcases mul_eq_zero.mp h1 with h2 | h2
      · exact hd₁0 h2
      · exact hd₂0 h2
    · exact hg3 h1
  exact ⟨n, e, a, b, t, d₁, d₂, hd₁m, hd₂m, he0, ha0, hb0, ht0, hxe, hna, hgb, ht.symm, hne⟩

/-! ### 5. Certificates -/

/-- A rank-0 certificate for `y² = x(x − e₂)(x − e₃)`. `tables` are triples `(p, M, L)`: a
prime `p`, a power `M` of it and the list `L` of squares modulo `M`. -/
structure Cert where
  e₂ : ℤ
  e₃ : ℤ
  S₁ : List ℕ
  S₂ : List ℕ
  good : List (ℤ × ℤ)
  tables : List (ℤ × ℕ × List ℤ)

/-- The class `(d₁, d₂)` is allowed or excluded. -/
def Cert.classOK (c : Cert) (d₁ d₂ : ℤ) : Bool :=
  decide ((d₁, d₂) ∈ c.good) || realExcluded c.e₂ c.e₃ d₁ d₂ ||
    c.tables.any fun T => excludedAt c.e₂ c.e₃ d₁ d₂ T.1 T.2.1 T.2.2

/-- **The certificate check**, evaluated by the kernel for each curve. -/
def Cert.check (c : Cert) : Bool :=
  (classesOf c.S₁).all fun d₁ => (classesOf c.S₂).all fun d₂ => c.classOK d₁ d₂

/-- The side conditions on the tables (each is checked by `decide` for each curve). -/
def Cert.TablesOK (c : Cert) : Prop :=
  ∀ T ∈ c.tables, Prime T.1 ∧ 0 < T.2.1 ∧ T.1 ∣ (T.2.1 : ℤ) ∧
    ∀ k : ℕ, k < T.2.1 → ((k : ℤ) * k) % (T.2.1 : ℤ) ∈ T.2.2

theorem Cert.good_of_check (c : Cert) (hc : c.check = true) (hT : c.TablesOK)
    {d₁ d₂ n e a b t : ℤ} (hd₁ : d₁ ∈ classesOf c.S₁) (hd₂ : d₂ ∈ classesOf c.S₂)
    (he : e ≠ 0) (ha : a ≠ 0) (hb : b ≠ 0) (ht : t ≠ 0) (hna : n = d₁ * a ^ 2)
    (hgb : n - c.e₂ * e ^ 2 = d₂ * b ^ 2) (hgt : d₁ * d₂ * (n - c.e₃ * e ^ 2) = t ^ 2)
    (hne : IsCoprime n e) : (d₁, d₂) ∈ c.good := by
  have hok := List.all_eq_true.mp (List.all_eq_true.mp hc d₁ hd₁) d₂ hd₂
  simp only [Cert.classOK, Bool.or_eq_true, decide_eq_true_eq] at hok
  rcases hok with (hg | hr) | hx
  · exact hg
  · rw [realExcluded_eq_false he ha hb ht hna hgb hgt] at hr
    exact absurd hr Bool.false_ne_true
  · exfalso
    obtain ⟨T, hTm, hTx⟩ := List.any_eq_true.mp hx
    obtain ⟨hp, hM, hpM, hL⟩ := hT T hTm
    have hprim : ¬ (T.1 ∣ d₁ * a ^ 2 ∧ T.1 ∣ e) := by
      rintro ⟨h1, h2⟩
      rw [← hna] at h1
      exact hp.not_isUnit (hne.isUnit_of_dvd' h1 h2)
    rw [excludedAt_eq_false hM hp hpM hL (by rw [← hna]; exact hgb) (by rw [← hna]; exact hgt)
      hprim] at hTx
    exact Bool.false_ne_true hTx

end LeanRank.RankZeroGen

#print axioms LeanRank.RankZeroGen.excludedAt_eq_false
#print axioms LeanRank.RankZeroGen.realExcluded_eq_false
#print axioms LeanRank.RankZeroGen.split_gen
#print axioms LeanRank.RankZeroGen.descent_gen
#print axioms LeanRank.RankZeroGen.Cert.good_of_check
