import MagicSquareOfSquares.ClassConic

/-!
# Classes of odd valuation when `p` divides `m` once, and the prime 5

* `pass_odd_valuation_gen` (section 4.3 of `COMPLETE.md`, for `m = p·m₁` with any `m₁`): at a
  point `(p·u₁, 1)` with `ω₀² ≡ m₁² + n²u₁² (mod p)`, the classes `p·d₁` with
  `d₁ · 2n(ω₀ + m₁)` a nonzero square modulo `p` pass modulo every power of `p`.
* `prime_mul_passes_at_zero_one`: the point `(0, 1)`.
* `prime_mul_unit_passes_seven_mod_eight`: at a prime `p ≡ 7 (mod 8)` with `p ∥ m`, every
  `d = p·d₁`, `p ∤ d₁`, passes. With `ClassLocal.lean` the row `p ≡ 7 (mod 8)` of the table of
  Theorem C is complete in Lean (lower bounds, modulo every power of `p`).
* `passes_at_five`: at `p = 5` dividing `m`, every `d ≡ ±2 (mod 5)` passes at `(1, 1)`.
* `prime_mul_unit_passes_one_mod_eight`: at a prime `p ≡ 1 (mod 8)` with `p ∥ m`, every
  `d = p·d₁`, `p ∤ d₁`, passes at some point `(p·u₁, 1)`.
-/

namespace MagicSquareOfSquares.CoverDescent

/-- **Classes of odd valuation, `m = p·m₁`.** Let `p` be an odd prime modulo which 2 is a square,
`p ∤ n`, and `ω₀` a unit with `ω₀² ≡ m₁² + (n·u₁)² (mod p)`. Then every `d = p·d₁` such that
`d₁ · 2n(ω₀ + m₁)` is a nonzero square modulo `p` passes at `(p·u₁, 1)` modulo every power of `p`. -/
theorem pass_odd_valuation_gen {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) {m1 n u1 ω0 r t d1 : ℤ}
    (hn : ¬(p : ℤ) ∣ n) (hω0 : ¬(p : ℤ) ∣ ω0) (hω : ω0 ^ 2 ≡ m1 ^ 2 + (n * u1) ^ 2 [ZMOD p])
    (hr : r ^ 2 ≡ 2 [ZMOD p]) (ht : t ^ 2 ≡ d1 * (2 * n * (ω0 + m1)) [ZMOD p])
    (ht0 : ¬(p : ℤ) ∣ t) (k : ℕ) :
    ∀ f : Bool × Bool × Bool,
      GSqM ((p : ℤ) ^ (k + 1)) ((p : ℤ) * d1) (Q ((p : ℤ) * m1) n ((p : ℤ) * u1) 1 f) := by
  have hp' : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  have h2 := not_dvd_two_of_odd hp hp2
  have hp0 : (p : ℤ) ≡ 0 [ZMOD p] := Int.modEq_zero_iff_dvd.mpr dvd_rfl
  have hr0 : ¬(p : ℤ) ∣ r := by
    intro h
    have h1 : (p : ℤ) ∣ r ^ 2 := Dvd.dvd.pow h two_ne_zero
    exact h2 (Int.modEq_zero_iff_dvd.mp (hr.symm.trans (Int.modEq_zero_iff_dvd.mpr h1)))
  have e8 : n ^ 2 ≡ nrm (P8 ((p : ℤ) * m1) n ((p : ℤ) * u1) 1) [ZMOD p] := by
    symm
    calc nrm (P8 ((p : ℤ) * m1) n ((p : ℤ) * u1) 1)
        = (n * 1) ^ 2 + ((p : ℤ) * m1 * ((p : ℤ) * u1)) ^ 2 := rfl
      _ ≡ (n * 1) ^ 2 + ((0 : ℤ) * m1 * (0 * u1)) ^ 2 [ZMOD p] :=
          (((hp0.mul_right m1).mul (hp0.mul_right u1)).pow 2).add_left _
      _ = n ^ 2 := by ring
  have eL : (n * r) ^ 2 ≡ nrm (PL ((p : ℤ) * m1) n ((p : ℤ) * u1) 1) [ZMOD p] := by
    symm
    calc nrm (PL ((p : ℤ) * m1) n ((p : ℤ) * u1) 1)
        = (((p : ℤ) * m1 + n) * ((p : ℤ) * u1 - 1)) ^ 2 +
            (((p : ℤ) * m1 - n) * ((p : ℤ) * u1 + 1)) ^ 2 := rfl
      _ ≡ (((0 : ℤ) * m1 + n) * (0 * u1 - 1)) ^ 2 +
            (((0 : ℤ) * m1 - n) * (0 * u1 + 1)) ^ 2 [ZMOD p] :=
          ((((hp0.mul_right m1).add_right n).mul ((hp0.mul_right u1).sub_right 1)).pow 2).add
            ((((hp0.mul_right m1).sub_right n).mul ((hp0.mul_right u1).add_right 1)).pow 2)
      _ = n ^ 2 * 2 := by ring
      _ ≡ n ^ 2 * r ^ 2 [ZMOD p] := hr.symm.mul_left _
      _ = (n * r) ^ 2 := by ring
  have eU : (n * r) ^ 2 ≡ nrm (PU ((p : ℤ) * m1) n ((p : ℤ) * u1) 1) [ZMOD p] := by
    symm
    calc nrm (PU ((p : ℤ) * m1) n ((p : ℤ) * u1) 1)
        = (((p : ℤ) * m1 - n) * ((p : ℤ) * u1 - 1)) ^ 2 +
            (((p : ℤ) * m1 + n) * ((p : ℤ) * u1 + 1)) ^ 2 := rfl
      _ ≡ (((0 : ℤ) * m1 - n) * (0 * u1 - 1)) ^ 2 +
            (((0 : ℤ) * m1 + n) * (0 * u1 + 1)) ^ 2 [ZMOD p] :=
          ((((hp0.mul_right m1).sub_right n).mul ((hp0.mul_right u1).sub_right 1)).pow 2).add
            ((((hp0.mul_right m1).add_right n).mul ((hp0.mul_right u1).add_right 1)).pow 2)
      _ = n ^ 2 * 2 := by ring
      _ ≡ n ^ 2 * r ^ 2 [ZMOD p] := hr.symm.mul_left _
      _ = (n * r) ^ 2 := by ring
  obtain ⟨ω, hωk, cω⟩ := exists_sq_modEq_pow hp hp2 hω0 hω k
  obtain ⟨w8, h8, c8⟩ := exists_sq_modEq_pow hp hp2 hn e8 k
  obtain ⟨wL, hL, cL⟩ := exists_sq_modEq_pow hp hp2 (not_dvd_mul_of_prime hp' hn hr0) eL k
  obtain ⟨wU, hU, cU⟩ := exists_sq_modEq_pow hp hp2 (not_dvd_mul_of_prime hp' hn hr0) eU k
  have h0 : ((p : ℤ) * ω) ^ 2 ≡ nrm (P0 ((p : ℤ) * m1) n ((p : ℤ) * u1) 1)
      [ZMOD (p : ℤ) ^ (k + 1)] := by
    calc ((p : ℤ) * ω) ^ 2 = (p : ℤ) ^ 2 * ω ^ 2 := by ring
      _ ≡ (p : ℤ) ^ 2 * (m1 ^ 2 + (n * u1) ^ 2) [ZMOD (p : ℤ) ^ (k + 1)] := hωk.mul_left _
      _ = nrm (P0 ((p : ℤ) * m1) n ((p : ℤ) * u1) 1) := by simp only [nrm, P0]; ring
  obtain ⟨D1, hD1⟩ : ∃ D1 : ℤ, D1 = (ω + m1) * (w8 + n * 1) *
      (wL + ((p : ℤ) * m1 + n) * ((p : ℤ) * u1 - 1)) *
      (wU + ((p : ℤ) * m1 - n) * ((p : ℤ) * u1 - 1)) := ⟨_, rfl⟩
  have eD : clsD ((p : ℤ) * m1) n ((p : ℤ) * u1) 1 ((p : ℤ) * ω) w8 wL wU = (p : ℤ) * D1 := by
    unfold clsD
    rw [hD1]
    ring
  have cD1 : D1 ≡ (ω0 + m1) * (n + n * 1) *
      (n * r + ((0 : ℤ) * m1 + n) * (0 * u1 - 1)) *
      (n * r + ((0 : ℤ) * m1 - n) * (0 * u1 - 1)) [ZMOD p] := by
    rw [hD1]
    exact (((cω.add_right m1).mul (c8.add_right _)).mul
      (cL.add (((hp0.mul_right m1).add_right n).mul ((hp0.mul_right u1).sub_right 1)))).mul
      (cU.add (((hp0.mul_right m1).sub_right n).mul ((hp0.mul_right u1).sub_right 1)))
  have hE : (4 * n * t) ^ 2 ≡ d1 * (16 * D1) [ZMOD p] := by
    symm
    calc d1 * (16 * D1)
        ≡ d1 * (16 * ((ω0 + m1) * (n + n * 1) *
            (n * r + ((0 : ℤ) * m1 + n) * (0 * u1 - 1)) *
            (n * r + ((0 : ℤ) * m1 - n) * (0 * u1 - 1)))) [ZMOD p] :=
          (cD1.mul_left 16).mul_left d1
      _ = 16 * n ^ 2 * (d1 * (2 * n * (ω0 + m1))) * (r ^ 2 - 1) := by ring
      _ ≡ 16 * n ^ 2 * t ^ 2 * (2 - 1) [ZMOD p] := (ht.symm.mul_left _).mul (hr.sub_right _)
      _ = (4 * n * t) ^ 2 := by ring
  have h4 : ¬(p : ℤ) ∣ 4 := by
    have e4 : (4 : ℤ) = 2 * 2 := by norm_num
    rw [e4]
    exact not_dvd_mul_of_prime hp' h2 h2
  have hg : ¬(p : ℤ) ∣ 4 * n * t :=
    not_dvd_mul_of_prime hp' (not_dvd_mul_of_prime hp' h4 hn) ht0
  obtain ⟨x, hx, hxg⟩ := exists_sq_modEq_pow hp hp2 hg hE k
  have hxp : ¬(p : ℤ) ∣ x := fun hd =>
    hg (Int.modEq_zero_iff_dvd.mp (hxg.symm.trans (Int.modEq_zero_iff_dvd.mpr hd)))
  obtain ⟨ρ, μ, hρ⟩ : IsCoprime x ((p : ℤ) ^ (k + 1)) :=
    (hp'.irreducible.coprime_iff_not_dvd.mpr hxp).symm.pow_right
  have hρ1 : ρ * x ≡ 1 [ZMOD (p : ℤ) ^ (k + 1)] :=
    Int.modEq_iff_dvd.mpr ⟨μ, by linear_combination -hρ⟩
  have hd : (p : ℤ) * d1 ≡ (ρ * d1) ^ 2 *
      (16 * clsD ((p : ℤ) * m1) n ((p : ℤ) * u1) 1 ((p : ℤ) * ω) w8 wL wU)
      [ZMOD (p : ℤ) ^ (k + 1)] := by
    rw [eD]
    calc (p : ℤ) * d1 = (p : ℤ) * d1 * 1 ^ 2 := by ring
      _ ≡ (p : ℤ) * d1 * (ρ * x) ^ 2 [ZMOD (p : ℤ) ^ (k + 1)] := ((hρ1.pow 2).symm).mul_left _
      _ = (p : ℤ) * d1 * ρ ^ 2 * x ^ 2 := by ring
      _ ≡ (p : ℤ) * d1 * ρ ^ 2 * (d1 * (16 * D1)) [ZMOD (p : ℤ) ^ (k + 1)] := hx.mul_left _
      _ = (ρ * d1) ^ 2 * (16 * ((p : ℤ) * D1)) := by ring
  exact fun f => gsqM_class h0.symm h8.symm hL.symm hU.symm hd f

/-- The point `(0, 1)`: for `m = p·m₁` with `p ∤ m₁`, the classes `p·d₁` with `d₁ · n·m₁` a
nonzero square modulo `p` pass. -/
theorem prime_mul_passes_at_zero_one {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) {m1 n r t d1 : ℤ}
    (hn : ¬(p : ℤ) ∣ n) (hm1 : ¬(p : ℤ) ∣ m1) (hr : r ^ 2 ≡ 2 [ZMOD p])
    (ht : t ^ 2 ≡ d1 * (n * m1) [ZMOD p]) (ht0 : ¬(p : ℤ) ∣ t) (k : ℕ) :
    ∀ f : Bool × Bool × Bool,
      GSqM ((p : ℤ) ^ (k + 1)) ((p : ℤ) * d1) (Q ((p : ℤ) * m1) n ((p : ℤ) * 0) 1 f) := by
  have hp' : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  have h2 := not_dvd_two_of_odd hp hp2
  have hω : m1 ^ 2 ≡ m1 ^ 2 + (n * 0) ^ 2 [ZMOD p] := by
    have e : m1 ^ 2 + (n * 0) ^ 2 = m1 ^ 2 := by ring
    rw [e]
  have ht' : (2 * t) ^ 2 ≡ d1 * (2 * n * (m1 + m1)) [ZMOD p] := by
    calc (2 * t) ^ 2 = 4 * t ^ 2 := by ring
      _ ≡ 4 * (d1 * (n * m1)) [ZMOD p] := ht.mul_left _
      _ = d1 * (2 * n * (m1 + m1)) := by ring
  exact pass_odd_valuation_gen hp hp2 hn hm1 hω hr ht' (not_dvd_mul_of_prime hp' h2 ht0) k

/-- **Every class `p · (unit)` passes at a prime `p ≡ 7 (mod 8)` dividing `m` exactly once**: at
`(0, 1)`, modulo every power of `p`; `d₁·n·m₁` or its negative is a square modulo `p`, and `d`,
`−d` pass together. -/
theorem prime_mul_unit_passes_seven_mod_eight {p : ℕ} (hp : p.Prime) (h7 : p % 8 = 7)
    {m1 n d1 : ℤ} (hn : ¬(p : ℤ) ∣ n) (hm1 : ¬(p : ℤ) ∣ m1) (hd1 : ¬(p : ℤ) ∣ d1) (k : ℕ) :
    ∀ f : Bool × Bool × Bool,
      GSqM ((p : ℤ) ^ (k + 1)) ((p : ℤ) * d1) (Q ((p : ℤ) * m1) n ((p : ℤ) * 0) 1 f) := by
  have : Fact p.Prime := ⟨hp⟩
  have hp2 : p ≠ 2 := by omega
  have hp' : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  obtain ⟨r, hr⟩ := exists_sq_modEq_two hp (Or.inr h7)
  have hc0 : (((d1 * (n * m1) : ℤ)) : ZMod p) ≠ 0 := by
    rw [Ne, ZMod.intCast_zmod_eq_zero_iff_dvd]
    exact not_dvd_mul_of_prime hp' hd1 (not_dvd_mul_of_prime hp' hn hm1)
  by_cases hsq : IsSquare ((((d1 * (n * m1) : ℤ)) : ℤ) : ZMod p)
  · obtain ⟨t, ht, ht0⟩ := exists_int_sq_of_isSquare (p := p) hsq hc0
    exact prime_mul_passes_at_zero_one hp hp2 hn hm1 hr ht ht0 k
  · have hm1' : quadraticChar (ZMod p) (-1) = -1 := by
      rw [quadraticChar_neg_one_iff_not_isSquare, ZMod.exists_sq_eq_neg_one_iff]
      omega
    have ha1 : quadraticChar (ZMod p) ((((d1 * (n * m1) : ℤ)) : ℤ) : ZMod p) = -1 :=
      quadraticChar_neg_one_iff_not_isSquare.mpr hsq
    have hneg : (((-d1 * (n * m1) : ℤ)) : ZMod p) = -1 * (((d1 * (n * m1) : ℤ)) : ZMod p) := by
      push_cast
      ring
    have hq : quadraticChar (ZMod p) (((-d1 * (n * m1) : ℤ)) : ZMod p) = 1 := by
      rw [hneg, map_mul, hm1', ha1]
      norm_num
    have hna0 : (((-d1 * (n * m1) : ℤ)) : ZMod p) ≠ 0 := by
      rw [hneg]
      exact mul_ne_zero (by norm_num) hc0
    obtain ⟨t, ht, ht0⟩ := exists_int_sq_of_isSquare (p := p)
      ((quadraticChar_one_iff_isSquare hna0).mp hq) hna0
    intro f
    have h5 := (prime_mul_passes_at_zero_one hp hp2 (d1 := -d1) hn hm1 hr ht ht0 k f).neg_left
    have e : -((p : ℤ) * -d1) = (p : ℤ) * d1 := by ring
    rwa [e] at h5

/-- **The prime 5**: at `p = 5` dividing `m`, every `d ≡ ±2 (mod 5)` passes at `(1, 1)` modulo
every power of 5. By the law at 5 (`class_mod_five`) these are the only classes of a solution when
25 does not divide `mn(m² − n²)`. -/
theorem passes_at_five {m n d : ℤ} (hm : (5 : ℤ) ∣ m) (hn : ¬(5 : ℤ) ∣ n)
    (hd : d % 5 = 2 ∨ d % 5 = 3) (k : ℕ) :
    ∀ f : Bool × Bool × Bool, GSqM ((5 : ℤ) ^ (k + 1)) d (Q m n 1 1 f) := by
  have hp : Nat.Prime 5 := by norm_num
  have h1 : ¬((5 : ℕ) : ℤ) ∣ 1 := by norm_num
  have h2 : ¬((5 : ℕ) : ℤ) ∣ 2 := by norm_num
  have hs : (2 : ℤ) ^ 2 ≡ 2 * (1 ^ 2 + 1 ^ 2) [ZMOD ((5 : ℕ) : ℤ)] := by norm_num
  have key : ∃ t : ℤ, t ^ 2 ≡ d * (2 * 1 * 1) [ZMOD ((5 : ℕ) : ℤ)] ∧ ¬((5 : ℕ) : ℤ) ∣ t := by
    rcases hd with h | h
    · refine ⟨2, ?_, h2⟩
      unfold Int.ModEq
      push_cast
      omega
    · refine ⟨1, ?_, h1⟩
      unfold Int.ModEq
      push_cast
      omega
  obtain ⟨t, ht, ht0⟩ := key
  obtain ⟨_, _, _, _, _, _, _, _, g⟩ := pass_first_order hp (by norm_num) (u := 1) (v := 1)
    (s := 2) (t := t) (m := m) (n := n) (d := d) (by exact_mod_cast hm) (by exact_mod_cast hn)
    h1 h1 (by norm_num) hs h2 ht ht0 k
  intro f
  have h5 := g f
  exact_mod_cast h5

/-- **Every class `p · (unit)` passes at a prime `p ≡ 1 (mod 8)` dividing `m` exactly once**, at
some point `(p·u₁, 1)`, modulo every power of `p`. If `λ = d₁·n·m₁` is a square modulo `p` the
point is `(0, 1)`; otherwise `ω₀ = m₁(λ² + 1)/(2λ)`, `n·u₁ = m₁(λ² − 1)/(2λ)`, and
`d₁ · 2n(ω₀ + m₁) = (λ + 1)²`. -/
theorem prime_mul_unit_passes_one_mod_eight {p : ℕ} (hp : p.Prime) (h18 : p % 8 = 1)
    {m1 n d1 : ℤ} (hn : ¬(p : ℤ) ∣ n) (hm1 : ¬(p : ℤ) ∣ m1) (hd1 : ¬(p : ℤ) ∣ d1) (k : ℕ) :
    ∃ u1 : ℤ, ∀ f : Bool × Bool × Bool,
      GSqM ((p : ℤ) ^ (k + 1)) ((p : ℤ) * d1) (Q ((p : ℤ) * m1) n ((p : ℤ) * u1) 1 f) := by
  have : Fact p.Prime := ⟨hp⟩
  have hp2 : p ≠ 2 := by omega
  have hp' : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  obtain ⟨r, hr⟩ := exists_sq_modEq_two hp (Or.inl h18)
  have h2ne := zmod_two_ne_zero (p := p) hp2
  have hn0 : ((n : ℤ) : ZMod p) ≠ 0 := by
    rw [Ne, ZMod.intCast_zmod_eq_zero_iff_dvd]
    exact hn
  have hm10 : ((m1 : ℤ) : ZMod p) ≠ 0 := by
    rw [Ne, ZMod.intCast_zmod_eq_zero_iff_dvd]
    exact hm1
  have hd10 : ((d1 : ℤ) : ZMod p) ≠ 0 := by
    rw [Ne, ZMod.intCast_zmod_eq_zero_iff_dvd]
    exact hd1
  obtain ⟨L, hLe⟩ : ∃ L : ZMod p,
      L = ((d1 : ℤ) : ZMod p) * (((n : ℤ) : ZMod p) * ((m1 : ℤ) : ZMod p)) := ⟨_, rfl⟩
  have hcast : (((d1 * (n * m1) : ℤ)) : ZMod p) = L := by
    rw [hLe, Int.cast_mul, Int.cast_mul]
  have hL0 : L ≠ 0 := by
    rw [hLe]
    exact mul_ne_zero hd10 (mul_ne_zero hn0 hm10)
  by_cases hsq : IsSquare L
  · obtain ⟨t, ht, ht0⟩ := exists_int_sq_of_isSquare (p := p) (c := d1 * (n * m1))
      (by rw [hcast]; exact hsq) (by rw [hcast]; exact hL0)
    exact ⟨0, prime_mul_passes_at_zero_one hp hp2 hn hm1 hr ht ht0 k⟩
  · have hL1 : L + 1 ≠ 0 := by
      intro h
      apply hsq
      have e : L = -1 := eq_neg_of_add_eq_zero_left h
      rw [e]
      exact (ZMod.exists_sq_eq_neg_one_iff (p := p)).mpr (by omega)
    have hL2 : L ^ 2 + 1 ≠ 0 := by
      intro h
      apply hsq
      have e : L ^ 2 = -1 := eq_neg_of_add_eq_zero_left h
      rw [ZMod.euler_criterion p hL0]
      have e2 : p / 2 = 2 * (2 * (p / 8)) := by omega
      rw [e2, pow_mul, e, pow_mul]
      norm_num
    obtain ⟨ω0, hω0e⟩ := ZMod.intCast_surjective (((m1 : ℤ) : ZMod p) * (L ^ 2 + 1) / (2 * L))
    obtain ⟨u1, hu1e⟩ := ZMod.intCast_surjective
      (((m1 : ℤ) : ZMod p) * (L ^ 2 - 1) / (2 * L * ((n : ℤ) : ZMod p)))
    obtain ⟨t, hte⟩ := ZMod.intCast_surjective (L + 1)
    have hω0 : ¬(p : ℤ) ∣ ω0 := by
      rw [← ZMod.intCast_zmod_eq_zero_iff_dvd, hω0e]
      exact div_ne_zero (mul_ne_zero hm10 hL2) (mul_ne_zero h2ne hL0)
    have ht0 : ¬(p : ℤ) ∣ t := by
      rw [← ZMod.intCast_zmod_eq_zero_iff_dvd, hte]
      exact hL1
    have hω : ω0 ^ 2 ≡ m1 ^ 2 + (n * u1) ^ 2 [ZMOD p] := by
      rw [← ZMod.intCast_eq_intCast_iff]
      push_cast
      rw [hω0e, hu1e]
      field_simp
      ring
    have ht : t ^ 2 ≡ d1 * (2 * n * (ω0 + m1)) [ZMOD p] := by
      rw [← ZMod.intCast_eq_intCast_iff]
      push_cast
      rw [hω0e, hte, hLe]
      field_simp
      ring
    exact ⟨u1, pass_odd_valuation_gen hp hp2 hn hω0 hω hr ht ht0 k⟩

end MagicSquareOfSquares.CoverDescent
