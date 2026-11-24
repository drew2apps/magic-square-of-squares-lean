import MagicSquareOfSquares.ClassOdd

/-!
# The class 1 at the prime 5 when 25 divides `m` (section 4.4 of `COMPLETE.md`)

When 5 divides `m` exactly once, only the classes `≡ ±2 (mod 5)` occur (`class_mod_five`). When 25
divides `m` the classes `≡ ±1 (mod 5)` occur too, at points with `u² + v² ≡ 0 (mod 25)`:

* `class_one_passes_at_five_of_125`: if 125 divides `m`, at the point `(1, 7)`.
* `class_one_passes_at_five_of_25`: if `m = 25·m₁` and `m₁·n` is a nonzero square modulo 5, at the
  point `(1, 57)`.

* `class_one_passes_at_five_of_25_nonres`: if `m = 25·m₁` and `m₁·n` is a non-residue modulo 5,
  at a point `(1, v)` that depends on `m, n` and on the power of 5.

With `passes_at_five` this is the row `p = 5` of the table of Theorem C, as lower bounds.
-/

namespace MagicSquareOfSquares.CoverDescent

/-- A point `(1, v)` at the prime 5 at which `GL` and `GU` are 100 times squares of units. -/
theorem five_pass_aux {m n v d AL AU c0 c8 cL cU g : ℤ}
    (hPL : nrm (PL m n 1 v) = 100 * AL) (hPU : nrm (PU m n 1 v) = 100 * AU)
    (e0 : c0 ^ 2 ≡ nrm (P0 m n 1 v) [ZMOD ((5 : ℕ) : ℤ)]) (h0 : ¬((5 : ℕ) : ℤ) ∣ c0)
    (e8 : c8 ^ 2 ≡ nrm (P8 m n 1 v) [ZMOD ((5 : ℕ) : ℤ)]) (h8 : ¬((5 : ℕ) : ℤ) ∣ c8)
    (eL : cL ^ 2 ≡ AL [ZMOD ((5 : ℕ) : ℤ)]) (hL : ¬((5 : ℕ) : ℤ) ∣ cL)
    (eU : cU ^ 2 ≡ AU [ZMOD ((5 : ℕ) : ℤ)]) (hU : ¬((5 : ℕ) : ℤ) ∣ cU)
    (hg : g ^ 2 ≡ d * (16 * ((c0 + m * v) * (c8 + n * v) * (10 * cL + (m + n) * (1 - v)) *
      (10 * cU + (m - n) * (1 - v)))) [ZMOD ((5 : ℕ) : ℤ)])
    (hgu : ¬((5 : ℕ) : ℤ) ∣ g) (k : ℕ) :
    ∀ f : Bool × Bool × Bool, GSqM (((5 : ℕ) : ℤ) ^ (k + 1)) d (Q m n 1 v f) := by
  have hp : Nat.Prime 5 := by norm_num
  have hp2 : (5 : ℕ) ≠ 2 := by norm_num
  have hp' : Prime (((5 : ℕ) : ℤ)) := Nat.prime_iff_prime_int.mp hp
  obtain ⟨w0, hw0, k0⟩ := exists_sq_modEq_pow hp hp2 h0 e0 k
  obtain ⟨w8, hw8, k8⟩ := exists_sq_modEq_pow hp hp2 h8 e8 k
  obtain ⟨ωL, hωL, kL⟩ := exists_sq_modEq_pow hp hp2 hL eL k
  obtain ⟨ωU, hωU, kU⟩ := exists_sq_modEq_pow hp hp2 hU eU k
  have hwL : (10 * ωL) ^ 2 ≡ nrm (PL m n 1 v) [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] := by
    rw [hPL]
    calc (10 * ωL) ^ 2 = 100 * ωL ^ 2 := by ring
      _ ≡ 100 * AL [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] := hωL.mul_left _
  have hwU : (10 * ωU) ^ 2 ≡ nrm (PU m n 1 v) [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] := by
    rw [hPU]
    calc (10 * ωU) ^ 2 = 100 * ωU ^ 2 := by ring
      _ ≡ 100 * AU [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] := hωU.mul_left _
  have cD : clsD m n 1 v w0 w8 (10 * ωL) (10 * ωU) ≡
      (c0 + m * v) * (c8 + n * v) * (10 * cL + (m + n) * (1 - v)) *
        (10 * cU + (m - n) * (1 - v)) [ZMOD ((5 : ℕ) : ℤ)] := by
    unfold clsD
    exact (((k0.add_right _).mul (k8.add_right _)).mul ((kL.mul_left 10).add_right _)).mul
      ((kU.mul_left 10).add_right _)
  have hE : g ^ 2 ≡ d * (16 * clsD m n 1 v w0 w8 (10 * ωL) (10 * ωU)) [ZMOD ((5 : ℕ) : ℤ)] :=
    hg.trans ((cD.symm.mul_left 16).mul_left d)
  obtain ⟨x, hx, hxg⟩ := exists_sq_modEq_pow hp hp2 hgu hE k
  have hxp : ¬((5 : ℕ) : ℤ) ∣ x := fun hd =>
    hgu (Int.modEq_zero_iff_dvd.mp (hxg.symm.trans (Int.modEq_zero_iff_dvd.mpr hd)))
  obtain ⟨ρ, μ, hρ⟩ : IsCoprime x (((5 : ℕ) : ℤ) ^ (k + 1)) :=
    (hp'.irreducible.coprime_iff_not_dvd.mpr hxp).symm.pow_right
  have hρ1 : ρ * x ≡ 1 [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] :=
    Int.modEq_iff_dvd.mpr ⟨μ, by linear_combination -hρ⟩
  have hd : d ≡ (ρ * d) ^ 2 * (16 * clsD m n 1 v w0 w8 (10 * ωL) (10 * ωU))
      [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] := by
    calc d = d * 1 ^ 2 := by ring
      _ ≡ d * (ρ * x) ^ 2 [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] := ((hρ1.pow 2).symm).mul_left d
      _ = d * ρ ^ 2 * x ^ 2 := by ring
      _ ≡ d * ρ ^ 2 * (d * (16 * clsD m n 1 v w0 w8 (10 * ωL) (10 * ωU)))
          [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] := hx.mul_left _
      _ = (ρ * d) ^ 2 * (16 * clsD m n 1 v w0 w8 (10 * ωL) (10 * ωU)) := by ring
  exact fun f => gsqM_class hw0.symm hw8.symm hwL.symm hwU.symm hd f

theorem exists_sqrt_mod_five {d : ℤ} (hd : d % 5 = 1 ∨ d % 5 = 4) :
    ∃ t : ℤ, t ^ 2 ≡ d [ZMOD ((5 : ℕ) : ℤ)] ∧ ¬((5 : ℕ) : ℤ) ∣ t := by
  rcases hd with h | h
  · refine ⟨1, ?_, by norm_num⟩
    unfold Int.ModEq
    push_cast
    omega
  · refine ⟨2, ?_, by norm_num⟩
    unfold Int.ModEq
    push_cast
    omega

/-- **The class 1 at 5 when 125 divides `m`**: every `d ≡ ±1 (mod 5)` passes at `(1, 7)` modulo
every power of 5. -/
theorem class_one_passes_at_five_of_125 {m3 n d : ℤ} (hn : ¬(5 : ℤ) ∣ n)
    (hd : d % 5 = 1 ∨ d % 5 = 4) (k : ℕ) :
    ∀ f : Bool × Bool × Bool, GSqM ((5 : ℤ) ^ (k + 1)) d (Q (125 * m3) n 1 7 f) := by
  have hp' : Prime (((5 : ℕ) : ℤ)) := Nat.prime_iff_prime_int.mp (by norm_num)
  have h5n : ¬((5 : ℕ) : ℤ) ∣ n := by exact_mod_cast hn
  obtain ⟨t, ht, ht0⟩ := exists_sqrt_mod_five hd
  have h7 : ¬((5 : ℕ) : ℤ) ∣ 7 * n := not_dvd_mul_of_prime hp' (by norm_num) h5n
  have hg : (n * n * t) ^ 2 ≡ d * (16 * ((n + 125 * m3 * 7) * (7 * n + n * 7) *
      (10 * n + (125 * m3 + n) * (1 - 7)) * (10 * n + (125 * m3 - n) * (1 - 7))))
      [ZMOD ((5 : ℕ) : ℤ)] := by
    have h1 : (n * n * t) ^ 2 ≡ n ^ 4 * d [ZMOD ((5 : ℕ) : ℤ)] := by
      calc (n * n * t) ^ 2 = n ^ 4 * t ^ 2 := by ring
        _ ≡ n ^ 4 * d [ZMOD ((5 : ℕ) : ℤ)] := ht.mul_left _
    refine h1.trans (Int.modEq_iff_dvd.mpr ⟨22050000000 * d * m3 ^ 3 * n -
      562800000 * d * m3 ^ 2 * n ^ 2 + 1836800 * d * m3 * n ^ 3 + 2867 * d * n ^ 4, ?_⟩)
    push_cast
    ring
  have key := five_pass_aux (m := 125 * m3) (n := n) (v := 7) (d := d)
    (AL := (125 * m3) ^ 2 + n ^ 2 - 70 * m3 * n) (AU := (125 * m3) ^ 2 + n ^ 2 + 70 * m3 * n)
    (c0 := n) (c8 := 7 * n) (cL := n) (cU := n) (g := n * n * t)
    (by simp only [nrm, PL]; ring) (by simp only [nrm, PU]; ring)
    (Int.modEq_iff_dvd.mpr ⟨153125 * m3 ^ 2, by simp only [nrm, P0]; push_cast; ring⟩) h5n
    (Int.modEq_iff_dvd.mpr ⟨3125 * m3 ^ 2, by simp only [nrm, P8]; push_cast; ring⟩) h7
    (Int.modEq_iff_dvd.mpr ⟨3125 * m3 ^ 2 - 14 * m3 * n, by push_cast; ring⟩) h5n
    (Int.modEq_iff_dvd.mpr ⟨3125 * m3 ^ 2 + 14 * m3 * n, by push_cast; ring⟩) h5n
    hg (not_dvd_mul_of_prime hp' (not_dvd_mul_of_prime hp' h5n h5n) ht0) k
  intro f
  have h := key f
  exact_mod_cast h

/-- **The class 1 at 5 when `m = 25·m₁` and `m₁·n` is a nonzero square modulo 5**: every
`d ≡ ±1 (mod 5)` passes at `(1, 57)` modulo every power of 5. -/
theorem class_one_passes_at_five_of_25 {m1 n τ d : ℤ} (hn : ¬(5 : ℤ) ∣ n)
    (hτ : τ ^ 2 ≡ m1 * n [ZMOD 5]) (hτ0 : ¬(5 : ℤ) ∣ τ) (hd : d % 5 = 1 ∨ d % 5 = 4) (k : ℕ) :
    ∀ f : Bool × Bool × Bool, GSqM ((5 : ℤ) ^ (k + 1)) d (Q (25 * m1) n 1 57 f) := by
  have hp' : Prime (((5 : ℕ) : ℤ)) := Nat.prime_iff_prime_int.mp (by norm_num)
  have h5n : ¬((5 : ℕ) : ℤ) ∣ n := by exact_mod_cast hn
  have h5τ : ¬((5 : ℕ) : ℤ) ∣ τ := by exact_mod_cast hτ0
  have hτ' : τ ^ 2 ≡ m1 * n [ZMOD ((5 : ℕ) : ℤ)] := by exact_mod_cast hτ
  obtain ⟨t, ht, ht0⟩ := exists_sqrt_mod_five hd
  have h57 : ¬((5 : ℕ) : ℤ) ∣ 57 * n := not_dvd_mul_of_prime hp' (by norm_num) h5n
  have h2τ : ¬((5 : ℕ) : ℤ) ∣ 2 * τ := not_dvd_mul_of_prime hp' (by norm_num) h5τ
  have eL : τ ^ 2 ≡ 65 * ((25 * m1) ^ 2 + n ^ 2) - 114 * m1 * n [ZMOD ((5 : ℕ) : ℤ)] :=
    hτ'.trans (Int.modEq_iff_dvd.mpr
      ⟨13 * ((25 * m1) ^ 2 + n ^ 2) - 23 * m1 * n, by push_cast; ring⟩)
  have eU : (2 * τ) ^ 2 ≡ 65 * ((25 * m1) ^ 2 + n ^ 2) + 114 * m1 * n [ZMOD ((5 : ℕ) : ℤ)] := by
    calc (2 * τ) ^ 2 = 4 * τ ^ 2 := by ring
      _ ≡ 4 * (m1 * n) [ZMOD ((5 : ℕ) : ℤ)] := hτ'.mul_left 4
      _ ≡ 65 * ((25 * m1) ^ 2 + n ^ 2) + 114 * m1 * n [ZMOD ((5 : ℕ) : ℤ)] :=
          Int.modEq_iff_dvd.mpr
            ⟨13 * ((25 * m1) ^ 2 + n ^ 2) + 22 * m1 * n, by push_cast; ring⟩
  have hg : (n * n * t) ^ 2 ≡ d * (16 * ((n + 25 * m1 * 57) * (57 * n + n * 57) *
      (10 * τ + (25 * m1 + n) * (1 - 57)) * (10 * (2 * τ) + (25 * m1 - n) * (1 - 57))))
      [ZMOD ((5 : ℕ) : ℤ)] := by
    have h1 : (n * n * t) ^ 2 ≡ n ^ 4 * d [ZMOD ((5 : ℕ) : ℤ)] := by
      calc (n * n * t) ^ 2 = n ^ 4 * t ^ 2 := by ring
        _ ≡ n ^ 4 * d [ZMOD ((5 : ℕ) : ℤ)] := ht.mul_left _
    refine h1.trans (Int.modEq_iff_dvd.mpr ⟨1018886400000 * d * m1 ^ 3 * n +
      715008000 * d * m1 ^ 2 * n ^ 2 - 21833280000 * d * m1 ^ 2 * n * τ -
      1630218240 * d * m1 * n ^ 3 - 306432000 * d * m1 * n ^ 2 * τ +
      103968000 * d * m1 * n * τ ^ 2 - 1144013 * d * n ^ 4 - 204288 * d * n ^ 3 * τ +
      72960 * d * n ^ 2 * τ ^ 2, ?_⟩)
    push_cast
    ring
  have key := five_pass_aux (m := 25 * m1) (n := n) (v := 57) (d := d)
    (AL := 65 * ((25 * m1) ^ 2 + n ^ 2) - 114 * m1 * n)
    (AU := 65 * ((25 * m1) ^ 2 + n ^ 2) + 114 * m1 * n)
    (c0 := n) (c8 := 57 * n) (cL := τ) (cU := 2 * τ) (g := n * n * t)
    (by simp only [nrm, PL]; ring) (by simp only [nrm, PU]; ring)
    (Int.modEq_iff_dvd.mpr ⟨406125 * m1 ^ 2, by simp only [nrm, P0]; push_cast; ring⟩) h5n
    (Int.modEq_iff_dvd.mpr ⟨125 * m1 ^ 2, by simp only [nrm, P8]; push_cast; ring⟩) h57
    eL h5τ eU h2τ
    hg (not_dvd_mul_of_prime hp' (not_dvd_mul_of_prime hp' h5n h5n) ht0) k
  intro f
  have h := key f
  exact_mod_cast h

/-- **The class 1 at 5 when `m = 25·m₁` and `m₁·n` is a non-residue modulo 5** (in the form
`3·m₁·n ≡ τ²`): every `d ≡ ±1 (mod 5)` passes modulo every power of 5 at a point `(1, v)`, where
`v` solves `GL(1, v) ≡ 625` (the quadratic formula, with a square root by Hensel's lemma). -/
theorem class_one_passes_at_five_of_25_nonres {m1 n τ d : ℤ} (hn : ¬(5 : ℤ) ∣ n)
    (hτ : τ ^ 2 ≡ 3 * (m1 * n) [ZMOD 5]) (hτ0 : ¬(5 : ℤ) ∣ τ) (hd : d % 5 = 1 ∨ d % 5 = 4)
    (k : ℕ) :
    ∃ v : ℤ, ∀ f : Bool × Bool × Bool, GSqM ((5 : ℤ) ^ (k + 1)) d (Q (25 * m1) n 1 v f) := by
  have hp : Nat.Prime 5 := by norm_num
  have hp2 : (5 : ℕ) ≠ 2 := by norm_num
  have hp' : Prime (((5 : ℕ) : ℤ)) := Nat.prime_iff_prime_int.mp hp
  have h5n : ¬((5 : ℕ) : ℤ) ∣ n := by exact_mod_cast hn
  have h5τ : ¬((5 : ℕ) : ℤ) ∣ τ := by exact_mod_cast hτ0
  have hτ' : τ ^ 2 ≡ 3 * (m1 * n) [ZMOD ((5 : ℕ) : ℤ)] := by exact_mod_cast hτ
  obtain ⟨t, ht, ht0⟩ := exists_sqrt_mod_five hd
  have hnn : ¬((5 : ℕ) : ℤ) ∣ n * n := not_dvd_mul_of_prime hp' h5n h5n
  -- the square root of the discriminant
  have hΔ : (n * n) ^ 2 ≡ 1250 * ((25 * m1) ^ 2 + n ^ 2) - 4 * ((25 * m1) ^ 2 - n ^ 2) ^ 2
      [ZMOD ((5 : ℕ) : ℤ)] :=
    Int.modEq_iff_dvd.mpr ⟨-312500 * m1 ^ 4 + 1000 * m1 ^ 2 * n ^ 2 + 156250 * m1 ^ 2 - n ^ 4 +
      250 * n ^ 2, by push_cast; ring⟩
  obtain ⟨s, hs, ks⟩ := exists_sq_modEq_pow hp hp2 hnn hΔ k
  have h2N : ¬((5 : ℕ) : ℤ) ∣ 2 * ((25 * m1) ^ 2 + n ^ 2) := by
    intro h
    have h1 : ((5 : ℕ) : ℤ) ∣ 2 * ((25 * m1) ^ 2 + n ^ 2) - ((5 : ℕ) : ℤ) * (250 * m1 ^ 2) :=
      dvd_sub h (dvd_mul_right _ _)
    have e : 2 * ((25 * m1) ^ 2 + n ^ 2) - ((5 : ℕ) : ℤ) * (250 * m1 ^ 2) = 2 * (n * n) := by
      push_cast
      ring
    rw [e] at h1
    exact not_dvd_mul_of_prime hp' (not_dvd_two_of_odd hp hp2) hnn h1
  obtain ⟨ι, μ, hι⟩ : IsCoprime (2 * ((25 * m1) ^ 2 + n ^ 2)) (((5 : ℕ) : ℤ) ^ (k + 1)) :=
    (hp'.irreducible.coprime_iff_not_dvd.mpr h2N).symm.pow_right
  have hι1 : ι * (2 * ((25 * m1) ^ 2 + n ^ 2)) ≡ 1 [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] :=
    Int.modEq_iff_dvd.mpr ⟨μ, by linear_combination -hι⟩
  obtain ⟨v, hv⟩ : ∃ v : ℤ, v = ι * (s + 4 * (25 * m1) * n) := ⟨_, rfl⟩
  have hlin : 2 * ((25 * m1) ^ 2 + n ^ 2) * v ≡ s + 4 * (25 * m1) * n
      [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] := by
    calc 2 * ((25 * m1) ^ 2 + n ^ 2) * v
        = ι * (2 * ((25 * m1) ^ 2 + n ^ 2)) * (s + 4 * (25 * m1) * n) := by rw [hv]; ring
      _ ≡ 1 * (s + 4 * (25 * m1) * n) [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] := hι1.mul_right _
      _ = s + 4 * (25 * m1) * n := by ring
  -- `GL(1, v) ≡ 625`
  have id1 : 2 * ((25 * m1) ^ 2 + n ^ 2) * (nrm (PL (25 * m1) n 1 v) - 625) =
      (2 * ((25 * m1) ^ 2 + n ^ 2) * v - 4 * (25 * m1) * n) ^ 2 -
        (1250 * ((25 * m1) ^ 2 + n ^ 2) - 4 * ((25 * m1) ^ 2 - n ^ 2) ^ 2) := by
    simp only [nrm, PL]
    ring
  have hA : (2 * ((25 * m1) ^ 2 + n ^ 2) * v - 4 * (25 * m1) * n) ^ 2 ≡
      1250 * ((25 * m1) ^ 2 + n ^ 2) - 4 * ((25 * m1) ^ 2 - n ^ 2) ^ 2
      [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] := by
    have h1 : 2 * ((25 * m1) ^ 2 + n ^ 2) * v - 4 * (25 * m1) * n ≡ s
        [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] := by
      have h2 := hlin.sub_right (4 * (25 * m1) * n)
      have e : s + 4 * (25 * m1) * n - 4 * (25 * m1) * n = s := by ring
      rwa [e] at h2
    exact (h1.pow 2).trans hs
  have hB : 2 * ((25 * m1) ^ 2 + n ^ 2) * (nrm (PL (25 * m1) n 1 v) - 625) ≡ 0
      [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] := by
    rw [id1]
    have h1 := hA.sub_right (1250 * ((25 * m1) ^ 2 + n ^ 2) - 4 * ((25 * m1) ^ 2 - n ^ 2) ^ 2)
    rwa [sub_self] at h1
  have hGL : nrm (PL (25 * m1) n 1 v) ≡ 25 ^ 2 [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] := by
    calc nrm (PL (25 * m1) n 1 v) = 625 + 1 * (nrm (PL (25 * m1) n 1 v) - 625) := by ring
      _ ≡ 625 + ι * (2 * ((25 * m1) ^ 2 + n ^ 2)) * (nrm (PL (25 * m1) n 1 v) - 625)
          [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] := (hι1.symm.mul_right _).add_left _
      _ = 625 + ι * (2 * ((25 * m1) ^ 2 + n ^ 2) * (nrm (PL (25 * m1) n 1 v) - 625)) := by ring
      _ ≡ 625 + ι * 0 [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] := (hB.mul_left ι).add_left _
      _ = 25 ^ 2 := by norm_num
  -- `v ≡ 3 (mod 5)`
  have hv3 : v ≡ 3 [ZMOD ((5 : ℕ) : ℤ)] := by
    have hlin5 := hlin.of_dvd (dvd_pow_self (((5 : ℕ) : ℤ)) (Nat.succ_ne_zero k))
    obtain ⟨a, ha⟩ := hlin5.dvd
    obtain ⟨b, hb⟩ := ks.dvd
    have h1 : ((5 : ℕ) : ℤ) ∣ n * n * (2 * v - 1) :=
      ⟨-a - b - 250 * m1 ^ 2 * v + 20 * m1 * n, by
        push_cast at ha hb ⊢
        linear_combination -ha - hb⟩
    have h2 : ((5 : ℕ) : ℤ) ∣ 2 * v - 1 := (hp'.dvd_or_dvd h1).resolve_left hnn
    have h3 : (5 : ℤ) ∣ 2 * v - 1 := by exact_mod_cast h2
    have h4 : (5 : ℤ) ∣ 3 - v := by omega
    exact Int.modEq_iff_dvd.mpr (by exact_mod_cast h4)
  have h5v : ¬((5 : ℕ) : ℤ) ∣ v := by
    intro h
    have h1 : (5 : ℤ) ∣ v := by exact_mod_cast h
    have h2 : (5 : ℤ) ∣ 3 - v := by exact_mod_cast hv3.dvd
    omega
  -- `GU(1, v) ≡ 25·(25 + 16 m₁ n v)`
  have hPU : nrm (PU (25 * m1) n 1 v) ≡ 25 * (25 + 16 * m1 * n * v)
      [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] := by
    calc nrm (PU (25 * m1) n 1 v) = nrm (PL (25 * m1) n 1 v) + 400 * m1 * n * v := by
          simp only [nrm, PL, PU]
          ring
      _ ≡ 25 ^ 2 + 400 * m1 * n * v [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] := hGL.add_right _
      _ = 25 * (25 + 16 * m1 * n * v) := by ring
  have eU : (4 * τ) ^ 2 ≡ 25 + 16 * m1 * n * v [ZMOD ((5 : ℕ) : ℤ)] := by
    calc (4 * τ) ^ 2 = 16 * τ ^ 2 := by ring
      _ ≡ 16 * (3 * (m1 * n)) [ZMOD ((5 : ℕ) : ℤ)] := hτ'.mul_left 16
      _ = 16 * m1 * n * 3 := by ring
      _ ≡ 16 * m1 * n * v [ZMOD ((5 : ℕ) : ℤ)] := hv3.symm.mul_left _
      _ ≡ 25 + 16 * m1 * n * v [ZMOD ((5 : ℕ) : ℤ)] :=
          Int.modEq_iff_dvd.mpr ⟨5, by push_cast; ring⟩
  obtain ⟨ωU, hωU, kU⟩ := exists_sq_modEq_pow hp hp2
    (not_dvd_mul_of_prime hp' (by norm_num) h5τ) eU k
  have hwU : (5 * ωU) ^ 2 ≡ nrm (PU (25 * m1) n 1 v) [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] := by
    calc (5 * ωU) ^ 2 = 25 * ωU ^ 2 := by ring
      _ ≡ 25 * (25 + 16 * m1 * n * v) [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] := hωU.mul_left _
      _ ≡ nrm (PU (25 * m1) n 1 v) [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] := hPU.symm
  have e0 : n ^ 2 ≡ nrm (P0 (25 * m1) n 1 v) [ZMOD ((5 : ℕ) : ℤ)] :=
    Int.modEq_iff_dvd.mpr ⟨125 * m1 ^ 2 * v ^ 2, by simp only [nrm, P0]; push_cast; ring⟩
  have e8 : (n * v) ^ 2 ≡ nrm (P8 (25 * m1) n 1 v) [ZMOD ((5 : ℕ) : ℤ)] :=
    Int.modEq_iff_dvd.mpr ⟨125 * m1 ^ 2, by simp only [nrm, P8]; push_cast; ring⟩
  obtain ⟨w0, hw0, k0⟩ := exists_sq_modEq_pow hp hp2 h5n e0 k
  obtain ⟨w8, hw8, k8⟩ := exists_sq_modEq_pow hp hp2 (not_dvd_mul_of_prime hp' h5n h5v) e8 k
  have one3 : (1 : ℤ) - v ≡ 1 - 3 [ZMOD ((5 : ℕ) : ℤ)] := (Int.ModEq.refl 1).sub hv3
  have cD : clsD (25 * m1) n 1 v w0 w8 25 (5 * ωU) ≡
      (n + 25 * m1 * 3) * (n * 3 + n * 3) * (25 + (25 * m1 + n) * (1 - 3)) *
        (5 * (4 * τ) + (25 * m1 - n) * (1 - 3)) [ZMOD ((5 : ℕ) : ℤ)] := by
    unfold clsD
    exact (((k0.add (hv3.mul_left _)).mul ((k8.trans (hv3.mul_left n)).add (hv3.mul_left n))).mul
      ((one3.mul_left _).add_left 25)).mul ((kU.mul_left 5).add (one3.mul_left _))
  have hE : (n * n * t) ^ 2 ≡ d * (16 * clsD (25 * m1) n 1 v w0 w8 25 (5 * ωU))
      [ZMOD ((5 : ℕ) : ℤ)] := by
    have h1 : (n * n * t) ^ 2 ≡ n ^ 4 * d [ZMOD ((5 : ℕ) : ℤ)] := by
      calc (n * n * t) ^ 2 = n ^ 4 * t ^ 2 := by ring
        _ ≡ n ^ 4 * d [ZMOD ((5 : ℕ) : ℤ)] := ht.mul_left _
    have h2 : n ^ 4 * d ≡ d * (16 * ((n + 25 * m1 * 3) * (n * 3 + n * 3) *
        (25 + (25 * m1 + n) * (1 - 3)) * (5 * (4 * τ) + (25 * m1 - n) * (1 - 3))))
        [ZMOD ((5 : ℕ) : ℤ)] :=
      Int.modEq_iff_dvd.mpr ⟨3600000 * d * m1 ^ 3 * n + 48000 * d * m1 ^ 2 * n ^ 2 -
        1440000 * d * m1 ^ 2 * n * τ - 1800000 * d * m1 ^ 2 * n - 5760 * d * m1 * n ^ 3 -
        76800 * d * m1 * n ^ 2 * τ + 48000 * d * m1 * n ^ 2 + 720000 * d * m1 * n * τ -
        77 * d * n ^ 4 - 768 * d * n ^ 3 * τ + 960 * d * n ^ 3 + 9600 * d * n ^ 2 * τ, by
          push_cast
          ring⟩
    exact (h1.trans h2).trans ((cD.symm.mul_left 16).mul_left d)
  have hgu : ¬((5 : ℕ) : ℤ) ∣ n * n * t := not_dvd_mul_of_prime hp' hnn ht0
  obtain ⟨x, hx, hxg⟩ := exists_sq_modEq_pow hp hp2 hgu hE k
  have hxp : ¬((5 : ℕ) : ℤ) ∣ x := fun hd =>
    hgu (Int.modEq_zero_iff_dvd.mp (hxg.symm.trans (Int.modEq_zero_iff_dvd.mpr hd)))
  obtain ⟨ρ, μ', hρ⟩ : IsCoprime x (((5 : ℕ) : ℤ) ^ (k + 1)) :=
    (hp'.irreducible.coprime_iff_not_dvd.mpr hxp).symm.pow_right
  have hρ1 : ρ * x ≡ 1 [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] :=
    Int.modEq_iff_dvd.mpr ⟨μ', by linear_combination -hρ⟩
  have hdd : d ≡ (ρ * d) ^ 2 * (16 * clsD (25 * m1) n 1 v w0 w8 25 (5 * ωU))
      [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] := by
    calc d = d * 1 ^ 2 := by ring
      _ ≡ d * (ρ * x) ^ 2 [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] := ((hρ1.pow 2).symm).mul_left d
      _ = d * ρ ^ 2 * x ^ 2 := by ring
      _ ≡ d * ρ ^ 2 * (d * (16 * clsD (25 * m1) n 1 v w0 w8 25 (5 * ωU)))
          [ZMOD ((5 : ℕ) : ℤ) ^ (k + 1)] := hx.mul_left _
      _ = (ρ * d) ^ 2 * (16 * clsD (25 * m1) n 1 v w0 w8 25 (5 * ωU)) := by ring
  refine ⟨v, fun f => ?_⟩
  have h := gsqM_class hw0.symm hw8.symm hGL hwU.symm hdd f
  exact_mod_cast h

end MagicSquareOfSquares.CoverDescent
