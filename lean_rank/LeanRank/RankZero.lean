import EllipticCurves.SelmerGroup
import EllipticCurves.MordellWeil

/-!
# `rank E'(ℚ) = 0` for `E' : y² = x(x − 12)(x − 16)`, through the `x − T` map, with the bad place `2`

Written by an AI model; not read by a person.  Experiment, not counted.

`E'` is the curve `y² = x(x² − 28x + 192)` of the proof of Lemma B of `UNIFORM_RANK.md`: it is
2-isogenous to `E_B : y² = x(x² + 14x + 1)`, the curve of Euler's equation
`X⁴ + 14X²Y² + Y⁴ = Z²`.  It has full rational 2-torsion (`x² − 28x + 192 = (x − 12)(x − 16)`), so
the étale algebra `A = ℚ[X]/(f)` of M. Stoll's `x − T` map `μ : E'(ℚ) → Aˣ/(Aˣ)²` is `ℚ³`.

Inputs from M. Stoll's library `EllipticCurves` (commit 2164fb5):
* `fg_point_of_numberField` (Mordell–Weil);
* `μ` and `card_range_μ`: `#(im μ) = 2 ^ rank · #E(ℚ)[2]`, which rests on `ker μ = 2E(ℚ)`
  (`ker_μ_eq`) and `[G : 2G] = 2 ^ rank · #G[2]`;
* `card_ker_nsmul_two`: `#E(ℚ)[2] = #(roots of f) + 1`.

Proved here:
* `Ψ : Aˣ/(Aˣ)² → (ℚˣ/(ℚˣ)²)³`, evaluation at the roots `0, 12, 16`, is injective (Lagrange
  interpolation);
* `descent`: for a rational point with `x(x − 12)(x − 16) ≠ 0`, the triple of square classes of
  `(x, x − 12, x − 16)` is one of `(1, 1, 1)`, `(3, −3, −1)`, `(2, −1, −2)`, `(6, 3, 2)`.  The proof
  writes `x = n/e²` (the denominator is a square), splits off the square classes of `n` and
  `n − 12e²` (only `2` and `3` can occur, by coprimality), uses the real place (`n > 0`), and
  excludes the other 28 classes by one certificate modulo `64` (`badClasses_excluded`, checked by
  `decide +kernel`; soundness `excluded_eq_false`): **the local condition at the bad place 2**.
* Hence `#(im μ) ≤ 4 = #E(ℚ)[2]`, so `rank E'(ℚ) = 0` (`finrank_eq_zero`).
-/

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

namespace LeanRank.RankZero

/-! ### 1. The certificate modulo `64` -/

/-- The squares modulo `64`. -/
def sq64 : List ℤ := [0, 1, 4, 9, 16, 17, 25, 33, 36, 41, 49, 57]

theorem sq64_mem_aux : ∀ k : ℕ, k < 64 → ((k : ℤ) * k) % 64 ∈ sq64 := by decide +kernel

theorem sq64_mem (z : ℤ) : z ^ 2 % 64 ∈ sq64 := by
  have h0 : 0 ≤ z % 64 := Int.emod_nonneg z (by norm_num)
  have h1 : z % 64 < 64 := Int.emod_lt_of_pos z (by norm_num)
  obtain ⟨k, hk⟩ := Int.eq_ofNat_of_zero_le h0
  have hk64 : k < 64 := by omega
  have := sq64_mem_aux k hk64
  rw [← hk] at this
  rw [sq, Int.mul_emod]
  exact this

/-- `pass d₁ d₂ s q` is `true` when the residues `s ≡ a²`, `q ≡ e² (mod 64)` are compatible with
`n = d₁a²`, `n − 12e² = d₂ · square`, `d₁d₂(n − 16e²) = square` (mod `64`) and `¬ (2 ∣ n ∧ 2 ∣ e)`. -/
def pass (d₁ d₂ s q : ℤ) : Bool :=
  ((d₁ * s) % 2 == 1 || q % 2 == 1) &&
    sq64.any (fun r => (d₁ * s - 12 * q - d₂ * r) % 64 == 0) &&
    sq64.any (fun r => (d₁ * d₂ * (d₁ * s - 16 * q) - r) % 64 == 0)

/-- The class `(d₁, d₂)` is excluded at `2` when no pair of residues passes. -/
def excluded (d₁ d₂ : ℤ) : Bool := !(sq64.any fun s => sq64.any fun q => pass d₁ d₂ s q)

/-- **Soundness of the certificate**: an integer solution with `¬ (2 ∣ d₁a² ∧ 2 ∣ e)` makes the
class pass, so a class with `excluded d₁ d₂ = true` has no such solution. -/
theorem excluded_eq_false {d₁ d₂ a e b t : ℤ} (h1 : d₁ * a ^ 2 - 12 * e ^ 2 = d₂ * b ^ 2)
    (h2 : d₁ * d₂ * (d₁ * a ^ 2 - 16 * e ^ 2) = t ^ 2) (hp : ¬ (2 ∣ d₁ * a ^ 2 ∧ 2 ∣ e)) :
    excluded d₁ d₂ = false := by
  have hA := Int.mod_modEq (a ^ 2) 64
  have hE := Int.mod_modEq (e ^ 2) 64
  have hB := Int.mod_modEq (b ^ 2) 64
  have hT := Int.mod_modEq (t ^ 2) 64
  have k1 : (d₁ * (a ^ 2 % 64) - 12 * (e ^ 2 % 64) - d₂ * (b ^ 2 % 64)) % 64 = 0 := by
    have := ((hA.mul_left d₁).sub (hE.mul_left 12)).sub (hB.mul_left d₂)
    have h0 : d₁ * a ^ 2 - 12 * e ^ 2 - d₂ * b ^ 2 = 0 := by rw [h1]; ring
    rw [h0] at this
    exact this
  have k2 : (d₁ * d₂ * (d₁ * (a ^ 2 % 64) - 16 * (e ^ 2 % 64)) - t ^ 2 % 64) % 64 = 0 := by
    have := (((hA.mul_left d₁).sub (hE.mul_left 16)).mul_left (d₁ * d₂)).sub hT
    have h0 : d₁ * d₂ * (d₁ * a ^ 2 - 16 * e ^ 2) - t ^ 2 = 0 := by rw [h2]; ring
    rw [h0] at this
    exact this
  have k3 : (d₁ * (a ^ 2 % 64)) % 2 = 1 ∨ (e ^ 2 % 64) % 2 = 1 := by
    have m1 : (d₁ * (a ^ 2 % 64)) % 2 = (d₁ * a ^ 2) % 2 :=
      (hA.mul_left d₁).of_dvd (by norm_num)
    have m2 : (e ^ 2 % 64) % 2 = (e ^ 2) % 2 := hE.of_dvd (by norm_num)
    by_cases hd : 2 ∣ d₁ * a ^ 2
    · right
      have he : ¬ 2 ∣ e := fun h => hp ⟨hd, h⟩
      have ho : Odd e := Int.not_even_iff_odd.mp (by rwa [even_iff_two_dvd])
      have : e ^ 2 % 2 = 1 := Int.odd_iff.mp ho.pow
      omega
    · left
      omega
  have hpass : pass d₁ d₂ (a ^ 2 % 64) (e ^ 2 % 64) = true := by
    simp only [pass, Bool.and_eq_true, Bool.or_eq_true, beq_iff_eq, List.any_eq_true]
    exact ⟨⟨k3, b ^ 2 % 64, sq64_mem b, k1⟩, t ^ 2 % 64, sq64_mem t, k2⟩
  have hany : (sq64.any fun s => sq64.any fun q => pass d₁ d₂ s q) = true :=
    List.any_eq_true.mpr ⟨_, sq64_mem a, List.any_eq_true.mpr ⟨_, sq64_mem e, hpass⟩⟩
  simp [excluded, hany]

/-- The four classes `(d₁, d₂)` that survive (they are the classes of the torsion points). -/
def goodClasses : List (ℤ × ℤ) := [(1, 1), (3, -3), (2, -1), (6, 3)]

/-- The other 28 classes with `d₁ > 0`. -/
def badClasses : List (ℤ × ℤ) :=
  [(1, -1), (1, 2), (1, -2), (1, 3), (1, -3), (1, 6), (1, -6),
   (2, 1), (2, 2), (2, -2), (2, 3), (2, -3), (2, 6), (2, -6),
   (3, 1), (3, -1), (3, 2), (3, -2), (3, 3), (3, 6), (3, -6),
   (6, 1), (6, -1), (6, 2), (6, -2), (6, -3), (6, 6), (6, -6)]

/-- **The 2-adic certificate**: each of the 28 classes is excluded modulo `64`. -/
theorem badClasses_excluded : badClasses.all (fun c => excluded c.1 c.2) = true := by
  decide +kernel

theorem classes_cover : ∀ d₁ ∈ ([1, 2, 3, 6] : List ℤ),
    ∀ d₂ ∈ ([1, -1, 2, -2, 3, -3, 6, -6] : List ℤ),
      (d₁, d₂) ∈ goodClasses ∨ (d₁, d₂) ∈ badClasses := by
  decide

/-! ### 2. Integer lemmas -/

theorem int_sq_of_rat_sq {r : ℚ} {N : ℤ} (h : r ^ 2 = N) : ∃ z : ℤ, z ^ 2 = N := by
  have : IsSquare (N : ℚ) := ⟨r, by rw [← h, sq]⟩
  obtain ⟨z, hz⟩ := Rat.isSquare_intCast_iff.mp this
  exact ⟨z, by rw [hz, sq]⟩

theorem prime_23 {p : ℕ} (hp : p.Prime) {i j : ℕ} (h : p ∣ 2 ^ i * 3 ^ j) : p = 2 ∨ p = 3 := by
  rcases (Nat.Prime.dvd_mul hp).mp h with h | h
  · exact Or.inl ((Nat.prime_dvd_prime_iff_eq hp Nat.prime_two).mp (hp.dvd_of_dvd_pow h))
  · exact Or.inr ((Nat.prime_dvd_prime_iff_eq hp Nat.prime_three).mp (hp.dvd_of_dvd_pow h))

/-- If `g · G` is a square and every prime dividing both `g` and `G` is `2` or `3`, then
`g = d · a²` with `d ∣ 6` up to sign. -/
theorem split6 {g G Z : ℤ} (hg : g ≠ 0) (hZ : g * G = Z ^ 2)
    (hcop : ∀ p : ℕ, p.Prime → (p : ℤ) ∣ g → (p : ℤ) ∣ G → p = 2 ∨ p = 3) :
    ∃ d a : ℤ, d ∈ ([1, -1, 2, -2, 3, -3, 6, -6] : List ℤ) ∧ g = d * a ^ 2 := by
  have hN : g.natAbs ≠ 0 := Int.natAbs_ne_zero.mpr hg
  obtain ⟨α, N₁, h2N₁, hN₁⟩ := Nat.exists_eq_pow_mul_and_not_dvd hN 2 (by norm_num)
  have hN₁0 : N₁ ≠ 0 := by
    rintro rfl
    exact hN (by rw [hN₁, mul_zero])
  obtain ⟨β, N₂, h3N₂, hN₂⟩ := Nat.exists_eq_pow_mul_and_not_dvd hN₁0 3 (by norm_num)
  have habs : g.natAbs = 2 ^ α * 3 ^ β * N₂ := by rw [hN₁, hN₂, mul_assoc]
  obtain ⟨σ, hσ, hgσ⟩ : ∃ σ : ℤ, (σ = 1 ∨ σ = -1) ∧ g = σ * (2 ^ α * 3 ^ β * (N₂ : ℤ)) := by
    have hc : ((g.natAbs : ℕ) : ℤ) = 2 ^ α * 3 ^ β * (N₂ : ℤ) := by rw [habs]; push_cast; ring
    rcases Int.natAbs_eq g with h | h
    · exact ⟨1, Or.inl rfl, by rw [one_mul, h, hc]⟩
    · exact ⟨-1, Or.inr rfl, by rw [h, hc]; ring⟩
  -- `N₂` is coprime to `σ · 2^α · 3^β · G`
  have hcopN : Int.gcd (N₂ : ℤ) (σ * 2 ^ α * 3 ^ β * G) = 1 := by
    rw [Int.gcd_eq_natAbs, Int.natAbs_natCast]
    apply Nat.coprime_of_dvd
    intro k hk hkN hkG
    have hk2 : k ≠ 2 := by rintro rfl; exact h2N₁ (by rw [hN₂]; exact Dvd.dvd.mul_left hkN _)
    have hk3 : k ≠ 3 := by rintro rfl; exact h3N₂ hkN
    rw [Int.natAbs_mul, Int.natAbs_mul, Int.natAbs_mul, Int.natAbs_pow, Int.natAbs_pow] at hkG
    have hσ1 : σ.natAbs = 1 := by rcases hσ with rfl | rfl <;> rfl
    rw [hσ1, one_mul] at hkG
    rcases (Nat.Prime.dvd_mul hk).mp hkG with h | h
    · rcases prime_23 hk (i := α) (j := β) (by simpa using h) with h | h
      · exact hk2 h
      · exact hk3 h
    · have hkg : (k : ℤ) ∣ g := by
        rw [Int.natCast_dvd, habs]
        exact Dvd.dvd.mul_left hkN _
      have hkG' : (k : ℤ) ∣ G := Int.natCast_dvd.mpr h
      rcases hcop k hk hkg hkG' with h | h
      · exact hk2 h
      · exact hk3 h
  have hprod : (N₂ : ℤ) * (σ * 2 ^ α * 3 ^ β * G) = Z ^ 2 := by
    rw [← hZ, hgσ]; ring
  obtain ⟨c, hc⟩ := Int.sq_of_gcd_eq_one hcopN hprod
  obtain ⟨s, hs, hsc⟩ : ∃ s : ℤ, (s = 1 ∨ s = -1) ∧ (N₂ : ℤ) = s * c ^ 2 := by
    rcases hc with h | h
    · exact ⟨1, Or.inl rfl, by rw [h, one_mul]⟩
    · exact ⟨-1, Or.inr rfl, by rw [h]; ring⟩
  have e2 : (2 : ℤ) ^ α = 2 ^ (α % 2) * (2 ^ (α / 2)) ^ 2 := by
    rw [← pow_mul, ← pow_add]; congr 1; omega
  have e3 : (3 : ℤ) ^ β = 3 ^ (β % 2) * (3 ^ (β / 2)) ^ 2 := by
    rw [← pow_mul, ← pow_add]; congr 1; omega
  refine ⟨σ * s * 2 ^ (α % 2) * 3 ^ (β % 2), 2 ^ (α / 2) * 3 ^ (β / 2) * c, ?_, ?_⟩
  · rcases hσ with rfl | rfl <;> rcases hs with rfl | rfl <;>
      rcases Nat.mod_two_eq_zero_or_one α with h | h <;>
      rcases Nat.mod_two_eq_zero_or_one β with h' | h' <;> simp [h, h']
  · rw [hgσ, hsc, e2, e3]; ring

/-! ### 3. The descent on `y² = x(x − 12)(x − 16)` -/

/-- **Descent**: for a rational point with `x(x − 12)(x − 16) ≠ 0`, the classes of `x`, `x − 12`,
`x − 16` modulo squares are `d₁`, `d₂`, `d₁d₂` with `(d₁, d₂)` one of the four good classes. -/
theorem descent {x y : ℚ} (h : y ^ 2 = x * (x - 12) * (x - 16)) (hx0 : x ≠ 0)
    (hx12 : x - 12 ≠ 0) (hx16 : x - 16 ≠ 0) :
    ∃ d₁ d₂ : ℤ, (d₁, d₂) ∈ goodClasses ∧ ∃ w₁ w₂ w₃ : ℚ, w₁ ≠ 0 ∧ w₂ ≠ 0 ∧ w₃ ≠ 0 ∧
      x = d₁ * w₁ ^ 2 ∧ x - 12 = d₂ * w₂ ^ 2 ∧ x - 16 = ((d₁ * d₂ : ℤ) : ℚ) * w₃ ^ 2 := by
  -- `x = n / q` in lowest terms
  obtain ⟨n, q, hq0, hxq, hcop⟩ : ∃ n q : ℤ, 0 < q ∧ x * q = n ∧ IsCoprime n q :=
    ⟨x.num, x.den, by exact_mod_cast x.den_pos, by exact_mod_cast Rat.mul_den_eq_num x,
      Int.isCoprime_iff_gcd_eq_one.mpr x.reduced⟩
  -- `q · n(n − 12q)(n − 16q)` is a square
  obtain ⟨Y, hY⟩ : ∃ Y : ℤ, Y ^ 2 = q * (n * (n - 12 * q) * (n - 16 * q)) := by
    apply int_sq_of_rat_sq (r := y * q ^ 2)
    push_cast
    linear_combination (q : ℚ) ^ 4 * h +
      (q : ℚ) * ((x * q) ^ 2 + x * q * n + n ^ 2 - 28 * q * (x * q + n) + 192 * q ^ 2) * hxq
  -- so `q` is a square
  obtain ⟨u, v, huv⟩ := hcop
  have hcop2 : IsCoprime q (n * (n - 12 * q) * (n - 16 * q)) := by
    have c1 : IsCoprime q n := ⟨v, u, by linear_combination huv⟩
    have c2 : IsCoprime q (n - 12 * q) := ⟨v + 12 * u, u, by linear_combination huv⟩
    have c3 : IsCoprime q (n - 16 * q) := ⟨v + 16 * u, u, by linear_combination huv⟩
    exact (c1.mul_right c2).mul_right c3
  obtain ⟨e, he⟩ := Int.sq_of_gcd_eq_one (Int.isCoprime_iff_gcd_eq_one.mp hcop2) hY.symm
  have hqe : q = e ^ 2 := by
    rcases he with he | he
    · exact he
    · nlinarith [sq_nonneg e]
  have he0 : e ≠ 0 := by rintro rfl; rw [hqe] at hq0; norm_num at hq0
  have hne : IsCoprime n e := ⟨u, v * e, by rw [hqe] at huv; linear_combination huv⟩
  have hxe : x * (e : ℚ) ^ 2 = n := by rw [← hxq, hqe]; push_cast; ring
  -- `n(n − 12e²)(n − 16e²) = Z²`
  obtain ⟨Z, hZ⟩ : ∃ Z : ℤ, Z ^ 2 = n * (n - 12 * e ^ 2) * (n - 16 * e ^ 2) := by
    apply int_sq_of_rat_sq (r := y * e ^ 3)
    push_cast
    linear_combination (e : ℚ) ^ 6 * h +
      ((x * e ^ 2) ^ 2 + x * e ^ 2 * n + n ^ 2 - 28 * e ^ 2 * (x * e ^ 2 + n) +
        192 * (e : ℚ) ^ 4) * hxe
  have hn0 : n ≠ 0 := by
    intro h0
    have : x * (e : ℚ) ^ 2 = 0 := by rw [hxe, h0]; simp
    rcases mul_eq_zero.mp this with h1 | h1
    · exact hx0 h1
    · exact he0 (by exact_mod_cast pow_eq_zero_iff (n := 2) (by norm_num) |>.mp h1)
  have hg2 : n - 12 * e ^ 2 ≠ 0 := by
    intro h0
    have : (x - 12) * (e : ℚ) ^ 2 = 0 := by
      have : ((n - 12 * e ^ 2 : ℤ) : ℚ) = 0 := by rw [h0]; simp
      push_cast at this
      linear_combination this + hxe
    rcases mul_eq_zero.mp this with h1 | h1
    · exact hx12 h1
    · exact he0 (by exact_mod_cast pow_eq_zero_iff (n := 2) (by norm_num) |>.mp h1)
  have hg3 : n - 16 * e ^ 2 ≠ 0 := by
    intro h0
    have : (x - 16) * (e : ℚ) ^ 2 = 0 := by
      have : ((n - 16 * e ^ 2 : ℤ) : ℚ) = 0 := by rw [h0]; simp
      push_cast at this
      linear_combination this + hxe
    rcases mul_eq_zero.mp this with h1 | h1
    · exact hx16 h1
    · exact he0 (by exact_mod_cast pow_eq_zero_iff (n := 2) (by norm_num) |>.mp h1)
  -- the real place: `n > 0`
  have hnpos : 0 < n := by
    rcases lt_or_gt_of_ne hn0 with hneg | hpos
    · exfalso
      have he2 : 0 < e ^ 2 := by positivity
      have h3 : 0 < (-n) * (12 * e ^ 2 - n) * (16 * e ^ 2 - n) := by
        apply mul_pos (mul_pos _ _) _ <;> linarith
      have h4 : (-n) * (12 * e ^ 2 - n) * (16 * e ^ 2 - n) = -(Z ^ 2) := by rw [hZ]; ring
      nlinarith [sq_nonneg Z]
    · exact hpos
  -- primes dividing `n` do not divide `e`
  have hpe : ∀ p : ℕ, p.Prime → (p : ℤ) ∣ n → ¬ (p : ℤ) ∣ e := by
    intro p hp hpn hpe'
    have h1 : (p : ℤ) ∣ 1 := (hne.isUnit_of_dvd' hpn hpe').dvd
    have h1' : p ∣ 1 := by exact_mod_cast h1
    exact hp.one_lt.ne' (Nat.dvd_one.mp h1')
  have hpe2 : ∀ p : ℕ, p.Prime → (p : ℤ) ∣ n → ∀ c : ℕ, (p : ℤ) ∣ (c : ℤ) * e ^ 2 →
      p ∣ c := by
    intro p hp hpn c hc
    have hpZ : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
    rcases hpZ.dvd_or_dvd hc with h | h
    · exact_mod_cast h
    · exact absurd (hpZ.dvd_of_dvd_pow h) (hpe p hp hpn)
  -- square classes of `n` and `n − 12e²`
  obtain ⟨d₁, a, hd₁, hna⟩ := split6 (G := (n - 12 * e ^ 2) * (n - 16 * e ^ 2)) hn0
    (by rw [hZ]; ring) (by
      intro p hp hpn hpG
      have hpZ : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
      rcases hpZ.dvd_or_dvd hpG with h | h
      · have : (p : ℤ) ∣ (12 : ℕ) * e ^ 2 := by
          have := dvd_sub hpn h; push_cast; convert this using 1; ring
        exact prime_23 hp (i := 2) (j := 1) (by simpa using hpe2 p hp hpn 12 this)
      · have : (p : ℤ) ∣ (16 : ℕ) * e ^ 2 := by
          have := dvd_sub hpn h; push_cast; convert this using 1; ring
        exact prime_23 hp (i := 4) (j := 0) (by simpa using hpe2 p hp hpn 16 this))
  obtain ⟨d₂, b, hd₂, hgb⟩ := split6 (G := n * (n - 16 * e ^ 2)) hg2
    (by rw [hZ]; ring) (by
      intro p hp hpg hpG
      have hpZ : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
      rcases hpZ.dvd_or_dvd hpG with h | h
      · have : (p : ℤ) ∣ (12 : ℕ) * e ^ 2 := by
          have := dvd_sub h hpg; push_cast; convert this using 1; ring
        exact prime_23 hp (i := 2) (j := 1) (by simpa using hpe2 p hp h 12 this)
      · have h4 : (p : ℤ) ∣ (4 : ℕ) * e ^ 2 := by
          have := dvd_sub hpg h; push_cast; convert this using 1; ring
        rcases hpZ.dvd_or_dvd h4 with h5 | h5
        · exact prime_23 hp (i := 2) (j := 0) (by simpa using (show p ∣ 4 by exact_mod_cast h5))
        · have hpe' := hpZ.dvd_of_dvd_pow h5
          have hpn : (p : ℤ) ∣ n := by
            have := dvd_add hpg (Dvd.dvd.mul_left (dvd_pow hpe' two_ne_zero) 12)
            convert this using 1; ring
          exact absurd hpe' (hpe p hp hpn))
  -- `d₁ > 0`
  have hd₁pos : 0 < d₁ := by
    by_contra hd
    have := mul_nonpos_of_nonpos_of_nonneg (not_lt.mp hd) (sq_nonneg a)
    linarith
  have hd₁' : d₁ ∈ ([1, 2, 3, 6] : List ℤ) := by
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hd₁ ⊢
    omega
  have ha0 : a ≠ 0 := by rintro rfl; apply hn0; rw [hna]; ring
  have hb0 : b ≠ 0 := by rintro rfl; apply hg2; rw [hgb]; ring
  -- the third class: `d₁d₂(n − 16e²)` is a square
  obtain ⟨t, ht⟩ : ∃ t : ℤ, t ^ 2 = d₁ * d₂ * (n - 16 * e ^ 2) := by
    apply int_sq_of_rat_sq (r := (Z : ℚ) / (a * b))
    have hint : Z ^ 2 = d₁ * d₂ * (n - 16 * e ^ 2) * (a * b) ^ 2 := by
      linear_combination hZ + (n - 16 * e ^ 2) * (n - 12 * e ^ 2) * hna +
        (n - 16 * e ^ 2) * d₁ * a ^ 2 * hgb
    have hab : ((a : ℚ) * b) ≠ 0 := by
      have : (a : ℚ) ≠ 0 := by exact_mod_cast ha0
      have : (b : ℚ) ≠ 0 := by exact_mod_cast hb0
      positivity
    rw [div_pow, div_eq_iff (pow_ne_zero 2 hab)]
    exact_mod_cast hint
  -- the certificate at `2`
  have hgood : (d₁, d₂) ∈ goodClasses := by
    rcases classes_cover d₁ hd₁' d₂ hd₂ with hg | hb
    · exact hg
    · exfalso
      have htrue := List.all_eq_true.mp badClasses_excluded _ hb
      have hfalse : excluded d₁ d₂ = false := by
        refine excluded_eq_false (a := a) (e := e) (b := b) (t := t) ?_ ?_ ?_
        · rw [← hna, hgb]
        · rw [← hna, ht]
        · rw [← hna]
          exact fun ⟨h1, h2⟩ => hpe 2 Nat.prime_two (by exact_mod_cast h1) (by exact_mod_cast h2)
      simp only at htrue
      rw [htrue] at hfalse
      exact Bool.noConfusion hfalse
  have hd₁0 : d₁ ≠ 0 := hd₁pos.ne'
  have hd₂0 : d₂ ≠ 0 := by rintro rfl; apply hg2; rw [hgb]; ring
  have ht0 : t ≠ 0 := by
    rintro rfl
    have : d₁ * d₂ * (n - 16 * e ^ 2) = 0 := by rw [← ht]; ring
    rcases mul_eq_zero.mp this with h1 | h1
    · rcases mul_eq_zero.mp h1 with h2 | h2
      · exact hd₁0 h2
      · exact hd₂0 h2
    · exact hg3 h1
  have heQ : (e : ℚ) ≠ 0 := by exact_mod_cast he0
  have hd12Q : ((d₁ * d₂ : ℤ) : ℚ) ≠ 0 := by exact_mod_cast mul_ne_zero hd₁0 hd₂0
  refine ⟨d₁, d₂, hgood, a / e, b / e, t / ((d₁ * d₂ : ℤ) * e), ?_, ?_, ?_, ?_, ?_, ?_⟩
  · exact div_ne_zero (by exact_mod_cast ha0) heQ
  · exact div_ne_zero (by exact_mod_cast hb0) heQ
  · exact div_ne_zero (by exact_mod_cast ht0) (mul_ne_zero hd12Q heQ)
  · have hnQ : (n : ℚ) = d₁ * a ^ 2 := by exact_mod_cast hna
    field_simp
    linear_combination hxe + hnQ
  · have hgQ : (n : ℚ) - 12 * e ^ 2 = d₂ * b ^ 2 := by exact_mod_cast hgb
    field_simp
    linear_combination hxe + hgQ
  · have htQ : (t : ℚ) ^ 2 = d₁ * d₂ * (n - 16 * e ^ 2) := by exact_mod_cast ht
    have hd1Q : (d₁ : ℚ) ≠ 0 := by exact_mod_cast hd₁0
    have hd2Q : (d₂ : ℚ) ≠ 0 := by exact_mod_cast hd₂0
    push_cast
    field_simp
    linear_combination (d₁ : ℚ) * d₂ * hxe - htQ

/-! ### 4. The curve `E'` and the map `Ψ` -/

/-- `E' : y² = x³ − 28x² + 192x = x(x − 12)(x − 16)`. -/
def E : Affine ℚ := ⟨0, -28, 0, 192, 0⟩

instance : E.IsCharNeTwoNF := ⟨rfl, rfl⟩

theorem Δ_E : E.Δ = 9437184 := by
  norm_num [E, WeierstrassCurve.Δ, WeierstrassCurve.b₂, WeierstrassCurve.b₄,
    WeierstrassCurve.b₆, WeierstrassCurve.b₈]

instance : E.IsElliptic := by
  rw [isElliptic_iff, isUnit_iff_ne_zero, Δ_E]
  norm_num

theorem f_eq : E.f = X * (X - C 12) * (X - C 16) := by
  simp only [Affine.f, E, map_neg, map_zero, map_ofNat]
  ring

theorem eval_f (x : ℚ) : E.f.eval x = x * (x - 12) * (x - 16) := by
  rw [f_eq]; simp

theorem equation_iff' (x y : ℚ) : E.Equation x y ↔ y ^ 2 = x * (x - 12) * (x - 16) := by
  rw [equation_iff_eval_f_eq_sq, eval_f, eq_comm]

theorem r0 : E.f.eval 0 = 0 := by rw [eval_f]; norm_num
theorem r12 : E.f.eval 12 = 0 := by rw [eval_f]; norm_num
theorem r16 : E.f.eval 16 = 0 := by rw [eval_f]; norm_num

/-- Evaluation at a root of `f`. -/
noncomputable def ev (r : ℚ) (hr : E.f.eval r = 0) : E.A →+* ℚ :=
  AdjoinRoot.lift (RingHom.id ℚ) r hr

theorem ev_mk (r : ℚ) (hr : E.f.eval r = 0) (p : ℚ[X]) :
    ev r hr (AdjoinRoot.mk E.f p) = p.eval r :=
  AdjoinRoot.lift_mk hr p

/-- `Ψ : Aˣ/(Aˣ)² → (ℚˣ/(ℚˣ)²)³`, evaluation at the roots `0, 12, 16`. -/
noncomputable def Ψ : E.M →* Units.modPow ℚ 2 × Units.modPow ℚ 2 × Units.modPow ℚ 2 :=
  (Units.modPow.map (ev 0 r0).toMonoidHom 2).prod
    ((Units.modPow.map (ev 12 r12).toMonoidHom 2).prod
      (Units.modPow.map (ev 16 r16).toMonoidHom 2))

theorem dvd_of_roots {p : ℚ[X]} (h0 : p.eval 0 = 0) (h12 : p.eval 12 = 0)
    (h16 : p.eval 16 = 0) : E.f ∣ p := by
  rw [f_eq]
  obtain ⟨p₁, rfl⟩ : X - C 0 ∣ p := dvd_iff_isRoot.mpr h0
  obtain ⟨p₂, rfl⟩ : X - C 12 ∣ p₁ := by
    apply dvd_iff_isRoot.mpr
    simp only [eval_mul, eval_sub, eval_X, eval_C] at h12
    simpa using h12
  obtain ⟨p₃, rfl⟩ : X - C 16 ∣ p₂ := by
    apply dvd_iff_isRoot.mpr
    simp only [eval_mul, eval_sub, eval_X, eval_C] at h16
    have : (16 - 0) * ((16 - 12) * eval 16 p₂) = 0 := h16
    have h' : (64 : ℚ) * eval 16 p₂ = 0 := by linear_combination this
    simpa using h'
  exact ⟨p₃, by simp only [map_zero, sub_zero]; ring⟩

/-- **`Ψ` is injective** (Lagrange interpolation at `0, 12, 16`). -/
theorem Ψ_injective : Function.Injective Ψ := by
  refine (injective_iff_map_eq_one Ψ).mpr fun m hm => ?_
  induction m using QuotientGroup.induction_on with
  | H u =>
    simp only [Ψ, MonoidHom.prod_apply, Units.modPow.map_mk, Prod.mk_eq_one] at hm
    obtain ⟨h0, h12, h16⟩ := hm
    rw [Units.modPow.mk_eq_one_iff_isSquare] at h0 h12 h16 ⊢
    simp only [Units.coe_map, RingHom.toMonoidHom_eq_coe, MonoidHom.coe_ofClass] at h0 h12 h16
    obtain ⟨g, hg⟩ := AdjoinRoot.mk_surjective (u : E.A)
    rw [← hg, ev_mk] at h0 h12 h16
    obtain ⟨s₀, hs₀⟩ := h0
    obtain ⟨s₁, hs₁⟩ := h12
    obtain ⟨s₂, hs₂⟩ := h16
    set L : ℚ[X] := C (s₀ / 192) * ((X - C 12) * (X - C 16)) +
      C (-s₁ / 48) * (X * (X - C 16)) + C (s₂ / 64) * (X * (X - C 12)) with hL
    have hL0 : L.eval 0 = s₀ := by simp [hL]; ring
    have hL12 : L.eval 12 = s₁ := by simp [hL]; ring
    have hL16 : L.eval 16 = s₂ := by simp [hL]; ring
    refine ⟨AdjoinRoot.mk E.f L, ?_⟩
    rw [← hg, ← map_mul, AdjoinRoot.mk_eq_mk]
    apply dvd_of_roots <;> simp [eval_mul, hL0, hL12, hL16, hs₀, hs₁, hs₂]

/-- The class of `d` in `ℚˣ/(ℚˣ)²` (and `1` for `d = 0`). -/
noncomputable def cls (d : ℚ) : Units.modPow ℚ 2 :=
  if hd : d = 0 then 1 else QuotientGroup.mk (Units.mk0 d hd)

theorem cls_eq {u : ℚˣ} {d w : ℚ} (hd : d ≠ 0) (hw : w ≠ 0) (h : (u : ℚ) = d * w ^ 2) :
    (QuotientGroup.mk u : Units.modPow ℚ 2) = cls d := by
  rw [cls, dif_neg hd, Units.modPow.mk_eq_mk_iff]
  refine ⟨(Units.mk0 w hw)⁻¹, Units.ext ?_⟩
  simp only [Units.val_mk0, Units.val_mul, Units.val_pow_eq_pow_val, Units.val_inv_eq_inv_val, h]
  field_simp

theorem cls_one : cls 1 = 1 := by
  rw [cls, dif_neg one_ne_zero]
  exact (Units.modPow.mk_eq_one_iff_isSquare).mpr ⟨1, by simp⟩

theorem comp_eq (r : ℚ) (hr : E.f.eval r = 0) {a : E.A} (ha : IsUnit a) {d w : ℚ} (hd : d ≠ 0)
    (hw : w ≠ 0) (h : ev r hr a = d * w ^ 2) :
    (QuotientGroup.mk (Units.map (ev r hr).toMonoidHom ha.unit) : Units.modPow ℚ 2) = cls d :=
  cls_eq hd hw (by simpa using h)

open Classical in
/-- The four triples of square classes. -/
noncomputable def T : Finset (Units.modPow ℚ 2 × Units.modPow ℚ 2 × Units.modPow ℚ 2) :=
  {(cls 1, cls 1, cls 1), (cls 3, cls (-3), cls (-9)), (cls 2, cls (-1), cls (-2)),
    (cls 6, cls 3, cls 18)}

open Classical in
theorem card_T : T.card ≤ 4 := Finset.card_le_four

theorem Ψ_unit {a : E.A} (ha : IsUnit a) :
    Ψ (ha.unit : E.M) = (QuotientGroup.mk (Units.map (ev 0 r0).toMonoidHom ha.unit),
      QuotientGroup.mk (Units.map (ev 12 r12).toMonoidHom ha.unit),
      QuotientGroup.mk (Units.map (ev 16 r16).toMonoidHom ha.unit)) := rfl

open Classical in
/-- **The image of every point lies in `T`.** -/
theorem mem_T (P : E.Point) : Ψ (μ (Multiplicative.ofAdd P)) ∈ T := by
  rw [μ_apply]
  rcases P with _ | ⟨x, y, h⟩
  · have : Ψ (μ₀ (W := E) 0) = (cls 1, cls 1, cls 1) := by
      rw [μ₀_zero, map_one, cls_one]; rfl
    rw [show (0 : E.Point) = .zero from rfl] at this
    rw [this]
    simp [T]
  · rw [μ₀_some]
    by_cases hx : E.f.eval x = 0
    · rw [μX_of_eval_f_eq_zero hx, Ψ_unit]
      rw [eval_f] at hx
      have hx' : x = 0 ∨ x = 12 ∨ x = 16 := by
        rcases mul_eq_zero.mp hx with h1 | h1
        · rcases mul_eq_zero.mp h1 with h2 | h2
          · exact Or.inl h2
          · exact Or.inr (Or.inl (by linarith))
        · exact Or.inr (Or.inr (by linarith))
      rcases hx' with rfl | rfl | rfl
      · rw [comp_eq 0 r0 _ (d := 3) (w := 8) (by norm_num) (by norm_num) (by
            rw [ev_mk]; simp [E]; norm_num),
          comp_eq 12 r12 _ (d := -3) (w := 2) (by norm_num) (by norm_num) (by
            rw [ev_mk]; simp [E]; norm_num),
          comp_eq 16 r16 _ (d := -9) (w := 4 / 3) (by norm_num) (by norm_num) (by
            rw [ev_mk]; simp [E]; norm_num)]
        simp [T]
      · rw [comp_eq 0 r0 _ (d := 3) (w := 2) (by norm_num) (by norm_num) (by
            rw [ev_mk]; simp [E]; norm_num),
          comp_eq 12 r12 _ (d := -3) (w := 4) (by norm_num) (by norm_num) (by
            rw [ev_mk]; simp [E]; norm_num),
          comp_eq 16 r16 _ (d := -9) (w := 2 / 3) (by norm_num) (by norm_num) (by
            rw [ev_mk]; simp [E]; norm_num)]
        simp [T]
      · rw [comp_eq 0 r0 _ (d := 1) (w := 4) (by norm_num) (by norm_num) (by
            rw [ev_mk]; simp [E]; norm_num),
          comp_eq 12 r12 _ (d := 1) (w := 2) (by norm_num) (by norm_num) (by
            rw [ev_mk]; simp [E]; norm_num),
          comp_eq 16 r16 _ (d := 1) (w := 8) (by norm_num) (by norm_num) (by
            rw [ev_mk]; simp [E]; norm_num)]
        simp [T]
    · rw [μX_of_eval_f_ne_zero hx, Ψ_unit]
      have hxe : y ^ 2 = x * (x - 12) * (x - 16) := (equation_iff' x y).mp h.1
      rw [eval_f] at hx
      have hx0 : x ≠ 0 := fun h0 => hx (by rw [h0]; ring)
      have hx12 : x - 12 ≠ 0 := fun h0 => hx (by rw [h0]; ring)
      have hx16 : x - 16 ≠ 0 := fun h0 => hx (by rw [h0]; ring)
      obtain ⟨d₁, d₂, hgood, w₁, w₂, w₃, hw₁, hw₂, hw₃, e1, e2, e3⟩ :=
        descent hxe hx0 hx12 hx16
      simp only [goodClasses, List.mem_cons, Prod.mk.injEq, List.not_mem_nil, or_false] at hgood
      rcases hgood with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ <;> norm_num at e1 e2 e3
      · rw [comp_eq 0 r0 _ (d := 1) (w := w₁) (by norm_num) hw₁ (by
            rw [ev_mk]; simp; linear_combination e1),
          comp_eq 12 r12 _ (d := 1) (w := w₂) (by norm_num) hw₂ (by
            rw [ev_mk]; simp; linear_combination e2),
          comp_eq 16 r16 _ (d := 1) (w := w₃) (by norm_num) hw₃ (by
            rw [ev_mk]; simp; linear_combination e3)]
        simp [T]
      · rw [comp_eq 0 r0 _ (d := 3) (w := w₁) (by norm_num) hw₁ (by
            rw [ev_mk]; simp; linear_combination e1),
          comp_eq 12 r12 _ (d := -3) (w := w₂) (by norm_num) hw₂ (by
            rw [ev_mk]; simp; linear_combination e2),
          comp_eq 16 r16 _ (d := -9) (w := w₃) (by norm_num) hw₃ (by
            rw [ev_mk]; simp; linear_combination e3)]
        simp [T]
      · rw [comp_eq 0 r0 _ (d := 2) (w := w₁) (by norm_num) hw₁ (by
            rw [ev_mk]; simp; linear_combination e1),
          comp_eq 12 r12 _ (d := -1) (w := w₂) (by norm_num) hw₂ (by
            rw [ev_mk]; simp; linear_combination e2),
          comp_eq 16 r16 _ (d := -2) (w := w₃) (by norm_num) hw₃ (by
            rw [ev_mk]; simp; linear_combination e3)]
        simp [T]
      · rw [comp_eq 0 r0 _ (d := 6) (w := w₁) (by norm_num) hw₁ (by
            rw [ev_mk]; simp; linear_combination e1),
          comp_eq 12 r12 _ (d := 3) (w := w₂) (by norm_num) hw₂ (by
            rw [ev_mk]; simp; linear_combination e2),
          comp_eq 16 r16 _ (d := 18) (w := w₃) (by norm_num) hw₃ (by
            rw [ev_mk]; simp; linear_combination e3)]
        simp [T]

/-! ### 5. Rank `0` -/

open Classical in
theorem card_range_μ_le : Nat.card (μ (W := E)).range ≤ 4 := by
  have hmem : ∀ m ∈ (μ (W := E)).range, Ψ m ∈ T := by
    rintro _ ⟨P, rfl⟩
    exact mem_T (Multiplicative.toAdd P)
  have hinj : Function.Injective
      (fun m : (μ (W := E)).range => (⟨Ψ m, hmem m m.2⟩ : T)) := by
    intro m m' hmm
    exact Subtype.ext (Ψ_injective (congrArg Subtype.val hmm))
  calc Nat.card (μ (W := E)).range ≤ Nat.card T := Nat.card_le_card_of_injective _ hinj
    _ = T.card := by rw [Nat.card_eq_fintype_card, Fintype.card_coe]
    _ ≤ 4 := card_T

theorem card_roots : Nat.card {x : ℚ // E.f.eval x = 0} = 3 := by
  have : {x : ℚ // E.f.eval x = 0} ≃ ({0, 12, 16} : Finset ℚ) := by
    refine Equiv.subtypeEquivRight fun x => ?_
    rw [eval_f]
    simp only [Finset.mem_insert, Finset.mem_singleton]
    constructor
    · intro hx
      rcases mul_eq_zero.mp hx with h1 | h1
      · rcases mul_eq_zero.mp h1 with h2 | h2
        · exact Or.inl h2
        · exact Or.inr (Or.inl (by linarith))
      · exact Or.inr (Or.inr (by linarith))
    · rintro (rfl | rfl | rfl) <;> norm_num
  rw [Nat.card_congr this, Nat.card_eq_fintype_card, Fintype.card_coe]
  rfl

/-- **`E'(ℚ)[2]` has `4` elements.** -/
theorem card_two_torsion : Nat.card (nsmulAddMonoidHom (α := E.Point) 2).ker = 4 := by
  rw [card_ker_nsmul_two, card_roots]

instance : IsAddFG E.Point := fg_point_of_numberField

/-- **`rank E'(ℚ) = 0`** for `E' : y² = x(x − 12)(x − 16)`: `2 ^ rank · 4 = #(im μ) ≤ 4`. -/
theorem finrank_eq_zero : Module.finrank ℤ E.Point = 0 := by
  have h := card_range_μ (W := E)
  rw [card_two_torsion] at h
  have h4 := card_range_μ_le
  rw [h] at h4
  by_contra hr
  have := Nat.le_self_pow hr 2
  omega

end LeanRank.RankZero

#print axioms LeanRank.RankZero.badClasses_excluded
#print axioms LeanRank.RankZero.excluded_eq_false
#print axioms LeanRank.RankZero.descent
#print axioms LeanRank.RankZero.Ψ_injective
#print axioms LeanRank.RankZero.mem_T
#print axioms LeanRank.RankZero.finrank_eq_zero
