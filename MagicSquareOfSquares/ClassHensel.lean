import MagicSquareOfSquares.ClassTable

/-!
# One power is enough (Proposition C of `COMPLETE.md`, section 5)

`ClassAt (p^K) (p³) d m n u v` at one power `p^K`, `K ≥ 18`, gives `ClassAt (p^k) (p³) d m n u v`
at **every** power `p^k`, at the same integer point `(u, v)`, for an odd prime `p` and `p² ∤ d`.
No `ℤ_p`: the proof is Hensel's lemma for integers (`exists_sq_modEq_pow`, `ClassLocal.lean`) and
for Gaussian integers modulo `p^k`, with valuations handled by explicit factorizations `p^e · c`.

* `sq_lift`: if `x ≡ y² (mod p^K)` and `p^j ∤ y` with `2(j − 1) < K`, then `x` is a square
  modulo every power of `p`, with a root not divisible by `p^j`.
* `classAt_norms_all`: **the norm half**, for `K ≥ 5`.
* `gsqM_lift_split` (`−1` a square modulo `p`) and `gsqM_lift_inert` (`−1` not a square modulo
  `p`): `d · Q_f` a Gaussian square modulo `p^K` gives the same modulo every `p^k`, when the norm
  of `Q_f` is not divisible by `p^17` and `p² ∤ d`; `K ≥ 18` and `K ≥ 10`.
* `classAt_all_of_classAt`: **the headline**, `K ≥ 18`. `passes_all_of_passes`: `Passes p 17 d m n`
  gives `Passes p k d m n` for every `k`.
* `exists_point_of_table`: every class of the table of Theorem C at an odd prime of `mn(m² − n²)`
  is the class modulo every power of `p` at **one** primitive integer point.
* Positive controls: `classAt_seven_all` and `classAt_seventeen_all`, from explicit points modulo
  `7^18` and `17^18` checked by the kernel (`scripts/class_hensel_witness.py` computes them).

The primitivity of `(u, v)` is not used. The hypothesis `p² ∤ d` is used only to bound the
valuation of `d · Q_f`; without some bound the statement is false (`d = p^19` passes trivially
modulo `p^18`). Statements about `ℤ_p` are still not in Lean: what is here is "modulo every power of
`p`, at one integer point".
-/

namespace MagicSquareOfSquares.CoverDescent

/-! ### Valuations, by factorization -/

/-- Every nonzero Gaussian integer `a + ib` is `p^r` times one that `p` does not divide. -/
theorem exists_content {p : ℤ} (hp : Prime p) :
    ∀ (N : ℕ) (a b : ℤ), a.natAbs + b.natAbs ≤ N → (a ≠ 0 ∨ b ≠ 0) →
      ∃ r : ℕ, ∃ a' b' : ℤ, a = p ^ r * a' ∧ b = p ^ r * b' ∧ ¬(p ∣ a' ∧ p ∣ b') := by
  intro N
  induction N with
  | zero =>
    intro a b hN hab
    have ha : a = 0 := Int.natAbs_eq_zero.mp (by omega)
    have hb : b = 0 := Int.natAbs_eq_zero.mp (by omega)
    rcases hab with h | h
    · exact absurd ha h
    · exact absurd hb h
  | succ N ih =>
    intro a b hN hab
    by_cases hd : p ∣ a ∧ p ∣ b
    · obtain ⟨⟨a₁, rfl⟩, ⟨b₁, rfl⟩⟩ := hd
      have hp2 : 2 ≤ p.natAbs := by
        have h0 : p.natAbs ≠ 0 := by
          rw [Ne, Int.natAbs_eq_zero]
          exact hp.ne_zero
        have h1 : p.natAbs ≠ 1 := fun h1 => hp.not_isUnit (Int.isUnit_iff_natAbs_eq.mpr h1)
        omega
      have hab₁ : a₁ ≠ 0 ∨ b₁ ≠ 0 := by
        rcases hab with h | h
        · left
          rintro rfl
          exact h (mul_zero p)
        · right
          rintro rfl
          exact h (mul_zero p)
      have hpos : 0 < a₁.natAbs + b₁.natAbs := by
        rcases hab₁ with h | h
        · have := Int.natAbs_pos.mpr h
          omega
        · have := Int.natAbs_pos.mpr h
          omega
      rw [Int.natAbs_mul, Int.natAbs_mul] at hN
      have e1 : 2 * a₁.natAbs ≤ p.natAbs * a₁.natAbs := Nat.mul_le_mul_right _ hp2
      have e2 : 2 * b₁.natAbs ≤ p.natAbs * b₁.natAbs := Nat.mul_le_mul_right _ hp2
      obtain ⟨r, a', b', ha, hb, hn⟩ := ih a₁ b₁ (by omega) hab₁
      exact ⟨r + 1, a', b', by rw [ha]; ring, by rw [hb]; ring, hn⟩
    · exact ⟨0, a, b, by ring, by ring, hd⟩

/-- Every nonzero integer is `p^e · c` with `p ∤ c`. -/
theorem exists_pow_mul_not_dvd {p : ℤ} (hp : Prime p) {y : ℤ} (hy : y ≠ 0) :
    ∃ e : ℕ, ∃ c : ℤ, y = p ^ e * c ∧ ¬p ∣ c := by
  obtain ⟨e, c, c', h1, h2, h3⟩ := exists_content hp _ y 0 le_rfl (Or.inl hy)
  have hc' : c' = 0 := by
    rcases mul_eq_zero.mp h2.symm with h | h
    · exact absurd h (pow_ne_zero _ hp.ne_zero)
    · exact h
  exact ⟨e, c, h1, fun hc => h3 ⟨hc, hc' ▸ dvd_zero p⟩⟩

/-- `p^j ∣ p^e · c` with `p ∤ c` forces `j ≤ e`. -/
theorem le_of_pow_dvd_pow_mul {p c : ℤ} (hp : Prime p) {j e : ℕ} (hc : ¬p ∣ c)
    (h : p ^ j ∣ p ^ e * c) : j ≤ e := by
  by_contra hje
  have h1 : p ^ (e + 1) ∣ p ^ e * c := dvd_trans (pow_dvd_pow p (by omega)) h
  rw [pow_succ] at h1
  exact hc ((mul_dvd_mul_iff_left (pow_ne_zero e hp.ne_zero)).mp h1)

theorem not_pow_dvd_of_lt {p c : ℤ} (hp : Prime p) {j e : ℕ} (hc : ¬p ∣ c) (hje : e < j) :
    ¬p ^ j ∣ p ^ e * c := fun h => absurd (le_of_pow_dvd_pow_mul hp hc h) (by omega)

/-- Valuations add: if `p^a ∤ x`, `p^b ∤ y` and `a + b ≤ c + 1`, then `p^c ∤ x y`. -/
theorem not_pow_dvd_mul {p x y : ℤ} (hp : Prime p) {a b c : ℕ} (hx : ¬p ^ a ∣ x)
    (hy : ¬p ^ b ∣ y) (habc : a + b ≤ c + 1) : ¬p ^ c ∣ x * y := by
  have hx0 : x ≠ 0 := by
    rintro rfl
    exact hx (dvd_zero _)
  have hy0 : y ≠ 0 := by
    rintro rfl
    exact hy (dvd_zero _)
  obtain ⟨e, c₁, rfl, hc₁⟩ := exists_pow_mul_not_dvd hp hx0
  obtain ⟨f, c₂, rfl, hc₂⟩ := exists_pow_mul_not_dvd hp hy0
  have he : e < a := by
    by_contra h
    exact hx (dvd_trans (pow_dvd_pow p (by omega)) (dvd_mul_right _ _))
  have hf : f < b := by
    by_contra h
    exact hy (dvd_trans (pow_dvd_pow p (by omega)) (dvd_mul_right _ _))
  have e1 : p ^ e * c₁ * (p ^ f * c₂) = p ^ (e + f) * (c₁ * c₂) := by ring
  rw [e1]
  exact not_pow_dvd_of_lt hp (not_dvd_mul_of_prime hp hc₁ hc₂) (by omega)

theorem dvd_of_modEq' {n a b : ℤ} (h : a ≡ b [ZMOD n]) (ha : n ∣ a) : n ∣ b :=
  Int.modEq_zero_iff_dvd.mp (h.symm.trans (Int.modEq_zero_iff_dvd.mpr ha))

theorem modEq_pow_iff_zmod {p k : ℕ} {a b : ℤ} :
    a ≡ b [ZMOD (p : ℤ) ^ k] ↔ (a : ZMod (p ^ k)) = (b : ZMod (p ^ k)) := by
  rw [ZMod.intCast_eq_intCast_iff, Nat.cast_pow]

/-! ### Squares of bounded valuation lift: the norm half -/

/-- **Square roots of bounded valuation lift** (`p` odd). If `x ≡ y² (mod p^K)` and `p^j ∤ y`
with `2(j − 1) < K`, then `x` is a square modulo every power of `p`, with a root not divisible by
`p^j`. Proof: `y = p^e c` with `p ∤ c` and `2e < K`, so `x = p^(2e) A` with `A ≡ c² (mod p)`, a
unit square modulo `p`; lift a root of `A` (`exists_sq_modEq_pow`) and multiply it by `p^e`. -/
theorem sq_lift {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) {K j : ℕ} (hjK : 2 * j ≤ K + 1) {x y : ℤ}
    (h : x ≡ y ^ 2 [ZMOD (p : ℤ) ^ K]) (hy : ¬(p : ℤ) ^ j ∣ y) (k : ℕ) :
    ∃ t : ℤ, x ≡ t ^ 2 [ZMOD (p : ℤ) ^ k] ∧ ¬(p : ℤ) ^ j ∣ t := by
  have hp' : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  have hy0 : y ≠ 0 := by
    rintro rfl
    exact hy (dvd_zero _)
  obtain ⟨e, c, rfl, hc⟩ := exists_pow_mul_not_dvd hp' hy0
  have hej : e < j := by
    by_contra hh
    exact hy (dvd_trans (pow_dvd_pow _ (by omega)) (dvd_mul_right _ _))
  obtain ⟨s, hs⟩ := h.symm.dvd
  have hK : (p : ℤ) ^ K = (p : ℤ) ^ (2 * e) * (p : ℤ) ^ (K - 2 * e) := by
    rw [← pow_add]
    congr 1
    omega
  have hxA : x = (p : ℤ) ^ (2 * e) * (c ^ 2 + (p : ℤ) ^ (K - 2 * e) * s) := by
    rw [hK] at hs
    linear_combination hs
  have hcA : c ^ 2 ≡ c ^ 2 + (p : ℤ) ^ (K - 2 * e) * s [ZMOD p] := by
    apply Int.modEq_iff_dvd.mpr
    rw [show c ^ 2 + (p : ℤ) ^ (K - 2 * e) * s - c ^ 2 = (p : ℤ) ^ (K - 2 * e) * s by ring]
    exact Dvd.dvd.mul_right (dvd_pow_self _ (by omega)) _
  obtain ⟨x', hx', hx'c⟩ := exists_sq_modEq_pow hp hp2 hc hcA k
  refine ⟨(p : ℤ) ^ e * x', ?_, not_pow_dvd_of_lt hp' (not_dvd_of_modEq hx'c hc) hej⟩
  have h1 := ((hx'.of_dvd (pow_dvd_pow _ (Nat.le_succ k))).mul_left ((p : ℤ) ^ (2 * e))).symm
  have e3 : ((p : ℤ) ^ e * x') ^ 2 = (p : ℤ) ^ (2 * e) * x' ^ 2 := by ring
  rw [hxA, e3]
  exact h1

/-- A norm with a root modulo `p^K` not divisible by `p³` (`K ≥ 5`) is not divisible by `p⁵`. -/
theorem norm_not_dvd {p : ℕ} (hp : p.Prime) {K : ℕ} (hK : 5 ≤ K) {N w : ℤ}
    (h : N ≡ w ^ 2 [ZMOD (p : ℤ) ^ K]) (hw : ¬(p : ℤ) ^ 3 ∣ w) : ¬(p : ℤ) ^ 5 ∣ N := by
  have hp' : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  intro hN
  have h5 : ¬(p : ℤ) ^ 5 ∣ w * w := not_pow_dvd_mul hp' hw hw (by norm_num)
  apply h5
  rw [← sq]
  exact dvd_of_modEq' (h.of_dvd (pow_dvd_pow _ hK)) hN

/-- **The norm half of Proposition C.** For an odd prime `p` and `K ≥ 5`, the roots in
`ClassAt (p^K) (p³)` give roots of the four norms modulo every power `p^k`, again not divisible by
`p³`, at the same point. (The valuation of a root is at most 2, so `v(N_j) = 2 v(w_j) ≤ 4 < K`,
and the unit part of `N_j` is a square modulo `p`.) -/
theorem classAt_norms_all {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) {K : ℕ} (hK : 5 ≤ K)
    {d m n u v : ℤ} (h : ClassAt ((p : ℤ) ^ K) ((p : ℤ) ^ 3) d m n u v) (k : ℕ) :
    ∃ w0 w8 wL wU : ℤ,
      (nrm (P0 m n u v) ≡ w0 ^ 2 [ZMOD (p : ℤ) ^ k] ∧ nrm (P8 m n u v) ≡ w8 ^ 2 [ZMOD (p : ℤ) ^ k] ∧
        nrm (PL m n u v) ≡ wL ^ 2 [ZMOD (p : ℤ) ^ k] ∧
          nrm (PU m n u v) ≡ wU ^ 2 [ZMOD (p : ℤ) ^ k]) ∧
      (¬(p : ℤ) ^ 3 ∣ w0 ∧ ¬(p : ℤ) ^ 3 ∣ w8 ∧ ¬(p : ℤ) ^ 3 ∣ wL ∧ ¬(p : ℤ) ^ 3 ∣ wU) := by
  obtain ⟨w0, w8, wL, wU, ⟨h0, h8, hL, hU⟩, ⟨q0, q8, qL, qU⟩, -⟩ := h
  obtain ⟨a0, e0, r0⟩ := sq_lift hp hp2 (j := 3) (by omega) h0 q0 k
  obtain ⟨a8, e8, r8⟩ := sq_lift hp hp2 (j := 3) (by omega) h8 q8 k
  obtain ⟨aL, eL, rL⟩ := sq_lift hp hp2 (j := 3) (by omega) hL qL k
  obtain ⟨aU, eU, rU⟩ := sq_lift hp hp2 (j := 3) (by omega) hU qU k
  exact ⟨a0, a8, aL, aU, ⟨e0, e8, eL, eU⟩, ⟨r0, r8, rL, rU⟩⟩

/-- At a `ClassAt (p^K) (p³)` point with `K ≥ 5`, the norm of each of the eight products,
`N_0 N_8 N_L N_U`, is not divisible by `p^17`. -/
theorem nrm_Q_not_dvd {p : ℕ} (hp : p.Prime) {K : ℕ} (hK : 5 ≤ K) {m n u v w0 w8 wL wU : ℤ}
    (h0 : nrm (P0 m n u v) ≡ w0 ^ 2 [ZMOD (p : ℤ) ^ K])
    (h8 : nrm (P8 m n u v) ≡ w8 ^ 2 [ZMOD (p : ℤ) ^ K])
    (hL : nrm (PL m n u v) ≡ wL ^ 2 [ZMOD (p : ℤ) ^ K])
    (hU : nrm (PU m n u v) ≡ wU ^ 2 [ZMOD (p : ℤ) ^ K]) (q0 : ¬(p : ℤ) ^ 3 ∣ w0)
    (q8 : ¬(p : ℤ) ^ 3 ∣ w8) (qL : ¬(p : ℤ) ^ 3 ∣ wL) (qU : ¬(p : ℤ) ^ 3 ∣ wU)
    (f : Bool × Bool × Bool) : ¬(p : ℤ) ^ 17 ∣ nrm (Q m n u v f) := by
  have hp' : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  rw [nrm_Q]
  have a0 := norm_not_dvd hp hK h0 q0
  have a8 := norm_not_dvd hp hK h8 q8
  have aL := norm_not_dvd hp hK hL qL
  have aU := norm_not_dvd hp hK hU qU
  have b1 : ¬(p : ℤ) ^ 9 ∣ nrm (P0 m n u v) * nrm (P8 m n u v) :=
    not_pow_dvd_mul hp' a0 a8 (by norm_num)
  have b2 : ¬(p : ℤ) ^ 13 ∣ nrm (P0 m n u v) * nrm (P8 m n u v) * nrm (PL m n u v) :=
    not_pow_dvd_mul hp' b1 aL (by norm_num)
  exact not_pow_dvd_mul hp' b2 aU (by norm_num)

/-! ### The Gaussian half when `−1` is a square modulo `p` -/

theorem ring_sq_split {R : Type*} [CommRing R] {ι X Y a b : R} (hι : ι ^ 2 = -1)
    (h1 : X = a ^ 2 - b ^ 2) (h2 : Y = 2 * a * b) : X + ι * Y = (a + ι * b) ^ 2 := by
  linear_combination h1 + ι * h2 - b ^ 2 * hι

theorem ring_norm_split {R : Type*} [CommRing R] {ι X Y : R} (hι : ι ^ 2 = -1) :
    (X + ι * Y) * (X + -ι * Y) = X ^ 2 + Y ^ 2 := by
  linear_combination (-Y ^ 2) * hι

theorem ring_two_roots {R : Type*} [CommRing R] {ι h X Y x₁ x₂ : R} (hι : ι ^ 2 = -1)
    (hh : 2 * h = 1) (h1 : X + ι * Y = x₁ ^ 2) (h2 : X + -ι * Y = x₂ ^ 2) :
    X = (h * (x₁ + x₂)) ^ 2 - (-ι * h * (x₁ - x₂)) ^ 2 ∧
      Y = 2 * (h * (x₁ + x₂)) * (-ι * h * (x₁ - x₂)) :=
  ⟨by linear_combination (-(2 * h + 1) * X) * hh + (2 * h ^ 2) * h1 + (2 * h ^ 2) * h2 +
      (h ^ 2 * (x₁ - x₂) ^ 2) * hι,
    by linear_combination (-(2 * h + 1) * Y) * hh + (4 * h ^ 2 * Y) * hι + (-2 * ι * h ^ 2) * h1 +
      (2 * ι * h ^ 2) * h2⟩

/-- If `ι² ≡ −1` and `X ± ιY` are squares modulo `p^k`, then `X + iY` is a Gaussian square
modulo `p^k` (`2` invertible). -/
theorem gsq_of_two_roots {p k : ℕ} {ι h X Y x₁ x₂ : ℤ} (hι : ι ^ 2 ≡ -1 [ZMOD (p : ℤ) ^ k])
    (hh : 2 * h ≡ 1 [ZMOD (p : ℤ) ^ k]) (h1 : X + ι * Y ≡ x₁ ^ 2 [ZMOD (p : ℤ) ^ k])
    (h2 : X + -ι * Y ≡ x₂ ^ 2 [ZMOD (p : ℤ) ^ k]) :
    ∃ a b : ℤ, X ≡ a ^ 2 - b ^ 2 [ZMOD (p : ℤ) ^ k] ∧ Y ≡ 2 * a * b [ZMOD (p : ℤ) ^ k] := by
  refine ⟨h * (x₁ + x₂), -ι * h * (x₁ - x₂), ?_, ?_⟩
  · rw [modEq_pow_iff_zmod] at hι hh h1 h2 ⊢
    push_cast at hι hh h1 h2 ⊢
    exact (ring_two_roots hι hh h1 h2).1
  · rw [modEq_pow_iff_zmod] at hι hh h1 h2 ⊢
    push_cast at hι hh h1 h2 ⊢
    exact (ring_two_roots hι hh h1 h2).2

/-- One of the two halves `X ± ιY` (split case): `d X + ι · d Y ≡ (a + ιb)²` has valuation at most
17, so it is a square modulo every power of `p`. -/
theorem split_half {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) {K : ℕ} (hK : 18 ≤ K)
    {d X Y a b ι : ℤ} (hd : ¬(p : ℤ) ^ 2 ∣ d) (hN : ¬(p : ℤ) ^ 17 ∣ X ^ 2 + Y ^ 2)
    (hι : ι ^ 2 ≡ -1 [ZMOD (p : ℤ) ^ K]) (h1 : d * X ≡ a ^ 2 - b ^ 2 [ZMOD (p : ℤ) ^ K])
    (h2 : d * Y ≡ 2 * a * b [ZMOD (p : ℤ) ^ K]) (k : ℕ) :
    ∃ x : ℤ, d * X + ι * (d * Y) ≡ x ^ 2 [ZMOD (p : ℤ) ^ k] := by
  have hp' : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  have hz : d * X + ι * (d * Y) ≡ (a + ι * b) ^ 2 [ZMOD (p : ℤ) ^ K] := by
    rw [modEq_pow_iff_zmod] at hι h1 h2 ⊢
    push_cast at hι h1 h2 ⊢
    exact ring_sq_split hι h1 h2
  have hq : (X + ι * Y) * (X + -ι * Y) ≡ X ^ 2 + Y ^ 2 [ZMOD (p : ℤ) ^ K] := by
    rw [modEq_pow_iff_zmod] at hι ⊢
    push_cast at hι ⊢
    exact ring_norm_split hι
  have h17 : ¬(p : ℤ) ^ 17 ∣ X + ι * Y := fun h =>
    hN (dvd_of_modEq' (hq.of_dvd (pow_dvd_pow _ (by omega))) (Dvd.dvd.mul_right h _))
  have h18 : ¬(p : ℤ) ^ 18 ∣ d * X + ι * (d * Y) := by
    rw [show d * X + ι * (d * Y) = d * (X + ι * Y) by ring]
    exact not_pow_dvd_mul hp' hd h17 (by norm_num)
  have h9 : ¬(p : ℤ) ^ 9 ∣ a + ι * b := by
    intro h
    apply h18
    have h' : (p : ℤ) ^ 18 ∣ (a + ι * b) ^ 2 := by
      have := pow_dvd_pow_of_dvd h 2
      rwa [← pow_mul] at this
    exact dvd_of_modEq' (hz.symm.of_dvd (pow_dvd_pow _ hK)) h'
  obtain ⟨t, ht, -⟩ := sq_lift hp hp2 (j := 9) (by omega) hz h9 k
  exact ⟨t, ht⟩

/-- **The Gaussian half, `−1` a square modulo `p`** (`p ≡ 1 (mod 4)`). With `ι² ≡ −1`,
`ℤ[i]/p^k ≅ ℤ/p^k × ℤ/p^k` by `x + iy ↦ (x + ιy, x − ιy)`. The two images of `d · Q_f` have
valuations adding up to that of `d² N(Q_f)`, so each is at most 17; each is a square modulo
`p^K`, `K ≥ 18`, hence modulo every power (`sq_lift`); and the two roots recombine. -/
theorem gsqM_lift_split {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2)
    (hs : ∃ c : ℤ, c ^ 2 ≡ -1 [ZMOD p]) {K : ℕ} (hK : 18 ≤ K) {d : ℤ} {q : ℤ × ℤ}
    (hd : ¬(p : ℤ) ^ 2 ∣ d) (hN : ¬(p : ℤ) ^ 17 ∣ nrm q) (h : GSqM ((p : ℤ) ^ K) d q) (k : ℕ) :
    GSqM ((p : ℤ) ^ k) d q := by
  have hp' : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  obtain ⟨c, hc⟩ := hs
  have hc0 : ¬(p : ℤ) ∣ c := by
    intro hpc
    have h1 : (p : ℤ) ∣ c ^ 2 := dvd_pow hpc two_ne_zero
    have h2 : (p : ℤ) ∣ -1 := dvd_of_modEq' hc h1
    exact hp'.not_isUnit (isUnit_of_dvd_one ((dvd_neg).mp h2))
  obtain ⟨ι, hι, -⟩ := exists_sq_modEq_pow hp hp2 hc0 hc (K + k)
  have hιK : ι ^ 2 ≡ -1 [ZMOD (p : ℤ) ^ K] := hι.of_dvd (pow_dvd_pow _ (by omega))
  have hιk : ι ^ 2 ≡ -1 [ZMOD (p : ℤ) ^ k] := hι.of_dvd (pow_dvd_pow _ (by omega))
  have hι' : (-ι) ^ 2 ≡ -1 [ZMOD (p : ℤ) ^ K] := by
    rw [neg_sq]
    exact hιK
  obtain ⟨X, Y⟩ := q
  change ¬(p : ℤ) ^ 17 ∣ X ^ 2 + Y ^ 2 at hN
  obtain ⟨a, b, h1, h2⟩ := h
  change d * X ≡ a ^ 2 - b ^ 2 [ZMOD (p : ℤ) ^ K] at h1
  change d * Y ≡ 2 * a * b [ZMOD (p : ℤ) ^ K] at h2
  show ∃ a b : ℤ, d * X ≡ a ^ 2 - b ^ 2 [ZMOD (p : ℤ) ^ k] ∧ d * Y ≡ 2 * a * b [ZMOD (p : ℤ) ^ k]
  obtain ⟨x₁, hx₁⟩ := split_half hp hp2 hK hd hN hιK h1 h2 k
  obtain ⟨x₂, hx₂⟩ := split_half hp hp2 hK hd hN hι' h1 h2 k
  have hcop : IsCoprime (2 : ℤ) ((p : ℤ) ^ k) :=
    ((hp'.irreducible.coprime_iff_not_dvd.mpr (not_dvd_two_of_odd hp hp2)).symm).pow_right
  obtain ⟨α, β, e⟩ := hcop
  have hh : 2 * α ≡ 1 [ZMOD (p : ℤ) ^ k] := Int.modEq_iff_dvd.mpr ⟨β, by linear_combination -e⟩
  exact gsq_of_two_roots hιk hh hx₁ hx₂

/-! ### The Gaussian half when `−1` is not a square modulo `p` -/

/-- If `−1` is not a square modulo `p` and `p ∤ x + iy`, then `p ∤ x² + y²`. -/
theorem not_dvd_nrm_of_inert {p : ℕ} (hp : p.Prime) (hs : ¬∃ c : ℤ, c ^ 2 ≡ -1 [ZMOD p])
    {x y : ℤ} (h : ¬((p : ℤ) ∣ x ∧ (p : ℤ) ∣ y)) : ¬(p : ℤ) ∣ x ^ 2 + y ^ 2 := by
  have hp' : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  rintro ⟨w, hw⟩
  by_cases hy : (p : ℤ) ∣ y
  · apply h
    refine ⟨hp'.dvd_of_dvd_pow (n := 2) ?_, hy⟩
    have e : x ^ 2 = (p : ℤ) * w - y ^ 2 := by linear_combination hw
    rw [e]
    exact dvd_sub (dvd_mul_right _ _) (dvd_pow hy two_ne_zero)
  · obtain ⟨η, β, hab⟩ : IsCoprime y (p : ℤ) :=
      (hp'.irreducible.coprime_iff_not_dvd.mpr hy).symm
    apply hs
    refine ⟨x * η, Int.modEq_iff_dvd.mpr ⟨-(η ^ 2 * w) - 2 * β + β ^ 2 * p, ?_⟩⟩
    linear_combination (-η ^ 2) * hw + (η * y + 1 - β * p) * hab

/-- **Hensel's lemma for squares in `ℤ[i]` modulo `p^k`, `−1` not a square modulo `p`.** A
Gaussian integer `A + iB` that is congruent modulo `p` to the square of some `c₁ + ic₂` not
divisible by `p` is a Gaussian square modulo every power of `p`. The Newton step: if
`A + iB − (x + iy)² = p^(k+1) (t₁ + it₂)`, add `p^(k+1) α (t₁ + it₂)(x − iy)`, where
`α · 2(x² + y²) ≡ 1 (mod p)`. -/
theorem gauss_hensel {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) (hs : ¬∃ c : ℤ, c ^ 2 ≡ -1 [ZMOD p])
    {A B c₁ c₂ : ℤ} (hc : ¬((p : ℤ) ∣ c₁ ∧ (p : ℤ) ∣ c₂)) (h1 : c₁ ^ 2 - c₂ ^ 2 ≡ A [ZMOD p])
    (h2 : 2 * c₁ * c₂ ≡ B [ZMOD p]) (k : ℕ) :
    ∃ x y : ℤ, x ^ 2 - y ^ 2 ≡ A [ZMOD (p : ℤ) ^ (k + 1)] ∧
      2 * x * y ≡ B [ZMOD (p : ℤ) ^ (k + 1)] ∧ x ≡ c₁ [ZMOD p] ∧ y ≡ c₂ [ZMOD p] := by
  have hp' : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  have h2p := not_dvd_two_of_odd hp hp2
  induction k with
  | zero =>
    refine ⟨c₁, c₂, ?_, ?_, Int.ModEq.refl _, Int.ModEq.refl _⟩
    · rw [zero_add, pow_one]
      exact h1
    · rw [zero_add, pow_one]
      exact h2
  | succ k ih =>
    obtain ⟨x, y, hx, hy, hxc, hyc⟩ := ih
    have hxy : ¬((p : ℤ) ∣ x ∧ (p : ℤ) ∣ y) := fun ⟨a, b⟩ =>
      hc ⟨dvd_of_modEq' hxc a, dvd_of_modEq' hyc b⟩
    have hN := not_dvd_nrm_of_inert hp hs hxy
    obtain ⟨α, β, hab⟩ : IsCoprime (2 * (x ^ 2 + y ^ 2)) (p : ℤ) :=
      (hp'.irreducible.coprime_iff_not_dvd.mpr (not_dvd_mul_of_prime hp' h2p hN)).symm
    obtain ⟨t₁, ht₁⟩ := hx.dvd
    obtain ⟨t₂, ht₂⟩ := hy.dvd
    have h0 : ∀ z : ℤ, (p : ℤ) ^ (k + 1) * z ≡ 0 [ZMOD p] := fun z =>
      Int.modEq_zero_iff_dvd.mpr (Dvd.dvd.mul_right (dvd_pow_self _ (Nat.succ_ne_zero k)) _)
    refine ⟨x + (p : ℤ) ^ (k + 1) * (α * (t₁ * x + t₂ * y)),
      y + (p : ℤ) ^ (k + 1) * (α * (t₂ * x - t₁ * y)), ?_, ?_, ?_, ?_⟩
    · apply Int.modEq_iff_dvd.mpr
      exact ⟨β * t₁ - (p : ℤ) ^ k *
          ((α * (t₁ * x + t₂ * y)) ^ 2 - (α * (t₂ * x - t₁ * y)) ^ 2), by
        linear_combination ht₁ - (p : ℤ) ^ (k + 1) * t₁ * hab⟩
    · apply Int.modEq_iff_dvd.mpr
      exact ⟨β * t₂ - (p : ℤ) ^ k * (2 * (α * (t₁ * x + t₂ * y)) * (α * (t₂ * x - t₁ * y))), by
        linear_combination ht₂ - (p : ℤ) ^ (k + 1) * t₂ * hab⟩
    · have := hxc.add (h0 (α * (t₁ * x + t₂ * y)))
      rwa [add_zero] at this
    · have := hyc.add (h0 (α * (t₂ * x - t₁ * y)))
      rwa [add_zero] at this

/-- **The Gaussian half, `−1` not a square modulo `p`** (`p ≡ 3 (mod 4)`, `K ≥ 10`). A root
`a + ib` of `d · Q_f` modulo `p^K` is `p^r (a' + ib')` with `p ∤ a' + ib'`; `2r < K` because
`d · Q_f` is not `0` modulo `p^K` (its norm `d² N(Q_f)` is not divisible by `p^19`). Then
`d · Q_f = p^(2r) A` with `A ≡ (a' + ib')² (mod p)`, and `gauss_hensel` lifts a root of `A`. -/
theorem gsqM_lift_inert {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2)
    (hs : ¬∃ c : ℤ, c ^ 2 ≡ -1 [ZMOD p]) {K : ℕ} (hK : 10 ≤ K) {d : ℤ} {q : ℤ × ℤ}
    (hd : ¬(p : ℤ) ^ 2 ∣ d) (hN : ¬(p : ℤ) ^ 17 ∣ nrm q) (h : GSqM ((p : ℤ) ^ K) d q) (k : ℕ) :
    GSqM ((p : ℤ) ^ k) d q := by
  have hp' : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  obtain ⟨X, Y⟩ := q
  change ¬(p : ℤ) ^ 17 ∣ X ^ 2 + Y ^ 2 at hN
  obtain ⟨a, b, h1, h2⟩ := h
  change d * X ≡ a ^ 2 - b ^ 2 [ZMOD (p : ℤ) ^ K] at h1
  change d * Y ≡ 2 * a * b [ZMOD (p : ℤ) ^ K] at h2
  show ∃ a b : ℤ, d * X ≡ a ^ 2 - b ^ 2 [ZMOD (p : ℤ) ^ k] ∧ d * Y ≡ 2 * a * b [ZMOD (p : ℤ) ^ k]
  have hz : ¬((p : ℤ) ^ K ∣ d * X ∧ (p : ℤ) ^ K ∣ d * Y) := by
    rintro ⟨⟨s, hs⟩, ⟨t, ht⟩⟩
    have h3 : ¬(p : ℤ) ^ 3 ∣ d * d := not_pow_dvd_mul hp' hd hd (by norm_num)
    have h19 : ¬(p : ℤ) ^ 19 ∣ d * d * (X ^ 2 + Y ^ 2) := not_pow_dvd_mul hp' h3 hN (by norm_num)
    apply h19
    have hpp : (p : ℤ) ^ 19 * (p : ℤ) ^ (2 * K - 19) = (p : ℤ) ^ K * (p : ℤ) ^ K := by
      rw [← pow_add, ← pow_add]
      congr 1
      omega
    have e : d * d * (X ^ 2 + Y ^ 2) =
        (p : ℤ) ^ 19 * ((p : ℤ) ^ (2 * K - 19) * (s ^ 2 + t ^ 2)) := by
      linear_combination (d * X + (p : ℤ) ^ K * s) * hs + (d * Y + (p : ℤ) ^ K * t) * ht -
        (s ^ 2 + t ^ 2) * hpp
    rw [e]
    exact dvd_mul_right _ _
  have hab : a ≠ 0 ∨ b ≠ 0 := by
    by_cases ha : a = 0
    · by_cases hb : b = 0
      · subst ha
        subst hb
        exact absurd ⟨dvd_of_modEq' h1.symm (by simp), dvd_of_modEq' h2.symm (by simp)⟩ hz
      · exact Or.inr hb
    · exact Or.inl ha
  obtain ⟨r, a', b', rfl, rfl, hc⟩ := exists_content hp' _ a b le_rfl hab
  have hr : 2 * r < K := by
    by_contra hr
    have hdv : (p : ℤ) ^ K ∣ (p : ℤ) ^ r * (p : ℤ) ^ r := by
      rw [← pow_add]
      exact pow_dvd_pow _ (by omega)
    exact hz ⟨dvd_of_modEq' h1.symm (dvd_trans hdv ⟨a' ^ 2 - b' ^ 2, by ring⟩),
      dvd_of_modEq' h2.symm (dvd_trans hdv ⟨2 * a' * b', by ring⟩)⟩
  obtain ⟨s₁, hs₁⟩ := h1.symm.dvd
  obtain ⟨s₂, hs₂⟩ := h2.symm.dvd
  have hKr : (p : ℤ) ^ K = (p : ℤ) ^ (2 * r) * (p : ℤ) ^ (K - 2 * r) := by
    rw [← pow_add]
    congr 1
    omega
  rw [hKr] at hs₁ hs₂
  have hpd : (p : ℤ) ∣ (p : ℤ) ^ (K - 2 * r) := dvd_pow_self _ (by omega)
  have c1 : a' ^ 2 - b' ^ 2 ≡ a' ^ 2 - b' ^ 2 + (p : ℤ) ^ (K - 2 * r) * s₁ [ZMOD p] := by
    apply Int.modEq_iff_dvd.mpr
    rw [show a' ^ 2 - b' ^ 2 + (p : ℤ) ^ (K - 2 * r) * s₁ - (a' ^ 2 - b' ^ 2) =
      (p : ℤ) ^ (K - 2 * r) * s₁ by ring]
    exact Dvd.dvd.mul_right hpd _
  have c2 : 2 * a' * b' ≡ 2 * a' * b' + (p : ℤ) ^ (K - 2 * r) * s₂ [ZMOD p] := by
    apply Int.modEq_iff_dvd.mpr
    rw [show 2 * a' * b' + (p : ℤ) ^ (K - 2 * r) * s₂ - 2 * a' * b' =
      (p : ℤ) ^ (K - 2 * r) * s₂ by ring]
    exact Dvd.dvd.mul_right hpd _
  obtain ⟨x, y, hx, hy, -, -⟩ := gauss_hensel hp hp2 hs hc c1 c2 k
  have hdk : (p : ℤ) ^ k ∣ (p : ℤ) ^ (k + 1) := pow_dvd_pow _ (Nat.le_succ k)
  refine ⟨(p : ℤ) ^ r * x, (p : ℤ) ^ r * y, ?_, ?_⟩
  · have h' := ((hx.of_dvd hdk).mul_left ((p : ℤ) ^ (2 * r))).symm
    have e1 : d * X = (p : ℤ) ^ (2 * r) * (a' ^ 2 - b' ^ 2 + (p : ℤ) ^ (K - 2 * r) * s₁) := by
      linear_combination hs₁
    have e2 : ((p : ℤ) ^ r * x) ^ 2 - ((p : ℤ) ^ r * y) ^ 2 =
        (p : ℤ) ^ (2 * r) * (x ^ 2 - y ^ 2) := by ring
    rw [e1, e2]
    exact h'
  · have h' := ((hy.of_dvd hdk).mul_left ((p : ℤ) ^ (2 * r))).symm
    have e1 : d * Y = (p : ℤ) ^ (2 * r) * (2 * a' * b' + (p : ℤ) ^ (K - 2 * r) * s₂) := by
      linear_combination hs₂
    have e2 : 2 * ((p : ℤ) ^ r * x) * ((p : ℤ) ^ r * y) = (p : ℤ) ^ (2 * r) * (2 * x * y) := by
      ring
    rw [e1, e2]
    exact h'

/-! ### Proposition C -/

/-- **The Gaussian half of Proposition C.** For an odd prime `p`, `K ≥ 18`, `p² ∤ d` and a
Gaussian integer `q` whose norm is not divisible by `p^17`: if `d · q` is a Gaussian square modulo
`p^K`, it is one modulo every power of `p`. -/
theorem gsqM_lift {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) {K : ℕ} (hK : 18 ≤ K) {d : ℤ}
    {q : ℤ × ℤ} (hd : ¬(p : ℤ) ^ 2 ∣ d) (hN : ¬(p : ℤ) ^ 17 ∣ nrm q)
    (h : GSqM ((p : ℤ) ^ K) d q) (k : ℕ) : GSqM ((p : ℤ) ^ k) d q := by
  by_cases hs : ∃ c : ℤ, c ^ 2 ≡ -1 [ZMOD p]
  · exact gsqM_lift_split hp hp2 hs hK hd hN h k
  · exact gsqM_lift_inert hp hp2 hs (by omega) hd hN h k

/-- **Proposition C of `COMPLETE.md` (one power is enough).** Let `p` be an odd prime, `p² ∤ d`,
and `K ≥ 18`. If `d` is the class of the integer point `(u, v)` modulo `p^K`
(`ClassAt (p^K) (p³) d m n u v`: roots of the four norms not divisible by `p³`, and `d` times each
of the eight products a Gaussian square), then the same holds modulo **every** power of `p`, with
the same bound `p³` on the roots, at the same point. The norm half needs only `K ≥ 5`
(`classAt_norms_all`); the eight products need `K ≥ 18` when `−1` is a square modulo `p` and
`K ≥ 10` when it is not. -/
theorem classAt_all_of_classAt {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) {K : ℕ} (hK : 18 ≤ K)
    {d m n u v : ℤ} (hd : ¬(p : ℤ) ^ 2 ∣ d) (h : ClassAt ((p : ℤ) ^ K) ((p : ℤ) ^ 3) d m n u v)
    (k : ℕ) : ClassAt ((p : ℤ) ^ k) ((p : ℤ) ^ 3) d m n u v := by
  obtain ⟨w0, w8, wL, wU, ⟨h0, h8, hL, hU⟩, ⟨q0, q8, qL, qU⟩, g⟩ := h
  obtain ⟨a0, e0, r0⟩ := sq_lift hp hp2 (j := 3) (by omega) h0 q0 k
  obtain ⟨a8, e8, r8⟩ := sq_lift hp hp2 (j := 3) (by omega) h8 q8 k
  obtain ⟨aL, eL, rL⟩ := sq_lift hp hp2 (j := 3) (by omega) hL qL k
  obtain ⟨aU, eU, rU⟩ := sq_lift hp hp2 (j := 3) (by omega) hU qU k
  refine ⟨a0, a8, aL, aU, ⟨e0, e8, eL, eU⟩, ⟨r0, r8, rL, rU⟩, fun f => ?_⟩
  exact gsqM_lift hp hp2 hK hd (nrm_Q_not_dvd hp (by omega) h0 h8 hL hU q0 q8 qL qU f) (g f) k

/-- **`Passes` at the single power `p^18` gives `Passes` at every power** (odd `p`, `p² ∤ d`), at
the same primitive point. -/
theorem passes_all_of_passes {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) {d m n : ℤ}
    (hd : ¬(p : ℤ) ^ 2 ∣ d) (h : Passes p 17 d m n) (k : ℕ) : Passes p k d m n := by
  obtain ⟨u, v, huv, hc⟩ := h
  exact ⟨u, v, huv, classAt_all_of_classAt hp hp2 (by norm_num) hd hc (k + 1)⟩

/-- A class of the table is not divisible by `p²`. -/
theorem not_sq_dvd_of_row {p : ℕ} (hp : p.Prime) {x d : ℤ} (h : Row p x d) :
    ¬(p : ℤ) ^ 2 ∣ d := by
  have hp' : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  rcases h with ⟨hd, -⟩ | ⟨d1, rfl, hd1, -⟩
  · exact fun h2 => hd (dvd_trans (dvd_pow_self _ two_ne_zero) h2)
  · intro h2
    rw [sq] at h2
    exact hd1 ((mul_dvd_mul_iff_left hp'.ne_zero).mp h2)

/-- **The lower bounds of Theorem C at one integer point.** For coprime `m`, `n`, an odd prime `p`
dividing one of `m`, `n`, `m + n`, `m − n`, and a class `d` of the row of that factor, there is
**one** integer point `(u, v)`, primitive modulo `p`, at which `d` is the class modulo every power
of `p`: `ClassAt (p^k) (p³) d m n u v` for all `k`. (`passes_of_table` gives, for each `k`, a point
that may depend on `k`; here the point of `k = 17` serves for every `k`, by
`classAt_all_of_classAt`.) -/
theorem exists_point_of_table {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) {m n d : ℤ}
    (hc : IsCoprime m n)
    (h : ((p : ℤ) ∣ m ∧ Row p m d) ∨ ((p : ℤ) ∣ n ∧ Row p n d) ∨
      ((p : ℤ) ∣ m + n ∧ Row p (m + n) d) ∨ ((p : ℤ) ∣ m - n ∧ Row p (m - n) d)) :
    ∃ u v : ℤ, (¬(p : ℤ) ∣ u ∨ ¬(p : ℤ) ∣ v) ∧
      ∀ k : ℕ, ClassAt ((p : ℤ) ^ k) ((p : ℤ) ^ 3) d m n u v := by
  have hd : ¬(p : ℤ) ^ 2 ∣ d := by
    rcases h with ⟨-, r⟩ | ⟨-, r⟩ | ⟨-, r⟩ | ⟨-, r⟩ <;> exact not_sq_dvd_of_row hp r
  obtain ⟨u, v, huv, hcl⟩ := passes_of_table hp hp2 hc h 17
  exact ⟨u, v, huv, fun k => classAt_all_of_classAt hp hp2 (by norm_num) hd hcl k⟩

/-! ### Positive controls

Two points modulo `p^18`, computed by `scripts/class_hensel_witness.py` and checked here by the
kernel, then lifted to every power by `classAt_all_of_classAt`. Both have `p ∣ d`, `p ∣ w_0`
(`N_0 = p² · 50` at 7, `p² · 290` at 17) and Gaussian content `p²` in every `d · Q_f`. -/

/-- `21` is the class of the point `(7, 1)` for the ratio `49 : 1` modulo `7^18` (`−1` is not a
square modulo 7). -/
theorem classAt_seven : ClassAt (((7 : ℕ) : ℤ) ^ 18) (((7 : ℕ) : ℤ) ^ 3) 21 49 1 7 1 := by
  refine ⟨1499564513154295, 1017756768591956,
    344394874057069, 206312675315127,
    by decide +kernel, by decide +kernel, ?_⟩
  rintro ⟨_ | _, _ | _, _ | _⟩
  · exact ⟨1407323163836282, 893882059475846,
      by decide +kernel, by decide +kernel⟩
  · exact ⟨910804828866527, 1398389747125448,
      by decide +kernel, by decide +kernel⟩
  · exact ⟨1399364424181556, 120436846507579,
      by decide +kernel, by decide +kernel⟩
  · exact ⟨1592663569525917, 1577849882151900,
      by decide +kernel, by decide +kernel⟩
  · exact ⟨1017249538956959, 1421836491000401,
      by decide +kernel, by decide +kernel⟩
  · exact ⟨998519696562934, 803891191955495,
      by decide +kernel, by decide +kernel⟩
  · exact ⟨1043173622143499, 1331215016735646,
      by decide +kernel, by decide +kernel⟩
  · exact ⟨553315069768713, 497842015910372,
      by decide +kernel, by decide +kernel⟩

/-- `34` is the class of the point `(17, 1)` for the ratio `289 : 1` modulo `17^18` (`−1` is a
square modulo 17). -/
theorem classAt_seventeen :
    ClassAt (((17 : ℕ) : ℤ) ^ 18) (((17 : ℕ) : ℤ) ^ 3) 34 289 1 17 1 := by
  refine ⟨9539564116624916567129, 8789427709714560534446,
    6837009740152637118450, 10983903513319529333223,
    by decide +kernel, by decide +kernel, ?_⟩
  rintro ⟨_ | _, _ | _, _ | _⟩
  · exact ⟨12602329142185559595963, 3906724355988192813374,
      by decide +kernel, by decide +kernel⟩
  · exact ⟨9003377763113330529692, 328711152694744620548,
      by decide +kernel, by decide +kernel⟩
  · exact ⟨7511789223806021553812, 11532397477959865351184,
      by decide +kernel, by decide +kernel⟩
  · exact ⟨6257406387504173158533, 1255545820459354662515,
      by decide +kernel, by decide +kernel⟩
  · exact ⟨12843151382803680433901, 9750330412946983507333,
      by decide +kernel, by decide +kernel⟩
  · exact ⟨3979113668546697350008, 1849728055357526439775,
      by decide +kernel, by decide +kernel⟩
  · exact ⟨185238514372281680285, 3626927607313203802973,
      by decide +kernel, by decide +kernel⟩
  · exact ⟨12305944333319741352002, 773139297583296184670,
      by decide +kernel, by decide +kernel⟩

/-- Positive control, `−1` not a square modulo `p`: the class `21` at `(7, 1)` for `49 : 1`, modulo
every power of 7. -/
theorem classAt_seven_all (k : ℕ) :
    ClassAt (((7 : ℕ) : ℤ) ^ k) (((7 : ℕ) : ℤ) ^ 3) 21 49 1 7 1 :=
  classAt_all_of_classAt (by norm_num) (by norm_num) le_rfl (by decide +kernel) classAt_seven k

/-- Positive control, `−1` a square modulo `p`: the class `34` at `(17, 1)` for `289 : 1`, modulo
every power of 17. -/
theorem classAt_seventeen_all (k : ℕ) :
    ClassAt (((17 : ℕ) : ℤ) ^ k) (((17 : ℕ) : ℤ) ^ 3) 34 289 1 17 1 :=
  classAt_all_of_classAt (by norm_num) (by norm_num) le_rfl (by decide +kernel) classAt_seventeen k

end MagicSquareOfSquares.CoverDescent

