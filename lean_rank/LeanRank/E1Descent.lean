import LeanRank.Zsqrt2
import LeanRank.E1RankGlue

/-!
# The two descent statements for `E1` over `ℚ(√2)`, hence `rank E1(ℚ(√2)) = 0`

Written by an AI model (Claude, made by Anthropic) on 26 November 2026 (UTC); not read by a person.
Experiment, not counted in the census.

`E1 : y² = x(x² + Ax + B)` with `A = √2⁷ε = 8(2 + √2)`, `B = −√2⁸ε² = −8(2 + √2)²`, and the isogenous
`E1' : Y² = X(X² + A'X + B')`, `A' = −√2⁹ε`, `B' = 3·√2¹²ε²`, over `ℤ√2` (`ε = 1 + √2`).
For `x = δw²` (`δ ∈ ℤ√2` squarefree) and `w = s/t` the equation gives `n² = r²·F(s, t)`,
`F = δ(δ²s⁴ + Aδs²t² + Bt⁴)`, over `ℤ√2` (`quartic_of_point`). Then:
* `prime_dvd_B`: a prime factor of a squarefree `δ` divides `B` (odd valuation, `not_sq_odd`);
* so `δ = unit · √2^i · 3^j` (`classify1`, `classify2`), and units are `±v²`, `±εv²` (`unit_classes`);
* local exclusions: odd `√2`-valuation (`V1`, `V2`), a mod-2 test for `δ = ±ε` on `E1` (`M1`,
  `not_sq_oddim`), the two real places on `E1'` (`real_excl`).
`descentE1`, `descentE1'` are the hypotheses `DescentW`, `DescentW'` of `E1RankGlue.lean`, and
**`E1K_finrank_eq_zero : Module.finrank ℤ E1K.Point = 0`**.
-/

namespace MagicSquareOfSquares.E1Descent

open Zsqrtd MagicSquareOfSquares.Zsqrt2 MagicSquareOfSquares.Family MagicSquareOfSquares.E1Rank

/-- `F(s, t) = c(c²s⁴ + Acs²t² + Bt⁴)`. -/
def Fq (A B c s t : ℤ√2) : ℤ√2 := c * (c ^ 2 * s ^ 4 + A * c * s ^ 2 * t ^ 2 + B * t ^ 4)

def A1 : ℤ√2 := pi2 ^ 7 * eps
def B1 : ℤ√2 := -(pi2 ^ 8 * eps ^ 2)
def A2 : ℤ√2 := -(pi2 ^ 9 * eps)
def B2 : ℤ√2 := three * pi2 ^ 12 * eps ^ 2

/-! ### Divisibility helpers -/

theorem norm_lt_of_mul {p r₁ : ℤ√2} (hp0 : p ≠ 0) (hpu : ¬IsUnit p) (hr₁ : r₁ ≠ 0) :
    r₁.norm.natAbs < (p * r₁).norm.natAbs := by
  rw [norm_mul, Int.natAbs_mul]
  have h1 := natAbs_norm_ge_two hp0 hpu
  have h2 := Int.natAbs_pos.mpr (norm_ne_zero hr₁)
  nlinarith

theorem nd_add {p X Y : ℤ√2} (hX : ¬p ∣ X) : ¬p ∣ X + p * Y :=
  fun h => hX ((dvd_add_left (dvd_mul_right _ _)).mp h)

theorem nd_mul {p X Y : ℤ√2} (hp : Prime p) (hX : ¬p ∣ X) (hY : ¬p ∣ Y) : ¬p ∣ X * Y :=
  fun h => (hp.dvd_or_dvd h).elim hX hY

theorem nd_pow {p X : ℤ√2} (hp : Prime p) (hX : ¬p ∣ X) (k : ℕ) : ¬p ∣ X ^ k :=
  fun h => hX (hp.dvd_of_dvd_pow h)

theorem nd_neg {p X : ℤ√2} (hX : ¬p ∣ X) : ¬p ∣ -X := fun h => hX ((dvd_neg).mp h)

theorem nd_unit {p u : ℤ√2} (hp : Prime p) (hu : IsUnit u) : ¬p ∣ u :=
  fun h => hp.not_unit (isUnit_of_dvd_unit h hu)

theorem pi2_nd_eps : ¬pi2 ∣ eps := by rw [pi2_dvd_iff]; decide

theorem pi2_nd_three : ¬pi2 ∣ three := by rw [pi2_dvd_iff]; decide

theorem pi2_dvd_or_sub_one (z : ℤ√2) : pi2 ∣ z ∨ pi2 ∣ z - 1 := by
  rw [pi2_dvd_iff, pi2_dvd_iff, re_sub, re_one]; omega

/-- A unit modulo `√2` has `s⁴ ≡ 1 (mod √2⁵)`. -/
theorem pow4_of_not_dvd {s : ℤ√2} (hs : ¬pi2 ∣ s) : ∃ k, s ^ 4 = 1 + pi2 ^ 5 * k := by
  obtain ⟨c, hc⟩ := (pi2_dvd_or_sub_one s).resolve_left hs
  have hs' : s = 1 + pi2 * c := by linear_combination hc
  have h2 : s ^ 2 = 1 + pi2 ^ 2 * (pi2 * c + c ^ 2) := by
    rw [hs']; linear_combination (-(pi2 * c)) * pi2_sq
  have h4 : s ^ 4 = 1 + pi2 ^ 4 * ((pi2 * c + c ^ 2) * (1 + (pi2 * c + c ^ 2))) := by
    have : s ^ 4 = (s ^ 2) ^ 2 := by ring
    rw [this, h2]; linear_combination (-(pi2 ^ 2 * (pi2 * c + c ^ 2))) * pi2_sq
  obtain ⟨m', hm'⟩ : pi2 ∣ (pi2 * c + c ^ 2) * (1 + (pi2 * c + c ^ 2)) := by
    rcases pi2_dvd_or_sub_one (pi2 * c + c ^ 2) with h | ⟨m, hm⟩
    · exact dvd_mul_of_dvd_left h _
    · exact dvd_mul_of_dvd_right ⟨m + pi2, by linear_combination hm - pi2_sq⟩ _
  exact ⟨m', by rw [h4, hm']; ring⟩

theorem im_add_pi2sq (c₀ W : ℤ√2) : (c₀ + pi2 ^ 2 * W).im = c₀.im + 2 * W.im := by
  rw [pi2_sq, im_add, im_mul]; simp

theorem odd_im_of {V c₀ W : ℤ√2} (h : V = c₀ + pi2 ^ 2 * W) (hc : Odd c₀.im) : Odd V.im := by
  rw [h, im_add_pi2sq]; exact hc.add_even ⟨W.im, by ring⟩

theorem odd_im_mul_eps2 {u : ℤ√2} (hu : Odd u.im) : Odd (u * eps ^ 2).im := by
  rw [show eps ^ 2 = (⟨3, 2⟩ : ℤ√2) by decide, im_mul]
  show Odd (u.re * 2 + u.im * 3)
  exact Even.add_odd ⟨u.re, by ring⟩ (hu.mul (by decide))

theorem odd_im_neg {z : ℤ√2} (h : Odd z.im) : Odd (-z).im := by rw [im_neg]; exact h.neg

/-! ### Two non-square criteria over `ℤ√2` -/

/-- **Odd valuation**: `n² ≠ r²·p·U` for a prime `p ∤ U`, `r ≠ 0`. -/
theorem not_sq_odd {p : ℤ√2} (hp : Prime p) (n r U : ℤ√2) (hr : r ≠ 0) (hU : ¬p ∣ U) :
    n ^ 2 ≠ r ^ 2 * (p * U) := by
  suffices H : ∀ k : ℕ, ∀ n r : ℤ√2, r.norm.natAbs = k → r ≠ 0 → n ^ 2 ≠ r ^ 2 * (p * U) from
    H _ n r rfl hr
  intro k
  induction k using Nat.strong_induction_on with
  | _ k ih =>
    intro n r hk hr hn
    have hp0 := hp.ne_zero
    have hpn : p ∣ n := hp.dvd_of_dvd_pow (n := 2) ⟨r ^ 2 * U, by rw [hn]; ring⟩
    obtain ⟨n₁, rfl⟩ := hpn
    have h1 : p * n₁ ^ 2 = r ^ 2 * U := by
      apply mul_left_cancel₀ hp0
      linear_combination hn
    have hpr : p ∣ r := by
      have : p ∣ r ^ 2 * U := ⟨n₁ ^ 2, h1.symm⟩
      rcases hp.dvd_or_dvd this with h2 | h2
      · exact hp.dvd_of_dvd_pow h2
      · exact absurd h2 hU
    obtain ⟨r₁, rfl⟩ := hpr
    have hr₁ : r₁ ≠ 0 := right_ne_zero_of_mul hr
    have hlt : r₁.norm.natAbs < k := by rw [← hk]; exact norm_lt_of_mul hp0 hp.not_unit hr₁
    have h2 : n₁ ^ 2 = r₁ ^ 2 * (p * U) := by
      apply mul_left_cancel₀ hp0
      linear_combination h1
    exact ih _ hlt n₁ r₁ rfl hr₁ h2

/-- **Mod 2**: `n² ≠ r²·V` if `V.im` is odd (`V` is not a square modulo `2` up to `√2`-powers). -/
theorem not_sq_oddim (V : ℤ√2) (hV : Odd V.im) (n r : ℤ√2) (hr : r ≠ 0) : n ^ 2 ≠ r ^ 2 * V := by
  suffices H : ∀ k : ℕ, ∀ n r : ℤ√2, r.norm.natAbs = k → r ≠ 0 → n ^ 2 ≠ r ^ 2 * V from
    H _ n r rfl hr
  intro k
  induction k using Nat.strong_induction_on with
  | _ k ih =>
    intro n r hk hr hn
    by_cases hpr : pi2 ∣ r
    · obtain ⟨r₁, rfl⟩ := hpr
      have hpn : pi2 ∣ n := by
        apply pi2_prime.dvd_of_dvd_pow (n := 2)
        exact ⟨pi2 * r₁ ^ 2 * V, by rw [hn]; ring⟩
      obtain ⟨n₁, rfl⟩ := hpn
      have h2 : n₁ ^ 2 = r₁ ^ 2 * V := by
        have h3 : pi2 ^ 2 * n₁ ^ 2 = pi2 ^ 2 * (r₁ ^ 2 * V) := by linear_combination hn
        exact mul_left_cancel₀ (pow_ne_zero 2 pi2_ne_zero) h3
      have hr₁ : r₁ ≠ 0 := right_ne_zero_of_mul hr
      have hlt : r₁.norm.natAbs < k := by
        rw [← hk]; exact norm_lt_of_mul pi2_ne_zero pi2_not_unit hr₁
      exact ih _ hlt n₁ r₁ rfl hr₁ h2
    · rw [pi2_dvd_iff] at hpr
      have hro : Odd r.re := Int.not_even_iff_odd.mp (fun he => hpr (even_iff_two_dvd.mp he))
      have e1 : (n ^ 2).im = 2 * (n.re * n.im) := by
        simp only [sq, im_mul]; ring
      have e2 : (r ^ 2 * V).im = (r.re * r.re + 2 * r.im * r.im) * V.im + 2 * (r.re * r.im * V.re) := by
        simp only [sq, im_mul, re_mul]; ring
      have him := congrArg Zsqrtd.im hn
      rw [e1, e2] at him
      have hodd : Odd ((r.re * r.re + 2 * r.im * r.im) * V.im + 2 * (r.re * r.im * V.re)) :=
        Odd.add_even (Odd.mul (Odd.add_even (hro.mul hro) ⟨r.im * r.im, by ring⟩) hV)
          ⟨r.re * r.im * V.re, by ring⟩
      rw [← him] at hodd
      exact (Int.not_even_iff_odd.mpr hodd) ⟨n.re * n.im, by ring⟩

/-- **Real place**: `n² ≠ r²·F` if `σ(F) < 0` for an injective `σ : ℤ√2 → ℝ`. -/
theorem not_sq_neg (σ : ℤ√2 →+* ℝ) (hσ : Function.Injective σ) (n r F : ℤ√2) (hr : r ≠ 0)
    (hF : σ F < 0) : n ^ 2 ≠ r ^ 2 * F := by
  intro h
  have h1 := congrArg σ h
  simp only [map_mul, map_pow] at h1
  have hr' : σ r ≠ 0 := by intro h0; apply hr; apply hσ; rw [h0, map_zero]
  have : 0 < σ r ^ 2 := by positivity
  nlinarith [sq_nonneg (σ n)]

/-! ### From a point to the quartic over `ℤ√2` -/

theorem prim {p : ℤ√2} (hp : Prime p) : ∀ (k : ℕ) (s t : ℤ√2), t.norm.natAbs = k → t ≠ 0 →
    ∃ s' t' : ℤ√2, ∃ m : ℕ, ¬(p ∣ s' ∧ p ∣ t') ∧ s = p ^ m * s' ∧ t = p ^ m * t' := by
  intro k
  induction k using Nat.strong_induction_on with
  | _ k ih =>
    intro s t hk ht
    by_cases h : p ∣ s ∧ p ∣ t
    · obtain ⟨⟨s₁, rfl⟩, ⟨t₁, rfl⟩⟩ := h
      have ht₁ : t₁ ≠ 0 := right_ne_zero_of_mul ht
      obtain ⟨s', t', m, hst, hs, ht'⟩ :=
        ih _ (by rw [← hk]; exact norm_lt_of_mul hp.ne_zero hp.not_unit ht₁) s₁ t₁ rfl ht₁
      exact ⟨s', t', m + 1, hst, by rw [hs]; ring, by rw [ht']; ring⟩
    · exact ⟨s, t, 0, h, by ring, by ring⟩

theorem quartic_of_point {A B c : ℤ√2} {x y w : Q2} (hx0 : x ≠ 0)
    (hxy : y ^ 2 = x * (x ^ 2 + iota2 A * x + iota2 B)) (hxw : x = iota2 c * w ^ 2) {p : ℤ√2}
    (hp : Prime p) :
    ∃ n r s t : ℤ√2, r ≠ 0 ∧ s ≠ 0 ∧ t ≠ 0 ∧ ¬(p ∣ s ∧ p ∣ t) ∧ n ^ 2 = r ^ 2 * Fq A B c s t := by
  have hw0 : w ≠ 0 := by rintro rfl; apply hx0; rw [hxw]; ring
  obtain ⟨s₀, t₀, ht₀, hst₀⟩ := Q2_frac w
  obtain ⟨s, t, m, hprim, hs, ht⟩ := prim hp _ s₀ t₀ rfl ht₀
  have hpm : iota2 (p ^ m) ≠ 0 := iota2_ne_zero (pow_ne_zero _ hp.ne_zero)
  have ht0 : t ≠ 0 := by rintro rfl; apply ht₀; rw [ht, mul_zero]
  have hT : iota2 t ≠ 0 := iota2_ne_zero ht0
  have hw : w = iota2 s / iota2 t := by
    rw [hst₀, hs, ht, map_mul, map_mul, mul_div_mul_left _ _ hpm]
  have hs0 : s ≠ 0 := by rintro rfl; apply hw0; rw [hw, map_zero, zero_div]
  have hS : iota2 s ≠ 0 := iota2_ne_zero hs0
  obtain ⟨n, r, hr, hnr⟩ := Q2_frac (y * iota2 t ^ 3 / iota2 s)
  have hR : iota2 r ≠ 0 := iota2_ne_zero hr
  refine ⟨n, r, s, t, hr, hs0, ht0, hprim, ?_⟩
  apply iota2_injective
  have key : (y * iota2 t ^ 3 / iota2 s) ^ 2 = iota2 (Fq A B c s t) := by
    unfold Fq
    simp only [map_mul, map_add, map_pow]
    rw [div_pow, mul_pow, hxy, hxw, hw]
    field_simp
    try ring
  rw [hnr, div_pow] at key
  simp only [map_mul, map_pow]
  rw [← key]
  field_simp

/-- Every `x ≠ 0` in `ℚ(√2)` is `ι(δ)·w²` with `δ` squarefree. -/
theorem exists_sqfree (x : Q2) (hx : x ≠ 0) :
    ∃ δ : ℤ√2, ∃ w : Q2, δ ≠ 0 ∧ Squarefree δ ∧ x = iota2 δ * w ^ 2 := by
  obtain ⟨p, q, hq, hpq⟩ := Q2_frac x
  have hp : p ≠ 0 := by rintro rfl; apply hx; rw [hpq, map_zero, zero_div]
  have hQ : iota2 q ≠ 0 := iota2_ne_zero hq
  suffices H : ∀ k : ℕ, ∀ δ : ℤ√2, δ.norm.natAbs = k → δ ≠ 0 → (∃ w : Q2, x = iota2 δ * w ^ 2) →
      ∃ δ : ℤ√2, ∃ w : Q2, δ ≠ 0 ∧ Squarefree δ ∧ x = iota2 δ * w ^ 2 from
    H _ (p * q) rfl (mul_ne_zero hp hq) ⟨1 / iota2 q, by rw [hpq, map_mul]; field_simp; try ring⟩
  intro k
  induction k using Nat.strong_induction_on with
  | _ k ih =>
    rintro δ hk hδ ⟨w, hw⟩
    by_cases hsq : Squarefree δ
    · exact ⟨δ, w, hδ, hsq, hw⟩
    · unfold Squarefree at hsq
      push_neg at hsq
      obtain ⟨c, ⟨δ', hδ'⟩, hc⟩ := hsq
      have hc0 : c ≠ 0 := by rintro rfl; apply hδ; rw [hδ']; ring
      have hδ'0 : δ' ≠ 0 := by rintro rfl; apply hδ; rw [hδ', mul_zero]
      have hlt : δ'.norm.natAbs < k := by
        rw [← hk, hδ', mul_assoc]
        exact (norm_lt_of_mul hc0 hc hδ'0).trans (norm_lt_of_mul hc0 hc (mul_ne_zero hc0 hδ'0))
      exact ih _ hlt δ' rfl hδ'0 ⟨iota2 c * w, by rw [hw, hδ']; simp only [map_mul]; ring⟩

/-! ### A prime factor of a squarefree `δ` divides `B` -/

theorem prime_dvd_B {A B c : ℤ√2} (hsq : Squarefree c) {p : ℤ√2} (hp : Prime p) (hpc : p ∣ c)
    (hpB : ¬p ∣ B) {n r s t : ℤ√2} (hr : r ≠ 0) (hst : ¬(p ∣ s ∧ p ∣ t))
    (h : n ^ 2 = r ^ 2 * Fq A B c s t) : False := by
  obtain ⟨c', rfl⟩ := hpc
  have hc' : ¬p ∣ c' := fun hd => hp.not_unit (hsq p (mul_dvd_mul_left p hd))
  unfold Fq at h
  by_cases ht : p ∣ t
  · have hs : ¬p ∣ s := fun hs => hst ⟨hs, ht⟩
    obtain ⟨t₁, rfl⟩ := ht
    exact not_sq_odd hp n (r * p) (c' ^ 3 * s ^ 4 + p * (A * c' ^ 2 * s ^ 2 * t₁ ^ 2 + p * B * c' * t₁ ^ 4))
      (mul_ne_zero hr hp.ne_zero) (nd_add (nd_mul hp (nd_pow hp hc' 3) (nd_pow hp hs 4)))
      (by rw [h]; ring)
  · refine not_sq_odd hp n r (c' * (p * (p * c' ^ 2 * s ^ 4 + A * c' * s ^ 2 * t ^ 2) + B * t ^ 4)) hr ?_
      (by rw [h]; ring)
    refine nd_mul hp hc' fun hd => ?_
    have : p ∣ B * t ^ 4 := (dvd_add_right (dvd_mul_right _ _)).mp hd
    exact nd_mul hp hpB (nd_pow hp ht 4) this

/-! ### Squarefree `δ` with prime factors among `√2`, `3` -/

theorem unit_of_no_prime {δ : ℤ√2} (hδ0 : δ ≠ 0) (h : ∀ p, Prime p → ¬p ∣ δ) : IsUnit δ := by
  by_contra hu
  obtain ⟨i, hi, hid⟩ := WfDvdMonoid.exists_irreducible_factor hu hδ0
  exact h i (UniqueFactorizationMonoid.irreducible_iff_prime.mp hi) hid

theorem classify1 {δ : ℤ√2} (hsq : Squarefree δ) (hδ0 : δ ≠ 0)
    (h : ∀ p, Prime p → p ∣ δ → Associated p pi2) : IsUnit δ ∨ ∃ u, IsUnit u ∧ δ = pi2 * u := by
  by_cases h2 : pi2 ∣ δ
  · obtain ⟨u, rfl⟩ := h2
    right
    refine ⟨u, unit_of_no_prime (right_ne_zero_of_mul hδ0) fun p hp hpu => ?_, rfl⟩
    have hpa := h p hp (dvd_mul_of_dvd_right hpu _)
    have : pi2 ∣ u := hpa.dvd_iff_dvd_left.mp hpu
    exact pi2_not_unit (hsq pi2 (mul_dvd_mul_left pi2 this))
  · left
    exact unit_of_no_prime hδ0 fun p hp hpd => h2 ((h p hp hpd).dvd_iff_dvd_left.mp hpd)

theorem classify2 {δ : ℤ√2} (hsq : Squarefree δ) (hδ0 : δ ≠ 0)
    (h : ∀ p, Prime p → p ∣ δ → Associated p pi2 ∨ Associated p three) :
    ∃ u : ℤ√2, ∃ i j : ℕ, IsUnit u ∧ i < 2 ∧ j < 2 ∧ δ = u * pi2 ^ i * three ^ j := by
  obtain ⟨δ₁, i, hi, hδ₁, hπ⟩ : ∃ δ₁ : ℤ√2, ∃ i : ℕ, i < 2 ∧ δ = pi2 ^ i * δ₁ ∧ ¬pi2 ∣ δ₁ := by
    by_cases h2 : pi2 ∣ δ
    · obtain ⟨d1, hd1⟩ := h2
      refine ⟨d1, 1, by norm_num, by rw [hd1]; ring, fun hd => pi2_not_unit (hsq pi2 ?_)⟩
      rw [hd1]; exact mul_dvd_mul_left pi2 hd
    · exact ⟨δ, 0, by norm_num, by ring, h2⟩
  obtain ⟨δ₂, j, hj, hδ₂, h3⟩ : ∃ δ₂ : ℤ√2, ∃ j : ℕ, j < 2 ∧ δ₁ = three ^ j * δ₂ ∧ ¬three ∣ δ₂ := by
    by_cases h3 : three ∣ δ₁
    · obtain ⟨d2, hd2⟩ := h3
      refine ⟨d2, 1, by norm_num, by rw [hd2]; ring, fun hd => three_not_unit (hsq three ?_)⟩
      rw [hδ₁, hd2]; exact dvd_mul_of_dvd_right (mul_dvd_mul_left three hd) _
    · exact ⟨δ₁, 0, by norm_num, by ring, h3⟩
  have hδ₂0 : δ₂ ≠ 0 := by rintro rfl; apply hδ0; rw [hδ₁, hδ₂]; ring
  refine ⟨δ₂, i, j, ?_, hi, hj, by rw [hδ₁, hδ₂]; ring⟩
  refine unit_of_no_prime hδ₂0 fun p hp hpd => ?_
  have hpδ : p ∣ δ := by
    rw [hδ₁, hδ₂]; exact dvd_mul_of_dvd_right (dvd_mul_of_dvd_right hpd _) _
  rcases h p hp hpδ with ha | ha
  · exact hπ (by rw [hδ₂]; exact dvd_mul_of_dvd_right (ha.dvd_iff_dvd_left.mp hpd) _)
  · exact h3 (ha.dvd_iff_dvd_left.mp hpd)

/-! ### Local exclusions -/

/-- `E1`, `δ = √2·w` with `√2 ∤ w`: odd `√2`-valuation. -/
theorem V1 {w n r s t : ℤ√2} (hw : ¬pi2 ∣ w) (hr : r ≠ 0) (hst : ¬(pi2 ∣ s ∧ pi2 ∣ t)) :
    n ^ 2 ≠ r ^ 2 * Fq A1 B1 (pi2 * w) s t := by
  intro h
  unfold Fq A1 B1 at h
  have hp := pi2_prime
  by_cases hs : pi2 ∣ s
  · have ht : ¬pi2 ∣ t := fun ht => hst ⟨hs, ht⟩
    obtain ⟨s₁, rfl⟩ := hs
    by_cases hs₁ : pi2 ∣ s₁
    · obtain ⟨s₂, rfl⟩ := hs₁
      exact not_sq_odd hp n (r * pi2 ^ 4)
        (-(eps ^ 2 * w * t ^ 4) + pi2 * (pi2 * (w ^ 3 * s₂ ^ 4 + pi2 ^ 2 * eps * w ^ 2 * s₂ ^ 2 * t ^ 2)))
        (mul_ne_zero hr (pow_ne_zero _ pi2_ne_zero))
        (nd_add (nd_neg (nd_mul hp (nd_mul hp (nd_pow hp pi2_nd_eps 2) hw) (nd_pow hp ht 4))))
        (by rw [h]; ring)
    · exact not_sq_odd hp n (r * pi2 ^ 3)
        (w ^ 3 * s₁ ^ 4 + pi2 * (pi2 * (pi2 ^ 2 * eps * w ^ 2 * s₁ ^ 2 * t ^ 2 - eps ^ 2 * w * t ^ 4)))
        (mul_ne_zero hr (pow_ne_zero _ pi2_ne_zero))
        (nd_add (nd_mul hp (nd_pow hp hw 3) (nd_pow hp hs₁ 4)))
        (by rw [h]; ring)
  · exact not_sq_odd hp n (r * pi2)
      (w ^ 3 * s ^ 4 + pi2 * (pi2 ^ 5 * (eps * w ^ 2 * s ^ 2 * t ^ 2 - eps ^ 2 * w * t ^ 4)))
      (mul_ne_zero hr pi2_ne_zero)
      (nd_add (nd_mul hp (nd_pow hp hw 3) (nd_pow hp hs 4)))
      (by rw [h]; ring)

/-- `E1`, `δ = √2⁴·u` with `u² = ε²`, `u.im` odd (that is, `δ ≡ ±ε` modulo squares): mod 2. -/
theorem M1 {u n r s t : ℤ√2} (hu2 : u ^ 2 = eps ^ 2) (hui : Odd u.im) (hpu : ¬pi2 ∣ u) (hr : r ≠ 0)
    (hst : ¬(pi2 ∣ s ∧ pi2 ∣ t)) : n ^ 2 ≠ r ^ 2 * Fq A1 B1 (pi2 ^ 4 * u) s t := by
  intro h
  have hp := pi2_prime
  have hF : Fq A1 B1 (pi2 ^ 4 * u) s t
      = (pi2 ^ 6) ^ 2 * (u * (u ^ 2 * s ^ 4 + pi2 ^ 3 * eps * u * s ^ 2 * t ^ 2 - eps ^ 2 * t ^ 4)) := by
    unfold Fq A1 B1; ring
  rw [hF] at h
  have hr6 : r * pi2 ^ 6 ≠ 0 := mul_ne_zero hr (pow_ne_zero _ pi2_ne_zero)
  by_cases hs : pi2 ∣ s
  · have ht : ¬pi2 ∣ t := fun ht => hst ⟨hs, ht⟩
    obtain ⟨s₁, rfl⟩ := hs
    obtain ⟨k₂, hk₂⟩ := pow4_of_not_dvd ht
    have hV : u * (u ^ 2 * (pi2 * s₁) ^ 4 + pi2 ^ 3 * eps * u * (pi2 * s₁) ^ 2 * t ^ 2 - eps ^ 2 * t ^ 4)
        = -(u * eps ^ 2) + pi2 ^ 2 * (u * (pi2 ^ 2 * u ^ 2 * s₁ ^ 4 + pi2 ^ 3 * eps * u * s₁ ^ 2 * t ^ 2
          - pi2 ^ 3 * eps ^ 2 * k₂)) := by
      linear_combination (-(u * eps ^ 2)) * hk₂
    exact not_sq_oddim _ (odd_im_of hV (odd_im_neg (odd_im_mul_eps2 hui))) n (r * pi2 ^ 6) hr6
      (by rw [h]; ring)
  · obtain ⟨k₁, hk₁⟩ := pow4_of_not_dvd hs
    by_cases ht : pi2 ∣ t
    · obtain ⟨t₁, rfl⟩ := ht
      have hV : u * (u ^ 2 * s ^ 4 + pi2 ^ 3 * eps * u * s ^ 2 * (pi2 * t₁) ^ 2 - eps ^ 2 * (pi2 * t₁) ^ 4)
          = u * eps ^ 2 + pi2 ^ 2 * (pi2 ^ 3 * u ^ 3 * k₁ + pi2 ^ 3 * eps * u ^ 2 * s ^ 2 * t₁ ^ 2
            - pi2 ^ 2 * u * eps ^ 2 * t₁ ^ 4) := by
        linear_combination u ^ 3 * hk₁ + u * hu2
      exact not_sq_oddim _ (odd_im_of hV (odd_im_mul_eps2 hui)) n (r * pi2 ^ 6) hr6
        (by rw [h]; ring)
    · obtain ⟨k₂, hk₂⟩ := pow4_of_not_dvd ht
      have hV : u * (u ^ 2 * s ^ 4 + pi2 ^ 3 * eps * u * s ^ 2 * t ^ 2 - eps ^ 2 * t ^ 4)
          = pi2 ^ 3 * (u * (eps * u * s ^ 2 * t ^ 2 + pi2 * (pi2 * eps ^ 2 * (k₁ - k₂)))) := by
        linear_combination u * s ^ 4 * hu2 + u * eps ^ 2 * hk₁ - u * eps ^ 2 * hk₂
      exact not_sq_odd hp n (r * pi2 ^ 6 * pi2) (u * (eps * u * s ^ 2 * t ^ 2 + pi2 * (pi2 * eps ^ 2 * (k₁ - k₂))))
        (mul_ne_zero hr6 pi2_ne_zero)
        (nd_mul hp hpu (nd_add (nd_mul hp (nd_mul hp (nd_mul hp pi2_nd_eps hpu) (nd_pow hp hs 2))
          (nd_pow hp ht 2))))
        (by rw [h, hV]; ring)

/-- `E1'`, `δ = √2·w` with `√2 ∤ w`: odd `√2`-valuation. -/
theorem V2 {w n r s t : ℤ√2} (hw : ¬pi2 ∣ w) (hr : r ≠ 0) (hst : ¬(pi2 ∣ s ∧ pi2 ∣ t)) :
    n ^ 2 ≠ r ^ 2 * Fq A2 B2 (pi2 * w) s t := by
  intro h
  unfold Fq A2 B2 at h
  have hp := pi2_prime
  by_cases hs : pi2 ∣ s
  · have ht : ¬pi2 ∣ t := fun ht => hst ⟨hs, ht⟩
    obtain ⟨s₁, rfl⟩ := hs
    by_cases hs₁ : pi2 ∣ s₁
    · obtain ⟨s₂, rfl⟩ := hs₁
      by_cases hs₂ : pi2 ∣ s₂
      · obtain ⟨s₃, rfl⟩ := hs₂
        exact not_sq_odd hp n (r * pi2 ^ 6)
          (three * eps ^ 2 * w * t ^ 4 + pi2 * (pi2 * (w ^ 3 * s₃ ^ 4 - pi2 ^ 2 * eps * w ^ 2 * s₃ ^ 2 * t ^ 2)))
          (mul_ne_zero hr (pow_ne_zero _ pi2_ne_zero))
          (nd_add (nd_mul hp (nd_mul hp (nd_mul hp pi2_nd_three (nd_pow hp pi2_nd_eps 2)) hw)
            (nd_pow hp ht 4)))
          (by rw [h]; ring)
      · exact not_sq_odd hp n (r * pi2 ^ 5)
          (w ^ 3 * s₂ ^ 4 + pi2 * (pi2 * (-(pi2 ^ 2 * eps * w ^ 2 * s₂ ^ 2 * t ^ 2) + three * eps ^ 2 * w * t ^ 4)))
          (mul_ne_zero hr (pow_ne_zero _ pi2_ne_zero))
          (nd_add (nd_mul hp (nd_pow hp hw 3) (nd_pow hp hs₂ 4)))
          (by rw [h]; ring)
    · exact not_sq_odd hp n (r * pi2 ^ 3)
        (w ^ 3 * s₁ ^ 4 + pi2 * (pi2 ^ 5 * (-(eps * w ^ 2 * s₁ ^ 2 * t ^ 2) + three * eps ^ 2 * w * t ^ 4)))
        (mul_ne_zero hr (pow_ne_zero _ pi2_ne_zero))
        (nd_add (nd_mul hp (nd_pow hp hw 3) (nd_pow hp hs₁ 4)))
        (by rw [h]; ring)
  · exact not_sq_odd hp n (r * pi2)
      (w ^ 3 * s ^ 4 + pi2 * (pi2 ^ 7 * (-(eps * w ^ 2 * s ^ 2 * t ^ 2) + pi2 ^ 2 * three * eps ^ 2 * w * t ^ 4)))
      (mul_ne_zero hr pi2_ne_zero)
      (nd_add (nd_mul hp (nd_pow hp hw 3) (nd_pow hp hs 4)))
      (by rw [h]; ring)

/-- The two real embeddings. -/
noncomputable def sig1 : ℤ√2 →+* ℝ := Zsqrtd.lift ⟨√2, by rw [Real.mul_self_sqrt (by norm_num)]; norm_num⟩
noncomputable def sig2 : ℤ√2 →+* ℝ :=
  Zsqrtd.lift ⟨-√2, by rw [neg_mul_neg, Real.mul_self_sqrt (by norm_num)]; norm_num⟩

theorem sig1_apply (z : ℤ√2) : sig1 z = z.re + z.im * √2 := by simp [sig1]
theorem sig2_apply (z : ℤ√2) : sig2 z = z.re + z.im * (-√2) := by simp [sig2]
theorem sig1_inj : Function.Injective sig1 := Zsqrtd.lift_injective _ nonsq
theorem sig2_inj : Function.Injective sig2 := Zsqrtd.lift_injective _ nonsq

theorem sqrt2_bounds : (1 : ℝ) < √2 ∧ √2 < 3 / 2 := by
  have h := Real.sq_sqrt (show (0 : ℝ) ≤ 2 by norm_num)
  have h0 := Real.sqrt_nonneg 2
  constructor <;> nlinarith

/-- `E1'` at a real place where `δ` is negative. -/
theorem real_excl (σ : ℤ√2 →+* ℝ) (hσ : Function.Injective σ) {c n r s t : ℤ√2} (hr : r ≠ 0) (ht : t ≠ 0)
    (hc : σ c < 0) (hq : σ (⟨768, 512⟩ : ℤ√2) > 0) : n ^ 2 ≠ r ^ 2 * Fq A2 B2 c s t := by
  refine not_sq_neg σ hσ n r _ hr ?_
  have key : 4 * (c ^ 2 * s ^ 4 + A2 * c * s ^ 2 * t ^ 2 + B2 * t ^ 4)
      = (2 * c * s ^ 2 + A2 * t ^ 2) ^ 2 + (4 * B2 - A2 ^ 2) * t ^ 4 := by ring
  have hC : 4 * B2 - A2 ^ 2 = (⟨768, 512⟩ : ℤ√2) := by unfold A2 B2 three pi2 eps; decide
  rw [hC] at key
  have h4 := congrArg σ key
  simp only [map_mul, map_add, map_pow, map_ofNat] at h4
  have ht' : σ t ≠ 0 := by intro h0; apply ht; apply hσ; rw [h0, map_zero]
  have hpos : 0 < σ (c ^ 2 * s ^ 4 + A2 * c * s ^ 2 * t ^ 2 + B2 * t ^ 4) := by
    simp only [map_mul, map_add, map_pow]
    have : 0 < σ t ^ 4 := by positivity
    nlinarith [sq_nonneg (2 * σ c * σ s ^ 2 + σ A2 * σ t ^ 2), mul_pos hq this]
  unfold Fq
  rw [map_mul]
  exact mul_neg_of_neg_of_pos hc hpos

theorem sig1_q : sig1 (⟨768, 512⟩ : ℤ√2) > 0 := by
  rw [sig1_apply]; dsimp only; push_cast; have := sqrt2_bounds; nlinarith

theorem sig2_q : sig2 (⟨768, 512⟩ : ℤ√2) > 0 := by
  rw [sig2_apply]; dsimp only; push_cast; have := sqrt2_bounds; nlinarith

theorem sig1_pi2 : sig1 pi2 = √2 := by rw [sig1_apply]; simp [pi2]
theorem sig1_three : sig1 three = 3 := by rw [sig1_apply]; simp [three]
theorem sig1_eps : sig1 eps = 1 + √2 := by rw [sig1_apply]; simp [eps]
theorem sig2_pi2 : sig2 pi2 = -√2 := by rw [sig2_apply]; simp [pi2]
theorem sig2_three : sig2 three = 3 := by rw [sig2_apply]; simp [three]
theorem sig2_eps : sig2 eps = 1 - √2 := by rw [sig2_apply]; simp [eps]; try ring

/-! ### The constants in `ℚ(√2)` -/

theorem iota2_A1 : iota2 A1 = 8 * (2 + s2) := by
  rw [show A1 = (⟨16, 8⟩ : ℤ√2) by unfold A1 pi2 eps; decide, iota2_apply]
  dsimp only; push_cast; ring

theorem iota2_B1 : iota2 B1 = -8 * (2 + s2) ^ 2 := by
  rw [show B1 = (⟨-48, -32⟩ : ℤ√2) by unfold B1 pi2 eps; decide, iota2_apply]
  dsimp only; push_cast; linear_combination 8 * s2_mul_s2

theorem iota2_A2 : iota2 A2 = -2 * (8 * (2 + s2)) := by
  rw [show A2 = (⟨-32, -16⟩ : ℤ√2) by unfold A2 pi2 eps; decide, iota2_apply]
  dsimp only; push_cast; ring

theorem iota2_B2 : iota2 B2 = (8 * (2 + s2)) ^ 2 - 4 * (-8 * (2 + s2) ^ 2) := by
  rw [show B2 = (⟨576, 384⟩ : ℤ√2) by unfold B2 three pi2 eps; decide, iota2_apply]
  dsimp only; push_cast; linear_combination (-96) * s2_mul_s2

theorem iota2_three : iota2 three = 3 := by rw [iota2_apply]; simp [three]

theorem isUnit_eps : IsUnit eps :=
  ⟨⟨eps, epsinv, eps_mul_epsinv, (mul_comm epsinv eps).trans eps_mul_epsinv⟩, rfl⟩

/-! ### The two descent statements -/

/-- **Descent on `E1`**: every affine point with `x ≠ 0` has `x` or `b/x` a square. -/
theorem descentE1 : DescentW (8 * (2 + s2) : Q2) (-8 * (2 + s2) ^ 2) := by
  intro x y hxy hx
  rw [← iota2_A1, ← iota2_B1] at hxy
  obtain ⟨δ, w, hδ0, hsq, hxw⟩ := exists_sqfree x hx
  have hprime : ∀ p, Prime p → p ∣ δ → Associated p pi2 := by
    intro p hp hpδ
    by_cases hpB : p ∣ B1
    · have h1 : p ∣ pi2 ^ 8 * (-(eps ^ 2)) := by
        rw [show pi2 ^ 8 * (-(eps ^ 2)) = B1 by unfold B1; ring]; exact hpB
      rcases hp.dvd_or_dvd h1 with h2 | h2
      · exact hp.irreducible.associated_of_dvd pi2_prime.irreducible (hp.dvd_of_dvd_pow h2)
      · exact absurd h2 (nd_unit hp (isUnit_eps.pow 2).neg)
    · exfalso
      obtain ⟨n, r, s, t, hr, -, -, hst, h⟩ := quartic_of_point hx hxy hxw hp
      exact prime_dvd_B hsq hp hpδ hpB hr hst h
  have hpi : iota2 pi2 ≠ 0 := iota2_ne_zero pi2_ne_zero
  rcases classify1 hsq hδ0 hprime with hu | ⟨u, hu, rfl⟩
  · obtain ⟨v, hv | hv | hv | hv⟩ := unit_classes hu
    · left; exact ⟨iota2 v * w, by rw [hxw, hv, map_pow]; ring⟩
    · right
      have hW0 : iota2 v * w ≠ 0 := by
        intro h0; apply hx; rw [hxw, hv, map_neg, map_pow]; linear_combination (-(iota2 v * w)) * h0
      refine ⟨2 * s2 * (2 + s2) / (iota2 v * w), ?_⟩
      rw [hxw, hv, map_neg, map_pow]
      have e : -iota2 v ^ 2 * w ^ 2 * (2 * s2 * (2 + s2) / (iota2 v * w)) ^ 2
          = -(4 * (s2 * s2) * (2 + s2) ^ 2) * ((iota2 v * w) ^ 2 / (iota2 v * w) ^ 2) := by ring
      rw [e, div_self (pow_ne_zero 2 hW0), s2_mul_s2]; ring
    · exfalso
      have hxw' : x = iota2 (pi2 ^ 4 * eps) * (iota2 v * w / iota2 (pi2 ^ 2)) ^ 2 := by
        rw [hxw, hv]; simp only [map_mul, map_pow]
        rw [div_pow, eq_comm, mul_div_assoc', div_eq_iff (pow_ne_zero 2 (pow_ne_zero 2 hpi))]; ring
      obtain ⟨n, r, s, t, hr, -, -, hst, h⟩ := quartic_of_point hx hxy hxw' pi2_prime
      exact M1 (u := eps) rfl ⟨0, rfl⟩ pi2_nd_eps hr hst h
    · exfalso
      have hxw' : x = iota2 (pi2 ^ 4 * (-eps)) * (iota2 v * w / iota2 (pi2 ^ 2)) ^ 2 := by
        rw [hxw, hv]; simp only [map_mul, map_pow, map_neg]
        rw [div_pow, eq_comm, mul_div_assoc', div_eq_iff (pow_ne_zero 2 (pow_ne_zero 2 hpi))]; ring
      obtain ⟨n, r, s, t, hr, -, -, hst, h⟩ := quartic_of_point hx hxy hxw' pi2_prime
      exact M1 (u := -eps) (by ring) ⟨-1, rfl⟩ (nd_neg pi2_nd_eps) hr hst h
  · exfalso
    obtain ⟨n, r, s, t, hr, -, -, hst, h⟩ := quartic_of_point hx hxy hxw pi2_prime
    exact V1 (nd_unit pi2_prime hu) hr hst h

/-- **Descent on `E1'`**: every affine point with `X ≠ 0` has `X` or `(a² − 4b)/X` a square. -/
theorem descentE1' : DescentW' (8 * (2 + s2) : Q2) (-8 * (2 + s2) ^ 2) := by
  intro X Y hXY hX
  rw [← iota2_A2, ← iota2_B2] at hXY
  obtain ⟨δ, w, hδ0, hsq, hxw⟩ := exists_sqfree X hX
  have hprime : ∀ p, Prime p → p ∣ δ → Associated p pi2 ∨ Associated p three := by
    intro p hp hpδ
    by_cases hpB : p ∣ B2
    · have hpB' := hpB
      unfold B2 at hpB'
      rcases hp.dvd_or_dvd hpB' with h1 | h1
      · rcases hp.dvd_or_dvd h1 with h2 | h2
        · exact Or.inr (hp.irreducible.associated_of_dvd three_prime.irreducible h2)
        · exact Or.inl (hp.irreducible.associated_of_dvd pi2_prime.irreducible (hp.dvd_of_dvd_pow h2))
      · exact absurd (hp.dvd_of_dvd_pow h1) (nd_unit hp isUnit_eps)
    · exfalso
      obtain ⟨n, r, s, t, hr, -, -, hst, h⟩ := quartic_of_point hX hXY hxw hp
      exact prime_dvd_B hsq hp hpδ hpB hr hst h
  obtain ⟨u, i, j, hu, hi, hj, rfl⟩ := classify2 hsq hδ0 hprime
  obtain ⟨v, hv | hv | hv | hv⟩ := unit_classes hu
  all_goals
    have hW : ∀ c : ℤ√2, u * pi2 ^ i * three ^ j = c * v ^ 2 →
        X = iota2 c * (iota2 v * w) ^ 2 := by
      intro c hc; rw [hxw, hc, map_mul, map_pow]; ring
    have hreal : ∀ c : ℤ√2, ∀ σ : ℤ√2 →+* ℝ, Function.Injective σ → σ (⟨768, 512⟩ : ℤ√2) > 0 →
        σ c < 0 → X = iota2 c * (iota2 v * w) ^ 2 → False := by
      intro c σ hσ hq hc hXc
      obtain ⟨n, r, s, t, hr, -, ht, -, h⟩ := quartic_of_point hX hXY hXc pi2_prime
      exact real_excl σ hσ hr ht hc hq h
    obtain ⟨hb1, hb2⟩ := sqrt2_bounds
  · interval_cases i
    · interval_cases j
      · left
        exact ⟨iota2 v * w, by rw [hW 1 (by rw [hv]; ring), map_one, one_mul]⟩
      · right
        have hXc := hW three (by rw [hv]; ring)
        have hW0 : iota2 v * w ≠ 0 := by
          intro h0; apply hX; rw [hXc, h0]; ring
        refine ⟨4 * s2 * (2 + s2) / (iota2 v * w), ?_⟩
        have e : 3 * (iota2 v * w) ^ 2 * (4 * s2 * (2 + s2) / (iota2 v * w)) ^ 2
            = 48 * (s2 * s2) * (2 + s2) ^ 2 * ((iota2 v * w) ^ 2 / (iota2 v * w) ^ 2) := by ring
        rw [hXc, iota2_three, e, div_self (pow_ne_zero 2 hW0), s2_mul_s2]; ring
    · exact False.elim <| hreal (pi2 * three ^ j) sig2 sig2_inj sig2_q
        (by rw [map_mul, map_pow, sig2_pi2, sig2_three]
            exact mul_neg_of_neg_of_pos (by linarith) (by positivity))
        (hW _ (by rw [hv]; ring))
  · exact False.elim <| hreal (-(pi2 ^ i * three ^ j)) sig1 sig1_inj sig1_q
      (by rw [map_neg, map_mul, map_pow, map_pow, sig1_pi2, sig1_three]
          have : 0 < √2 ^ i * (3 : ℝ) ^ j := by positivity
          linarith)
      (hW _ (by rw [hv]; ring))
  · interval_cases i
    · exact False.elim <| hreal (eps * three ^ j) sig2 sig2_inj sig2_q
        (by rw [map_mul, map_pow, sig2_eps, sig2_three]
            exact mul_neg_of_neg_of_pos (by linarith) (by positivity))
        (hW _ (by rw [hv]; ring))
    · obtain ⟨n, r, s, t, hr, -, -, hst, h⟩ :=
        quartic_of_point hX hXY (hW (pi2 * (eps * three ^ j)) (by rw [hv]; ring)) pi2_prime
      exact (V2 (nd_mul pi2_prime pi2_nd_eps (nd_pow pi2_prime pi2_nd_three j)) hr hst h).elim
  · exact False.elim <| hreal (-(eps * (pi2 ^ i * three ^ j))) sig1 sig1_inj sig1_q
      (by rw [map_neg, map_mul, map_mul, map_pow, map_pow, sig1_pi2, sig1_three, sig1_eps]
          have : 0 < (1 + √2) * (√2 ^ i * (3 : ℝ) ^ j) := by positivity
          linarith)
      (hW _ (by rw [hv]; ring))

/-- **`rank E1(ℚ(√2)) = 0`**: `E1 : y² = x(x² + 8(2 + √2)x − 8(2 + √2)²)` over `ℚ⟮√2⟯`. -/
theorem E1K_finrank_eq_zero : Module.finrank ℤ E1K.Point = 0 :=
  E1K_finrank_eq_zero_of_descent descentE1 descentE1'

end MagicSquareOfSquares.E1Descent

#print axioms MagicSquareOfSquares.E1Descent.not_sq_odd
#print axioms MagicSquareOfSquares.E1Descent.not_sq_oddim
#print axioms MagicSquareOfSquares.E1Descent.quartic_of_point
#print axioms MagicSquareOfSquares.E1Descent.exists_sqfree
#print axioms MagicSquareOfSquares.E1Descent.prime_dvd_B
#print axioms MagicSquareOfSquares.E1Descent.classify2
#print axioms MagicSquareOfSquares.E1Descent.V1
#print axioms MagicSquareOfSquares.E1Descent.M1
#print axioms MagicSquareOfSquares.E1Descent.V2
#print axioms MagicSquareOfSquares.E1Descent.real_excl
#print axioms MagicSquareOfSquares.E1Descent.descentE1
#print axioms MagicSquareOfSquares.E1Descent.descentE1'
#print axioms MagicSquareOfSquares.E1Descent.E1K_finrank_eq_zero
