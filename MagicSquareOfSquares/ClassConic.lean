import MagicSquareOfSquares.ClassLocal

/-!
# The conic of section 4.2 of `COMPLETE.md`, without the bound of Hasse

For a prime `p ≡ 1 (mod 4)` the proof by hand of Theorem C needs, for each sign `c`, a `v` modulo `p`
with `2(1 + v²)` a nonzero square and `χ(2v) = c`. `COMPLETE.md` gets it from the bound of Hasse for
`y² = v³ + v`. Here it is proved for `p ≥ 37` from an elementary identity of Jacobsthal:

* `sum_qc_pair`: `∑_k χ((k + a)(k + b))` is `p − 1` if `a = b` and `−1` otherwise.
* `jac_sum_sq`: for `S(k) = ∑_v χ(v(k + v²))`, `∑_k S(k)² = 2p(p − 1)` when `−1` is a square.
* `jac_sq`: `S(t²)² = S(1)²` for `t ≠ 0`; so `3·S(1)² ≤ 2p(p − 1)` (the squares 1, 4, 9).
* `exists_conic_point`: the point of the conic, for `p ≡ 1 (mod 4)`, `p ≥ 37`.
* `unit_passes_of_one_mod_four`: at such a prime dividing `m`, **every** `d` prime to `p` passes
  at some point `(1, v)` modulo every power of `p` (with `pass_first_order`).

* `unit_passes_of_one_mod_four_ne_five`: the same for every prime `≡ 1 (mod 4)` except 5
  (13, 17, 29 by explicit points; 5 is the exception of the law at 5).
-/

namespace MagicSquareOfSquares.CoverDescent

section Jacobsthal

variable {p : ℕ} [Fact p.Prime]

local notation "χ" => quadraticChar (ZMod p)

theorem zmod_intCast_ne_zero {n : ℤ} (h0 : 0 < n) (hn : n < p) : ((n : ℤ) : ZMod p) ≠ 0 := by
  intro h
  have h1 := (ZMod.intCast_zmod_eq_zero_iff_dvd n p).mp h
  have h2 := Int.le_of_dvd h0 h1
  omega

theorem zmod_two_ne_zero (hp2 : p ≠ 2) : (2 : ZMod p) ≠ 0 := by
  intro h
  have h1 : (((2 : ℤ) : ℤ) : ZMod p) = 0 := by exact_mod_cast h
  exact not_dvd_two_of_odd Fact.out hp2 ((ZMod.intCast_zmod_eq_zero_iff_dvd 2 p).mp h1)

theorem sum_qc (hp2 : p ≠ 2) : ∑ a : ZMod p, χ a = 0 :=
  quadraticChar_sum_zero (by rw [ZMod.ringChar_zmod_n]; exact hp2)

theorem sum_qc_shift (hp2 : p ≠ 2) {c : ZMod p} (hc : c ≠ 0) :
    ∑ k : ZMod p, χ (k * (k + c)) = -1 := by
  have key : ∀ k : ZMod p,
      χ (k * (k + c)) = χ (1 + c * k⁻¹) - (if k⁻¹ = 0 then (1 : ℤ) else 0) := by
    intro k
    by_cases hk : k = 0
    · subst hk
      simp
    · have e : k * (k + c) = k ^ 2 * (1 + c * k⁻¹) := by
        have h := mul_inv_cancel₀ hk
        linear_combination (-(k * c)) * h
      rw [e, map_mul, quadraticChar_sq_one' hk, one_mul, if_neg (inv_ne_zero hk), sub_zero]
  have h1 : ∑ j : ZMod p, χ (1 + c * j) = ∑ x : ZMod p, χ x :=
    Equiv.sum_comp ((Equiv.mulLeft₀ c hc).trans (Equiv.addLeft (1 : ZMod p))) (fun x => χ x)
  calc ∑ k : ZMod p, χ (k * (k + c))
      = ∑ k : ZMod p, (χ (1 + c * k⁻¹) - (if k⁻¹ = 0 then (1 : ℤ) else 0)) :=
        Finset.sum_congr rfl (fun k _ => key k)
    _ = ∑ j : ZMod p, (χ (1 + c * j) - (if j = 0 then (1 : ℤ) else 0)) :=
        Equiv.sum_comp (Equiv.inv (ZMod p))
          (fun j => χ (1 + c * j) - (if j = 0 then (1 : ℤ) else 0))
    _ = -1 := by
        rw [Finset.sum_sub_distrib, h1, sum_qc hp2]
        simp

theorem sum_qc_pair (hp2 : p ≠ 2) (a b : ZMod p) :
    ∑ k : ZMod p, χ ((k + a) * (k + b)) = if a = b then (p : ℤ) - 1 else -1 := by
  by_cases h : a = b
  · subst h
    rw [if_pos rfl]
    have h1 : ∑ k : ZMod p, χ ((k + a) * (k + a)) = ∑ x : ZMod p, χ (x * x) :=
      Equiv.sum_comp (Equiv.addRight a) (fun x => χ (x * x))
    have h2 : ∀ x : ZMod p, χ (x * x) = 1 - (if x = 0 then (1 : ℤ) else 0) := by
      intro x
      by_cases hx : x = 0
      · subst hx
        simp
      · rw [if_neg hx, sub_zero, ← pow_two, quadraticChar_sq_one' hx]
    calc ∑ k : ZMod p, χ ((k + a) * (k + a)) = ∑ x : ZMod p, χ (x * x) := h1
      _ = ∑ x : ZMod p, ((1 : ℤ) - (if x = 0 then (1 : ℤ) else 0)) :=
          Finset.sum_congr rfl (fun x _ => h2 x)
      _ = (p : ℤ) - 1 := by
          rw [Finset.sum_sub_distrib]
          simp [ZMod.card]
  · rw [if_neg h]
    have hc : b - a ≠ 0 := sub_ne_zero.mpr (Ne.symm h)
    have h1 : ∑ k : ZMod p, χ ((k + a) * (k + a + (b - a))) =
        ∑ x : ZMod p, χ (x * (x + (b - a))) :=
      Equiv.sum_comp (Equiv.addRight a) (fun x => χ (x * (x + (b - a))))
    have h2 : ∀ k : ZMod p, (k + a) * (k + b) = (k + a) * (k + a + (b - a)) := by
      intro k
      ring
    calc ∑ k : ZMod p, χ ((k + a) * (k + b))
        = ∑ k : ZMod p, χ ((k + a) * (k + a + (b - a))) :=
          Finset.sum_congr rfl (fun k _ => by rw [h2 k])
      _ = ∑ x : ZMod p, χ (x * (x + (b - a))) := h1
      _ = -1 := sum_qc_shift hp2 hc

/-- The sum of Jacobsthal, `S(k) = ∑_v χ(v(k + v²))`. -/
def jac (p : ℕ) [Fact p.Prime] (k : ZMod p) : ℤ :=
  ∑ v : ZMod p, quadraticChar (ZMod p) (v * (k + v ^ 2))

theorem qc_neg (ι : ZMod p) (hι : ι ^ 2 = -1) (x : ZMod p) : χ (-x) = χ x := by
  have hι0 : ι ≠ 0 := by
    rintro rfl
    simp at hι
  have h1 : χ (-1) = 1 := by
    rw [← hι]
    exact quadraticChar_sq_one' hι0
  rw [← neg_one_mul, map_mul, h1, one_mul]

theorem jac_sum_sq (hp2 : p ≠ 2) (ι : ZMod p) (hι : ι ^ 2 = -1) :
    ∑ k : ZMod p, jac p k ^ 2 = 2 * (p : ℤ) * ((p : ℤ) - 1) := by
  have h2ne := zmod_two_ne_zero (p := p) hp2
  have hneg := qc_neg ι hι
  have hS : ∀ k : ZMod p, jac p k ^ 2 =
      ∑ v : ZMod p, ∑ w : ZMod p, χ v * χ w * χ ((k + v ^ 2) * (k + w ^ 2)) := by
    intro k
    unfold jac
    rw [pow_two, Finset.sum_mul_sum]
    refine Finset.sum_congr rfl (fun v _ => Finset.sum_congr rfl (fun w _ => ?_))
    simp only [map_mul]
    ring
  have pt : ∀ v w : ZMod p, χ v * χ w * (if v ^ 2 = w ^ 2 then (p : ℤ) - 1 else -1) =
      (p : ℤ) * (χ v * ((if w = v then χ w else 0) + (if w = -v then χ w else 0))) -
        χ v * χ w := by
    intro v w
    by_cases h : v ^ 2 = w ^ 2
    · rw [if_pos h]
      have hor : w = v ∨ w = -v := by
        have h0 : (w - v) * (w + v) = 0 := by linear_combination -h
        rcases mul_eq_zero.mp h0 with h1 | h1
        · exact Or.inl (sub_eq_zero.mp h1)
        · exact Or.inr (eq_neg_of_add_eq_zero_left h1)
      have key : χ v * ((if w = v then χ w else 0) + (if w = -v then χ w else 0)) =
          χ v * χ w := by
        by_cases h1 : w = v
        · by_cases h2 : w = -v
          · have hv0 : v = 0 := by
              have h4 : v = -v := h1.symm.trans h2
              have h3 : (2 : ZMod p) * v = 0 := by linear_combination h4
              exact (mul_eq_zero.mp h3).resolve_left h2ne
            subst hv0
            simp
          · rw [if_pos h1, if_neg h2, add_zero]
        · rw [if_neg h1, if_pos (hor.resolve_left h1), zero_add]
      rw [key]
      ring
    · have h1 : w ≠ v := fun e => h (by rw [e])
      have h2 : w ≠ -v := fun e => h (by rw [e]; ring)
      rw [if_neg h, if_neg h1, if_neg h2]
      ring
  have inner : ∀ v : ZMod p, ∑ w : ZMod p,
      ((p : ℤ) * (χ v * ((if w = v then χ w else 0) + (if w = -v then χ w else 0))) -
        χ v * χ w) = 2 * (p : ℤ) * χ v ^ 2 - χ v * ∑ w : ZMod p, χ w := by
    intro v
    rw [Finset.sum_sub_distrib, ← Finset.mul_sum, ← Finset.mul_sum, ← Finset.mul_sum,
      Finset.sum_add_distrib]
    simp only [Finset.sum_ite_eq', Finset.mem_univ, if_true]
    rw [hneg v]
    ring
  have h3 : ∑ v : ZMod p, χ v ^ 2 = (p : ℤ) - 1 := by
    have h4 : ∀ x : ZMod p, χ x ^ 2 = 1 - (if x = 0 then (1 : ℤ) else 0) := by
      intro x
      by_cases hx : x = 0
      · subst hx
        simp
      · rw [if_neg hx, sub_zero, quadraticChar_sq_one hx]
    calc ∑ v : ZMod p, χ v ^ 2 = ∑ x : ZMod p, ((1 : ℤ) - (if x = 0 then (1 : ℤ) else 0)) :=
          Finset.sum_congr rfl (fun x _ => h4 x)
      _ = (p : ℤ) - 1 := by
          rw [Finset.sum_sub_distrib]
          simp [ZMod.card]
  calc ∑ k : ZMod p, jac p k ^ 2
      = ∑ k : ZMod p, ∑ v : ZMod p, ∑ w : ZMod p,
          χ v * χ w * χ ((k + v ^ 2) * (k + w ^ 2)) :=
        Finset.sum_congr rfl (fun k _ => hS k)
    _ = ∑ v : ZMod p, ∑ k : ZMod p, ∑ w : ZMod p,
          χ v * χ w * χ ((k + v ^ 2) * (k + w ^ 2)) := Finset.sum_comm
    _ = ∑ v : ZMod p, ∑ w : ZMod p, ∑ k : ZMod p,
          χ v * χ w * χ ((k + v ^ 2) * (k + w ^ 2)) :=
        Finset.sum_congr rfl (fun v _ => Finset.sum_comm)
    _ = ∑ v : ZMod p, ∑ w : ZMod p,
          χ v * χ w * (if v ^ 2 = w ^ 2 then (p : ℤ) - 1 else -1) := by
        refine Finset.sum_congr rfl (fun v _ => Finset.sum_congr rfl (fun w _ => ?_))
        rw [← Finset.mul_sum, sum_qc_pair hp2]
    _ = ∑ v : ZMod p, (2 * (p : ℤ) * χ v ^ 2 - χ v * ∑ w : ZMod p, χ w) :=
        Finset.sum_congr rfl
          (fun v _ => (Finset.sum_congr rfl (fun w _ => pt v w)).trans (inner v))
    _ = 2 * (p : ℤ) * ((p : ℤ) - 1) := by
        rw [sum_qc hp2]
        simp only [mul_zero, sub_zero]
        rw [← Finset.mul_sum, h3]

theorem jac_sq (t : ZMod p) (ht : t ≠ 0) : jac p (t ^ 2) ^ 2 = jac p 1 ^ 2 := by
  have h1 : jac p (t ^ 2) = χ t * jac p 1 := by
    unfold jac
    rw [Finset.mul_sum]
    calc ∑ v : ZMod p, χ (v * (t ^ 2 + v ^ 2))
        = ∑ v : ZMod p, χ ((t * v) * (t ^ 2 + (t * v) ^ 2)) :=
          (Equiv.sum_comp (Equiv.mulLeft₀ t ht) (fun v => χ (v * (t ^ 2 + v ^ 2)))).symm
      _ = ∑ v : ZMod p, χ t * χ (v * (1 + v ^ 2)) := by
          refine Finset.sum_congr rfl (fun v _ => ?_)
          have e : t * v * (t ^ 2 + (t * v) ^ 2) = t ^ 2 * (t * (v * (1 + v ^ 2))) := by ring
          rw [e, map_mul, map_mul, quadraticChar_sq_one' ht, one_mul]
  rw [h1, mul_pow, quadraticChar_sq_one ht, one_mul]

theorem qc_le_one (x : ZMod p) : χ x ≤ 1 ∧ -1 ≤ χ x := by
  by_cases hx : x = 0
  · subst hx
    simp
  · rcases quadraticChar_dichotomy hx with h | h <;> rw [h] <;> norm_num

/-- **The point of the conic.** For a prime `p ≡ 1 (mod 4)`, `p ≥ 37`, and each sign `c`, there is
`v` modulo `p` with `2(1 + v²)` a nonzero square and `χ(2v) = c`. -/
theorem exists_conic_point (h14 : p % 4 = 1) (hbig : 37 ≤ p) {c : ℤ} (hc : c = 1 ∨ c = -1) :
    ∃ v : ZMod p, χ (2 * (1 + v ^ 2)) = 1 ∧ χ (2 * v) = c := by
  have hp2 : p ≠ 2 := by omega
  obtain ⟨ι, hι'⟩ := (ZMod.exists_sq_eq_neg_one_iff (p := p)).mpr (by omega)
  have hι : ι ^ 2 = -1 := by
    rw [pow_two]
    exact hι'.symm
  have h2ne := zmod_two_ne_zero (p := p) hp2
  have hιne : ι ≠ -ι := by
    intro h
    have h3 : (2 : ZMod p) * ι = 0 := by linear_combination h
    have h4 : ι = 0 := (mul_eq_zero.mp h3).resolve_left h2ne
    rw [h4] at hι
    simp at hι
  by_contra hno
  simp only [not_exists, not_and] at hno
  have hsum1 : ∑ v : ZMod p, χ (1 + v ^ 2) = -1 := by
    have h := sum_qc_pair hp2 ι (-ι)
    rw [if_neg hιne] at h
    rw [← h]
    refine Finset.sum_congr rfl (fun v _ => ?_)
    have e : 1 + v ^ 2 = (v + ι) * (v + -ι) := by linear_combination hι
    rw [e]
  have hcard : ∑ _v : ZMod p, (1 : ℤ) = p := by simp [ZMod.card]
  have ha : ∑ v : ZMod p, (1 + χ (2 * (1 + v ^ 2))) * (1 + c * χ (2 * v)) =
      (p : ℤ) - χ 2 + c * jac p 1 := by
    have e : ∀ v : ZMod p, (1 + χ (2 * (1 + v ^ 2))) * (1 + c * χ (2 * v)) =
        1 + χ 2 * χ (1 + v ^ 2) + c * (χ 2 * χ v) + c * χ (v * (1 + v ^ 2)) := by
      intro v
      have e1 : χ (2 * (1 + v ^ 2)) * χ (2 * v) = χ (v * (1 + v ^ 2)) := by
        rw [← map_mul]
        have e2 : 2 * (1 + v ^ 2) * (2 * v) = 2 ^ 2 * (v * (1 + v ^ 2)) := by ring
        rw [e2, map_mul, quadraticChar_sq_one' h2ne, one_mul]
      rw [← e1, map_mul, map_mul]
      ring
    calc ∑ v : ZMod p, (1 + χ (2 * (1 + v ^ 2))) * (1 + c * χ (2 * v))
        = ∑ v : ZMod p,
            (1 + χ 2 * χ (1 + v ^ 2) + c * (χ 2 * χ v) + c * χ (v * (1 + v ^ 2))) :=
          Finset.sum_congr rfl (fun v _ => e v)
      _ = ∑ _v : ZMod p, (1 : ℤ) + χ 2 * ∑ v : ZMod p, χ (1 + v ^ 2) +
            c * (χ 2 * ∑ v : ZMod p, χ v) + c * ∑ v : ZMod p, χ (v * (1 + v ^ 2)) := by
          simp only [Finset.sum_add_distrib, Finset.mul_sum]
      _ = (p : ℤ) - χ 2 + c * jac p 1 := by
          rw [hsum1, sum_qc hp2, hcard]
          unfold jac
          ring
  have hb : ∑ v : ZMod p, (1 + χ (2 * (1 + v ^ 2))) * (1 + c * χ (2 * v)) ≤ 6 := by
    have hpt : ∀ v : ZMod p, (1 + χ (2 * (1 + v ^ 2))) * (1 + c * χ (2 * v)) ≤
        if v ∈ ({0, ι, -ι} : Finset (ZMod p)) then (2 : ℤ) else 0 := by
      intro v
      have ha' := qc_le_one (p := p) (2 * (1 + v ^ 2))
      have hb' := qc_le_one (p := p) (2 * v)
      by_cases hv : v ∈ ({0, ι, -ι} : Finset (ZMod p))
      · rw [if_pos hv]
        have hz : χ (2 * (1 + v ^ 2)) = 0 ∨ χ (2 * v) = 0 := by
          simp only [Finset.mem_insert, Finset.mem_singleton] at hv
          rcases hv with h | h | h
          · right
            rw [h]
            simp
          · left
            have e : (2 : ZMod p) * (1 + v ^ 2) = 0 := by rw [h, hι]; ring
            rw [e]
            simp
          · left
            have e : (2 : ZMod p) * (1 + v ^ 2) = 0 := by rw [h, neg_sq, hι]; ring
            rw [e]
            simp
        rcases hz with h | h
        · rw [h]
          rcases hc with rfl | rfl <;> nlinarith [hb'.1, hb'.2]
        · rw [h]
          rcases hc with rfl | rfl <;> nlinarith [ha'.1, ha'.2]
      · rw [if_neg hv]
        simp only [Finset.mem_insert, Finset.mem_singleton, not_or] at hv
        obtain ⟨hv0, hv1, hv2⟩ := hv
        have hne1 : (2 : ZMod p) * (1 + v ^ 2) ≠ 0 := by
          refine mul_ne_zero h2ne (fun h => ?_)
          have h0 : (v - ι) * (v + ι) = 0 := by linear_combination h - hι
          rcases mul_eq_zero.mp h0 with h1 | h1
          · exact hv1 (sub_eq_zero.mp h1)
          · exact hv2 (eq_neg_of_add_eq_zero_left h1)
        have hne2 : (2 : ZMod p) * v ≠ 0 := mul_ne_zero h2ne hv0
        rcases quadraticChar_dichotomy hne1 with h1 | h1
        · have h2 := hno v h1
          rcases quadraticChar_dichotomy hne2 with h3 | h3
          · rcases hc with rfl | rfl
            · exact absurd h3 h2
            · rw [h1, h3]
              norm_num
          · rcases hc with rfl | rfl
            · rw [h1, h3]
              norm_num
            · exact absurd h3 h2
        · rw [h1]
          simp
    calc ∑ v : ZMod p, (1 + χ (2 * (1 + v ^ 2))) * (1 + c * χ (2 * v))
        ≤ ∑ v : ZMod p, (if v ∈ ({0, ι, -ι} : Finset (ZMod p)) then (2 : ℤ) else 0) :=
          Finset.sum_le_sum (fun v _ => hpt v)
      _ = ∑ _v ∈ ({0, ι, -ι} : Finset (ZMod p)), (2 : ℤ) := by
          rw [Finset.sum_ite_mem, Finset.univ_inter]
      _ ≤ 6 := by
          rw [Finset.sum_const, nsmul_eq_mul]
          have h3 := Finset.card_le_three (a := (0 : ZMod p)) (b := ι) (c := -ι)
          have h4 : ((({0, ι, -ι} : Finset (ZMod p)).card : ℕ) : ℤ) ≤ 3 := by exact_mod_cast h3
          linarith
  have hχ2 := (qc_le_one (p := p) 2).1
  have hS : c * jac p 1 ≤ 7 - (p : ℤ) := by linarith
  have hp37 : (37 : ℤ) ≤ p := by exact_mod_cast hbig
  have hsq : ((p : ℤ) - 7) ^ 2 ≤ jac p 1 ^ 2 := by
    rcases hc with rfl | rfl
    · nlinarith [mul_nonneg (by linarith : (0 : ℤ) ≤ 7 - (p : ℤ) - jac p 1)
        (by linarith : (0 : ℤ) ≤ (p : ℤ) - 7 - jac p 1)]
    · nlinarith [mul_nonneg (by linarith : (0 : ℤ) ≤ jac p 1 - ((p : ℤ) - 7))
        (by linarith : (0 : ℤ) ≤ jac p 1 + ((p : ℤ) - 7))]
  have h3 : 3 * jac p 1 ^ 2 ≤ 2 * (p : ℤ) * ((p : ℤ) - 1) := by
    rw [← jac_sum_sq hp2 ι hι]
    have n3 : (3 : ZMod p) ≠ 0 := by
      have h := zmod_intCast_ne_zero (p := p) (n := 3) (by norm_num) (by omega)
      simpa using h
    have n5 : (5 : ZMod p) ≠ 0 := by
      have h := zmod_intCast_ne_zero (p := p) (n := 5) (by norm_num) (by omega)
      simpa using h
    have n8 : (8 : ZMod p) ≠ 0 := by
      have h := zmod_intCast_ne_zero (p := p) (n := 8) (by norm_num) (by omega)
      simpa using h
    have d14 : (1 : ZMod p) ≠ 4 := fun h => n3 (by linear_combination -h)
    have d19 : (1 : ZMod p) ≠ 9 := fun h => n8 (by linear_combination -h)
    have d49 : (4 : ZMod p) ≠ 9 := fun h => n5 (by linear_combination -h)
    have e4 : jac p 4 ^ 2 = jac p 1 ^ 2 := by
      have e : (4 : ZMod p) = 2 ^ 2 := by norm_num
      rw [e]
      exact jac_sq 2 h2ne
    have e9 : jac p 9 ^ 2 = jac p 1 ^ 2 := by
      have e : (9 : ZMod p) = 3 ^ 2 := by norm_num
      rw [e]
      exact jac_sq 3 n3
    have hsub : ∑ k ∈ ({1, 4, 9} : Finset (ZMod p)), jac p k ^ 2 ≤ ∑ k : ZMod p, jac p k ^ 2 :=
      Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ _) (fun k _ _ => sq_nonneg _)
    have heq : ∑ k ∈ ({1, 4, 9} : Finset (ZMod p)), jac p k ^ 2 = 3 * jac p 1 ^ 2 := by
      rw [Finset.sum_insert (by
          rw [Finset.mem_insert, Finset.mem_singleton]
          exact fun h => h.elim d14 d19),
        Finset.sum_insert (by
          rw [Finset.mem_singleton]
          exact d49), Finset.sum_singleton, e4, e9]
      ring
    linarith
  nlinarith [hsq, h3, mul_nonneg (by linarith : (0 : ℤ) ≤ (p : ℤ) - 37)
    (by linarith : (0 : ℤ) ≤ (p : ℤ) - 3)]

theorem exists_int_sq_of_isSquare {c : ℤ} (hc : IsSquare ((c : ℤ) : ZMod p))
    (hc0 : ((c : ℤ) : ZMod p) ≠ 0) : ∃ t : ℤ, t ^ 2 ≡ c [ZMOD p] ∧ ¬(p : ℤ) ∣ t := by
  obtain ⟨r, hr⟩ := hc
  obtain ⟨t, rfl⟩ := ZMod.intCast_surjective r
  refine ⟨t, ?_, ?_⟩
  · rw [← ZMod.intCast_eq_intCast_iff]
    push_cast
    rw [hr]
    ring
  · intro h
    apply hc0
    rw [hr, (ZMod.intCast_zmod_eq_zero_iff_dvd t p).mpr h]
    ring

/-- From the two points of the conic to the unit classes: if for each sign `c` there is `v` with
`2(1 + v²)` a nonzero square and `χ(2v) = c`, every `d` prime to `p` passes at some `(1, v)`. -/
theorem unit_passes_of_conic (h14 : p % 4 = 1)
    (hcon : ∀ {c : ℤ}, c = 1 ∨ c = -1 →
      ∃ v : ZMod p, χ (2 * (1 + v ^ 2)) = 1 ∧ χ (2 * v) = c) {m n d : ℤ}
    (hm : (p : ℤ) ∣ m) (hn : ¬(p : ℤ) ∣ n) (hd : ¬(p : ℤ) ∣ d) (k : ℕ) :
    ∃ v : ℤ, ∀ f : Bool × Bool × Bool, GSqM ((p : ℤ) ^ (k + 1)) d (Q m n 1 v f) := by
  have hp : p.Prime := Fact.out
  have hp2 : p ≠ 2 := by omega
  have hp' : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  have h2 := not_dvd_two_of_odd hp hp2
  have h1 : ¬(p : ℤ) ∣ 1 := fun h => hp'.not_isUnit (isUnit_of_dvd_one h)
  have h2ne := zmod_two_ne_zero (p := p) hp2
  obtain ⟨ι, hι'⟩ := (ZMod.exists_sq_eq_neg_one_iff (p := p)).mpr (by omega)
  have hι : ι ^ 2 = -1 := by
    rw [pow_two]
    exact hι'.symm
  have hd0 : ((d : ℤ) : ZMod p) ≠ 0 := by
    rw [Ne, ZMod.intCast_zmod_eq_zero_iff_dvd]
    exact hd
  obtain ⟨v1, hv1', hv2'⟩ := hcon (quadraticChar_dichotomy hd0)
  obtain ⟨v0, hv1, hv2, hv3⟩ : ∃ v0 : ZMod p, χ (2 * (1 + v0 ^ 2)) = 1 ∧
      χ (2 * v0) = χ ((d : ℤ) : ZMod p) ∧ 1 + v0 ≠ 0 := by
    by_cases h : 1 + v1 = 0
    · refine ⟨1, ?_, ?_, ?_⟩
      · have e : (2 : ZMod p) * (1 + 1 ^ 2) = 2 ^ 2 := by ring
        rw [e]
        exact quadraticChar_sq_one' h2ne
      · have e : (2 : ZMod p) * v1 = -(2 * 1) := by linear_combination 2 * h
        rw [← hv2', e, qc_neg ι hι]
      · exact fun h' => h2ne (by linear_combination h')
    · exact ⟨v1, hv1', hv2', h⟩
  obtain ⟨v, rfl⟩ := ZMod.intCast_surjective v0
  have hv0 : ¬(p : ℤ) ∣ v := by
    intro h
    have e : ((v : ℤ) : ZMod p) = 0 := (ZMod.intCast_zmod_eq_zero_iff_dvd v p).mpr h
    rw [e, mul_zero, MulChar.map_zero] at hv2
    rcases quadraticChar_dichotomy hd0 with h5 | h5 <;> rw [h5] at hv2 <;> norm_num at hv2
  have huv : ¬(p : ℤ) ∣ 1 + v := by
    intro h
    apply hv3
    have e := (ZMod.intCast_zmod_eq_zero_iff_dvd (1 + v) p).mpr h
    push_cast at e
    exact e
  have hne : (2 : ZMod p) * (1 + ((v : ℤ) : ZMod p) ^ 2) ≠ 0 := by
    intro h
    rw [h, MulChar.map_zero] at hv1
    norm_num at hv1
  have ec : (((2 * (1 ^ 2 + v ^ 2) : ℤ)) : ZMod p) = 2 * (1 + ((v : ℤ) : ZMod p) ^ 2) := by
    push_cast
    ring
  obtain ⟨s, hs, hs0⟩ := exists_int_sq_of_isSquare (p := p) (c := 2 * (1 ^ 2 + v ^ 2))
    (by rw [ec]; exact (quadraticChar_one_iff_isSquare hne).mp hv1) (by rw [ec]; exact hne)
  have hq : χ (((d : ℤ) : ZMod p) * (2 * ((v : ℤ) : ZMod p))) = 1 := by
    rw [map_mul, hv2, ← pow_two, quadraticChar_sq_one hd0]
  have hne' : ((d : ℤ) : ZMod p) * (2 * ((v : ℤ) : ZMod p)) ≠ 0 := by
    intro h
    rw [h, MulChar.map_zero] at hq
    norm_num at hq
  have ed : (((d * (2 * 1 * v) : ℤ)) : ZMod p) = ((d : ℤ) : ZMod p) * (2 * ((v : ℤ) : ZMod p)) := by
    push_cast
    ring
  obtain ⟨t, ht, ht0⟩ := exists_int_sq_of_isSquare (p := p) (c := d * (2 * 1 * v))
    (by rw [ed]; exact (quadraticChar_one_iff_isSquare hne').mp hq) (by rw [ed]; exact hne')
  obtain ⟨_, _, _, _, _, _, _, _, g⟩ := pass_first_order hp hp2 (u := 1) (v := v) (s := s)
    (t := t) hm hn h1 hv0 huv hs hs0 ht ht0 k
  exact ⟨v, g⟩

/-- **Every unit class passes at a prime `p ≡ 1 (mod 4)`, `p ≥ 37`, dividing `m`**: a `d` not
divisible by `p` passes at some point `(1, v)` modulo every power of `p`. This is section 4.2 of
`COMPLETE.md` with `pass_first_order`. -/
theorem unit_passes_of_one_mod_four (h14 : p % 4 = 1) (hbig : 37 ≤ p) {m n d : ℤ}
    (hm : (p : ℤ) ∣ m) (hn : ¬(p : ℤ) ∣ n) (hd : ¬(p : ℤ) ∣ d) (k : ℕ) :
    ∃ v : ℤ, ∀ f : Bool × Bool × Bool, GSqM ((p : ℤ) ^ (k + 1)) d (Q m n 1 v f) :=
  unit_passes_of_conic h14 (fun hc => exists_conic_point h14 hbig hc) hm hn hd k

end Jacobsthal

/-! ### The primes 13, 17, 29, and the statement for every prime `≡ 1 (mod 4)` except 5 -/

local instance fact13 : Fact (Nat.Prime 13) := ⟨by norm_num⟩
local instance fact17 : Fact (Nat.Prime 17) := ⟨by norm_num⟩
local instance fact29 : Fact (Nat.Prime 29) := ⟨by norm_num⟩

theorem conic13 {c : ℤ} (hc : c = 1 ∨ c = -1) : ∃ v : ZMod 13,
    quadraticChar (ZMod 13) (2 * (1 + v ^ 2)) = 1 ∧ quadraticChar (ZMod 13) (2 * v) = c := by
  rcases hc with rfl | rfl
  · exact ⟨2, (quadraticChar_one_iff_isSquare (by decide)).mpr ⟨6, by decide⟩,
      (quadraticChar_one_iff_isSquare (by decide)).mpr ⟨2, by decide⟩⟩
  · refine ⟨1, (quadraticChar_one_iff_isSquare (by decide)).mpr ⟨2, by decide⟩,
      quadraticChar_neg_one_iff_not_isSquare.mpr ?_⟩
    rintro ⟨r, hr⟩
    revert r
    decide

theorem conic17 {c : ℤ} (hc : c = 1 ∨ c = -1) : ∃ v : ZMod 17,
    quadraticChar (ZMod 17) (2 * (1 + v ^ 2)) = 1 ∧ quadraticChar (ZMod 17) (2 * v) = c := by
  rcases hc with rfl | rfl
  · exact ⟨1, (quadraticChar_one_iff_isSquare (by decide)).mpr ⟨2, by decide⟩,
      (quadraticChar_one_iff_isSquare (by decide)).mpr ⟨6, by decide⟩⟩
  · refine ⟨5, (quadraticChar_one_iff_isSquare (by decide)).mpr ⟨1, by decide⟩,
      quadraticChar_neg_one_iff_not_isSquare.mpr ?_⟩
    rintro ⟨r, hr⟩
    revert r
    decide

theorem conic29 {c : ℤ} (hc : c = 1 ∨ c = -1) : ∃ v : ZMod 29,
    quadraticChar (ZMod 29) (2 * (1 + v ^ 2)) = 1 ∧ quadraticChar (ZMod 29) (2 * v) = c := by
  rcases hc with rfl | rfl
  · exact ⟨3, (quadraticChar_one_iff_isSquare (by decide)).mpr ⟨7, by decide⟩,
      (quadraticChar_one_iff_isSquare (by decide)).mpr ⟨8, by decide⟩⟩
  · refine ⟨1, (quadraticChar_one_iff_isSquare (by decide)).mpr ⟨2, by decide⟩,
      quadraticChar_neg_one_iff_not_isSquare.mpr ?_⟩
    rintro ⟨r, hr⟩
    revert r
    decide

/-- **Every unit class passes at every prime `p ≡ 1 (mod 4)` dividing `m`, except `p = 5`**: a `d`
not divisible by `p` passes at some point `(1, v)` modulo every power of `p`. For `p = 5` this is
false when 25 does not divide `m` (the law at 5, `class_mod_five`). With
`unit_passes_of_three_mod_four`: the unit classes pass at every odd prime of `m` except 5. -/
theorem unit_passes_of_one_mod_four_ne_five {p : ℕ} (hp : p.Prime) (h14 : p % 4 = 1)
    (h5 : p ≠ 5) {m n d : ℤ} (hm : (p : ℤ) ∣ m) (hn : ¬(p : ℤ) ∣ n) (hd : ¬(p : ℤ) ∣ d)
    (k : ℕ) :
    ∃ v : ℤ, ∀ f : Bool × Bool × Bool, GSqM ((p : ℤ) ^ (k + 1)) d (Q m n 1 v f) := by
  by_cases hbig : 37 ≤ p
  · have : Fact p.Prime := ⟨hp⟩
    exact unit_passes_of_one_mod_four h14 hbig hm hn hd k
  · have h : p = 13 ∨ p = 17 ∨ p = 29 := by
      have hlt : p < 37 := by omega
      have h2 : 2 ≤ p := hp.two_le
      interval_cases p <;> first | omega | exact absurd hp (by decide)
    rcases h with rfl | rfl | rfl
    · exact unit_passes_of_conic h14 conic13 hm hn hd k
    · exact unit_passes_of_conic h14 conic17 hm hn hd k
    · exact unit_passes_of_conic h14 conic29 hm hn hd k

end MagicSquareOfSquares.CoverDescent
