import LeanRank.LemmaDBridge

/-!
# Arithmetic of `ℤ[√2]` for the descent over `ℚ(√2)`

Written by an AI model (Claude, made by Anthropic) on 26 November 2026 (UTC); not read by a person.
Experiment, not counted in the census.

* `ℤ√2` (Mathlib's `Zsqrtd 2`) is a Euclidean domain (norm-Euclidean, rounding as in Mathlib's
  `GaussianInt`), hence a PID and a UFD;
* `pi2 = √2` and `three = 3` are prime, with `pi2 ∣ z ↔ 2 ∣ z.re`;
* every unit is `± v²` or `± ε v²`, `ε = 1 + √2` (`unit_classes`, from Mathlib's Pell theory with the
  fundamental solution `3 + 2√2 = ε²`);
* the ring hom `iota2 : ℤ√2 → ℚ(√2)` (`Q2` of `LemmaDBridge.lean`) is injective, and every element of `Q2`
  is `iota2 p / iota2 q` (`Q2_frac`).
-/

namespace MagicSquareOfSquares.Zsqrt2

open Zsqrtd MagicSquareOfSquares.Family

theorem nonsq : ∀ n : ℤ, (2 : ℤ) ≠ n * n := by
  intro n hn
  have h1 : n < 2 := by nlinarith
  have h2 : -2 < n := by nlinarith
  interval_cases n <;> omega

theorem norm_ne_zero {z : ℤ√2} (hz : z ≠ 0) : z.norm ≠ 0 :=
  fun h => hz ((norm_eq_zero nonsq z).mp h)

/-! ### Euclidean division -/

instance : Div (ℤ√2) :=
  ⟨fun x y => ⟨round (((x * star y).re : ℚ) / (y.norm : ℚ)), round (((x * star y).im : ℚ) / (y.norm : ℚ))⟩⟩

instance : Mod (ℤ√2) := ⟨fun x y => x - y * (x / y)⟩

theorem div_def (x y : ℤ√2) : x / y =
    ⟨round (((x * star y).re : ℚ) / (y.norm : ℚ)), round (((x * star y).im : ℚ) / (y.norm : ℚ))⟩ := rfl

theorem mod_def (x y : ℤ√2) : x % y = x - y * (x / y) := rfl

theorem round_bound (P n : ℤ) (hn : n ≠ 0) : (2 * (P - n * round ((P : ℚ) / n))) ^ 2 ≤ n ^ 2 := by
  have hn' : (n : ℚ) ≠ 0 := by exact_mod_cast hn
  have h := abs_le.mp (abs_sub_round ((P : ℚ) / n))
  have key : (((2 * (P - n * round ((P : ℚ) / n))) ^ 2 : ℤ) : ℚ)
      = 4 * (n : ℚ) ^ 2 * ((P : ℚ) / n - round ((P : ℚ) / n)) ^ 2 := by
    push_cast; field_simp; ring
  have h2 : ((P : ℚ) / n - round ((P : ℚ) / n)) ^ 2 ≤ 1 / 4 := by
    nlinarith [mul_nonneg (sub_nonneg.mpr h.1) (sub_nonneg.mpr h.2)]
  have : (((2 * (P - n * round ((P : ℚ) / n))) ^ 2 : ℤ) : ℚ) ≤ ((n ^ 2 : ℤ) : ℚ) := by
    rw [key]; push_cast; nlinarith [sq_nonneg (n : ℚ)]
  exact_mod_cast this

theorem norm_mod_lt (x : ℤ√2) {y : ℤ√2} (hy : y ≠ 0) : (x % y).norm.natAbs < y.norm.natAbs := by
  have hn := norm_ne_zero hy
  have hrel : (x % y) * star y = x * star y - (y.norm : ℤ√2) * (x / y) := by
    rw [mod_def, norm_eq_mul_conj]; ring
  have hq1 : (x / y).re = round (((x * star y).re : ℚ) / (y.norm : ℚ)) := rfl
  have hq2 : (x / y).im = round (((x * star y).im : ℚ) / (y.norm : ℚ)) := rfl
  have hnorm : (x % y).norm * y.norm = ((x * star y).re - y.norm * (x / y).re) ^ 2
      - 2 * ((x * star y).im - y.norm * (x / y).im) ^ 2 := by
    have := congrArg Zsqrtd.norm hrel
    rw [norm_mul, norm_conj] at this
    rw [this, norm_def]
    simp only [re_sub, im_sub, re_smul, im_smul]
    ring
  have hb1 := round_bound (x * star y).re y.norm hn
  have hb2 := round_bound (x * star y).im y.norm hn
  rw [← hq1] at hb1
  rw [← hq2] at hb2
  have h3 : |2 * ((x % y).norm * y.norm)| ≤ y.norm ^ 2 := by
    rw [abs_le, hnorm]
    constructor
    · nlinarith [sq_nonneg ((x * star y).re - y.norm * (x / y).re),
        sq_nonneg ((x * star y).im - y.norm * (x / y).im)]
    · nlinarith [sq_nonneg ((x * star y).re - y.norm * (x / y).re),
        sq_nonneg ((x * star y).im - y.norm * (x / y).im)]
  rw [abs_mul, abs_mul, abs_two, ← sq_abs y.norm] at h3
  have hpos : 0 < |y.norm| := abs_pos.mpr hn
  have h4 : 2 * |(x % y).norm| ≤ |y.norm| := by
    have h5 : (2 * |(x % y).norm|) * |y.norm| ≤ |y.norm| * |y.norm| := by nlinarith
    exact le_of_mul_le_mul_right h5 hpos
  have h6 : |(x % y).norm| < |y.norm| := by
    have := abs_nonneg (x % y).norm
    linarith
  have h7 : ((x % y).norm.natAbs : ℤ) < (y.norm.natAbs : ℤ) := by
    rw [Int.natCast_natAbs, Int.natCast_natAbs]; exact h6
  exact_mod_cast h7

theorem norm_le_norm_mul_left (x : ℤ√2) {y : ℤ√2} (hy : y ≠ 0) :
    x.norm.natAbs ≤ (x * y).norm.natAbs := by
  rw [norm_mul, Int.natAbs_mul]
  exact Nat.le_mul_of_pos_right _ (Int.natAbs_pos.mpr (norm_ne_zero hy))

instance : EuclideanDomain (ℤ√2) :=
  { (inferInstance : CommRing (ℤ√2)), (inferInstance : Nontrivial (ℤ√2)) with
    quotient := (· / ·)
    remainder := (· % ·)
    quotient_zero := fun a => by
      rw [div_def]
      ext <;> simp
    quotient_mul_add_remainder_eq := fun x y => by rw [mod_def]; ring
    r := _
    r_wellFounded := (measure (Int.natAbs ∘ Zsqrtd.norm)).wf
    remainder_lt := fun x _ hy => norm_mod_lt x hy
    mul_left_not_lt := fun a _ hb => not_lt_of_ge (norm_le_norm_mul_left a hb) }

theorem natAbs_norm_ge_two {c : ℤ√2} (hc0 : c ≠ 0) (hc : ¬IsUnit c) : 2 ≤ c.norm.natAbs := by
  have h1 : c.norm.natAbs ≠ 1 := fun h => hc (norm_eq_one_iff.mp h)
  have h0 : c.norm.natAbs ≠ 0 := by
    intro h; exact norm_ne_zero hc0 (Int.natAbs_eq_zero.mp h)
  omega

/-! ### The primes `√2` and `3` -/

/-- `√2`. -/
def pi2 : ℤ√2 := ⟨0, 1⟩

/-- `3`. -/
def three : ℤ√2 := ⟨3, 0⟩

/-- `ε = 1 + √2`. -/
def eps : ℤ√2 := ⟨1, 1⟩

/-- `ε⁻¹ = −1 + √2`. -/
def epsinv : ℤ√2 := ⟨-1, 1⟩

theorem eps_mul_epsinv : eps * epsinv = 1 := by decide

theorem pi2_sq : pi2 ^ 2 = 2 := by decide

theorem three_eq : three = 3 := by decide

theorem pi2_mul (z : ℤ√2) : pi2 * z = ⟨2 * z.im, z.re⟩ := by
  ext <;> simp [pi2]

theorem pi2_dvd_iff (z : ℤ√2) : pi2 ∣ z ↔ 2 ∣ z.re := by
  constructor
  · rintro ⟨w, rfl⟩
    rw [pi2_mul]
    exact ⟨w.im, rfl⟩
  · rintro ⟨k, hk⟩
    refine ⟨⟨z.im, k⟩, ?_⟩
    rw [pi2_mul]
    ext
    · exact hk
    · rfl

theorem pi2_ne_zero : pi2 ≠ 0 := by decide

theorem pi2_not_unit : ¬IsUnit pi2 := by
  intro h
  have := norm_eq_one_iff.mpr h
  revert this
  decide

theorem pi2_prime : Prime pi2 := by
  refine ⟨pi2_ne_zero, pi2_not_unit, ?_⟩
  intro a b h
  rw [pi2_dvd_iff, re_mul] at h
  rw [pi2_dvd_iff, pi2_dvd_iff]
  have : (2 : ℤ) ∣ a.re * b.re := by
    rcases h with ⟨k, hk⟩
    exact ⟨k - a.im * b.im, by linarith⟩
  exact Int.prime_two.dvd_mul.mp this

theorem three_of_norm {a b : ℤ} (h : (3 : ℤ) ∣ a * a - 2 * b * b) : (3 : ℤ) ∣ a ∧ (3 : ℤ) ∣ b := by
  have key : ∀ x y : ZMod 3, x * x - 2 * y * y = 0 → x = 0 ∧ y = 0 := by decide
  have h' : ((a * a - 2 * b * b : ℤ) : ZMod 3) = 0 :=
    (ZMod.intCast_zmod_eq_zero_iff_dvd _ 3).mpr (by exact_mod_cast h)
  push_cast at h'
  obtain ⟨ha, hb⟩ := key _ _ h'
  exact ⟨by exact_mod_cast (ZMod.intCast_zmod_eq_zero_iff_dvd a 3).mp ha,
    by exact_mod_cast (ZMod.intCast_zmod_eq_zero_iff_dvd b 3).mp hb⟩

theorem three_dvd_iff (z : ℤ√2) : three ∣ z ↔ (3 : ℤ) ∣ z.re ∧ (3 : ℤ) ∣ z.im := by
  have : three = ((3 : ℤ) : ℤ√2) := by decide
  rw [this, intCast_dvd]

theorem three_ne_zero : three ≠ 0 := by decide

theorem three_not_unit : ¬IsUnit three := by
  intro h
  have := norm_eq_one_iff.mpr h
  revert this
  decide

theorem three_prime : Prime three := by
  refine ⟨three_ne_zero, three_not_unit, ?_⟩
  intro a b h
  rw [three_dvd_iff] at h
  have hN : (3 : ℤ) ∣ (a * b).norm := by
    rw [norm_def]
    obtain ⟨h1, h2⟩ := h
    exact dvd_sub (dvd_mul_of_dvd_left h1 _) (dvd_mul_of_dvd_right h2 _)
  rw [norm_mul] at hN
  rcases Int.prime_three.dvd_mul.mp hN with h3 | h3
  · left; rw [three_dvd_iff]; rw [norm_def] at h3; exact three_of_norm h3
  · right; rw [three_dvd_iff]; rw [norm_def] at h3; exact three_of_norm h3

/-! ### Units: `±v²` or `±εv²` -/

/-- The fundamental solution `3 + 2√2 = ε²` of `x² − 2y² = 1`. -/
def fund : Pell.Solution₁ 2 := Pell.Solution₁.mk 3 2 (by norm_num)

theorem fund_coe : ((fund : Pell.Solution₁ 2) : ℤ√2) = eps ^ 2 := by
  rw [fund, Pell.Solution₁.coe_mk]; decide

theorem fund_isFundamental : Pell.IsFundamental fund := by
  refine ⟨by simp [fund], by simp [fund], ?_⟩
  intro b hb
  simp only [fund, Pell.Solution₁.x_mk]
  by_contra h
  push_neg at h
  have hx : b.x = 2 := by omega
  have := b.prop
  rw [hx] at this
  have : 2 * b.y ^ 2 = 3 := by linarith
  omega

theorem coe_sol_mul (a b : Pell.Solution₁ 2) : ((a * b : Pell.Solution₁ 2) : ℤ√2) = a * b := rfl

theorem fund_zpow_sq (n : ℤ) : ∃ v : ℤ√2, ((fund ^ n : Pell.Solution₁ 2) : ℤ√2) = v ^ 2 := by
  obtain ⟨m, hm | hm⟩ := Int.even_or_odd' n
  · refine ⟨((fund ^ m : Pell.Solution₁ 2) : ℤ√2), ?_⟩
    rw [hm, mul_comm, zpow_mul, zpow_two, coe_sol_mul, sq]
  · refine ⟨((fund ^ m : Pell.Solution₁ 2) : ℤ√2) * eps, ?_⟩
    rw [hm, zpow_add, zpow_one, mul_comm 2 m, zpow_mul, zpow_two, coe_sol_mul, coe_sol_mul,
      fund_coe]
    ring

theorem sol_sq (a : Pell.Solution₁ 2) : ∃ v : ℤ√2, (a : ℤ√2) = v ^ 2 ∨ (a : ℤ√2) = -v ^ 2 := by
  obtain ⟨n, hn | hn⟩ := fund_isFundamental.eq_zpow_or_neg_zpow a
  · obtain ⟨v, hv⟩ := fund_zpow_sq n
    exact ⟨v, Or.inl (by rw [hn, hv])⟩
  · obtain ⟨v, hv⟩ := fund_zpow_sq n
    refine ⟨v, Or.inr ?_⟩
    rw [hn, ← hv]
    rfl

/-- **Units of `ℤ√2` modulo squares**: `±1`, `±ε`. -/
theorem unit_classes {u : ℤ√2} (hu : IsUnit u) :
    ∃ v : ℤ√2, u = v ^ 2 ∨ u = -v ^ 2 ∨ u = eps * v ^ 2 ∨ u = -(eps * v ^ 2) := by
  have h1 : u.norm.natAbs = 1 := norm_eq_one_iff.mpr hu
  rcases Int.natAbs_eq_iff.mp h1 with h | h
  · obtain ⟨v, hv | hv⟩ := sol_sq ⟨u, norm_eq_one_iff_mem_unitary.mp (by exact_mod_cast h)⟩
    · exact ⟨v, Or.inl hv⟩
    · exact ⟨v, Or.inr (Or.inl hv)⟩
  · have hn : (u * eps).norm = 1 := by
      rw [norm_mul, h]; decide
    obtain ⟨v, hv | hv⟩ := sol_sq ⟨u * eps, norm_eq_one_iff_mem_unitary.mp hn⟩
    · have hv' : u * eps = v ^ 2 := hv
      refine ⟨v * epsinv, Or.inr (Or.inr (Or.inl ?_))⟩
      calc u = u * eps * epsinv := by rw [mul_assoc, eps_mul_epsinv, mul_one]
        _ = v ^ 2 * epsinv := by rw [hv']
        _ = eps * (v * epsinv) ^ 2 := by linear_combination (-(v ^ 2 * epsinv)) * eps_mul_epsinv
    · have hv' : u * eps = -v ^ 2 := hv
      refine ⟨v * epsinv, Or.inr (Or.inr (Or.inr ?_))⟩
      calc u = u * eps * epsinv := by rw [mul_assoc, eps_mul_epsinv, mul_one]
        _ = -v ^ 2 * epsinv := by rw [hv']
        _ = -(eps * (v * epsinv) ^ 2) := by linear_combination (v ^ 2 * epsinv) * eps_mul_epsinv

/-! ### The embedding `ℤ√2 → ℚ(√2)` -/

theorem s2_mul_s2 : s2 * s2 = 2 := by
  apply (algebraMap Q2 ℝ).injective
  rw [map_mul, algebraMap_s2, map_ofNat]
  exact Real.mul_self_sqrt (by norm_num)

/-- `ℤ√2 → ℚ(√2)`, `√2 ↦ s2`. -/
noncomputable def iota2 : ℤ√2 →+* Q2 := Zsqrtd.lift ⟨s2, by rw [s2_mul_s2]; norm_num⟩

theorem iota2_apply (z : ℤ√2) : iota2 z = (z.re : Q2) + (z.im : Q2) * s2 := by
  simp [iota2]

theorem iota2_injective : Function.Injective iota2 := Zsqrtd.lift_injective _ nonsq

theorem iota2_ne_zero {z : ℤ√2} (hz : z ≠ 0) : iota2 z ≠ 0 := by
  intro h; apply hz; apply iota2_injective; rw [h, map_zero]

/-- `{a + b√2}` as a subalgebra of `ℝ`. -/
noncomputable def S2 : Subalgebra ℚ ℝ where
  carrier := {z | ∃ a b : ℚ, z = a + b * √2}
  mul_mem' := by
    rintro _ _ ⟨a, b, rfl⟩ ⟨c, e, rfl⟩
    refine ⟨a * c + 2 * b * e, a * e + b * c, ?_⟩
    have h : √2 * √2 = (2 : ℝ) := Real.mul_self_sqrt (by norm_num)
    push_cast
    linear_combination ((b : ℝ) * e) * h
  one_mem' := ⟨1, 0, by simp⟩
  add_mem' := by
    rintro _ _ ⟨a, b, rfl⟩ ⟨c, e, rfl⟩
    exact ⟨a + c, b + e, by push_cast; ring⟩
  zero_mem' := ⟨0, 0, by simp⟩
  algebraMap_mem' := fun r => ⟨r, 0, by simp⟩

theorem mem_S2 (z : Q2) : ∃ a b : ℚ, (z : ℝ) = a + b * √2 := by
  have h1 : (z : ℝ) ∈ Q2.toSubalgebra := z.2
  rw [IntermediateField.adjoin_simple_toSubalgebra_of_isAlgebraic isIntegral_sqrt_two.isAlgebraic]
    at h1
  have h2 : Algebra.adjoin ℚ {(√2 : ℝ)} ≤ S2 :=
    Algebra.adjoin_le (Set.singleton_subset_iff.mpr ⟨0, 1, by simp⟩)
  exact h2 h1

/-- Every element of `ℚ(√2)` is a quotient of two elements of `ℤ√2`. -/
theorem Q2_frac (z : Q2) : ∃ p q : ℤ√2, q ≠ 0 ∧ z = iota2 p / iota2 q := by
  obtain ⟨a, b, h⟩ := mem_S2 z
  refine ⟨⟨a.num * b.den, b.num * a.den⟩, ⟨a.den * b.den, 0⟩, ?_, ?_⟩
  · intro h0
    have := congrArg Zsqrtd.re h0
    simp only [Zsqrtd.re_zero] at this
    exact absurd this (by positivity)
  · apply (algebraMap Q2 ℝ).injective
    rw [map_div₀, iota2_apply, iota2_apply]
    simp only [map_add, map_mul, map_intCast, algebraMap_s2]
    have hz : algebraMap Q2 ℝ z = a + b * √2 := h
    rw [hz, Rat.cast_def a, Rat.cast_def b]
    have hda : (a.den : ℝ) ≠ 0 := by positivity
    have hdb : (b.den : ℝ) ≠ 0 := by positivity
    push_cast
    field_simp
    ring

end MagicSquareOfSquares.Zsqrt2

#print axioms MagicSquareOfSquares.Zsqrt2.norm_mod_lt
#print axioms MagicSquareOfSquares.Zsqrt2.pi2_prime
#print axioms MagicSquareOfSquares.Zsqrt2.three_prime
#print axioms MagicSquareOfSquares.Zsqrt2.unit_classes
#print axioms MagicSquareOfSquares.Zsqrt2.Q2_frac
