import LeanRank.RankZero
import EllipticCurves.ReductionAtPrime
import EllipticCurves.Mathlib.EllipticCurvePoint

/-!
# `E'(ℚ)` has exactly 8 points; Euler's `X⁴ + 14X²Y² + Y⁴ = Z²`; `E_B(ℚ)` has 4 points

Written by an AI model; not read by a person.  Experiment, not counted.

From `rank E'(ℚ) = 0` (`LeanRank.RankZero.finrank_eq_zero`) every point of
`E' : y² = x(x − 12)(x − 16)` has finite order.  M. Stoll's reduction map at `5`
(`EllipticCurves.ReductionAtPrime`: a group homomorphism, injective on torsion at a prime of good
reduction with `p` odd and unramified) is then injective on `E'(ℚ)`, and the reduction has `8`
points (`decide` over `ZMod 5`).  The 8 known points `O, (0,0), (12,0), (16,0), (8,±16), (24,±48)`
are therefore all of `E'(ℚ)` (`points_E`).

Consequences, all for rational numbers:
* `euler14`: `X⁴ + 14X²Y² + Y⁴ = Z²` implies `XY(X² − Y²) = 0` (Euler; the case `d = 14` of
  `x⁴ + dx²y² + y⁴ = z²`), through the map `(u, v) ↦ (14 + 2u² − 2v, −2u(14 + 2u² − 2v))` from
  `v² = u⁴ + 14u² + 1` to `E'`;
* `lemmaB`: `m⁸ + 14m⁴n⁴ + n⁸ = k²` implies `mn(m² − n²) = 0` (Lemma B of `UNIFORM_RANK.md`);
* `points_EB`, `finrank_EB`: `E_B : y² = x(x² + 14x + 1)` has the 4 points `O, (0, 0), (1, ±4)`,
  so its rank is `0` too (through the 2-isogeny `(x, y) ↦ (y²/x², y(1 − x²)/x²)`, used only as a
  map of sets: it sends affine points with `x ≠ 0` to affine points of `E'`).
-/

open WeierstrassCurve WeierstrassCurve.Affine IsDedekindDomain IsDedekindDomain.HeightOneSpectrum

namespace LeanRank.RankZero

/-! ### Every point has finite order -/

theorem isOfFinAddOrder_all (P : E.Point) : IsOfFinAddOrder P := by
  have hfin : Module.Finite ℤ E.Point := by
    exact inferInstance
  have ht : Module.IsTorsion ℤ E.Point := Module.finrank_eq_zero_iff_isTorsion.mp finrank_eq_zero
  exact (isAddTorsion_iff_isTorsion_int.mpr ht) P

/-! ### Reduction modulo `5` -/

/-- `E'` over `ℤ`. -/
def W₀ : WeierstrassCurve ℤ := ⟨0, -28, 0, 192, 0⟩

theorem map_W₀ : W₀.map (algebraMap ℤ ℚ) = E := by
  ext <;> simp [W₀, E]

theorem Δ_W₀ : W₀.Δ = 9437184 := by
  norm_num [W₀, WeierstrassCurve.Δ, WeierstrassCurve.b₂, WeierstrassCurve.b₄,
    WeierstrassCurve.b₆, WeierstrassCurve.b₈]

/-- The height-one prime `(p)` of `ℤ` (as in Stoll's `InfiniteOrderExample`). -/
def intPrime (p : ℕ) [Fact p.Prime] : HeightOneSpectrum ℤ :=
  .ofPrime (p := Ideal.span {(p : ℤ)})
    (Ideal.prime_span_singleton_iff.mpr (Nat.prime_iff_prime_int.mp Fact.out))

@[simp] theorem intPrime_asIdeal (p : ℕ) [Fact p.Prime] :
    (intPrime p).asIdeal = Ideal.span {(p : ℤ)} := rfl

instance (p : ℕ) [Fact p.Prime] : DecidableEq (ℤ ⧸ (intPrime p).asIdeal) :=
  (Int.quotientSpanNatEquivZMod p).toEquiv.decidableEq

instance : Fact (Nat.Prime 5) := ⟨by norm_num⟩

instance : (redCurve (intPrime 5) W₀).IsElliptic := by
  rw [isElliptic_iff, isUnit_iff_ne_zero]
  change (W₀.map (algebraMap ℤ (ℤ ⧸ (intPrime 5).asIdeal))).Δ ≠ 0
  rw [Ne, map_Δ, Ideal.Quotient.algebraMap_eq, Ideal.Quotient.eq_zero_iff_mem,
    intPrime_asIdeal, Ideal.mem_span_singleton, Δ_W₀]
  norm_num

/-- `E'` modulo `5`: `y² = x³ + 2x² + 2x`. -/
def E5 : WeierstrassCurve (ZMod 5) := ⟨0, 2, 0, 2, 0⟩

theorem baseChange_zmod_five :
    ((W₀.toAffine ⁄ (ZMod 5)) : WeierstrassCurve _).toAffine = E5.toAffine := by
  ext <;> decide +kernel

/-- The reduction modulo `5` has `8` points. -/
theorem card_point_E5 : Fintype.card E5.toAffine.Point = 8 := by decide

/-- `ℤ ⧸ (5) ≃ ZMod 5`. -/
def residueZModAlgEquiv (p : ℕ) [Fact p.Prime] : (ℤ ⧸ (intPrime p).asIdeal) ≃ₐ[ℤ] ZMod p :=
  (Int.quotientSpanNatEquivZMod p).toIntAlgEquiv

theorem natCard_red5 : Nat.card (redCurve (intPrime 5) W₀).Point = 8 := by
  have e : (redCurve (intPrime 5) W₀).Point ≃+ E5.toAffine.Point :=
    (Point.mapEquiv (W' := W₀.toAffine) (residueZModAlgEquiv 5)).trans
      (Point.congr baseChange_zmod_five)
  rw [Nat.card_congr e.toEquiv, Nat.card_eq_fintype_card, card_point_E5]

/-- The reduction homomorphism `E'(ℚ) → E'(𝔽₅)`. -/
noncomputable def red5 : E.Point →+ (redCurve (intPrime 5) W₀).Point :=
  redHom (intPrime 5) map_W₀

/-- **Reduction modulo `5` is injective on `E'(ℚ)`** (every point is torsion). -/
theorem red5_injective : Function.Injective red5 := by
  intro P Q hPQ
  have h0 : red (intPrime 5) map_W₀ (P - Q) = 0 := by
    rw [← coe_redHom, map_sub]
    exact sub_eq_zero.mpr hPQ
  exact sub_eq_zero.mp (eq_zero_of_isOfFinAddOrder_of_red_eq_zero (intPrime 5) map_W₀ (p := 5)
    (by norm_num) (Ideal.mem_span_singleton_self _)
    (by rw [intPrime_asIdeal, Ideal.span_singleton_pow, Ideal.mem_span_singleton]; norm_num)
    (isOfFinAddOrder_all _) h0)

instance : Finite E.Point := by
  have : Finite (redCurve (intPrime 5) W₀).Point :=
    Nat.finite_of_card_ne_zero (by rw [natCard_red5]; norm_num)
  exact Finite.of_injective _ red5_injective

theorem card_le_eight : Nat.card E.Point ≤ 8 := by
  have : Finite (redCurve (intPrime 5) W₀).Point :=
    Nat.finite_of_card_ne_zero (by rw [natCard_red5]; norm_num)
  calc Nat.card E.Point ≤ Nat.card (redCurve (intPrime 5) W₀).Point :=
        Nat.card_le_card_of_injective _ red5_injective
    _ = 8 := natCard_red5

/-! ### The eight points -/

theorem ns {x y : ℚ} (h : y ^ 2 = x * (x - 12) * (x - 16)) : E.Nonsingular x y :=
  E.equation_iff_nonsingular.mp ((equation_iff' x y).mpr h)

/-- The eight known points of `E'(ℚ)`. -/
noncomputable def pts : Fin 8 → E.Point :=
  ![0, .some 0 0 (ns (by norm_num)), .some 12 0 (ns (by norm_num)),
    .some 16 0 (ns (by norm_num)), .some 8 16 (ns (by norm_num)),
    .some 8 (-16) (ns (by norm_num)), .some 24 48 (ns (by norm_num)),
    .some 24 (-48) (ns (by norm_num))]

/-- Coordinates of a point (`none` for `O`). -/
def key : E.Point → Option (ℚ × ℚ)
  | .zero => none
  | .some x y _ => some (x, y)

theorem key_injective : Function.Injective key := by
  intro P Q h
  rcases P with _ | ⟨x, y, _⟩ <;> rcases Q with _ | ⟨x', y', _⟩
  · rfl
  · simp [key] at h
  · simp [key] at h
  · simp only [key, Option.some.injEq, Prod.mk.injEq] at h
    obtain ⟨rfl, rfl⟩ := h
    rfl

theorem pts_injective : Function.Injective pts := by
  have : Function.Injective (key ∘ pts) := by decide +kernel
  exact this.of_comp

/-- **`E'(ℚ)` is exactly the eight points.** -/
theorem points_E (P : E.Point) : ∃ i, pts i = P :=
  (pts_injective.bijective_of_nat_card_le (by rw [Nat.card_fin]; exact card_le_eight)).2 P

theorem card_E : Nat.card E.Point = 8 :=
  le_antisymm card_le_eight (by
    have := Nat.card_le_card_of_injective _ pts_injective
    rwa [Nat.card_fin] at this)

/-- **The affine rational points of `y² = x(x − 12)(x − 16)`.** -/
theorem affine_points_E {x y : ℚ} (h : y ^ 2 = x * (x - 12) * (x - 16)) :
    (x, y) ∈ ([(0, 0), (12, 0), (16, 0), (8, 16), (8, -16), (24, 48), (24, -48)] :
      List (ℚ × ℚ)) := by
  obtain ⟨i, hi⟩ := points_E (.some x y (ns h))
  have hk := congrArg key hi
  fin_cases i <;> simp [pts, key] at hk <;> obtain ⟨rfl, rfl⟩ := hk <;> simp

/-! ### Euler's equation `X⁴ + 14X²Y² + Y⁴ = Z²` and Lemma B -/

/-- **`v² = u⁴ + 14u² + 1` has only the rational solutions with `u ∈ {0, 1, −1}`.** -/
theorem quartic {u v : ℚ} (h : v ^ 2 = u ^ 4 + 14 * u ^ 2 + 1) : u = 0 ∨ u = 1 ∨ u = -1 := by
  have hE : (-2 * u * (14 + 2 * u ^ 2 - 2 * v)) ^ 2 =
      (14 + 2 * u ^ 2 - 2 * v) * ((14 + 2 * u ^ 2 - 2 * v) - 12) *
        ((14 + 2 * u ^ 2 - 2 * v) - 16) := by
    linear_combination (-4 * (14 + 2 * u ^ 2 - 2 * v)) * h
  have hm := affine_points_E hE
  simp only [List.mem_cons, Prod.mk.injEq, List.not_mem_nil, or_false] at hm
  rcases hm with ⟨h1, -⟩ | ⟨h1, h2⟩ | ⟨h1, h2⟩ | ⟨h1, h2⟩ | ⟨h1, h2⟩ | ⟨h1, h2⟩ | ⟨h1, h2⟩
  · exfalso
    have hv : v = 7 + u ^ 2 := by linarith
    rw [hv] at h
    have h48 : (48 : ℚ) = 0 := by linear_combination h
    norm_num at h48
  all_goals rw [h1] at h2
  · left; linarith
  · left; linarith
  · right; right; linarith
  · right; left; linarith
  · right; right; linarith
  · right; left; linarith

/-- **Euler (`d = 14`)**: if `X⁴ + 14X²Y² + Y⁴ = Z²` in rationals, then `XY(X² − Y²) = 0`. -/
theorem euler14 {X Y Z : ℚ} (h : X ^ 4 + 14 * X ^ 2 * Y ^ 2 + Y ^ 4 = Z ^ 2) :
    X * Y * (X ^ 2 - Y ^ 2) = 0 := by
  by_cases hY : Y = 0
  · simp [hY]
  have key : (Z / Y ^ 2) ^ 2 - ((X / Y) ^ 4 + 14 * (X / Y) ^ 2 + 1) =
      (Z ^ 2 - (X ^ 4 + 14 * X ^ 2 * Y ^ 2 + Y ^ 4)) / Y ^ 4 := by
    field_simp
  rw [h, sub_self, zero_div, sub_eq_zero] at key
  rcases quartic key with h0 | h1 | h2
  · rw [div_eq_zero_iff] at h0
    rcases h0 with h0 | h0
    · simp [h0]
    · exact absurd h0 hY
  · rw [div_eq_one_iff_eq hY] at h1
    simp [h1]
  · rw [div_eq_iff hY] at h2
    rw [h2]; ring

/-- **Lemma B** of `UNIFORM_RANK.md`: if `m⁸ + 14m⁴n⁴ + n⁸` is the square of a rational, then
`mn(m² − n²) = 0`. -/
theorem lemmaB {m n k : ℚ} (h : m ^ 8 + 14 * m ^ 4 * n ^ 4 + n ^ 8 = k ^ 2) :
    m * n * (m ^ 2 - n ^ 2) = 0 := by
  have h' : (m ^ 2) ^ 4 + 14 * (m ^ 2) ^ 2 * (n ^ 2) ^ 2 + (n ^ 2) ^ 4 = k ^ 2 := by
    rw [← h]; ring
  have he := euler14 h'
  rcases mul_eq_zero.mp he with h1 | h1
  · rcases mul_eq_zero.mp h1 with h2 | h2
    · simp [pow_eq_zero_iff (two_ne_zero) |>.mp h2]
    · simp [pow_eq_zero_iff (two_ne_zero) |>.mp h2]
  · have h3 : (m ^ 2 - n ^ 2) * (m ^ 2 + n ^ 2) = 0 := by linear_combination h1
    rcases mul_eq_zero.mp h3 with h4 | h4
    · rw [h4, mul_zero]
    · have hm2 : m ^ 2 = 0 := by nlinarith [sq_nonneg m, sq_nonneg n]
      simp [pow_eq_zero_iff (two_ne_zero) |>.mp hm2]

/-! ### `E_B : y² = x(x² + 14x + 1)` -/

theorem not_sq_nat {N : ℕ} (hN : 2 ≤ N) (hb : ∀ k : ℕ, k < N → k * k ≠ N) {r : ℚ}
    (h : r ^ 2 = (N : ℚ)) : False := by
  obtain ⟨z, hz⟩ := int_sq_of_rat_sq (N := (N : ℤ)) (by exact_mod_cast h)
  have hk : z.natAbs * z.natAbs = N := by
    have : ((z.natAbs * z.natAbs : ℕ) : ℤ) = (N : ℤ) := by rw [Int.natAbs_mul_self, ← sq, hz]
    exact_mod_cast this
  refine hb z.natAbs ?_ hk
  by_contra hc
  have hc' := not_lt.mp hc
  nlinarith

theorem not_sq_48 {r : ℚ} (h : r ^ 2 = 48) : False :=
  not_sq_nat (N := 48) (by norm_num) (by decide) (by exact_mod_cast h)
theorem not_sq_12 {r : ℚ} (h : r ^ 2 = 12) : False :=
  not_sq_nat (N := 12) (by norm_num) (by decide) (by exact_mod_cast h)
theorem not_sq_8 {r : ℚ} (h : r ^ 2 = 8) : False :=
  not_sq_nat (N := 8) (by norm_num) (by decide) (by exact_mod_cast h)
theorem not_sq_24 {r : ℚ} (h : r ^ 2 = 24) : False :=
  not_sq_nat (N := 24) (by norm_num) (by decide) (by exact_mod_cast h)

/-- **The affine rational points of `y² = x(x² + 14x + 1)`**: `(0, 0)` and `(1, ±4)`. -/
theorem affine_points_EB {x y : ℚ} (h : y ^ 2 = x * (x ^ 2 + 14 * x + 1)) :
    (x = 0 ∧ y = 0) ∨ (x = 1 ∧ (y = 4 ∨ y = -4)) := by
  by_cases hx : x = 0
  · left
    refine ⟨hx, ?_⟩
    rw [hx] at h
    exact pow_eq_zero_iff (two_ne_zero) |>.mp (by linear_combination h)
  right
  -- the image under the 2-isogeny
  have hE : (y * (1 - x ^ 2) / x ^ 2) ^ 2 = ((x ^ 2 + 14 * x + 1) / x) *
      ((x ^ 2 + 14 * x + 1) / x - 12) * ((x ^ 2 + 14 * x + 1) / x - 16) := by
    have key : (y * (1 - x ^ 2) / x ^ 2) ^ 2 - ((x ^ 2 + 14 * x + 1) / x) *
        ((x ^ 2 + 14 * x + 1) / x - 12) * ((x ^ 2 + 14 * x + 1) / x - 16) =
        (1 - x ^ 2) ^ 2 / x ^ 4 * (y ^ 2 - x * (x ^ 2 + 14 * x + 1)) := by
      field_simp
      ring
    rw [h, sub_self, mul_zero, sub_eq_zero] at key
    exact key
  have hm := affine_points_E hE
  simp only [List.mem_cons, Prod.mk.injEq, List.not_mem_nil, or_false] at hm
  have hX : ∀ c : ℚ, (x ^ 2 + 14 * x + 1) / x = c → x ^ 2 + (14 - c) * x + 1 = 0 := by
    intro c hc
    rw [div_eq_iff hx] at hc
    linear_combination hc
  rcases hm with ⟨h1, -⟩ | ⟨h1, -⟩ | ⟨h1, -⟩ | ⟨h1, -⟩ | ⟨h1, -⟩ | ⟨h1, -⟩ | ⟨h1, -⟩
  · exact (not_sq_48 (r := x + 7) (by linear_combination hX 0 h1)).elim
  · -- x = -1, then y² = 12
    have hx1 : x = -1 := by
      have : (x + 1) ^ 2 = 0 := by linear_combination hX 12 h1
      linarith [pow_eq_zero_iff (two_ne_zero) |>.mp this]
    exact (not_sq_12 (r := y) (by rw [h, hx1]; norm_num)).elim
  · have hx1 : x = 1 := by
      have : (x - 1) ^ 2 = 0 := by linear_combination hX 16 h1
      linarith [pow_eq_zero_iff (two_ne_zero) |>.mp this]
    refine ⟨hx1, ?_⟩
    have : (y - 4) * (y + 4) = 0 := by rw [hx1] at h; linear_combination h
    rcases mul_eq_zero.mp this with h4 | h4
    · left; linarith
    · right; linarith
  · exact (not_sq_8 (r := x + 3) (by linear_combination hX 8 h1)).elim
  · exact (not_sq_8 (r := x + 3) (by linear_combination hX 8 h1)).elim
  · exact (not_sq_24 (r := x - 5) (by linear_combination hX 24 h1)).elim
  · exact (not_sq_24 (r := x - 5) (by linear_combination hX 24 h1)).elim

/-- `E_B : y² = x³ + 14x² + x`. -/
def EB : Affine ℚ := ⟨0, 14, 0, 1, 0⟩

theorem EB_equation {x y : ℚ} (h : EB.Equation x y) : y ^ 2 = x * (x ^ 2 + 14 * x + 1) := by
  rw [Affine.equation_iff] at h
  simp only [EB] at h
  linear_combination h

/-- Coordinates of a point of `E_B` (`none` for `O`). -/
def keyB : EB.Point → Option (ℚ × ℚ)
  | .zero => none
  | .some x y _ => some (x, y)

theorem keyB_injective : Function.Injective keyB := by
  intro P Q h
  rcases P with _ | ⟨x, y, _⟩ <;> rcases Q with _ | ⟨x', y', _⟩
  · rfl
  · simp [keyB] at h
  · simp [keyB] at h
  · simp only [keyB, Option.some.injEq, Prod.mk.injEq] at h
    obtain ⟨rfl, rfl⟩ := h
    rfl

/-- **`E_B(ℚ) = {O, (0, 0), (1, 4), (1, −4)}`**, so `E_B(ℚ)` is finite. -/
theorem points_EB (P : EB.Point) :
    keyB P ∈ ({none, some (0, 0), some (1, 4), some (1, -4)} : Set (Option (ℚ × ℚ))) := by
  rcases P with _ | ⟨x, y, hP⟩
  · simp [keyB]
  · rcases affine_points_EB (EB_equation hP.1) with ⟨rfl, rfl⟩ | ⟨rfl, rfl | rfl⟩ <;>
      simp [keyB]

instance : Finite EB.Point :=
  Finite.of_injective (fun P => (⟨keyB P, points_EB P⟩ :
    ({none, some (0, 0), some (1, 4), some (1, -4)} : Set (Option (ℚ × ℚ)))))
    (fun P Q h => keyB_injective (congrArg Subtype.val h))

/-- **`rank E_B(ℚ) = 0`** for `E_B : y² = x(x² + 14x + 1)`, the curve of Lemma B. -/
theorem finrank_EB : Module.finrank ℤ EB.Point = 0 := by
  have : Module.Finite ℤ EB.Point := Module.Finite.of_finite
  exact Module.finrank_eq_zero_iff_isTorsion.mpr
    (isAddTorsion_iff_isTorsion_int.mp isAddTorsion_of_finite)

end LeanRank.RankZero

#print axioms LeanRank.RankZero.red5_injective
#print axioms LeanRank.RankZero.points_E
#print axioms LeanRank.RankZero.card_E
#print axioms LeanRank.RankZero.affine_points_E
#print axioms LeanRank.RankZero.quartic
#print axioms LeanRank.RankZero.euler14
#print axioms LeanRank.RankZero.lemmaB
#print axioms LeanRank.RankZero.affine_points_EB
#print axioms LeanRank.RankZero.finrank_EB
