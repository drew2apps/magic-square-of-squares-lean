import Mathlib.AlgebraicGeometry.EllipticCurve.Affine.Point
import Mathlib.Data.Rat.Lemmas
import Mathlib.Tactic
import MagicSquareOfSquares.G1SigmaGeneric

/-!
# The 2-descent map for `y² = (X − e)(X² + aX + b)`, and `UNIFORM_RANK.md` Theorem 2 at given ratios

Written by an AI model (Claude, made by Anthropic) on 25 November 2026 (UTC); not read by a person.
Experiment, not counted in the census.

## Setting

`K` is any field with decidable equality; `curve e a b` is the Weierstrass curve with `a₁ = a₃ = 0`, `a₂ = a − e`,
`a₄ = b − ae`, `a₆ = −eb`, that is `y² = f(X) = (X − e)(X² + aX + b)`, and `q(e) = e² + ae + b` (`= f′(e)`). Points
are Mathlib's `WeierstrassCurve.Affine.Point` with Mathlib's group law. No hypothesis on the characteristic is used
(in characteristic 2 every point of this model satisfies `P + P = O`, and the statements below are still true).

## What is proved (all on Mathlib's points unless marked "coordinates")

1. Chord (`chord_coord`, coordinates; `chord_point`, on points): if `P₁ + P₂ = P₃` with `x₁ ≠ x₂` and `ℓ` is the
   slope, then `(x₁ − e)(x₂ − e)(x₃ − e) = (ℓ(e − x₁) + y₁)²`. Also `chord_T_coord`, `chord_third_T_coord` (lines
   through `T = (e, 0)`).
2. Tangent (`tangent_coord`, coordinates): `4y²(x(2P) − e) = ((x − e)² − q(e))²`.
3. Doubles: if `W + W = (x, y)` then `x − e` is a square in `K` (`double_isSquare`, no hypothesis); if moreover
   `q(e) ≠ 0`, the same through `α` (`double_x_sub_isSquare`), and `q(e)` is a square when `x = e`
   (`double_T_isSquare`). The halving test at the other
   two roots (`gauss_double`): if `a = −2r`, `b = r² + s²` (roots `r ± is` in `K(i)`) and `W + W = (x, y)`, then
   `x − r = u² − v²` and `−s = 2uv` for some `u, v ∈ K`, i.e. `x − (r + is)` is a square in `K(i)`; stated over `K`,
   so no field `K(i)` is constructed.
4. The 2-descent map: `alpha` (`α(O) = 1`, `α(e, 0) = q(e)`, `α(x, y) = x − e`), `alpha_add`
   (`α(P₁ + P₂)·α(P₁)·α(P₂)` is a square for all `P₁, P₂`, if `q(e) ≠ 0`), and the group homomorphism
   `alphaHom : E(K) →+ Additive (Kˣ/(Kˣ)²)` with `2E(K)` in its kernel (`alphaHom_two_smul`).
5. `UNIFORM_RANK.md` Theorem 2 at given ratios. `Eplus m n` is `E⁺` at `(m, n)` (its cubic is `genCub` of
   `G1SigmaGeneric.lean`, `Eplus_cubic`); `E⁻` at `(m, n)` is `Eplus m (−n)`. `theorem2_of`: if
   `m² + n² = 2k²`, `m, n, m ± n ≠ 0`, and none of `H(m, n)`, `H(m, −n)`, `H(m, n)H(m, −n)`, `2mn`, `−2mn` is a
   square in `ℚ`, then `R`, `Q`, `T` are independent in `E(ℚ)/2E(ℚ)` (`Theorem2At`: `c₁R + c₂Q + c₃T = W + W`
   with `cᵢ ∈ ℤ` forces all `cᵢ` even). Certificates: `theorem2_7_1`, `theorem2_7_m1` here, and all 127 ratios
   with `m ≤ 1000`, both curves, in `TwoDescentCertificates.lean`. `two_torsion_Q`, `two_torsion_Eplus`: over `ℚ`,
   if `s ≠ 0`, the only point of order 2 is `T`. `Eplus_isElliptic`: `Δ = −64·q(e)²·s² ≠ 0` when `m, n, m ± n ≠ 0`.
   `R_add_R`: `R + R = P = (−(m + n)²S², S³(m − n)(m + n)²)` on Mathlib's points, for every such ratio.
6. Negative controls: `control_7_1` (Mathlib's `R + R` at `(7, 1)` is the point `P` of `z = 2`, `α(P) = 300²`, and
   the halving test at `θ` passes on it) and `control_bracket` (`300²` has no square-root bracket).

## Not in this file (out of scope)

* Lemma D (the genus-2 curve `D2`) and Lemma B (Euler, `x⁴ + 14x²y² + y⁴ = z²`) of `UNIFORM_RANK.md`, which give
  the five non-square hypotheses of `theorem2_of` for every twice-square ratio at once and so make Theorem 2
  uniform. Here they are replaced, ratio by ratio, by square-root brackets.
* "rank ≥ 2": that needs Mordell–Weil and `dim E(ℚ)/2E(ℚ) = rank + dim E(ℚ)[2]`, which are not in Mathlib.
* The rest of Theorem 1 of `UNIFORM_RANK.md`: that `R` has infinite order (it needs Lemma D); `P = 2R` is
  `R_add_R`.
* That the Weierstrass model of `E⁻` at `(m, n)` equals `Eplus m (−n)` is `UNIFORM_RANK.md` A2 (sympy); on the
  quartics it is `gen_minus_plus`.
-/

namespace MagicSquareOfSquares.TwoDescent

open WeierstrassCurve WeierstrassCurve.Affine

section General

variable {K : Type*} [Field K]

/-- The Weierstrass curve `y² = (X − e)(X² + aX + b)`: `a₁ = a₃ = 0`, `a₂ = a − e`, `a₄ = b − ae`,
`a₆ = −eb`. -/
def curve (e a b : K) : WeierstrassCurve K where
  a₁ := 0
  a₂ := a - e
  a₃ := 0
  a₄ := b - a * e
  a₆ := -(e * b)

/-- `q(e) = e² + ae + b`, the value of the quadratic factor at `e` (it is `f′(e)`). -/
def qe (e a b : K) : K := e ^ 2 + a * e + b

variable {e a b : K}

@[simp] theorem curve_a₁ : (curve e a b).toAffine.a₁ = 0 := rfl
@[simp] theorem curve_a₂ : (curve e a b).toAffine.a₂ = a - e := rfl
@[simp] theorem curve_a₃ : (curve e a b).toAffine.a₃ = 0 := rfl
@[simp] theorem curve_a₄ : (curve e a b).toAffine.a₄ = b - a * e := rfl
@[simp] theorem curve_a₆ : (curve e a b).toAffine.a₆ = -(e * b) := rfl

theorem curve_equation_iff (x y : K) :
    (curve e a b).toAffine.Equation x y ↔ y ^ 2 = (x - e) * (x ^ 2 + a * x + b) := by
  rw [Affine.equation_iff]
  simp only [curve_a₁, curve_a₂, curve_a₃, curve_a₄, curve_a₆]
  constructor <;> intro h <;> linear_combination h

theorem curve_negY (x y : K) : (curve e a b).toAffine.negY x y = -y := by
  simp only [Affine.negY, curve_a₁, curve_a₃]
  ring

/-- The discriminant: `Δ = 16·q(e)²·(a² − 4b)` (`16` times the discriminant of the cubic). -/
theorem curve_Δ : (curve e a b).Δ = 16 * qe e a b ^ 2 * (a ^ 2 - 4 * b) := by
  simp only [WeierstrassCurve.Δ, WeierstrassCurve.b₂, WeierstrassCurve.b₄, WeierstrassCurve.b₆,
    WeierstrassCurve.b₈, curve]
  unfold qe
  ring

/-! ### Step 1 and 2 in coordinates -/

/-- **Chord identity.** If `(x₁, y₁)`, `(x₂, y₂)` are on the curve, `x₁ ≠ x₂`, and `ℓ` is the slope of the
line through them, then with `x₃ = ℓ² − a₂ − x₁ − x₂` (the abscissa of the third point on the line) and
`ν = y₁ − ℓx₁`: `(x₁ − e)(x₂ − e)(x₃ − e) = (ℓe + ν)²`. -/
theorem chord_coord {x₁ x₂ y₁ y₂ ℓ : K} (h₁ : y₁ ^ 2 = (x₁ - e) * (x₁ ^ 2 + a * x₁ + b))
    (h₂ : y₂ ^ 2 = (x₂ - e) * (x₂ ^ 2 + a * x₂ + b)) (hx : x₁ ≠ x₂)
    (hℓ : ℓ * (x₁ - x₂) = y₁ - y₂) :
    (x₁ - e) * (x₂ - e) * (ℓ ^ 2 - (a - e) - x₁ - x₂ - e) = (ℓ * (e - x₁) + y₁) ^ 2 := by
  have hG : x₁ ^ 2 + x₁ * x₂ + x₂ ^ 2 + (a - e) * (x₁ + x₂) + (b - a * e) - 2 * y₁ * ℓ
      + ℓ ^ 2 * (x₁ - x₂) = 0 := by
    have h : (x₁ - x₂) * (x₁ ^ 2 + x₁ * x₂ + x₂ ^ 2 + (a - e) * (x₁ + x₂) + (b - a * e)
        - 2 * y₁ * ℓ + ℓ ^ 2 * (x₁ - x₂)) = 0 := by
      linear_combination -h₁ + h₂ + (ℓ * (x₁ - x₂) - y₁ - y₂) * hℓ
    exact (mul_eq_zero.mp h).resolve_left (sub_ne_zero.mpr hx)
  linear_combination (e - x₁) * hG - h₁

/-- **Translation by `T = (e, 0)`.** If `(x₂, y₂)` is on the curve, `x₂ ≠ e`, and `ℓ` is the slope of the line
through `T` and `(x₂, y₂)`, then the third point `x₃ = ℓ² − a₂ − e − x₂` satisfies `(x₂ − e)(x₃ − e) = q(e)`. -/
theorem chord_T_coord {x₂ y₂ ℓ : K} (h₂ : y₂ ^ 2 = (x₂ - e) * (x₂ ^ 2 + a * x₂ + b)) (hx : x₂ ≠ e)
    (hℓ : ℓ * (x₂ - e) = y₂) : (x₂ - e) * (ℓ ^ 2 - (a - e) - e - x₂ - e) = qe e a b := by
  have h : (x₂ - e) * (ℓ ^ 2 * (x₂ - e) - (x₂ ^ 2 + a * x₂ + b)) = 0 := by
    linear_combination (ℓ * (x₂ - e) + y₂) * hℓ + h₂
  have h' := (mul_eq_zero.mp h).resolve_left (sub_ne_zero.mpr hx)
  unfold qe
  linear_combination h'

/-- **The chord through `T`.** If the line of slope `ℓ` through `(x₁, y₁)` (on the curve, `x₁ ≠ e`) passes through
`T` (`ℓ(e − x₁) + y₁ = 0`) and its third point has abscissa `x₂`, i.e. `ℓ² − a₂ − x₁ − x₂ = e`, then
`(x₁ − e)(x₂ − e) = q(e)`. -/
theorem chord_third_T_coord {x₁ x₂ y₁ ℓ : K} (h₁ : y₁ ^ 2 = (x₁ - e) * (x₁ ^ 2 + a * x₁ + b))
    (hx : x₁ ≠ e) (hL : ℓ * (e - x₁) + y₁ = 0) (h3 : ℓ ^ 2 - (a - e) - x₁ - x₂ = e) :
    (x₁ - e) * (x₂ - e) = qe e a b := by
  have h : (x₁ - e) * (ℓ ^ 2 * (x₁ - e) - (x₁ ^ 2 + a * x₁ + b)) = 0 := by
    linear_combination (-(ℓ * (x₁ - e) + y₁)) * hL + h₁
  have h' := (mul_eq_zero.mp h).resolve_left (sub_ne_zero.mpr hx)
  unfold qe
  linear_combination h' - (x₁ - e) * h3

/-- **Tangent identity.** If `(x, y)` is on the curve and `ℓ · 2y = f′(x)` (the tangent slope), then with
`x(2P) = ℓ² − a₂ − 2x`: `4y²·(x(2P) − e) = ((x − e)² − q(e))²`. -/
theorem tangent_coord {x y ℓ : K} (h : y ^ 2 = (x - e) * (x ^ 2 + a * x + b))
    (hℓ : ℓ * (2 * y) = 3 * x ^ 2 + 2 * (a - e) * x + (b - a * e)) :
    4 * y ^ 2 * (ℓ ^ 2 - (a - e) - x - x - e) = ((x - e) ^ 2 - qe e a b) ^ 2 := by
  unfold qe
  linear_combination (ℓ * (2 * y) + (3 * x ^ 2 + 2 * (a - e) * x + (b - a * e))) * hℓ
    - 4 * (a + 2 * x) * h

/-- **The tangent identity at the other two roots, written over `K`.** If `a = −2r` and `b = r² + s²`, so that
the quadratic factor is `(X − r)² + s²` with roots `r ± is` in `K(i)`, and `(x, y)` is on the curve with `2y ≠ 0`
and tangent slope `ℓ`, then `x(2P) − (r + is) = (u + iv)²` for some `u, v ∈ K`; as two equations over `K`:
`x(2P) − r = u² − v²` and `−s = 2uv`. (`u + iv = ((x − r)² + s² − 2is(x − e))/(2y)`.) -/
theorem gauss_coord {r s x y ℓ : K} (ha : a = -2 * r) (hb : b = r ^ 2 + s ^ 2)
    (h : y ^ 2 = (x - e) * (x ^ 2 + a * x + b)) (hy : 2 * y ≠ 0)
    (hℓ : ℓ * (2 * y) = 3 * x ^ 2 + 2 * (a - e) * x + (b - a * e)) :
    ∃ u v : K, ℓ ^ 2 - (a - e) - x - x - r = u ^ 2 - v ^ 2 ∧ -s = 2 * u * v := by
  subst ha hb
  have k1 : 4 * y ^ 2 * (ℓ ^ 2 - (-2 * r - e) - x - x - r) =
      ((x - r) ^ 2 + s ^ 2) ^ 2 - (-2 * s * (x - e)) ^ 2 := by
    linear_combination (ℓ * (2 * y) + (3 * x ^ 2 + 2 * (-2 * r - e) * x
      + (r ^ 2 + s ^ 2 - -2 * r * e))) * hℓ - 4 * (-2 * r - e + 2 * x + r) * h
  have k2 : 4 * y ^ 2 * (-s) = 2 * ((x - r) ^ 2 + s ^ 2) * (-2 * s * (x - e)) := by
    linear_combination (-4 * s) * h
  obtain ⟨w, hw⟩ : ∃ w : K, 2 * y * w = 1 := ⟨(2 * y)⁻¹, mul_inv_cancel₀ hy⟩
  refine ⟨((x - r) ^ 2 + s ^ 2) * w, (-2 * s * (x - e)) * w, ?_, ?_⟩
  · linear_combination (-(ℓ ^ 2 - (-2 * r - e) - x - x - r) * (1 + 2 * y * w)) * hw + w ^ 2 * k1
  · linear_combination (s * (1 + 2 * y * w)) * hw + w ^ 2 * k2

/-- In a field: if `y` is a non-zero square and `x·y` is a square, then `x` is a square. -/
theorem isSquare_cancel {x y : K} (hy : IsSquare y) (hy0 : y ≠ 0) (h : IsSquare (x * y)) :
    IsSquare x := by
  obtain ⟨r, hr⟩ := hy
  obtain ⟨t, ht⟩ := h
  have hr0 : r ≠ 0 := by rintro rfl; exact hy0 (by rw [hr, mul_zero])
  obtain ⟨w, hw⟩ : ∃ w : K, r * w = 1 := ⟨r⁻¹, mul_inv_cancel₀ hr0⟩
  refine ⟨t * w, ?_⟩
  linear_combination (-x * (1 + r * w)) * hw + w ^ 2 * ht - w ^ 2 * x * hr

theorem isSquare_cancel_left {x y : K} (hx : IsSquare x) (hx0 : x ≠ 0) (h : IsSquare (x * y)) :
    IsSquare y :=
  isSquare_cancel hx hx0 (by rwa [mul_comm] at h)

/-! ### Points (Mathlib's `WeierstrassCurve.Affine.Point`) -/

theorem nonsingular_of_y {x y : K} (h : y ^ 2 = (x - e) * (x ^ 2 + a * x + b)) (hy : 2 * y ≠ 0) :
    (curve e a b).toAffine.Nonsingular x y := by
  rw [Affine.nonsingular_iff']
  refine ⟨(curve_equation_iff x y).mpr h, Or.inr ?_⟩
  simp only [curve_a₁, curve_a₃]
  intro h0; apply hy; linear_combination h0

theorem nonsingular_T (hq : qe e a b ≠ 0) : (curve e a b).toAffine.Nonsingular e 0 := by
  rw [Affine.nonsingular_iff']
  refine ⟨(curve_equation_iff e 0).mpr (by ring), Or.inl ?_⟩
  simp only [curve_a₁, curve_a₂, curve_a₄]
  intro h0; apply hq; unfold qe; linear_combination -h0

variable [DecidableEq K]

theorem slope_double {x y : K} (hy : y ≠ (curve e a b).toAffine.negY x y) :
    (curve e a b).toAffine.slope x x y y * (2 * y) = 3 * x ^ 2 + 2 * (a - e) * x + (b - a * e) ∧
      2 * y ≠ 0 := by
  have h2 : 2 * y ≠ 0 := by
    intro h; apply hy; rw [curve_negY]; linear_combination h
  refine ⟨?_, h2⟩
  rw [Affine.slope_of_Y_ne rfl hy, curve_negY]
  simp only [curve_a₁, curve_a₂, curve_a₄]
  rw [show y - -y = 2 * y by ring, div_mul_cancel₀ _ h2]
  ring

/-- **The 2-descent map** (as a function to `K`): `α(O) = 1`, `α(T) = q(e)` for `T = (e, 0)`, and
`α(x, y) = x − e` otherwise. Its class in `Kˣ/(Kˣ)²` is the usual map (Silverman–Tate III.5). -/
def alpha (e a b : K) : (curve e a b).toAffine.Point → K
  | .zero => 1
  | .some x _ _ => if x = e then qe e a b else x - e

@[simp] theorem alpha_zero : alpha e a b 0 = 1 := rfl

theorem alpha_some {x y : K} (h : (curve e a b).toAffine.Nonsingular x y) :
    alpha e a b (.some x y h) = if x = e then qe e a b else x - e := rfl

theorem alpha_ne_zero (hq : qe e a b ≠ 0) (P : (curve e a b).toAffine.Point) :
    alpha e a b P ≠ 0 := by
  rcases P with _ | ⟨x, y, h⟩
  · exact one_ne_zero
  · rw [alpha_some]
    split_ifs with hx
    · exact hq
    · exact sub_ne_zero.mpr hx

/-- **`α` is multiplicative up to squares**: for all points `P₁, P₂`, `α(P₁ + P₂)·α(P₁)·α(P₂)` is a square in
`K`. Hypothesis: `q(e) ≠ 0` (the point `T = (e, 0)` is nonsingular). No hypothesis on the characteristic. -/
theorem alpha_add (hq : qe e a b ≠ 0) (P₁ P₂ : (curve e a b).toAffine.Point) :
    IsSquare (alpha e a b (P₁ + P₂) * alpha e a b P₁ * alpha e a b P₂) := by
  rcases P₁ with _ | ⟨x₁, y₁, h₁⟩
  · show IsSquare (alpha e a b (0 + P₂) * alpha e a b 0 * alpha e a b P₂)
    rw [zero_add, alpha_zero]
    exact ⟨alpha e a b P₂, by ring⟩
  rcases P₂ with _ | ⟨x₂, y₂, h₂⟩
  · show IsSquare (alpha e a b (Point.some x₁ y₁ h₁ + 0) * alpha e a b (Point.some x₁ y₁ h₁) *
      alpha e a b 0)
    rw [add_zero, alpha_zero]
    exact ⟨alpha e a b (Point.some x₁ y₁ h₁), by ring⟩
  have e₁ := (curve_equation_iff x₁ y₁).mp h₁.1
  have e₂ := (curve_equation_iff x₂ y₂).mp h₂.1
  by_cases hxy : x₁ = x₂ ∧ y₁ = (curve e a b).toAffine.negY x₂ y₂
  · rw [Point.add_of_Y_eq hxy.1 hxy.2, alpha_zero, alpha_some, alpha_some, hxy.1]
    exact ⟨if x₂ = e then qe e a b else x₂ - e, by ring⟩
  rw [Point.add_some hxy]
  simp only [alpha_some]
  generalize hℓ : (curve e a b).toAffine.slope x₁ x₂ y₁ y₂ = ℓ
  have hX : (curve e a b).toAffine.addX x₁ x₂ ℓ = ℓ ^ 2 - (a - e) - x₁ - x₂ := by
    simp only [Affine.addX, curve_a₁, curve_a₂]; ring
  rw [hX]
  by_cases hx : x₁ = x₂
  · -- doubling
    have hy : y₁ ≠ (curve e a b).toAffine.negY x₂ y₂ := fun h => hxy ⟨hx, h⟩
    have hyy : y₁ = y₂ := Affine.Y_eq_of_Y_ne h₁.1 h₂.1 hx hy
    subst hx hyy
    obtain ⟨hs, h2⟩ := slope_double hy
    rw [hℓ] at hs
    have key := tangent_coord e₁ hs
    have hA : IsSquare (if ℓ ^ 2 - (a - e) - x₁ - x₁ = e then qe e a b
        else ℓ ^ 2 - (a - e) - x₁ - x₁ - e) := by
      obtain ⟨w, hw⟩ : ∃ w : K, 2 * y₁ * w = 1 := ⟨(2 * y₁)⁻¹, mul_inv_cancel₀ h2⟩
      split_ifs with h3
      · refine ⟨x₁ - e, ?_⟩
        have h0 : ((x₁ - e) ^ 2 - qe e a b) ^ 2 = 0 := by
          linear_combination -key + 4 * y₁ ^ 2 * h3
        have h00 := (pow_eq_zero_iff two_ne_zero).mp h0
        linear_combination -h00
      · refine ⟨((x₁ - e) ^ 2 - qe e a b) * w, ?_⟩
        linear_combination (-(ℓ ^ 2 - (a - e) - x₁ - x₁ - e) * (1 + 2 * y₁ * w)) * hw + w ^ 2 * key
    rw [mul_assoc]
    exact hA.mul ⟨_, rfl⟩
  · -- chord
    have hs : ℓ * (x₁ - x₂) = y₁ - y₂ := by
      rw [← hℓ, Affine.slope_of_X_ne hx, div_mul_cancel₀ _ (sub_ne_zero.mpr hx)]
    have hC := chord_coord e₁ e₂ hx hs
    by_cases hx1 : x₁ = e
    · have hy1 : y₁ = 0 := by
        have h0 : y₁ ^ 2 = 0 := by rw [e₁, hx1, sub_self, zero_mul]
        exact (pow_eq_zero_iff two_ne_zero).mp h0
      have hx2 : x₂ ≠ e := fun h => hx (hx1.trans h.symm)
      have hs' : ℓ * (x₂ - e) = y₂ := by linear_combination -hs + ℓ * hx1 - hy1
      have hT := chord_T_coord e₂ hx2 hs'
      have h3 : ℓ ^ 2 - (a - e) - x₁ - x₂ ≠ e := by
        intro h3; apply hq
        linear_combination -hT + (x₂ - e) * h3 + (x₂ - e) * hx1
      rw [if_pos hx1, if_neg hx2, if_neg h3]
      exact ⟨qe e a b, by linear_combination (qe e a b) * hT - (qe e a b) * (x₂ - e) * hx1⟩
    by_cases hx2 : x₂ = e
    · have hy2 : y₂ = 0 := by
        have h0 : y₂ ^ 2 = 0 := by rw [e₂, hx2, sub_self, zero_mul]
        exact (pow_eq_zero_iff two_ne_zero).mp h0
      have hs' : ℓ * (x₁ - e) = y₁ := by linear_combination hs + ℓ * hx2 - hy2
      have hT := chord_T_coord e₁ hx1 hs'
      have h3 : ℓ ^ 2 - (a - e) - x₁ - x₂ ≠ e := by
        intro h3; apply hq
        linear_combination -hT + (x₁ - e) * h3 + (x₁ - e) * hx2
      rw [if_neg hx1, if_pos hx2, if_neg h3]
      exact ⟨qe e a b, by linear_combination (qe e a b) * hT - (qe e a b) * (x₁ - e) * hx2⟩
    rw [if_neg hx1, if_neg hx2]
    by_cases h3 : ℓ ^ 2 - (a - e) - x₁ - x₂ = e
    · rw [if_pos h3]
      have hL : ℓ * (e - x₁) + y₁ = 0 := by
        have h0 : (ℓ * (e - x₁) + y₁) ^ 2 = 0 := by
          linear_combination -hC + (x₁ - e) * (x₂ - e) * h3
        exact (pow_eq_zero_iff two_ne_zero).mp h0
      have h := chord_third_T_coord e₁ hx1 hL h3
      exact ⟨qe e a b, by linear_combination (qe e a b) * h⟩
    · rw [if_neg h3]
      exact ⟨ℓ * (e - x₁) + y₁, by linear_combination hC⟩

/-- **The kernel of `α` contains `2E(K)`**: `α(W + W)` is a square for every point `W`. -/
theorem alpha_double (hq : qe e a b ≠ 0) (W₀ : (curve e a b).toAffine.Point) :
    IsSquare (alpha e a b (W₀ + W₀)) := by
  have h := alpha_add hq W₀ W₀
  rw [mul_assoc] at h
  exact isSquare_cancel ⟨_, rfl⟩ (mul_ne_zero (alpha_ne_zero hq W₀) (alpha_ne_zero hq W₀)) h

/-- **Step 3, doubles.** If `P = W + W` is affine with `x(P) ≠ e` then `x(P) − e` is a square in `K`. -/
theorem double_x_sub_isSquare (hq : qe e a b ≠ 0) {W₀ : (curve e a b).toAffine.Point} {x y : K}
    {h : (curve e a b).toAffine.Nonsingular x y} (hP : W₀ + W₀ = .some x y h) (hx : x ≠ e) :
    IsSquare (x - e) := by
  have := alpha_double hq W₀
  rwa [hP, alpha_some, if_neg hx] at this

/-- **Step 3, `T`.** If `T = (e, 0)` is `W + W` then `q(e)` is a square in `K`. -/
theorem double_T_isSquare (hq : qe e a b ≠ 0) {W₀ : (curve e a b).toAffine.Point} {y : K}
    {h : (curve e a b).toAffine.Nonsingular e y} (hP : W₀ + W₀ = .some e y h) :
    IsSquare (qe e a b) := by
  have := alpha_double hq W₀
  rwa [hP, alpha_some, if_pos rfl] at this

/-- **Step 3, the product rule for three collinear points**, on Mathlib's points: if `P₁ + P₂ = P₃` with
`x₁ ≠ x₂` (so `P₁, P₂, −P₃` lie on the line `y = ℓ(x − x₁) + y₁`), then
`(x₁ − e)(x₂ − e)(x₃ − e) = (ℓ(e − x₁) + y₁)²`. -/
theorem chord_point {x₁ y₁ x₂ y₂ x₃ y₃ : K} {h₁ : (curve e a b).toAffine.Nonsingular x₁ y₁}
    {h₂ : (curve e a b).toAffine.Nonsingular x₂ y₂} {h₃ : (curve e a b).toAffine.Nonsingular x₃ y₃}
    (hx : x₁ ≠ x₂) (hP : Point.some x₁ y₁ h₁ + Point.some x₂ y₂ h₂ = Point.some x₃ y₃ h₃) :
    (x₁ - e) * (x₂ - e) * (x₃ - e) =
      ((curve e a b).toAffine.slope x₁ x₂ y₁ y₂ * (e - x₁) + y₁) ^ 2 := by
  rw [Point.add_of_X_ne hx] at hP
  obtain ⟨hx3, -⟩ := Point.some.inj hP
  rw [← hx3]
  have e₁ := (curve_equation_iff x₁ y₁).mp h₁.1
  have e₂ := (curve_equation_iff x₂ y₂).mp h₂.1
  have hs : (curve e a b).toAffine.slope x₁ x₂ y₁ y₂ * (x₁ - x₂) = y₁ - y₂ := by
    rw [Affine.slope_of_X_ne hx, div_mul_cancel₀ _ (sub_ne_zero.mpr hx)]
  have hC := chord_coord e₁ e₂ hx hs
  simp only [Affine.addX, curve_a₁, curve_a₂]
  linear_combination hC

/-- **Step 3, doubles, without hypotheses**: if `W + W = (x, y)` is affine then `x − e` is a square in `K`
(`x − e = (((x₀ − e)² − q(e))/(2y₀))²` for `W = (x₀, y₀)`). No condition on `q(e)` or on the characteristic. -/
theorem double_isSquare {W₀ : (curve e a b).toAffine.Point} {x y : K}
    {h : (curve e a b).toAffine.Nonsingular x y} (hP : W₀ + W₀ = .some x y h) : IsSquare (x - e) := by
  rcases W₀ with _ | ⟨x₀, y₀, h₀⟩
  · exact absurd hP.symm (Point.some_ne_zero h)
  by_cases hy : y₀ = (curve e a b).toAffine.negY x₀ y₀
  · rw [Point.add_self_of_Y_eq hy] at hP
    exact absurd hP.symm (Point.some_ne_zero h)
  rw [Point.add_self_of_Y_ne hy] at hP
  obtain ⟨hx3, -⟩ := Point.some.inj hP
  obtain ⟨hs, h2⟩ := slope_double hy
  have e₀ := (curve_equation_iff x₀ y₀).mp h₀.1
  have key := tangent_coord e₀ hs
  obtain ⟨w, hw⟩ : ∃ w : K, 2 * y₀ * w = 1 := ⟨(2 * y₀)⁻¹, mul_inv_cancel₀ h2⟩
  refine ⟨((x₀ - e) ^ 2 - qe e a b) * w, ?_⟩
  rw [← hx3]
  simp only [Affine.addX, curve_a₁, curve_a₂]
  linear_combination (-((curve e a b).toAffine.slope x₀ x₀ y₀ y₀ ^ 2 - (a - e) - x₀ - x₀ - e) *
    (1 + 2 * y₀ * w)) * hw + w ^ 2 * key

/-- **Step 3, the halving test at the conjugate roots**, over `K`: if `a = −2r`, `b = r² + s²` and `P = W + W`
is affine, then `x(P) − (r + is)` is a square in `K(i)`, i.e. `x(P) − r = u² − v²`, `−s = 2uv` for some
`u, v ∈ K`. -/
theorem gauss_double {r s : K} (ha : a = -2 * r) (hb : b = r ^ 2 + s ^ 2)
    {W₀ : (curve e a b).toAffine.Point} {x y : K} {h : (curve e a b).toAffine.Nonsingular x y}
    (hP : W₀ + W₀ = .some x y h) : ∃ u v : K, x - r = u ^ 2 - v ^ 2 ∧ -s = 2 * u * v := by
  rcases W₀ with _ | ⟨x₀, y₀, h₀⟩
  · exact absurd hP.symm (Point.some_ne_zero h)
  by_cases hy : y₀ = (curve e a b).toAffine.negY x₀ y₀
  · rw [Point.add_self_of_Y_eq hy] at hP
    exact absurd hP.symm (Point.some_ne_zero h)
  rw [Point.add_self_of_Y_ne hy] at hP
  obtain ⟨hx3, -⟩ := Point.some.inj hP
  obtain ⟨hs, h2⟩ := slope_double hy
  have e₀ := (curve_equation_iff x₀ y₀).mp h₀.1
  obtain ⟨u, v, h1, h2'⟩ := gauss_coord ha hb e₀ h2 hs
  refine ⟨u, v, ?_, h2'⟩
  rw [← hx3]
  simp only [Affine.addX, curve_a₁, curve_a₂]
  linear_combination h1

/-! ### Step 4: `α` as a homomorphism `E(K) → Kˣ/(Kˣ)²` -/

/-- The group of square classes `Kˣ/(Kˣ)²`. -/
abbrev SqClass (K : Type*) [Field K] := Kˣ ⧸ (powMonoidHom 2 : Kˣ →* Kˣ).range

omit [DecidableEq K] in
theorem sqClass_mk_eq_one {u : Kˣ} : (QuotientGroup.mk u : SqClass K) = 1 ↔ IsSquare (u : K) := by
  rw [QuotientGroup.eq_one_iff, MonoidHom.mem_range]
  constructor
  · rintro ⟨v, hv⟩
    exact ⟨v, by rw [← hv]; simp [sq]⟩
  · rintro ⟨r, hr⟩
    have hr0 : r ≠ 0 := by rintro rfl; exact u.ne_zero (by rw [hr, mul_zero])
    exact ⟨Units.mk0 r hr0, by ext; simp [sq, hr]⟩

/-- **The 2-descent map as a group homomorphism** `α : E(K) → Kˣ/(Kˣ)²` (written additively), for
`y² = (X − e)(X² + aX + b)` with `q(e) ≠ 0`: `α(O) = 1`, `α(e, 0) = q(e)`, `α(x, y) = x − e` otherwise. -/
noncomputable def alphaHom (hq : qe e a b ≠ 0) :
    (curve e a b).toAffine.Point →+ Additive (SqClass K) where
  toFun P := Additive.ofMul (QuotientGroup.mk (Units.mk0 (alpha e a b P) (alpha_ne_zero hq P)))
  map_zero' := by
    rw [ofMul_eq_zero, sqClass_mk_eq_one]
    exact ⟨1, by simp [alpha_zero]⟩
  map_add' P₁ P₂ := by
    have hsq : ∀ u : Kˣ, Additive.ofMul (QuotientGroup.mk u : SqClass K) +
        Additive.ofMul (QuotientGroup.mk u : SqClass K) = 0 := by
      intro u
      rw [← ofMul_mul, ← QuotientGroup.mk_mul, ofMul_eq_zero, sqClass_mk_eq_one]
      exact ⟨u, by simp⟩
    have h3 : Additive.ofMul (QuotientGroup.mk (Units.mk0 (alpha e a b (P₁ + P₂))
          (alpha_ne_zero hq _)) : SqClass K) +
        Additive.ofMul (QuotientGroup.mk (Units.mk0 (alpha e a b P₁) (alpha_ne_zero hq _))) +
        Additive.ofMul (QuotientGroup.mk (Units.mk0 (alpha e a b P₂) (alpha_ne_zero hq _))) = 0 := by
      rw [← ofMul_mul, ← ofMul_mul, ← QuotientGroup.mk_mul, ← QuotientGroup.mk_mul, ofMul_eq_zero,
        sqClass_mk_eq_one]
      simpa using alpha_add hq P₁ P₂
    calc _ = (Additive.ofMul (QuotientGroup.mk (Units.mk0 (alpha e a b (P₁ + P₂))
              (alpha_ne_zero hq _)) : SqClass K) +
            Additive.ofMul (QuotientGroup.mk (Units.mk0 (alpha e a b P₁) (alpha_ne_zero hq _))) +
            Additive.ofMul (QuotientGroup.mk (Units.mk0 (alpha e a b P₂) (alpha_ne_zero hq _)))) +
          (Additive.ofMul (QuotientGroup.mk (Units.mk0 (alpha e a b P₁) (alpha_ne_zero hq _))) +
            Additive.ofMul (QuotientGroup.mk (Units.mk0 (alpha e a b P₂) (alpha_ne_zero hq _)))) -
          (Additive.ofMul (QuotientGroup.mk (Units.mk0 (alpha e a b P₁) (alpha_ne_zero hq _))) +
            Additive.ofMul (QuotientGroup.mk (Units.mk0 (alpha e a b P₁) (alpha_ne_zero hq _)))) -
          (Additive.ofMul (QuotientGroup.mk (Units.mk0 (alpha e a b P₂) (alpha_ne_zero hq _))) +
            Additive.ofMul (QuotientGroup.mk (Units.mk0 (alpha e a b P₂) (alpha_ne_zero hq _)))) := by
          abel
      _ = _ := by rw [h3, hsq, hsq]; abel

theorem alphaHom_eq_zero_iff (hq : qe e a b ≠ 0) (P : (curve e a b).toAffine.Point) :
    alphaHom hq P = 0 ↔ IsSquare (alpha e a b P) := by
  change Additive.ofMul _ = 0 ↔ _
  rw [ofMul_eq_zero, sqClass_mk_eq_one, Units.val_mk0]

/-- **`2E(K)` lies in the kernel of the homomorphism `α`.** -/
theorem alphaHom_two_smul (hq : qe e a b ≠ 0) (W₀ : (curve e a b).toAffine.Point) :
    alphaHom hq ((2 : ℤ) • W₀) = 0 := by
  rw [two_zsmul, alphaHom_eq_zero_iff]
  exact alpha_double hq W₀

/-- **The independence argument of `UNIFORM_RANK.md`, Theorem 2, abstractly.** Let `q(e) ≠ 0` and let `R, Q, T`
be points with `α(R)`, `α(T)`, `α(R)·α(T)` not squares in `K`, `α(Q)` a square, and `Q` not of the form `W + W`.
Then no non-trivial combination `ε₁R + ε₂Q + ε₃T` (`εᵢ ∈ {0, 1}`) is of the form `W + W`: the three points are
independent in `E(K)/2E(K)`. -/
theorem indep_of_alpha (hq : qe e a b ≠ 0) (R Q T : (curve e a b).toAffine.Point)
    (nR : ¬ IsSquare (alpha e a b R)) (nT : ¬ IsSquare (alpha e a b T))
    (nRT : ¬ IsSquare (alpha e a b R * alpha e a b T)) (sQ : IsSquare (alpha e a b Q))
    (hQ : ∀ W : (curve e a b).toAffine.Point, Q ≠ W + W) (ε₁ ε₂ ε₃ : ℕ) (h₁ : ε₁ < 2) (h₂ : ε₂ < 2)
    (h₃ : ε₃ < 2) (W : (curve e a b).toAffine.Point) (hW : ε₁ • R + ε₂ • Q + ε₃ • T = W + W) :
    ε₁ = 0 ∧ ε₂ = 0 ∧ ε₃ = 0 := by
  have ne := alpha_ne_zero hq
  have hd : ∀ P : (curve e a b).toAffine.Point, P = W + W → IsSquare (alpha e a b P) := by
    rintro P rfl; exact alpha_double hq W
  have c₁ : ε₁ = 0 ∨ ε₁ = 1 := by omega
  have c₂ : ε₂ = 0 ∨ ε₂ = 1 := by omega
  have c₃ : ε₃ = 0 ∨ ε₃ = 1 := by omega
  rcases c₁ with rfl | rfl <;> rcases c₂ with rfl | rfl <;> rcases c₃ with rfl | rfl <;>
    simp only [zero_nsmul, one_nsmul, zero_add, add_zero] at hW
  · exact ⟨rfl, rfl, rfl⟩
  · exact absurd (hd T hW) nT
  · exact absurd hW (hQ W)
  · have m := alpha_add hq Q T
    rw [mul_assoc] at m
    exact absurd (isSquare_cancel_left sQ (ne Q) (isSquare_cancel_left (hd _ hW) (ne _) m)) nT
  · exact absurd (hd R hW) nR
  · have m := alpha_add hq R T
    rw [mul_assoc] at m
    exact absurd (isSquare_cancel_left (hd _ hW) (ne _) m) nRT
  · have m := alpha_add hq R Q
    rw [mul_assoc] at m
    have m' := isSquare_cancel_left (hd _ hW) (ne _) m
    rw [mul_comm] at m'
    exact absurd (isSquare_cancel_left sQ (ne Q) m') nR
  · have m₁ := alpha_add hq (R + Q) T
    rw [mul_assoc] at m₁
    have A : IsSquare (alpha e a b (R + Q) * alpha e a b T) :=
      isSquare_cancel_left (hd _ hW) (ne _) m₁
    have B : IsSquare (alpha e a b (R + Q) * alpha e a b R * alpha e a b Q) := alpha_add hq R Q
    have C : IsSquare ((alpha e a b (R + Q) * alpha e a b (R + Q)) *
        (alpha e a b Q * (alpha e a b R * alpha e a b T))) := by
      obtain ⟨r, hr⟩ := A
      obtain ⟨t, ht⟩ := B
      refine ⟨r * t, ?_⟩
      linear_combination (alpha e a b (R + Q) * alpha e a b R * alpha e a b Q) * hr + (r * r) * ht
    have D := isSquare_cancel_left ⟨_, rfl⟩ (mul_ne_zero (ne _) (ne _)) C
    exact absurd (isSquare_cancel_left sQ (ne Q) D) nRT

/-- The same with integer coefficients: if `aR + bQ + cT = W + W` then `a`, `b`, `c` are even. -/
theorem indep_of_alpha_int (hq : qe e a b ≠ 0) (R Q T : (curve e a b).toAffine.Point)
    (nR : ¬ IsSquare (alpha e a b R)) (nT : ¬ IsSquare (alpha e a b T))
    (nRT : ¬ IsSquare (alpha e a b R * alpha e a b T)) (sQ : IsSquare (alpha e a b Q))
    (hQ : ∀ W : (curve e a b).toAffine.Point, Q ≠ W + W) (c₁ c₂ c₃ : ℤ)
    (W : (curve e a b).toAffine.Point) (hW : c₁ • R + c₂ • Q + c₃ • T = W + W) :
    Even c₁ ∧ Even c₂ ∧ Even c₃ := by
  have key : ∀ c : ℤ, ∃ ε : ℕ, ε < 2 ∧ ∃ d : ℤ, c = (ε : ℤ) + d + d ∧ (Even c ↔ ε = 0) := by
    intro c
    rcases Int.even_or_odd' c with ⟨d, hd | hd⟩
    · exact ⟨0, by norm_num, d, by rw [hd]; push_cast; ring,
        ⟨fun _ => rfl, fun _ => ⟨d, by rw [hd]; ring⟩⟩⟩
    · exact ⟨1, by norm_num, d, by rw [hd]; push_cast; ring,
        ⟨fun h => absurd h (Int.not_even_iff_odd.mpr ⟨d, hd⟩), fun h => absurd h one_ne_zero⟩⟩
  obtain ⟨ε₁, h₁, d₁, e₁, i₁⟩ := key c₁
  obtain ⟨ε₂, h₂, d₂, e₂, i₂⟩ := key c₂
  obtain ⟨ε₃, h₃, d₃, e₃, i₃⟩ := key c₃
  have hW' : ε₁ • R + ε₂ • Q + ε₃ • T = (W - (d₁ • R + d₂ • Q + d₃ • T)) +
      (W - (d₁ • R + d₂ • Q + d₃ • T)) := by
    rw [e₁, e₂, e₃] at hW
    simp only [add_zsmul, natCast_zsmul] at hW
    rw [← sub_eq_zero] at hW ⊢
    rw [← hW]
    abel
  obtain ⟨z₁, z₂, z₃⟩ := indep_of_alpha hq R Q T nR nT nRT sQ hQ ε₁ ε₂ ε₃ h₁ h₂ h₃ _ hW'
  exact ⟨i₁.mpr z₁, i₂.mpr z₂, i₃.mpr z₃⟩

end General

/-! ## Non-squares in `ℚ` -/

/-- **Square-root bracket.** If `s² < N < (s + 1)²` then `N` is not a square in `ℚ`. -/
theorem not_isSquare_of_bracket (N s : ℕ) (h1 : s * s < N) (h2 : N < (s + 1) * (s + 1)) :
    ¬ IsSquare (N : ℚ) := by
  rw [Rat.isSquare_natCast_iff]
  rintro ⟨t, rfl⟩
  rcases le_or_gt t s with h | h
  · exact absurd h1 (not_lt.mpr (Nat.mul_self_le_mul_self h))
  · exact absurd h2 (not_lt.mpr (Nat.mul_self_le_mul_self h))

theorem not_isSquare_of_neg {x : ℚ} (hx : x < 0) : ¬ IsSquare x := by
  rintro ⟨r, rfl⟩
  nlinarith [mul_self_nonneg r]

/-- The bracket test for a rational number given as an expression `x` equal to the natural number `N`. -/
theorem nsq {x : ℚ} (N s : ℕ) (h : x = (N : ℚ)) (h1 : s * s < N) (h2 : N < (s + 1) * (s + 1)) :
    ¬ IsSquare x := by
  rw [h]; exact not_isSquare_of_bracket N s h1 h2

/-! ## The curves `E±` of `UNIFORM_RANK.md` -/

section Family

open MagicSquareOfSquares.G1

/-- `e = −2S³`, `S = m² + n²`: the rational root of the cubic of `E⁺`. -/
def eP (m n : ℚ) : ℚ := -2 * (m ^ 2 + n ^ 2) ^ 3
/-- `a = 16m²n²(m + n)²`. -/
def aP (m n : ℚ) : ℚ := 16 * m ^ 2 * n ^ 2 * (m + n) ^ 2
/-- `b = 16m²n²(m + n)⁴S²`. -/
def bP (m n : ℚ) : ℚ := 16 * m ^ 2 * n ^ 2 * (m + n) ^ 4 * (m ^ 2 + n ^ 2) ^ 2
/-- `E⁺` at `(m, n)`: `y² = (X + 2S³)(X² + 16m²n²(m+n)²X + 16m²n²(m+n)⁴S²)`. `E⁻` at `(m, n)` is `Eplus m (−n)`
(`UNIFORM_RANK.md`, A2, checked by sympy; on the quartics this is `gen_minus_plus`). -/
abbrev Eplus (m n : ℚ) : WeierstrassCurve ℚ := curve (eP m n) (aP m n) (bP m n)

/-- The cubic of `Eplus m n` is the cubic `genCub m n` of `G1SigmaGeneric.lean`, which is the model
`sigma_exact.model(m, n, 2)`. -/
theorem Eplus_cubic (m n X : ℚ) :
    (X - eP m n) * (X ^ 2 + aP m n * X + bP m n) = genCub m n X := by
  rw [gen_cubic_factor]; unfold eP aP bP; ring

/-- `x(R) = −4mnS(m² − n²)`. -/
def xR (m n : ℚ) : ℚ := -4 * m * n * (m ^ 2 + n ^ 2) * (m ^ 2 - n ^ 2)
/-- `y(R) = −8mn(m + n)S·H(m, n)`. -/
def yR (m n : ℚ) : ℚ := -8 * m * n * (m + n) * (m ^ 2 + n ^ 2) * genH m n
/-- `y(Q) = 8kmn(m + n)²S²` for `S = 2k²` (`Q = (0, y(Q))`). -/
def yQ (m n k : ℚ) : ℚ := 8 * k * m * n * (m + n) ^ 2 * (m ^ 2 + n ^ 2) ^ 2
/-- The quadratic factor is `(X − r)² + s²` with `r = −8m²n²(m + n)²` ... -/
def rP (m n : ℚ) : ℚ := -8 * m ^ 2 * n ^ 2 * (m + n) ^ 2
/-- ... and `s = 4mn(m + n)³(m − n)`; its roots are `θ, θ̄ = r ± is`. -/
def sP (m n : ℚ) : ℚ := 4 * m * n * (m + n) ^ 3 * (m - n)

theorem aP_eq (m n : ℚ) : aP m n = -2 * rP m n := by unfold aP rP; ring

theorem bP_eq (m n : ℚ) : bP m n = rP m n ^ 2 + sP m n ^ 2 := by unfold bP rP sP; ring

theorem qe_P (m n : ℚ) :
    qe (eP m n) (aP m n) (bP m n) = 4 * (m ^ 2 + n ^ 2) ^ 2 * (genH m n * genH m (-n)) := by
  unfold qe eP aP bP genH; ring

theorem R_on (m n : ℚ) :
    yR m n ^ 2 = (xR m n - eP m n) * (xR m n ^ 2 + aP m n * xR m n + bP m n) := by
  unfold yR xR eP aP bP genH; ring

theorem alphaR_val (m n : ℚ) : xR m n - eP m n = 2 * (m ^ 2 + n ^ 2) * genH m n := by
  unfold xR eP genH; ring

theorem Q_on {m n k : ℚ} (hk : m ^ 2 + n ^ 2 = 2 * k ^ 2) :
    yQ m n k ^ 2 = (0 - eP m n) * (0 ^ 2 + aP m n * 0 + bP m n) := by
  unfold yQ eP aP bP
  linear_combination (-32 * m ^ 2 * n ^ 2 * (m + n) ^ 4 * (m ^ 2 + n ^ 2) ^ 4) * hk

/-- `H(m, n) = (m² − n² − mn)² + 3(mn)² > 0` for `m ≠ 0`. -/
theorem genH_pos {m n : ℚ} (hm : m ≠ 0) : 0 < genH m n := by
  have h : genH m n = (m ^ 2 - n ^ 2 - m * n) ^ 2 + 3 * (m * n) ^ 2 := by unfold genH; ring
  rw [h]
  have hm2 : 0 < m ^ 2 := lt_of_le_of_ne (sq_nonneg m) (Ne.symm (pow_ne_zero 2 hm))
  rcases eq_or_ne n 0 with hn | hn
  · subst hn
    have h2 : (m ^ 2 - 0 ^ 2 - m * 0) ^ 2 + 3 * (m * 0) ^ 2 = (m ^ 2) ^ 2 := by ring
    rw [h2]; exact pow_pos hm2 2
  · have : 0 < (m * n) ^ 2 :=
      lt_of_le_of_ne (sq_nonneg _) (Ne.symm (pow_ne_zero 2 (mul_ne_zero hm hn)))
    nlinarith [sq_nonneg (m ^ 2 - n ^ 2 - m * n)]

theorem Eplus_Δ (m n : ℚ) :
    (Eplus m n).Δ = -64 * qe (eP m n) (aP m n) (bP m n) ^ 2 * sP m n ^ 2 := by
  rw [curve_Δ]
  have h : aP m n ^ 2 - 4 * bP m n = -4 * sP m n ^ 2 := by rw [aP_eq, bP_eq]; ring
  rw [h]; ring

/-- **`E⁺` is an elliptic curve** (`Δ ≠ 0`) whenever `m, n, m ± n ≠ 0`. -/
theorem Eplus_isElliptic {m n : ℚ} (hm : m ≠ 0) (hn : n ≠ 0) (hp : m + n ≠ 0) (hd : m - n ≠ 0) :
    (Eplus m n).IsElliptic := by
  have hS : m ^ 2 + n ^ 2 ≠ 0 := by
    have h1 : 0 < m ^ 2 := lt_of_le_of_ne (sq_nonneg m) (Ne.symm (pow_ne_zero 2 hm))
    have h2 := sq_nonneg n
    exact (by linarith : 0 < m ^ 2 + n ^ 2).ne'
  have hm' : m ≠ 0 := hm
  have hq : qe (eP m n) (aP m n) (bP m n) ≠ 0 := by
    rw [qe_P]
    exact mul_ne_zero (mul_ne_zero (by norm_num) (pow_ne_zero 2 hS))
      (mul_ne_zero (genH_pos hm).ne' (genH_pos hm').ne')
  have hs : sP m n ≠ 0 := by
    unfold sP
    exact mul_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero (by norm_num) hm) hn)
      (pow_ne_zero 3 hp)) hd
  refine ⟨isUnit_iff_ne_zero.mpr ?_⟩
  rw [Eplus_Δ]
  exact mul_ne_zero (mul_ne_zero (by norm_num) (pow_ne_zero 2 hq)) (pow_ne_zero 2 hs)

/-- `x(P) = −(m + n)²S²`, the point of `z = 2`. -/
def xP (m n : ℚ) : ℚ := -((m + n) ^ 2 * (m ^ 2 + n ^ 2) ^ 2)
/-- `y(P) = S³(m − n)(m + n)²`. -/
def yP (m n : ℚ) : ℚ := (m ^ 2 + n ^ 2) ^ 3 * (m - n) * (m + n) ^ 2

theorem P_on (m n : ℚ) :
    yP m n ^ 2 = (xP m n - eP m n) * (xP m n ^ 2 + aP m n * xP m n + bP m n) := by
  unfold yP xP eP aP bP; ring

/-- **`P = 2R` on Mathlib's points, for every ratio** with `m, n, m ± n ≠ 0` (`UNIFORM_RANK.md` section 2, A5;
the identities are `gen_tangent_x`, `gen_tangent_y` of `G1SigmaGeneric.lean`). Then `α(P) = x(P) − e =
(S(m − n))²` is a square, as `alpha_double` says it must be. -/
theorem R_add_R {m n : ℚ} (hm : m ≠ 0) (hn : n ≠ 0) (hp : m + n ≠ 0) (hd : m - n ≠ 0) :
    ∃ (hR : (Eplus m n).toAffine.Nonsingular (xR m n) (yR m n))
      (hP : (Eplus m n).toAffine.Nonsingular (xP m n) (yP m n)),
      Point.some _ _ hR + Point.some _ _ hR = Point.some _ _ hP := by
  have hS : m ^ 2 + n ^ 2 ≠ 0 := by
    have h1 : 0 < m ^ 2 := lt_of_le_of_ne (sq_nonneg m) (Ne.symm (pow_ne_zero 2 hm))
    have h2 := sq_nonneg n
    exact (by linarith : 0 < m ^ 2 + n ^ 2).ne'
  have hyR : 2 * yR m n ≠ 0 := by
    unfold yR
    exact mul_ne_zero two_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero
      (mul_ne_zero (by norm_num) hm) hn) hp) hS) (genH_pos hm).ne')
  have hyP : 2 * yP m n ≠ 0 := by
    unfold yP
    exact mul_ne_zero two_ne_zero (mul_ne_zero (mul_ne_zero (pow_ne_zero 3 hS) hd) (pow_ne_zero 2 hp))
  have hR : (Eplus m n).toAffine.Nonsingular (xR m n) (yR m n) := nonsingular_of_y (R_on m n) hyR
  have hP : (Eplus m n).toAffine.Nonsingular (xP m n) (yP m n) := nonsingular_of_y (P_on m n) hyP
  refine ⟨hR, hP, ?_⟩
  have hy : yR m n ≠ (Eplus m n).toAffine.negY (xR m n) (yR m n) := by
    rw [curve_negY]; intro h; apply hyR; linear_combination h
  have hF : ∀ x : ℚ, 3 * x ^ 2 + 2 * (aP m n - eP m n) * x + (bP m n - aP m n * eP m n) =
      genCubD m n x := by
    intro x; unfold genCubD genA2 genA genC3 aP eP bP; ring
  have hA2 : genA2 m n = aP m n - eP m n := by unfold genA2 aP eP; ring
  have key : (3 * xR m n ^ 2 + 2 * (aP m n - eP m n) * xR m n + (bP m n - aP m n * eP m n)) ^ 2 =
      4 * yR m n ^ 2 * (xP m n + (aP m n - eP m n) + 2 * xR m n) := by
    have h := gen_tangent_x m n
    rw [← hF, hA2] at h
    exact h
  have key2 : 2 * yR m n * yP m n + 2 * yR m n ^ 2 +
      (3 * xR m n ^ 2 + 2 * (aP m n - eP m n) * xR m n + (bP m n - aP m n * eP m n)) *
        (xP m n - xR m n) = 0 := by
    have h := gen_tangent_y m n
    rw [← hF] at h
    exact h
  obtain ⟨hs, -⟩ := slope_double hy
  rw [Point.add_self_of_Y_ne hy, Point.some.injEq]
  have hX : (Eplus m n).toAffine.addX (xR m n) (xR m n)
      ((Eplus m n).toAffine.slope (xR m n) (xR m n) (yR m n) (yR m n)) = xP m n := by
    simp only [Affine.addX, curve_a₁, curve_a₂]
    have h4 : (2 * yR m n) ^ 2 * ((Eplus m n).toAffine.slope (xR m n) (xR m n) (yR m n) (yR m n) ^ 2
        - (aP m n - eP m n) - 2 * xR m n - xP m n) = 0 := by
      linear_combination ((Eplus m n).toAffine.slope (xR m n) (xR m n) (yR m n) (yR m n) * (2 * yR m n)
        + (3 * xR m n ^ 2 + 2 * (aP m n - eP m n) * xR m n + (bP m n - aP m n * eP m n))) * hs + key
    have h5 := (mul_eq_zero.mp h4).resolve_left (pow_ne_zero 2 hyR)
    linear_combination h5
  refine ⟨hX, ?_⟩
  simp only [Affine.addY, curve_negY, Affine.negAddY]
  rw [hX]
  have h2 : (2 * yR m n) * (-((Eplus m n).toAffine.slope (xR m n) (xR m n) (yR m n) (yR m n) *
      (xP m n - xR m n) + yR m n) - yP m n) = 0 := by
    linear_combination (-(xP m n - xR m n)) * hs - key2
  linear_combination (mul_eq_zero.mp h2).resolve_left hyR

/-- The halving test at `θ` for the point `Q` (`x(Q) = 0`): if `−θ = (u + iv)²` with `u, v ∈ ℚ` then `2mn` or
`−2mn` is a square in `ℚ` (`UNIFORM_RANK.md`: `−θ = (2i(m + n)(m + in))²·(i·mn)`). -/
theorem gauss_family {m n u v : ℚ} (hp : m + n ≠ 0) (hd : m - n ≠ 0)
    (h1 : 0 - rP m n = u ^ 2 - v ^ 2) (h2 : -sP m n = 2 * u * v) :
    IsSquare (2 * m * n) ∨ IsSquare (-2 * m * n) := by
  unfold rP at h1
  unfold sP at h2
  have h : (u ^ 2 + v ^ 2) ^ 2 = (4 * m * n * (m + n) ^ 2 * (m ^ 2 + n ^ 2)) ^ 2 := by
    linear_combination (-(u ^ 2 - v ^ 2 + 8 * m ^ 2 * n ^ 2 * (m + n) ^ 2)) * h1
      + (-(2 * u * v - 4 * m * n * (m + n) ^ 3 * (m - n))) * h2
  rcases sq_eq_sq_iff_eq_or_eq_neg.mp h with h3 | h3
  · left
    have hu : u ^ 2 = 2 * m * n * (m + n) ^ 4 := by
      linear_combination (-1 / 2 : ℚ) * h1 + (1 / 2 : ℚ) * h3
    obtain ⟨w, hw⟩ : ∃ w : ℚ, (m + n) ^ 2 * w = 1 := ⟨((m + n) ^ 2)⁻¹, mul_inv_cancel₀ (pow_ne_zero 2 hp)⟩
    exact ⟨u * w, by linear_combination (-2 * m * n * (1 + (m + n) ^ 2 * w)) * hw - w ^ 2 * hu⟩
  · right
    have hu : u ^ 2 = -2 * m * n * (m + n) ^ 2 * (m - n) ^ 2 := by
      linear_combination (-1 / 2 : ℚ) * h1 + (1 / 2 : ℚ) * h3
    obtain ⟨w, hw⟩ : ∃ w : ℚ, ((m + n) * (m - n)) * w = 1 :=
      ⟨((m + n) * (m - n))⁻¹, mul_inv_cancel₀ (mul_ne_zero hp hd)⟩
    exact ⟨u * w, by
      linear_combination (2 * m * n * (1 + ((m + n) * (m - n)) * w)) * hw - w ^ 2 * hu⟩

/-- **2-torsion over `ℚ`.** If `a = −2r`, `b = r² + s²` with `s ≠ 0` (the quadratic factor has no rational root),
then an affine point `P` with `P + P = O` is `T = (e, 0)`. For `Eplus m n`, `s = sP m n ≠ 0` when
`m, n, m ± n ≠ 0` (`two_torsion_Eplus`), so `E⁺(ℚ)[2] ⊆ {O, T}` (and `T` is in it). -/
theorem two_torsion_Q {e a b r s : ℚ} (ha : a = -2 * r) (hb : b = r ^ 2 + s ^ 2) (hs : s ≠ 0) {x y : ℚ}
    {h : (curve e a b).toAffine.Nonsingular x y} (hP : Point.some x y h + Point.some x y h = 0) :
    x = e ∧ y = 0 := by
  have hy : y = (curve e a b).toAffine.negY x y := by
    by_contra hy
    rw [Point.add_self_of_Y_ne hy] at hP
    exact Point.some_ne_zero _ hP
  rw [curve_negY] at hy
  have hy0 : y = 0 := by linarith
  have e₀ := (curve_equation_iff x y).mp h.1
  have hpos : 0 < (x - r) ^ 2 + s ^ 2 := by
    have : 0 < s ^ 2 := lt_of_le_of_ne (sq_nonneg s) (Ne.symm (pow_ne_zero 2 hs))
    nlinarith [sq_nonneg (x - r)]
  have h0 : (x - e) * ((x - r) ^ 2 + s ^ 2) = 0 := by
    rw [hy0, ha, hb] at e₀; linear_combination -e₀
  refine ⟨?_, hy0⟩
  rcases mul_eq_zero.mp h0 with h1 | h1
  · linarith
  · linarith

theorem two_torsion_Eplus {m n : ℚ} (hm : m ≠ 0) (hn : n ≠ 0) (hp : m + n ≠ 0) (hd : m - n ≠ 0) {x y : ℚ}
    {h : (Eplus m n).toAffine.Nonsingular x y} (hP : Point.some x y h + Point.some x y h = 0) :
    x = eP m n ∧ y = 0 :=
  two_torsion_Q (aP_eq m n) (bP_eq m n)
    (by unfold sP; exact mul_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero (by norm_num) hm) hn)
      (pow_ne_zero 3 hp)) hd) hP

/-- The conclusion of `UNIFORM_RANK.md` Theorem 2 on `Eplus m n` (`k` gives `y(Q)`): the three affine points
`R`, `Q`, `T` are nonsingular, and `c₁R + c₂Q + c₃T = W + W` with `cᵢ ∈ ℤ`, `W ∈ E(ℚ)` forces `c₁, c₂, c₃`
even. -/
def Theorem2At (m n k : ℚ) : Prop :=
  ∃ (hR : (Eplus m n).toAffine.Nonsingular (xR m n) (yR m n))
    (hQ : (Eplus m n).toAffine.Nonsingular 0 (yQ m n k))
    (hT : (Eplus m n).toAffine.Nonsingular (eP m n) 0),
    ∀ c₁ c₂ c₃ : ℤ, ∀ W : (Eplus m n).toAffine.Point,
      c₁ • Point.some _ _ hR + c₂ • Point.some _ _ hQ + c₃ • Point.some _ _ hT = W + W →
        Even c₁ ∧ Even c₂ ∧ Even c₃

/-- **`UNIFORM_RANK.md` Theorem 2 at one ratio, from five non-square facts.**
Hypotheses, all on rational numbers `m, n, k`: `m² + n² = 2k²`; `m, n, m + n, m − n ≠ 0`; none of `H(m, n)`,
`H(m, −n)`, `H(m, n)·H(m, −n)`, `2mn`, `−2mn` is a square in `ℚ`. Conclusion: on `E⁺` at `(m, n)` (the curve
`y² = (X + 2S³)(X² + 16m²n²(m+n)²X + 16m²n²(m+n)⁴S²)`, `S = m² + n²`), the points
`R = (−4mnS(m² − n²), −8mn(m+n)S·H(m, n))`, `Q = (0, 8kmn(m+n)²S²)` and `T = (−2S³, 0)` are nonsingular points,
and if `c₁R + c₂Q + c₃T = W + W` for integers `cᵢ` and a rational point `W`, then `c₁, c₂, c₃` are even: `R, Q, T`
are independent in `E⁺(ℚ)/2E⁺(ℚ)`. (`E⁻` at `(m, n)` is the case `(m, −n)`.)
`E⁺(ℚ)[2] ⊆ {O, T}` is `two_torsion_Eplus`; "rank `≥ 2`" follows by Mordell–Weil and
`dim E(ℚ)/2E(ℚ) = rank + dim E(ℚ)[2]`, which this library does not have; the sub-project `lean_rank/`
(Lean v4.36.0, on M. Stoll's `EllipticCurves`) proves it for the 254 curves of `TwoDescentCertificates`. Not included: Lemmas B and D, which give the
non-square hypotheses for every ratio (`UNIFORM_RANK.md`, section 4). -/
theorem theorem2_of (m n k : ℚ) (hk : m ^ 2 + n ^ 2 = 2 * k ^ 2) (hm : m ≠ 0) (hn : n ≠ 0)
    (hp : m + n ≠ 0) (hd : m - n ≠ 0)
    (hH : ¬ IsSquare (genH m n)) (hH' : ¬ IsSquare (genH m (-n)))
    (hHH : ¬ IsSquare (genH m n * genH m (-n)))
    (hmn : ¬ IsSquare (2 * m * n)) (hmn' : ¬ IsSquare (-2 * m * n)) :
    Theorem2At m n k := by
  unfold Theorem2At
  have hS : m ^ 2 + n ^ 2 ≠ 0 := by
    have h1 : 0 < m ^ 2 := lt_of_le_of_ne (sq_nonneg m) (Ne.symm (pow_ne_zero 2 hm))
    have h2 := sq_nonneg n
    exact (by linarith : 0 < m ^ 2 + n ^ 2).ne'
  have hk0 : k ≠ 0 := by rintro rfl; apply hS; linear_combination hk
  have hH0 : genH m n ≠ 0 := fun h => hH ⟨0, by rw [h]; ring⟩
  have hH0' : genH m (-n) ≠ 0 := fun h => hH' ⟨0, by rw [h]; ring⟩
  have hq : qe (eP m n) (aP m n) (bP m n) ≠ 0 := by
    rw [qe_P]
    exact mul_ne_zero (mul_ne_zero (by norm_num) (pow_ne_zero 2 hS)) (mul_ne_zero hH0 hH0')
  have hyR : 2 * yR m n ≠ 0 := by
    unfold yR
    exact mul_ne_zero two_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero
      (mul_ne_zero (by norm_num) hm) hn) hp) hS) hH0)
  have hyQ : 2 * yQ m n k ≠ 0 := by
    unfold yQ
    exact mul_ne_zero two_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero
      (mul_ne_zero (by norm_num) hk0) hm) hn) (pow_ne_zero 2 hp)) (pow_ne_zero 2 hS))
  have hR : (Eplus m n).toAffine.Nonsingular (xR m n) (yR m n) := nonsingular_of_y (R_on m n) hyR
  have hQ : (Eplus m n).toAffine.Nonsingular 0 (yQ m n k) := nonsingular_of_y (Q_on hk) hyQ
  have hT : (Eplus m n).toAffine.Nonsingular (eP m n) 0 := nonsingular_T hq
  refine ⟨hR, hQ, hT, ?_⟩
  have h2S : IsSquare (2 * (m ^ 2 + n ^ 2)) := ⟨2 * k, by linear_combination 2 * hk⟩
  have h2S0 : 2 * (m ^ 2 + n ^ 2) ≠ 0 := mul_ne_zero two_ne_zero hS
  have hxR : xR m n ≠ eP m n := by
    intro h; apply mul_ne_zero h2S0 hH0; rw [← alphaR_val, h, sub_self]
  have hxQ : (0 : ℚ) ≠ eP m n := by
    unfold eP; intro h; apply hS
    have h3 : (m ^ 2 + n ^ 2) ^ 3 = 0 := by linear_combination h / 2
    exact (pow_eq_zero_iff (by norm_num)).mp h3
  have aR : alpha (eP m n) (aP m n) (bP m n) (Point.some _ _ hR) = 2 * (m ^ 2 + n ^ 2) * genH m n := by
    rw [alpha_some, if_neg hxR, alphaR_val]
  have aQ : alpha (eP m n) (aP m n) (bP m n) (Point.some _ _ hQ) = 2 * (m ^ 2 + n ^ 2) ^ 3 := by
    rw [alpha_some, if_neg hxQ]; unfold eP; ring
  have aT : alpha (eP m n) (aP m n) (bP m n) (Point.some _ _ hT) = qe (eP m n) (aP m n) (bP m n) := by
    rw [alpha_some, if_pos rfl]
  intro c₁ c₂ c₃ W hW
  refine indep_of_alpha_int hq _ _ _ ?_ ?_ ?_ ?_ ?_ c₁ c₂ c₃ W hW
  · rw [aR]; intro h; exact hH (isSquare_cancel_left h2S h2S0 h)
  · rw [aT, qe_P]; intro h
    exact hHH (isSquare_cancel_left ⟨2 * (m ^ 2 + n ^ 2), by ring⟩
      (mul_ne_zero (by norm_num) (pow_ne_zero 2 hS)) h)
  · rw [aR, aT, qe_P]; intro h
    apply hH'
    refine isSquare_cancel_left (x := 2 * (m ^ 2 + n ^ 2) * (2 * (m ^ 2 + n ^ 2) * genH m n) ^ 2)
      ⟨2 * k * (2 * (m ^ 2 + n ^ 2) * genH m n), ?_⟩
      (mul_ne_zero h2S0 (pow_ne_zero 2 (mul_ne_zero h2S0 hH0))) ?_
    · linear_combination 2 * (2 * (m ^ 2 + n ^ 2) * genH m n) ^ 2 * hk
    · convert h using 1; ring
  · rw [aQ]; exact ⟨2 * k * (m ^ 2 + n ^ 2), by linear_combination 2 * (m ^ 2 + n ^ 2) ^ 2 * hk⟩
  · intro W' hW'
    obtain ⟨u, v, h1, h2⟩ := gauss_double (aP_eq m n) (bP_eq m n) hW'.symm
    rcases gauss_family hp hd h1 h2 with h | h
    · exact hmn h
    · exact hmn' h

end Family

/-! ## Certificates: the ratio `(7, 1)`, `m² + n² = 50 = 2·5²`

`H(7, 1) = 1828` (`42² < 1828 < 43²`), `H(7, −1) = 3172` (`56² < 3172 < 57²`), their product `5798416`
(`2407² < 5798416 < 2408²`), `2·7·1 = 14` (`3² < 14 < 4²`) and `−14 < 0`. -/

section Certificates

open MagicSquareOfSquares.G1

/-- `UNIFORM_RANK.md` Theorem 2 on `E⁺` at `(7, 1)`. -/
theorem theorem2_7_1 :
    ∃ (hR : (Eplus 7 1).toAffine.Nonsingular (xR 7 1) (yR 7 1))
      (hQ : (Eplus 7 1).toAffine.Nonsingular 0 (yQ 7 1 5))
      (hT : (Eplus 7 1).toAffine.Nonsingular (eP 7 1) 0),
      ∀ c₁ c₂ c₃ : ℤ, ∀ W : (Eplus 7 1).toAffine.Point,
        c₁ • Point.some _ _ hR + c₂ • Point.some _ _ hQ + c₃ • Point.some _ _ hT = W + W →
          Even c₁ ∧ Even c₂ ∧ Even c₃ :=
  theorem2_of 7 1 5 (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by rw [show genH 7 1 = ((1828 : ℕ) : ℚ) by unfold genH; norm_num]
        exact not_isSquare_of_bracket 1828 42 (by norm_num) (by norm_num))
    (by rw [show genH 7 (-1) = ((3172 : ℕ) : ℚ) by unfold genH; norm_num]
        exact not_isSquare_of_bracket 3172 56 (by norm_num) (by norm_num))
    (by rw [show genH 7 1 * genH 7 (-1) = ((5798416 : ℕ) : ℚ) by unfold genH; norm_num]
        exact not_isSquare_of_bracket 5798416 2407 (by norm_num) (by norm_num))
    (by rw [show (2 : ℚ) * 7 * 1 = ((14 : ℕ) : ℚ) by norm_num]
        exact not_isSquare_of_bracket 14 3 (by norm_num) (by norm_num))
    (not_isSquare_of_neg (by norm_num))

/-- `UNIFORM_RANK.md` Theorem 2 on `E⁻` at `(7, 1)`, which is `Eplus 7 (−1)`. -/
theorem theorem2_7_m1 :
    ∃ (hR : (Eplus 7 (-1)).toAffine.Nonsingular (xR 7 (-1)) (yR 7 (-1)))
      (hQ : (Eplus 7 (-1)).toAffine.Nonsingular 0 (yQ 7 (-1) 5))
      (hT : (Eplus 7 (-1)).toAffine.Nonsingular (eP 7 (-1)) 0),
      ∀ c₁ c₂ c₃ : ℤ, ∀ W : (Eplus 7 (-1)).toAffine.Point,
        c₁ • Point.some _ _ hR + c₂ • Point.some _ _ hQ + c₃ • Point.some _ _ hT = W + W →
          Even c₁ ∧ Even c₂ ∧ Even c₃ :=
  theorem2_of 7 (-1) 5 (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by rw [show genH 7 (-1) = ((3172 : ℕ) : ℚ) by unfold genH; norm_num]
        exact not_isSquare_of_bracket 3172 56 (by norm_num) (by norm_num))
    (by rw [show genH 7 (-(-1)) = ((1828 : ℕ) : ℚ) by unfold genH; norm_num]
        exact not_isSquare_of_bracket 1828 42 (by norm_num) (by norm_num))
    (by rw [show genH 7 (-1) * genH 7 (-(-1)) = ((5798416 : ℕ) : ℚ) by unfold genH; norm_num]
        exact not_isSquare_of_bracket 5798416 2407 (by norm_num) (by norm_num))
    (not_isSquare_of_neg (by norm_num))
    (by rw [show (-2 : ℚ) * 7 * (-1) = ((14 : ℕ) : ℚ) by norm_num]
        exact not_isSquare_of_bracket 14 3 (by norm_num) (by norm_num))

/-- **Negative control.** On `E⁺` at `(7, 1)`, Mathlib's `R + R` is `P = (−160000, 48000000)`, the point of
`z = 2`; `α(P) = x(P) − e = 90000 = 300²` is a square, and the halving test at `θ` passes as well:
`x(P) − r = u² − v²`, `−s = 2uv` with `u = −112`, `v = 384` (`u + iv = i(m + n)(m + in)²`). So the tests do
not refuse a point that is a double. -/
theorem control_7_1 (hR : (Eplus 7 1).toAffine.Nonsingular (xR 7 1) (yR 7 1)) :
    ∃ hP : (Eplus 7 1).toAffine.Nonsingular (-160000) 48000000,
      Point.some _ _ hR + Point.some _ _ hR = Point.some _ _ hP ∧
      alpha (eP 7 1) (aP 7 1) (bP 7 1) (Point.some _ _ hP) = 300 * 300 ∧
      (-160000 : ℚ) - rP 7 1 = (-112) ^ 2 - 384 ^ 2 ∧ -sP 7 1 = 2 * (-112) * 384 := by
  have hy : yR 7 1 ≠ (Eplus 7 1).toAffine.negY (xR 7 1) (yR 7 1) := by
    rw [curve_negY]; unfold yR genH; norm_num
  have hP : (Eplus 7 1).toAffine.Nonsingular (-160000) 48000000 :=
    nonsingular_of_y (by unfold eP aP bP; norm_num) (by norm_num)
  refine ⟨hP, ?_, ?_, ?_, ?_⟩
  · rw [Point.add_self_of_Y_ne hy, Point.some.injEq]
    have hℓ : (Eplus 7 1).toAffine.slope (xR 7 1) (xR 7 1) (yR 7 1) (yR 7 1) = 76 := by
      rw [Affine.slope_of_Y_ne rfl hy, curve_negY]
      simp only [curve_a₁, curve_a₂, curve_a₄]
      unfold xR yR genH eP aP bP; norm_num
    rw [hℓ]
    simp only [Affine.addY, Affine.negAddY, Affine.addX, curve_negY, curve_a₁, curve_a₂]
    unfold xR yR genH eP aP; norm_num
  · rw [alpha_some, if_neg (by unfold eP; norm_num)]; unfold eP; norm_num
  · unfold rP; norm_num
  · unfold sP; norm_num

/-- **Negative control for the bracket**: `90000 = 300²` has no square-root bracket. -/
theorem control_bracket : ¬ ∃ s : ℕ, s * s < 90000 ∧ 90000 < (s + 1) * (s + 1) := by
  rintro ⟨s, h1, h2⟩
  rcases le_or_gt s 299 with h | h
  · have := Nat.mul_self_le_mul_self (show s + 1 ≤ 300 by omega)
    omega
  · have := Nat.mul_self_le_mul_self (show 300 ≤ s from h)
    omega

end Certificates

end MagicSquareOfSquares.TwoDescent
