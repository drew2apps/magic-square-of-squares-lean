import Mathlib

/-!
# Squares in arithmetic progression and doubles on a congruent-number curve

For a rational `n ≠ 0`, let `E_n` be the curve `y² = x³ - n²x`. A rational number `x` has
`x - n`, `x`, `x + n` all squares of rationals if and only if `x` is the `x`-coordinate of a
rational point of `E_n` that is the double of a rational point of `E_n`.

## Main results

* `congruentCurve`: the curve `E_n` over `ℚ`.
* `IsDoubleX`: `x` is the `x`-coordinate of the double of a rational point of `E_n`.
* `isSquare_iff_isDoubleX`: the equivalence above.
-/

namespace MagicSquareOfSquares

open WeierstrassCurve.Affine

/-- The curve `y² = x³ - n²x` over `ℚ`. -/
def congruentCurve (n : ℚ) : WeierstrassCurve.Affine ℚ := ⟨0, 0, 0, -n ^ 2, 0⟩

/-- `x` is the `x`-coordinate of an affine rational point of `y² = x³ - n²x` that is the
double of a rational point of that curve. -/
def IsDoubleX (n x : ℚ) : Prop :=
  ∃ (y : ℚ) (h : (congruentCurve n).Nonsingular x y) (Q : (congruentCurve n).Point),
    Point.some x y h = Q + Q

private theorem sq_helper {v A B c : ℚ} (hv : v ≠ 0) (key : A ^ 2 - c * v ^ 2 = B ^ 2) :
    (A / (2 * v)) ^ 2 - c / 4 = (B / (2 * v)) ^ 2 := by
  have h2 : (2 * v) ^ 2 ≠ 0 := pow_ne_zero 2 (mul_ne_zero two_ne_zero hv)
  rw [div_pow, div_pow, ← key, eq_div_iff h2, sub_mul, div_mul_cancel₀ _ h2]
  ring

private theorem curve_eq {x n α β γ : ℚ} (hα : x - n = α * α) (hγ : x = γ * γ)
    (hβ : x + n = β * β) :
    ((α + γ) * (γ + β) * (α + β)) ^ 2 =
      ((α + γ) * (γ + β)) ^ 3 - n ^ 2 * ((α + γ) * (γ + β)) := by
  linear_combination ((n + γ ^ 2 - α ^ 2) * ((α + γ) * (γ + β))) * (hγ - hα) -
    ((α + γ) * ((α + γ) * (γ + β) * (α + β))) * (hα + hβ - 2 * hγ)

private theorem half_eq {x n α β γ : ℚ} (hα : x - n = α * α) (hγ : x = γ * γ)
    (hβ : x + n = β * β) :
    ((α + γ) * (γ + β)) ^ 2 + n ^ 2 = 2 * ((α + γ) * (γ + β) * (α + β)) * γ := by
  linear_combination (n + γ ^ 2 - α ^ 2) * (hγ - hα) - (α ^ 2 - γ ^ 2) * (hα + hβ - 2 * hγ)

private theorem slope_eq {n u v : ℚ} (hY : v ≠ (congruentCurve n).negY u v) :
    (congruentCurve n).slope u u v v = (3 * u ^ 2 - n ^ 2) / (2 * v) := by
  rw [slope_of_Y_ne rfl hY]
  simp only [negY, congruentCurve]
  congr 1 <;> ring

private theorem addX_eq {n u v : ℚ} (hY : v ≠ (congruentCurve n).negY u v) :
    (congruentCurve n).addX u u ((congruentCurve n).slope u u v v) =
      ((3 * u ^ 2 - n ^ 2) / (2 * v)) ^ 2 - 2 * u := by
  rw [slope_eq hY]
  simp only [addX, congruentCurve]
  ring

/-- For `n ≠ 0`: the rationals `x - n`, `x`, `x + n` are all squares if and only if `x` is
the `x`-coordinate of the double of a rational point of `y² = x³ - n²x`. -/
theorem isSquare_iff_isDoubleX {n x : ℚ} (hn : n ≠ 0) :
    (IsSquare (x - n) ∧ IsSquare x ∧ IsSquare (x + n)) ↔ IsDoubleX n x := by
  constructor
  · rintro ⟨⟨α, hα⟩, ⟨γ, hγ⟩, ⟨β, hβ⟩⟩
    have h1 : α + γ ≠ 0 := by
      intro h0
      have hn' : n = (γ - α) * (α + γ) := by linear_combination hγ - hα
      rw [h0, mul_zero] at hn'
      exact hn hn'
    have h2 : γ + β ≠ 0 := by
      intro h0
      have hn' : n = (β - γ) * (γ + β) := by linear_combination hβ - hγ
      rw [h0, mul_zero] at hn'
      exact hn hn'
    have h3 : α + β ≠ 0 := by
      intro h0
      have hn' : 2 * n = (β - α) * (α + β) := by linear_combination hβ - hα
      rw [h0, mul_zero] at hn'
      exact hn (by linarith)
    obtain ⟨u, hu⟩ : ∃ u, u = (α + γ) * (γ + β) := ⟨_, rfl⟩
    obtain ⟨v, hv⟩ : ∃ v, v = (α + γ) * (γ + β) * (α + β) := ⟨_, rfl⟩
    have hv0 : v ≠ 0 := by
      rw [hv]
      exact mul_ne_zero (mul_ne_zero h1 h2) h3
    have hE : v ^ 2 = u ^ 3 - n ^ 2 * u := by
      rw [hu, hv]
      exact curve_eq hα hγ hβ
    have hD : u ^ 2 + n ^ 2 = 2 * v * γ := by
      rw [hu, hv]
      exact half_eq hα hγ hβ
    have hEq : (congruentCurve n).Equation u v := by
      rw [equation_iff]
      simp only [congruentCurve]
      linear_combination hE
    have hNs : (congruentCurve n).Nonsingular u v := by
      refine (nonsingular_iff' u v).mpr ⟨hEq, Or.inr ?_⟩
      simp only [congruentCurve]
      intro h0
      apply hv0
      linarith
    have hY : v ≠ (congruentCurve n).negY u v := by
      simp only [negY, congruentCurve]
      intro h0
      apply hv0
      linarith
    have hXval : ((3 * u ^ 2 - n ^ 2) / (2 * v)) ^ 2 - 2 * u = γ * γ := by
      have k0 : (3 * u ^ 2 - n ^ 2) ^ 2 - (8 * u) * v ^ 2 = (u ^ 2 + n ^ 2) ^ 2 := by
        rw [hE]
        ring
      have hs := sq_helper hv0 k0
      have hγ' : (u ^ 2 + n ^ 2) / (2 * v) = γ := by
        rw [hD]
        field_simp
      rw [hγ'] at hs
      linear_combination hs
    have hX : (congruentCurve n).addX u u ((congruentCurve n).slope u u v v) = x := by
      rw [addX_eq hY, hγ]
      exact hXval
    subst hX
    exact ⟨_, _, Point.some u v hNs, (Point.add_self_of_Y_ne hY).symm⟩
  · rintro ⟨y, h, Q, hQ⟩
    rcases Q with _ | ⟨u, v, hNs⟩
    · rw [← Point.zero_def, add_zero] at hQ
      exact absurd hQ (Point.some_ne_zero h)
    · by_cases hY : v = (congruentCurve n).negY u v
      · rw [Point.add_self_of_Y_eq hY] at hQ
        exact absurd hQ (Point.some_ne_zero h)
      · rw [Point.add_self_of_Y_ne hY] at hQ
        have hx : x = ((3 * u ^ 2 - n ^ 2) / (2 * v)) ^ 2 - 2 * u := by
          rw [← addX_eq hY]
          exact (Point.some.inj hQ).1
        have hv0 : v ≠ 0 := by
          intro h0
          apply hY
          simp [h0, negY, congruentCurve]
        have hE : v ^ 2 = u ^ 3 - n ^ 2 * u := by
          have hEq := (equation_iff u v).mp ((nonsingular_iff' u v).mp hNs).1
          simp only [congruentCurve] at hEq
          linear_combination hEq
        have k0 : (3 * u ^ 2 - n ^ 2) ^ 2 - (8 * u) * v ^ 2 = (u ^ 2 + n ^ 2) ^ 2 := by
          rw [hE]
          ring
        have kp : (3 * u ^ 2 - n ^ 2) ^ 2 - (8 * u + 4 * n) * v ^ 2 =
            (u ^ 2 - 2 * n * u - n ^ 2) ^ 2 := by
          rw [hE]
          ring
        have km : (3 * u ^ 2 - n ^ 2) ^ 2 - (8 * u - 4 * n) * v ^ 2 =
            (u ^ 2 + 2 * n * u - n ^ 2) ^ 2 := by
          rw [hE]
          ring
        refine ⟨⟨(u ^ 2 - 2 * n * u - n ^ 2) / (2 * v), ?_⟩, ⟨(u ^ 2 + n ^ 2) / (2 * v), ?_⟩,
          ⟨(u ^ 2 + 2 * n * u - n ^ 2) / (2 * v), ?_⟩⟩
        · rw [hx]
          linear_combination sq_helper hv0 kp
        · rw [hx]
          linear_combination sq_helper hv0 k0
        · rw [hx]
          linear_combination sq_helper hv0 km

end MagicSquareOfSquares
