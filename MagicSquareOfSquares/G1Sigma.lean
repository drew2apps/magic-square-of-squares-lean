import Mathlib.Tactic

/-!
# The curves of `F0·F8·GL`: the algebra behind `scripts/g1_full/dm_sigma.py`

`G1_FULL_SYSTEM.md`, Result 4, "The curves of `F0·F8·GL`". For the four corner forms
`F0 = m² + n²t²`, `F8 = n² + m²t²`, `GL = 2(m² + n²)(t² + 1) − 8mnt`, `GU = 2(m² + n²)(t² + 1) + 8mnt`
of the family G1, and `z = t + 1/t`:

* `F0·F8·GL = t³·g(z)` with `g(z) = (2Sz − 8mn)·Q(z)`, `S = m² + n²`, `Q(z) = m²n²z² + (m² − n²)²`
  (`sigma_palin`); `F0·F8·GU = −t³·g(−z)` (`sigma_palin_U`).
* So `(z + 2)·g(z)` and `(z − 2)·g(z)` at `z = t + 1/t` are squares times `F0·F8·GL`, and at
  `z = −(t + 1/t)` squares times `F0·F8·GU` (`sigma_plus`, `sigma_minus`, `sigma_plus_U`,
  `sigma_minus_U`). The points of the `U` side lie on the same two curves, not on twists: the
  factors are `((t ± 1)/t²)²`.
* A quartic `Y² = c₁h + c₂h² + c₃h³ + c₄h⁴` in `h = z − z₀` goes to the Weierstrass model
  `y² = X³ + c₂X² + c₁c₃X + c₁²c₄` by `X = c₁/h`, `y = c₁Y/h²` (`quartic_root_weierstrass`); and
  `z ↦ −z` is `X ↦ −c₁X/(2z₀X + c₁)` (`sigma_involution`).
* The heights: for that Möbius map, `|c₁p|` and `|2z₀p + c₁q|` are at most `(|c₁| + 2|z₀|)·max(|p|, |q|)`
  (`moebius_height`), which gives `H(X') ≤ (|A| + 4)·H(X)` for `z₀ = ±2`.
* Translation by a point of order 2: on `y² = (x − e)(x² + bx + c)`, the abscissa of `P + (e, 0)` is
  `e + (e² + be + c)/(x − e)` (`two_torsion_translate`).

Ranks, generators and heights are not in Lean (`sigma_exact.py`, `sigma_certify.py`).
-/

namespace MagicSquareOfSquares.G1

/-- The cubic `g(z) = (2(m² + n²)z − 8mn)(m²n²z² + (m² − n²)²)`. -/
def sigmaG (m n z : ℚ) : ℚ :=
  (2 * (m ^ 2 + n ^ 2) * z - 8 * m * n) * (m ^ 2 * n ^ 2 * z ^ 2 + (m ^ 2 - n ^ 2) ^ 2)

theorem sigma_palin (m n t : ℚ) (ht : t ≠ 0) :
    (m ^ 2 + n ^ 2 * t ^ 2) * (n ^ 2 + m ^ 2 * t ^ 2) *
        (2 * (m ^ 2 + n ^ 2) * (t ^ 2 + 1) - 8 * m * n * t) =
      t ^ 3 * sigmaG m n (t + 1 / t) := by
  unfold sigmaG
  field_simp
  ring

theorem sigma_palin_U (m n t : ℚ) (ht : t ≠ 0) :
    (m ^ 2 + n ^ 2 * t ^ 2) * (n ^ 2 + m ^ 2 * t ^ 2) *
        (2 * (m ^ 2 + n ^ 2) * (t ^ 2 + 1) + 8 * m * n * t) =
      -(t ^ 3 * sigmaG m n (-(t + 1 / t))) := by
  unfold sigmaG
  field_simp
  ring

theorem sigma_plus (m n t : ℚ) (ht : t ≠ 0) :
    ((t + 1 / t) + 2) * sigmaG m n (t + 1 / t) =
      ((t + 1) / t ^ 2) ^ 2 *
        ((m ^ 2 + n ^ 2 * t ^ 2) * (n ^ 2 + m ^ 2 * t ^ 2) *
          (2 * (m ^ 2 + n ^ 2) * (t ^ 2 + 1) - 8 * m * n * t)) := by
  unfold sigmaG
  field_simp
  ring

theorem sigma_minus (m n t : ℚ) (ht : t ≠ 0) :
    ((t + 1 / t) - 2) * sigmaG m n (t + 1 / t) =
      ((t - 1) / t ^ 2) ^ 2 *
        ((m ^ 2 + n ^ 2 * t ^ 2) * (n ^ 2 + m ^ 2 * t ^ 2) *
          (2 * (m ^ 2 + n ^ 2) * (t ^ 2 + 1) - 8 * m * n * t)) := by
  unfold sigmaG
  field_simp
  ring

theorem sigma_plus_U (m n t : ℚ) (ht : t ≠ 0) :
    (-(t + 1 / t) + 2) * sigmaG m n (-(t + 1 / t)) =
      ((t - 1) / t ^ 2) ^ 2 *
        ((m ^ 2 + n ^ 2 * t ^ 2) * (n ^ 2 + m ^ 2 * t ^ 2) *
          (2 * (m ^ 2 + n ^ 2) * (t ^ 2 + 1) + 8 * m * n * t)) := by
  unfold sigmaG
  field_simp
  ring

theorem sigma_minus_U (m n t : ℚ) (ht : t ≠ 0) :
    (-(t + 1 / t) - 2) * sigmaG m n (-(t + 1 / t)) =
      ((t + 1) / t ^ 2) ^ 2 *
        ((m ^ 2 + n ^ 2 * t ^ 2) * (n ^ 2 + m ^ 2 * t ^ 2) *
          (2 * (m ^ 2 + n ^ 2) * (t ^ 2 + 1) + 8 * m * n * t)) := by
  unfold sigmaG
  field_simp
  ring

/-- A quartic with the root `z₀`, written in `h = z − z₀`, and its Weierstrass model. -/
theorem quartic_root_weierstrass (c1 c2 c3 c4 h Y : ℚ) (hh : h ≠ 0)
    (hY : Y ^ 2 = c1 * h + c2 * h ^ 2 + c3 * h ^ 3 + c4 * h ^ 4) :
    (c1 * Y / h ^ 2) ^ 2 =
      (c1 / h) ^ 3 + c2 * (c1 / h) ^ 2 + c1 * c3 * (c1 / h) + c1 ^ 2 * c4 := by
  rw [div_pow, mul_pow, hY]
  field_simp

/-- `z ↦ −z` in the coordinate `X = c₁/(z − z₀)`. -/
theorem sigma_involution (c1 z z0 : ℚ) (hc : c1 ≠ 0) (h1 : z - z0 ≠ 0) (h2 : -z - z0 ≠ 0) :
    c1 / (-z - z0) = -c1 * (c1 / (z - z0)) / (2 * z0 * (c1 / (z - z0)) + c1) := by
  have h3 : 2 * z0 * (c1 / (z - z0)) + c1 = c1 * (z + z0) / (z - z0) := by
    field_simp
    ring
  have h4 : z + z0 ≠ 0 := by
    intro h
    apply h2
    linarith
  rw [h3]
  field_simp
  ring

/-- The height bound for the Möbius map `X ↦ −c₁X/(2z₀X + c₁)` at `X = p/q`. -/
theorem moebius_height (c1 z0 p q : ℤ) :
    |c1 * p| ≤ (|c1| + 2 * |z0|) * max |p| |q| ∧
      |2 * z0 * p + c1 * q| ≤ (|c1| + 2 * |z0|) * max |p| |q| := by
  have hp : |p| ≤ max |p| |q| := le_max_left _ _
  have hq : |q| ≤ max |p| |q| := le_max_right _ _
  have hz : 0 ≤ |z0| := abs_nonneg _
  have hc : 0 ≤ |c1| := abs_nonneg _
  have hm : 0 ≤ max |p| |q| := le_trans (abs_nonneg p) hp
  constructor
  · rw [abs_mul]
    nlinarith [mul_le_mul_of_nonneg_left hp hc]
  · calc |2 * z0 * p + c1 * q| ≤ |2 * z0 * p| + |c1 * q| := abs_add_le _ _
      _ = 2 * |z0| * |p| + |c1| * |q| := by
          rw [abs_mul, abs_mul, abs_mul]
          norm_num
      _ ≤ 2 * |z0| * max |p| |q| + |c1| * max |p| |q| := by
          gcongr
      _ = (|c1| + 2 * |z0|) * max |p| |q| := by ring

/-- Translation by the point `(e, 0)` of order 2 on `y² = (x − e)(x² + bx + c)`. -/
theorem two_torsion_translate (b c e x y : ℚ) (hx : x - e ≠ 0)
    (hy : y ^ 2 = (x - e) * (x ^ 2 + b * x + c)) :
    (y / (x - e)) ^ 2 - (b - e) - x - e = e + (e ^ 2 + b * e + c) / (x - e) := by
  rw [div_pow, hy]
  field_simp
  ring

end MagicSquareOfSquares.G1
