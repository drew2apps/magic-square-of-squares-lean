import LeanRank.LemmaDReduction

/-!
# The 2-isogeny pair of `y² = x(x² + ax + b)` as explicit maps

Written by an AI model (Claude, made by Anthropic) on 25 November 2026 (UTC); not read by a person.
Experiment, not counted in the census. A first piece of step 1d of the plan for Lemma D
(a plan kept outside this repository); field identities only, over any field of characteristic 0.

`E : y² = x(x² + ax + b)`, `E′ : Y² = X(X² − 2aX + a² − 4b)` (for `E1`: `a = 8d`, `b = −8d²`).
* `isog_on`: `φ(x, y) = (y²/x², y(b − x²)/x²)` lies on `E′` (`x ≠ 0`).
* `isog_dual_on`: `φ̂(X, Y) = (Y²/(4X²), Y(a² − 4b − X²)/(8X²))` lies on `E` (`X ≠ 0`).
* `isog_comp_x`: `x(φ̂(φ(P))) = (x² − b)²/(4y²)`, and `double_x`: this is `λ² − a − 2x` with the tangent slope
  `λ = (3x² + 2ax + b)/(2y)`, i.e. the abscissa of `2P` (Mathlib's `addX x x (slope x x y y)` for `a₁ = a₃ = 0`;
  that last identification is not made here).
-/

namespace MagicSquareOfSquares.Family

variable {K : Type*} [Field K] [CharZero K]

omit [CharZero K] in
theorem isog_on {a b x y : K} (hx : x ≠ 0) (h : y ^ 2 = x * (x ^ 2 + a * x + b)) :
    (y * (b - x ^ 2) / x ^ 2) ^ 2 =
      y ^ 2 / x ^ 2 * ((y ^ 2 / x ^ 2) ^ 2 - 2 * a * (y ^ 2 / x ^ 2) + (a ^ 2 - 4 * b)) := by
  rw [← sub_eq_zero]
  have e : (y * (b - x ^ 2) / x ^ 2) ^ 2 -
      y ^ 2 / x ^ 2 * ((y ^ 2 / x ^ 2) ^ 2 - 2 * a * (y ^ 2 / x ^ 2) + (a ^ 2 - 4 * b))
      = -y ^ 2 * (-a * x ^ 2 + b * x + x ^ 3 + y ^ 2) * (y ^ 2 - x * (x ^ 2 + a * x + b)) / x ^ 6 := by
    field_simp
    ring
  rw [e, h, sub_self, mul_zero, zero_div]

theorem isog_dual_on {a b X Y : K} (hX : X ≠ 0)
    (h : Y ^ 2 = X * (X ^ 2 - 2 * a * X + (a ^ 2 - 4 * b))) :
    (Y * (a ^ 2 - 4 * b - X ^ 2) / (8 * X ^ 2)) ^ 2 =
      Y ^ 2 / (4 * X ^ 2) * ((Y ^ 2 / (4 * X ^ 2)) ^ 2 + a * (Y ^ 2 / (4 * X ^ 2)) + b) := by
  rw [← sub_eq_zero]
  have e : (Y * (a ^ 2 - 4 * b - X ^ 2) / (8 * X ^ 2)) ^ 2 -
      Y ^ 2 / (4 * X ^ 2) * ((Y ^ 2 / (4 * X ^ 2)) ^ 2 + a * (Y ^ 2 / (4 * X ^ 2)) + b)
      = -Y ^ 2 * (X ^ 3 + 2 * X ^ 2 * a + X * a ^ 2 - 4 * X * b + Y ^ 2)
        * (Y ^ 2 - X * (X ^ 2 - 2 * a * X + (a ^ 2 - 4 * b))) / (64 * X ^ 6) := by
    field_simp
    ring
  rw [e, h, sub_self, mul_zero, zero_div]

theorem isog_comp_x {b x y : K} (hx : x ≠ 0) (hy : y ≠ 0) :
    (y * (b - x ^ 2) / x ^ 2) ^ 2 / (4 * (y ^ 2 / x ^ 2) ^ 2) = (x ^ 2 - b) ^ 2 / (4 * y ^ 2) := by
  field_simp
  ring

theorem double_x {a b x y : K} (hy : y ≠ 0) (h : y ^ 2 = x * (x ^ 2 + a * x + b)) :
    ((3 * x ^ 2 + 2 * a * x + b) / (2 * y)) ^ 2 - a - 2 * x = (x ^ 2 - b) ^ 2 / (4 * y ^ 2) := by
  rw [← sub_eq_zero]
  have e : ((3 * x ^ 2 + 2 * a * x + b) / (2 * y)) ^ 2 - a - 2 * x - (x ^ 2 - b) ^ 2 / (4 * y ^ 2)
      = -(a + 2 * x) * (y ^ 2 - x * (x ^ 2 + a * x + b)) / y ^ 2 := by
    field_simp
    ring
  rw [e, h, sub_self, mul_zero, zero_div]

end MagicSquareOfSquares.Family

#print axioms MagicSquareOfSquares.Family.isog_on
#print axioms MagicSquareOfSquares.Family.isog_dual_on
#print axioms MagicSquareOfSquares.Family.isog_comp_x
#print axioms MagicSquareOfSquares.Family.double_x
