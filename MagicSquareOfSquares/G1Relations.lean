import Mathlib.Tactic

/-!
# The two relations used by the Demjanenko–Manin computation

`G1_FULL_SYSTEM.md`, Result 4, uses two elliptic curves, each with two rational points coming
from one solution `t` of the five-form system. The rank and height computations of that result
are outside Lean. This file checks only the algebra: the points lie on the stated curves, on the
same curve and not on a twist.

Integer forms: `F0 = m² + n²t²`, `F8 = n² + m²t²`, `G4 = 2(m² + n²)(t² + 1)`,
`GL = G4 − 8mnt`, `GU = G4 + 8mnt`.

* Curve `E04: y² = x(x + 2Sm²)(x + 2Sn²)`, `S = m² + n²`, `κ = 2Sn²`. At `x = κt²` the right side
  is a square times `F0·G4` (`e04_point`); at `x = κ/t²` it is a square times `F8·G4`
  (`e04_point'`). The two abscissas have product `κ²`.
* Curve `D1: y² = d1 m n x`. With `c0 = −16Smn(m² − n²)²` and `z = t + 1/t`, at `x = c0/z` the
  right side is a square times `F0·F8·G4·GL` (`d1_point`), and at `x = −c0/z` it is a square times
  `F0·F8·G4·GU` (`d1_point'`). The two abscissas are opposite.
-/

namespace MagicSquareOfSquares.G1

theorem e04_point (m n t : ℚ) :
    (2 * (m ^ 2 + n ^ 2) * n ^ 2 * t ^ 2) *
        (2 * (m ^ 2 + n ^ 2) * n ^ 2 * t ^ 2 + 2 * (m ^ 2 + n ^ 2) * m ^ 2) *
        (2 * (m ^ 2 + n ^ 2) * n ^ 2 * t ^ 2 + 2 * (m ^ 2 + n ^ 2) * n ^ 2) =
      (2 * (m ^ 2 + n ^ 2) * n ^ 2 * t) ^ 2 *
        ((m ^ 2 + n ^ 2 * t ^ 2) * (2 * (m ^ 2 + n ^ 2) * (t ^ 2 + 1))) := by
  ring

theorem e04_point' (m n t : ℚ) (ht : t ≠ 0) :
    (2 * (m ^ 2 + n ^ 2) * n ^ 2 / t ^ 2) *
        (2 * (m ^ 2 + n ^ 2) * n ^ 2 / t ^ 2 + 2 * (m ^ 2 + n ^ 2) * m ^ 2) *
        (2 * (m ^ 2 + n ^ 2) * n ^ 2 / t ^ 2 + 2 * (m ^ 2 + n ^ 2) * n ^ 2) =
      (2 * (m ^ 2 + n ^ 2) * n ^ 2 / t ^ 3) ^ 2 *
        ((n ^ 2 + m ^ 2 * t ^ 2) * (2 * (m ^ 2 + n ^ 2) * (t ^ 2 + 1))) := by
  field_simp
  ring

/-- The cubic of the curve `D1`: `x³ + c1·x² + c0·c2·x + c0²·c3`. -/
def d1 (m n x : ℚ) : ℚ :=
  x ^ 3 + 4 * (m ^ 2 + n ^ 2) ^ 2 * (m ^ 2 - n ^ 2) ^ 2 * x ^ 2 +
    (-16 * (m ^ 2 + n ^ 2) * m * n * (m ^ 2 - n ^ 2) ^ 2) *
      (-16 * (m ^ 2 + n ^ 2) * m ^ 3 * n ^ 3) * x +
    (-16 * (m ^ 2 + n ^ 2) * m * n * (m ^ 2 - n ^ 2) ^ 2) ^ 2 *
      (4 * (m ^ 2 + n ^ 2) ^ 2 * m ^ 2 * n ^ 2)

theorem d1_point (m n t : ℚ) :
    d1 m n (-16 * (m ^ 2 + n ^ 2) * m * n * (m ^ 2 - n ^ 2) ^ 2 * t / (t ^ 2 + 1)) =
      (-16 * (m ^ 2 + n ^ 2) * m * n * (m ^ 2 - n ^ 2) ^ 2 / (t ^ 2 + 1) ^ 2) ^ 2 *
        ((m ^ 2 + n ^ 2 * t ^ 2) * (n ^ 2 + m ^ 2 * t ^ 2) * (2 * (m ^ 2 + n ^ 2) * (t ^ 2 + 1)) *
          (2 * (m ^ 2 + n ^ 2) * (t ^ 2 + 1) - 8 * m * n * t)) := by
  have h1 : t ^ 2 + 1 ≠ 0 := by positivity
  unfold d1
  field_simp
  ring

theorem d1_point' (m n t : ℚ) :
    d1 m n (-(-16 * (m ^ 2 + n ^ 2) * m * n * (m ^ 2 - n ^ 2) ^ 2 * t / (t ^ 2 + 1))) =
      (-16 * (m ^ 2 + n ^ 2) * m * n * (m ^ 2 - n ^ 2) ^ 2 / (t ^ 2 + 1) ^ 2) ^ 2 *
        ((m ^ 2 + n ^ 2 * t ^ 2) * (n ^ 2 + m ^ 2 * t ^ 2) * (2 * (m ^ 2 + n ^ 2) * (t ^ 2 + 1)) *
          (2 * (m ^ 2 + n ^ 2) * (t ^ 2 + 1) + 8 * m * n * t)) := by
  have h1 : t ^ 2 + 1 ≠ 0 := by positivity
  unfold d1
  field_simp
  ring

end MagicSquareOfSquares.G1
