import Mathlib.Tactic

/-!
# Curve-membership algebra for the Demjanenko–Manin computation

`G1_FULL_SYSTEM.md`, Result 4, uses three elliptic curves (`scripts/g1_full/dm.py` for `E04`,
`dm2.py` for `D1`, `dm3.py` for `E0L`), each with two rational points coming from one solution
`t` of the five-form system. The rank and height computations of that result
are outside Lean. This file checks only the algebra of curve membership: for arbitrary rational
`m`, `n`, `t` the stated abscissas give a square times the stated product of forms, so when the
forms are squares the points are rational points of the stated curves, of the same curve and not
of a twist. It does not prove anything about ranks, generators or heights.

Integer forms: `F0 = m² + n²t²`, `F8 = n² + m²t²`, `G4 = 2(m² + n²)(t² + 1)`,
`GL = G4 − 8mnt`, `GU = G4 + 8mnt`.

* Curve `E04: y² = x(x + 2Sm²)(x + 2Sn²)`, `S = m² + n²`, `κ = 2Sn²`. At `x = κt²` the right side
  is a square times `F0·G4` (`e04_point`); at `x = κ/t²` it is a square times `F8·G4`
  (`e04_point'`). The two abscissas have product `κ²`.
* Curve `D1: y² = d1 m n x`. With `c0 = −16Smn(m² − n²)²` and `z = t + 1/t`, at `x = c0/z` the
  right side is a square times `F0·F8·G4·GL` (`d1_point`), and at `x = −c0/z` it is a square times
  `F0·F8·G4·GU` (`d1_point'`). The two abscissas are opposite.
* Curve `E0L: y² = F0(t)·GL(t)`. Its value at `−t` is `F0(t)·GU(t)` (`e0l_neg`; a symmetry of the
  factorisation, not a point), and `t = m/n` is a rational point of the quartic (`e0l_base`).
  The passage from this quartic to the Weierstrass model used by `dm3.py` is not in Lean.
* The partner pair: with `M = m + n`, `N = m − n` and `t' = (t + 1)/(t − 1)`, the five forms of
  `(M, N)` at `t'` are the five forms of `(m, n)` at `t`, in another order, times squares
  (`partner_forms`). When `m`, `n` are both odd the scripts use `(M/2, N/2)`, which divides every
  form by `4`.
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

theorem e0l_neg (m n t : ℚ) :
    (m ^ 2 + n ^ 2 * (-t) ^ 2) * (2 * (m ^ 2 + n ^ 2) * ((-t) ^ 2 + 1) - 8 * m * n * (-t)) =
      (m ^ 2 + n ^ 2 * t ^ 2) * (2 * (m ^ 2 + n ^ 2) * (t ^ 2 + 1) + 8 * m * n * t) := by
  ring

theorem e0l_base (m n : ℚ) (hn : n ≠ 0) :
    (m ^ 2 + n ^ 2 * (m / n) ^ 2) *
        (2 * (m ^ 2 + n ^ 2) * ((m / n) ^ 2 + 1) - 8 * m * n * (m / n)) =
      (2 * m * (m ^ 2 - n ^ 2) / n) ^ 2 := by
  field_simp
  ring

/-- The five forms of the partner pair `(m + n, m − n)` at `(t + 1)/(t − 1)`. -/
theorem partner_forms (m n t : ℚ) (ht : t - 1 ≠ 0) :
    (m + n) ^ 2 + (m - n) ^ 2 * ((t + 1) / (t - 1)) ^ 2 =
        (2 * (m ^ 2 + n ^ 2) * (t ^ 2 + 1) - 8 * m * n * t) / (t - 1) ^ 2 ∧
      (m - n) ^ 2 + (m + n) ^ 2 * ((t + 1) / (t - 1)) ^ 2 =
        (2 * (m ^ 2 + n ^ 2) * (t ^ 2 + 1) + 8 * m * n * t) / (t - 1) ^ 2 ∧
      2 * ((m + n) ^ 2 + (m - n) ^ 2) * (((t + 1) / (t - 1)) ^ 2 + 1) =
        2 ^ 2 * (2 * (m ^ 2 + n ^ 2) * (t ^ 2 + 1)) / (t - 1) ^ 2 ∧
      2 * ((m + n) ^ 2 + (m - n) ^ 2) * (((t + 1) / (t - 1)) ^ 2 + 1) -
          8 * (m + n) * (m - n) * ((t + 1) / (t - 1)) =
        4 ^ 2 * (m ^ 2 + n ^ 2 * t ^ 2) / (t - 1) ^ 2 ∧
      2 * ((m + n) ^ 2 + (m - n) ^ 2) * (((t + 1) / (t - 1)) ^ 2 + 1) +
          8 * (m + n) * (m - n) * ((t + 1) / (t - 1)) =
        4 ^ 2 * (n ^ 2 + m ^ 2 * t ^ 2) / (t - 1) ^ 2 := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩ <;> field_simp <;> ring

end MagicSquareOfSquares.G1
