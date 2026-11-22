import Mathlib.Tactic
import Mathlib.RingTheory.Coprime.Basic

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
  form by `4` (`forms_double`). The partner pair is again coprime (`partner_coprime` when `m + n`
  is odd, `partner_half_coprime` for the halved pair).
* The quartic `v² = aτ⁴ + bτ³ + cτ² + dτ + q²` and its Weierstrass model, by the formulas
  `dm3.py` uses (`quartic_to_weierstrass`). For `E0L`, `τ = t − m/n` and `q = 2m(m² − n²)/n`.
  The further change to a minimal model is made by Sage and is not in Lean.
* The octic `m⁸ + 14m⁴n⁴ + n⁸` of `G1_FULL_SYSTEM.md`, Result 1: it is a quarter of the
  resultant of `F0` and `GL` (`resultant_F0_GL`), and modulo it `S²·F0·F8·GL·GU` is a square
  (`octic_square`).
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

/-- The five forms of `(2m, 2n)` are `4` times the five forms of `(m, n)`. -/
theorem forms_double (m n t : ℚ) :
    (2 * m) ^ 2 + (2 * n) ^ 2 * t ^ 2 = 4 * (m ^ 2 + n ^ 2 * t ^ 2) ∧
    (2 * n) ^ 2 + (2 * m) ^ 2 * t ^ 2 = 4 * (n ^ 2 + m ^ 2 * t ^ 2) ∧
    2 * ((2 * m) ^ 2 + (2 * n) ^ 2) * (t ^ 2 + 1) = 4 * (2 * (m ^ 2 + n ^ 2) * (t ^ 2 + 1)) ∧
    2 * ((2 * m) ^ 2 + (2 * n) ^ 2) * (t ^ 2 + 1) - 8 * (2 * m) * (2 * n) * t =
      4 * (2 * (m ^ 2 + n ^ 2) * (t ^ 2 + 1) - 8 * m * n * t) ∧
    2 * ((2 * m) ^ 2 + (2 * n) ^ 2) * (t ^ 2 + 1) + 8 * (2 * m) * (2 * n) * t =
      4 * (2 * (m ^ 2 + n ^ 2) * (t ^ 2 + 1) + 8 * m * n * t) := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩ <;> ring

/-- If `m`, `n` are coprime and `m + n` is odd, the partner pair `(m + n, m − n)` is coprime. -/
theorem partner_coprime (m n : ℤ) (h : IsCoprime m n) (hodd : Odd (m + n)) :
    IsCoprime (m + n) (m - n) := by
  obtain ⟨a, b, hab⟩ := h
  obtain ⟨k, hk⟩ := hodd
  exact ⟨1 - k * (a + b), -(k * (a - b)), by linear_combination hk - 2 * k * hab⟩

/-- If `M + N` and `M − N` are coprime, so are `M` and `N`. With `M = (m + n)/2`, `N = (m − n)/2`
for odd coprime `m`, `n` this says that the halved partner pair is coprime. -/
theorem partner_half_coprime (M N : ℤ) (h : IsCoprime (M + N) (M - N)) : IsCoprime M N := by
  obtain ⟨a, b, hab⟩ := h
  exact ⟨a + b, a - b, by linear_combination hab⟩

/-- A rational point `(τ, v)`, `τ ≠ 0`, of the quartic `v² = aτ⁴ + bτ³ + cτ² + dτ + q²` gives the
point `(X, Y)` of the Weierstrass curve
`Y² + (d/q)XY + 2qbY = X³ + (c − d²/(4q²))X² − 4q²aX − 4q²a(c − d²/(4q²))`. -/
theorem quartic_to_weierstrass (a b c d q τ v X Y : ℚ) (hq : q ≠ 0) (hτ : τ ≠ 0)
    (h : v ^ 2 = a * τ ^ 4 + b * τ ^ 3 + c * τ ^ 2 + d * τ + q ^ 2)
    (hX : X = (2 * q * (v + q) + d * τ) / τ ^ 2)
    (hY : Y = (4 * q ^ 2 * (v + q) + 2 * q * (d * τ + c * τ ^ 2) - d ^ 2 * τ ^ 2 / (2 * q)) / τ ^ 3) :
    Y ^ 2 + d / q * X * Y + 2 * q * b * Y =
      X ^ 3 + (c - d ^ 2 / (4 * q ^ 2)) * X ^ 2 - 4 * q ^ 2 * a * X -
        4 * q ^ 2 * a * (c - d ^ 2 / (4 * q ^ 2)) := by
  have h0 : v ^ 2 - (a * τ ^ 4 + b * τ ^ 3 + c * τ ^ 2 + d * τ + q ^ 2) = 0 := by rw [h]; ring
  have hK : Y ^ 2 + d / q * X * Y + 2 * q * b * Y -
      (X ^ 3 + (c - d ^ 2 / (4 * q ^ 2)) * X ^ 2 - 4 * q ^ 2 * a * X -
        4 * q ^ 2 * a * (c - d ^ 2 / (4 * q ^ 2))) =
      -(4 * c * q ^ 2 * τ ^ 2 - d ^ 2 * τ ^ 2 + 4 * d * q ^ 2 * τ + 8 * q ^ 4 + 8 * q ^ 3 * v) / τ ^ 6 *
        (v ^ 2 - (a * τ ^ 4 + b * τ ^ 3 + c * τ ^ 2 + d * τ + q ^ 2)) := by
    rw [hX, hY]
    field_simp
    ring
  rw [h0, mul_zero] at hK
  exact sub_eq_zero.mp hK

/-- The resultant of `F0 = n²t² + m²` and `GL = 2St² − 8mnt + 2S` (for quadratics
`at² + bt + c`, `At² + Bt + C` it is `(aC − cA)² − (aB − bA)(bC − cB)`) is four times the octic. -/
theorem resultant_F0_GL (m n : ℚ) :
    (n ^ 2 * (2 * (m ^ 2 + n ^ 2)) - m ^ 2 * (2 * (m ^ 2 + n ^ 2))) ^ 2 -
      (n ^ 2 * (-8 * m * n) - 0 * (2 * (m ^ 2 + n ^ 2))) *
        (0 * (2 * (m ^ 2 + n ^ 2)) - m ^ 2 * (-8 * m * n)) =
    4 * (m ^ 8 + 14 * m ^ 4 * n ^ 4 + n ^ 8) := by
  ring

/-- `S²·F0·F8·GL·GU` is a square plus a multiple of the octic. So modulo a prime dividing the
octic and not dividing `S`, the product `F0·F8·GL·GU` is a square for every `t`. -/
theorem octic_square (m n t : ℚ) :
    (m ^ 2 + n ^ 2) ^ 2 * ((m ^ 2 + n ^ 2 * t ^ 2) * (n ^ 2 + m ^ 2 * t ^ 2) *
      (2 * (m ^ 2 + n ^ 2) * (t ^ 2 + 1) - 8 * m * n * t) *
      (2 * (m ^ 2 + n ^ 2) * (t ^ 2 + 1) + 8 * m * n * t)) =
    (2 * m * n * ((m ^ 2 + n ^ 2) ^ 2 * (t ^ 4 + 1) +
      (2 * (m ^ 2 + n ^ 2) ^ 2 - 16 * m ^ 2 * n ^ 2) * t ^ 2)) ^ 2 +
    4 * (m ^ 8 + 14 * m ^ 4 * n ^ 4 + n ^ 8) * t ^ 2 * ((m ^ 2 + n ^ 2) ^ 2 * (t ^ 4 + 1) +
      (2 * (m ^ 2 + n ^ 2) ^ 2 - 16 * m ^ 2 * n ^ 2) * t ^ 2) := by
  ring

end MagicSquareOfSquares.G1
