import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Ring
import Mathlib.NumberTheory.Zsqrtd.GaussianInt

/-!
# The polynomial identities of the second descent (`SECOND_DESCENT.md`, sections 1-3)

Written by an AI model (Claude); not read by a person. Only identities are proved here: no class,
no rank, no torsion and no elimination of a ratio is in Lean.

Notation of `SECOND_DESCENT.md`: a ratio `m : n`, `t = u / v`, `S = m² + n²`, `D = m² - n²`,
`s = -i t`, and over `ℤ[i]`

`α0 = m v + i n u`, `α8 = n v + i m u`, `γ4 = u + i v`,
`γL = (m + i n) u - i (m - i n) v`, `γU = (m + i n) u + i (m - i n) v`,

`L1 = m - n s`, `L2 = m + n s`, `L3 = n + m s`, `L4 = n - m s`, `L5 = 1 - s`, `L6 = 1 + s`,
`q₊ = S s² - 2 D s + S`, `q₋ = S s² + 2 D s + S`.

Everything in the first sections is stated in a commutative ring `R` with an element `I`,
`I * I = -1`, and an element `s` with `s * v = -I * u` (that is, `s = -i u / v` when `v` is
invertible), so it applies to `ℚ(i)` and, after clearing `v`, to `ℤ[i]`.

* `lin_L1` … `lin_L6`, `lin_L5L6`, `lin_qplus`, `lin_qminus`: `v · L1 = α0`, `v · L2 = conj α0`,
  `v · L3 = conj α8`, `v · L4 = α8`, `v · L5 = i · conj γ4`, `v · L6 = -i · γ4`,
  `v² L5 L6 = N(γ4)`, `v² q₊ = -γL · conj γU`, `v² q₋ = -conj γL · γU`.
* `quartic_b`, `quartic_p14`, `quartic_p23`, `quartic_m14`, `quartic_m23`: `v⁴` times the quartics
  `L1 L3 L5 L6`, `q₊ L1 L4`, `q₊ L2 L3`, `q₋ L1 L4`, `q₋ L2 L3` as products of the Gaussian linear
  forms and their conjugates. These are the facts behind "the class of `Q(s0)` is the product of the
  classes of the tuple and of their conjugates" (section 2).
* `coords_quartic_b`, `coords_quartic_p14`, `coords_quartic_p23`: the same three identities in
  `GaussianInt`, in integer coordinates: with `Q(s) = Σ q_k s^k`, the Gaussian integer
  `v⁴ Q(-i u / v) = (q0 v⁴ - q2 u² v² + q4 u⁴) + i (-q1 u v³ + q3 u³ v)` is the product.
* `degenerate_p14`, `degenerate_p23`, `degenerate_b_plus`, `class_p14_square_type`, …: the values at
  the degenerate solutions `s = ∓i` (`t = ±1`) and the classes `d` they give: `2D` for `q₊ L1 L4`
  and `q₊ L2 L3` when `S` is a square, `2mn` when `S` is twice a square (at `t = 0`; `t = ∞` is the
  case `(u, v) = (1, 0)` of the coordinate identities), `mn` for `L1 L3 L5 L6` when `S` is twice a
  square (read off `coeffs_b`), `1` when `S` is a square (`quartic_b_degenerate` in `GaussDescent.lean`).
  Each `class_*` is an identity about the value; that the value is the quartic's at the degenerate
  solution is `degenerate_*` or `coeffs_*`, and the two are composed by hand.
* The chord formulas for `x(P ± T)` are hypotheses of `add_sum`, `add_prod`, `sigma_equation` (the
  group law is not in Lean); `frey_twist` is the `K`-isomorphism `(x, y) ↦ (-x, i y)`.
* `frey_twist`: the twist by `-1` of the Frey curve `y² = x (x + a)(x - b)` is `y² = x (x - a)(x + b)`,
  by `(x, y) ↦ (-x, i y)`.
* `sigma_x`: the action of complex conjugation, `s ↦ -s`, in the coordinate `X = A/(s - r)`:
  `X' (2 r X + A) = -A X` (the map `g` of section 3; the version with divisions is
  `sigma_involution` in `G1Sigma.lean`).
* `add_sum`, `add_prod`, `sigma_equation`: for points `P = (x1, y1)`, `T = (x2, y2)` of
  `y² = x³ + a2 x² + a4 x + a6`, the sum and the product of `x(P + T)` and `x(P - T)` as rational
  functions of `x1` (`y1`, `y2` eliminated), and hence: if `x(σP) = g` is `x(P + T)` or `x(P - T)`,
  then `g² (x1 - x2)² - g · 2((x1 + x2)(x1 x2 + a4) + 2 a2 x1 x2 + 2 a6) + (x1 x2 - a4)² - 4 a6 (x1 + x2 + a2) = 0`.
  This is the polynomial equation of section 3 ("squared to remove `y`"), for `T` not of order 2
  and not `O` (`x1 ≠ x2`); `second_control.py` compared it with Sage's group law.
-/

namespace MagicSquareOfSquares.SecondDescentIdentities

variable {R : Type*} [CommRing R]

/-! ### The six linear forms and the two quadratics at `s = -i u / v` -/

theorem lin_L1 (I m n u v s : R) (hs : s * v = -I * u) : v * (m - n * s) = m * v + I * n * u := by
  linear_combination (-n) * hs

theorem lin_L2 (I m n u v s : R) (hs : s * v = -I * u) : v * (m + n * s) = m * v - I * n * u := by
  linear_combination n * hs

theorem lin_L3 (I m n u v s : R) (hs : s * v = -I * u) : v * (n + m * s) = n * v - I * m * u := by
  linear_combination m * hs

theorem lin_L4 (I m n u v s : R) (hs : s * v = -I * u) : v * (n - m * s) = n * v + I * m * u := by
  linear_combination (-m) * hs

theorem lin_L5 (I u v s : R) (hI : I * I = -1) (hs : s * v = -I * u) :
    v * (1 - s) = I * (u - I * v) := by
  linear_combination (-1 : R) * hs + v * hI

theorem lin_L6 (I u v s : R) (hI : I * I = -1) (hs : s * v = -I * u) :
    v * (1 + s) = -I * (u + I * v) := by
  linear_combination hs + v * hI

theorem lin_L5L6 (I u v s : R) (hI : I * I = -1) (hs : s * v = -I * u) :
    v ^ 2 * ((1 - s) * (1 + s)) = (u + I * v) * (u - I * v) := by
  linear_combination (I * u - s * v) * hs + (-u ^ 2 + v ^ 2) * hI

theorem lin_qplus (I m n u v s : R) (hI : I * I = -1) (hs : s * v = -I * u) :
    v ^ 2 * ((m ^ 2 + n ^ 2) * s ^ 2 - 2 * (m ^ 2 - n ^ 2) * s + (m ^ 2 + n ^ 2))
      = -(((m + I * n) * u - I * (m - I * n) * v) * ((m - I * n) * u - I * (m + I * n) * v)) := by
  linear_combination
    (-I * m ^ 2 * u - I * n ^ 2 * u + m ^ 2 * s * v - 2 * m ^ 2 * v + n ^ 2 * s * v + 2 * n ^ 2 * v) * hs
    + (-I ^ 2 * n ^ 2 * v ^ 2 - 2 * I * n ^ 2 * u * v + m ^ 2 * u ^ 2 + m ^ 2 * v ^ 2 + n ^ 2 * v ^ 2) * hI

theorem lin_qminus (I m n u v s : R) (hI : I * I = -1) (hs : s * v = -I * u) :
    v ^ 2 * ((m ^ 2 + n ^ 2) * s ^ 2 + 2 * (m ^ 2 - n ^ 2) * s + (m ^ 2 + n ^ 2))
      = -(((m - I * n) * u + I * (m + I * n) * v) * ((m + I * n) * u + I * (m - I * n) * v)) := by
  linear_combination
    (-I * m ^ 2 * u - I * n ^ 2 * u + m ^ 2 * s * v + 2 * m ^ 2 * v + n ^ 2 * s * v - 2 * n ^ 2 * v) * hs
    + (-I ^ 2 * n ^ 2 * v ^ 2 + 2 * I * n ^ 2 * u * v + m ^ 2 * u ^ 2 + m ^ 2 * v ^ 2 + n ^ 2 * v ^ 2) * hI

/-! ### The quartics of the second descent as products of the forms -/

/-- `v⁴ · L1 L3 L5 L6 = α0 · conj α8 · N(γ4)`. -/
theorem quartic_b (I m n u v s : R) (hI : I * I = -1) (hs : s * v = -I * u) :
    v ^ 4 * ((m - n * s) * (n + m * s) * ((1 - s) * (1 + s)))
      = (m * v + I * n * u) * (n * v - I * m * u) * ((u + I * v) * (u - I * v)) := by
  have e : v ^ 4 * ((m - n * s) * (n + m * s) * ((1 - s) * (1 + s)))
      = (v * (m - n * s)) * (v * (n + m * s)) * (v ^ 2 * ((1 - s) * (1 + s))) := by ring
  rw [e, lin_L1 I m n u v s hs, lin_L3 I m n u v s hs, lin_L5L6 I u v s hI hs]

/-- `v⁴ · q₊ L1 L4 = -γL · conj γU · α0 · α8`. -/
theorem quartic_p14 (I m n u v s : R) (hI : I * I = -1) (hs : s * v = -I * u) :
    v ^ 4 * (((m ^ 2 + n ^ 2) * s ^ 2 - 2 * (m ^ 2 - n ^ 2) * s + (m ^ 2 + n ^ 2)) * (m - n * s) * (n - m * s))
      = -(((m + I * n) * u - I * (m - I * n) * v) * ((m - I * n) * u - I * (m + I * n) * v))
          * (m * v + I * n * u) * (n * v + I * m * u) := by
  have e : v ^ 4 * (((m ^ 2 + n ^ 2) * s ^ 2 - 2 * (m ^ 2 - n ^ 2) * s + (m ^ 2 + n ^ 2)) * (m - n * s) * (n - m * s))
      = (v ^ 2 * ((m ^ 2 + n ^ 2) * s ^ 2 - 2 * (m ^ 2 - n ^ 2) * s + (m ^ 2 + n ^ 2)))
          * (v * (m - n * s)) * (v * (n - m * s)) := by ring
  rw [e, lin_qplus I m n u v s hI hs, lin_L1 I m n u v s hs, lin_L4 I m n u v s hs]

/-- `v⁴ · q₊ L2 L3 = -γL · conj γU · conj α0 · conj α8`. -/
theorem quartic_p23 (I m n u v s : R) (hI : I * I = -1) (hs : s * v = -I * u) :
    v ^ 4 * (((m ^ 2 + n ^ 2) * s ^ 2 - 2 * (m ^ 2 - n ^ 2) * s + (m ^ 2 + n ^ 2)) * (m + n * s) * (n + m * s))
      = -(((m + I * n) * u - I * (m - I * n) * v) * ((m - I * n) * u - I * (m + I * n) * v))
          * (m * v - I * n * u) * (n * v - I * m * u) := by
  have e : v ^ 4 * (((m ^ 2 + n ^ 2) * s ^ 2 - 2 * (m ^ 2 - n ^ 2) * s + (m ^ 2 + n ^ 2)) * (m + n * s) * (n + m * s))
      = (v ^ 2 * ((m ^ 2 + n ^ 2) * s ^ 2 - 2 * (m ^ 2 - n ^ 2) * s + (m ^ 2 + n ^ 2)))
          * (v * (m + n * s)) * (v * (n + m * s)) := by ring
  rw [e, lin_qplus I m n u v s hI hs, lin_L2 I m n u v s hs, lin_L3 I m n u v s hs]

/-- `v⁴ · q₋ L1 L4 = -conj γL · γU · α0 · α8`. -/
theorem quartic_m14 (I m n u v s : R) (hI : I * I = -1) (hs : s * v = -I * u) :
    v ^ 4 * (((m ^ 2 + n ^ 2) * s ^ 2 + 2 * (m ^ 2 - n ^ 2) * s + (m ^ 2 + n ^ 2)) * (m - n * s) * (n - m * s))
      = -(((m - I * n) * u + I * (m + I * n) * v) * ((m + I * n) * u + I * (m - I * n) * v))
          * (m * v + I * n * u) * (n * v + I * m * u) := by
  have e : v ^ 4 * (((m ^ 2 + n ^ 2) * s ^ 2 + 2 * (m ^ 2 - n ^ 2) * s + (m ^ 2 + n ^ 2)) * (m - n * s) * (n - m * s))
      = (v ^ 2 * ((m ^ 2 + n ^ 2) * s ^ 2 + 2 * (m ^ 2 - n ^ 2) * s + (m ^ 2 + n ^ 2)))
          * (v * (m - n * s)) * (v * (n - m * s)) := by ring
  rw [e, lin_qminus I m n u v s hI hs, lin_L1 I m n u v s hs, lin_L4 I m n u v s hs]

/-- `v⁴ · q₋ L2 L3 = -conj γL · γU · conj α0 · conj α8`. -/
theorem quartic_m23 (I m n u v s : R) (hI : I * I = -1) (hs : s * v = -I * u) :
    v ^ 4 * (((m ^ 2 + n ^ 2) * s ^ 2 + 2 * (m ^ 2 - n ^ 2) * s + (m ^ 2 + n ^ 2)) * (m + n * s) * (n + m * s))
      = -(((m - I * n) * u + I * (m + I * n) * v) * ((m + I * n) * u + I * (m - I * n) * v))
          * (m * v - I * n * u) * (n * v - I * m * u) := by
  have e : v ^ 4 * (((m ^ 2 + n ^ 2) * s ^ 2 + 2 * (m ^ 2 - n ^ 2) * s + (m ^ 2 + n ^ 2)) * (m + n * s) * (n + m * s))
      = (v ^ 2 * ((m ^ 2 + n ^ 2) * s ^ 2 + 2 * (m ^ 2 - n ^ 2) * s + (m ^ 2 + n ^ 2)))
          * (v * (m + n * s)) * (v * (n + m * s)) := by ring
  rw [e, lin_qminus I m n u v s hI hs, lin_L2 I m n u v s hs, lin_L3 I m n u v s hs]

/-! ### The same in `ℤ[i]`, in integer coordinates

With `Q(s) = Σ q_k s^k`, `v⁴ Q(-i u / v) = (q0 v⁴ - q2 u² v² + q4 u⁴) + i (-q1 u v³ + q3 u³ v)`.
For `L1 L3 L5 L6`: `q = (mn, m² - n², -2mn, n² - m², mn)`. -/

theorem coords_quartic_b (m n u v : ℤ) :
    (⟨m * v, n * u⟩ : GaussianInt) * ⟨n * v, -(m * u)⟩ * ⟨u ^ 2 + v ^ 2, 0⟩
      = ⟨(m * n) * v ^ 4 - (-2 * m * n) * u ^ 2 * v ^ 2 + (m * n) * u ^ 4,
          -((m ^ 2 - n ^ 2) * u * v ^ 3) + (n ^ 2 - m ^ 2) * u ^ 3 * v⟩ := by
  apply Zsqrtd.ext <;> simp only [Zsqrtd.re_mul, Zsqrtd.im_mul] <;> ring

/-- For `q₊ L1 L4`: `q0 = q4 = mn S`, `q1 = q3 = -(m⁴ + 2m³n + 2m²n² - 2mn³ + n⁴)`,
`q2 = 2(m⁴ + m³n + mn³ - n⁴)`. -/
theorem coords_quartic_p14 (m n u v : ℤ) :
    -((⟨m * u - n * v, n * u - m * v⟩ : GaussianInt) * ⟨m * u + n * v, -(n * u + m * v)⟩)
        * ⟨m * v, n * u⟩ * ⟨n * v, m * u⟩
      = ⟨(m * n * (m ^ 2 + n ^ 2)) * v ^ 4 - (2 * (m ^ 4 + m ^ 3 * n + m * n ^ 3 - n ^ 4)) * u ^ 2 * v ^ 2
            + (m * n * (m ^ 2 + n ^ 2)) * u ^ 4,
          -((-(m ^ 4 + 2 * m ^ 3 * n + 2 * m ^ 2 * n ^ 2 - 2 * m * n ^ 3 + n ^ 4)) * u * v ^ 3)
            + (-(m ^ 4 + 2 * m ^ 3 * n + 2 * m ^ 2 * n ^ 2 - 2 * m * n ^ 3 + n ^ 4)) * u ^ 3 * v⟩ := by
  apply Zsqrtd.ext <;> simp only [Zsqrtd.re_mul, Zsqrtd.im_mul, Zsqrtd.re_neg, Zsqrtd.im_neg] <;> ring

/-- For `q₊ L2 L3`: `q0 = q4 = mn S`, `q1 = q3 = m⁴ - 2m³n + 2m²n² + 2mn³ + n⁴`,
`q2 = 2(-m⁴ + m³n + mn³ + n⁴)`. -/
theorem coords_quartic_p23 (m n u v : ℤ) :
    -((⟨m * u - n * v, n * u - m * v⟩ : GaussianInt) * ⟨m * u + n * v, -(n * u + m * v)⟩)
        * ⟨m * v, -(n * u)⟩ * ⟨n * v, -(m * u)⟩
      = ⟨(m * n * (m ^ 2 + n ^ 2)) * v ^ 4 - (2 * (-m ^ 4 + m ^ 3 * n + m * n ^ 3 + n ^ 4)) * u ^ 2 * v ^ 2
            + (m * n * (m ^ 2 + n ^ 2)) * u ^ 4,
          -((m ^ 4 - 2 * m ^ 3 * n + 2 * m ^ 2 * n ^ 2 + 2 * m * n ^ 3 + n ^ 4) * u * v ^ 3)
            + (m ^ 4 - 2 * m ^ 3 * n + 2 * m ^ 2 * n ^ 2 + 2 * m * n ^ 3 + n ^ 4) * u ^ 3 * v⟩ := by
  apply Zsqrtd.ext <;> simp only [Zsqrtd.re_mul, Zsqrtd.im_mul, Zsqrtd.re_neg, Zsqrtd.im_neg] <;> ring

/-- The coefficients used above are those of the quartics (a check of the transcription). -/
theorem coeffs_p14 (m n s : R) :
    ((m ^ 2 + n ^ 2) * s ^ 2 - 2 * (m ^ 2 - n ^ 2) * s + (m ^ 2 + n ^ 2)) * (m - n * s) * (n - m * s)
      = m * n * (m ^ 2 + n ^ 2)
        + (-(m ^ 4 + 2 * m ^ 3 * n + 2 * m ^ 2 * n ^ 2 - 2 * m * n ^ 3 + n ^ 4)) * s
        + (2 * (m ^ 4 + m ^ 3 * n + m * n ^ 3 - n ^ 4)) * s ^ 2
        + (-(m ^ 4 + 2 * m ^ 3 * n + 2 * m ^ 2 * n ^ 2 - 2 * m * n ^ 3 + n ^ 4)) * s ^ 3
        + m * n * (m ^ 2 + n ^ 2) * s ^ 4 := by ring

theorem coeffs_p23 (m n s : R) :
    ((m ^ 2 + n ^ 2) * s ^ 2 - 2 * (m ^ 2 - n ^ 2) * s + (m ^ 2 + n ^ 2)) * (m + n * s) * (n + m * s)
      = m * n * (m ^ 2 + n ^ 2)
        + (m ^ 4 - 2 * m ^ 3 * n + 2 * m ^ 2 * n ^ 2 + 2 * m * n ^ 3 + n ^ 4) * s
        + (2 * (-m ^ 4 + m ^ 3 * n + m * n ^ 3 + n ^ 4)) * s ^ 2
        + (m ^ 4 - 2 * m ^ 3 * n + 2 * m ^ 2 * n ^ 2 + 2 * m * n ^ 3 + n ^ 4) * s ^ 3
        + m * n * (m ^ 2 + n ^ 2) * s ^ 4 := by ring

theorem coeffs_b (m n s : R) :
    (m - n * s) * (n + m * s) * ((1 - s) * (1 + s))
      = m * n + (m ^ 2 - n ^ 2) * s + (-2 * m * n) * s ^ 2 + (n ^ 2 - m ^ 2) * s ^ 3 + m * n * s ^ 4 := by
  ring

/-! ### Values at the degenerate solutions, and the class `d` -/

theorem degenerate_p14 (I m n : R) (hI : I * I = -1) :
    ((m ^ 2 + n ^ 2) * I ^ 2 - 2 * (m ^ 2 - n ^ 2) * I + (m ^ 2 + n ^ 2)) * (m - n * I) * (n - m * I)
        = -2 * (m ^ 2 - n ^ 2) * (m ^ 2 + n ^ 2)
      ∧ ((m ^ 2 + n ^ 2) * (-I) ^ 2 - 2 * (m ^ 2 - n ^ 2) * (-I) + (m ^ 2 + n ^ 2)) * (m - n * (-I)) * (n - m * (-I))
        = -2 * (m ^ 2 - n ^ 2) * (m ^ 2 + n ^ 2) := by
  constructor
  · linear_combination (I ^ 2 * m ^ 3 * n + I ^ 2 * m * n ^ 3 - I * m ^ 4 - 2 * I * m ^ 3 * n - 2 * I * m ^ 2 * n ^ 2
      + 2 * I * m * n ^ 3 - I * n ^ 4 + 2 * m ^ 4 + m ^ 3 * n + m * n ^ 3 - 2 * n ^ 4) * hI
  · linear_combination (I ^ 2 * m ^ 3 * n + I ^ 2 * m * n ^ 3 + I * m ^ 4 + 2 * I * m ^ 3 * n + 2 * I * m ^ 2 * n ^ 2
      - 2 * I * m * n ^ 3 + I * n ^ 4 + 2 * m ^ 4 + m ^ 3 * n + m * n ^ 3 - 2 * n ^ 4) * hI

theorem degenerate_p23 (I m n : R) (hI : I * I = -1) :
    ((m ^ 2 + n ^ 2) * I ^ 2 - 2 * (m ^ 2 - n ^ 2) * I + (m ^ 2 + n ^ 2)) * (m + n * I) * (n + m * I)
        = 2 * (m ^ 2 - n ^ 2) * (m ^ 2 + n ^ 2)
      ∧ ((m ^ 2 + n ^ 2) * (-I) ^ 2 - 2 * (m ^ 2 - n ^ 2) * (-I) + (m ^ 2 + n ^ 2)) * (m + n * (-I)) * (n + m * (-I))
        = 2 * (m ^ 2 - n ^ 2) * (m ^ 2 + n ^ 2) := by
  constructor
  · linear_combination (I ^ 2 * m ^ 3 * n + I ^ 2 * m * n ^ 3 + I * m ^ 4 - 2 * I * m ^ 3 * n + 2 * I * m ^ 2 * n ^ 2
      + 2 * I * m * n ^ 3 + I * n ^ 4 - 2 * m ^ 4 + m ^ 3 * n + m * n ^ 3 + 2 * n ^ 4) * hI
  · linear_combination (I ^ 2 * m ^ 3 * n + I ^ 2 * m * n ^ 3 - I * m ^ 4 + 2 * I * m ^ 3 * n - 2 * I * m ^ 2 * n ^ 2
      - 2 * I * m * n ^ 3 - I * n ^ 4 - 2 * m ^ 4 + m ^ 3 * n + m * n ^ 3 + 2 * n ^ 4) * hI

/-- The other degenerate solution of square type, `t = -1`, `s = i` (the conjugate of
`quartic_b_degenerate` of `GaussDescent.lean`): there too `L1 L3 L5 L6` is a square. -/
theorem degenerate_b_plus (I m n : R) (hI : I * I = -1) :
    (m - n * I) * (n + m * I) * (1 - I ^ 2) = ((1 + I) * (m - I * n)) ^ 2 := by
  linear_combination (I ^ 2 * m * n - I ^ 2 * n ^ 2 - I * m ^ 2 + 2 * I * m * n - I * n ^ 2 - m ^ 2 + m * n) * hI

/-- Square type, `m² + n² = k²`: at both degenerate solutions `2 (m² - n²) · q₊ L1 L4` is the
square of `2 i (m² - n²) k`, so the class of `q₊ L1 L4` there is that of `d = 2 (m² - n²)`. -/
theorem class_p14_square_type (I m n k : R) (hI : I * I = -1) (hk : m ^ 2 + n ^ 2 = k ^ 2) :
    2 * (m ^ 2 - n ^ 2) * (-2 * (m ^ 2 - n ^ 2) * (m ^ 2 + n ^ 2)) = (2 * I * (m ^ 2 - n ^ 2) * k) ^ 2 := by
  linear_combination (-4 * (m ^ 2 - n ^ 2) ^ 2) * hk + (-4 * (m ^ 2 - n ^ 2) ^ 2 * k ^ 2) * hI

/-- Square type: `2 (m² - n²) · q₊ L2 L3 = (2 (m² - n²) k)²` at both degenerate solutions. -/
theorem class_p23_square_type (m n k : R) (hk : m ^ 2 + n ^ 2 = k ^ 2) :
    2 * (m ^ 2 - n ^ 2) * (2 * (m ^ 2 - n ^ 2) * (m ^ 2 + n ^ 2)) = (2 * (m ^ 2 - n ^ 2) * k) ^ 2 := by
  linear_combination (4 * (m ^ 2 - n ^ 2) ^ 2) * hk

/-- Twice-square type, `m² + n² = 2k²`: at `t = 0` (`s = 0`) and at `t = ∞` (leading
coefficient) both `q₊ L1 L4` and `q₊ L2 L3` take the value `m n (m² + n²)` (`coeffs_p14`,
`coeffs_p23`), and `2 m n` times it is the square of `2 k m n`: the class is that of `d = 2 m n`. -/
theorem class_p_twice_square_type (m n k : R) (hk : m ^ 2 + n ^ 2 = 2 * k ^ 2) :
    2 * m * n * (m * n * (m ^ 2 + n ^ 2)) = (2 * k * m * n) ^ 2 := by
  linear_combination (2 * m ^ 2 * n ^ 2) * hk

-- Twice-square type, `L1 L3 L5 L6`: the value at `s = 0` and the leading coefficient are both `m n`
-- (`coeffs_b`, `coords_quartic_b` at `(u, v) = (0, 1), (1, 0)`), so the class there is that of `d = m n`.
-- (An earlier version stated this as the tautology `m n · m n = (m n)²`; removed after review.)

/-! ### The twist by `-1` of the Frey curve -/

theorem frey_twist (I a b x y : R) (hI : I * I = -1) (h : y ^ 2 = x * (x + a) * (x - b)) :
    (I * y) ^ 2 = (-x) * ((-x) - a) * ((-x) + b) := by
  linear_combination y ^ 2 * hI - h

/-! ### Complex conjugation in the coordinate `X = A / (s - r)` -/

/-- If `X (s - r) = A` and `X' (-s - r) = A` (the point and its conjugate, `σ s = -s`), then
`X' (2 r X + A) = -A X`: `X' = g(X)`, `g(x) = -A x / (2 r x + A)`. -/
theorem sigma_x (A r s X X' : R) (h1 : X * (s - r) = A) (h2 : X' * (-s - r) = A) :
    X' * (2 * r * X + A) = -A * X := by
  linear_combination (-X') * h1 - X * h2

/-! ### The equation of section 3: `x(σP) = x(P ± T)` with `y` removed

For `P = (x1, y1)`, `T = (x2, y2)` on `y² = x³ + a2 x² + a4 x + a6` with `x1 ≠ x2`, the addition
law gives `x(P + T) (x2 - x1)² = (y2 - y1)² - (a2 + x1 + x2)(x2 - x1)²` and
`x(P - T) (x2 - x1)² = (y2 + y1)² - (a2 + x1 + x2)(x2 - x1)²` (slope `(∓y2 - y1)/(x2 - x1)`). -/

theorem add_sum (a2 a4 a6 x1 y1 x2 y2 x3 x4 : R)
    (h1 : y1 ^ 2 = x1 ^ 3 + a2 * x1 ^ 2 + a4 * x1 + a6) (h2 : y2 ^ 2 = x2 ^ 3 + a2 * x2 ^ 2 + a4 * x2 + a6)
    (h3 : x3 * (x2 - x1) ^ 2 = (y2 - y1) ^ 2 - (a2 + x1 + x2) * (x2 - x1) ^ 2)
    (h4 : x4 * (x2 - x1) ^ 2 = (y2 + y1) ^ 2 - (a2 + x1 + x2) * (x2 - x1) ^ 2) :
    (x3 + x4) * (x2 - x1) ^ 2 = 2 * ((x1 + x2) * (x1 * x2 + a4) + 2 * a2 * x1 * x2 + 2 * a6) := by
  linear_combination h3 + h4 + 2 * h1 + 2 * h2

theorem add_prod (a2 a4 a6 x1 y1 x2 y2 x3 x4 : R)
    (h1 : y1 ^ 2 = x1 ^ 3 + a2 * x1 ^ 2 + a4 * x1 + a6) (h2 : y2 ^ 2 = x2 ^ 3 + a2 * x2 ^ 2 + a4 * x2 + a6)
    (h3 : x3 * (x2 - x1) ^ 2 = (y2 - y1) ^ 2 - (a2 + x1 + x2) * (x2 - x1) ^ 2)
    (h4 : x4 * (x2 - x1) ^ 2 = (y2 + y1) ^ 2 - (a2 + x1 + x2) * (x2 - x1) ^ 2) :
    x3 * x4 * (x2 - x1) ^ 4 = ((x1 * x2 - a4) ^ 2 - 4 * a6 * (x1 + x2 + a2)) * (x2 - x1) ^ 2 := by
  linear_combination (x4 * (x2 - x1) ^ 2) * h3 + ((y2 - y1) ^ 2 - (a2 + x1 + x2) * (x2 - x1) ^ 2) * h4
    + (-a2 * x1 ^ 2 + 4 * a2 * x1 * x2 - 2 * a2 * x2 ^ 2 + a4 * x1 + a6 - x1 ^ 3 + 2 * x1 ^ 2 * x2
        + 2 * x1 * x2 ^ 2 - 2 * x2 ^ 3 + y1 ^ 2 - 2 * y2 ^ 2) * h1
    + (-4 * a2 * x1 ^ 2 + 4 * a2 * x1 * x2 - a2 * x2 ^ 2 - 2 * a4 * x1 + a4 * x2 - a6 - 4 * x1 ^ 3
        + 2 * x1 ^ 2 * x2 + 2 * x1 * x2 ^ 2 - x2 ^ 3 + y2 ^ 2) * h2

/-- The polynomial equation of section 3: if `g = x(σP)` equals `x(P + T)` or `x(P - T)`, then
`g` and `x1 = x(P)` satisfy an equation with no `y` in it. -/
theorem sigma_equation {F : Type*} [Field F] (a2 a4 a6 x1 y1 x2 y2 x3 x4 g : F)
    (h1 : y1 ^ 2 = x1 ^ 3 + a2 * x1 ^ 2 + a4 * x1 + a6) (h2 : y2 ^ 2 = x2 ^ 3 + a2 * x2 ^ 2 + a4 * x2 + a6)
    (h3 : x3 * (x2 - x1) ^ 2 = (y2 - y1) ^ 2 - (a2 + x1 + x2) * (x2 - x1) ^ 2)
    (h4 : x4 * (x2 - x1) ^ 2 = (y2 + y1) ^ 2 - (a2 + x1 + x2) * (x2 - x1) ^ 2)
    (hx : x2 - x1 ≠ 0) (hg : g = x3 ∨ g = x4) :
    g ^ 2 * (x2 - x1) ^ 2 - g * (2 * ((x1 + x2) * (x1 * x2 + a4) + 2 * a2 * x1 * x2 + 2 * a6))
      + ((x1 * x2 - a4) ^ 2 - 4 * a6 * (x1 + x2 + a2)) = 0 := by
  have hs := add_sum a2 a4 a6 x1 y1 x2 y2 x3 x4 h1 h2 h3 h4
  have hp := add_prod a2 a4 a6 x1 y1 x2 y2 x3 x4 h1 h2 h3 h4
  have h0 : (g - x3) * (g - x4) = 0 := by
    rcases hg with rfl | rfl <;> ring
  have hd : (x2 - x1) ^ 2 ≠ 0 := pow_ne_zero 2 hx
  have key : (x2 - x1) ^ 2 * (g ^ 2 * (x2 - x1) ^ 2 - g * (2 * ((x1 + x2) * (x1 * x2 + a4) + 2 * a2 * x1 * x2 + 2 * a6))
      + ((x1 * x2 - a4) ^ 2 - 4 * a6 * (x1 + x2 + a2))) = (g - x3) * (g - x4) * (x2 - x1) ^ 4 := by
    linear_combination (g * (x2 - x1) ^ 2) * hs - hp
  rw [h0, zero_mul] at key
  exact (mul_eq_zero.mp key).resolve_left hd

end MagicSquareOfSquares.SecondDescentIdentities

