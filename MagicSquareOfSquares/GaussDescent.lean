import Mathlib

/-!
# Identities behind the Gaussian descent (`GAUSS_DESCENT.md`, `SECOND_DESCENT.md`)

The two descents are experiments and are not proved in Lean. This file checks only the polynomial
identities they start from. Everything is stated in a commutative ring `R` with an element `I`
such that `I * I = -1`, so it applies to `ℤ[i]` and to `ℚ(i)`; the norms are also stated for
Mathlib's `GaussianInt`.

With `t = u/v` the five forms of the G1 family are norms of the linear forms

`α0 = m v + i n u`, `α8 = n v + i m u`, `γ4 = u + i v`,
`γL = (m + i n) u - i (m - i n) v`, `γU = (m + i n) u + i (m - i n) v`:

`v² F0 = N(α0)`, `v² F8 = N(α8)`, `v² G4 = 2 (m² + n²) N(γ4)`, `v² GL = 2 N(γL)`, `v² GU = 2 N(γU)`.

* `norm_*`: the five norms.
* `support_*`: if `u, v` are coprime (`a u + b v = 1`), a common divisor of a form and its
  conjugate divides `2 m n` (for `α0`, `α8`), `2` (for `γ4`), `2 (m² - n²)` (for `γL`, `γU`).
  This is the algebraic half of the support lemma. The other half (a prime with odd exponent in a
  form whose norm is a square divides the conjugate too) uses unique factorisation in `ℤ[i]` and
  is not in Lean.
* `quartic_b_value`, `quad_plus_value`, `quartic_b_degenerate`: the identities of the second
  descent for `s = -i t`: the quartic `(m - n s)(n + m s)(1 - s²)` is `α0 · conj α8 · N(γ4) / v⁴`,
  the quadratic `q₊` is `-γL · conj γU / v²`, and at the degenerate solution `s = -i` the quartic is
  the square of `(1 - i)(m + i n)`.
* `quartic_b_frey`: the curve `Y² = (m - n s)(n + m s)(1 - s²)` maps to `y² = x (x + (m+n)²)(x - (m-n)²)`.

No local computation, no rank and no elimination of a ratio is in Lean.
-/

namespace MagicSquareOfSquares.GaussDescent

variable {R : Type*} [CommRing R]

/-! ### The five norms -/

theorem norm_alpha0 (I m n u v : R) (hI : I * I = -1) :
    (m * v + I * n * u) * (m * v - I * n * u) = m ^ 2 * v ^ 2 + n ^ 2 * u ^ 2 := by
  linear_combination (-n ^ 2 * u ^ 2) * hI

theorem norm_alpha8 (I m n u v : R) (hI : I * I = -1) :
    (n * v + I * m * u) * (n * v - I * m * u) = n ^ 2 * v ^ 2 + m ^ 2 * u ^ 2 := by
  linear_combination (-m ^ 2 * u ^ 2) * hI

theorem norm_gamma4 (I u v : R) (hI : I * I = -1) :
    (u + I * v) * (u - I * v) = u ^ 2 + v ^ 2 := by
  linear_combination (-v ^ 2) * hI

theorem norm_gammaL (I m n u v : R) (hI : I * I = -1) :
    ((m + I * n) * u - I * (m - I * n) * v) * ((m - I * n) * u + I * (m + I * n) * v)
      = (m ^ 2 + n ^ 2) * (u ^ 2 + v ^ 2) - 4 * m * n * u * v := by
  linear_combination
    (I ^ 2 * n ^ 2 * v ^ 2 - m ^ 2 * v ^ 2 + 4 * m * n * u * v - n ^ 2 * u ^ 2 - n ^ 2 * v ^ 2) * hI

theorem norm_gammaU (I m n u v : R) (hI : I * I = -1) :
    ((m + I * n) * u + I * (m - I * n) * v) * ((m - I * n) * u - I * (m + I * n) * v)
      = (m ^ 2 + n ^ 2) * (u ^ 2 + v ^ 2) + 4 * m * n * u * v := by
  linear_combination
    (I ^ 2 * n ^ 2 * v ^ 2 - m ^ 2 * v ^ 2 - 4 * m * n * u * v - n ^ 2 * u ^ 2 - n ^ 2 * v ^ 2) * hI

/-! ### The same norms in `GaussianInt` -/

theorem gaussianInt_norm_alpha0 (m n u v : ℤ) :
    Zsqrtd.norm (⟨m * v, n * u⟩ : GaussianInt) = m ^ 2 * v ^ 2 + n ^ 2 * u ^ 2 := by
  simp only [Zsqrtd.norm_def]; ring

theorem gaussianInt_norm_alpha8 (m n u v : ℤ) :
    Zsqrtd.norm (⟨n * v, m * u⟩ : GaussianInt) = n ^ 2 * v ^ 2 + m ^ 2 * u ^ 2 := by
  simp only [Zsqrtd.norm_def]; ring

theorem gaussianInt_norm_gamma4 (u v : ℤ) :
    Zsqrtd.norm (⟨u, v⟩ : GaussianInt) = u ^ 2 + v ^ 2 := by
  simp only [Zsqrtd.norm_def]; ring

theorem gaussianInt_norm_gammaL (m n u v : ℤ) :
    Zsqrtd.norm (⟨m * u - n * v, n * u - m * v⟩ : GaussianInt)
      = (m ^ 2 + n ^ 2) * (u ^ 2 + v ^ 2) - 4 * m * n * u * v := by
  simp only [Zsqrtd.norm_def]; ring

theorem gaussianInt_norm_gammaU (m n u v : ℤ) :
    Zsqrtd.norm (⟨m * u + n * v, n * u + m * v⟩ : GaussianInt)
      = (m ^ 2 + n ^ 2) * (u ^ 2 + v ^ 2) + 4 * m * n * u * v := by
  simp only [Zsqrtd.norm_def]; ring

/-! ### A common divisor of a form and its conjugate -/

theorem support_alpha0 (I m n u v a b d : R) (hI : I * I = -1) (hab : a * u + b * v = 1)
    (h1 : d ∣ m * v + I * n * u) (h2 : d ∣ m * v - I * n * u) : d ∣ 2 * m * n := by
  have h : 2 * m * n = (b * n - I * a * m) * (m * v + I * n * u)
      + (b * n + I * a * m) * (m * v - I * n * u) := by
    linear_combination (-2 * m * n) * hab + (2 * a * m * n * u) * hI
  rw [h]
  exact dvd_add (Dvd.dvd.mul_left h1 _) (Dvd.dvd.mul_left h2 _)

theorem support_alpha8 (I m n u v a b d : R) (hI : I * I = -1) (hab : a * u + b * v = 1)
    (h1 : d ∣ n * v + I * m * u) (h2 : d ∣ n * v - I * m * u) : d ∣ 2 * m * n := by
  have h := support_alpha0 I n m u v a b d hI hab h1 h2
  have e : 2 * m * n = 2 * n * m := by ring
  rw [e]; exact h

theorem support_gamma4 (I u v a b d : R) (hI : I * I = -1) (hab : a * u + b * v = 1)
    (h1 : d ∣ u + I * v) (h2 : d ∣ u - I * v) : d ∣ 2 := by
  have h : (2 : R) = (a - I * b) * (u + I * v) + (a + I * b) * (u - I * v) := by
    linear_combination (-2) * hab + (2 * b * v) * hI
  rw [h]
  exact dvd_add (Dvd.dvd.mul_left h1 _) (Dvd.dvd.mul_left h2 _)

theorem support_gammaL (I m n u v a b d : R) (hI : I * I = -1) (hab : a * u + b * v = 1)
    (h1 : d ∣ (m + I * n) * u - I * (m - I * n) * v)
    (h2 : d ∣ (m - I * n) * u + I * (m + I * n) * v) : d ∣ 2 * (m ^ 2 - n ^ 2) := by
  have h : 2 * (m ^ 2 - n ^ 2)
      = (a * (m + I * n) + b * I * (m - I * n)) * ((m + I * n) * u - I * (m - I * n) * v)
        + (a * (m - I * n) - b * I * (m + I * n)) * ((m - I * n) * u + I * (m + I * n) * v) := by
    linear_combination (-2 * m ^ 2 + 2 * n ^ 2) * hab
      + (2 * I ^ 2 * b * n ^ 2 * v - 2 * a * n ^ 2 * u + 2 * b * m ^ 2 * v - 2 * b * n ^ 2 * v) * hI
  rw [h]
  exact dvd_add (Dvd.dvd.mul_left h1 _) (Dvd.dvd.mul_left h2 _)

theorem support_gammaU (I m n u v a b d : R) (hI : I * I = -1) (hab : a * u + b * v = 1)
    (h1 : d ∣ (m + I * n) * u + I * (m - I * n) * v)
    (h2 : d ∣ (m - I * n) * u - I * (m + I * n) * v) : d ∣ 2 * (m ^ 2 - n ^ 2) := by
  have h : 2 * (m ^ 2 - n ^ 2)
      = (a * (m + I * n) - b * I * (m - I * n)) * ((m + I * n) * u + I * (m - I * n) * v)
        + (a * (m - I * n) + b * I * (m + I * n)) * ((m - I * n) * u - I * (m + I * n) * v) := by
    linear_combination (-2 * m ^ 2 + 2 * n ^ 2) * hab
      + (2 * I ^ 2 * b * n ^ 2 * v - 2 * a * n ^ 2 * u + 2 * b * m ^ 2 * v - 2 * b * n ^ 2 * v) * hI
  rw [h]
  exact dvd_add (Dvd.dvd.mul_left h1 _) (Dvd.dvd.mul_left h2 _)

/-! ### The second descent: `s = -i t`, written with `s v = -i u` -/

/-- `v⁴ Q_b(s) = α0 · conj α8 · N(γ4)` for `Q_b(s) = (m - n s)(n + m s)(1 - s²)`, `s = -i u / v`. -/
theorem quartic_b_value (I m n u v : R) (hI : I * I = -1) :
    (m * v - n * (-I * u)) * (n * v + m * (-I * u)) * (v ^ 2 - (-I * u) ^ 2)
      = (m * v + I * n * u) * (n * v - I * m * u) * (u ^ 2 + v ^ 2) := by
  linear_combination (I^2*m*n*u^4 + I*m^2*u^3*v - I*n^2*u^3*v - m*n*u^2*v^2) * hI

/-- `v² q₊(s) = -γL · conj γU` for `q₊(s) = (m² + n²) s² - 2 (m² - n²) s + (m² + n²)`. -/
theorem quad_plus_value (I m n u v : R) (hI : I * I = -1) :
    (m ^ 2 + n ^ 2) * (-I * u) ^ 2 - 2 * (m ^ 2 - n ^ 2) * (-I * u) * v + (m ^ 2 + n ^ 2) * v ^ 2
      = -(((m + I * n) * u - I * (m - I * n) * v) * ((m - I * n) * u - I * (m + I * n) * v)) := by
  linear_combination (-I^2*n^2*v^2 - 2*I*n^2*u*v + m^2*u^2 + m^2*v^2 + n^2*v^2) * hI

/-- At the degenerate solution `t = 1`, `s = -i`, the quartic `Q_b` is a square. -/
theorem quartic_b_degenerate (I m n : R) (hI : I * I = -1) :
    (m - n * (-I)) * (n + m * (-I)) * (1 - (-I) ^ 2) = ((1 - I) * (m + I * n)) ^ 2 := by
  linear_combination (I^2*m*n - I^2*n^2 + I*m^2 - 2*I*m*n + I*n^2 - m^2 + m*n) * hI

/-! ### The main curve of the second descent is a Frey curve -/

/-- If `Y² = (m - n s)(n + m s)(1 - s²)`, then `x = (n² - m²)(s + 1)/(s - 1)`, `y = 2 (n² - m²) Y/(s - 1)²` is a point of
`y² = x (x + (m + n)²)(x - (m - n)²)`. Stated with the denominators cleared. -/
theorem quartic_b_frey (m n s Y : R) (hY : Y ^ 2 = (m - n * s) * (n + m * s) * (1 - s ^ 2)) :
    (2 * (n ^ 2 - m ^ 2)) ^ 2 * Y ^ 2
      = ((n ^ 2 - m ^ 2) * (s + 1)) * ((n ^ 2 - m ^ 2) * (s + 1) + (m + n) ^ 2 * (s - 1))
        * ((n ^ 2 - m ^ 2) * (s + 1) - (m - n) ^ 2 * (s - 1)) * (s - 1) := by
  linear_combination (4 * (n ^ 2 - m ^ 2) ^ 2) * hY

end MagicSquareOfSquares.GaussDescent
