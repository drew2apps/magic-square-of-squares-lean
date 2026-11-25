import MagicSquareOfSquares.G1Sigma

/-!
# The curves of `F0·F8·GL` at every ratio: the algebra behind `UNIFORM_RANK.md`

Written by an AI model; not read by a person; not compiled by its author (no Lean installation was available there); compiled by the conductor of the run, after
removing a `genC4` from two `unfold` lists. Every identity below is also checked by sympy in
`scripts/g1_full/sigma_generic.py`, Part A.

For arbitrary rational `m, n`, with `S = m² + n²` and `H(m, n) = m⁴ − 2m³n + 2m²n² + 2mn³ + n⁴`:

* the quartic `(z + 2)·g(z)` of `E⁺` in `h = z + 2`, and its Weierstrass model `y² = X³ + a₂X² + a₄X + a₆`
  (`gen_quartic`, by `quartic_root_weierstrass` of `G1Sigma.lean`);
* the cubic factors as `(X + 2S³)(X² + 16m²n²(m+n)²X + 16m²n²(m+n)⁴S²)` (`gen_cubic_factor`), and the quadratic
  factor has discriminant `−64m²n²(m − n)²(m + n)⁶` (`gen_disc`);
* `P = (−(m+n)²S², S³(m − n)(m + n)²)` and `R = (−4mnS(m² − n²), −8mn(m+n)S·H)` lie on it (`gen_P_on`, `gen_R_on`),
  and the tangent at `R` meets the curve again at `(x(P), −y(P))` (`gen_tangent_x`, `gen_tangent_y`), so `2R = P`;
* the values of the 2-isogeny map: `x(R) + 2S³ = 2S·H` (`gen_alpha_R`), `f′(−2S³) = 4S²(m⁸ + 14m⁴n⁴ + n⁸)`
  (`gen_alpha_T`), `H(m, n)·H(m, −n) = m⁸ + 14m⁴n⁴ + n⁸` (`gen_H_mul`), and the value at `X = 0` is
  `2S·(4mn(m+n)²S²)²` (`gen_Q_on`);
* the quartic of `E⁻` at `(m, n)` is the quartic of `E⁺` at `(m, −n)` after `z ↦ −z` (`gen_minus_plus`);
* the duplication identity `f′² − 4f·(2x − e₂ − e₃) = ((x − e₁)² − (e₁ − e₂)(e₁ − e₃))²` (`dup_identity`), which says
  that `x(2P) − e₁` is a square in the field of `e₁`;
* the map of Lemma B: a solution of `X⁴ + 14X²Y² + Y⁴ = Z²` with `Y ≠ 0` gives the point `(X²/Y², XZ/Y³)` of
  `y² = x³ + 14x² + x` (`lemmaB_identity`, in homogeneous form).

The descents (Lemma B, Lemma D), the torsion bounds and the local tests are not in Lean.
-/

namespace MagicSquareOfSquares.G1

/-- `a₂` of the model of `E⁺` (`sigma_exact.model(m, n, 2)`). -/
def genA2 (m n : ℚ) : ℚ := 2 * (m ^ 6 + 11 * m ^ 4 * n ^ 2 + 16 * m ^ 3 * n ^ 3 + 11 * m ^ 2 * n ^ 4 + n ^ 6)
/-- `c₁ = A = q′(−2)`. -/
def genA (m n : ℚ) : ℚ := -4 * (m + n) ^ 2 * (m ^ 2 + n ^ 2) ^ 2
/-- `c₃` of the quartic in `h = z + 2`. -/
def genC3 (m n : ℚ) : ℚ := -4 * m ^ 2 * n ^ 2 * (3 * m ^ 2 + 2 * m * n + 3 * n ^ 2)
/-- `c₄`, the leading coefficient `2m²n²S` of the quartic. -/
def genC4 (m n : ℚ) : ℚ := 2 * m ^ 2 * n ^ 2 * (m ^ 2 + n ^ 2)
/-- the cubic of the model: `X³ + c₂X² + c₁c₃X + c₁²c₄`. -/
def genCub (m n X : ℚ) : ℚ :=
  X ^ 3 + genA2 m n * X ^ 2 + genA m n * genC3 m n * X + genA m n ^ 2 * genC4 m n
/-- its derivative in `X`. -/
def genCubD (m n X : ℚ) : ℚ := 3 * X ^ 2 + 2 * genA2 m n * X + genA m n * genC3 m n
/-- `H(m, n) = m⁴ − 2m³n + 2m²n² + 2mn³ + n⁴`. -/
def genH (m n : ℚ) : ℚ := m ^ 4 - 2 * m ^ 3 * n + 2 * m ^ 2 * n ^ 2 + 2 * m * n ^ 3 + n ^ 4

theorem gen_quartic (m n z : ℚ) :
    (z + 2) * sigmaG m n z =
      genA m n * (z + 2) + genA2 m n * (z + 2) ^ 2 + genC3 m n * (z + 2) ^ 3 + genC4 m n * (z + 2) ^ 4 := by
  unfold sigmaG genA genA2 genC3 genC4
  ring

theorem gen_cubic_factor (m n X : ℚ) :
    genCub m n X = (X + 2 * (m ^ 2 + n ^ 2) ^ 3) *
      (X ^ 2 + 16 * m ^ 2 * n ^ 2 * (m + n) ^ 2 * X + 16 * m ^ 2 * n ^ 2 * (m + n) ^ 4 * (m ^ 2 + n ^ 2) ^ 2) := by
  unfold genCub genA2 genA genC3 genC4
  ring

theorem gen_disc (m n : ℚ) :
    (16 * m ^ 2 * n ^ 2 * (m + n) ^ 2) ^ 2 - 4 * (16 * m ^ 2 * n ^ 2 * (m + n) ^ 4 * (m ^ 2 + n ^ 2) ^ 2) =
      -64 * m ^ 2 * n ^ 2 * (m - n) ^ 2 * (m + n) ^ 6 := by
  ring

theorem gen_P_on (m n : ℚ) :
    ((m ^ 2 + n ^ 2) ^ 3 * (m - n) * (m + n) ^ 2) ^ 2 = genCub m n (-((m + n) ^ 2 * (m ^ 2 + n ^ 2) ^ 2)) := by
  unfold genCub genA2 genA genC3 genC4
  ring

theorem gen_R_on (m n : ℚ) :
    (-8 * m * n * (m + n) * (m ^ 2 + n ^ 2) * genH m n) ^ 2 =
      genCub m n (-4 * m * n * (m ^ 2 + n ^ 2) * (m ^ 2 - n ^ 2)) := by
  unfold genCub genA2 genA genC3 genC4 genH
  ring

set_option maxHeartbeats 1000000 in
/-- The tangent at `R` (slope `f′(x_R)/(2y_R)`) meets the curve again at abscissa `x(P)`:
`λ² = x_P + a₂ + 2x_R`, cleared of denominators. -/
theorem gen_tangent_x (m n : ℚ) :
    genCubD m n (-4 * m * n * (m ^ 2 + n ^ 2) * (m ^ 2 - n ^ 2)) ^ 2 =
      4 * (-8 * m * n * (m + n) * (m ^ 2 + n ^ 2) * genH m n) ^ 2 *
        (-((m + n) ^ 2 * (m ^ 2 + n ^ 2) ^ 2) + genA2 m n + 2 * (-4 * m * n * (m ^ 2 + n ^ 2) * (m ^ 2 - n ^ 2))) := by
  unfold genCubD genA2 genA genC3 genH
  ring

set_option maxHeartbeats 1000000 in
/-- ... and at the point `(x_P, −y_P)`: `y(2R) = −(y_R + λ(x_P − x_R)) = y_P`, cleared of denominators. -/
theorem gen_tangent_y (m n : ℚ) :
    2 * (-8 * m * n * (m + n) * (m ^ 2 + n ^ 2) * genH m n) * ((m ^ 2 + n ^ 2) ^ 3 * (m - n) * (m + n) ^ 2) +
      2 * (-8 * m * n * (m + n) * (m ^ 2 + n ^ 2) * genH m n) ^ 2 +
      genCubD m n (-4 * m * n * (m ^ 2 + n ^ 2) * (m ^ 2 - n ^ 2)) *
        (-((m + n) ^ 2 * (m ^ 2 + n ^ 2) ^ 2) - (-4 * m * n * (m ^ 2 + n ^ 2) * (m ^ 2 - n ^ 2))) = 0 := by
  unfold genCubD genA2 genA genC3 genH
  ring

theorem gen_alpha_R (m n : ℚ) :
    -4 * m * n * (m ^ 2 + n ^ 2) * (m ^ 2 - n ^ 2) + 2 * (m ^ 2 + n ^ 2) ^ 3 = 2 * (m ^ 2 + n ^ 2) * genH m n := by
  unfold genH
  ring

theorem gen_alpha_T (m n : ℚ) :
    genCubD m n (-2 * (m ^ 2 + n ^ 2) ^ 3) = 4 * (m ^ 2 + n ^ 2) ^ 2 * (m ^ 8 + 14 * m ^ 4 * n ^ 4 + n ^ 8) := by
  unfold genCubD genA2 genA genC3
  ring

theorem gen_H_mul (m n : ℚ) : genH m n * genH m (-n) = m ^ 8 + 14 * m ^ 4 * n ^ 4 + n ^ 8 := by
  unfold genH
  ring

theorem gen_Q_on (m n : ℚ) :
    genCub m n 0 = 2 * (m ^ 2 + n ^ 2) * (4 * m * n * (m + n) ^ 2 * (m ^ 2 + n ^ 2) ^ 2) ^ 2 := by
  unfold genCub genA genC4
  ring

theorem gen_minus_plus (m n z : ℚ) : (z - 2) * sigmaG m n z = (-z + 2) * sigmaG m (-n) (-z) := by
  unfold sigmaG
  ring

theorem dup_identity (x e₁ e₂ e₃ : ℚ) :
    ((x - e₂) * (x - e₃) + (x - e₁) * (x - e₃) + (x - e₁) * (x - e₂)) ^ 2 -
        4 * ((x - e₁) * (x - e₂) * (x - e₃)) * (2 * x - e₂ - e₃) =
      ((x - e₁) ^ 2 - (e₁ - e₂) * (e₁ - e₃)) ^ 2 := by
  ring

/-- The map of Lemma B in homogeneous form: `(x : y : 1) = (X²Y : XZ : Y³)` on `y² = x³ + 14x² + x`. -/
theorem lemmaB_identity (X Y Z : ℚ) :
    (X * Z) ^ 2 - ((X ^ 2) ^ 3 + 14 * (X ^ 2) ^ 2 * Y ^ 2 + X ^ 2 * Y ^ 4) =
      X ^ 2 * (Z ^ 2 - (X ^ 4 + 14 * X ^ 2 * Y ^ 2 + Y ^ 4)) := by
  ring

end MagicSquareOfSquares.G1

