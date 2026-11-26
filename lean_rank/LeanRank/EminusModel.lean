import LeanRank.UniformRank

/-!
# `E⁻` is `Eplus m (−n)`: the Weierstrass models of `UNIFORM_RANK.md` §1, built from the quartic

Written by an AI model; not read by a person.

`UNIFORM_RANK.md` §1 defines `E^±` as the quotients `Y² = q(z)` of `y² = F0·F8·GL`, with
`q(z) = (z ± 2)·g(z)`, `g(z) = (2Sz − 8mn)(m²n²z² + (m² − n²)²)`, `S = m² + n²`, and puts them in Weierstrass
form by the recipe: `z0 = ∓2` (a root of `q`), `A = q′(z0)`,

    a2 = q″(z0)/2,   a4 = A·q‴(z0)/6,   a6 = A²·q⁗(z0)/24,

with the map `(z, Y) ↦ (A/(z − z0), A·Y/(z − z0)²)`. Until now the library used `TwoDescent.Eplus m n` for `E⁺`
and `TwoDescent.Eplus m (−n)` for `E⁻`, the second identification being a sympy check (A2). Here:

* `weierstrassOf q z0`: the recipe, for an arbitrary polynomial `q` over `ℚ` and a point `z0`, with Mathlib's
  `Polynomial.derivative` and `Polynomial.eval` (nothing about the answer is pasted in);
* `qPoly s m n = (X + s)·g`, so `E⁺` is `weierstrassOf (qPoly 2 m n) (−2)` (`EplusModel`) and `E⁻` is
  `weierstrassOf (qPoly (−2) m n) 2` (`Eminus`);
* `EplusModel_eq : EplusModel m n = Eplus m n` and `Eminus_eq : Eminus m n = Eplus m (−n)`, equalities of
  Weierstrass curves over `ℚ` (all five coefficients), for all rational `m, n`;
* `EplusModel_map`, `Eminus_map`: the map sends every point `(z, Y)` of `Y² = q(z)` with `z ≠ z0` to a point of
  the model; `sextic_to_Eplus`, `sextic_to_Eminus`: composed with `G1.sigma_plus`, `G1.sigma_minus`
  (`G1Sigma.lean`), every point `(t, y)` of `y² = F0·F8·GL` with `t ≠ 0, ∓1` lands on `E^±`;
* `theorem1_Eminus`, `theorem2_Eminus`: Theorems 1 and 2 of `UniformRank.lean` restated with these names
  (coprime pairs; for all rational pairs see `UniformRankAll.lean`).
-/

namespace MagicSquareOfSquares.EminusModel

open Polynomial MagicSquareOfSquares.TwoDescent

/-- **The recipe of `UNIFORM_RANK.md` §1**, for any polynomial `q` over `ℚ` and any `z0`:
`y² = x³ + a2·x² + a4·x + a6` with `A = q′(z0)`, `a2 = q″(z0)/2`, `a4 = A·q‴(z0)/6`, `a6 = A²·q⁗(z0)/24`
(`a1 = a3 = 0`). It is meant for a quartic `q` (linear times cubic) with `q(z0) = 0`. -/
noncomputable def weierstrassOf (q : ℚ[X]) (z0 : ℚ) : WeierstrassCurve ℚ where
  a₁ := 0
  a₂ := (derivative (derivative q)).eval z0 / 2
  a₃ := 0
  a₄ := (derivative q).eval z0 * (derivative (derivative (derivative q))).eval z0 / 6
  a₆ := (derivative q).eval z0 ^ 2 *
    (derivative (derivative (derivative (derivative q)))).eval z0 / 24

/-- `A = q′(z0)`, the constant of the map. -/
noncomputable def mapA (q : ℚ[X]) (z0 : ℚ) : ℚ := (derivative q).eval z0

/-- `g(z) = (2Sz − 8mn)(m²n²z² + (m² − n²)²)`, `S = m² + n²`. -/
noncomputable def gPoly (m n : ℚ) : ℚ[X] :=
  (C (2 * (m ^ 2 + n ^ 2)) * X - C (8 * m * n)) * (C (m ^ 2 * n ^ 2) * X ^ 2 + C ((m ^ 2 - n ^ 2) ^ 2))

/-- `q(z) = (z + s)·g(z)`; `s = 2` for `E⁺`, `s = −2` for `E⁻`. -/
noncomputable def qPoly (s m n : ℚ) : ℚ[X] := (X + C s) * gPoly m n

/-- `E⁺` at `(m, n)`, by the recipe: `q = (z + 2)·g`, `z0 = −2`. -/
noncomputable def EplusModel (m n : ℚ) : WeierstrassCurve ℚ := weierstrassOf (qPoly 2 m n) (-2)

/-- `E⁻` at `(m, n)`, by the recipe: `q = (z − 2)·g`, `z0 = 2`. -/
noncomputable def Eminus (m n : ℚ) : WeierstrassCurve ℚ := weierstrassOf (qPoly (-2) m n) 2

/-- `z0 = −s` is a root of `q = (z + s)·g`. -/
theorem qPoly_root (s m n : ℚ) : (qPoly s m n).eval (-s) = 0 := by
  simp [qPoly]

/-- The value of `q` at `z`. -/
theorem qPoly_eval (s m n z : ℚ) : (qPoly s m n).eval z =
    (z + s) * ((2 * (m ^ 2 + n ^ 2) * z - 8 * m * n) * (m ^ 2 * n ^ 2 * z ^ 2 + (m ^ 2 - n ^ 2) ^ 2)) := by
  simp [qPoly, gPoly]

/-- The four derivatives of `q` at `z0 = −s`, computed by Mathlib's `derivative` (product rule; no formula
for the answer is assumed). -/
theorem qPoly_d1 (s m n : ℚ) : (derivative (qPoly s m n)).eval (-s) =
    (2 * (m ^ 2 + n ^ 2) * (-s) - 8 * m * n) * (m ^ 2 * n ^ 2 * s ^ 2 + (m ^ 2 - n ^ 2) ^ 2) := by
  simp only [qPoly, gPoly, derivative_mul, derivative_add, derivative_sub, derivative_X,
    derivative_C, derivative_X_sq, eval_mul, eval_add, eval_sub, eval_C, eval_X,
    eval_pow, eval_one, eval_zero]
  ring

theorem qPoly_d2 (s m n : ℚ) : (derivative (derivative (qPoly s m n))).eval (-s) =
    2 * (2 * (m ^ 2 + n ^ 2) * (m ^ 2 * n ^ 2 * s ^ 2 + (m ^ 2 - n ^ 2) ^ 2) +
      (2 * (m ^ 2 + n ^ 2) * (-s) - 8 * m * n) * (2 * m ^ 2 * n ^ 2 * (-s))) := by
  simp only [qPoly, gPoly, derivative_mul, derivative_add, derivative_sub, derivative_X,
    derivative_C, derivative_X_sq, derivative_one, derivative_zero, eval_mul, eval_add, eval_sub, eval_C, eval_X,
    eval_pow, eval_one, eval_zero]
  ring

theorem qPoly_d3 (s m n : ℚ) : (derivative (derivative (derivative (qPoly s m n)))).eval (-s) =
    3 * (2 * (2 * (m ^ 2 + n ^ 2) * (2 * m ^ 2 * n ^ 2 * (-s))) +
      (2 * (m ^ 2 + n ^ 2) * (-s) - 8 * m * n) * (2 * m ^ 2 * n ^ 2)) := by
  simp only [qPoly, gPoly, derivative_mul, derivative_add, derivative_sub, derivative_X,
    derivative_C, derivative_X_sq, derivative_one, derivative_zero, eval_mul, eval_add, eval_sub, eval_C, eval_X,
    eval_pow, eval_one, eval_zero]
  ring

theorem qPoly_d4 (s m n : ℚ) :
    (derivative (derivative (derivative (derivative (qPoly s m n))))).eval (-s) =
      4 * 3 * 2 * (2 * (m ^ 2 + n ^ 2) * (m ^ 2 * n ^ 2)) := by
  simp only [qPoly, gPoly, derivative_mul, derivative_add, derivative_sub, derivative_X,
    derivative_C, derivative_X_sq, derivative_one, derivative_zero, eval_mul, eval_add, eval_sub, eval_C, eval_X,
    eval_pow, eval_one, eval_zero]
  ring

@[simp] theorem weierstrassOf_a₁ (q : ℚ[X]) (z0 : ℚ) : (weierstrassOf q z0).toAffine.a₁ = 0 := rfl
@[simp] theorem weierstrassOf_a₃ (q : ℚ[X]) (z0 : ℚ) : (weierstrassOf q z0).toAffine.a₃ = 0 := rfl
theorem weierstrassOf_a₂ (q : ℚ[X]) (z0 : ℚ) :
    (weierstrassOf q z0).toAffine.a₂ = (derivative (derivative q)).eval z0 / 2 := rfl
theorem weierstrassOf_a₄ (q : ℚ[X]) (z0 : ℚ) : (weierstrassOf q z0).toAffine.a₄ =
    (derivative q).eval z0 * (derivative (derivative (derivative q))).eval z0 / 6 := rfl
theorem weierstrassOf_a₆ (q : ℚ[X]) (z0 : ℚ) : (weierstrassOf q z0).toAffine.a₆ =
    (derivative q).eval z0 ^ 2 * (derivative (derivative (derivative (derivative q)))).eval z0 / 24 := rfl

/-- The recipe applied to `q = (z + s)·g` at its root `−s`, coefficient by coefficient (`q′(−s) = g(−s)`,
`q″(−s) = 2g′(−s)`, `q‴(−s) = 3g″(−s)`, `q⁗ = 4g‴ = 48·S·m²n²`), for every `s`. -/
theorem weierstrassOf_qPoly (s m n : ℚ) : weierstrassOf (qPoly s m n) (-s) =
    { a₁ := 0
      a₂ := 2 * (m ^ 2 + n ^ 2) * (m ^ 2 * n ^ 2 * s ^ 2 + (m ^ 2 - n ^ 2) ^ 2) +
        (2 * (m ^ 2 + n ^ 2) * (-s) - 8 * m * n) * (2 * m ^ 2 * n ^ 2 * (-s))
      a₃ := 0
      a₄ := (2 * (m ^ 2 + n ^ 2) * (-s) - 8 * m * n) * (m ^ 2 * n ^ 2 * s ^ 2 + (m ^ 2 - n ^ 2) ^ 2) *
        (2 * (2 * (m ^ 2 + n ^ 2) * (2 * m ^ 2 * n ^ 2 * (-s))) +
          (2 * (m ^ 2 + n ^ 2) * (-s) - 8 * m * n) * (2 * m ^ 2 * n ^ 2)) / 2
      a₆ := ((2 * (m ^ 2 + n ^ 2) * (-s) - 8 * m * n) * (m ^ 2 * n ^ 2 * s ^ 2 + (m ^ 2 - n ^ 2) ^ 2)) ^ 2 *
        (2 * (m ^ 2 + n ^ 2) * (m ^ 2 * n ^ 2)) } := by
  ext
  · rfl
  · show (derivative (derivative (qPoly s m n))).eval (-s) / 2 = _
    rw [qPoly_d2]; ring
  · rfl
  · show (derivative (qPoly s m n)).eval (-s) * (derivative (derivative (derivative (qPoly s m n)))).eval (-s) / 6 = _
    rw [qPoly_d1, qPoly_d3]; ring
  · show (derivative (qPoly s m n)).eval (-s) ^ 2 *
      (derivative (derivative (derivative (derivative (qPoly s m n))))).eval (-s) / 24 = _
    rw [qPoly_d1, qPoly_d4]; ring

/-- **`E⁺` by the recipe is `Eplus m n`**, coefficient for coefficient, for all rational `m, n`. -/
theorem EplusModel_eq (m n : ℚ) : EplusModel m n = Eplus m n := by
  rw [EplusModel, weierstrassOf_qPoly]
  ext <;> simp only [Eplus, curve, aP, bP, eP] <;> ring

/-- `Eminus m n` is the recipe at `z0 = −(−2)`. -/
theorem Eminus_def' (m n : ℚ) : Eminus m n = weierstrassOf (qPoly (-2) m n) (-(-2)) := by
  rw [neg_neg]; rfl

/-- **`E⁻` by the recipe is `Eplus m (−n)`**, coefficient for coefficient, for all rational `m, n` (this was
line A2 of `sigma_generic.txt`, a sympy check). -/
theorem Eminus_eq (m n : ℚ) : Eminus m n = Eplus m (-n) := by
  rw [Eminus_def', weierstrassOf_qPoly]
  ext <;> simp only [Eplus, curve, aP, bP, eP] <;> ring

/-- `A = q′(−s) = g(−s)`. -/
theorem mapA_qPoly (s m n : ℚ) : mapA (qPoly s m n) (-s) =
    (2 * (m ^ 2 + n ^ 2) * (-s) - 8 * m * n) * (m ^ 2 * n ^ 2 * s ^ 2 + (m ^ 2 - n ^ 2) ^ 2) :=
  qPoly_d1 s m n

/-- `A = −4(m + n)²S²` on `E⁺` and `A = 4(m − n)²S²` on `E⁻` (`UNIFORM_RANK.md` §2, A1). -/
theorem mapA_plus (m n : ℚ) : mapA (qPoly 2 m n) (-2) = -4 * (m + n) ^ 2 * (m ^ 2 + n ^ 2) ^ 2 := by
  rw [mapA_qPoly]; ring

theorem mapA_minus (m n : ℚ) : mapA (qPoly (-2) m n) 2 = 4 * (m - n) ^ 2 * (m ^ 2 + n ^ 2) ^ 2 := by
  have h := mapA_qPoly (-2) m n
  rw [neg_neg] at h
  rw [h]; ring

/-- **The map of the recipe**: if `Y² = q(z)` with `q = (z + s)·g` and `z ≠ −s`, then
`(A/(z + s), A·Y/(z + s)²)` is on the model `weierstrassOf q (−s)` (for every `s`, `m`, `n`). -/
theorem qPoly_map {s m n z Y : ℚ} (hz : z + s ≠ 0) (h : Y ^ 2 = (qPoly s m n).eval z) :
    (weierstrassOf (qPoly s m n) (-s)).toAffine.Equation (mapA (qPoly s m n) (-s) / (z + s))
      (mapA (qPoly s m n) (-s) * Y / (z + s) ^ 2) := by
  rw [WeierstrassCurve.Affine.equation_iff, weierstrassOf_a₁, weierstrassOf_a₃, weierstrassOf_a₂,
    weierstrassOf_a₄, weierstrassOf_a₆, qPoly_d1, qPoly_d2, qPoly_d3, qPoly_d4, mapA_qPoly]
  rw [qPoly_eval] at h
  obtain ⟨t, rfl⟩ : ∃ t, z = t - s := ⟨z + s, by ring⟩
  rw [show t - s + s = t by ring] at hz h ⊢
  have key : ∀ A : ℚ, (A * Y / t ^ 2) ^ 2 = A ^ 2 * Y ^ 2 / t ^ 4 := fun A => by ring
  rw [key, h]
  field_simp
  ring

/-- **The map sends `Y² = q⁺(z)`, `z ≠ −2`, into `E⁺ = Eplus m n`**: `(A/(z + 2), A·Y/(z + 2)²)`,
`A = −4(m + n)²S²`. -/
theorem EplusModel_map {m n z Y : ℚ} (hz : z ≠ -2) (h : Y ^ 2 = (qPoly 2 m n).eval z) :
    (Eplus m n).toAffine.Equation (mapA (qPoly 2 m n) (-2) / (z - (-2)))
      (mapA (qPoly 2 m n) (-2) * Y / (z - (-2)) ^ 2) := by
  rw [← EplusModel_eq, sub_neg_eq_add]
  exact qPoly_map (by intro h0; apply hz; linear_combination h0) h

/-- **The map sends `Y² = q⁻(z)`, `z ≠ 2`, into `E⁻ = Eminus m n = Eplus m (−n)`**: `(A/(z − 2), A·Y/(z − 2)²)`,
`A = 4(m − n)²S²`. -/
theorem Eminus_map {m n z Y : ℚ} (hz : z ≠ 2) (h : Y ^ 2 = (qPoly (-2) m n).eval z) :
    (Eminus m n).toAffine.Equation (mapA (qPoly (-2) m n) 2 / (z - 2))
      (mapA (qPoly (-2) m n) 2 * Y / (z - 2) ^ 2) := by
  have h1 := qPoly_map (s := -2) (z := z) (by intro h0; apply hz; linear_combination h0) h
  rw [neg_neg, ← sub_eq_add_neg] at h1
  exact h1

/-- On `E⁻` the point of `z = −2` (`t = −1`) goes to `P` of `Eplus m (−n)`: its abscissa is
`A/(−2 − 2) = −(m − n)²S² = xP m (−n)`. -/
theorem Eminus_xP (m n : ℚ) : mapA (qPoly (-2) m n) 2 / (-2 - 2) = xP m (-n) := by
  rw [mapA_minus]; unfold xP; ring

/-- On `E⁺` the point of `z = 2` (`t = 1`) goes to `P`: `A/(2 + 2) = −(m + n)²S² = xP m n`. -/
theorem EplusModel_xP (m n : ℚ) : mapA (qPoly 2 m n) (-2) / (2 - (-2)) = xP m n := by
  rw [mapA_plus]; unfold xP; ring

/-! ## From the curve `y² = F0·F8·GL` to the models

`g` is the main library's `G1.sigmaG` (`G1Sigma.lean`), and `G1.sigma_plus`, `G1.sigma_minus` say that at
`z = t + 1/t`, `(z ± 2)·g(z) = ((t ± 1)/t²)²·F0·F8·GL(t)`. So a point `(t, y)` of `y² = F0·F8·GL` with `t ≠ 0`
gives the point `(z, (t ± 1)y/t²)` of `Y² = q^±(z)`, and the recipe's map takes it to `E^±`. -/

theorem gPoly_eval (m n z : ℚ) : (gPoly m n).eval z = MagicSquareOfSquares.G1.sigmaG m n z := by
  simp [gPoly, MagicSquareOfSquares.G1.sigmaG]

/-- **From `y² = F0·F8·GL` to `E⁺ = Eplus m n`**, for `t ≠ 0, −1` (`t = −1` is `z = z0 = −2`). -/
theorem sextic_to_Eplus {m n t y : ℚ} (ht : t ≠ 0) (ht1 : t ≠ -1)
    (h : y ^ 2 = (m ^ 2 + n ^ 2 * t ^ 2) * (n ^ 2 + m ^ 2 * t ^ 2) *
      (2 * (m ^ 2 + n ^ 2) * (t ^ 2 + 1) - 8 * m * n * t)) :
    (Eplus m n).toAffine.Equation (mapA (qPoly 2 m n) (-2) / ((t + 1 / t) - (-2)))
      (mapA (qPoly 2 m n) (-2) * ((t + 1) / t ^ 2 * y) / ((t + 1 / t) - (-2)) ^ 2) := by
  refine EplusModel_map ?_ ?_
  · intro h0
    have e : (t + 1) ^ 2 = t * ((t + 1 / t) - (-2)) := by
      rw [mul_sub, mul_add, mul_one_div_cancel ht]; ring
    rw [h0, sub_self, mul_zero] at e
    exact ht1 (by linear_combination (pow_eq_zero_iff two_ne_zero).mp e)
  · have hs := MagicSquareOfSquares.G1.sigma_plus m n t ht
    rw [qPoly_eval]
    unfold MagicSquareOfSquares.G1.sigmaG at hs
    linear_combination ((t + 1) / t ^ 2) ^ 2 * h - hs

/-- **From `y² = F0·F8·GL` to `E⁻ = Eminus m n = Eplus m (−n)`**, for `t ≠ 0, 1` (`t = 1` is `z = z0 = 2`). -/
theorem sextic_to_Eminus {m n t y : ℚ} (ht : t ≠ 0) (ht1 : t ≠ 1)
    (h : y ^ 2 = (m ^ 2 + n ^ 2 * t ^ 2) * (n ^ 2 + m ^ 2 * t ^ 2) *
      (2 * (m ^ 2 + n ^ 2) * (t ^ 2 + 1) - 8 * m * n * t)) :
    (Eminus m n).toAffine.Equation (mapA (qPoly (-2) m n) 2 / ((t + 1 / t) - 2))
      (mapA (qPoly (-2) m n) 2 * ((t - 1) / t ^ 2 * y) / ((t + 1 / t) - 2) ^ 2) := by
  refine Eminus_map ?_ ?_
  · intro h0
    have e : (t - 1) ^ 2 = t * ((t + 1 / t) - 2) := by
      rw [mul_sub, mul_add, mul_one_div_cancel ht]; ring
    rw [h0, sub_self, mul_zero] at e
    exact ht1 (by linear_combination (pow_eq_zero_iff two_ne_zero).mp e)
  · have hs := MagicSquareOfSquares.G1.sigma_minus m n t ht
    rw [qPoly_eval]
    unfold MagicSquareOfSquares.G1.sigmaG at hs
    linear_combination ((t - 1) / t ^ 2) ^ 2 * h - hs

/-! ## Theorems 1 and 2 of `UniformRank.lean` with the names `EplusModel`, `Eminus` -/

/-- **Theorem 1** (rank `≥ 1`) on `E⁺ = EplusModel m n` and `E⁻ = Eminus m n`, the curves of the recipe of
`UNIFORM_RANK.md` §1, for coprime `m, n` with `mn(m² − n²) ≠ 0`. -/
theorem theorem1_Eminus {m n : ℤ} (hmn : IsCoprime m n) (hm : m ≠ 0) (hn : n ≠ 0) (hne : m ^ 2 ≠ n ^ 2) :
    1 ≤ Module.finrank ℤ (EplusModel (m : ℚ) (n : ℚ)).toAffine.Point ∧
      1 ≤ Module.finrank ℤ (Eminus (m : ℚ) (n : ℚ)).toAffine.Point := by
  rw [EplusModel_eq, Eminus_eq]
  exact UniformRank.theorem1 hmn hm hn hne

/-- **Theorem 2** (rank `≥ 2`) on `E⁺ = EplusModel m n` and `E⁻ = Eminus m n`, for coprime `m, n` with
`m² + n² = 2k²`, `m² ≠ n²`. (The independence of `R, Q, T` is `UniformRank.theorem2`, stated on
`Eplus m (±n)`, which `EplusModel_eq` and `Eminus_eq` identify with these curves.) -/
theorem theorem2_Eminus {m n k : ℤ} (hmn : IsCoprime m n) (hk : m ^ 2 + n ^ 2 = 2 * k ^ 2)
    (hne : m ^ 2 ≠ n ^ 2) :
    2 ≤ Module.finrank ℤ (EplusModel (m : ℚ) (n : ℚ)).toAffine.Point ∧
      2 ≤ Module.finrank ℤ (Eminus (m : ℚ) (n : ℚ)).toAffine.Point := by
  rw [EplusModel_eq, Eminus_eq]
  obtain ⟨-, -, h1, h2⟩ := UniformRank.theorem2 hmn hk hne
  exact ⟨h1, h2⟩

/-! ## Controls -/

/-- The recipe is not vacuous: at `(7, 1)` the model of `E⁻` differs from that of `E⁺`. -/
theorem control_Eminus_ne_Eplus : Eminus 7 1 ≠ Eplus 7 1 := by
  rw [Eminus_eq]
  intro h
  have := congrArg WeierstrassCurve.a₂ h
  simp only [Eplus, curve] at this
  unfold aP eP at this
  norm_num at this

/-- Using the wrong root (`z0 = 2` for `q⁺`) gives another curve. -/
theorem control_wrong_root : weierstrassOf (qPoly 2 7 1) 2 ≠ Eplus 7 1 := by
  intro h
  have := congrArg WeierstrassCurve.a₂ h
  simp only [weierstrassOf, Eplus, curve] at this
  unfold aP eP at this
  simp only [qPoly, gPoly, derivative_mul, derivative_add, derivative_sub, derivative_X,
    derivative_C, derivative_X_sq, derivative_one, derivative_zero, eval_mul, eval_add, eval_sub, eval_C, eval_X,
    eval_pow, eval_one, eval_zero] at this
  norm_num at this

end MagicSquareOfSquares.EminusModel

#print axioms MagicSquareOfSquares.EminusModel.EplusModel_eq
#print axioms MagicSquareOfSquares.EminusModel.Eminus_eq
#print axioms MagicSquareOfSquares.EminusModel.qPoly_map
#print axioms MagicSquareOfSquares.EminusModel.EplusModel_map
#print axioms MagicSquareOfSquares.EminusModel.Eminus_map
#print axioms MagicSquareOfSquares.EminusModel.sextic_to_Eplus
#print axioms MagicSquareOfSquares.EminusModel.sextic_to_Eminus
#print axioms MagicSquareOfSquares.EminusModel.theorem1_Eminus
#print axioms MagicSquareOfSquares.EminusModel.theorem2_Eminus
#print axioms MagicSquareOfSquares.EminusModel.control_Eminus_ne_Eplus
