import LeanRank.LemmaDBridge
import LeanRank.TwoIsogenyMaps

/-!
# Rank 0 from the two descent statements (2-isogeny descent, glue part)

Written by an AI model (Claude, made by Anthropic) on 26 November 2026 (UTC); not read by a person.
Experiment, not counted in the census.

For `W : y² = x(x² + ax + b)` over a field `F` of characteristic 0 (`b ≠ 0`, `a² − 4b ≠ 0`) and the
isogenous `W' : Y² = X(X² − 2aX + (a² − 4b))`, with the explicit maps
`φ(x, y) = (y²/x², y(b − x²)/x²)` and `ψ(X, Y) = (Y²/(4X²), Y(a² − 4b − X²)/(8X²))` on Mathlib's points
(no homomorphism property is used):
* `psi_phi`: `ψ(φ(R)) = R + R`;
* `psi_add_T`: `ψ(Q + T') = ψ(Q)`;
* `exists_psi_eq`: a point with `x = u²` is in the image of `ψ`; `exists_phi_eq`: the same for `φ`;
* `add_T_some`: `x(P + T) · x(P) = b`.
`cover`: if every affine point of `W` with `x ≠ 0` has `x` or `b/x` a square (`DescentW`) and every affine
point of `W'` with `X ≠ 0` has `X` or `(a² − 4b)/X` a square (`DescentW'`), then every point is `2R` or
`2R + T`. `finrank_eq_zero_of_cover`: in a finitely generated abelian group with an element `T₀ ≠ 0` of
order 2, this forces `finrank ℤ G = 0` (Stoll's `AddSubgroup.index_range_nsmul_of_fg`).
`E1K_finrank_eq_zero_of_descent`: the case of `E1K` over `ℚ(√2)`.
-/

namespace MagicSquareOfSquares.E1Rank

open WeierstrassCurve WeierstrassCurve.Affine MagicSquareOfSquares.Family Module

section Generic

variable {F : Type*} [Field F] [CharZero F]

/-- `y² = x(x² + ax + b)`. -/
def Wc (a b : F) : WeierstrassCurve.Affine F := ⟨0, a, 0, b, 0⟩

variable {a b : F}

@[simp] lemma Wc_a₁ : (Wc a b).a₁ = 0 := rfl
@[simp] lemma Wc_a₂ : (Wc a b).a₂ = a := rfl
@[simp] lemma Wc_a₃ : (Wc a b).a₃ = 0 := rfl
@[simp] lemma Wc_a₄ : (Wc a b).a₄ = b := rfl
@[simp] lemma Wc_a₆ : (Wc a b).a₆ = 0 := rfl

lemma Wc_eq {x y : F} : (Wc a b).Equation x y ↔ y ^ 2 = x * (x ^ 2 + a * x + b) := by
  rw [Affine.equation_iff]
  simp only [Wc_a₁, Wc_a₂, Wc_a₃, Wc_a₄, Wc_a₆]
  constructor <;> intro h <;> linear_combination h

lemma Wc_Δ : (Wc a b).Δ = 16 * b ^ 2 * (a ^ 2 - 4 * b) := by
  simp only [Wc, WeierstrassCurve.Δ, WeierstrassCurve.b₂, WeierstrassCurve.b₄, WeierstrassCurve.b₆,
    WeierstrassCurve.b₈]
  ring

lemma Wc_ns (hb : b ≠ 0) (hd : a ^ 2 - 4 * b ≠ 0) {x y : F} (h : y ^ 2 = x * (x ^ 2 + a * x + b)) :
    (Wc a b).Nonsingular x y := by
  have hΔ : (Wc a b).Δ ≠ 0 := by
    rw [Wc_Δ]; exact mul_ne_zero (mul_ne_zero (by norm_num) (pow_ne_zero _ hb)) hd
  exact (equation_iff_nonsingular_of_Δ_ne_zero hΔ).mp (Wc_eq.mpr h)

lemma Wc_negY (x y : F) : (Wc a b).negY x y = -y := by
  simp only [negY, Wc_a₁, Wc_a₃]; ring

lemma Wc_addX (x₁ x₂ ℓ : F) : (Wc a b).addX x₁ x₂ ℓ = ℓ ^ 2 - a - x₁ - x₂ := by
  simp only [addX, Wc_a₁, Wc_a₂]; ring

lemma Wc_addY (x₁ x₂ y₁ ℓ : F) :
    (Wc a b).addY x₁ x₂ y₁ ℓ = -(ℓ * (ℓ ^ 2 - a - x₁ - x₂ - x₁) + y₁) := by
  simp only [addY, negAddY, Wc_negY, Wc_addX]

lemma y_eq_zero_of_x_eq_zero {y : F} (h : (Wc a b).Equation 0 y) : y = 0 := by
  have h' := Wc_eq.mp h
  exact pow_eq_zero_iff (n := 2) two_ne_zero |>.mp (by rw [h']; ring)

theorem some_eq {W : WeierstrassCurve.Affine F} {x y x' y' : F} {h : W.Nonsingular x y}
    {h' : W.Nonsingular x' y'} (hx : x = x') (hy : y = y') :
    Point.some x y h = Point.some x' y' h' := by
  subst hx; subst hy; rfl

/-- The hypotheses `b ≠ 0`, `a² − 4b ≠ 0`. -/
structure Good (a b : F) : Prop where
  hb : b ≠ 0
  hd : a ^ 2 - 4 * b ≠ 0

lemma Good.ns (G : Good a b) {x y : F} (h : y ^ 2 = x * (x ^ 2 + a * x + b)) :
    (Wc a b).Nonsingular x y := Wc_ns G.hb G.hd h

lemma Good.ns' (G : Good a b) {X Y : F}
    (h : Y ^ 2 = X * (X ^ 2 + (-2 * a) * X + (a ^ 2 - 4 * b))) :
    (Wc (-2 * a) (a ^ 2 - 4 * b)).Nonsingular X Y :=
  Wc_ns G.hd (by intro h0; apply G.hb; linear_combination h0 / 16) h

lemma phi_ns (G : Good a b) {x y : F} (h : (Wc a b).Nonsingular x y) (hx : x ≠ 0) :
    (Wc (-2 * a) (a ^ 2 - 4 * b)).Nonsingular (y ^ 2 / x ^ 2) (y * (b - x ^ 2) / x ^ 2) :=
  G.ns' (by rw [isog_on hx (Wc_eq.mp h.1)]; ring)

lemma psi_ns (G : Good a b) {X Y : F} (h : (Wc (-2 * a) (a ^ 2 - 4 * b)).Nonsingular X Y)
    (hX : X ≠ 0) :
    (Wc a b).Nonsingular (Y ^ 2 / (4 * X ^ 2)) (Y * (a ^ 2 - 4 * b - X ^ 2) / (8 * X ^ 2)) :=
  G.ns (isog_dual_on hX (by rw [Wc_eq.mp h.1]; ring))

variable [DecidableEq F]

/-- The 2-isogeny `φ : W → W'`. -/
noncomputable def phi (G : Good a b) : (Wc a b).Point → (Wc (-2 * a) (a ^ 2 - 4 * b)).Point
  | .zero => 0
  | .some x y h => if hx : x = 0 then 0 else .some _ _ (phi_ns G h hx)

/-- The dual 2-isogeny `ψ : W' → W`. -/
noncomputable def psi (G : Good a b) : (Wc (-2 * a) (a ^ 2 - 4 * b)).Point → (Wc a b).Point
  | .zero => 0
  | .some X Y h => if hX : X = 0 then 0 else .some _ _ (psi_ns G h hX)

lemma phi_some (G : Good a b) {x y : F} (h : (Wc a b).Nonsingular x y) (hx : x ≠ 0) :
    phi G (.some x y h) = .some _ _ (phi_ns G h hx) := by
  simp only [phi, dif_neg hx]

lemma phi_some_zero (G : Good a b) {y : F} (h : (Wc a b).Nonsingular 0 y) :
    phi G (.some 0 y h) = 0 := by
  simp [phi]

lemma psi_some (G : Good a b) {X Y : F} (h : (Wc (-2 * a) (a ^ 2 - 4 * b)).Nonsingular X Y)
    (hX : X ≠ 0) : psi G (.some X Y h) = .some _ _ (psi_ns G h hX) := by
  simp only [psi, dif_neg hX]

lemma psi_some_zero (G : Good a b) {X Y : F} (h : (Wc (-2 * a) (a ^ 2 - 4 * b)).Nonsingular X Y)
    (hX : X = 0) : psi G (.some X Y h) = 0 := by
  simp only [psi, dif_pos hX]

/-- `T = (0, 0)` on `W`. -/
noncomputable def T0 (G : Good a b) : (Wc a b).Point := .some 0 0 (G.ns (by ring))

/-- `T' = (0, 0)` on `W'`. -/
noncomputable def T1 (G : Good a b) : (Wc (-2 * a) (a ^ 2 - 4 * b)).Point :=
  .some 0 0 (G.ns' (by ring))

theorem T0_ne_zero (G : Good a b) : T0 G ≠ 0 := Point.some_ne_zero _

theorem T0_add_T0 (G : Good a b) : T0 G + T0 G = 0 := by
  unfold T0
  exact Point.add_self_of_Y_eq (by rw [Wc_negY, neg_zero])

/-- **`ψ(φ(R)) = 2R`.** -/
theorem psi_phi (G : Good a b) (R : (Wc a b).Point) : psi G (phi G R) = R + R := by
  rcases R with _ | ⟨x, y, h⟩
  · show psi G (phi G 0) = 0 + 0
    rw [add_zero]; rfl
  · by_cases hx : x = 0
    · subst hx
      have hy : y = 0 := y_eq_zero_of_x_eq_zero h.1
      subst hy
      rw [Point.add_self_of_Y_eq (h₁ := h) (by rw [Wc_negY, neg_zero]), phi_some_zero]
      rfl
    · have hx' : x ≠ 0 := hx
      have heq := Wc_eq.mp h.1
      by_cases hy : y = 0
      · subst hy
        rw [Point.add_self_of_Y_eq (h₁ := h) (by rw [Wc_negY, neg_zero]), phi_some G h hx,
          psi_some_zero G _ (show (0 : F) ^ 2 / x ^ 2 = 0 by ring)]
      · have hy0 : y ≠ 0 := hy
        have hX : y ^ 2 / x ^ 2 ≠ 0 := div_ne_zero (pow_ne_zero _ hy) (pow_ne_zero _ hx)
        have hy' : y ≠ (Wc a b).negY x y := by
          rw [Wc_negY]; intro h'; apply hy; linear_combination h' / 2
        rw [phi_some G h hx, psi_some G _ hX, Point.add_self_of_Y_ne hy']
        have hℓ : (Wc a b).slope x x y y = (3 * x ^ 2 + 2 * a * x + b) / (2 * y) := by
          rw [slope_of_Y_ne rfl hy', Wc_negY]
          simp only [Wc_a₁, Wc_a₂, Wc_a₄]
          congr 1 <;> ring
        apply some_eq
        · rw [Wc_addX, hℓ, ← sub_eq_zero]
          have e : (y * (b - x ^ 2) / x ^ 2) ^ 2 / (4 * (y ^ 2 / x ^ 2) ^ 2) -
              (((3 * x ^ 2 + 2 * a * x + b) / (2 * y)) ^ 2 - a - x - x)
              = (y ^ 2 - x * (x ^ 2 + a * x + b)) * (a + 2 * x) / y ^ 2 := by
            field_simp
            ring
          rw [e, heq, sub_self, zero_mul, zero_div]
        · rw [Wc_addY, hℓ, ← sub_eq_zero]
          have e : y * (b - x ^ 2) / x ^ 2 * (a ^ 2 - 4 * b - (y ^ 2 / x ^ 2) ^ 2) /
                (8 * (y ^ 2 / x ^ 2) ^ 2) -
              -((3 * x ^ 2 + 2 * a * x + b) / (2 * y) *
                (((3 * x ^ 2 + 2 * a * x + b) / (2 * y)) ^ 2 - a - x - x - x) + y)
              = (y ^ 2 - x * (x ^ 2 + a * x + b)) *
                (-(8 * a ^ 2 * x ^ 3 + 5 * a * b * x ^ 2 + 27 * a * x ^ 4 + b ^ 2 * x + 4 * b * x ^ 3
                  + b * y ^ 2 + 27 * x ^ 5 - 9 * x ^ 2 * y ^ 2) / (8 * x ^ 2 * y ^ 3)) := by
            field_simp
            ring
          rw [e, heq, sub_self, zero_mul]

/-- **`x(P + T) · x(P) = β`** for an affine point `P` with `x ≠ 0` on `y² = x(x² + αx + β)`. -/
theorem add_T_some {α β : F} (hT : (Wc α β).Nonsingular 0 0) {x y : F}
    (h : (Wc α β).Nonsingular x y) (hx : x ≠ 0) :
    ∃ x' y' h', Point.some x y h + Point.some 0 0 hT = Point.some x' y' h' ∧ x' * x = β := by
  refine ⟨_, _, _, Point.add_of_X_ne hx, ?_⟩
  rw [slope_of_X_ne hx, Wc_addX]
  have heq := Wc_eq.mp h.1
  have hx' : x ≠ 0 := hx
  have e : (((y - 0) / (x - 0)) ^ 2 - α - x - 0) * x - β
      = (y ^ 2 - x * (x ^ 2 + α * x + β)) / x := by
    field_simp
    ring
  rw [← sub_eq_zero, e, heq, sub_self, zero_div]

/-- **`ψ(Q + T') = ψ(Q)`.** -/
theorem psi_add_T (G : Good a b) (Q : (Wc (-2 * a) (a ^ 2 - 4 * b)).Point) :
    psi G (Q + T1 G) = psi G Q := by
  rcases Q with _ | ⟨X, Y, h⟩
  · show psi G (0 + T1 G) = psi G 0
    rw [zero_add, T1, psi_some_zero G _ rfl]; rfl
  · by_cases hX : X = 0
    · subst hX
      have hY : Y = 0 := y_eq_zero_of_x_eq_zero h.1
      subst hY
      rw [T1, Point.add_self_of_Y_eq (h₁ := h) (by rw [Wc_negY, neg_zero]), psi_some_zero G _ rfl]
      rfl
    · have heq := Wc_eq.mp h.1
      have hbp : a ^ 2 - 4 * b ≠ 0 := G.hd
      have hX' : X ≠ 0 := hX
      have hℓ : (Wc (-2 * a) (a ^ 2 - 4 * b)).slope X 0 Y 0 = Y / X := by
        rw [slope_of_X_ne hX]; simp
      have hX3 : (Y / X) ^ 2 - -2 * a - X - 0 = (a ^ 2 - 4 * b) / X := by
        rw [← sub_eq_zero]
        have e : (Y / X) ^ 2 - -2 * a - X - 0 - (a ^ 2 - 4 * b) / X
            = (Y ^ 2 - X * (X ^ 2 + -2 * a * X + (a ^ 2 - 4 * b))) / X ^ 2 := by
          field_simp
          ring
        rw [e, heq, sub_self, zero_div]
      have hX3' : (Wc (-2 * a) (a ^ 2 - 4 * b)).addX X 0 ((Wc (-2 * a) (a ^ 2 - 4 * b)).slope X 0 Y 0)
          = (a ^ 2 - 4 * b) / X := by
        rw [Wc_addX, hℓ, hX3]
      have hX3ne : (Wc (-2 * a) (a ^ 2 - 4 * b)).addX X 0
          ((Wc (-2 * a) (a ^ 2 - 4 * b)).slope X 0 Y 0) ≠ 0 := by
        rw [hX3']; exact div_ne_zero hbp hX
      rw [T1, Point.add_of_X_ne hX, psi_some G h hX, psi_some G _ hX3ne]
      apply some_eq
      · rw [hX3', Wc_addY, hℓ, hX3]
        field_simp
        ring
      · rw [hX3', Wc_addY, hℓ, hX3]
        field_simp
        ring

/-- A point with `x = u²` is `ψ` of a point of `W'`. -/
theorem exists_psi_eq (G : Good a b) {x y : F} (h : (Wc a b).Nonsingular x y) {u : F} (hu : u ≠ 0)
    (hxu : x = u ^ 2) : ∃ Q, psi G Q = .some x y h := by
  have heq := Wc_eq.mp h.1
  obtain ⟨X, hXdef⟩ : ∃ X : F, X = a + 2 * u ^ 2 + 2 * y / u := ⟨_, rfl⟩
  have hy' : y = u * (X - a - 2 * u ^ 2) / 2 := by
    rw [hXdef]; field_simp; ring
  have R6 : X ^ 2 - 2 * a * X - 4 * u ^ 2 * X + (a ^ 2 - 4 * b) = 0 := by
    have h4 : u ^ 2 * (X ^ 2 - 2 * a * X - 4 * u ^ 2 * X + (a ^ 2 - 4 * b)) = 0 := by
      rw [hxu, hy'] at heq
      linear_combination 4 * heq
    rcases mul_eq_zero.mp h4 with h0 | h0
    · exact absurd (pow_eq_zero_iff two_ne_zero |>.mp h0) hu
    · exact h0
  have hX : X ≠ 0 := by
    intro h0; apply G.hd; rw [h0] at R6; linear_combination R6
  have hQ : (-2 * u * X) ^ 2 = X * (X ^ 2 + (-2 * a) * X + (a ^ 2 - 4 * b)) := by
    linear_combination (-X) * R6
  refine ⟨.some X (-2 * u * X) (G.ns' hQ), ?_⟩
  rw [psi_some G _ hX]
  apply some_eq
  · rw [hxu]; field_simp; ring
  · rw [hy', div_eq_div_iff (mul_ne_zero (by norm_num) (pow_ne_zero 2 hX)) (by norm_num)]
    linear_combination (-4 * u * X) * R6

/-- A point of `W'` with `X = v²` is `φ` of a point of `W`. -/
theorem exists_phi_eq (G : Good a b) {X Y : F} (h : (Wc (-2 * a) (a ^ 2 - 4 * b)).Nonsingular X Y)
    {v : F} (hv : v ≠ 0) (hXv : X = v ^ 2) : ∃ R, phi G R = .some X Y h := by
  have heq := Wc_eq.mp h.1
  obtain ⟨x, hxdef⟩ : ∃ x : F, x = (v ^ 2 - a + Y / v) / 2 := ⟨_, rfl⟩
  have hY' : Y = v * (2 * x - v ^ 2 + a) := by rw [hxdef]; field_simp; ring
  have R7 : x ^ 2 + a * x - v ^ 2 * x + b = 0 := by
    have h4 : v ^ 2 * (x ^ 2 + a * x - v ^ 2 * x + b) = 0 := by
      rw [hXv, hY'] at heq
      linear_combination heq / 4
    rcases mul_eq_zero.mp h4 with h0 | h0
    · exact absurd (pow_eq_zero_iff two_ne_zero |>.mp h0) hv
    · exact h0
  have hx : x ≠ 0 := by intro h0; apply G.hb; rw [h0] at R7; linear_combination R7
  have hR : (-v * x) ^ 2 = x * (x ^ 2 + a * x + b) := by linear_combination (-x) * R7
  refine ⟨.some x (-v * x) (G.ns hR), ?_⟩
  rw [phi_some G _ hx]
  apply some_eq
  · rw [hXv]; field_simp
  · rw [hY', div_eq_iff (pow_ne_zero 2 hx)]
    linear_combination (-v * x) * R7

/-- Descent statement on `W`: `x` or `b/x` is a square. -/
def DescentW (a b : F) : Prop := ∀ x y : F, y ^ 2 = x * (x ^ 2 + a * x + b) → x ≠ 0 →
  (∃ u, x = u ^ 2) ∨ (∃ u, b = x * u ^ 2)

/-- Descent statement on `W'`: `X` or `(a² − 4b)/X` is a square. -/
def DescentW' (a b : F) : Prop := ∀ X Y : F,
  Y ^ 2 = X * (X ^ 2 + (-2 * a) * X + (a ^ 2 - 4 * b)) → X ≠ 0 →
  (∃ v, X = v ^ 2) ∨ (∃ v, a ^ 2 - 4 * b = X * v ^ 2)

/-- Every `ψ(Q)` is a double. -/
theorem psi_double (G : Good a b) (hD' : DescentW' a b) (Q : (Wc (-2 * a) (a ^ 2 - 4 * b)).Point) :
    ∃ R, psi G Q = R + R := by
  rcases Q with _ | ⟨X, Y, h⟩
  · exact ⟨0, by rw [add_zero]; rfl⟩
  · by_cases hX : X = 0
    · exact ⟨0, by rw [psi_some_zero G h hX, add_zero]⟩
    · rcases hD' X Y (Wc_eq.mp h.1) hX with ⟨v, hv⟩ | ⟨v, hv⟩
      · have hv0 : v ≠ 0 := by rintro rfl; apply hX; rw [hv]; ring
        obtain ⟨R, hR⟩ := exists_phi_eq G h hv0 hv
        exact ⟨R, by rw [← hR, psi_phi]⟩
      · obtain ⟨X', Y', h', hQ, hX'⟩ := add_T_some (G.ns' (by ring)) h hX
        have hv0 : v ≠ 0 := by rintro rfl; apply G.hd; rw [hv]; ring
        have hX'v : X' = v ^ 2 := by
          apply mul_right_cancel₀ hX; rw [hX', hv]; ring
        obtain ⟨R, hR⟩ := exists_phi_eq G h' hv0 hX'v
        refine ⟨R, ?_⟩
        rw [← psi_add_T G, show Point.some X Y h + T1 G = Point.some X' Y' h' from hQ, ← hR, psi_phi]

/-- **Every point of `W` is `2R` or `2R + T`**, given the two descent statements. -/
theorem cover (G : Good a b) (hD : DescentW a b) (hD' : DescentW' a b) (P : (Wc a b).Point) :
    ∃ R, P = R + R ∨ P = R + R + T0 G := by
  rcases P with _ | ⟨x, y, h⟩
  · exact ⟨0, Or.inl (by rw [add_zero]; rfl)⟩
  · by_cases hx : x = 0
    · subst hx
      have hy : y = 0 := y_eq_zero_of_x_eq_zero h.1
      subst hy
      exact ⟨0, Or.inr (by rw [add_zero, zero_add]; rfl)⟩
    · rcases hD x y (Wc_eq.mp h.1) hx with ⟨u, hu⟩ | ⟨u, hu⟩
      · have hu0 : u ≠ 0 := by rintro rfl; apply hx; rw [hu]; ring
        obtain ⟨Q, hQ⟩ := exists_psi_eq G h hu0 hu
        obtain ⟨R, hR⟩ := psi_double G hD' Q
        exact ⟨R, Or.inl (by rw [← hQ, hR])⟩
      · obtain ⟨x', y', h', hP, hx'⟩ := add_T_some (G.ns (by ring)) h hx
        have hu0 : u ≠ 0 := by rintro rfl; apply G.hb; rw [hu]; ring
        have hx'u : x' = u ^ 2 := by apply mul_right_cancel₀ hx; rw [hx', hu]; ring
        obtain ⟨Q, hQ⟩ := exists_psi_eq G h' hu0 hx'u
        obtain ⟨R, hR⟩ := psi_double G hD' Q
        refine ⟨R, Or.inr ?_⟩
        have hP' : Point.some x y h + T0 G = R + R := by
          rw [show Point.some x y h + T0 G = Point.some x' y' h' from hP, ← hQ, hR]
        rw [← hP', add_assoc, T0_add_T0, add_zero]

end Generic

/-- **Group theory.** A finitely generated abelian group with an element `T₀ ≠ 0` of order 2 in which
every element is `2R` or `2R + T₀` has `finrank ℤ G = 0`. -/
theorem finrank_eq_zero_of_cover {G : Type*} [AddCommGroup G] [IsAddFG G] {T₀ : G} (hT₀ : T₀ ≠ 0)
    (hTT : T₀ + T₀ = 0) (hcov : ∀ P : G, ∃ R : G, P = R + R ∨ P = R + R + T₀) :
    finrank ℤ G = 0 := by
  set H := (nsmulAddMonoidHom (α := G) 2).range with hH
  let f : Bool → G ⧸ H := fun c => if c then ((T₀ : G) : G ⧸ H) else ((0 : G) : G ⧸ H)
  have hf : Function.Surjective f := by
    intro q
    induction q using QuotientAddGroup.induction_on with
    | H P =>
      obtain ⟨R, hR | hR⟩ := hcov P
      · refine ⟨false, ?_⟩
        show ((0 : G) : G ⧸ H) = (P : G ⧸ H)
        rw [QuotientAddGroup.eq]
        exact AddMonoidHom.mem_range.mpr ⟨R, by show 2 • R = _; rw [two_nsmul, hR]; abel⟩
      · refine ⟨true, ?_⟩
        show ((T₀ : G) : G ⧸ H) = (P : G ⧸ H)
        rw [QuotientAddGroup.eq]
        exact AddMonoidHom.mem_range.mpr ⟨R, by show 2 • R = _; rw [two_nsmul, hR]; abel⟩
  have hfin : Finite (G ⧸ H) := Finite.of_surjective f hf
  have hidx_le : H.index ≤ 2 := by
    have := Nat.card_le_card_of_surjective f hf
    rw [Nat.card_eq_fintype_card (α := Bool), Fintype.card_bool] at this
    exact this
  have hidx0 : H.index ≠ 0 := AddSubgroup.index_ne_zero_of_finite
  have hform : H.index = 2 ^ finrank ℤ G * Nat.card (nsmulAddMonoidHom (α := G) 2).ker :=
    AddSubgroup.index_range_nsmul_of_fg G two_ne_zero
  set K := (nsmulAddMonoidHom (α := G) 2).ker with hK
  have hK0 : Nat.card K ≠ 0 := by
    intro h0; rw [h0, mul_zero] at hform; exact hidx0 hform
  haveI hKfin : Finite K := Nat.finite_of_card_ne_zero hK0
  have hTmem : T₀ ∈ K := by
    show 2 • T₀ = 0
    rw [two_nsmul, hTT]
  haveI : Nontrivial K := ⟨⟨0, ⟨T₀, hTmem⟩, fun h => hT₀ (congrArg Subtype.val h).symm⟩⟩
  have hK2 : 2 ≤ Nat.card K := Finite.one_lt_card
  by_contra hr
  have h2r : 2 ≤ 2 ^ finrank ℤ G := Nat.le_self_pow hr 2
  have := Nat.mul_le_mul h2r hK2
  omega

section E1

theorem good_E1 : Good (8 * (2 + s2) : Q2) (-8 * (2 + s2) ^ 2) := by
  have h := two_add_s2_ne_zero
  constructor
  · intro h0
    apply h
    have h1 : (2 + s2) ^ 2 = 0 := by linear_combination (-1 / 8 : Q2) * h0
    exact pow_eq_zero_iff two_ne_zero |>.mp h1
  · intro h0
    apply h
    have h1 : (2 + s2) ^ 2 = 0 := by linear_combination (1 / 96 : Q2) * h0
    exact pow_eq_zero_iff two_ne_zero |>.mp h1

/-- **Rank 0 of `E1(ℚ(√2))` from the two descent statements** (`DescentW`, `DescentW'` for
`a = 8(2 + √2)`, `b = −8(2 + √2)²`). -/
theorem E1K_finrank_eq_zero_of_descent
    (hD : DescentW (8 * (2 + s2) : Q2) (-8 * (2 + s2) ^ 2))
    (hD' : DescentW' (8 * (2 + s2) : Q2) (-8 * (2 + s2) ^ 2)) :
    finrank ℤ E1K.Point = 0 := by
  exact @finrank_eq_zero_of_cover _ _ E1K_fg _ (T0_ne_zero good_E1) (T0_add_T0 good_E1)
    (cover good_E1 hD hD')

end E1

end MagicSquareOfSquares.E1Rank

#print axioms MagicSquareOfSquares.E1Rank.psi_phi
#print axioms MagicSquareOfSquares.E1Rank.psi_add_T
#print axioms MagicSquareOfSquares.E1Rank.cover
#print axioms MagicSquareOfSquares.E1Rank.finrank_eq_zero_of_cover
#print axioms MagicSquareOfSquares.E1Rank.E1K_finrank_eq_zero_of_descent
