import MagicSquareOfSquares.CoverDescent

/-!
# The class of a solution, explicitly; and local points

Toward Theorem C of `COMPLETE.md`, whose proof is by hand. What is here, in Lean:

* `class_iff` (Lemma 1 of `COMPLETE.md`, over the integers): at a solution on the four forms with
  roots `w0, w8, wL, wU` of the four norms, a `d` has the property of the eight products **if and
  only if** `d · D` or `−d · D` is a square, where `D = clsD = ∏ (w_j + x_j)`. So the class of a
  solution is given by a formula.
* `gsqM_class`: the same identity modulo any `M`: roots of the norms modulo `M` and
  `d ≡ σ² · 16 D (mod M)` give `d · Q_f ≡ (a + ib)² (mod M)` for the eight products.
* `exists_sq_modEq_pow`: a unit that is a square modulo an odd prime `p` is a square modulo
  every power of `p` (Hensel, by induction).
* `pass_first_order` (section 4.1 of `COMPLETE.md`): for an odd prime `p` dividing `m`, at a point
  `(u, v)` of first order every `d` in the class of `2uv` modulo `p` passes modulo **every** power
  of `p`. `two_passes`: in particular `d = 2` at `(1, 1)`.

* `pass_odd_valuation` (section 4.3, the case `p² | m`): for an odd prime `p` modulo which 2 is a
  square, the classes `p · d₁` pass at `(p·u₁, 1)` modulo every power of `p`.
* With no hypothesis on a point: `unit_passes_of_three_mod_four` (every unit class at a prime
  `p ≡ 3 (mod 4)` dividing `m`) and `prime_mul_unit_passes_mod_eight` (every class `p · unit` at
  a prime `p ≡ ±1 (mod 8)` with `p² | m`).

"Passes modulo every power of `p`" is what a `p`-adic point with that class amounts to; the step
from all powers to a `p`-adic point is not in Lean. In later files: the conic of 4.2, by an identity
of Jacobsthal in place of the bound of Hasse (`ClassConic.lean`); the case `p ∥ m` of 4.3
(`ClassOdd.lean`); the prime 5 (`ClassFive.lean`); the reduction from `n, m ± n` to `m`
(`ClassTable.lean`). Not in Lean: the upper bounds for `p`-adic points.
-/

namespace MagicSquareOfSquares.CoverDescent

/-! ### The class formula -/

/-- `2(w + x) · (x + iy) = (w + x + iy)²` when `x² + y² = w²`. -/
theorem GSq.of_norm {x y w : ℤ} (h : x ^ 2 + y ^ 2 = w ^ 2) : GSq (2 * (w + x)) x y :=
  ⟨w + x, y, by linear_combination h, by ring⟩

theorem GSq.scale2 {d x y : ℤ} (s : ℤ) (h : GSq d x y) : GSq (s ^ 2 * d) x y := by
  obtain ⟨a, b, ha, hb⟩ := h
  exact ⟨s * a, s * b, by linear_combination s ^ 2 * ha, by linear_combination s ^ 2 * hb⟩

/-- The product of the four numbers `w_j + x_j`: the class of the point, up to squares and sign. -/
def clsD (m n u v w0 w8 wL wU : ℤ) : ℤ :=
  (w0 + m * v) * (w8 + n * v) * (wL + (m + n) * (u - v)) * (wU + (m - n) * (u - v))

/-- **The class of a solution is `∏ (w_j + x_j)`**: `16 D` times each of the eight products is a
Gaussian square. -/
theorem gsqP_class {m n u v w0 w8 wL wU : ℤ} (h0 : nrm (P0 m n u v) = w0 ^ 2)
    (h8 : nrm (P8 m n u v) = w8 ^ 2) (hL : nrm (PL m n u v) = wL ^ 2)
    (hU : nrm (PU m n u v) = wU ^ 2) (f : Bool × Bool × Bool) :
    GSqP (16 * clsD m n u v w0 w8 wL wU) (Q m n u v f) := by
  have g0 : GSqP (2 * (w0 + m * v)) (P0 m n u v) := GSq.of_norm h0
  have g8 : GSqP (2 * (w8 + n * v)) (P8 m n u v) := GSq.of_norm h8
  have gL : GSqP (2 * (wL + (m + n) * (u - v))) (PL m n u v) := GSq.of_norm hL
  have gU : GSqP (2 * (wU + (m - n) * (u - v))) (PU m n u v) := GSq.of_norm hU
  have e : 16 * clsD m n u v w0 w8 wL wU = 2 * (w0 + m * v) * (2 * (w8 + n * v)) *
      (2 * (wL + (m + n) * (u - v))) * (2 * (wU + (m - n) * (u - v))) := by
    unfold clsD
    ring
  rw [e]
  exact ((g0.gm (g8.cjIf f.1)).gm (gL.cjIf f.2.1)).gm (gU.cjIf f.2.2)

theorem nrm_Q (m n u v : ℤ) (f : Bool × Bool × Bool) :
    nrm (Q m n u v f) =
      nrm (P0 m n u v) * nrm (P8 m n u v) * nrm (PL m n u v) * nrm (PU m n u v) := by
  obtain ⟨a, b, c⟩ := f
  cases a <;> cases b <;> cases c <;>
    simp only [Q, gm, cj, cjIf, nrm, P0, P8, PL, PU, if_true, if_false, Bool.false_eq_true] <;> ring

theorem isSquare_of_mul_sq_eq {x y z : ℤ} (hy : y ≠ 0) (h : x * y ^ 2 = z ^ 2) : IsSquare x := by
  have hd : y ∣ z := (Int.pow_dvd_pow_iff two_ne_zero).mp ⟨x, by linear_combination -h⟩
  obtain ⟨t, rfl⟩ := hd
  refine ⟨t, ?_⟩
  have h1 : y ^ 2 * x = y ^ 2 * (t * t) := by linear_combination h
  exact mul_left_cancel₀ (pow_ne_zero 2 hy) h1

/-- If `d` times one of the eight products is a Gaussian square, `d · D` or `−d · D` is a square. -/
theorem class_unique {m n u v w0 w8 wL wU d : ℤ} (h0 : nrm (P0 m n u v) = w0 ^ 2)
    (h8 : nrm (P8 m n u v) = w8 ^ 2) (hL : nrm (PL m n u v) = wL ^ 2)
    (hU : nrm (PU m n u v) = wU ^ 2) (hW : w0 * w8 * wL * wU ≠ 0) {f : Bool × Bool × Bool}
    (g : GSqP d (Q m n u v f)) :
    IsSquare (d * clsD m n u v w0 w8 wL wU) ∨ IsSquare (-(d * clsD m n u v w0 w8 wL wU)) := by
  obtain ⟨a, b, ha, hb⟩ := g
  obtain ⟨c, e, hc, he⟩ := gsqP_class h0 h8 hL hU f
  have hN : (Q m n u v f).1 ^ 2 + (Q m n u v f).2 ^ 2 = (w0 * w8 * wL * wU) ^ 2 := by
    have h1 : (Q m n u v f).1 ^ 2 + (Q m n u v f).2 ^ 2 = nrm (Q m n u v f) := rfl
    rw [h1, nrm_Q, h0, h8, hL, hU]
    ring
  generalize (Q m n u v f).1 = X at ha hb hc he hN
  generalize (Q m n u v f).2 = Y at ha hb hc he hN
  generalize clsD m n u v w0 w8 wL wU = D at hc he ⊢
  generalize w0 * w8 * wL * wU = W at hW hN
  have hz2 : 2 * ((a * c + b * e) * (b * c - a * e)) = 0 := by
    linear_combination (-(c ^ 2 - e ^ 2)) * hb - (d * Y) * hc + (2 * c * e) * ha + (d * X) * he
  have hz : (a * c + b * e) * (b * c - a * e) = 0 := by linarith
  have hS : 16 * (d * D) * W ^ 2 = (a * c + b * e) ^ 2 - (b * c - a * e) ^ 2 := by
    linear_combination (-(16 * d * D)) * hN + (16 * D * X) * ha + (a ^ 2 - b ^ 2) * hc +
      (16 * D * Y) * hb + (2 * a * b) * he
  have h4W : 4 * W ≠ 0 := mul_ne_zero (by norm_num) hW
  rcases mul_eq_zero.mp hz with h1 | h2
  · right
    exact isSquare_of_mul_sq_eq (y := 4 * W) (z := b * c - a * e) h4W
      (by linear_combination -hS - (a * c + b * e) * h1)
  · left
    exact isSquare_of_mul_sq_eq (y := 4 * W) (z := a * c + b * e) h4W
      (by linear_combination hS - (b * c - a * e) * h2)

/-- **The class formula** (Lemma 1 of `COMPLETE.md`, over the integers). At a solution on the four
forms with roots `w0, w8, wL, wU` of the norms, none zero, and `D = ∏ (w_j + x_j) ≠ 0`: a `d ≠ 0`
has the property of the eight products if and only if `d · D` or `−d · D` is a square. -/
theorem class_iff {m n u v w0 w8 wL wU d : ℤ} (h0 : nrm (P0 m n u v) = w0 ^ 2)
    (h8 : nrm (P8 m n u v) = w8 ^ 2) (hL : nrm (PL m n u v) = wL ^ 2)
    (hU : nrm (PU m n u v) = wU ^ 2) (hW : w0 * w8 * wL * wU ≠ 0)
    (hdD : d * clsD m n u v w0 w8 wL wU ≠ 0) :
    (∀ f : Bool × Bool × Bool, GSqP d (Q m n u v f)) ↔
      IsSquare (d * clsD m n u v w0 w8 wL wU) ∨ IsSquare (-(d * clsD m n u v w0 w8 wL wU)) := by
  constructor
  · intro g
    exact class_unique h0 h8 hL hU hW (g (false, false, false))
  · intro hsq f
    have base : GSq (16 * clsD m n u v w0 w8 wL wU) (Q m n u v f).1 (Q m n u v f).2 :=
      gsqP_class h0 h8 hL hU f
    show GSq d (Q m n u v f).1 (Q m n u v f).2
    generalize clsD m n u v w0 w8 wL wU = D at hdD hsq base
    rcases hsq with ⟨t, ht⟩ | ⟨t, ht⟩
    · have ht0 : t ≠ 0 := by
        rintro rfl
        exact hdD (by linear_combination ht)
      have hpos : 0 < (4 * t).natAbs := Int.natAbs_pos.mpr (mul_ne_zero (by norm_num) ht0)
      have e1 : d * (((4 * t).natAbs : ℕ) : ℤ) ^ 2 = d ^ 2 * (16 * D) := by
        rw [Int.natAbs_sq]
        linear_combination (-16 * d) * ht
      refine GSq.reduce (4 * t).natAbs hpos d _ _ ?_
      rw [e1]
      exact GSq.scale2 d base
    · have ht0 : t ≠ 0 := by
        rintro rfl
        exact hdD (by linear_combination -ht)
      have hpos : 0 < (4 * t).natAbs := Int.natAbs_pos.mpr (mul_ne_zero (by norm_num) ht0)
      have e1 : (-d) * (((4 * t).natAbs : ℕ) : ℤ) ^ 2 = d ^ 2 * (16 * D) := by
        rw [Int.natAbs_sq]
        linear_combination (16 * d) * ht
      have h1 : GSq (-d) (Q m n u v f).1 (Q m n u v f).2 := by
        refine GSq.reduce (4 * t).natAbs hpos (-d) _ _ ?_
        rw [e1]
        exact GSq.scale2 d base
      have h2 := GSq.neg_left h1
      rwa [neg_neg] at h2

/-! ### The same modulo `M` -/

/-- `d (x + i y) ≡ (a + i b)² (mod M)` for some integers `a, b`. -/
def GSqM (M d : ℤ) (p : ℤ × ℤ) : Prop :=
  ∃ a b : ℤ, d * p.1 ≡ a ^ 2 - b ^ 2 [ZMOD M] ∧ d * p.2 ≡ 2 * a * b [ZMOD M]

theorem GSqM.of_norm {M w : ℤ} {p : ℤ × ℤ} (h : nrm p ≡ w ^ 2 [ZMOD M]) :
    GSqM M (2 * (w + p.1)) p := by
  refine ⟨w + p.1, p.2, ?_, Int.ModEq.refl _⟩
  calc 2 * (w + p.1) * p.1 = (w + p.1) ^ 2 - p.2 ^ 2 + (nrm p - w ^ 2) := by unfold nrm; ring
    _ ≡ (w + p.1) ^ 2 - p.2 ^ 2 + (w ^ 2 - w ^ 2) [ZMOD M] := (h.sub_right _).add_left _
    _ = (w + p.1) ^ 2 - p.2 ^ 2 := by ring

theorem GSqM.gm {M d₁ d₂ : ℤ} {p q : ℤ × ℤ} (h₁ : GSqM M d₁ p) (h₂ : GSqM M d₂ q) :
    GSqM M (d₁ * d₂) (CoverDescent.gm p q) := by
  obtain ⟨a₁, b₁, e₁, f₁⟩ := h₁
  obtain ⟨a₂, b₂, e₂, f₂⟩ := h₂
  refine ⟨a₁ * a₂ - b₁ * b₂, a₁ * b₂ + b₁ * a₂, ?_, ?_⟩
  · calc d₁ * d₂ * (CoverDescent.gm p q).1 = d₁ * p.1 * (d₂ * q.1) - d₁ * p.2 * (d₂ * q.2) := by
          simp only [CoverDescent.gm]; ring
      _ ≡ (a₁ ^ 2 - b₁ ^ 2) * (a₂ ^ 2 - b₂ ^ 2) - 2 * a₁ * b₁ * (2 * a₂ * b₂) [ZMOD M] :=
          (e₁.mul e₂).sub (f₁.mul f₂)
      _ = (a₁ * a₂ - b₁ * b₂) ^ 2 - (a₁ * b₂ + b₁ * a₂) ^ 2 := by ring
  · calc d₁ * d₂ * (CoverDescent.gm p q).2 = d₁ * p.1 * (d₂ * q.2) + d₁ * p.2 * (d₂ * q.1) := by
          simp only [CoverDescent.gm]; ring
      _ ≡ (a₁ ^ 2 - b₁ ^ 2) * (2 * a₂ * b₂) + 2 * a₁ * b₁ * (a₂ ^ 2 - b₂ ^ 2) [ZMOD M] :=
          (e₁.mul f₂).add (f₁.mul e₂)
      _ = 2 * (a₁ * a₂ - b₁ * b₂) * (a₁ * b₂ + b₁ * a₂) := by ring

theorem GSqM.cjIf {M d : ℤ} {p : ℤ × ℤ} (h : GSqM M d p) (c : Bool) :
    GSqM M d (CoverDescent.cjIf c p) := by
  cases c
  · exact h
  · obtain ⟨a, b, e, f⟩ := h
    refine ⟨a, -b, ?_, ?_⟩
    · calc d * (CoverDescent.cjIf true p).1 = d * p.1 := by simp [CoverDescent.cjIf, CoverDescent.cj]
        _ ≡ a ^ 2 - b ^ 2 [ZMOD M] := e
        _ = a ^ 2 - (-b) ^ 2 := by ring
    · calc d * (CoverDescent.cjIf true p).2 = -(d * p.2) := by simp [CoverDescent.cjIf, CoverDescent.cj]
        _ ≡ -(2 * a * b) [ZMOD M] := f.neg
        _ = 2 * a * -b := by ring

theorem GSqM.scale2 {M d : ℤ} {p : ℤ × ℤ} (s : ℤ) (h : GSqM M d p) : GSqM M (s ^ 2 * d) p := by
  obtain ⟨a, b, e, f⟩ := h
  refine ⟨s * a, s * b, ?_, ?_⟩
  · calc s ^ 2 * d * p.1 = s ^ 2 * (d * p.1) := by ring
      _ ≡ s ^ 2 * (a ^ 2 - b ^ 2) [ZMOD M] := e.mul_left _
      _ = (s * a) ^ 2 - (s * b) ^ 2 := by ring
  · calc s ^ 2 * d * p.2 = s ^ 2 * (d * p.2) := by ring
      _ ≡ s ^ 2 * (2 * a * b) [ZMOD M] := f.mul_left _
      _ = 2 * (s * a) * (s * b) := by ring

theorem GSqM.congr {M d d' : ℤ} {p : ℤ × ℤ} (hd : d ≡ d' [ZMOD M]) (h : GSqM M d' p) :
    GSqM M d p := by
  obtain ⟨a, b, e, f⟩ := h
  exact ⟨a, b, (hd.mul_right _).trans e, (hd.mul_right _).trans f⟩

/-- **The class criterion modulo `M`.** If `w0, w8, wL, wU` are square roots of the four norms
modulo `M` and `d ≡ σ² · 16 D (mod M)` with `D = ∏ (w_j + x_j)`, then `d` times each of the eight
products is a Gaussian square modulo `M`. -/
theorem gsqM_class {M m n u v w0 w8 wL wU σ d : ℤ} (h0 : nrm (P0 m n u v) ≡ w0 ^ 2 [ZMOD M])
    (h8 : nrm (P8 m n u v) ≡ w8 ^ 2 [ZMOD M]) (hL : nrm (PL m n u v) ≡ wL ^ 2 [ZMOD M])
    (hU : nrm (PU m n u v) ≡ wU ^ 2 [ZMOD M])
    (hd : d ≡ σ ^ 2 * (16 * clsD m n u v w0 w8 wL wU) [ZMOD M]) (f : Bool × Bool × Bool) :
    GSqM M d (Q m n u v f) := by
  have g0 : GSqM M (2 * (w0 + m * v)) (P0 m n u v) := GSqM.of_norm h0
  have g8 : GSqM M (2 * (w8 + n * v)) (P8 m n u v) := GSqM.of_norm h8
  have gL : GSqM M (2 * (wL + (m + n) * (u - v))) (PL m n u v) := GSqM.of_norm hL
  have gU : GSqM M (2 * (wU + (m - n) * (u - v))) (PU m n u v) := GSqM.of_norm hU
  have e : 16 * clsD m n u v w0 w8 wL wU = 2 * (w0 + m * v) * (2 * (w8 + n * v)) *
      (2 * (wL + (m + n) * (u - v))) * (2 * (wU + (m - n) * (u - v))) := by
    unfold clsD
    ring
  have g : GSqM M (16 * clsD m n u v w0 w8 wL wU) (Q m n u v f) := by
    rw [e]
    exact ((g0.gm (g8.cjIf f.1)).gm (gL.cjIf f.2.1)).gm (gU.cjIf f.2.2)
  exact GSqM.congr hd (g.scale2 σ)

/-! ### Square roots modulo powers of an odd prime -/

theorem not_dvd_mul_of_prime {q a b : ℤ} (hq : Prime q) (ha : ¬q ∣ a) (hb : ¬q ∣ b) :
    ¬q ∣ a * b :=
  fun h => (hq.dvd_or_dvd h).elim ha hb

theorem not_dvd_two_of_odd {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) : ¬(p : ℤ) ∣ 2 := by
  intro h
  have h1 : p ∣ 2 := by exact_mod_cast h
  have h3 := Nat.le_of_dvd (by norm_num) h1
  have h4 := hp.two_le
  omega

/-- A unit that is a square modulo an odd prime `p` is a square modulo every power of `p`, with a
root in the given class modulo `p`. -/
theorem exists_sq_modEq_pow {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) {a c : ℤ} (hc : ¬(p : ℤ) ∣ c)
    (h : c ^ 2 ≡ a [ZMOD p]) (k : ℕ) :
    ∃ x : ℤ, x ^ 2 ≡ a [ZMOD (p : ℤ) ^ (k + 1)] ∧ x ≡ c [ZMOD p] := by
  have hp' : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  have h2 := not_dvd_two_of_odd hp hp2
  induction k with
  | zero =>
    refine ⟨c, ?_, Int.ModEq.refl _⟩
    rw [zero_add, pow_one]
    exact h
  | succ k ih =>
    obtain ⟨x, hx, hxc⟩ := ih
    have hxp : ¬(p : ℤ) ∣ x := fun hd =>
      hc (Int.modEq_zero_iff_dvd.mp (hxc.symm.trans (Int.modEq_zero_iff_dvd.mpr hd)))
    obtain ⟨α, β, hab⟩ : IsCoprime (2 * x) (p : ℤ) :=
      (hp'.irreducible.coprime_iff_not_dvd.mpr (not_dvd_mul_of_prime hp' h2 hxp)).symm
    obtain ⟨t, ht⟩ := hx.dvd
    refine ⟨x + (p : ℤ) ^ (k + 1) * (α * t), ?_, ?_⟩
    · apply Int.modEq_iff_dvd.mpr
      exact ⟨β * t - (p : ℤ) ^ k * (α * t) ^ 2, by
        linear_combination ht - (p : ℤ) ^ (k + 1) * t * hab⟩
    · have h0 : (p : ℤ) ^ (k + 1) * (α * t) ≡ 0 [ZMOD p] :=
        Int.modEq_zero_iff_dvd.mpr (Dvd.dvd.mul_right (dvd_pow_self _ (Nat.succ_ne_zero k)) _)
      have h1 := hxc.add h0
      rwa [add_zero] at h1

/-! ### Points of `X_p` modulo a power of `p`, with their class -/

/-- `d` is the class of the point `(u, v)` modulo `M`: the four norms have square roots modulo
`M`, none of them divisible by `q`, and `d` times each of the eight products is a Gaussian square
modulo `M`. It is used with `M = p^(k+1)`, for every `k`, and `q = p³`. The roots, with the bound
on their valuation, are what the last condition alone does not give: at a prime `≡ 3 (mod 4)` the
eight products can pass at a point where a norm is not a square (`COMPLETE.md`, section 5; the
first version of `Passes` in `ClassTable.lean` had only the last condition). -/
def ClassAt (M q d m n u v : ℤ) : Prop :=
  ∃ w0 w8 wL wU : ℤ,
    (nrm (P0 m n u v) ≡ w0 ^ 2 [ZMOD M] ∧ nrm (P8 m n u v) ≡ w8 ^ 2 [ZMOD M] ∧
      nrm (PL m n u v) ≡ wL ^ 2 [ZMOD M] ∧ nrm (PU m n u v) ≡ wU ^ 2 [ZMOD M]) ∧
    (¬q ∣ w0 ∧ ¬q ∣ w8 ∧ ¬q ∣ wL ∧ ¬q ∣ wU) ∧
    ∀ f : Bool × Bool × Bool, GSqM M d (Q m n u v f)

theorem ClassAt.gsq {M q d m n u v : ℤ} (h : ClassAt M q d m n u v) (f : Bool × Bool × Bool) :
    GSqM M d (Q m n u v f) := by
  obtain ⟨_, _, _, _, _, _, g⟩ := h
  exact g f

theorem not_dvd_of_modEq {q x c : ℤ} (h : x ≡ c [ZMOD q]) (hc : ¬q ∣ c) : ¬q ∣ x :=
  fun hd => hc (Int.modEq_zero_iff_dvd.mp (h.symm.trans (Int.modEq_zero_iff_dvd.mpr hd)))

theorem not_cube_dvd {q x : ℤ} (hx : ¬q ∣ x) : ¬q ^ 3 ∣ x :=
  fun h => hx (dvd_trans (dvd_pow_self q (by norm_num)) h)

theorem not_cube_dvd_mul {q x : ℤ} (hq : Prime q) (hx : ¬q ∣ x) : ¬q ^ 3 ∣ q * x := by
  intro h
  have e : q ^ 3 = q * q ^ 2 := by ring
  rw [e] at h
  exact hx (dvd_trans (dvd_pow_self q two_ne_zero) ((mul_dvd_mul_iff_left hq.ne_zero).mp h))

/-! ### Points of first order (section 4.1 of `COMPLETE.md`) -/

/-- **Points of first order.** Let `p` be an odd prime dividing `m` and not `n`, and `u, v` units
modulo `p` with `u + v` a unit and `2(u² + v²)` a nonzero square modulo `p`. Then every `d` such
that `d · 2uv` is a nonzero square modulo `p` passes at `(u, v)` modulo every power of `p`: the
four norms have square roots, and `d` times each of the eight products is a Gaussian square. -/
theorem pass_first_order {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) {m n u v s t d : ℤ}
    (hm : (p : ℤ) ∣ m) (hn : ¬(p : ℤ) ∣ n) (hu : ¬(p : ℤ) ∣ u) (hv : ¬(p : ℤ) ∣ v)
    (huv : ¬(p : ℤ) ∣ u + v) (hs : s ^ 2 ≡ 2 * (u ^ 2 + v ^ 2) [ZMOD p]) (hs0 : ¬(p : ℤ) ∣ s)
    (ht : t ^ 2 ≡ d * (2 * u * v) [ZMOD p]) (ht0 : ¬(p : ℤ) ∣ t) (k : ℕ) :
    ClassAt ((p : ℤ) ^ (k + 1)) ((p : ℤ) ^ 3) d m n u v := by
  have hp' : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  have h2 := not_dvd_two_of_odd hp hp2
  have hm0 : m ≡ 0 [ZMOD p] := Int.modEq_zero_iff_dvd.mpr hm
  -- the four norms modulo `p`
  have e0 : (n * u) ^ 2 ≡ nrm (P0 m n u v) [ZMOD p] := by
    symm
    calc nrm (P0 m n u v) = (m * v) ^ 2 + (n * u) ^ 2 := rfl
      _ ≡ (0 * v) ^ 2 + (n * u) ^ 2 [ZMOD p] := ((hm0.mul_right v).pow 2).add_right _
      _ = (n * u) ^ 2 := by ring
  have e8 : (n * v) ^ 2 ≡ nrm (P8 m n u v) [ZMOD p] := by
    symm
    calc nrm (P8 m n u v) = (n * v) ^ 2 + (m * u) ^ 2 := rfl
      _ ≡ (n * v) ^ 2 + (0 * u) ^ 2 [ZMOD p] := ((hm0.mul_right u).pow 2).add_left _
      _ = (n * v) ^ 2 := by ring
  have eL : (n * s) ^ 2 ≡ nrm (PL m n u v) [ZMOD p] := by
    symm
    calc nrm (PL m n u v) = ((m + n) * (u - v)) ^ 2 + ((m - n) * (u + v)) ^ 2 := rfl
      _ ≡ ((0 + n) * (u - v)) ^ 2 + ((0 - n) * (u + v)) ^ 2 [ZMOD p] :=
          (((hm0.add_right n).mul_right _).pow 2).add (((hm0.sub_right n).mul_right _).pow 2)
      _ = n ^ 2 * (2 * (u ^ 2 + v ^ 2)) := by ring
      _ ≡ n ^ 2 * s ^ 2 [ZMOD p] := hs.symm.mul_left _
      _ = (n * s) ^ 2 := by ring
  have eU : (n * s) ^ 2 ≡ nrm (PU m n u v) [ZMOD p] := by
    symm
    calc nrm (PU m n u v) = ((m - n) * (u - v)) ^ 2 + ((m + n) * (u + v)) ^ 2 := rfl
      _ ≡ ((0 - n) * (u - v)) ^ 2 + ((0 + n) * (u + v)) ^ 2 [ZMOD p] :=
          (((hm0.sub_right n).mul_right _).pow 2).add (((hm0.add_right n).mul_right _).pow 2)
      _ = n ^ 2 * (2 * (u ^ 2 + v ^ 2)) := by ring
      _ ≡ n ^ 2 * s ^ 2 [ZMOD p] := hs.symm.mul_left _
      _ = (n * s) ^ 2 := by ring
  obtain ⟨w0, h0, c0⟩ := exists_sq_modEq_pow hp hp2 (not_dvd_mul_of_prime hp' hn hu) e0 k
  obtain ⟨w8, h8, c8⟩ := exists_sq_modEq_pow hp hp2 (not_dvd_mul_of_prime hp' hn hv) e8 k
  obtain ⟨wL, hL, cL⟩ := exists_sq_modEq_pow hp hp2 (not_dvd_mul_of_prime hp' hn hs0) eL k
  obtain ⟨wU, hU, cU⟩ := exists_sq_modEq_pow hp hp2 (not_dvd_mul_of_prime hp' hn hs0) eU k
  -- the class modulo `p`
  have hD : clsD m n u v w0 w8 wL wU ≡
      (n * u + 0 * v) * (n * v + n * v) * (n * s + (0 + n) * (u - v)) *
        (n * s + (0 - n) * (u - v)) [ZMOD p] := by
    unfold clsD
    exact (((c0.add (hm0.mul_right v)).mul (c8.add_right _)).mul
      (cL.add ((hm0.add_right n).mul_right _))).mul (cU.add ((hm0.sub_right n).mul_right _))
  have hE : (4 * n * n * (u + v) * t) ^ 2 ≡ d * (16 * clsD m n u v w0 w8 wL wU) [ZMOD p] := by
    symm
    calc d * (16 * clsD m n u v w0 w8 wL wU)
        ≡ d * (16 * ((n * u + 0 * v) * (n * v + n * v) * (n * s + (0 + n) * (u - v)) *
            (n * s + (0 - n) * (u - v)))) [ZMOD p] := (hD.mul_left 16).mul_left d
      _ = 16 * n ^ 4 * (d * (2 * u * v)) * (s ^ 2 - (u - v) ^ 2) := by ring
      _ ≡ 16 * n ^ 4 * t ^ 2 * (2 * (u ^ 2 + v ^ 2) - (u - v) ^ 2) [ZMOD p] :=
          (ht.symm.mul_left _).mul (hs.sub_right _)
      _ = (4 * n * n * (u + v) * t) ^ 2 := by ring
  have h4 : ¬(p : ℤ) ∣ 4 := by
    have e4 : (4 : ℤ) = 2 * 2 := by norm_num
    rw [e4]
    exact not_dvd_mul_of_prime hp' h2 h2
  have hg : ¬(p : ℤ) ∣ 4 * n * n * (u + v) * t :=
    not_dvd_mul_of_prime hp' (not_dvd_mul_of_prime hp'
      (not_dvd_mul_of_prime hp' (not_dvd_mul_of_prime hp' h4 hn) hn) huv) ht0
  obtain ⟨x, hx, hxg⟩ := exists_sq_modEq_pow hp hp2 hg hE k
  have hxp : ¬(p : ℤ) ∣ x := fun hd =>
    hg (Int.modEq_zero_iff_dvd.mp (hxg.symm.trans (Int.modEq_zero_iff_dvd.mpr hd)))
  obtain ⟨ρ, μ, hρ⟩ : IsCoprime x ((p : ℤ) ^ (k + 1)) :=
    (hp'.irreducible.coprime_iff_not_dvd.mpr hxp).symm.pow_right
  have hρ1 : ρ * x ≡ 1 [ZMOD (p : ℤ) ^ (k + 1)] :=
    Int.modEq_iff_dvd.mpr ⟨μ, by linear_combination -hρ⟩
  have hd : d ≡ (ρ * d) ^ 2 * (16 * clsD m n u v w0 w8 wL wU) [ZMOD (p : ℤ) ^ (k + 1)] := by
    calc d = d * 1 ^ 2 := by ring
      _ ≡ d * (ρ * x) ^ 2 [ZMOD (p : ℤ) ^ (k + 1)] := ((hρ1.pow 2).symm).mul_left d
      _ = d * ρ ^ 2 * x ^ 2 := by ring
      _ ≡ d * ρ ^ 2 * (d * (16 * clsD m n u v w0 w8 wL wU)) [ZMOD (p : ℤ) ^ (k + 1)] :=
          hx.mul_left _
      _ = (ρ * d) ^ 2 * (16 * clsD m n u v w0 w8 wL wU) := by ring
  exact ⟨w0, w8, wL, wU, ⟨h0.symm, h8.symm, hL.symm, hU.symm⟩,
    ⟨not_cube_dvd (not_dvd_of_modEq c0 (not_dvd_mul_of_prime hp' hn hu)),
      not_cube_dvd (not_dvd_of_modEq c8 (not_dvd_mul_of_prime hp' hn hv)),
      not_cube_dvd (not_dvd_of_modEq cL (not_dvd_mul_of_prime hp' hn hs0)),
      not_cube_dvd (not_dvd_of_modEq cU (not_dvd_mul_of_prime hp' hn hs0))⟩,
    fun f => gsqM_class h0.symm h8.symm hL.symm hU.symm hd f⟩

/-- **The class 2 is never excluded at an odd prime of `m`**: at `(u, v) = (1, 1)`, modulo every
power of `p`. -/
theorem two_passes {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) {m n : ℤ} (hm : (p : ℤ) ∣ m)
    (hn : ¬(p : ℤ) ∣ n) (k : ℕ) :
    ClassAt ((p : ℤ) ^ (k + 1)) ((p : ℤ) ^ 3) 2 m n 1 1 := by
  have hp' : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  have h2 := not_dvd_two_of_odd hp hp2
  have h1 : ¬(p : ℤ) ∣ 1 := fun h => hp'.not_isUnit (isUnit_of_dvd_one h)
  exact pass_first_order hp hp2 (u := 1) (v := 1) (s := 2) (t := 2)
    (d := 2) hm hn h1 h1 (by norm_num; exact h2) (by norm_num) h2
    (by norm_num) h2 k

/-! ### Classes of odd valuation (section 4.3 of `COMPLETE.md`, the case `p² | m`) -/

/-- **Classes of odd valuation.** Let `p` be an odd prime modulo which 2 is a square, `m = p²·m₂`,
`p ∤ n`, and `u = p·u₁` with `u₁` a unit, `v = 1`. Then every `d = p·d₁` such that `d₁ · 2u₁` is a
nonzero square modulo `p` passes at `(u, v)` modulo every power of `p`. This is Theorem E read
backwards: a prime `≡ ±1 (mod 8)` can divide the class. -/
theorem pass_odd_valuation {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) {m2 n u1 r t d1 : ℤ}
    (hn : ¬(p : ℤ) ∣ n) (hu : ¬(p : ℤ) ∣ u1) (hr : r ^ 2 ≡ 2 [ZMOD p])
    (ht : t ^ 2 ≡ d1 * (2 * u1) [ZMOD p]) (ht0 : ¬(p : ℤ) ∣ t) (k : ℕ) :
    ClassAt ((p : ℤ) ^ (k + 1)) ((p : ℤ) ^ 3) ((p : ℤ) * d1) ((p : ℤ) ^ 2 * m2) n ((p : ℤ) * u1) 1 := by
  have hp' : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  have h2 := not_dvd_two_of_odd hp hp2
  have hp0 : (p : ℤ) ≡ 0 [ZMOD p] := Int.modEq_zero_iff_dvd.mpr dvd_rfl
  have hr0 : ¬(p : ℤ) ∣ r := by
    intro h
    have h1 : (p : ℤ) ∣ r ^ 2 := Dvd.dvd.pow h two_ne_zero
    exact h2 (Int.modEq_zero_iff_dvd.mp (hr.symm.trans (Int.modEq_zero_iff_dvd.mpr h1)))
  -- the four norms modulo `p`; the first is `p²` times a unit square
  have e0 : (n * u1) ^ 2 ≡ (n * u1) ^ 2 + ((p : ℤ) * m2) ^ 2 [ZMOD p] := by
    calc (n * u1) ^ 2 = (n * u1) ^ 2 + (0 * m2) ^ 2 := by ring
      _ ≡ (n * u1) ^ 2 + ((p : ℤ) * m2) ^ 2 [ZMOD p] := ((hp0.symm.mul_right m2).pow 2).add_left _
  have e8 : n ^ 2 ≡ nrm (P8 ((p : ℤ) ^ 2 * m2) n ((p : ℤ) * u1) 1) [ZMOD p] := by
    symm
    calc nrm (P8 ((p : ℤ) ^ 2 * m2) n ((p : ℤ) * u1) 1)
        = (n * 1) ^ 2 + ((p : ℤ) ^ 2 * m2 * ((p : ℤ) * u1)) ^ 2 := rfl
      _ ≡ (n * 1) ^ 2 + ((0 : ℤ) ^ 2 * m2 * (0 * u1)) ^ 2 [ZMOD p] :=
          ((((hp0.pow 2).mul_right m2).mul (hp0.mul_right u1)).pow 2).add_left _
      _ = n ^ 2 := by ring
  have eL : (n * r) ^ 2 ≡ nrm (PL ((p : ℤ) ^ 2 * m2) n ((p : ℤ) * u1) 1) [ZMOD p] := by
    symm
    calc nrm (PL ((p : ℤ) ^ 2 * m2) n ((p : ℤ) * u1) 1)
        = (((p : ℤ) ^ 2 * m2 + n) * ((p : ℤ) * u1 - 1)) ^ 2 +
            (((p : ℤ) ^ 2 * m2 - n) * ((p : ℤ) * u1 + 1)) ^ 2 := rfl
      _ ≡ (((0 : ℤ) ^ 2 * m2 + n) * (0 * u1 - 1)) ^ 2 +
            (((0 : ℤ) ^ 2 * m2 - n) * (0 * u1 + 1)) ^ 2 [ZMOD p] :=
          (((((hp0.pow 2).mul_right m2).add_right n).mul ((hp0.mul_right u1).sub_right 1)).pow 2).add
            (((((hp0.pow 2).mul_right m2).sub_right n).mul ((hp0.mul_right u1).add_right 1)).pow 2)
      _ = n ^ 2 * 2 := by ring
      _ ≡ n ^ 2 * r ^ 2 [ZMOD p] := hr.symm.mul_left _
      _ = (n * r) ^ 2 := by ring
  have eU : (n * r) ^ 2 ≡ nrm (PU ((p : ℤ) ^ 2 * m2) n ((p : ℤ) * u1) 1) [ZMOD p] := by
    symm
    calc nrm (PU ((p : ℤ) ^ 2 * m2) n ((p : ℤ) * u1) 1)
        = (((p : ℤ) ^ 2 * m2 - n) * ((p : ℤ) * u1 - 1)) ^ 2 +
            (((p : ℤ) ^ 2 * m2 + n) * ((p : ℤ) * u1 + 1)) ^ 2 := rfl
      _ ≡ (((0 : ℤ) ^ 2 * m2 - n) * (0 * u1 - 1)) ^ 2 +
            (((0 : ℤ) ^ 2 * m2 + n) * (0 * u1 + 1)) ^ 2 [ZMOD p] :=
          (((((hp0.pow 2).mul_right m2).sub_right n).mul ((hp0.mul_right u1).sub_right 1)).pow 2).add
            (((((hp0.pow 2).mul_right m2).add_right n).mul ((hp0.mul_right u1).add_right 1)).pow 2)
      _ = n ^ 2 * 2 := by ring
      _ ≡ n ^ 2 * r ^ 2 [ZMOD p] := hr.symm.mul_left _
      _ = (n * r) ^ 2 := by ring
  obtain ⟨ω0, hω, cω⟩ := exists_sq_modEq_pow hp hp2 (not_dvd_mul_of_prime hp' hn hu) e0 k
  obtain ⟨w8, h8, c8⟩ := exists_sq_modEq_pow hp hp2 hn e8 k
  obtain ⟨wL, hL, cL⟩ := exists_sq_modEq_pow hp hp2 (not_dvd_mul_of_prime hp' hn hr0) eL k
  obtain ⟨wU, hU, cU⟩ := exists_sq_modEq_pow hp hp2 (not_dvd_mul_of_prime hp' hn hr0) eU k
  have h0 : ((p : ℤ) * ω0) ^ 2 ≡ nrm (P0 ((p : ℤ) ^ 2 * m2) n ((p : ℤ) * u1) 1)
      [ZMOD (p : ℤ) ^ (k + 1)] := by
    calc ((p : ℤ) * ω0) ^ 2 = (p : ℤ) ^ 2 * ω0 ^ 2 := by ring
      _ ≡ (p : ℤ) ^ 2 * ((n * u1) ^ 2 + ((p : ℤ) * m2) ^ 2) [ZMOD (p : ℤ) ^ (k + 1)] :=
          hω.mul_left _
      _ = nrm (P0 ((p : ℤ) ^ 2 * m2) n ((p : ℤ) * u1) 1) := by simp only [nrm, P0]; ring
  -- the class is `p` times a unit `D1`
  obtain ⟨D1, hD1⟩ : ∃ D1 : ℤ, D1 = (ω0 + (p : ℤ) * m2) * (w8 + n * 1) *
      (wL + ((p : ℤ) ^ 2 * m2 + n) * ((p : ℤ) * u1 - 1)) *
      (wU + ((p : ℤ) ^ 2 * m2 - n) * ((p : ℤ) * u1 - 1)) := ⟨_, rfl⟩
  have eD : clsD ((p : ℤ) ^ 2 * m2) n ((p : ℤ) * u1) 1 ((p : ℤ) * ω0) w8 wL wU = (p : ℤ) * D1 := by
    unfold clsD
    rw [hD1]
    ring
  have cD1 : D1 ≡ (n * u1 + 0 * m2) * (n + n * 1) *
      (n * r + ((0 : ℤ) ^ 2 * m2 + n) * (0 * u1 - 1)) *
      (n * r + ((0 : ℤ) ^ 2 * m2 - n) * (0 * u1 - 1)) [ZMOD p] := by
    rw [hD1]
    exact (((cω.add (hp0.mul_right m2)).mul (c8.add_right _)).mul
      (cL.add ((((hp0.pow 2).mul_right m2).add_right n).mul ((hp0.mul_right u1).sub_right 1)))).mul
      (cU.add ((((hp0.pow 2).mul_right m2).sub_right n).mul ((hp0.mul_right u1).sub_right 1)))
  have hE : (4 * n * n * t) ^ 2 ≡ d1 * (16 * D1) [ZMOD p] := by
    symm
    calc d1 * (16 * D1)
        ≡ d1 * (16 * ((n * u1 + 0 * m2) * (n + n * 1) *
            (n * r + ((0 : ℤ) ^ 2 * m2 + n) * (0 * u1 - 1)) *
            (n * r + ((0 : ℤ) ^ 2 * m2 - n) * (0 * u1 - 1)))) [ZMOD p] :=
          (cD1.mul_left 16).mul_left d1
      _ = 16 * n ^ 4 * (d1 * (2 * u1)) * (r ^ 2 - 1) := by ring
      _ ≡ 16 * n ^ 4 * t ^ 2 * (2 - 1) [ZMOD p] := (ht.symm.mul_left _).mul (hr.sub_right _)
      _ = (4 * n * n * t) ^ 2 := by ring
  have h4 : ¬(p : ℤ) ∣ 4 := by
    have e4 : (4 : ℤ) = 2 * 2 := by norm_num
    rw [e4]
    exact not_dvd_mul_of_prime hp' h2 h2
  have hg : ¬(p : ℤ) ∣ 4 * n * n * t :=
    not_dvd_mul_of_prime hp' (not_dvd_mul_of_prime hp' (not_dvd_mul_of_prime hp' h4 hn) hn) ht0
  obtain ⟨x, hx, hxg⟩ := exists_sq_modEq_pow hp hp2 hg hE k
  have hxp : ¬(p : ℤ) ∣ x := fun hd =>
    hg (Int.modEq_zero_iff_dvd.mp (hxg.symm.trans (Int.modEq_zero_iff_dvd.mpr hd)))
  obtain ⟨ρ, μ, hρ⟩ : IsCoprime x ((p : ℤ) ^ (k + 1)) :=
    (hp'.irreducible.coprime_iff_not_dvd.mpr hxp).symm.pow_right
  have hρ1 : ρ * x ≡ 1 [ZMOD (p : ℤ) ^ (k + 1)] :=
    Int.modEq_iff_dvd.mpr ⟨μ, by linear_combination -hρ⟩
  have hd : (p : ℤ) * d1 ≡ (ρ * d1) ^ 2 *
      (16 * clsD ((p : ℤ) ^ 2 * m2) n ((p : ℤ) * u1) 1 ((p : ℤ) * ω0) w8 wL wU)
      [ZMOD (p : ℤ) ^ (k + 1)] := by
    rw [eD]
    calc (p : ℤ) * d1 = (p : ℤ) * d1 * 1 ^ 2 := by ring
      _ ≡ (p : ℤ) * d1 * (ρ * x) ^ 2 [ZMOD (p : ℤ) ^ (k + 1)] := ((hρ1.pow 2).symm).mul_left _
      _ = (p : ℤ) * d1 * ρ ^ 2 * x ^ 2 := by ring
      _ ≡ (p : ℤ) * d1 * ρ ^ 2 * (d1 * (16 * D1)) [ZMOD (p : ℤ) ^ (k + 1)] := hx.mul_left _
      _ = (ρ * d1) ^ 2 * (16 * ((p : ℤ) * D1)) := by ring
  exact ⟨(p : ℤ) * ω0, w8, wL, wU, ⟨h0.symm, h8.symm, hL.symm, hU.symm⟩,
    ⟨not_cube_dvd_mul hp' (not_dvd_of_modEq cω (not_dvd_mul_of_prime hp' hn hu)),
      not_cube_dvd (not_dvd_of_modEq c8 hn),
      not_cube_dvd (not_dvd_of_modEq cL (not_dvd_mul_of_prime hp' hn hr0)),
      not_cube_dvd (not_dvd_of_modEq cU (not_dvd_mul_of_prime hp' hn hr0))⟩,
    fun f => gsqM_class h0.symm h8.symm hL.symm hU.symm hd f⟩

/-- **The class `2p` is never excluded at an odd prime `p` with `p² | m` modulo which 2 is a
square**: at `(u, v) = (p, 1)`, modulo every power of `p`. -/
theorem two_mul_prime_passes {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) {m2 n r : ℤ}
    (hn : ¬(p : ℤ) ∣ n) (hr : r ^ 2 ≡ 2 [ZMOD p]) (k : ℕ) :
    ClassAt ((p : ℤ) ^ (k + 1)) ((p : ℤ) ^ 3) ((p : ℤ) * 2) ((p : ℤ) ^ 2 * m2) n ((p : ℤ) * 1) 1 := by
  have hp' : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  have h2 := not_dvd_two_of_odd hp hp2
  have h1 : ¬(p : ℤ) ∣ 1 := fun h => hp'.not_isUnit (isUnit_of_dvd_one h)
  exact pass_odd_valuation hp hp2 (u1 := 1) (d1 := 2) (t := 2)
    (m2 := m2) hn h1 hr (by norm_num) h2 k

/-! ### Two rows of the table of Theorem C, lower bounds, with no further hypothesis -/

theorem GSqM.neg_left {M d : ℤ} {p : ℤ × ℤ} (h : GSqM M d p) : GSqM M (-d) p := by
  obtain ⟨a, b, e, f⟩ := h
  refine ⟨b, -a, ?_, ?_⟩
  · calc -d * p.1 = -(d * p.1) := by ring
      _ ≡ -(a ^ 2 - b ^ 2) [ZMOD M] := e.neg
      _ = b ^ 2 - (-a) ^ 2 := by ring
  · calc -d * p.2 = -(d * p.2) := by ring
      _ ≡ -(2 * a * b) [ZMOD M] := f.neg
      _ = 2 * b * -a := by ring

theorem ClassAt.neg {M q d m n u v : ℤ} (h : ClassAt M q d m n u v) :
    ClassAt M q (-d) m n u v := by
  obtain ⟨w0, w8, wL, wU, hr, hn, g⟩ := h
  exact ⟨w0, w8, wL, wU, hr, hn, fun f => (g f).neg_left⟩

/-- **Every class `p · (unit)` passes** at an odd prime `p` modulo which 2 is a square, when
`p² | m`: `d = p·d₁` passes at `(u, v) = (2p·d₁, 1)` modulo every power of `p`. -/
theorem prime_mul_unit_passes {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) {m2 n r d1 : ℤ}
    (hn : ¬(p : ℤ) ∣ n) (hd1 : ¬(p : ℤ) ∣ d1) (hr : r ^ 2 ≡ 2 [ZMOD p]) (k : ℕ) :
    ClassAt ((p : ℤ) ^ (k + 1)) ((p : ℤ) ^ 3) ((p : ℤ) * d1)
      ((p : ℤ) ^ 2 * m2) n ((p : ℤ) * (2 * d1)) 1 := by
  have hp' : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  have h2 := not_dvd_two_of_odd hp hp2
  have hu : ¬(p : ℤ) ∣ 2 * d1 := not_dvd_mul_of_prime hp' h2 hd1
  have ht : (2 * d1) ^ 2 ≡ d1 * (2 * (2 * d1)) [ZMOD p] := by
    have e : (2 * d1) ^ 2 = d1 * (2 * (2 * d1)) := by ring
    rw [e]
  exact pass_odd_valuation hp hp2 (u1 := 2 * d1) (d1 := d1)
    (t := 2 * d1) (m2 := m2) hn hu hr ht hu k

/-- **Every unit class passes** at a prime `p ≡ 3 (mod 4)` dividing `m`: a `d` not divisible by `p`
passes at `(u, v) = (1, 1)` modulo every power of `p` (`d` or `−d` is twice a square modulo `p`,
and `d`, `−d` pass together). With Theorem E this is the row `p ≡ 3 (mod 8)` of the table. -/
theorem unit_passes_of_three_mod_four {p : ℕ} (hp : p.Prime) (h3 : p % 4 = 3) {m n d : ℤ}
    (hm : (p : ℤ) ∣ m) (hn : ¬(p : ℤ) ∣ n) (hd : ¬(p : ℤ) ∣ d) (k : ℕ) :
    ClassAt ((p : ℤ) ^ (k + 1)) ((p : ℤ) ^ 3) d m n 1 1 := by
  have : Fact p.Prime := ⟨hp⟩
  have hp2 : p ≠ 2 := by omega
  have hp' : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  have h2 := not_dvd_two_of_odd hp hp2
  have h1 : ¬(p : ℤ) ∣ 1 := fun h => hp'.not_isUnit (isUnit_of_dvd_one h)
  have h11 : ¬(p : ℤ) ∣ 1 + 1 := by norm_num; exact h2
  have hs : (2 : ℤ) ^ 2 ≡ 2 * (1 ^ 2 + 1 ^ 2) [ZMOD p] := by norm_num
  -- `2d` or `−2d` is a square modulo `p`
  have ha0 : ((2 * d : ℤ) : ZMod p) ≠ 0 := by
    rw [Ne, ZMod.intCast_zmod_eq_zero_iff_dvd]
    exact not_dvd_mul_of_prime hp' h2 hd
  have key : ∀ c : ℤ, IsSquare ((c : ℤ) : ZMod p) → ((c : ℤ) : ZMod p) ≠ 0 →
      ∃ t : ℤ, t ^ 2 ≡ c [ZMOD p] ∧ ¬(p : ℤ) ∣ t := by
    intro c hc hc0
    obtain ⟨r, hr⟩ := hc
    obtain ⟨t, rfl⟩ := ZMod.intCast_surjective r
    refine ⟨t, ?_, ?_⟩
    · rw [← ZMod.intCast_eq_intCast_iff]
      push_cast
      rw [hr]
      ring
    · intro h
      apply hc0
      rw [hr, (ZMod.intCast_zmod_eq_zero_iff_dvd t p).mpr h]
      ring
  by_cases hsq : IsSquare (((2 * d : ℤ) : ℤ) : ZMod p)
  · obtain ⟨t, ht, ht0⟩ := key (2 * d) hsq ha0
    have ht' : t ^ 2 ≡ d * (2 * 1 * 1) [ZMOD p] := by
      have e : d * (2 * 1 * 1) = 2 * d := by ring
      rw [e]
      exact ht
    exact pass_first_order hp hp2 (u := 1) (v := 1) (s := 2)
      hm hn h1 h1 h11 hs h2 ht' ht0 k
  · have hm1 : quadraticChar (ZMod p) (-1) = -1 := by
      rw [quadraticChar_neg_one_iff_not_isSquare, ZMod.exists_sq_eq_neg_one_iff]
      omega
    have ha1 : quadraticChar (ZMod p) (((2 * d : ℤ) : ℤ) : ZMod p) = -1 :=
      quadraticChar_neg_one_iff_not_isSquare.mpr hsq
    have hneg : ((-(2 * d) : ℤ) : ZMod p) = -1 * ((2 * d : ℤ) : ZMod p) := by
      push_cast
      ring
    have hq : quadraticChar (ZMod p) ((-(2 * d) : ℤ) : ZMod p) = 1 := by
      rw [hneg, map_mul, hm1, ha1]
      norm_num
    have hna0 : ((-(2 * d) : ℤ) : ZMod p) ≠ 0 := by
      rw [hneg]
      exact mul_ne_zero (by norm_num) ha0
    obtain ⟨t, ht, ht0⟩ := key (-(2 * d)) ((quadraticChar_one_iff_isSquare hna0).mp hq) hna0
    have ht' : t ^ 2 ≡ -d * (2 * 1 * 1) [ZMOD p] := by
      have e : -d * (2 * 1 * 1) = -(2 * d) := by ring
      rw [e]
      exact ht
    have h5 := (pass_first_order hp hp2 (u := 1) (v := 1) (s := 2)
      (d := -d) hm hn h1 h1 h11 hs h2 ht' ht0 k).neg
    rwa [neg_neg] at h5

theorem exists_sq_modEq_two {p : ℕ} (hp : p.Prime) (h8 : p % 8 = 1 ∨ p % 8 = 7) :
    ∃ r : ℤ, r ^ 2 ≡ 2 [ZMOD p] := by
  have : Fact p.Prime := ⟨hp⟩
  have hp2 : p ≠ 2 := by omega
  obtain ⟨r, hr⟩ := (ZMod.exists_sq_eq_two_iff hp2).mpr h8
  obtain ⟨t, rfl⟩ := ZMod.intCast_surjective r
  refine ⟨t, ?_⟩
  rw [← ZMod.intCast_eq_intCast_iff]
  push_cast
  rw [hr]
  ring

/-- **The classes of odd valuation at a prime `p ≡ ±1 (mod 8)` with `p² | m`**: every `d = p·d₁`
with `p ∤ d₁` passes modulo every power of `p`. With `unit_passes_of_three_mod_four` this is the
row `p ≡ 7 (mod 8)` of the table of Theorem C when `p² | m`: both classes occur. -/
theorem prime_mul_unit_passes_mod_eight {p : ℕ} (hp : p.Prime) (h8 : p % 8 = 1 ∨ p % 8 = 7)
    {m2 n d1 : ℤ} (hn : ¬(p : ℤ) ∣ n) (hd1 : ¬(p : ℤ) ∣ d1) (k : ℕ) :
    ClassAt ((p : ℤ) ^ (k + 1)) ((p : ℤ) ^ 3) ((p : ℤ) * d1)
      ((p : ℤ) ^ 2 * m2) n ((p : ℤ) * (2 * d1)) 1 := by
  obtain ⟨r, hr⟩ := exists_sq_modEq_two hp h8
  exact prime_mul_unit_passes hp (by omega) hn hd1 hr k

end MagicSquareOfSquares.CoverDescent
