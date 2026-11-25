import EllipticCurves.MordellWeil
import EllipticCurves.Mathlib.SelmerGroup
import MagicSquareOfSquares.TwoDescentCertificates

/-!
# `rank E(ℚ) ≥ 2` for the twice-square curves

Written by an AI model; not read by a person.  Experiment, not counted.

`rank` here is `Module.finrank ℤ E(ℚ)`, the notion used by Michael Stoll's `EllipticCurves`
(`card_range_μ`, `finrank_point_le_one`): for a finitely generated abelian group it is the rank of
its free part.  Inputs:

* Mordell–Weil: `WeierstrassCurve.Affine.fg_point_of_numberField` (Stoll, `EllipticCurves`
  2164fb5): `E(K)` is finitely generated for an elliptic curve over a number field.
* `AddSubgroup.index_range_nsmul_of_fg` (Stoll, `EllipticCurves.Mathlib.SelmerGroup`):
  `[G : 2G] = 2 ^ finrank ℤ G · #G[2]` for a finitely generated abelian group `G`.
* `Eplus_isElliptic`, `two_torsion_Eplus` (`E⁺(ℚ)[2] ⊆ {O, T}`) and `Theorem2At` (independence of
  `R`, `Q`, `T` in `E(ℚ)/2E(ℚ)`) from `MagicSquareOfSquares.TwoDescentMap`, with the 254
  instances of `MagicSquareOfSquares.TwoDescentCertificates`.

Then `8 ≤ [E : 2E] = 2 ^ rank · #E[2] ≤ 2 ^ rank · 2`, so `rank ≥ 2`.
-/

namespace MagicSquareOfSquares.Rank

open Module

section GroupTheory

variable {G : Type*} [AddCommGroup G]

/-- Two elements of `Fin 2` whose difference in `ℤ` is even are equal. -/
theorem fin2_eq_of_even {a b : Fin 2} (h : Even (((b : ℕ) : ℤ) - ((a : ℕ) : ℤ))) : a = b := by
  obtain ⟨k, hk⟩ := h
  have ha := a.isLt
  have hb := b.isLt
  apply Fin.ext
  omega

/-- The eight sums `ε₁R + ε₂Q + ε₃T`, `εᵢ ∈ {0, 1}`. -/
def eightSums (R Q T : G) (c : Fin 2 × Fin 2 × Fin 2) : G :=
  (c.1 : ℕ) • R + (c.2.1 : ℕ) • Q + (c.2.2 : ℕ) • T

/-- Independence modulo `2G` (with integer coefficients) makes the eight sums pairwise incongruent
modulo `2G`. -/
theorem eightSums_injective {R Q T : G}
    (hind : ∀ c₁ c₂ c₃ : ℤ, ∀ W : G, c₁ • R + c₂ • Q + c₃ • T = W + W →
      Even c₁ ∧ Even c₂ ∧ Even c₃) :
    Function.Injective
      (fun c => (eightSums R Q T c : G ⧸ (nsmulAddMonoidHom (α := G) 2).range)) := by
  intro c c' h
  have h' : -(eightSums R Q T c) + eightSums R Q T c' ∈ (nsmulAddMonoidHom (α := G) 2).range :=
    QuotientAddGroup.eq.mp h
  obtain ⟨W, hW⟩ := AddMonoidHom.mem_range.mp h'
  have hW' : 2 • W = -(eightSums R Q T c) + eightSums R Q T c' := hW
  have key : (((c'.1 : ℕ) : ℤ) - ((c.1 : ℕ) : ℤ)) • R + (((c'.2.1 : ℕ) : ℤ) - ((c.2.1 : ℕ) : ℤ)) • Q
      + (((c'.2.2 : ℕ) : ℤ) - ((c.2.2 : ℕ) : ℤ)) • T = W + W := by
    rw [show W + W = 2 • W from (two_nsmul W).symm, hW']
    simp only [eightSums, sub_zsmul, natCast_zsmul]
    abel
  obtain ⟨h1, h2, h3⟩ := hind _ _ _ W key
  exact Prod.ext (fin2_eq_of_even h1) (Prod.ext (fin2_eq_of_even h2) (fin2_eq_of_even h3))

/-- **`#(G/2G) ≥ 8`** when `R`, `Q`, `T` are independent modulo `2G` (and `G/2G` is finite). -/
theorem eight_le_card_quotient {R Q T : G}
    (hind : ∀ c₁ c₂ c₃ : ℤ, ∀ W : G, c₁ • R + c₂ • Q + c₃ • T = W + W →
      Even c₁ ∧ Even c₂ ∧ Even c₃)
    [Finite (G ⧸ (nsmulAddMonoidHom (α := G) 2).range)] :
    8 ≤ Nat.card (G ⧸ (nsmulAddMonoidHom (α := G) 2).range) := by
  have h := Nat.card_le_card_of_injective _ (eightSums_injective hind)
  have h8 : Nat.card (Fin 2 × Fin 2 × Fin 2) = 8 := by
    rw [Nat.card_eq_fintype_card, Fintype.card_prod, Fintype.card_prod, Fintype.card_fin]
  rwa [h8] at h

/-- **Group theory.** In a finitely generated abelian group `G` whose 2-torsion is contained in
`{0, T₀}`, three elements independent modulo `2G` force `finrank ℤ G ≥ 2`. -/
theorem two_le_finrank_of_indep [IsAddFG G] {R Q T T₀ : G}
    (htors : ∀ P : G, P + P = 0 → P = 0 ∨ P = T₀)
    (hind : ∀ c₁ c₂ c₃ : ℤ, ∀ W : G, c₁ • R + c₂ • Q + c₃ • T = W + W →
      Even c₁ ∧ Even c₂ ∧ Even c₃) :
    2 ≤ finrank ℤ G := by
  set K := (nsmulAddMonoidHom (α := G) 2).ker with hKdef
  have hsub : (K : Set G) ⊆ {0, T₀} := by
    intro P hP
    have hP' : (nsmulAddMonoidHom (α := G) 2) P = 0 := AddMonoidHom.mem_ker.mp hP
    have hPP : P + P = 0 := by rw [← two_nsmul]; exact hP'
    rcases htors P hPP with h | h
    · exact Or.inl h
    · exact Or.inr h
  have hKfin : (K : Set G).Finite := (Set.toFinite _).subset hsub
  have hK2 : Nat.card K ≤ 2 := by
    have h1 : Nat.card K = (K : Set G).ncard := Nat.card_coe_set_eq _
    rw [h1]
    calc (K : Set G).ncard ≤ ({0, T₀} : Set G).ncard := Set.ncard_le_ncard hsub (Set.toFinite _)
      _ ≤ ({T₀} : Set G).ncard + 1 := Set.ncard_insert_le _ _
      _ = 2 := by rw [Set.ncard_singleton]
  have hK0 : Nat.card K ≠ 0 := by
    have : Finite K := hKfin.to_subtype
    exact Nat.card_ne_zero.mpr ⟨⟨⟨0, K.zero_mem⟩⟩, inferInstance⟩
  have hidx : (nsmulAddMonoidHom (α := G) 2).range.index = 2 ^ finrank ℤ G * Nat.card K :=
    AddSubgroup.index_range_nsmul_of_fg G (n := 2) two_ne_zero
  have hne : (nsmulAddMonoidHom (α := G) 2).range.index ≠ 0 := by
    rw [hidx]; exact mul_ne_zero (pow_ne_zero _ two_ne_zero) hK0
  have : Finite (G ⧸ (nsmulAddMonoidHom (α := G) 2).range) := Nat.finite_of_card_ne_zero hne
  have h8 := eight_le_card_quotient hind
  rw [← AddSubgroup.index_eq_card, hidx] at h8
  by_contra hr
  have hr' : finrank ℤ G ≤ 1 := by omega
  have hp : 2 ^ finrank ℤ G ≤ 2 := by
    calc 2 ^ finrank ℤ G ≤ 2 ^ 1 := Nat.pow_le_pow_right (by norm_num) hr'
      _ = 2 := by norm_num
  have := Nat.mul_le_mul hp hK2
  omega

end GroupTheory

open WeierstrassCurve WeierstrassCurve.Affine MagicSquareOfSquares.TwoDescent

/-- `E⁺(ℚ)[2] ⊆ {O, T}` on Mathlib's points (`two_torsion_Eplus` handles the affine points). -/
theorem two_torsion_point {m n : ℚ} (hm : m ≠ 0) (hn : n ≠ 0) (hp : m + n ≠ 0) (hd : m - n ≠ 0)
    (hT : (Eplus m n).toAffine.Nonsingular (eP m n) 0) (P : (Eplus m n).toAffine.Point)
    (hP : P + P = 0) : P = 0 ∨ P = Point.some _ _ hT := by
  rcases P with _ | ⟨x, y, h⟩
  · exact Or.inl rfl
  · obtain ⟨rfl, rfl⟩ := two_torsion_Eplus hm hn hp hd (h := h) hP
    exact Or.inr rfl

/-- **Mordell–Weil for `E⁺`** (Stoll's `fg_point_of_numberField`, with `Δ ≠ 0` from
`Eplus_isElliptic`). -/
theorem fg_Eplus {m n : ℚ} (hm : m ≠ 0) (hn : n ≠ 0) (hp : m + n ≠ 0) (hd : m - n ≠ 0) :
    IsAddFG (Eplus m n).toAffine.Point := by
  have := Eplus_isElliptic hm hn hp hd
  exact fg_point_of_numberField

/-- **`rank E⁺(ℚ) ≥ 2` from `Theorem2At`** (independence of `R`, `Q`, `T` in `E(ℚ)/2E(ℚ)`). -/
theorem two_le_rank_of_theorem2At {m n k : ℚ} (hm : m ≠ 0) (hn : n ≠ 0) (hp : m + n ≠ 0)
    (hd : m - n ≠ 0) (h : Theorem2At m n k) :
    2 ≤ finrank ℤ (Eplus m n).toAffine.Point := by
  have := fg_Eplus hm hn hp hd
  obtain ⟨hR, hQ, hT, hind⟩ := h
  exact two_le_finrank_of_indep (two_torsion_point hm hn hp hd hT) hind

/-- **`rank ≥ 2` at `(7, 1)`, both curves** (`E⁺ = Eplus 7 1`, `E⁻ = Eplus 7 (−1)`). -/
theorem rank_7_1 : 2 ≤ finrank ℤ (Eplus 7 1).toAffine.Point ∧
    2 ≤ finrank ℤ (Eplus 7 (-1)).toAffine.Point :=
  ⟨two_le_rank_of_theorem2At (by norm_num) (by norm_num) (by norm_num) (by norm_num) t2_7_1_p,
    two_le_rank_of_theorem2At (by norm_num) (by norm_num) (by norm_num) (by norm_num) t2_7_1_m⟩

/-- **`rank E(ℚ) ≥ 2` on all 254 twice-square curves `m ≤ 1000`**: for each `(m, n, k)` of
`twiceSquareRatios` (127 coprime ratios, `0 < n < m ≤ 1000`, `m² + n² = 2k²`), both `E⁺ = Eplus m n`
and `E⁻ = Eplus m (−n)` have `Module.finrank ℤ E(ℚ) ≥ 2`. -/
theorem rank_twiceSquare : ∀ p ∈ twiceSquareRatios,
    2 ≤ finrank ℤ (Eplus (p.1 : ℚ) (p.2.1 : ℚ)).toAffine.Point ∧
      2 ≤ finrank ℤ (Eplus (p.1 : ℚ) (-(p.2.1 : ℚ))).toAffine.Point := by
  intro p hp
  obtain ⟨hplus, hminus⟩ := theorem2_twiceSquare p hp
  obtain ⟨-, hn0, hnm, -, -⟩ := twiceSquareRatios_ok.2 p hp
  have hn : (0 : ℚ) < (p.2.1 : ℚ) := by exact_mod_cast hn0
  have hnm' : (p.2.1 : ℚ) < (p.1 : ℚ) := by exact_mod_cast hnm
  exact ⟨two_le_rank_of_theorem2At (by intro h; linarith) (by intro h; linarith)
      (by intro h; linarith) (by intro h; linarith) hplus,
    two_le_rank_of_theorem2At (by intro h; linarith) (by intro h; linarith)
      (by intro h; linarith) (by intro h; linarith) hminus⟩

end MagicSquareOfSquares.Rank

#print axioms MagicSquareOfSquares.Rank.two_le_finrank_of_indep
#print axioms MagicSquareOfSquares.Rank.fg_Eplus
#print axioms MagicSquareOfSquares.Rank.two_le_rank_of_theorem2At
#print axioms MagicSquareOfSquares.Rank.rank_7_1
#print axioms MagicSquareOfSquares.Rank.rank_twiceSquare
