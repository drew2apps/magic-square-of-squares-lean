import LeanRank.Rank

/-!
# `[E(ℚ) : 2E(ℚ)] = 2 ^ (rank + 1)` and `#(E(ℚ)/2E(ℚ)) ≥ 8` for the twice-square curves

Written by an AI model; not read by a person.  Experiment, not counted.
-/

namespace MagicSquareOfSquares.Rank

open Module

section GroupTheory

variable {G : Type*} [AddCommGroup G]

/-- If `G[2] = {0, T₀}` with `T₀ ≠ 0`, then `#G[2] = 2`. -/
theorem card_ker_two_eq {T₀ : G} (hT₀ : T₀ ≠ 0) (hTT : T₀ + T₀ = 0)
    (htors : ∀ P : G, P + P = 0 → P = 0 ∨ P = T₀) :
    Nat.card (nsmulAddMonoidHom (α := G) 2).ker = 2 := by
  have hset : ((nsmulAddMonoidHom (α := G) 2).ker : Set G) = {0, T₀} := by
    ext P
    simp only [SetLike.mem_coe, AddMonoidHom.mem_ker, Set.mem_insert_iff, Set.mem_singleton_iff]
    change 2 • P = 0 ↔ _
    rw [two_nsmul]
    constructor
    · exact htors P
    · rintro (rfl | rfl)
      · exact add_zero 0
      · exact hTT
  have h1 : Nat.card (nsmulAddMonoidHom (α := G) 2).ker =
      ((nsmulAddMonoidHom (α := G) 2).ker : Set G).ncard := Nat.card_coe_set_eq _
  rw [h1, hset, Set.ncard_pair (Ne.symm hT₀)]

/-- **`[G : 2G] = 2 ^ (finrank ℤ G + 1)`** for a finitely generated abelian group with
`G[2] = {0, T₀}`, `T₀ ≠ 0` (Stoll's `AddSubgroup.index_range_nsmul_of_fg`). -/
theorem index_two_eq [IsAddFG G] {T₀ : G} (hT₀ : T₀ ≠ 0) (hTT : T₀ + T₀ = 0)
    (htors : ∀ P : G, P + P = 0 → P = 0 ∨ P = T₀) :
    (nsmulAddMonoidHom (α := G) 2).range.index = 2 ^ (finrank ℤ G + 1) := by
  rw [AddSubgroup.index_range_nsmul_of_fg G two_ne_zero, card_ker_two_eq hT₀ hTT htors, pow_succ]

end GroupTheory

open WeierstrassCurve WeierstrassCurve.Affine MagicSquareOfSquares.TwoDescent

/-- `T = (e, 0)` has order 2. -/
theorem T_add_T {m n : ℚ} (hT : (Eplus m n).toAffine.Nonsingular (eP m n) 0) :
    Point.some _ _ hT + Point.some _ _ hT = 0 :=
  Point.add_self_of_Y_eq (by rw [curve_negY, neg_zero])

/-- **`[E⁺(ℚ) : 2E⁺(ℚ)] = 2 ^ (rank + 1)`**, from `E⁺(ℚ)[2] = {O, T}`, at every ratio with
`m, n, m ± n ≠ 0`. -/
theorem index_Eplus {m n : ℚ} (hm : m ≠ 0) (hn : n ≠ 0) (hp : m + n ≠ 0) (hd : m - n ≠ 0)
    (hT : (Eplus m n).toAffine.Nonsingular (eP m n) 0) :
    (nsmulAddMonoidHom (α := (Eplus m n).toAffine.Point) 2).range.index =
      2 ^ (finrank ℤ (Eplus m n).toAffine.Point + 1) := by
  have := fg_Eplus hm hn hp hd
  exact index_two_eq (Point.some_ne_zero hT) (T_add_T hT) (two_torsion_point hm hn hp hd hT)

/-- **`#(E⁺(ℚ)/2E⁺(ℚ)) ≥ 8` from `Theorem2At`.** -/
theorem eight_le_index_of_theorem2At {m n k : ℚ} (hm : m ≠ 0) (hn : n ≠ 0) (hp : m + n ≠ 0)
    (hd : m - n ≠ 0) (h : Theorem2At m n k) :
    8 ≤ (nsmulAddMonoidHom (α := (Eplus m n).toAffine.Point) 2).range.index := by
  obtain ⟨hR, hQ, hT, hind⟩ := h
  have hidx := index_Eplus hm hn hp hd hT
  have : Finite ((Eplus m n).toAffine.Point ⧸
      (nsmulAddMonoidHom (α := (Eplus m n).toAffine.Point) 2).range) :=
    Nat.finite_of_card_ne_zero (by rw [← AddSubgroup.index_eq_card, hidx]; positivity)
  rw [AddSubgroup.index_eq_card]
  exact eight_le_card_quotient hind

/-- **`#(E(ℚ)/2E(ℚ)) ≥ 8` on the same 254 curves** (and `[E : 2E] = 2 ^ (rank + 1)` by
`index_Eplus`). -/
theorem index_twiceSquare : ∀ p ∈ twiceSquareRatios,
    8 ≤ (nsmulAddMonoidHom (α := (Eplus (p.1 : ℚ) (p.2.1 : ℚ)).toAffine.Point) 2).range.index ∧
      8 ≤ (nsmulAddMonoidHom (α := (Eplus (p.1 : ℚ) (-(p.2.1 : ℚ))).toAffine.Point) 2).range.index := by
  intro p hp
  obtain ⟨hplus, hminus⟩ := theorem2_twiceSquare p hp
  obtain ⟨-, hn0, hnm, -, -⟩ := twiceSquareRatios_ok.2 p hp
  have hn : (0 : ℚ) < (p.2.1 : ℚ) := by exact_mod_cast hn0
  have hnm' : (p.2.1 : ℚ) < (p.1 : ℚ) := by exact_mod_cast hnm
  exact ⟨eight_le_index_of_theorem2At (by intro h; linarith) (by intro h; linarith)
      (by intro h; linarith) (by intro h; linarith) hplus,
    eight_le_index_of_theorem2At (by intro h; linarith) (by intro h; linarith)
      (by intro h; linarith) (by intro h; linarith) hminus⟩

end MagicSquareOfSquares.Rank

#print axioms MagicSquareOfSquares.Rank.index_two_eq
#print axioms MagicSquareOfSquares.Rank.index_Eplus
#print axioms MagicSquareOfSquares.Rank.eight_le_index_of_theorem2At
#print axioms MagicSquareOfSquares.Rank.index_twiceSquare
