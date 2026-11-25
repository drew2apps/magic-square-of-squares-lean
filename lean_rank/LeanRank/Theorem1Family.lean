import LeanRank.Theorem2Family

/-!
# `UNIFORM_RANK.md` Theorem 1 (`rank ≥ 1`) for every ratio, conditional only on Lemma D

Written by an AI model (Claude, made by Anthropic) on 25 November 2026 (UTC); not read by a person.
Experiment, not counted in the census.

For coprime integers `m, n` with `mn ≠ 0`, `m² ≠ n²` (every ratio, not only the twice-square ones):
`α(R) = 2S·H(m, n)`, `α(T) = 4S²·H(m, n)H(m, −n)` and `α(R)α(T) ≡ 2S·H(m, −n)` are not squares (Lemma D twice,
Lemma B), so `R` and `T` are independent in `E⁺(ℚ)/2E⁺(ℚ)`; with `E⁺(ℚ)[2] = {O, T}` and Stoll's Mordell–Weil
theorem, `[E : 2E] = 2^rank·2 ≥ 4`, so `rank E⁺(ℚ) ≥ 1` (`rank_one_of_lemmaD`; both curves: `rank_one_family`).
The only hypothesis beyond the ratio is `LemmaD` (`Theorem2Family.lean`).
-/

namespace MagicSquareOfSquares.Family

open Module WeierstrassCurve WeierstrassCurve.Affine
open MagicSquareOfSquares.G1 MagicSquareOfSquares.TwoDescent MagicSquareOfSquares.Rank
open MagicSquareOfSquares.LemmaB

section GroupTheory

variable {G : Type*} [AddCommGroup G]

/-- The four sums `ε₁R + ε₂T`, `εᵢ ∈ {0, 1}`. -/
def fourSums (R T : G) (c : Fin 2 × Fin 2) : G := (c.1 : ℕ) • R + (c.2 : ℕ) • T

theorem fourSums_injective {R T : G}
    (hind : ∀ c₁ c₂ : ℤ, ∀ W : G, c₁ • R + c₂ • T = W + W → Even c₁ ∧ Even c₂) :
    Function.Injective
      (fun c => (fourSums R T c : G ⧸ (nsmulAddMonoidHom (α := G) 2).range)) := by
  intro c c' h
  have h' : -(fourSums R T c) + fourSums R T c' ∈ (nsmulAddMonoidHom (α := G) 2).range :=
    QuotientAddGroup.eq.mp h
  obtain ⟨W, hW⟩ := AddMonoidHom.mem_range.mp h'
  have hW' : 2 • W = -(fourSums R T c) + fourSums R T c' := hW
  have key : (((c'.1 : ℕ) : ℤ) - ((c.1 : ℕ) : ℤ)) • R + (((c'.2 : ℕ) : ℤ) - ((c.2 : ℕ) : ℤ)) • T
      = W + W := by
    rw [show W + W = 2 • W from (two_nsmul W).symm, hW']
    simp only [fourSums, sub_zsmul, natCast_zsmul]
    abel
  obtain ⟨h1, h2⟩ := hind _ _ W key
  exact Prod.ext (fin2_eq_of_even h1) (fin2_eq_of_even h2)

/-- **Group theory.** In a finitely generated abelian group `G` whose 2-torsion is contained in
`{0, T₀}`, two elements independent modulo `2G` force `finrank ℤ G ≥ 1`. -/
theorem one_le_finrank_of_indep [IsAddFG G] {R T T₀ : G}
    (htors : ∀ P : G, P + P = 0 → P = 0 ∨ P = T₀)
    (hind : ∀ c₁ c₂ : ℤ, ∀ W : G, c₁ • R + c₂ • T = W + W → Even c₁ ∧ Even c₂) :
    1 ≤ finrank ℤ G := by
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
  have h4 : 4 ≤ Nat.card (G ⧸ (nsmulAddMonoidHom (α := G) 2).range) := by
    have h := Nat.card_le_card_of_injective _ (fourSums_injective hind)
    have h4 : Nat.card (Fin 2 × Fin 2) = 4 := by
      rw [Nat.card_eq_fintype_card, Fintype.card_prod, Fintype.card_fin]
    rwa [h4] at h
  rw [← AddSubgroup.index_eq_card, hidx] at h4
  by_contra hr
  have hr' : finrank ℤ G = 0 := by omega
  rw [hr', pow_zero, one_mul] at h4
  omega

end GroupTheory

/-- **Two points independent modulo `2E(K)`** from the 2-descent map: if `α(R)`, `α(T)`, `α(R)α(T)` are not
squares then `c₁R + c₂T = W + W` forces `c₁, c₂` even. -/
theorem indep2_of_alpha {K : Type*} [Field K] [DecidableEq K] {e a b : K} (hq : qe e a b ≠ 0)
    (R T : (curve e a b).toAffine.Point)
    (nR : ¬ IsSquare (alpha e a b R)) (nT : ¬ IsSquare (alpha e a b T))
    (nRT : ¬ IsSquare (alpha e a b R * alpha e a b T)) (c₁ c₂ : ℤ)
    (W : (curve e a b).toAffine.Point) (hW : c₁ • R + c₂ • T = W + W) : Even c₁ ∧ Even c₂ := by
  have key : ∀ c : ℤ, ∃ d : ℤ, c = d + d ∨ c = 1 + d + d := by
    intro c
    rcases Int.even_or_odd' c with ⟨d, hd | hd⟩
    · exact ⟨d, Or.inl (by rw [hd]; ring)⟩
    · exact ⟨d, Or.inr (by rw [hd]; ring)⟩
  obtain ⟨d₁, e₁ | e₁⟩ := key c₁ <;> obtain ⟨d₂, e₂ | e₂⟩ := key c₂ <;> subst e₁ e₂ <;>
    simp only [add_zsmul, one_zsmul] at hW
  · exact ⟨⟨d₁, rfl⟩, ⟨d₂, rfl⟩⟩
  · exfalso; apply nT
    have hT : T = (W - d₁ • R - d₂ • T) + (W - d₁ • R - d₂ • T) := by
      rw [← sub_eq_zero] at hW ⊢; rw [← hW]; abel
    rw [hT]; exact alpha_double hq _
  · exfalso; apply nR
    have hR : R = (W - d₁ • R - d₂ • T) + (W - d₁ • R - d₂ • T) := by
      rw [← sub_eq_zero] at hW ⊢; rw [← hW]; abel
    rw [hR]; exact alpha_double hq _
  · exfalso; apply nRT
    have hRT : R + T = (W - d₁ • R - d₂ • T) + (W - d₁ • R - d₂ • T) := by
      rw [← sub_eq_zero] at hW ⊢; rw [← hW]; abel
    have m := alpha_add hq R T
    rw [mul_assoc] at m
    have s := alpha_double hq (W - d₁ • R - d₂ • T)
    rw [← hRT] at s
    exact isSquare_cancel_left s (alpha_ne_zero hq _) m

/-- **Lemma D, general form**: if `2(m² + n²)·H(m, n)` is a square in `ℚ` and `n ≠ 0`, then `m = ±n`. -/
theorem lemmaD_bridge (hD : LemmaD) {m n : ℚ} (hn : n ≠ 0)
    (h : IsSquare (2 * (m ^ 2 + n ^ 2) * genH m n)) : m = n ∨ m = -n := by
  obtain ⟨z, hz⟩ := h
  obtain ⟨x, rfl⟩ : ∃ x, m = x * n := ⟨m / n, (div_mul_cancel₀ m hn).symm⟩
  have hw : (z / n ^ 3) ^ 2 = 2 * (x ^ 2 + 1) * (x ^ 4 - 2 * x ^ 3 + 2 * x ^ 2 + 2 * x + 1) := by
    rw [div_pow, div_eq_iff (pow_ne_zero 2 (pow_ne_zero 3 hn))]
    unfold genH at hz
    linear_combination -hz
  rcases hD x _ hw with rfl | rfl
  · left; ring
  · right; ring

/-- **`UNIFORM_RANK.md` Theorem 1 (`rank ≥ 1`) for every ratio, conditional only on Lemma D.** For coprime
integers `m, n` with `mn ≠ 0` and `m² ≠ n²`, `R` and `T` are independent in `E⁺(ℚ)/2E⁺(ℚ)` and
`Module.finrank ℤ E⁺(ℚ) ≥ 1`. -/
theorem rank_one_of_lemmaD (hD : LemmaD) {m n : ℤ} (hmn : IsCoprime m n) (hm : m ≠ 0) (hn : n ≠ 0)
    (hne : m ^ 2 ≠ n ^ 2) : 1 ≤ finrank ℤ (Eplus (m : ℚ) (n : ℚ)).toAffine.Point := by
  have hm' : (m : ℚ) ≠ 0 := by exact_mod_cast hm
  have hn' : (n : ℚ) ≠ 0 := by exact_mod_cast hn
  have hmn1 : (m : ℚ) ≠ n := by
    intro h; apply hne; have : m = n := by exact_mod_cast h
    rw [this]
  have hmn2 : (m : ℚ) ≠ -n := by
    intro h; apply hne; have : m = -n := by exact_mod_cast h
    rw [this]; ring
  have hp : (m : ℚ) + n ≠ 0 := by intro h; apply hmn2; linear_combination h
  have hd : (m : ℚ) - n ≠ 0 := by intro h; apply hmn1; linear_combination h
  have hq := qe_ne_zero (n := (n : ℚ)) hm'
  have hS := S_ne_zero (n := (n : ℚ)) hm'
  have hH0 : genH (m : ℚ) n ≠ 0 := (genH_pos hm').ne'
  have hyR : 2 * yR (m : ℚ) n ≠ 0 := by
    unfold yR
    exact mul_ne_zero two_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero
      (mul_ne_zero (by norm_num) hm') hn') hp) hS) hH0)
  have hR : (Eplus (m : ℚ) n).toAffine.Nonsingular (xR m n) (yR m n) := nonsingular_of_y (R_on _ _) hyR
  have hT : (Eplus (m : ℚ) n).toAffine.Nonsingular (eP m n) 0 := nonsingular_T hq
  have hxR : xR (m : ℚ) n ≠ eP m n := by
    intro h; apply mul_ne_zero (mul_ne_zero two_ne_zero hS) hH0; rw [← alphaR_val, h, sub_self]
  have aR : alpha (eP (m : ℚ) n) (aP m n) (bP m n) (Point.some _ _ hR) =
      2 * ((m : ℚ) ^ 2 + (n : ℚ) ^ 2) * genH m n := by
    rw [alpha_some, if_neg hxR, alphaR_val]
  have aT : alpha (eP (m : ℚ) n) (aP m n) (bP m n) (Point.some _ _ hT) = qe (eP m n) (aP m n) (bP m n) := by
    rw [alpha_some, if_pos rfl]
  have nR : ¬ IsSquare (alpha (eP (m : ℚ) n) (aP m n) (bP m n) (Point.some _ _ hR)) := by
    rw [aR]; intro h
    rcases lemmaD_bridge hD hn' h with h1 | h1
    · exact hmn1 h1
    · exact hmn2 h1
  have nT : ¬ IsSquare (alpha (eP (m : ℚ) n) (aP m n) (bP m n) (Point.some _ _ hT)) := by
    rw [aT, qe_P]; intro h
    have := isSquare_cancel_left ⟨2 * ((m : ℚ) ^ 2 + (n : ℚ) ^ 2), by ring⟩
      (mul_ne_zero (by norm_num) (pow_ne_zero 2 hS)) h
    rw [gen_H_mul] at this
    exact lemmaB_rat hmn hm hn hne this
  have nRT : ¬ IsSquare (alpha (eP (m : ℚ) n) (aP m n) (bP m n) (Point.some _ _ hR) *
      alpha (eP (m : ℚ) n) (aP m n) (bP m n) (Point.some _ _ hT)) := by
    rw [aR, aT, qe_P]; intro h
    have h' : IsSquare (2 * ((m : ℚ) ^ 2 + (-(n : ℚ)) ^ 2) * genH m (-n)) := by
      refine isSquare_cancel_left (x := (2 * ((m : ℚ) ^ 2 + (n : ℚ) ^ 2) * genH m n) ^ 2)
        ⟨_, sq _⟩ (pow_ne_zero 2 (mul_ne_zero (mul_ne_zero two_ne_zero hS) hH0)) ?_
      convert h using 1; ring
    rcases lemmaD_bridge hD (neg_ne_zero.mpr hn') h' with h1 | h1
    · exact hmn2 h1
    · exact hmn1 (by rw [h1, neg_neg])
  have hind := indep2_of_alpha hq _ _ nR nT nRT
  have := fg_Eplus hm' hn' hp hd
  exact one_le_finrank_of_indep (two_torsion_point hm' hn' hp hd hT) hind

/-- **`rank E±(ℚ) ≥ 1` for every ratio, conditional only on Lemma D**: `E⁺ = Eplus m n` and
`E⁻ = Eplus m (−n)`. -/
theorem rank_one_family (hD : LemmaD) {m n : ℤ} (hmn : IsCoprime m n) (hm : m ≠ 0) (hn : n ≠ 0)
    (hne : m ^ 2 ≠ n ^ 2) :
    1 ≤ finrank ℤ (Eplus (m : ℚ) (n : ℚ)).toAffine.Point ∧
      1 ≤ finrank ℤ (Eplus (m : ℚ) (-(n : ℚ))).toAffine.Point := by
  refine ⟨rank_one_of_lemmaD hD hmn hm hn hne, ?_⟩
  have := rank_one_of_lemmaD hD hmn.neg_right hm (neg_ne_zero.mpr hn) (by rwa [neg_sq])
  rwa [Int.cast_neg] at this

end MagicSquareOfSquares.Family

#print axioms MagicSquareOfSquares.Family.one_le_finrank_of_indep
#print axioms MagicSquareOfSquares.Family.indep2_of_alpha
#print axioms MagicSquareOfSquares.Family.lemmaD_bridge
#print axioms MagicSquareOfSquares.Family.rank_one_of_lemmaD
#print axioms MagicSquareOfSquares.Family.rank_one_family
