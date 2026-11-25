import LeanRank.Theorem1Family

/-!
# `UNIFORM_RANK.md` Theorem 1 in full, for every ratio: `R` (and `P = 2R`) has infinite order

Written by an AI model (Claude, made by Anthropic) on 25 November 2026 (UTC); not read by a person.
Experiment, not counted in the census.

`R_infinite_order` (hypothesis `LemmaD` only): for coprime integers `m, n` with `mn ≠ 0`, `m² ≠ n²`, no
multiple `N•R` with `N ≥ 1` is `O` on `E⁺ = Eplus m n`. Proof: write `N = 2^a·d` with `d` odd; `V = d•R` has
`2^a•V = O`, so `2V = O` (`no_order_four`, by induction on `a`), so `V ∈ {O, T}`; and `R − V ∈ 2E(ℚ)`, so
`α(R)` or `α(R)α(T)` would be a square, which Lemma D excludes. With `R_add_R` (`R + R = P`) this gives that
`P` has infinite order too (`P_infinite_order`). Together with `two_torsion_point` (`E⁺(ℚ)[2] = {O, T}`),
`no_order_four` and `rank_one_family` this is all of Theorem 1, for every ratio, conditional on Lemma D.
-/

namespace MagicSquareOfSquares.Family

open Module WeierstrassCurve WeierstrassCurve.Affine
open MagicSquareOfSquares.G1 MagicSquareOfSquares.TwoDescent MagicSquareOfSquares.Rank
open MagicSquareOfSquares.LemmaB

/-- `2^a•W = O` implies `2W = O` on `E⁺(ℚ)` (no point of order 4). -/
theorem two_pow_torsion {m n : ℤ} (hmn : IsCoprime m n) (hm : m ≠ 0) (hn : n ≠ 0) (hne : m ^ 2 ≠ n ^ 2)
    (a : ℕ) : ∀ W : (Eplus (m : ℚ) (n : ℚ)).toAffine.Point, (2 ^ a) • W = 0 → W + W = 0 := by
  induction a with
  | zero =>
    intro W h
    rw [pow_zero, one_smul] at h
    rw [h, add_zero]
  | succ a ih =>
    intro W h
    rw [pow_succ, mul_smul, two_nsmul] at h
    exact no_order_four hmn hm hn hne W (ih (W + W) h)

/-- **`R` has infinite order, for every ratio** (hypothesis: Lemma D). -/
theorem R_infinite_order (hD : LemmaD) {m n : ℤ} (hmn : IsCoprime m n) (hm : m ≠ 0) (hn : n ≠ 0)
    (hne : m ^ 2 ≠ n ^ 2) (hR : (Eplus (m : ℚ) (n : ℚ)).toAffine.Nonsingular (xR m n) (yR m n))
    (N : ℕ) (hN : 0 < N) : N • Point.some _ _ hR ≠ 0 := by
  intro h0
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
  have hT : (Eplus (m : ℚ) n).toAffine.Nonsingular (eP m n) 0 := nonsingular_T hq
  have hxR : xR (m : ℚ) n ≠ eP m n := by
    intro h; apply mul_ne_zero (mul_ne_zero two_ne_zero hS) hH0; rw [← alphaR_val, h, sub_self]
  have aR : alpha (eP (m : ℚ) n) (aP m n) (bP m n) (Point.some _ _ hR) =
      2 * ((m : ℚ) ^ 2 + (n : ℚ) ^ 2) * genH m n := by
    rw [alpha_some, if_neg hxR, alphaR_val]
  have aT : alpha (eP (m : ℚ) n) (aP m n) (bP m n) (Point.some _ _ hT) =
      qe (eP m n) (aP m n) (bP m n) := by
    rw [alpha_some, if_pos rfl]
  -- `N = 2^a·d`, `d` odd
  obtain ⟨a, d, hdo, rfl⟩ := Nat.exists_eq_two_pow_mul_odd hN.ne'
  set R := Point.some _ _ hR with hRdef
  have hV : (2 ^ a) • (d • R) = 0 := by rw [← mul_smul]; exact h0
  have hVV := two_pow_torsion hmn hm hn hne a (d • R) hV
  obtain ⟨e, rfl⟩ := hdo
  have hdR : (2 * e + 1) • R = (e • R + e • R) + R := by
    rw [add_smul, one_smul, mul_comm, mul_smul, two_nsmul, nsmul_add]
  rcases two_torsion_point hm' hn' hp hd hT _ hVV with hV0 | hVT
  · -- `R = W + W` with `W = −(e•R)`
    have hR2 : R = (-(e • R)) + (-(e • R)) := by
      rw [hdR] at hV0
      rw [← sub_eq_zero, ← hV0]; abel
    have := alpha_double hq (-(e • R))
    rw [← hR2, aR] at this
    rcases lemmaD_bridge hD hn' this with h1 | h1
    · exact hmn1 h1
    · exact hmn2 h1
  · -- `R + T = W + W`
    have hTT : Point.some _ _ hT + Point.some _ _ hT = 0 := by rw [← hVT]; exact hVV
    have hRT : R + Point.some _ _ hT = (-(e • R)) + (-(e • R)) := by
      rw [hdR] at hVT
      have e1 : R + Point.some _ _ hT =
          (e • R + e • R + R) + Point.some _ _ hT - (e • R + e • R) := by abel
      rw [e1, hVT, hTT]; abel
    have s := alpha_double hq (-(e • R))
    rw [← hRT] at s
    have m₁ := alpha_add hq R (Point.some _ _ hT)
    rw [mul_assoc] at m₁
    have h' := isSquare_cancel_left s (alpha_ne_zero hq _) m₁
    rw [aR, aT, qe_P] at h'
    have h'' : IsSquare (2 * ((m : ℚ) ^ 2 + (-(n : ℚ)) ^ 2) * genH m (-n)) := by
      refine isSquare_cancel_left (x := (2 * ((m : ℚ) ^ 2 + (n : ℚ) ^ 2) * genH m n) ^ 2)
        ⟨_, sq _⟩ (pow_ne_zero 2 (mul_ne_zero (mul_ne_zero two_ne_zero hS) hH0)) ?_
      convert h' using 1; ring
    rcases lemmaD_bridge hD (neg_ne_zero.mpr hn') h'' with h1 | h1
    · exact hmn2 h1
    · exact hmn1 (by rw [h1, neg_neg])

/-- **`P = 2R` (the point of `z = 2`) has infinite order, for every ratio** (hypothesis: Lemma D). -/
theorem P_infinite_order (hD : LemmaD) {m n : ℤ} (hmn : IsCoprime m n) (hm : m ≠ 0) (hn : n ≠ 0)
    (hne : m ^ 2 ≠ n ^ 2) :
    ∃ hP : (Eplus (m : ℚ) (n : ℚ)).toAffine.Nonsingular (xP m n) (yP m n),
      ∀ N : ℕ, 0 < N → N • Point.some _ _ hP ≠ 0 := by
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
  obtain ⟨hR, hP, hRR⟩ := R_add_R hm' hn' hp hd
  refine ⟨hP, fun N hN h0 => R_infinite_order hD hmn hm hn hne hR (2 * N) (by omega) ?_⟩
  rw [mul_comm, mul_smul, two_nsmul, hRR]; exact h0

end MagicSquareOfSquares.Family

#print axioms MagicSquareOfSquares.Family.two_pow_torsion
#print axioms MagicSquareOfSquares.Family.R_infinite_order
#print axioms MagicSquareOfSquares.Family.P_infinite_order
