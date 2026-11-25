import LeanRank.LemmaB
import LeanRank.Rank

/-!
# `UNIFORM_RANK.md` Theorem 2 for every twice-square ratio, conditional only on Lemma D

Written by an AI model (Claude, made by Anthropic) on 25 November 2026 (UTC); not read by a person.
Experiment, not counted in the census.

Of the five non-square hypotheses of `theorem2_of` (`TwoDescentMap.lean`), three are proved here for every
twice-square ratio: `H(m, n)·H(m, −n)` (Lemma B, `LemmaB.lean`) and `±2mn` (parity: `m, n` are odd). The other
two, `H(m, ±n)`, are Lemma D of `UNIFORM_RANK.md`, kept as a hypothesis:

* `LemmaD`: every rational point `(x, w)` of `D2 : w² = 2(x² + 1)(x⁴ − 2x³ + 2x² + 2x + 1)` has `x = ±1`.
  (`UNIFORM_RANK.md` section 4; LMFDB genus-2 curve 576.a.576.1; not proved in Lean.)

Results.
* `theorem2_family_of_H`: for coprime integers `1 ≤ n < m` with `m² + n² = 2k²`, if `H(m, n)` and `H(m, −n)` are
  not squares in `ℚ`, then `R, Q, T` are independent in `E(ℚ)/2E(ℚ)` on both `E⁺` and `E⁻`
  (`Theorem2At`), and `rank E±(ℚ) ≥ 2` (`rank_family_of_H`, Stoll's Mordell–Weil theorem).
* `theorem2_family`, `rank_family`: the same with the single hypothesis `LemmaD` (for all ratios at once).
* Unconditional and uniform (no Lemma D): `T_not_double` (`T ∉ 2E⁺(ℚ)`) and `no_order_four` (`E⁺(ℚ)` has no
  point of order 4) for every coprime `m, n` with `mn ≠ 0`, `m² ≠ n²`; with `E⁻ = Eplus m (−n)` this is the
  torsion part of Theorem 1 of `UNIFORM_RANK.md` for every ratio.
-/

namespace MagicSquareOfSquares.Family

open Module WeierstrassCurve WeierstrassCurve.Affine
open MagicSquareOfSquares.G1 MagicSquareOfSquares.TwoDescent MagicSquareOfSquares.Rank
open MagicSquareOfSquares.LemmaB

/-- **Lemma D of `UNIFORM_RANK.md`** (as a hypothesis): the affine rational points of the genus-2 curve
`D2 : w² = 2(x² + 1)(x⁴ − 2x³ + 2x² + 2x + 1)` all have `x = ±1`. -/
def LemmaD : Prop :=
  ∀ x w : ℚ, w ^ 2 = 2 * (x ^ 2 + 1) * (x ^ 4 - 2 * x ^ 3 + 2 * x ^ 2 + 2 * x + 1) → x = 1 ∨ x = -1

/-- `D2` does have the points `(±1, ±4)`, so `LemmaD` is the strongest statement of its kind. -/
theorem D2_points (x w : ℚ) (hx : x = 1 ∨ x = -1) (hw : w = 4 ∨ w = -4) :
    w ^ 2 = 2 * (x ^ 2 + 1) * (x ^ 4 - 2 * x ^ 3 + 2 * x ^ 2 + 2 * x + 1) := by
  rcases hx with rfl | rfl <;> rcases hw with rfl | rfl <;> norm_num

/-- **Lemma D at a twice-square ratio**: if `m² + n² = 2k²`, `n ≠ 0`, `m ≠ ±n`, then `H(m, n)` is not a square
(else `(m/n, 2kh/n³)` is a point of `D2` with `m/n ≠ ±1`). -/
theorem not_isSquare_genH (hD : LemmaD) {m n k : ℚ} (hn : n ≠ 0) (hk : m ^ 2 + n ^ 2 = 2 * k ^ 2)
    (h1 : m ≠ n) (h2 : m ≠ -n) : ¬ IsSquare (genH m n) := by
  rintro ⟨h, hh⟩
  obtain ⟨x, rfl⟩ : ∃ x, m = x * n := ⟨m / n, (div_mul_cancel₀ m hn).symm⟩
  have hw : (2 * k * h / n ^ 3) ^ 2 =
      2 * (x ^ 2 + 1) * (x ^ 4 - 2 * x ^ 3 + 2 * x ^ 2 + 2 * x + 1) := by
    rw [div_pow, div_eq_iff (pow_ne_zero 2 (pow_ne_zero 3 hn))]
    unfold genH at hh
    linear_combination (-2 * h ^ 2) * hk + (-2 * n ^ 2 * (x ^ 2 + 1)) * hh
  rcases hD x _ hw with rfl | rfl
  · exact h1 (by ring)
  · exact h2 (by ring)

/-- At a twice-square ratio with `m, n` coprime, `m` and `n` are odd. -/
theorem odd_of_twice_square {m n k : ℤ} (hmn : IsCoprime m n) (hk : m ^ 2 + n ^ 2 = 2 * k ^ 2) :
    Odd m ∧ Odd n := by
  rcases Int.even_or_odd m with ⟨i, rfl⟩ | hm <;> rcases Int.even_or_odd n with ⟨j, rfl⟩ | hn
  · obtain ⟨u, v, huv⟩ := hmn
    have : (2 : ℤ) * (u * i + v * j) = 1 := by linear_combination huv
    generalize u * i + v * j = X at this
    omega
  · exfalso
    obtain ⟨j, rfl⟩ := hn
    have : (2 : ℤ) * (2 * i ^ 2 + 2 * j ^ 2 + 2 * j - k ^ 2) = -1 := by linear_combination hk
    generalize 2 * i ^ 2 + 2 * j ^ 2 + 2 * j - k ^ 2 = X at this
    omega
  · exfalso
    obtain ⟨i, rfl⟩ := hm
    have : (2 : ℤ) * (2 * i ^ 2 + 2 * i + 2 * j ^ 2 - k ^ 2) = -1 := by linear_combination hk
    generalize 2 * i ^ 2 + 2 * i + 2 * j ^ 2 - k ^ 2 = X at this
    omega
  · exact ⟨hm, hn⟩

/-- **The parity fact**: for odd `m, n`, `2mn` is not a square in `ℚ`. -/
theorem not_isSquare_two_mul {m n : ℤ} (hm : Odd m) (hn : Odd n) : ¬ IsSquare (2 * (m : ℚ) * (n : ℚ)) := by
  rw [show (2 : ℚ) * m * n = ((2 * m * n : ℤ) : ℚ) by push_cast; ring, Rat.isSquare_intCast_iff]
  rintro ⟨z, hz⟩
  have hze : Even (z * z) := ⟨m * n, by rw [← hz]; ring⟩
  obtain ⟨w, rfl⟩ : Even z := (Int.even_mul.mp hze).elim id id
  have h2 : m * n = 2 * (w * w) :=
    mul_left_cancel₀ (two_ne_zero : (2 : ℤ) ≠ 0) (by linear_combination hz)
  exact (Int.not_even_iff_odd.mpr (hm.mul hn)) ⟨w * w, by rw [h2]; ring⟩

/-- **Theorem 2 for every twice-square ratio, from `H(m, ±n)` not squares.** For coprime integers
`1 ≤ n < m` with `m² + n² = 2k²`: if `H(m, n)` and `H(m, −n)` are not squares in `ℚ`, then `R, Q, T` are
independent in `E(ℚ)/2E(ℚ)` on `E⁺ = Eplus m n` and on `E⁻ = Eplus m (−n)`. Lemma B and the parity fact are
proved, not assumed. -/
theorem theorem2_family_of_H {m n k : ℤ} (hmn : IsCoprime m n) (hn : 1 ≤ n) (hnm : n < m)
    (hk : m ^ 2 + n ^ 2 = 2 * k ^ 2) (hH : ¬ IsSquare (genH (m : ℚ) (n : ℚ)))
    (hH' : ¬ IsSquare (genH (m : ℚ) (-(n : ℚ)))) :
    Theorem2At (m : ℚ) (n : ℚ) (k : ℚ) ∧ Theorem2At (m : ℚ) (-(n : ℚ)) (k : ℚ) := by
  obtain ⟨hmo, hno⟩ := odd_of_twice_square hmn hk
  have hn' : (1 : ℚ) ≤ n := by exact_mod_cast hn
  have hnm' : (n : ℚ) < m := by exact_mod_cast hnm
  have hkQ : (m : ℚ) ^ 2 + (n : ℚ) ^ 2 = 2 * (k : ℚ) ^ 2 := by exact_mod_cast hk
  have hkQ' : (m : ℚ) ^ 2 + (-(n : ℚ)) ^ 2 = 2 * (k : ℚ) ^ 2 := by rw [neg_sq]; exact hkQ
  have hm0 : (m : ℚ) ≠ 0 := by intro h; linarith
  have hn0 : (n : ℚ) ≠ 0 := by intro h; linarith
  have hne : m ^ 2 ≠ n ^ 2 := by intro h; nlinarith
  have hHH : ¬ IsSquare (genH (m : ℚ) n * genH m (-n)) := by
    rw [gen_H_mul]
    exact lemmaB_rat hmn (by exact_mod_cast hm0) (by exact_mod_cast hn0) hne
  have h2 : ¬ IsSquare (2 * (m : ℚ) * n) := not_isSquare_two_mul hmo hno
  have h2' : ¬ IsSquare (2 * (m : ℚ) * -(n : ℚ)) := by
    have := not_isSquare_two_mul hmo hno.neg
    rwa [Int.cast_neg] at this
  refine ⟨theorem2_of _ _ _ hkQ hm0 hn0 (by intro h; linarith) (by intro h; linarith) hH hH' hHH h2
    (by rw [show (-2 : ℚ) * m * n = 2 * m * -n by ring]; exact h2'), ?_⟩
  refine theorem2_of _ _ _ hkQ' hm0 (neg_ne_zero.mpr hn0) (by intro h; linarith) (by intro h; linarith)
    hH' (by rwa [neg_neg]) (by rwa [neg_neg, mul_comm]) h2'
    (by rw [show (-2 : ℚ) * m * -n = 2 * m * n by ring]; exact h2)

/-- **`rank E±(ℚ) ≥ 2` for every twice-square ratio, from `H(m, ±n)` not squares** (Stoll's Mordell–Weil
theorem through `two_le_rank_of_theorem2At`). -/
theorem rank_family_of_H {m n k : ℤ} (hmn : IsCoprime m n) (hn : 1 ≤ n) (hnm : n < m)
    (hk : m ^ 2 + n ^ 2 = 2 * k ^ 2) (hH : ¬ IsSquare (genH (m : ℚ) (n : ℚ)))
    (hH' : ¬ IsSquare (genH (m : ℚ) (-(n : ℚ)))) :
    2 ≤ finrank ℤ (Eplus (m : ℚ) (n : ℚ)).toAffine.Point ∧
      2 ≤ finrank ℤ (Eplus (m : ℚ) (-(n : ℚ))).toAffine.Point := by
  obtain ⟨h1, h2⟩ := theorem2_family_of_H hmn hn hnm hk hH hH'
  have hn' : (1 : ℚ) ≤ n := by exact_mod_cast hn
  have hnm' : (n : ℚ) < m := by exact_mod_cast hnm
  exact ⟨two_le_rank_of_theorem2At (by intro h; linarith) (by intro h; linarith) (by intro h; linarith)
      (by intro h; linarith) h1,
    two_le_rank_of_theorem2At (by intro h; linarith) (by intro h; linarith) (by intro h; linarith)
      (by intro h; linarith) h2⟩

/-- **Theorem 2 of `UNIFORM_RANK.md` for every twice-square ratio, conditional only on Lemma D.** -/
theorem theorem2_family (hD : LemmaD) {m n k : ℤ} (hmn : IsCoprime m n) (hn : 1 ≤ n) (hnm : n < m)
    (hk : m ^ 2 + n ^ 2 = 2 * k ^ 2) :
    Theorem2At (m : ℚ) (n : ℚ) (k : ℚ) ∧ Theorem2At (m : ℚ) (-(n : ℚ)) (k : ℚ) := by
  have hn' : (1 : ℚ) ≤ n := by exact_mod_cast hn
  have hnm' : (n : ℚ) < m := by exact_mod_cast hnm
  have hkQ : (m : ℚ) ^ 2 + (n : ℚ) ^ 2 = 2 * (k : ℚ) ^ 2 := by exact_mod_cast hk
  have hkQ' : (m : ℚ) ^ 2 + (-(n : ℚ)) ^ 2 = 2 * (k : ℚ) ^ 2 := by rw [neg_sq]; exact hkQ
  have hn0 : (n : ℚ) ≠ 0 := by intro h; linarith
  exact theorem2_family_of_H hmn hn hnm hk
    (not_isSquare_genH hD hn0 hkQ (by intro h; linarith) (by intro h; linarith))
    (not_isSquare_genH hD (neg_ne_zero.mpr hn0) hkQ' (by intro h; linarith)
      (by intro h; rw [neg_neg] at h; linarith))

/-- **`rank E±(ℚ) ≥ 2` for every twice-square ratio, conditional only on Lemma D.** For coprime integers
`1 ≤ n < m` with `m² + n² = 2k²`, both `E⁺ = Eplus m n` and `E⁻ = Eplus m (−n)` have
`Module.finrank ℤ E(ℚ) ≥ 2`. -/
theorem rank_family (hD : LemmaD) {m n k : ℤ} (hmn : IsCoprime m n) (hn : 1 ≤ n) (hnm : n < m)
    (hk : m ^ 2 + n ^ 2 = 2 * k ^ 2) :
    2 ≤ finrank ℤ (Eplus (m : ℚ) (n : ℚ)).toAffine.Point ∧
      2 ≤ finrank ℤ (Eplus (m : ℚ) (-(n : ℚ))).toAffine.Point := by
  obtain ⟨h1, h2⟩ := theorem2_family hD hmn hn hnm hk
  have hn' : (1 : ℚ) ≤ n := by exact_mod_cast hn
  have hnm' : (n : ℚ) < m := by exact_mod_cast hnm
  exact ⟨two_le_rank_of_theorem2At (by intro h; linarith) (by intro h; linarith) (by intro h; linarith)
      (by intro h; linarith) h1,
    two_le_rank_of_theorem2At (by intro h; linarith) (by intro h; linarith) (by intro h; linarith)
      (by intro h; linarith) h2⟩

/-! ## Unconditional and uniform: the torsion part of Theorem 1 -/

theorem S_ne_zero {m n : ℚ} (hm : m ≠ 0) : m ^ 2 + n ^ 2 ≠ 0 := by
  have h1 : 0 < m ^ 2 := lt_of_le_of_ne (sq_nonneg m) (Ne.symm (pow_ne_zero 2 hm))
  have h2 := sq_nonneg n
  exact (by linarith : 0 < m ^ 2 + n ^ 2).ne'

theorem qe_ne_zero {m n : ℚ} (hm : m ≠ 0) : qe (eP m n) (aP m n) (bP m n) ≠ 0 := by
  rw [qe_P]
  exact mul_ne_zero (mul_ne_zero (by norm_num) (pow_ne_zero 2 (S_ne_zero hm)))
    (mul_ne_zero (genH_pos hm).ne' (genH_pos hm).ne')

/-- **`T ∉ 2E⁺(ℚ)` for every ratio** (Lemma B): for coprime integers `m, n` with `mn ≠ 0`, `m² ≠ n²`, no
rational point `W` of `E⁺` has `W + W = T = (−2S³, 0)`. (For `E⁻` apply it to `(m, −n)`.) -/
theorem T_not_double {m n : ℤ} (hmn : IsCoprime m n) (hm : m ≠ 0) (hn : n ≠ 0) (hne : m ^ 2 ≠ n ^ 2)
    (hT : (Eplus (m : ℚ) (n : ℚ)).toAffine.Nonsingular (eP m n) 0)
    (W : (Eplus (m : ℚ) (n : ℚ)).toAffine.Point) : W + W ≠ Point.some _ _ hT := by
  intro h
  have hm' : (m : ℚ) ≠ 0 := by exact_mod_cast hm
  have hsq := double_T_isSquare (qe_ne_zero (n := (n : ℚ)) hm') h
  rw [qe_P] at hsq
  have := isSquare_cancel_left ⟨2 * ((m : ℚ) ^ 2 + (n : ℚ) ^ 2), by ring⟩
    (mul_ne_zero (by norm_num) (pow_ne_zero 2 (S_ne_zero (n := (n : ℚ)) hm'))) hsq
  rw [gen_H_mul] at this
  exact lemmaB_rat hmn hm hn hne this

/-- **No point of order 4, for every ratio**: if `4W = O` on `E⁺(ℚ)` then `2W = O` (`E⁺(ℚ)[2] = {O, T}` and
`T ∉ 2E⁺(ℚ)`). -/
theorem no_order_four {m n : ℤ} (hmn : IsCoprime m n) (hm : m ≠ 0) (hn : n ≠ 0) (hne : m ^ 2 ≠ n ^ 2)
    (W : (Eplus (m : ℚ) (n : ℚ)).toAffine.Point) (h4 : (W + W) + (W + W) = 0) : W + W = 0 := by
  have hm' : (m : ℚ) ≠ 0 := by exact_mod_cast hm
  have hn' : (n : ℚ) ≠ 0 := by exact_mod_cast hn
  have hp : (m : ℚ) + n ≠ 0 := by
    intro h; apply hne
    have : ((m + n : ℤ) : ℚ) = 0 := by push_cast; exact h
    have h' : m + n = 0 := by exact_mod_cast this
    linear_combination (m - n) * h'
  have hd : (m : ℚ) - n ≠ 0 := by
    intro h; apply hne
    have : ((m - n : ℤ) : ℚ) = 0 := by push_cast; exact h
    have h' : m - n = 0 := by exact_mod_cast this
    linear_combination (m + n) * h'
  have hT := nonsingular_T (qe_ne_zero (n := (n : ℚ)) hm')
  rcases two_torsion_point hm' hn' hp hd hT (W + W) h4 with h | h
  · exact h
  · exact absurd h (T_not_double hmn hm hn hne hT W)

/-! ## Negative controls -/

/-- Lemma B needs `m² ≠ n²`: at `m = n = 1`, `m⁸ + 14m⁴n⁴ + n⁸ = 16 = 4²`. -/
theorem control_lemmaB_trivial : IsSquare ((1 : ℤ) ^ 8 + 14 * 1 ^ 4 * 1 ^ 4 + 1 ^ 8) := ⟨4, by norm_num⟩

/-- So the statement of Lemma B without `m² ≠ n²` is false (and Lean proves its negation). -/
theorem control_lemmaB_false :
    ¬ ∀ m n : ℤ, IsCoprime m n → m ≠ 0 → n ≠ 0 → ¬ IsSquare (m ^ 8 + 14 * m ^ 4 * n ^ 4 + n ^ 8) := by
  intro h
  exact h 1 1 isCoprime_one_left one_ne_zero one_ne_zero control_lemmaB_trivial

/-- Pocklington's equation needs `a + b` odd: `1⁴ − 1²·1² + 1⁴ = 1²`. -/
theorem control_pocklington_parity : (1 : ℤ) ^ 4 - 1 ^ 2 * 1 ^ 2 + 1 ^ 4 = 1 ^ 2 := by norm_num

/-- The hypothesis `¬ IsSquare H(m, n)` has content: off the twice-square type it can fail,
`H(4, 1) = 169 = 13²` (`y² = H(x, 1)` has rank 1, `UNIFORM_RANK.md` section 6). -/
theorem control_H_4_1 : IsSquare (genH 4 1) := ⟨13, by unfold genH; norm_num⟩

end MagicSquareOfSquares.Family

#print axioms MagicSquareOfSquares.Family.not_isSquare_genH
#print axioms MagicSquareOfSquares.Family.not_isSquare_two_mul
#print axioms MagicSquareOfSquares.Family.theorem2_family_of_H
#print axioms MagicSquareOfSquares.Family.rank_family_of_H
#print axioms MagicSquareOfSquares.Family.theorem2_family
#print axioms MagicSquareOfSquares.Family.rank_family
#print axioms MagicSquareOfSquares.Family.T_not_double
#print axioms MagicSquareOfSquares.Family.no_order_four
#print axioms MagicSquareOfSquares.Family.control_lemmaB_false
