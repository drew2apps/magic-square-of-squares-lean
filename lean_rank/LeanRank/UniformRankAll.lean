import LeanRank.EminusModel

/-!
# Theorems 1 and 2 of `UNIFORM_RANK.md` for all rational `m, n` (no coprimality)

Written by an AI model; not read by a person.

`UniformRank.lean` proves Theorems 1 and 2 for coprime integer pairs. Coprimality enters those proofs at two
places only: Lemma B (`LemmaB.lemmaB`, proved by descent on coprime integers) and, in Theorem 2, the parity fact
"`m, n` odd, so `±2mn` is not a square". Lemma D (`lemmaD_bridge`, `not_isSquare_genH`), the 2-descent map, the
torsion and Mordell–Weil steps are already stated over `ℚ`. Both facts are about homogeneous forms of even degree,
so they pass to every rational pair by scaling the numbers (not the curves): if `m = c·p`, `n = c·q` with `p, q`
coprime integers and `c ∈ ℚˣ`, then `m⁸ + 14m⁴n⁴ + n⁸ = (c⁴)²·(p⁸ + 14p⁴q⁴ + q⁸)` and `2mn = c²·2pq`, and
`m² + n² = 2k²` gives `p² + q² = 2u²` for an integer `u` (`exists_coprime_scale`, `lemmaB_rat_all`,
`not_isSquare_two_mul_all`). No isomorphism of curves is needed.

Results, for all `m, n : ℚ` (hence all integers):
* `theorem1_all`: `mn(m² − n²) ≠ 0` ⇒ `rank E⁺(ℚ) ≥ 1` and `rank E⁻(ℚ) ≥ 1` (`E⁺ = Eplus m n`,
  `E⁻ = Eplus m (−n) = Eminus m n`), and `theorem1_P_all`: `P` has infinite order on `E⁺` (on `E⁻` apply it to
  `(m, −n)`); also `T_not_double_all`, `no_order_four_all`;
* `theorem2_all`: `m² + n² = 2k²`, `m² ≠ n²` ⇒ `R, Q, T` independent in `E(ℚ)/2E(ℚ)` on both curves
  (`Theorem2At`), and both ranks `≥ 2`;
* `theorem1_all_Eminus`, `theorem2_all_Eminus`: the same rank statements with the curves of the recipe of §1
  (`EminusModel.lean`).
-/

namespace MagicSquareOfSquares.UniformRankAll

open Module WeierstrassCurve WeierstrassCurve.Affine
open MagicSquareOfSquares.G1 MagicSquareOfSquares.TwoDescent MagicSquareOfSquares.Rank
open MagicSquareOfSquares.Family MagicSquareOfSquares.EminusModel

/-! ## Scaling the numbers -/

/-- **Every rational pair with `n ≠ 0` is a non-zero rational multiple of a coprime integer pair**:
`m = c·p`, `n = c·q` with `p = num(m/n)`, `q = den(m/n)`, `c = n/q`. -/
theorem exists_coprime_scale {m n : ℚ} (hn : n ≠ 0) :
    ∃ (p q : ℤ) (c : ℚ), IsCoprime p q ∧ q ≠ 0 ∧ c ≠ 0 ∧ m = c * p ∧ n = c * q := by
  have hd : ((m / n).den : ℚ) ≠ 0 := by exact_mod_cast (m / n).den_nz
  refine ⟨(m / n).num, (m / n).den, n / (m / n).den, ?_, ?_, div_ne_zero hn hd, ?_, ?_⟩
  · rw [Int.isCoprime_iff_gcd_eq_one]
    exact (m / n).reduced
  · exact_mod_cast (m / n).den_nz
  · have h := Rat.num_div_den (m / n)
    calc m = (m / n) * n := (div_mul_cancel₀ m hn).symm
      _ = ((m / n).num / (m / n).den) * n := by rw [h]
      _ = n / (m / n).den * (m / n).num := by ring
  · push_cast
    field_simp

/-- `2` is not a square in `ℚ` (from the parity fact at `m = n = 1`). -/
theorem not_isSquare_two : ¬ IsSquare (2 : ℚ) := by
  have := not_isSquare_two_mul (m := 1) (n := 1) odd_one odd_one
  simpa using this

/-- **Lemma B for every rational pair**: if `mn(m² − n²) ≠ 0` then `m⁸ + 14m⁴n⁴ + n⁸` is not a square in `ℚ`
(from `LemmaB.lemmaB_rat`, for coprime integers, by scaling). -/
theorem lemmaB_rat_all {m n : ℚ} (hm : m ≠ 0) (hn : n ≠ 0) (hne : m ^ 2 ≠ n ^ 2) :
    ¬ IsSquare (m ^ 8 + 14 * m ^ 4 * n ^ 4 + n ^ 8) := by
  obtain ⟨p, q, c, hpq, hq, hc, rfl, rfl⟩ := exists_coprime_scale (m := m) hn
  have hp : p ≠ 0 := by rintro rfl; simp at hm
  have hne' : p ^ 2 ≠ q ^ 2 := by
    intro h; apply hne
    have h' : (p : ℚ) ^ 2 = (q : ℚ) ^ 2 := by exact_mod_cast h
    rw [mul_pow, mul_pow, h']
  intro hsq
  apply LemmaB.lemmaB_rat hpq hp hq hne'
  have e : (c * p) ^ 8 + 14 * (c * p) ^ 4 * (c * q) ^ 4 + (c * q) ^ 8 =
      (c ^ 4) ^ 2 * ((p : ℚ) ^ 8 + 14 * (p : ℚ) ^ 4 * (q : ℚ) ^ 4 + (q : ℚ) ^ 8) := by ring
  rw [e] at hsq
  exact isSquare_cancel_left ⟨c ^ 4, sq _⟩ (pow_ne_zero _ (pow_ne_zero _ hc)) hsq

/-- **The parity fact for every rational pair**: if `m² + n² = 2k²` with `n ≠ 0`, then `2mn` is not a square in
`ℚ`. (Scale to coprime `p, q`; `p² + q² = 2u²` with `u` an integer forces `p, q` odd.) -/
theorem not_isSquare_two_mul_all {m n k : ℚ} (hn : n ≠ 0) (hk : m ^ 2 + n ^ 2 = 2 * k ^ 2) :
    ¬ IsSquare (2 * m * n) := by
  obtain ⟨p, q, c, hpq, hq, hc, rfl, rfl⟩ := exists_coprime_scale (m := m) hn
  -- `2(p² + q²) = (2k/c)²` is an integer square, so `p² + q² = 2u²`
  have h1 : (((2 * (p ^ 2 + q ^ 2) : ℤ)) : ℚ) = (2 * k / c) * (2 * k / c) := by
    rw [show (2 * k / c) * (2 * k / c) = 4 * k ^ 2 / c ^ 2 by ring, eq_div_iff (pow_ne_zero 2 hc)]
    push_cast
    linear_combination 2 * hk
  obtain ⟨t, ht⟩ := Rat.isSquare_intCast_iff.mp ⟨_, h1⟩
  obtain ⟨u, rfl⟩ : Even t := by
    have : Even (t * t) := ⟨p ^ 2 + q ^ 2, by rw [← ht]; ring⟩
    exact (Int.even_mul.mp this).elim id id
  have hu : p ^ 2 + q ^ 2 = 2 * u ^ 2 :=
    mul_left_cancel₀ (two_ne_zero : (2 : ℤ) ≠ 0) (by linear_combination ht)
  obtain ⟨hpo, hqo⟩ := odd_of_twice_square hpq hu
  intro hsq
  apply not_isSquare_two_mul hpo hqo
  have e : 2 * (c * p) * (c * q) = c ^ 2 * (2 * (p : ℚ) * (q : ℚ)) := by ring
  rw [e] at hsq
  exact isSquare_cancel_left ⟨c, sq c⟩ (pow_ne_zero 2 hc) hsq

/-- At a twice-square pair with `m² ≠ n²`, neither `m` nor `n` is `0` (else `2` would be a square). -/
theorem ne_zero_of_twice_square {m n k : ℚ} (hk : m ^ 2 + n ^ 2 = 2 * k ^ 2) (hne : m ^ 2 ≠ n ^ 2) :
    m ≠ 0 ∧ n ≠ 0 := by
  have key : ∀ x y : ℚ, x ^ 2 + y ^ 2 = 2 * k ^ 2 → x ^ 2 ≠ y ^ 2 → y ≠ 0 := by
    intro x y hxy hne' hy
    subst hy
    have hx : x ≠ 0 := by rintro rfl; exact hne' rfl
    have hk0 : k ≠ 0 := by rintro rfl; apply hx; simpa using hxy
    apply not_isSquare_two
    refine ⟨x / k, ?_⟩
    field_simp
    linear_combination -hxy
  exact ⟨key n m (by linear_combination hk) (Ne.symm hne), key m n hk hne⟩

/-! ## Theorem 1 for every rational pair -/

/-- `T ∉ 2E⁺(ℚ)` for every rational pair with `mn(m² − n²) ≠ 0`. -/
theorem T_not_double_all {m n : ℚ} (hm : m ≠ 0) (hn : n ≠ 0) (hne : m ^ 2 ≠ n ^ 2)
    (hT : (Eplus m n).toAffine.Nonsingular (eP m n) 0)
    (W : (Eplus m n).toAffine.Point) : W + W ≠ Point.some _ _ hT := by
  intro h
  have hsq := double_T_isSquare (qe_ne_zero (n := n) hm) h
  rw [qe_P] at hsq
  have := isSquare_cancel_left ⟨2 * (m ^ 2 + n ^ 2), by ring⟩
    (mul_ne_zero (by norm_num) (pow_ne_zero 2 (S_ne_zero (n := n) hm))) hsq
  rw [gen_H_mul] at this
  exact lemmaB_rat_all hm hn hne this

theorem ne_of_sq_ne {m n : ℚ} (hne : m ^ 2 ≠ n ^ 2) : m + n ≠ 0 ∧ m - n ≠ 0 := by
  refine ⟨fun h => hne ?_, fun h => hne ?_⟩
  · have : m = -n := by linear_combination h
    rw [this]; ring
  · have : m = n := by linear_combination h
    rw [this]

/-- No point of order 4 on `E⁺(ℚ)`, for every rational pair with `mn(m² − n²) ≠ 0`. -/
theorem no_order_four_all {m n : ℚ} (hm : m ≠ 0) (hn : n ≠ 0) (hne : m ^ 2 ≠ n ^ 2)
    (W : (Eplus m n).toAffine.Point) (h4 : (W + W) + (W + W) = 0) : W + W = 0 := by
  obtain ⟨hp, hd⟩ := ne_of_sq_ne hne
  have hT := nonsingular_T (qe_ne_zero (n := n) hm)
  rcases two_torsion_point hm hn hp hd hT (W + W) h4 with h | h
  · exact h
  · exact absurd h (T_not_double_all hm hn hne hT W)

theorem two_pow_torsion_all {m n : ℚ} (hm : m ≠ 0) (hn : n ≠ 0) (hne : m ^ 2 ≠ n ^ 2)
    (a : ℕ) : ∀ W : (Eplus m n).toAffine.Point, (2 ^ a) • W = 0 → W + W = 0 := by
  induction a with
  | zero =>
    intro W h
    rw [pow_zero, one_smul] at h
    rw [h, add_zero]
  | succ a ih =>
    intro W h
    rw [pow_succ, mul_smul, two_nsmul] at h
    exact no_order_four_all hm hn hne W (ih (W + W) h)

/-- **Theorem 1, rank**: for all rational `m, n` with `mn(m² − n²) ≠ 0`, `rank E⁺(ℚ) ≥ 1` (`E⁺ = Eplus m n`). -/
theorem rank_one_all {m n : ℚ} (hm : m ≠ 0) (hn : n ≠ 0) (hne : m ^ 2 ≠ n ^ 2) :
    1 ≤ finrank ℤ (Eplus m n).toAffine.Point := by
  have hD := UniformRank.lemmaD
  have hmn1 : m ≠ n := by intro h; apply hne; rw [h]
  have hmn2 : m ≠ -n := by intro h; apply hne; rw [h]; ring
  obtain ⟨hp, hd⟩ := ne_of_sq_ne hne
  have hq := qe_ne_zero (n := n) hm
  have hS := S_ne_zero (n := n) hm
  have hH0 : genH m n ≠ 0 := (genH_pos hm).ne'
  have hyR : 2 * yR m n ≠ 0 := by
    unfold yR
    exact mul_ne_zero two_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero
      (mul_ne_zero (by norm_num) hm) hn) hp) hS) hH0)
  have hR : (Eplus m n).toAffine.Nonsingular (xR m n) (yR m n) := nonsingular_of_y (R_on _ _) hyR
  have hT : (Eplus m n).toAffine.Nonsingular (eP m n) 0 := nonsingular_T hq
  have hxR : xR m n ≠ eP m n := by
    intro h; apply mul_ne_zero (mul_ne_zero two_ne_zero hS) hH0; rw [← alphaR_val, h, sub_self]
  have aR : alpha (eP m n) (aP m n) (bP m n) (Point.some _ _ hR) = 2 * (m ^ 2 + n ^ 2) * genH m n := by
    rw [alpha_some, if_neg hxR, alphaR_val]
  have aT : alpha (eP m n) (aP m n) (bP m n) (Point.some _ _ hT) = qe (eP m n) (aP m n) (bP m n) := by
    rw [alpha_some, if_pos rfl]
  have nR : ¬ IsSquare (alpha (eP m n) (aP m n) (bP m n) (Point.some _ _ hR)) := by
    rw [aR]; intro h
    rcases lemmaD_bridge hD hn h with h1 | h1
    · exact hmn1 h1
    · exact hmn2 h1
  have nT : ¬ IsSquare (alpha (eP m n) (aP m n) (bP m n) (Point.some _ _ hT)) := by
    rw [aT, qe_P]; intro h
    have := isSquare_cancel_left ⟨2 * (m ^ 2 + n ^ 2), by ring⟩
      (mul_ne_zero (by norm_num) (pow_ne_zero 2 hS)) h
    rw [gen_H_mul] at this
    exact lemmaB_rat_all hm hn hne this
  have nRT : ¬ IsSquare (alpha (eP m n) (aP m n) (bP m n) (Point.some _ _ hR) *
      alpha (eP m n) (aP m n) (bP m n) (Point.some _ _ hT)) := by
    rw [aR, aT, qe_P]; intro h
    have h' : IsSquare (2 * (m ^ 2 + (-n) ^ 2) * genH m (-n)) := by
      refine isSquare_cancel_left (x := (2 * (m ^ 2 + n ^ 2) * genH m n) ^ 2)
        ⟨_, sq _⟩ (pow_ne_zero 2 (mul_ne_zero (mul_ne_zero two_ne_zero hS) hH0)) ?_
      convert h using 1; ring
    rcases lemmaD_bridge hD (neg_ne_zero.mpr hn) h' with h1 | h1
    · exact hmn2 h1
    · exact hmn1 (by rw [h1, neg_neg])
  have hind := indep2_of_alpha hq _ _ nR nT nRT
  have := fg_Eplus hm hn hp hd
  exact one_le_finrank_of_indep (two_torsion_point hm hn hp hd hT) hind

/-- `R` has infinite order on `E⁺(ℚ)`, for every rational pair with `mn(m² − n²) ≠ 0`. -/
theorem R_infinite_order_all {m n : ℚ} (hm : m ≠ 0) (hn : n ≠ 0) (hne : m ^ 2 ≠ n ^ 2)
    (hR : (Eplus m n).toAffine.Nonsingular (xR m n) (yR m n))
    (N : ℕ) (hN : 0 < N) : N • Point.some _ _ hR ≠ 0 := by
  intro h0
  have hD := UniformRank.lemmaD
  have hmn1 : m ≠ n := by intro h; apply hne; rw [h]
  have hmn2 : m ≠ -n := by intro h; apply hne; rw [h]; ring
  obtain ⟨hp, hd⟩ := ne_of_sq_ne hne
  have hq := qe_ne_zero (n := n) hm
  have hS := S_ne_zero (n := n) hm
  have hH0 : genH m n ≠ 0 := (genH_pos hm).ne'
  have hT : (Eplus m n).toAffine.Nonsingular (eP m n) 0 := nonsingular_T hq
  have hxR : xR m n ≠ eP m n := by
    intro h; apply mul_ne_zero (mul_ne_zero two_ne_zero hS) hH0; rw [← alphaR_val, h, sub_self]
  have aR : alpha (eP m n) (aP m n) (bP m n) (Point.some _ _ hR) = 2 * (m ^ 2 + n ^ 2) * genH m n := by
    rw [alpha_some, if_neg hxR, alphaR_val]
  have aT : alpha (eP m n) (aP m n) (bP m n) (Point.some _ _ hT) = qe (eP m n) (aP m n) (bP m n) := by
    rw [alpha_some, if_pos rfl]
  obtain ⟨a, d, hdo, rfl⟩ := Nat.exists_eq_two_pow_mul_odd hN.ne'
  set R := Point.some _ _ hR with hRdef
  have hV : (2 ^ a) • (d • R) = 0 := by rw [← mul_smul]; exact h0
  have hVV := two_pow_torsion_all hm hn hne a (d • R) hV
  obtain ⟨e, rfl⟩ := hdo
  have hdR : (2 * e + 1) • R = (e • R + e • R) + R := by
    rw [add_smul, one_smul, mul_comm, mul_smul, two_nsmul, nsmul_add]
  rcases two_torsion_point hm hn hp hd hT _ hVV with hV0 | hVT
  · have hR2 : R = (-(e • R)) + (-(e • R)) := by
      rw [hdR] at hV0
      rw [← sub_eq_zero, ← hV0]; abel
    have := alpha_double hq (-(e • R))
    rw [← hR2, aR] at this
    rcases lemmaD_bridge hD hn this with h1 | h1
    · exact hmn1 h1
    · exact hmn2 h1
  · have hTT : Point.some _ _ hT + Point.some _ _ hT = 0 := by rw [← hVT]; exact hVV
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
    have h'' : IsSquare (2 * (m ^ 2 + (-n) ^ 2) * genH m (-n)) := by
      refine isSquare_cancel_left (x := (2 * (m ^ 2 + n ^ 2) * genH m n) ^ 2)
        ⟨_, sq _⟩ (pow_ne_zero 2 (mul_ne_zero (mul_ne_zero two_ne_zero hS) hH0)) ?_
      convert h' using 1; ring
    rcases lemmaD_bridge hD (neg_ne_zero.mpr hn) h'' with h1 | h1
    · exact hmn2 h1
    · exact hmn1 (by rw [h1, neg_neg])

/-- **Theorem 1, the point**: for all rational `m, n` with `mn(m² − n²) ≠ 0`, the point `P = 2R` of `E⁺`
(`z = 2`) is nonsingular and of infinite order. On `E⁻` (the point of `z = −2`) apply it to `(m, −n)`. -/
theorem theorem1_P_all {m n : ℚ} (hm : m ≠ 0) (hn : n ≠ 0) (hne : m ^ 2 ≠ n ^ 2) :
    ∃ hP : (Eplus m n).toAffine.Nonsingular (xP m n) (yP m n),
      ∀ N : ℕ, 0 < N → N • Point.some _ _ hP ≠ 0 := by
  obtain ⟨hp, hd⟩ := ne_of_sq_ne hne
  obtain ⟨hR, hP, hRR⟩ := R_add_R hm hn hp hd
  refine ⟨hP, fun N hN h0 => R_infinite_order_all hm hn hne hR (2 * N) (by omega) ?_⟩
  rw [mul_comm, mul_smul, two_nsmul, hRR]; exact h0

/-- Theorem 1, the point, on `E⁻ = Eplus m (−n)`: its `P` (the image of `z = −2`, `EminusModel.Eminus_xP`) is
nonsingular and of infinite order. -/
theorem theorem1_P_all_minus {m n : ℚ} (hm : m ≠ 0) (hn : n ≠ 0) (hne : m ^ 2 ≠ n ^ 2) :
    ∃ hP : (Eplus m (-n)).toAffine.Nonsingular (xP m (-n)) (yP m (-n)),
      ∀ N : ℕ, 0 < N → N • Point.some _ _ hP ≠ 0 :=
  theorem1_P_all hm (neg_ne_zero.mpr hn) (by rwa [neg_sq])

/-- **Theorem 1 for every rational pair** (rank `≥ 1`): for all `m, n : ℚ` with `mn(m² − n²) ≠ 0`, both
`E⁺ = Eplus m n` and `E⁻ = Eplus m (−n)` have rank at least 1. -/
theorem theorem1_all {m n : ℚ} (hm : m ≠ 0) (hn : n ≠ 0) (hne : m ^ 2 ≠ n ^ 2) :
    1 ≤ finrank ℤ (Eplus m n).toAffine.Point ∧ 1 ≤ finrank ℤ (Eplus m (-n)).toAffine.Point :=
  ⟨rank_one_all hm hn hne, rank_one_all hm (neg_ne_zero.mpr hn) (by rwa [neg_sq])⟩

/-! ## Theorem 2 for every rational pair -/

/-- **Theorem 2 for every rational pair** (rank `≥ 2`): for all `m, n, k : ℚ` with `m² + n² = 2k²` and
`m² ≠ n²`, the points `R, Q, T` are independent in `E(ℚ)/2E(ℚ)` on `E⁺ = Eplus m n` and on
`E⁻ = Eplus m (−n)`, and both have rank at least 2. -/
theorem theorem2_all {m n k : ℚ} (hk : m ^ 2 + n ^ 2 = 2 * k ^ 2) (hne : m ^ 2 ≠ n ^ 2) :
    Theorem2At m n k ∧ Theorem2At m (-n) k ∧
      2 ≤ finrank ℤ (Eplus m n).toAffine.Point ∧ 2 ≤ finrank ℤ (Eplus m (-n)).toAffine.Point := by
  have hD := UniformRank.lemmaD
  obtain ⟨hm0, hn0⟩ := ne_zero_of_twice_square hk hne
  have h1 : m ≠ n := by intro h; apply hne; rw [h]
  have h2 : m ≠ -n := by intro h; apply hne; rw [h]; ring
  obtain ⟨hp, hd⟩ := ne_of_sq_ne hne
  have hp' : m + -n ≠ 0 := by intro h; apply h1; linear_combination h
  have hd' : m - -n ≠ 0 := by intro h; apply h2; linear_combination h
  have hk' : m ^ 2 + (-n) ^ 2 = 2 * k ^ 2 := by rw [neg_sq]; exact hk
  have hH : ¬ IsSquare (genH m n) := not_isSquare_genH hD hn0 hk h1 h2
  have hH' : ¬ IsSquare (genH m (-n)) :=
    not_isSquare_genH hD (neg_ne_zero.mpr hn0) hk' h2 (by rw [neg_neg]; exact h1)
  have hHH : ¬ IsSquare (genH m n * genH m (-n)) := by
    rw [gen_H_mul]; exact lemmaB_rat_all hm0 hn0 hne
  have s2 : ¬ IsSquare (2 * m * n) := not_isSquare_two_mul_all hn0 hk
  have s2' : ¬ IsSquare (2 * m * -n) := not_isSquare_two_mul_all (neg_ne_zero.mpr hn0) hk'
  have t1 : Theorem2At m n k :=
    theorem2_of _ _ _ hk hm0 hn0 hp hd hH hH' hHH s2
      (by rw [show (-2 : ℚ) * m * n = 2 * m * -n by ring]; exact s2')
  have t2 : Theorem2At m (-n) k :=
    theorem2_of _ _ _ hk' hm0 (neg_ne_zero.mpr hn0) hp' hd' hH' (by rwa [neg_neg])
      (by rwa [neg_neg, mul_comm]) s2' (by rw [show (-2 : ℚ) * m * -n = 2 * m * n by ring]; exact s2)
  exact ⟨t1, t2, two_le_rank_of_theorem2At hm0 hn0 hp hd t1,
    two_le_rank_of_theorem2At hm0 (neg_ne_zero.mpr hn0) hp' hd' t2⟩

/-! ## The same, with the curves of the recipe of §1 -/

/-- **Theorem 1** on `E⁺ = EplusModel m n`, `E⁻ = Eminus m n` (`EminusModel.lean`), all rational pairs with
`mn(m² − n²) ≠ 0`. -/
theorem theorem1_all_Eminus {m n : ℚ} (hm : m ≠ 0) (hn : n ≠ 0) (hne : m ^ 2 ≠ n ^ 2) :
    1 ≤ finrank ℤ (EplusModel m n).toAffine.Point ∧ 1 ≤ finrank ℤ (Eminus m n).toAffine.Point := by
  rw [EplusModel_eq, Eminus_eq]; exact theorem1_all hm hn hne

/-- **Theorem 2** on `E⁺ = EplusModel m n`, `E⁻ = Eminus m n`, all rational pairs with `m² + n² = 2k²`,
`m² ≠ n²`. -/
theorem theorem2_all_Eminus {m n k : ℚ} (hk : m ^ 2 + n ^ 2 = 2 * k ^ 2) (hne : m ^ 2 ≠ n ^ 2) :
    2 ≤ finrank ℤ (EplusModel m n).toAffine.Point ∧ 2 ≤ finrank ℤ (Eminus m n).toAffine.Point := by
  rw [EplusModel_eq, Eminus_eq]
  obtain ⟨-, -, h1, h2⟩ := theorem2_all hk hne
  exact ⟨h1, h2⟩

/-- Integer form of Theorem 1: every pair of integers with `mn(m² − n²) ≠ 0`, coprime or not. -/
theorem theorem1_int {m n : ℤ} (hm : m ≠ 0) (hn : n ≠ 0) (hne : m ^ 2 ≠ n ^ 2) :
    1 ≤ finrank ℤ (Eminus (m : ℚ) (n : ℚ)).toAffine.Point ∧
      1 ≤ finrank ℤ (EplusModel (m : ℚ) (n : ℚ)).toAffine.Point := by
  have h := theorem1_all_Eminus (m := (m : ℚ)) (n := (n : ℚ)) (by exact_mod_cast hm) (by exact_mod_cast hn)
    (by exact_mod_cast hne)
  exact ⟨h.2, h.1⟩

/-! ## Controls -/

/-- A non-coprime pair is covered: `(14, 2)`, twice `(7, 1)`, has rank `≥ 2` on both curves. -/
theorem control_14_2 : 2 ≤ finrank ℤ (Eplus 14 2).toAffine.Point ∧
    2 ≤ finrank ℤ (Eplus 14 (-2)).toAffine.Point := by
  obtain ⟨-, -, h1, h2⟩ := theorem2_all (m := 14) (n := 2) (k := 10) (by norm_num) (by norm_num)
  exact ⟨h1, h2⟩

/-- A non-integral rational pair is covered: `(7/3, 1/3)`. -/
theorem control_7_3_1_3 : 2 ≤ finrank ℤ (Eplus (7 / 3) (1 / 3)).toAffine.Point :=
  (theorem2_all (m := 7 / 3) (n := 1 / 3) (k := 5 / 3) (by norm_num) (by norm_num)).2.2.1

/-- The hypothesis `m² ≠ n²` of Lemma B is needed over `ℚ` too: at `(1/2, 1/2)` the form is `(1/4)²·4²`. -/
theorem control_lemmaB_all_sharp : IsSquare ((1 / 2 : ℚ) ^ 8 + 14 * (1 / 2) ^ 4 * (1 / 2) ^ 4 + (1 / 2) ^ 8) :=
  ⟨1 / 4, by norm_num⟩

/-- The parity fact needs `m² + n² = 2k²`: `2·1·2 = 2²` at the pair `(1, 2)` (`1 + 4 = 5` is not twice a
square). -/
theorem control_parity_sharp : IsSquare (2 * (1 : ℚ) * 2) := ⟨2, by norm_num⟩

end MagicSquareOfSquares.UniformRankAll

#print axioms MagicSquareOfSquares.UniformRankAll.exists_coprime_scale
#print axioms MagicSquareOfSquares.UniformRankAll.lemmaB_rat_all
#print axioms MagicSquareOfSquares.UniformRankAll.not_isSquare_two_mul_all
#print axioms MagicSquareOfSquares.UniformRankAll.theorem1_all
#print axioms MagicSquareOfSquares.UniformRankAll.theorem1_P_all
#print axioms MagicSquareOfSquares.UniformRankAll.theorem1_P_all_minus
#print axioms MagicSquareOfSquares.UniformRankAll.theorem2_all
#print axioms MagicSquareOfSquares.UniformRankAll.theorem1_all_Eminus
#print axioms MagicSquareOfSquares.UniformRankAll.theorem2_all_Eminus
#print axioms MagicSquareOfSquares.UniformRankAll.theorem1_int
#print axioms MagicSquareOfSquares.UniformRankAll.control_14_2
