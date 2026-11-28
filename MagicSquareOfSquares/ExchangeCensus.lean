import MagicSquareOfSquares.Exchange
import MagicSquareOfSquares.CoverDescentTable

/-!
# The second ratio of a solution is in the census list too

`G1.fullSquares_exchange` (`Exchange.lean`) with the census of this library: if the five cells of the G1 family are
squares at the ratio `m : n` and the parameter `t = u/v`, where `u : v` is a coprime pair with `v < u ≤ 1000`, then
`(u, v)` is one of the 425 pairs of `G1.exceptions1000` (`no_solution_le_1000`) and, more precisely, one of the 306
pairs of `CoverDescent.exceptions` (`CoverDescent.no_fullSquares_of_mem_ratios`). So the parameter of a solution
obeys the same restrictions as the ratio. Nothing is assumed about `m : n`.

The five forms depend on `t` only up to sign and inverse: `G1.fullSquares_neg` (`F0, F4, F8` are even in `t`, and
`L, U` are exchanged) and `G1.fullSquares_inv` (for `t ≠ 0`, each form at `1/t` is a form at `t` divided by `t²`,
with `F0` and `F8` exchanged). So the corollary covers every parameter `t ∉ {0, 1, −1}`, written as
`t = ±(u/v)^{±1}` with `0 < v < u` coprime (`exists_ratio_of_ne`): `fullSquares_param_mem_exceptions`.

Written by an AI model (Claude); not read by a person.
-/

namespace MagicSquareOfSquares.CoverDescent

/-- **The parameter `u/v` of a solution, as a ratio, is one of the 306 pairs left by the census in Lean**
(for `v < u ≤ 1000` coprime; any `m`, any `n ≠ 0`). -/
theorem mem_exceptions_of_fullSquares {m n u v : ℕ} (hn : n ≠ 0) (hv : 0 < v) (hvu : v < u)
    (hu : u ≤ 1000) (hcop : Nat.Coprime u v) (h : G1.FullSquares m n ((u : ℚ) / v)) :
    (u, v) ∈ exceptions := by
  have h' := (G1.fullSquares_exchange hn hv.ne').mp h
  have hex : (u, v) ∈ G1.exceptions1000 := by
    by_contra hex
    exact G1.no_solution_le_1000 hv hvu hu hcop hex _ h'
  unfold exceptions
  rw [List.mem_filter]
  refine ⟨hex, decide_eq_true ?_⟩
  intro hr
  exact no_fullSquares_of_mem_ratios hr _ h'

end MagicSquareOfSquares.CoverDescent

/-! ## Every parameter: the sign and the inverse of `t` -/

namespace MagicSquareOfSquares.G1

theorem isSquare_transport {a b : ℚ} (h : IsSquare a) (e : a = b) : IsSquare b := e ▸ h

/-- The five forms at `−t` are those at `t`, with `L` and `U` exchanged. -/
theorem fullSquares_neg {m n : ℕ} {t : ℚ} : FullSquares m n (-t) ↔ FullSquares m n t := by
  unfold FullSquares DiagonalSquares
  constructor
  · rintro ⟨⟨h0, h4, h8⟩, hL, hU⟩
    exact ⟨⟨isSquare_transport h0 (by ring), isSquare_transport h4 (by ring), isSquare_transport h8 (by ring)⟩,
      isSquare_transport hU (by ring), isSquare_transport hL (by ring)⟩
  · rintro ⟨⟨h0, h4, h8⟩, hL, hU⟩
    exact ⟨⟨isSquare_transport h0 (by ring), isSquare_transport h4 (by ring), isSquare_transport h8 (by ring)⟩,
      isSquare_transport hU (by ring), isSquare_transport hL (by ring)⟩

/-- The five forms at `1/t` are those at `t` divided by `t²`, with `F0` and `F8` exchanged. -/
theorem fullSquares_inv {m n : ℕ} {t : ℚ} (ht : t ≠ 0) : FullSquares m n t⁻¹ ↔ FullSquares m n t := by
  have h1 : (1 : ℚ) ≠ 0 := one_ne_zero
  have i0 : IsSquare ((m : ℚ) ^ 2 + (n : ℚ) ^ 2 * t⁻¹ ^ 2) ↔ IsSquare ((n : ℚ) ^ 2 + (m : ℚ) ^ 2 * t ^ 2) :=
    isSquare_congr_of_mul_sq ht h1 (by field_simp; ring)
  have i8 : IsSquare ((n : ℚ) ^ 2 + (m : ℚ) ^ 2 * t⁻¹ ^ 2) ↔ IsSquare ((m : ℚ) ^ 2 + (n : ℚ) ^ 2 * t ^ 2) :=
    isSquare_congr_of_mul_sq ht h1 (by field_simp; ring)
  have i4 : IsSquare (((m : ℚ) ^ 2 + (n : ℚ) ^ 2) * (t⁻¹ ^ 2 + 1) / 2) ↔
      IsSquare (((m : ℚ) ^ 2 + (n : ℚ) ^ 2) * (t ^ 2 + 1) / 2) :=
    isSquare_congr_of_mul_sq ht h1 (by field_simp; ring)
  have iL : IsSquare (((m : ℚ) ^ 2 + (n : ℚ) ^ 2) * (t⁻¹ ^ 2 + 1) / 2 - 2 * m * n * t⁻¹) ↔
      IsSquare (((m : ℚ) ^ 2 + (n : ℚ) ^ 2) * (t ^ 2 + 1) / 2 - 2 * m * n * t) :=
    isSquare_congr_of_mul_sq ht h1 (by field_simp; ring)
  have iU : IsSquare (((m : ℚ) ^ 2 + (n : ℚ) ^ 2) * (t⁻¹ ^ 2 + 1) / 2 + 2 * m * n * t⁻¹) ↔
      IsSquare (((m : ℚ) ^ 2 + (n : ℚ) ^ 2) * (t ^ 2 + 1) / 2 + 2 * m * n * t) :=
    isSquare_congr_of_mul_sq ht h1 (by field_simp; ring)
  unfold FullSquares DiagonalSquares
  constructor
  · rintro ⟨⟨h0, h4, h8⟩, hL, hU⟩
    exact ⟨⟨i8.mp h8, i4.mp h4, i0.mp h0⟩, iL.mp hL, iU.mp hU⟩
  · rintro ⟨⟨h0, h4, h8⟩, hL, hU⟩
    exact ⟨⟨i0.mpr h8, i4.mpr h4, i8.mpr h0⟩, iL.mpr hL, iU.mpr hU⟩

end MagicSquareOfSquares.G1

namespace MagicSquareOfSquares.CoverDescent

/-- Every rational `t ∉ {0, 1, −1}` is `±u/v` or `±v/u` for a coprime pair `0 < v < u`. -/
theorem exists_ratio_of_ne {t : ℚ} (h0 : t ≠ 0) (h1 : t ≠ 1) (h2 : t ≠ -1) :
    ∃ u v : ℕ, 0 < v ∧ v < u ∧ Nat.Coprime u v ∧
      (t = (u : ℚ) / v ∨ t = -((u : ℚ) / v) ∨ t = (v : ℚ) / u ∨ t = -((v : ℚ) / u)) := by
  have hcop : Nat.Coprime t.num.natAbs t.den := t.reduced
  have hden : 0 < t.den := t.den_pos
  have hnum : t.num ≠ 0 := Rat.num_ne_zero.mpr h0
  have hapos : 0 < t.num.natAbs := Int.natAbs_pos.mpr hnum
  have e : (t.num : ℚ) / t.den = t := Rat.num_div_den t
  have hsign : t = (t.num.natAbs : ℚ) / t.den ∨ t = -((t.num.natAbs : ℚ) / t.den) := by
    rcases le_or_gt 0 t.num with hle | hlt
    · left
      have h3 : ((t.num.natAbs : ℕ) : ℤ) = t.num := Int.natAbs_of_nonneg hle
      have hc : ((t.num.natAbs : ℕ) : ℚ) = (t.num : ℚ) := by
        rw [← Int.cast_natCast, h3]
      calc t = (t.num : ℚ) / t.den := e.symm
        _ = (t.num.natAbs : ℚ) / t.den := by rw [hc]
    · right
      have h3 : ((t.num.natAbs : ℕ) : ℤ) = -t.num := by
        rw [← Int.natAbs_neg]
        exact Int.natAbs_of_nonneg (by omega)
      have hc : ((t.num.natAbs : ℕ) : ℚ) = -(t.num : ℚ) := by
        rw [← Int.cast_natCast, h3, Int.cast_neg]
      calc t = (t.num : ℚ) / t.den := e.symm
        _ = -((t.num.natAbs : ℚ) / t.den) := by rw [hc]; ring
  have hne : t.num.natAbs ≠ t.den := by
    intro heq
    have hc1 : Nat.Coprime t.den t.den := by rw [heq] at hcop; exact hcop
    have hd1 : t.den = 1 := (Nat.coprime_self _).mp hc1
    rcases hsign with hs | hs
    · exact h1 (by rw [hs, heq, hd1]; simp)
    · exact h2 (by rw [hs, heq, hd1]; simp)
  rcases lt_or_gt_of_ne hne with hlt | hgt
  · refine ⟨t.den, t.num.natAbs, hapos, hlt, hcop.symm, ?_⟩
    rcases hsign with hs | hs
    · exact Or.inr (Or.inr (Or.inl hs))
    · exact Or.inr (Or.inr (Or.inr hs))
  · refine ⟨t.num.natAbs, t.den, hden, hgt, hcop, ?_⟩
    rcases hsign with hs | hs
    · exact Or.inl hs
    · exact Or.inr (Or.inl hs)

/-- The corollary for `t = ±(u/v)^{±1}`: the coprime pair `0 < v < u ≤ 1000` is one of the 306 pairs. -/
theorem mem_exceptions_of_fullSquares_any {m n : ℕ} (hn : n ≠ 0) {t : ℚ} (h : G1.FullSquares m n t)
    {u v : ℕ} (hv : 0 < v) (hvu : v < u) (hu : u ≤ 1000) (hcop : Nat.Coprime u v)
    (ht : t = (u : ℚ) / v ∨ t = -((u : ℚ) / v) ∨ t = (v : ℚ) / u ∨ t = -((v : ℚ) / u)) :
    (u, v) ∈ exceptions := by
  have hu0 : (u : ℚ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  have hv0 : (v : ℚ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  have hvu0 : (v : ℚ) / u ≠ 0 := div_ne_zero hv0 hu0
  have huv : G1.FullSquares m n ((u : ℚ) / v) := by
    rcases ht with rfl | rfl | rfl | rfl
    · exact h
    · exact G1.fullSquares_neg.mp h
    · have h' := (G1.fullSquares_inv hvu0).mpr h
      rwa [inv_div] at h'
    · have h' := (G1.fullSquares_inv hvu0).mpr (G1.fullSquares_neg.mp h)
      rwa [inv_div] at h'
  exact mem_exceptions_of_fullSquares hn hv hvu hu hcop huv

/-- **Every parameter.** If the five cells of the G1 family are squares at the ratio `m : n` (any `m`, any `n ≠ 0`)
and a parameter `t ∉ {0, 1, −1}`, then `t = ±(u/v)^{±1}` for a coprime pair `0 < v < u`, and if `u ≤ 1000` the pair
`(u, v)` is one of the 306 pairs left by the census in Lean. -/
theorem fullSquares_param_mem_exceptions {m n : ℕ} (hn : n ≠ 0) {t : ℚ} (h : G1.FullSquares m n t)
    (h0 : t ≠ 0) (h1 : t ≠ 1) (h2 : t ≠ -1) :
    ∃ u v : ℕ, 0 < v ∧ v < u ∧ Nat.Coprime u v ∧
      (t = (u : ℚ) / v ∨ t = -((u : ℚ) / v) ∨ t = (v : ℚ) / u ∨ t = -((v : ℚ) / u)) ∧
      (u ≤ 1000 → (u, v) ∈ exceptions) := by
  obtain ⟨u, v, hv, hvu, hcop, ht⟩ := exists_ratio_of_ne h0 h1 h2
  exact ⟨u, v, hv, hvu, hcop, ht, fun hu => mem_exceptions_of_fullSquares_any hn h hv hvu hu hcop ht⟩

end MagicSquareOfSquares.CoverDescent
