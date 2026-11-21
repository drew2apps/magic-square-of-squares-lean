import MagicSquareOfSquares.Statement
import MagicSquareOfSquares.Curve

/-!
# The open problem in terms of rational squares and of an elliptic curve

## Main results

* `openStatement_iff_rat_progressions`: the arithmetic-progression reformulation with
  rational squares instead of integer squares.
* `openStatement_iff_three_doubles`: Robertson's elliptic-curve reformulation. A 3×3 magic
  square of distinct squares exists if and only if, for some rational `n`, the curve
  `y² = x³ - n²x` has three rational points, each the double of a rational point, whose
  `x`-coordinates `x, y, z` are in arithmetic progression (and the nine numbers `x ± n, x`,
  `y ± n, y`, `z ± n, z` are distinct).
-/

namespace MagicSquareOfSquares

private theorem exists_common_den (l : List ℚ) :
    ∃ D : ℕ, D ≠ 0 ∧ ∀ r ∈ l, ∃ m : ℤ, (m : ℚ) = r * D := by
  induction l with
  | nil => exact ⟨1, one_ne_zero, by simp⟩
  | cons r l ih =>
    obtain ⟨D, hD, h⟩ := ih
    refine ⟨r.den * D, mul_ne_zero r.den_nz hD, ?_⟩
    intro s hs
    rcases List.mem_cons.mp hs with h1 | h1
    · refine ⟨r.num * D, ?_⟩
      rw [h1]
      push_cast
      rw [← mul_assoc, Rat.mul_den_eq_num]
    · obtain ⟨m, hm⟩ := h s h1
      refine ⟨m * r.den, ?_⟩
      push_cast
      rw [hm]
      ring

private theorem isSquare_cast_sub {a b : ℤ} (h : IsSquare (a - b)) : IsSquare ((a : ℚ) - b) := by
  obtain ⟨r, hr⟩ := h
  exact ⟨r, by rw [← Int.cast_sub, hr, Int.cast_mul]⟩

private theorem isSquare_cast_add {a b : ℤ} (h : IsSquare (a + b)) : IsSquare ((a : ℚ) + b) := by
  obtain ⟨r, hr⟩ := h
  exact ⟨r, by rw [← Int.cast_add, hr, Int.cast_mul]⟩

private theorem isSquare_cast {a : ℤ} (h : IsSquare a) : IsSquare (a : ℚ) := by
  obtain ⟨r, hr⟩ := h
  exact ⟨r, by rw [hr, Int.cast_mul]⟩

/-- The arithmetic-progression reformulation with rational squares. -/
theorem openStatement_iff_rat_progressions :
    OpenStatement ↔ ∃ x y z n : ℚ, x + z = 2 * y ∧
      IsSquare (x - n) ∧ IsSquare x ∧ IsSquare (x + n) ∧
      IsSquare (y - n) ∧ IsSquare y ∧ IsSquare (y + n) ∧
      IsSquare (z - n) ∧ IsSquare z ∧ IsSquare (z + n) ∧
      [x - n, x, x + n, y - n, y, y + n, z - n, z, z + n].Nodup := by
  rw [openStatement_iff_three_progressions]
  constructor
  · rintro ⟨x, y, z, n, hxz, s1, s2, s3, s4, s5, s6, s7, s8, s9, hnd⟩
    refine ⟨x, y, z, n, by exact_mod_cast hxz, isSquare_cast_sub s1, isSquare_cast s2,
      isSquare_cast_add s3, isSquare_cast_sub s4, isSquare_cast s5, isSquare_cast_add s6,
      isSquare_cast_sub s7, isSquare_cast s8, isSquare_cast_add s9, ?_⟩
    have hmap : [(x : ℚ) - n, x, x + n, y - n, y, y + n, z - n, z, z + n] =
        [x - n, x, x + n, y - n, y, y + n, z - n, z, z + n].map (Int.cast : ℤ → ℚ) := by
      simp only [List.map_cons, List.map_nil, Int.cast_sub, Int.cast_add]
    rw [hmap]
    exact hnd.map Int.cast_injective
  · rintro ⟨x, y, z, n, hxz, ⟨a1, h1⟩, ⟨a2, h2⟩, ⟨a3, h3⟩, ⟨a4, h4⟩, ⟨a5, h5⟩, ⟨a6, h6⟩,
      ⟨a7, h7⟩, ⟨a8, h8⟩, ⟨a9, h9⟩, hnd⟩
    obtain ⟨D, hD, hall⟩ := exists_common_den [a1, a2, a3, a4, a5, a6, a7, a8, a9]
    obtain ⟨m1, hm1⟩ := hall a1 (by simp)
    obtain ⟨m2, hm2⟩ := hall a2 (by simp)
    obtain ⟨m3, hm3⟩ := hall a3 (by simp)
    obtain ⟨m4, hm4⟩ := hall a4 (by simp)
    obtain ⟨m5, hm5⟩ := hall a5 (by simp)
    obtain ⟨m6, hm6⟩ := hall a6 (by simp)
    obtain ⟨m7, hm7⟩ := hall a7 (by simp)
    obtain ⟨m8, hm8⟩ := hall a8 (by simp)
    obtain ⟨m9, hm9⟩ := hall a9 (by simp)
    have hDq : ((D : ℚ)) ^ 2 ≠ 0 := pow_ne_zero 2 (by exact_mod_cast hD)
    obtain ⟨X, hXdef⟩ : ∃ X : ℤ, X = m2 * m2 := ⟨_, rfl⟩
    obtain ⟨Y, hYdef⟩ : ∃ Y : ℤ, Y = m5 * m5 := ⟨_, rfl⟩
    obtain ⟨Z, hZdef⟩ : ∃ Z : ℤ, Z = m8 * m8 := ⟨_, rfl⟩
    obtain ⟨N, hNdef⟩ : ∃ N : ℤ, N = m3 * m3 - m2 * m2 := ⟨_, rfl⟩
    have hX : (X : ℚ) = (D : ℚ) ^ 2 * x := by
      rw [hXdef]
      push_cast
      rw [hm2]
      linear_combination (-(D : ℚ) ^ 2) * h2
    have hY : (Y : ℚ) = (D : ℚ) ^ 2 * y := by
      rw [hYdef]
      push_cast
      rw [hm5]
      linear_combination (-(D : ℚ) ^ 2) * h5
    have hZ : (Z : ℚ) = (D : ℚ) ^ 2 * z := by
      rw [hZdef]
      push_cast
      rw [hm8]
      linear_combination (-(D : ℚ) ^ 2) * h8
    have hN : (N : ℚ) = (D : ℚ) ^ 2 * n := by
      rw [hNdef]
      push_cast
      rw [hm3, hm2]
      linear_combination (-(D : ℚ) ^ 2) * (h3 - h2)
    have e1 : X - N = m1 * m1 := by
      apply Int.cast_injective (α := ℚ)
      push_cast
      rw [hX, hN, hm1]
      linear_combination (D : ℚ) ^ 2 * h1
    have e3 : X + N = m3 * m3 := by
      apply Int.cast_injective (α := ℚ)
      push_cast
      rw [hX, hN, hm3]
      linear_combination (D : ℚ) ^ 2 * h3
    have e4 : Y - N = m4 * m4 := by
      apply Int.cast_injective (α := ℚ)
      push_cast
      rw [hY, hN, hm4]
      linear_combination (D : ℚ) ^ 2 * h4
    have e6 : Y + N = m6 * m6 := by
      apply Int.cast_injective (α := ℚ)
      push_cast
      rw [hY, hN, hm6]
      linear_combination (D : ℚ) ^ 2 * h6
    have e7 : Z - N = m7 * m7 := by
      apply Int.cast_injective (α := ℚ)
      push_cast
      rw [hZ, hN, hm7]
      linear_combination (D : ℚ) ^ 2 * h7
    have e9 : Z + N = m9 * m9 := by
      apply Int.cast_injective (α := ℚ)
      push_cast
      rw [hZ, hN, hm9]
      linear_combination (D : ℚ) ^ 2 * h9
    have hXZ : X + Z = 2 * Y := by
      apply Int.cast_injective (α := ℚ)
      push_cast
      rw [hX, hY, hZ]
      linear_combination (D : ℚ) ^ 2 * hxz
    refine ⟨X, Y, Z, N, hXZ, ⟨m1, e1⟩, ⟨m2, hXdef⟩, ⟨m3, e3⟩, ⟨m4, e4⟩, ⟨m5, hYdef⟩, ⟨m6, e6⟩,
      ⟨m7, e7⟩, ⟨m8, hZdef⟩, ⟨m9, e9⟩, ?_⟩
    have hmap : [X - N, X, X + N, Y - N, Y, Y + N, Z - N, Z, Z + N].map (Int.cast : ℤ → ℚ) =
        [x - n, x, x + n, y - n, y, y + n, z - n, z, z + n].map (fun t => (D : ℚ) ^ 2 * t) := by
      simp only [List.map_cons, List.map_nil, Int.cast_sub, Int.cast_add, hX, hY, hZ, hN,
        mul_sub, mul_add]
    exact List.Nodup.of_map _ (hmap ▸ hnd.map (mul_right_injective₀ hDq))

private theorem ne_zero_of_nodup {x y z n : ℚ}
    (hnd : [x - n, x, x + n, y - n, y, y + n, z - n, z, z + n].Nodup) : n ≠ 0 := by
  intro h0
  have hnot := (List.nodup_cons.mp hnd).1
  apply hnot
  rw [h0, sub_zero]
  simp

/-- Robertson's elliptic-curve reformulation: a 3×3 magic square of distinct squares exists
if and only if, for some rational `n`, the curve `y² = x³ - n²x` has three rational points,
each the double of a rational point, whose `x`-coordinates `x, y, z` form an arithmetic
progression, with the nine numbers `x - n, x, x + n, y - n, y, y + n, z - n, z, z + n`
distinct. -/
theorem openStatement_iff_three_doubles :
    OpenStatement ↔ ∃ n x y z : ℚ, x + z = 2 * y ∧
      IsDoubleX n x ∧ IsDoubleX n y ∧ IsDoubleX n z ∧
      [x - n, x, x + n, y - n, y, y + n, z - n, z, z + n].Nodup := by
  rw [openStatement_iff_rat_progressions]
  constructor
  · rintro ⟨x, y, z, n, hxz, s1, s2, s3, s4, s5, s6, s7, s8, s9, hnd⟩
    have hn := ne_zero_of_nodup hnd
    exact ⟨n, x, y, z, hxz, (isSquare_iff_isDoubleX hn).mp ⟨s1, s2, s3⟩,
      (isSquare_iff_isDoubleX hn).mp ⟨s4, s5, s6⟩, (isSquare_iff_isDoubleX hn).mp ⟨s7, s8, s9⟩,
      hnd⟩
  · rintro ⟨n, x, y, z, hxz, dx, dy, dz, hnd⟩
    have hn := ne_zero_of_nodup hnd
    obtain ⟨s1, s2, s3⟩ := (isSquare_iff_isDoubleX hn).mpr dx
    obtain ⟨s4, s5, s6⟩ := (isSquare_iff_isDoubleX hn).mpr dy
    obtain ⟨s7, s8, s9⟩ := (isSquare_iff_isDoubleX hn).mpr dz
    exact ⟨x, y, z, n, hxz, s1, s2, s3, s4, s5, s6, s7, s8, s9, hnd⟩

end MagicSquareOfSquares
