import MagicSquareOfSquares.Elliptic

/-!
# The open problem in terms of right triangles

Robertson's second reformulation: a 3×3 magic square of distinct squares exists if and only
if there are three rational right triangles with the same area whose squared hypotenuses
form an arithmetic progression.

A right triangle with legs `p, q` and hypotenuse `r` gives the three squares
`(q - p)², r², (q + p)²` in arithmetic progression with common difference `2pq`, which is
four times the area.

## Main result

* `openStatement_iff_three_triangles`
-/

namespace MagicSquareOfSquares

private theorem triangle_of_progression {x n a b c : ℚ} (hn : 0 < n) (ha : x - n = a * a)
    (hb : x = b * b) (hc : x + n = c * c) :
    ∃ p q r : ℚ, 0 < p ∧ 0 < q ∧ 0 < r ∧ p ^ 2 + q ^ 2 = r ^ 2 ∧ 2 * (p * q) = n ∧
      (q - p) ^ 2 = x - n ∧ r ^ 2 = x ∧ (q + p) ^ 2 = x + n := by
  have hA : |a| ^ 2 = x - n := by rw [sq_abs, ha]; ring
  have hB : |b| ^ 2 = x := by rw [sq_abs, hb]; ring
  have hC : |c| ^ 2 = x + n := by rw [sq_abs, hc]; ring
  have hac : |a| < |c| := by
    apply sq_lt_sq.mp
    nlinarith
  have ha0 : 0 ≤ |a| := abs_nonneg a
  have hb0 : 0 < |b| := by
    apply abs_pos.mpr
    rintro rfl
    nlinarith [mul_self_nonneg a]
  refine ⟨(|c| - |a|) / 2, (|c| + |a|) / 2, |b|, by linarith, by linarith, hb0, ?_, ?_, ?_, hB,
    ?_⟩
  · linear_combination (1 / 2) * hC + (1 / 2) * hA - hB
  · linear_combination (1 / 2) * hC - (1 / 2) * hA
  · linear_combination hA
  · linear_combination hC

private theorem nodup_neg {x y z n : ℚ}
    (hnd : [x - n, x, x + n, y - n, y, y + n, z - n, z, z + n].Nodup) :
    [x - -n, x, x + -n, y - -n, y, y + -n, z - -n, z, z + -n].Nodup := by
  simp only [sub_neg_eq_add, ← sub_eq_add_neg]
  let f : Fin 9 → ℚ := ![x - n, x, x + n, y - n, y, y + n, z - n, z, z + n]
  have h1 : [x + n, x, x - n, y + n, y, y - n, z + n, z, z - n] =
      ([2, 1, 0, 5, 4, 3, 8, 7, 6] : List (Fin 9)).map f := rfl
  have h2 : [x - n, x, x + n, y - n, y, y + n, z - n, z, z + n] =
      ([0, 1, 2, 3, 4, 5, 6, 7, 8] : List (Fin 9)).map f := rfl
  rw [h1]
  rw [h2] at hnd
  exact (List.Perm.nodup_iff (List.Perm.map f (by decide))).mpr hnd

private theorem triangles_of_pos {x y z n a1 a2 a3 a4 a5 a6 a7 a8 a9 : ℚ} (hn : 0 < n)
    (hxz : x + z = 2 * y) (h1 : x - n = a1 * a1) (h2 : x = a2 * a2) (h3 : x + n = a3 * a3)
    (h4 : y - n = a4 * a4) (h5 : y = a5 * a5) (h6 : y + n = a6 * a6)
    (h7 : z - n = a7 * a7) (h8 : z = a8 * a8) (h9 : z + n = a9 * a9)
    (hnd : [x - n, x, x + n, y - n, y, y + n, z - n, z, z + n].Nodup) :
    ∃ p₁ q₁ r₁ p₂ q₂ r₂ p₃ q₃ r₃ : ℚ,
      (0 < p₁ ∧ 0 < q₁ ∧ 0 < r₁) ∧ (0 < p₂ ∧ 0 < q₂ ∧ 0 < r₂) ∧ (0 < p₃ ∧ 0 < q₃ ∧ 0 < r₃) ∧
      p₁ ^ 2 + q₁ ^ 2 = r₁ ^ 2 ∧ p₂ ^ 2 + q₂ ^ 2 = r₂ ^ 2 ∧ p₃ ^ 2 + q₃ ^ 2 = r₃ ^ 2 ∧
      p₁ * q₁ = p₂ * q₂ ∧ p₂ * q₂ = p₃ * q₃ ∧
      r₁ ^ 2 + r₃ ^ 2 = 2 * r₂ ^ 2 ∧
      [(q₁ - p₁) ^ 2, r₁ ^ 2, (q₁ + p₁) ^ 2, (q₂ - p₂) ^ 2, r₂ ^ 2, (q₂ + p₂) ^ 2,
        (q₃ - p₃) ^ 2, r₃ ^ 2, (q₃ + p₃) ^ 2].Nodup := by
  obtain ⟨p₁, q₁, r₁, hp1, hq1, hr1, t1, ar1, lo1, mid1, hi1⟩ :=
    triangle_of_progression hn h1 h2 h3
  obtain ⟨p₂, q₂, r₂, hp2, hq2, hr2, t2, ar2, lo2, mid2, hi2⟩ :=
    triangle_of_progression hn h4 h5 h6
  obtain ⟨p₃, q₃, r₃, hp3, hq3, hr3, t3, ar3, lo3, mid3, hi3⟩ :=
    triangle_of_progression hn h7 h8 h9
  refine ⟨p₁, q₁, r₁, p₂, q₂, r₂, p₃, q₃, r₃, ⟨hp1, hq1, hr1⟩, ⟨hp2, hq2, hr2⟩,
    ⟨hp3, hq3, hr3⟩, t1, t2, t3, by linarith, by linarith, ?_, ?_⟩
  · rw [mid1, mid2, mid3]
    exact hxz
  · rw [lo1, mid1, hi1, lo2, mid2, hi2, lo3, mid3, hi3]
    exact hnd

/-- Robertson's right-triangle reformulation: a 3×3 magic square of distinct squares exists
if and only if there are three rational right triangles, with legs `pₖ, qₖ` and hypotenuse
`rₖ`, that have the same area and whose squared hypotenuses form an arithmetic progression,
with the nine squares `(qₖ - pₖ)², rₖ², (qₖ + pₖ)²` distinct. -/
theorem openStatement_iff_three_triangles :
    OpenStatement ↔ ∃ p₁ q₁ r₁ p₂ q₂ r₂ p₃ q₃ r₃ : ℚ,
      (0 < p₁ ∧ 0 < q₁ ∧ 0 < r₁) ∧ (0 < p₂ ∧ 0 < q₂ ∧ 0 < r₂) ∧ (0 < p₃ ∧ 0 < q₃ ∧ 0 < r₃) ∧
      p₁ ^ 2 + q₁ ^ 2 = r₁ ^ 2 ∧ p₂ ^ 2 + q₂ ^ 2 = r₂ ^ 2 ∧ p₃ ^ 2 + q₃ ^ 2 = r₃ ^ 2 ∧
      p₁ * q₁ = p₂ * q₂ ∧ p₂ * q₂ = p₃ * q₃ ∧
      r₁ ^ 2 + r₃ ^ 2 = 2 * r₂ ^ 2 ∧
      [(q₁ - p₁) ^ 2, r₁ ^ 2, (q₁ + p₁) ^ 2, (q₂ - p₂) ^ 2, r₂ ^ 2, (q₂ + p₂) ^ 2,
        (q₃ - p₃) ^ 2, r₃ ^ 2, (q₃ + p₃) ^ 2].Nodup := by
  rw [openStatement_iff_rat_progressions]
  constructor
  · rintro ⟨x, y, z, n, hxz, ⟨a1, h1⟩, ⟨a2, h2⟩, ⟨a3, h3⟩, ⟨a4, h4⟩, ⟨a5, h5⟩, ⟨a6, h6⟩,
      ⟨a7, h7⟩, ⟨a8, h8⟩, ⟨a9, h9⟩, hnd⟩
    rcases lt_trichotomy n 0 with hneg | h0 | hpos
    · exact triangles_of_pos (n := -n) (by linarith) hxz
        (by rw [sub_neg_eq_add]; exact h3) h2 (by rw [← sub_eq_add_neg]; exact h1)
        (by rw [sub_neg_eq_add]; exact h6) h5 (by rw [← sub_eq_add_neg]; exact h4)
        (by rw [sub_neg_eq_add]; exact h9) h8 (by rw [← sub_eq_add_neg]; exact h7)
        (nodup_neg hnd)
    · exfalso
      have hnot := (List.nodup_cons.mp hnd).1
      apply hnot
      rw [h0, sub_zero]
      simp
    · exact triangles_of_pos hpos hxz h1 h2 h3 h4 h5 h6 h7 h8 h9 hnd
  · rintro ⟨p₁, q₁, r₁, p₂, q₂, r₂, p₃, q₃, r₃, -, -, -, t1, t2, t3, ar12, ar23, hap, hnd⟩
    have e1 : r₁ ^ 2 - 2 * (p₁ * q₁) = (q₁ - p₁) ^ 2 := by linear_combination -t1
    have e3 : r₁ ^ 2 + 2 * (p₁ * q₁) = (q₁ + p₁) ^ 2 := by linear_combination -t1
    have e4 : r₂ ^ 2 - 2 * (p₁ * q₁) = (q₂ - p₂) ^ 2 := by linear_combination -t2 - 2 * ar12
    have e6 : r₂ ^ 2 + 2 * (p₁ * q₁) = (q₂ + p₂) ^ 2 := by linear_combination -t2 + 2 * ar12
    have e7 : r₃ ^ 2 - 2 * (p₁ * q₁) = (q₃ - p₃) ^ 2 := by
      linear_combination -t3 - 2 * ar12 - 2 * ar23
    have e9 : r₃ ^ 2 + 2 * (p₁ * q₁) = (q₃ + p₃) ^ 2 := by
      linear_combination -t3 + 2 * ar12 + 2 * ar23
    refine ⟨r₁ ^ 2, r₂ ^ 2, r₃ ^ 2, 2 * (p₁ * q₁), hap, ⟨q₁ - p₁, by rw [e1, sq]⟩, ⟨r₁, sq r₁⟩,
      ⟨q₁ + p₁, by rw [e3, sq]⟩, ⟨q₂ - p₂, by rw [e4, sq]⟩, ⟨r₂, sq r₂⟩,
      ⟨q₂ + p₂, by rw [e6, sq]⟩, ⟨q₃ - p₃, by rw [e7, sq]⟩, ⟨r₃, sq r₃⟩,
      ⟨q₃ + p₃, by rw [e9, sq]⟩, ?_⟩
    rw [e1, e3, e4, e6, e7, e9]
    exact hnd

end MagicSquareOfSquares
