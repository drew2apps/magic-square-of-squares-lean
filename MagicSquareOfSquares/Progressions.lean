import MagicSquareOfSquares.Basic

/-!
# Magic squares of squares and arithmetic progressions of squares

Robertson (1996) showed that a 3×3 magic square of distinct squares exists if and only if
there are three arithmetic progressions of three squares each, all with the same common
difference, whose middle terms themselves form an arithmetic progression.

## Main result

* `exists_isMSS_iff_three_progressions`: the equivalence. The progressions are
  `x - n, x, x + n`, `y - n, y, y + n`, `z - n, z, z + n` with `x, y, z` in progression.
-/

namespace MagicSquareOfSquares

/-- A 3×3 magic square of distinct squares exists if and only if there are three 3-term
arithmetic progressions of squares with the same common difference `n` whose middle terms
`x, y, z` form an arithmetic progression (and the nine squares are distinct). -/
theorem exists_isMSS_iff_three_progressions :
    (∃ a b c d e f g h i T : ℤ, IsMSS a b c d e f g h i T ∧
        [a ^ 2, b ^ 2, c ^ 2, d ^ 2, e ^ 2, f ^ 2, g ^ 2, h ^ 2, i ^ 2].Nodup) ↔
    ∃ x y z n : ℤ, x + z = 2 * y ∧
      IsSquare (x - n) ∧ IsSquare x ∧ IsSquare (x + n) ∧
      IsSquare (y - n) ∧ IsSquare y ∧ IsSquare (y + n) ∧
      IsSquare (z - n) ∧ IsSquare z ∧ IsSquare (z + n) ∧
      [y + n, x - n, z, z - n, y, x + n, x, z + n, y - n].Nodup := by
  constructor
  · rintro ⟨a, b, c, d, e, f, g, h, i, T, H, hnd⟩
    have hT := H.sum_eq
    have HM : IsMagic (a ^ 2) (b ^ 2) (c ^ 2) (d ^ 2) (e ^ 2) (f ^ 2) (g ^ 2) (h ^ 2) (i ^ 2) T :=
      H
    obtain ⟨r1, r2, r3, c1, c2, c3, d1, d2⟩ := HM
    have e1 : e ^ 2 + (a ^ 2 - e ^ 2) = a ^ 2 := by ring
    have e2 : g ^ 2 - (a ^ 2 - e ^ 2) = b ^ 2 := by linarith
    have e4 : c ^ 2 - (a ^ 2 - e ^ 2) = d ^ 2 := by linarith
    have e6 : g ^ 2 + (a ^ 2 - e ^ 2) = f ^ 2 := by linarith
    have e8 : c ^ 2 + (a ^ 2 - e ^ 2) = h ^ 2 := by linarith
    have e9 : e ^ 2 - (a ^ 2 - e ^ 2) = i ^ 2 := by linarith
    refine ⟨g ^ 2, e ^ 2, c ^ 2, a ^ 2 - e ^ 2, by linarith, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_,
      ?_⟩
    · rw [e2]; exact ⟨b, sq b⟩
    · exact ⟨g, sq g⟩
    · rw [e6]; exact ⟨f, sq f⟩
    · rw [e9]; exact ⟨i, sq i⟩
    · exact ⟨e, sq e⟩
    · rw [e1]; exact ⟨a, sq a⟩
    · rw [e4]; exact ⟨d, sq d⟩
    · exact ⟨c, sq c⟩
    · rw [e8]; exact ⟨h, sq h⟩
    · rw [e1, e2, e4, e6, e8, e9]; exact hnd
  · rintro ⟨x, y, z, n, hxz, ⟨b, hb⟩, ⟨g, hg⟩, ⟨f, hf⟩, ⟨i, hi⟩, ⟨e, he⟩, ⟨a, ha⟩, ⟨d, hd⟩,
      ⟨c, hc⟩, ⟨h, hh⟩, hnd⟩
    refine ⟨a, b, c, d, e, f, g, h, i, 3 * y, ?_, ?_⟩
    · unfold IsMSS
      constructor <;> (simp only [sq]; linarith)
    · simpa only [sq, ← ha, ← hb, ← hc, ← hd, ← he, ← hf, ← hg, ← hh, ← hi] using hnd

end MagicSquareOfSquares
