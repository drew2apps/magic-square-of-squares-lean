import MagicSquareOfSquares.G1Size
import MagicSquareOfSquares.Statement

/-!
# The size bound for a square with distinct entries

`IsMSS.edge_root_gt_500` assumes three inequalities between squares of edge roots. For a magic
square of squares with nine distinct entries, as in the open problem, they hold automatically.
This file states the bound in that form, for the definition used in formal-conjectures.
-/

namespace MagicSquareOfSquares

/-- **In a 3×3 magic square of nine distinct squares, one of the four edge entries is larger
than `250000 = 500²`.** Published computer searches exclude far larger squares; this is only
the bound that follows from the theorems of this repository. -/
theorem IsMagicSquareOfSquares.edge_gt {m : Fin 3 → Fin 3 → ℕ} {t : ℕ}
    (H : IsMagicSquareOfSquares m t) (hinj : m.Injective2) :
    250000 < m 0 1 ∨ 250000 < m 1 0 ∨ 250000 < m 1 2 ∨ 250000 < m 2 1 := by
  obtain ⟨r, hr, hs⟩ := H.exists_roots
  have ne : ∀ i j i' j' : Fin 3, (i, j) ≠ (i', j') → (r i j : ℤ) ^ 2 ≠ (r i' j' : ℤ) ^ 2 := by
    intro i j i' j' hij h
    have h' : m i j = m i' j' := by
      rw [hr i j, hr i' j']
      exact_mod_cast h
    obtain ⟨h1, h2⟩ := hinj h'
    exact hij (by rw [h1, h2])
  have key : ∀ i j, 500 < |(r i j : ℤ)| → 250000 < m i j := by
    intro i j h
    rw [abs_of_nonneg (Int.natCast_nonneg _)] at h
    have h' : 500 < r i j := by exact_mod_cast h
    rw [hr i j]
    nlinarith
  rcases IsMSS.edge_root_gt_500 hs (ne 0 1 1 0 (by decide)) (ne 2 1 1 2 (by decide))
    (ne 0 1 1 2 (by decide)) with h | h | h | h
  · exact Or.inl (key _ _ h)
  · exact Or.inr (Or.inl (key _ _ h))
  · exact Or.inr (Or.inr (Or.inl (key _ _ h)))
  · exact Or.inr (Or.inr (Or.inr (key _ _ h)))

end MagicSquareOfSquares
