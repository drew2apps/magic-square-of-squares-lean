import MagicSquareOfSquares.G1Size
import MagicSquareOfSquares.G1Norm
import MagicSquareOfSquares.Statement

/-!
# The G1 results for the matrix form of the problem

The theorems on the edge roots are stated in `G1Size.lean` and `G1Norm.lean` for nine integer
roots. This file restates them for a matrix `m : Fin 3 → Fin 3 → ℕ` of squares, the definition
used in formal-conjectures. The edge entries are `m 0 1`, `m 1 0`, `m 1 2`, `m 2 1`.

* `IsMagicSquareOfSquares.edge_gt_of_ne`, `edge_gt`: one edge entry exceeds `250000 = 500²`.
  Only three inequalities between edge entries are used; they hold when the nine entries are
  distinct.
* `IsMagicSquareOfSquares.edge_norm`: the conditions on `m² + n²` of `IsMSS.edge_norm`.
-/

namespace MagicSquareOfSquares

/-- If `m 0 1 ≠ m 1 0`, `m 2 1 ≠ m 1 2` and `m 0 1 ≠ m 1 2`, one of the four edge entries is
larger than `250000`. -/
theorem IsMagicSquareOfSquares.edge_gt_of_ne {m : Fin 3 → Fin 3 → ℕ} {t : ℕ}
    (H : IsMagicSquareOfSquares m t) (h1 : m 0 1 ≠ m 1 0) (h2 : m 2 1 ≠ m 1 2)
    (h3 : m 0 1 ≠ m 1 2) :
    250000 < m 0 1 ∨ 250000 < m 1 0 ∨ 250000 < m 1 2 ∨ 250000 < m 2 1 := by
  obtain ⟨r, hr, hs⟩ := H.exists_roots
  have ne : ∀ i j i' j' : Fin 3, m i j ≠ m i' j' → (r i j : ℤ) ^ 2 ≠ (r i' j' : ℤ) ^ 2 := by
    intro i j i' j' hij h
    apply hij
    rw [hr i j, hr i' j']
    exact_mod_cast h
  have key : ∀ i j, 500 < |(r i j : ℤ)| → 250000 < m i j := by
    intro i j h
    rw [abs_of_nonneg (Int.natCast_nonneg _)] at h
    have h' : 500 < r i j := by exact_mod_cast h
    rw [hr i j]
    nlinarith
  rcases IsMSS.edge_root_gt_500 hs (ne 0 1 1 0 h1) (ne 2 1 1 2 h2) (ne 0 1 1 2 h3) with
    h | h | h | h
  · exact Or.inl (key _ _ h)
  · exact Or.inr (Or.inl (key _ _ h))
  · exact Or.inr (Or.inr (Or.inl (key _ _ h)))
  · exact Or.inr (Or.inr (Or.inr (key _ _ h)))

/-- **In a 3×3 magic square of nine distinct squares, one of the four edge entries is larger
than `250000 = 500²`.** Published computer searches exclude far larger squares; this is only
the bound that follows from the theorems of this repository. -/
theorem IsMagicSquareOfSquares.edge_gt {m : Fin 3 → Fin 3 → ℕ} {t : ℕ}
    (H : IsMagicSquareOfSquares m t) (hinj : m.Injective2) :
    250000 < m 0 1 ∨ 250000 < m 1 0 ∨ 250000 < m 1 2 ∨ 250000 < m 2 1 :=
  H.edge_gt_of_ne (fun h => absurd (hinj h).1 (by decide)) (fun h => absurd (hinj h).1 (by decide))
    (fun h => absurd (hinj h).1 (by decide))

/-- The conditions on `μ² + ν²` for the matrix form: if `m 0 1 ≠ m 1 0` and `m 2 1 ≠ m 1 2`, the
entries have roots `r` such that, with `μ : ν` the ratio `(r 0 1 + r 1 0) : (r 2 1 + r 1 2)` in
lowest terms, `μ² + ν²` is a square modulo `7`, `17`, `31`, `47` and every prime `≡ 5 (mod 8)`
divides it to an even power. -/
theorem IsMagicSquareOfSquares.edge_norm {m : Fin 3 → Fin 3 → ℕ} {t : ℕ}
    (H : IsMagicSquareOfSquares m t) (h1 : m 0 1 ≠ m 1 0) (h2 : m 2 1 ≠ m 1 2) :
    ∃ r : Fin 3 → Fin 3 → ℕ, (∀ i j, m i j = r i j ^ 2) ∧ ∃ μ ν : ℕ, Nat.Coprime μ ν ∧
      r 0 1 + r 1 0 = μ * Nat.gcd (r 0 1 + r 1 0) (r 2 1 + r 1 2) ∧
      r 2 1 + r 1 2 = ν * Nat.gcd (r 0 1 + r 1 0) (r 2 1 + r 1 2) ∧
      (∀ p : ℕ, p = 7 ∨ p = 17 ∨ p = 31 ∨ p = 47 → IsSquare (((μ ^ 2 + ν ^ 2 : ℕ) : ZMod p))) ∧
      ∀ p : ℕ, p.Prime → p % 8 = 5 → Even (padicValNat p (μ * μ + ν * ν)) := by
  obtain ⟨r, hr, hs⟩ := H.exists_roots
  have hbd : (r 0 1 : ℤ) ≠ (r 1 0 : ℤ) := by
    intro h
    apply h1
    have h' : r 0 1 = r 1 0 := by exact_mod_cast h
    rw [hr 0 1, hr 1 0, h']
  have hne : (r 2 1 : ℤ) + (r 1 2 : ℤ) ≠ 0 := by
    intro h
    have h' : r 2 1 + r 1 2 = 0 := by exact_mod_cast h
    have ha : r 2 1 = 0 := by omega
    have hb : r 1 2 = 0 := by omega
    apply h2
    rw [hr 2 1, hr 1 2, ha, hb]
  obtain ⟨μ, ν, hcop, e1, e2, hq, h5⟩ := hs.edge_norm hne hbd
  have n1 : ((r 0 1 : ℤ) + (r 1 0 : ℤ)).natAbs = r 0 1 + r 1 0 := by omega
  have n2 : ((r 2 1 : ℤ) + (r 1 2 : ℤ)).natAbs = r 2 1 + r 1 2 := by omega
  rw [n1, n2] at e1 e2
  exact ⟨r, hr, μ, ν, hcop, e1, e2, hq, h5⟩

end MagicSquareOfSquares
