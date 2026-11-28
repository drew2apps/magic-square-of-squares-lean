import MagicSquareOfSquares.Exchange
import MagicSquareOfSquares.CoverDescentTable

/-!
# The second ratio of a solution is in the census list too

`G1.fullSquares_exchange` (`Exchange.lean`) with the census of this library: if the five cells of the G1 family are
squares at the ratio `m : n` and the parameter `t = u/v`, where `u : v` is a coprime pair with `v < u ≤ 1000`, then
`(u, v)` is one of the 425 pairs of `G1.exceptions1000` (`no_solution_le_1000`) and, more precisely, one of the 306
pairs of `CoverDescent.exceptions` (`CoverDescent.no_fullSquares_of_mem_ratios`). So the parameter of a solution
obeys the same restrictions as the ratio. Nothing is assumed about `m : n`.

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
