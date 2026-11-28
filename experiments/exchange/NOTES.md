# The exchange symmetry `(m, n, u, v) ↔ (u, v, m, n)` (run 10, 27 November 2026, UTC)

Written by an AI model (Claude); not read by a person. `MagicSquareOfSquares/Exchange.lean` is from a separate AI
research session (Claude, Fable 5.2), committed unchanged; `ExchangeCensus.lean` and the check below are the
conductor's. The session's note (27 November 2026, section 1) reported the symmetry as possibly new; it is not
committed (its other sections rest on scripts this repository did not receive). This part of run 10 has had no
separate review apart from the pre-publication referee (`reviews/referee_run10.md`). **No count changes.**

## What was already here

For a magic square of squares the second ratio was already in the library: with roots `a … i`,
`(b + d) : (h + f) = m : n` and `(b + d) : (f − h) = u : v` (`IsMSS.edge_params`, `IsMSS.centre_two_ratios`;
`G1_FULL_SYSTEM.md`, "The centre" and "The two ratios together"; `PAPER.md`, Theorem 3), and negating `h` gives a
magic square whose ratio is `u : v`, so the census theorems apply to the second ratio too. In the G1 family `u : v` is
the parameter `t`.

## What is new (tier A)

- `G1.fullSquares_exchange` (`Exchange.lean`): for natural numbers with `n, v ≠ 0`,
  `FullSquares m n (u/v) ↔ FullSquares u v (m/n)`, for the library's rational five-form system (`G1Full.lean`), with
  no hypothesis of `IsMSS`. Each of `F0, F4, F8, L, U` at `(m, n, u/v)` times `v²` is the same form at `(u, v, m/n)`
  times `n²`. The file also has the corner pairs (`P0_exchange`, `PL_exchange`: coordinates swapped; `P8_exchange`,
  `PU_exchange`: unchanged), their norms, `F0 + F8 = (m² + n²)(u² + v²)`, `L + U = 4(m² + n²)(u² + v²)`, and
  `fiveSquaresAt_exchange` for the integer forms (a separate statement, not used for `fullSquares_exchange`).
- `CoverDescent.mem_exceptions_of_fullSquares` (`ExchangeCensus.lean`): for any `m` and `n ≠ 0`, if
  `FullSquares m n (u/v)` with `u : v` coprime and `0 < v < u ≤ 1000`, then `(u, v)` is one of the 306 pairs left by the
  census in Lean. It covers `t = u/v > 1` only. Since run 12 (after the release candidate `c8b677a`, on an outside
  review's suggestion) the rest is in Lean too: `G1.fullSquares_neg` (`t ↦ −t` swaps `L` and `U`), `G1.fullSquares_inv`
  (for `t ≠ 0`, `t ↦ 1/t` swaps `F0` and `F8`, each form divided by `t²`), `CoverDescent.exists_ratio_of_ne` (every
  `t ∉ {0, ±1}` is `±(u/v)^{±1}` with `0 < v < u` coprime) and `CoverDescent.fullSquares_param_mem_exceptions` (that
  pair, if `u ≤ 1000`, is one of the 306). Every instance of its hypotheses found by the referee has
  `m/n ∈ {0, 1}`: the exchanges of the trivial points `t = 0, 1` of a ratio, e.g. `(m, n, u, v) = (1, 1, 4, 3)`.

Ring identities and the census theorems; standard axioms (the axiom audit, `scripts/Audit.lean`, listed 1,681
theorems at run 10 and lists 1,697 since run 12).

The nine cells are the same numbers up to the factor `v²/n²`, in the same arrangement (the referee's check in
homogeneous coordinates: the matrices of entries are equal, the roots differ by `h ↦ −h`). So the exchange reads the
same square with the parameter and the edge ratio swapped.

## A check at the nine open ratios (tier D; `open_pairs.py`, output `open_pairs.txt`)

At each open ratio `m : n` and `t = ±(u/v)^{±1}` with `u : v` one of the nine open ratios (9 · 9 · 4 = 324 cases),
exact arithmetic: **at most one of the five forms is a square** (164 cases none, 160 cases exactly one, always `F4`:
one ratio of square type and the other of twice-square type, `F4 = (m² + n²)(u² + v²)/(2v²)`; 5 square-type and 4
twice-square-type open ratios, 2 · 5 · 4 · 4 = 160).

**Consequence (prose).** Every open ratio has a degenerate solution (`t = ±1` or `t = 0, ∞`), whose second ratio is
`1 : 1` or `0 : 1`. For a solution with `t ∉ {0, ±1, ∞}` at one of the nine open ratios, normalise `t = ±(u/v)^{±1}`
with `u > v ≥ 1` coprime (Lean since run 12, as above); if `u ≤ 1000`, then `(u, v)` is not outside the 425 pairs nor
among the 119 (Lean, as above), not among the 297 closed by the census as counted (tier B, including the 8 counted apart: there "closed" means
that every solution has `t' ∈ {0, ±1, ∞}`, and `t' = m/n` is not), and not open (the check). So `u > 1000`.

Rerun: `python3 -I experiments/exchange/open_pairs.py | diff - experiments/exchange/open_pairs.txt` (CI job `quick`).
