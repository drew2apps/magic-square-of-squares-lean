# The Cassels–Tate pairing on the pair curves, in closed form where it has one

> **Read `NOTES.md` first.** This is the session's report, committed unchanged apart from this line; where they differ, `NOTES.md` (after review) holds.

27 November 2026 (UTC). Written by an AI model (Claude, Fable 5.2); not read by a person; not reviewed
by a separate agent. Tier D; no count changes. Input: the pilot `conductor-t10` (its four output
files, 8,930 kill lines). PARI is used for `p`-adic square roots and Hilbert symbols.

## Result in one paragraph

Two of the six pair curves, `E_08` and `E_LU`, are Frey curves of the Pythagorean triple
`(m² − n², 2mn, m² + n²)`: all three differences of their 2-torsion are squares. For these two the
class of a tuple can be moved to the form `(1, δ, δ)`, the three conics of Cassels' formula have
rational points that do not depend on the ratio, and the pairing becomes a product of Hilbert
symbols with a closed form at every odd prime. It reproduces all 2,154 kill lines of the pilot on
these curves, and its closed forms agree with the numerical local terms at all 21,838 odd local
terms. The symbols are Legendre symbols, and quartic residue symbols at the primes where the
first component of `b` has odd valuation. For the four mixed curves (6,776 kill lines) there is no
such formula: one conic depends on the tuple and has to be solved, a norm equation of the type that
gives Rédei symbols.

## 1. The curves and the classes (checked on every line of the pilot)

`S = m² + n²`, `J± = m² ± 2mn − n²`, `R± = (J±² + S²)/2`, so `R₊R₋ = m⁸ + 14m⁴n⁴ + n⁸`.
Coordinates of a class: `(x − e_r, x − e_m, x − e_*)`, product a square; the pilot prints the first
two (found by testing local solubility, `assign.py`; for `E_LU` its header names `e* = e3` but the
printed second component is at `e2`).

| curve | model | class of a tuple `τ = (d0, d8, dL, dU)` |
|---|---|---|
| `E_08` | `y² = x(x − A²)(x − C²)`, `A = m² − n²`, `C = S`, `C² − A² = (2mn)²` | `(1, δ, δ)`, `δ = mn·d0·d8` |
| `E_LU` | the same with `A = 2mn`, `C = S`, `C² − A² = (m² − n²)²` | `(2S, 2·dL·dU, ·)`; times the image `(2S, 2(m² − n²), ·)` of a point of order 4 it is `(1, δ, δ)`, `δ = (m² − n²)·dL·dU` |
| `E_0L, E_0U, E_8L, E_8U` | `y² = x(x − 2R₊)(x − 2R₋)` (one curve) | `(2, κ·d_j·d_k, ·)`, `κ = mn(m² − n²)` up to a factor 2 that depends on the curve and on parity |

1,280 + 874 lines confirm the first two rows exactly; the models of the mixed curves are confirmed
on all 804 curves, their class constant on a sample only. The sign of `δ` is free: `(1, −1, −1)` is
the image of a 2-torsion point.

## 2. The formula for `E_08` and `E_LU`

Let `a = (1, δ, δ)` and `b = (b_r, b_m, b_*) ∈ Sel²`. The 2-covering of `a` is
`u² = x`, `δ·u_m² = x − A²`, `δ·u_*² = x − C²`. Its three conics

    u² − δu_m² = A²,      u² − δu_*² = C²,      δ(u_m² − u_*²) = C² − A²

have the rational points `(A, 0)`, `(C, 0)` and the point at infinity `u_m = u_*`, with tangent
forms `u − A`, `u − C`, `u_m − u_*`. Cassels' formula (1998; equal to the Cassels–Tate pairing by
Fisher–Schaefer–Stoll) then reads

    ⟨a, b⟩ = ∏_v (u − A, b_*)_v · (u − C, b_m)_v · (u_m − u_*, b_r)_v ,

at any local point of the 2-covering at each place `v`. The real place contributes 1. Which
component of `b` goes with which form was fixed by the data (the other assignment fails on 222
lines).

**Closed form at an odd prime `p`** (derived by choosing a local point in each case; then checked):

| case | local term |
|---|---|
| `p` divides `C − A` | 1 |
| `p` divides `C + A` | `(δ, b_*)_p` |
| otherwise, `δ` a square in `ℚ_p` | `(2√δ, b_r)_p` |
| otherwise, `−δ` a square in `ℚ_p` | `(−A, b_*)_p · (−C, b_m)_p · ((A − C)·√(−δ), b_r)_p` |

For `E_08`: `C − A = 2n²`, `C + A = 2m²`. For `E_LU`: `C − A = (m − n)²`, `C + A = (m + n)²`.
Every odd bad prime of every kill line fell in one of the four cases. By hand: at `p | A` a class
`(1, δ, δ)` is locally soluble only if `p ∤ δ` and `δ` or `−δ` is a square; at `p | C` only if `δ` is
a square.

When `b_r = 1` (so `b = (1, β, β)`; 668 of the 1,280 lines of `E_08`) no square root is left:

    ⟨(1, δ, δ), (1, β, β)⟩ = t₂ · ∏_{p | C + A} (δ, β)_p · ∏_{p | A, δ not a square mod p} (AC, β)_p .

When `b_r` has odd valuation at `p` (possible only at primes of `A·C`), the term `(2√δ, b_r)_p` is
the quartic residue symbol of `δ` at `p` times `(2/p)`. This is the same shape as for the congruent
number curves (Monsky's matrix for Sel², quartic symbols for the pairing).

**Not in closed form here:** the term `t₂` at 2 (computed numerically; it is −1 on 156 + 122
lines), and Sel² itself as a Monsky-type matrix (only the three necessary conditions above).

## 3. Check against the pilot (`ct08.py`, `check_local.py`; outputs `check_local_08.txt`, `_LU.txt`)

- `E_08`: 1,280 of 1,280 kill lines give −1. Odd local terms against the table: 2,644 + 2,972 +
  4,906 + 2,462 = 12,984, all agree.
- `E_LU`: 874 of 874 give −1. Odd local terms: 1,676 + 2,528 + 3,210 + 1,440 = 8,854, all agree.
- Controls (`controls.py`, 150 lines drawn at random from the two curves): the formula gives 1 for
  `b = a` (alternating), for the image `(1, −1, −1)` of the 2-torsion and for the images of both
  points of order 4, and still −1 for the killing `b` times the image of a point of order 4.
- Not tested: tuples that are not killed (the files hold kill lines only), and `b` over a whole
  basis of Sel². So the check shows that the formula finds every −1 of the pilot and is 1 where it
  must be on torsion, not yet that it agrees with `rank_second` on every pair of classes.

## 4. The four mixed curves

Class `(2, c, 2c)` with `c = κ·D`, `D = d_j·d_k`; differences `2R₊`, `2R₋`, `8mn(m² − n²)`. The conic
`2u² − c·u_*² = 2R` has the rational point `(J, 2)` when `c = mn(m² − n²)`, because
`J² − S² = 4mn(m² − n²)`: that is the class of the generic point over `t = ±m/n, ±n/m`. For a
general tuple the conics are `2u² − κD·w² = 2R±` and `κD(2u_m² − u_*²) = ±8mn(m² − n²)`, and their
points depend on `D`. The last one asks to write `D` (times a unit) as `x² − 2y²`, which is where
the class of a prime modulo 8 enters. The pairing on these curves is therefore of Rédei-symbol
type (Smith; Morgan–Smith), not a product of Legendre symbols of the primes of the ratio. This is
what I expect, not something I proved.

## 5. The question about tuples with trivial pairing

Not answered. What the formula gives: on `E_08` the tuple enters only through `δ = mn·d0·d8`, on
`E_LU` only through `(m² − n²)·dL·dU`, so these two curves see two of the six products. A tuple
escapes `E_08` when every `b` of Sel² has an even number of −1 terms, which for `b_r = 1` is a
linear condition on Legendre symbols of `δ` against the primes of `m` and of `m² − n²`. At
(5623,1017) Sel² of `E_08` is the torsion image, so there is no `b` at all; that is a rank-0
statement, not a pairing statement, as the pilot says.

## 6. Literature (searched 27 November 2026; none gives this formula for this family)

- J. W. S. Cassels, Second descents for elliptic curves, J. reine angew. Math. 494 (1998): the
  pairing by conics and Hilbert symbols.
- T. Fisher, E. Schaefer, M. Stoll, The yoga of the Cassels–Tate pairing, LMS J. Comput. Math. 13
  (2010): Cassels' pairing is the Cassels–Tate pairing. H. Shukla, M. Stoll, arXiv:2302.01640.
- T. Fisher, a formula for the pairing on Sel² from the invariant theory of binary quartics that
  solves no conics (Research in Number Theory, 2022): the tool to try on the mixed curves.
- P. Monsky, appendix to D. R. Heath-Brown, The size of Selmer groups for the congruent number
  problem II (1994): Sel² as the kernel of a matrix of Legendre symbols.
- Z. Wang, arXiv:1511.03810 and arXiv:1703.05909: for congruent number curves the pairing is a
  quartic residue symbol, tied to the 8-rank of a class group.
- arXiv:2303.05059: compares Sel² and the Cassels pairing of twists of `y² = x(x − e1)(x + e2)` and
  of `y² = x(x − e1·a²)(x + e2·b²)` when `e1a² + e2b² + e3c² = 0`. With `e1 = e2 = 1` the second curve
  is our Frey curve and the first is the congruent number curve. Its abstract asks for odd `a, b, c`,
  which a Pythagorean triple is not; it has to be read before it is used. I read the abstract only.
- B. Naskręcki, arXiv:1210.6933: ranks in the family `y² = x(x − a²)(x − b²)` of Pythagorean triples.

## Reproduce

    python3 ct08.py 08 ; python3 ct08.py LU      # pairing by the explicit formula on the kill lines
    python3 check_local.py 08 ; python3 check_local.py LU   # the table of section 2
    python3 assign.py                            # which roots the printed components belong to

(`parse.py` reads `t10_*_[ab].txt` from the current directory.)
