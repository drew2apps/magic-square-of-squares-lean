# A closed formula for the Cassels–Tate pairing on two pair curves (run 10, 27 November 2026, UTC)

Written by an AI model (Claude); not read by a person. `PAIRING_FORMULA.md` and its scripts (`ct08.py`,
`check_local.py`, `controls.py`, `assign.py`, `parse.py`, `natural_points_627.py`, outputs `check_local_08.txt`,
`check_local_LU.txt`) are from a separate AI research session (Claude, Fable 5.2), committed unchanged apart from
the banner line of the report. Reviewed by a separate hostile AI agent (z12; `reviews/review_z12.md`, its checks in
`../referee_z12/`), not by a person. **Experiment, not counted; tier D of `LEDGER.md`.** Nothing here changes the
census.

## The claim, in one paragraph

Of the six pair curves of the pilot `../conductor-t10/`, two are Frey curves of the Pythagorean triple
`(m² − n², 2mn, m² + n²)`: `E_08: y² = x(x − A²)(x − C²)` with `A = m² − n²`, `C = m² + n²`, and `E_LU`, the same
with `A = 2mn`. The class of a tuple `τ` is `(1, δ, δ)` with `δ = mn·d0·d8` on `E_08` and, after the image of a point of
order 4, `δ = (m² − n²)·dL·dU` on `E_LU`. The three conics of Cassels' formula then have points that do not depend
on the tuple, so the pairing is a product of Hilbert symbols, `∏_v (u − A, b_*)_v (u − C, b_m)_v (u_m − u_*, b_r)_v`,
with a closed form at each odd prime (four cases, `PAIRING_FORMULA.md` §2); the term at 2 is computed numerically.
On these two curves the tuple enters only through `δ`, by Legendre symbols, and by quartic residue symbols at primes
where `b_r` has odd valuation.

## What the review found (z12)

- **The formula agrees with the repository's complete 2-descent and pairing** (`rank_second.full_selmer2`,
  `DFull.pairing`, through `t10_pairct.py`) **on every case tested**: 28,918 distinct (ratio, curve, `a`, `b`), 0
  disagreements, on 146 of the pilot's 201 ratios for `E_08` and 133 for `E_LU`, (5623,1017) and (3320,2303)
  included; `a` = the class of every residual tuple, killed or not, both signs (1,736 + 1,540), and every `(1, δ, δ)`
  in `Sel²` (16,692 + 8,950); `b` over a basis of `Sel²`; over all of `Sel²` the formula is the linear extension of
  rank_second's values (117,016 checks). The session had tested only the pilot's kill lines. The sample is not all
  201 ratios because the run was stopped at about 90 minutes, not by choice.
- **The 2-adic term** does not depend on the random local point (8 seeds; PARI refuses an under-precise Hilbert symbol
  rather than return a sign).
- **The assignment of the components of `b` to the three tangent forms** is not a parameter fitted to the data, as the
  report says ("fixed by the data"): it follows from the Kummer map (the reviewer's derivation, not compared with
  Cassels' or Fisher–Schaefer–Stoll's indexing), it is the one `rank_second` uses, and the other assignment is not
  locally constant when `b_r ≠ 1`, so it is not a pairing at all.

## Corrections to `PAIRING_FORMULA.md`

1. **Hypotheses.** The formula needs `gcd(m, n) = 1` (with a common factor the case "`p | C − A`: 1" fails on 34 of 156
   tests, `../referee_z12/zlocal.txt`) and `A > 0`, i.e. `m > n` on `E_08`: "the real place contributes 1" is true
   only then (with the conic point `(−A, 0)` the formula is wrong on the 323 of 3,912 triples with `δ < 0`,
   `b_m < 0`). Every pilot ratio is coprime with `m > n`, so no number of the report changes.
2. **Counts are repetitions.** The 1,280 + 874 = 2,154 kill lines are 133 + 111 = **244 distinct** (ratio, class, `b`)
   on 121 + 94 ratios; the 21,838 odd local terms are **2,415 distinct**; "the other assignment fails on 222 lines" is
   `E_08` only (15 of 133 distinct; on `E_LU` 112 lines, 15 of 111); the 150 controls are drawn from repeated lines.
3. **Cases 1 and 2 of the odd-prime table** (`p | C ∓ A`) are checked numerically only (8,476 pairs `(a, b)` at
   coprime inputs in the committed `../referee_z12/zlocal.txt`, 6,828 with `m > n` and 1,648 with `m < n`; the review's
   text says 9,412; all four cases, `p = 3` included), not derived by hand; cases 3 and 4 and the `b_r = 1` product
   formula were re-derived by the reviewer. So "closed form at every odd prime" is a tested statement, not a proved one.
4. **The four mixed curves** (6,776 of the 8,930 kill lines): "there is no such formula" (the report's first paragraph)
   and "therefore of Rédei-symbol type" (§4) are expectations. That one conic depends on the tuple does not rule out a
   closed form (Fisher 2022 gives the pairing without solving conics).
5. **Literature.** arXiv:2303.05059 assumes `e1 + e2 + e3 = 0`; with `e1 = e2 = 1` its second curve has root
   differences `a², b², 2c²`, which is **not** our Frey curve (the report says it is). Naskręcki (arXiv:1210.6933)
   treats `y² = x(x − a²)(x − b²)` with `a, b` the two legs, a different family (a leg and the hypotenuse here). One
   search by the reviewer found no source for this formula; that is not conclusive. It is a direct specialisation of
   Cassels' formula; what is new here is the four local cases.

## What it means for the pilot

On `E_08` and `E_LU` the kills of `../conductor-t10/` (2,154 of 8,930 lines) are reproduced by symbols of `δ`
against the primes of `m`, `n`, `m² ± n²`: tested, with the gaps above. The other 6,776 lines, on the four mixed
curves, have no known formula. A tuple escapes `E_08` when every `b` of `Sel²` pairs to 1 with `(1, δ, δ)`, a linear
condition on these symbols; which tuples escape all six curves is not answered.

## Rerun

From `experiments/conductor-t10/` (the scripts read the pilot's `t10_*_[ab].txt` from the current folder and import
each other from their own folder, so run them without `-I`):

    python3 ../research-pairing/ct08.py 08 ; python3 ../research-pairing/ct08.py LU      # the kill lines, about 3 s each
    python3 ../research-pairing/check_local.py 08 | diff - ../research-pairing/check_local_08.txt    # CI: descents
    python3 ../research-pairing/check_local.py LU | diff - ../research-pairing/check_local_LU.txt    # CI: descents
    python3 ../research-pairing/controls.py

`natural_points_627.py` is the session's search for points of small height on the four quartics `L·L·M·M` of
(627,364) over `ℚ(i)` (1,910 candidates; square values only at `s = ±i`; not reviewed; not used here; the session's
suggested search on explicit 4-coverings was not done).
