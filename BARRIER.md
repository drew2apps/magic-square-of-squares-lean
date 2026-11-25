# Why no local test can exclude a ratio of the two degenerate types

Written by an AI model on 25 November 2026 (UTC), after an outside review that Drew passed on (its author is not
recorded here). Read by nobody. Elementary; the facts
are classical, and the repository already says that "no congruence can help" for these ratios. What this file
adds is the consequence for two methods suggested by that review: a Brauer–Manin obstruction and a sieve with
all three line-pair pencils at once.

**Proposition (computed for m ≤ 1000, with a one-line proof for every ratio).** Let `m > n ≥ 1` be coprime.

1. If `m² + n² = k²`, then `u = v` (that is, `t = 1`) gives `b = f`, `d = −h` in Theorem 1 of `PAPER.md`, and the
   nine numbers of the square are all squares: a 3×3 magic square of squares with edge ratio `m : n` and
   repeated entries.
2. If `m² + n² = 2k²`, then `u = 0` (`t = 0`) gives `b = −d`, `f = −h`, and again a magic square of squares with
   repeated entries.

*Proof.* The squares are, for `m² + n² = k²` and for `m² + n² = 2k²`,

    k²        (m + n)²  (m − n)²            m²  n²  k²
    (m − n)²  k²        (m + n)²            n²  k²  m²
    (m + n)²  (m − n)²  k²                  k²  m²  n²

Each is a cyclic arrangement of three squares in arithmetic progression (`(m − n)² + (m + n)² = 2k²`, resp.
`n² + m² = 2k²`), so every row, column and diagonal sums to `3k²`. With the roots `b = m + n`, `d = m − n`,
`f = m + n`, `h = n − m`, resp. `b = n`, `d = −n`, `f = m`, `h = −m`, both lie on the curve `C(m : n)` of
`G1_FULL_SYSTEM.md` (cut out on the surface of magic squares of squares by `n(b + d) = m(h + f)` and
`m(b − d) = n(f − h)`): in the coordinates `b + d = mu`, `h + f = nu`, `b − d = nv`, `f − h = mv` of Theorem 1 of
`PAPER.md` they are the points `u = v = 2` and `u = 0`, `v = 2`. Theorem 1 assumes `h + f ≠ 0`, and the second
point has `b + d = h + f = 0`: it is a rational point of `C(m : n)`, not a square of ratio `m : n` in the sense of
Theorem 7 of `CLASS.md`. ∎

`scripts/degenerate_squares.py` builds the square for every such ratio with `m ≤ 1000` and checks exactly the
eight line sums, the nine squares and the four equations above: 179 of the first type and 127 of the
second, the 306 pairs that Theorem 7 of `CLASS.md` leaves, no failure (`degenerate_squares.txt`, compared in the CI job `quick`).

**Consequence.** For such a ratio the curve `C(m : n)` has a rational point, and so does the fibre of G1 over
`m : n` (at `t = 1`, resp. `t = 0`). So:

- every congruence and every local condition, at any finite set of primes and to any precision, is satisfied;
- the Brauer–Manin set of that curve is not empty, since it contains the rational point; a Brauer class can
  only help together with the Mordell–Weil group of a Jacobian (a Mordell–Weil sieve, ratio by ratio), not on
  its own;
- the three line-pair pencils of the degenerate square are those of one square, so a joint condition on the
  three edge ratios is satisfied too; a joint sieve cannot exclude these ratios either.

What can exclude such a ratio is a global argument that separates the degenerate points from the others on
one curve: Demjanenko–Manin (`G1_FULL_SYSTEM.md`, Result 4), Chabauty, or a Mordell–Weil sieve with a descent
(`LAST_PAIR.md`, Theorem A). Each needs the rank of an elliptic curve that depends on the ratio, so none of
them is uniform on the two families as it stands. A uniform statement on an infinite subfamily would need a
rank bound proved for the whole subfamily, for instance a 2-Selmer group that is provably small there.

This does not touch the ratios that are not of the two types (Theorem 7 of `CLASS.md` removes them to 1000).

## Why the curves of `F0·F8·GL` close no twice-square ratio

`check_sigma.py` reports that none of the 127 ratios to 1000 with `m² + n²` twice a square is closed by its own
curves `E⁺`, `E⁻`. The reason is computed with exact rationals in `scripts/g1_full/sigma_twice.py` (output
`sigma_twice.txt`): on both curves of every one of the 127, the point with `z = ±2` (from `t = ±1`, which exists for
every ratio) and the point `X = 0` (from `t = 0`, which exists exactly when `2(m² + n²)` is a square) are
independent, so both curves have rank at least 2 and Demjanenko–Manin with one generator cannot apply. This covers
the 127 ratios to 1000 only. For four of them, (313,103), (503,329), (593,257) and (713,503),
`dm_sigma_306.jsonl` records `E⁻` with the rank bounds `[1, 3]`; the rank is in fact 2 or 3.

For the square type the point at `t = 1` is the same as the universal one (`z = 2`), so only one point is forced,
and Demjanenko–Manin works whenever the rank is exactly 1. A proof for every twice-square ratio would show that
the two sections are independent over the function field of the conic `x² + y² = 2z²` and specialise
injectively; that is not done.
