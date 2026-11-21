# Audit of the G1 closures in Kepler71/magic-square-of-squares

[Kepler71/magic-square-of-squares](https://github.com/Kepler71/magic-square-of-squares) studies a two-parameter family of 3×3 magic squares called G1. For coprime `1 ≤ n < m ≤ 20` (127 pairs) it reports that 120 pairs cannot give a magic square of nine distinct squares: 62 because an elliptic curve has rank 0, and 58 because of congruences. Its handoff note says these 120 closures were inherited without a consolidated audit.

This is an independent check of those 120 closures. It was made from the statements in `RESULTS_FOR_CLAUDE_2026-09-12_G1.md` and `SUMMARY_2026-09-12.md` (read on 21 November 2026). The repository's code and logs were not read or run.

## The claim that is checked

The three diagonal cells of the family are

```
F0 = m² + n²t²,   F4 = (m² + n²)(t² + 1)/2,   F8 = n² + m²t².
```

A magic square of nine squares in the family needs a rational `t` for which all three are squares. With `t = U/V` in lowest terms this needs the integers `m²V² + n²U²`, `2(m² + n²)(U² + V²)`, `n²V² + m²U²` to be squares.

## Results

**Reproduced.** Every number below matches the Kepler71 report.

| Item | Kepler71 | This audit |
|---|---|---|
| Torsion of `E_B: Y² = X(X − (m²−n²)²)(X − (m²+n²)²)` | `Z/2 × Z/4` for all 127 pairs | same |
| Rank bounds from PARI `ellrank` | 62 × [0,0], 60 × [1,1], 4 × [2,2], 1 × [0,2] | same |
| Pairs without rank 0 that are excluded by a congruence | 58 | 58 |
| First obstruction among those 58, in the order 32, 25, 7, 23 | 32, 7, 15, 4 | 32, 7, 15, 4 |
| Pairs left open | (11,4), (13,8), (15,1), (15,8), (16,5), (19,5), (19,16) | the same seven |
| The 2-adic criterion (lemma 1) | stated for all coprime `m, n` | agrees with the scan on all 127 pairs |

The ranks were recomputed with the same software (PARI 2.17.2), so that row is a reproduction, not an independent proof. The congruence scan is independent code.

**New.**

1. **56 of the 62 rank-0 pairs are also excluded by a congruence.** So 114 of the 120 closures need no rank computation and no computer algebra system. Only six closures depend on the rank being 0: (4,3), (7,1), (12,5), (17,7), (20,1), (20,7).
2. This covers (13,3) and (19,11). For these two the report says that eclib gave only the bound [0,2] and that the closure rested on PARI alone. Both are excluded modulo 7.
3. The pair (12,7) is excluded modulo 47, and by no smaller modulus that was tested. The report scanned odd primes up to 31.
4. **All 114 congruence exclusions are checked by the Lean kernel.** `MagicSquareOfSquares/G1.lean` proves that a finite check modulo `M` rules out every rational `t`. `MagicSquareOfSquares/G1Certificates.lean` lists the 114 triples `(m, n, M)` and the kernel evaluates each check. The theorems depend only on `propext`, `Classical.choice` and `Quot.sound`.

The moduli used are 7 (68 pairs), 8 (23), 23 (8), 25 (9), 32 (5) and 47 (1).

## Beyond `m, n ≤ 20`: whole remainder classes

The check depends only on the remainders of `m` and `n`, so each obstruction rules out infinitely many pairs at once. These statements hold for all natural numbers `m`, `n` in the class and are checked by Lean (`G1Classes.lean`, `G1Primes.lean`). `DiagonalSquares m n t` says that `F0`, `F4`, `F8` are squares of rationals; each theorem concludes `¬DiagonalSquares m n t` for every rational `t`.

| Theorem | Condition on `m`, `n` | How |
|---|---|---|
| `no_point_of_two_mod_four` | one is `2` modulo `4`, the other is odd | table modulo 8 |
| `no_point_of_odd` | both odd, one is `±1` and the other `±3` modulo `8` | table modulo 32 |
| `no_point_of_five_mod_eight` | a prime `p ≡ 5 (mod 8)` divides `m² + n²` to an odd power, `p ∤ n` | general proof, any such `p` |
| `no_point_mod_seven` | `m ≡ ±2n` or `n ≡ ±2m` modulo 7, not both divisible by 7 | table modulo 7 |
| `no_point_mod_twenty_three` | the same modulo 23 | table modulo 23 |
| `no_point_mod_forty_seven` | `m ≡ ±5n` or `n ≡ ±5m` modulo 47, not both divisible by 47 | table modulo 47 |

The first three are the two general lemmas of the Kepler71 notes (the 2-adic criterion and the primes `5 mod 8`), now with Lean proofs. The classes modulo 7, 23 and 47 appear pair by pair in the notes; I did not find them stated as classes there.

For a table modulo a prime power `p^k`, the Lean check uses that every admissible `(U, V)` is a unit multiple of `(1, v)` or of `(u, 1)` with `p | u`, so about `2p^k` cases per class are enough.

## Why no further primes appear

**Claim.** Let `m`, `n` be coprime and let `p` be a prime not dividing `2mn(m² − n²)(m² + n²)`. Then there is an obstruction at `p` exactly when `p` is 7, 23 or 47 and `(m : n)` is in one of the classes in the table.

**Argument (not formalized, not refereed).** For such `p` the equations `a² = m²V² + n²U²`, `b² = 2(m² + n²)(U² + V²)`, `c² = n²V² + m²U²` define a smooth curve over `F_p`. It covers the `(U : V)` line with group `(Z/2)³`, branched over six points, so its genus is 5. By the Weil bound it has at least `p + 1 − 10√p` points over `F_p`, which is positive for `p ≥ 101`. A point modulo `p` lifts to a `p`-adic point by Hensel's lemma, so there is no obstruction modulo any power of `p`. For the odd primes below 101, `scripts/g1_good_primes.py` counts the points for every ratio: there is no point only for the classes in the table (it runs to 400, and the Weil bound holds in every case, with worst ratio 0.854).

**Primes dividing `mn(m² − n²)`, odd.** No obstruction: if `p | mn` take `(U, V) = (1, 1)`, and if `p | m² − n²` take `(U, V) = (1, 0)`. The three numbers are then nonzero squares modulo `p`.

**What is left unproved.** Two converse statements: that the 2-adic table is complete (lemma 1 of the Kepler71 notes says so), and that an odd prime dividing `m² + n²` obstructs only when it is `5 mod 8` with odd exponent. Both agree with every direct test in `scripts/g1_three_kinds.py` (prime powers up to 2809, no disagreement).

If those two converses hold, a coprime pair has a local obstruction exactly when one of the three kinds applies.

## Census for `m, n ≤ 200` and `≤ 1000`

`scripts/g1_census.py` applies the three kinds to all coprime pairs.

| Range | Coprime pairs | Excluded by a congruence | Not excluded |
|---|---|---|---|
| `m ≤ 20` | 127 | 114 (89.8%) | 13 |
| `m ≤ 200` | 12,231 | 11,051 (90.4%) | 1,180 |
| `m ≤ 1000` | 304,191 | 278,873 (91.7%) | 25,318 |

For the 1,180 pairs up to 200, PARI `ellrank` gives these bounds for the rank of `E_B`: [0,0] 252 times, [0,2] 190, [1,1] 550, [1,3] 24, [2,2] 155, [3,3] 9. By the Kepler71 criterion (rank 0 excludes nine distinct squares) 252 more pairs are closed. That leaves **928 pairs up to 200 that no method here excludes**, against 7 up to 20. The congruence exclusions in this table are computed by the script; only the classes in the previous section are Lean theorems.

## Limits

- Everything here concerns the family G1. It says nothing about the open problem outside that family.
- The rank statements rely on PARI and on the Kepler71 criterion that rank 0 of `E_B` excludes nine distinct squares. Neither is checked in Lean.
- The formulas for `F0`, `F4`, `F8` are taken from the Kepler71 notes. The family itself, and the fact that it is magic, are not formalized here.
- The Weil-bound argument above is a standard one written by an AI model. It has not been checked by a number theorist.
- The seven open pairs were not examined.

## Reproduce

```
python3 scripts/g1_audit.py            # scan; recomputes ranks if cypari2 is installed
python3 scripts/g1_good_primes.py      # primes of good reduction, Weil bound
python3 scripts/g1_three_kinds.py      # direct tests at prime powers
python3 scripts/g1_census.py 200 rank  # census; "rank" needs cypari2
lake build                             # includes the certificates and the class theorems
```
