# The full G1 system: congruences, trivial points, and what is left

This note replaces an earlier audit, now in [archive/KEPLER71_G1_AUDIT.md](archive/KEPLER71_G1_AUDIT.md). That file looked only at the three diagonal cells of the G1 family of [Kepler71/magic-square-of-squares](https://github.com/Kepler71/magic-square-of-squares). This one uses all nine cells.

Status: the congruence census for `m ≤ 1000` and the table of obstructed classes are **Lean theorems** (see "What is checked in Lean"). Everything else is **computed, not formalized, not refereed.** The scripts are in `scripts/g1_full/`. Prepared with Claude, an AI model.

## Which numbers are checked by Lean

| Number | Status |
|---|---|
| The 484 obstructed classes at 34 primes | Lean: each class is obstructed (`fullTable_valid`, `maskTable_valid`). That there are no others is computed |
| 425 exceptions for `m ≤ 1000`, 60 for `m ≤ 200` | Lean (`exceptions1000_length`, `exceptions200_length`): every coprime pair outside these lists is excluded |
| 304,191 coprime pairs, hence 303,766 excluded | The count of coprime pairs is by script; the Lean theorem does not need it |
| 306 of the 425 have a trivial point | Computed by the generator; for `m ≤ 200` all 60 do, and that is Lean (`exceptions200_spec`) |
| Some edge root exceeds 500 | Lean (`IsMSS.edge_root_gt_500`) |
| 44 of 60, 183 of 306, 85 of 119 closed by rank 0; 157 open | Computed with PARI. Not in Lean |
| Exactly 34 primes; the Weil-bound argument | Computed and argued. Not in Lean |

## Correction to the earlier audit

The audit was written from two Kepler71 notes dated 12 September 2026. Later notes in the same repository go further, and I had not read them:

- Their later notes report all 127 pairs with `m ≤ 20` as excluded, not 120. For the last one, (15,8), a note of 13 September gives a computational proof and lists its audit as pending; their README still names (15,8) and (11,4) as open. I have not checked either claim.
- The five-form local test below is theirs (`PARALLEL_RESULT_FULL_NINE_126_2026-09-12.md`). They state the forbidden classes modulo 17 and 41 for all `(m : n)`.
- They also use other elliptic quotients of the same curves, and report 1,715 of the 12,231 pairs with `m ≤ 200` closed by a rank-0 criterion.

So the audit's phrase "seven open pairs" describes their state on 12 September only, and class-wide congruence obstructions were not new. What follows was checked against their notes by text search; I did not find it there, but I have not read all 250 files.

## The system

For coprime `m > n ≥ 1` and `t = a/b`, four cells of the square are squares identically: `(ma ± nb)²`, `(na ± mb)²`. The other five must be squares:

```
F0 = m²b² + n²a²        F8 = n²b² + m²a²        F4 = s(a² + b²),  s = (m² + n²)/2
L  = F4 − 2mnab         U  = F4 + 2mnab
```

Two facts used below:

- Every 3×3 magic square of rational squares arises this way (my derivation): the four identically square cells are the edge cells, and any two ways of writing a number as a sum of two squares have the form `(ma+nb)² + (na−mb)² = (ma−nb)² + (na+mb)²`. So G1, taken over all `(m : n)`, is the whole problem, cut into curves.
- `(m, n, t)` and `(m+n, m−n, (t+1)/(t−1))` give the same square, scaled by 4, with `F0, F8` and `L, U` exchanged. Pairs therefore come in couples (divide by 2 when `m`, `n` are both odd).

A pair has a **trivial point** when `m² + n²` is a square (`t = ±1`) or `(m² + n²)/2` is a square (`t = 0, ∞`). These are degenerate squares with repeated cells.

## Result 1: the primes at which the full system can fail modulo `p`

For `p ∤ 2mn(m² − n²)(m² + n²)` the five conditions define a smooth curve of genus 49 over `F_p` (a `(Z/2)⁵` cover of the line branched at ten points; the Kepler71 notes give the same genus and cite Auel–Singer, arXiv:2609.09351). The Weil bound gives at least `p + 1 − 98√p` points, which is positive for `p ≥ 9613`. `classify.c` tests every ratio `(m : n)` for every odd prime below 9700.

**There are exactly 34 primes with an obstructed class, the largest is 607:**

7, 17, 23, 31, 41, 47, 73, 79, 89, 103, 113, 127, 137, 151, 167, 191, 193, 199, 223, 233, 239, 257, 263, 271, 281, 311, 313, 337, 353, 439, 457, 521, 599, 607.

There are 484 obstructed classes in total (`modp_classes.txt`), all at ratios of good reduction. All 34 primes are `±1 mod 8`; I have no proof of why. The classes for 17 and 41 agree with the Kepler71 lists.

## Result 2: census

Local tests used: the 34 primes above, the 2-adic test modulo 32, 64, 128, and powers of the odd primes dividing `mn(m² − n²)(m² + n²)`.

| Range | Coprime pairs | Excluded by a congruence | Not excluded |
|---|---|---|---|
| `m ≤ 20` | 127 | 122 | 5 |
| `m ≤ 200` | 12,231 | 12,171 | 60 |
| `m ≤ 1000` | 304,191 | 303,766 | 425 |

**For `m ≤ 200`, the pairs not excluded by a congruence are exactly the 60 pairs with a trivial point.** Every pair without one fails at some prime.

For `m ≤ 1000` the 425 survivors are the 306 pairs with a trivial point and 119 pairs without. The first of the 119 is (231, 16). For these 119 there is no congruence obstruction among the tests above and no rational point is known.

799 closures were re-checked by a separate brute-force program (`verify_bruteforce.py`), with no mismatch, and no pair with a trivial point is ever reported as obstructed.

## Result 3: the survivors, by elliptic quotients of rank 0

With `w = t²`, a solution gives a rational point with first coordinate `w` on `y² = P(w)` for each product `P` of `w`, `f0 = m² + n²w`, `f8 = n² + m²w`, `g4 = 2(m² + n²)(w + 1)`, `K = (m² + n²)²(w + 1)² − 16m²n²w` (`K` is `L·U` up to a square). Twelve of these that are elliptic curves are used. If one has rank 0, `w` is the coordinate of a torsion point, and the finitely many candidates are tested directly. `w·f0·f8` is the curve `E_B` of the Kepler71 notes up to isogeny; `w·K` is its partner under the symmetry above.

Rank 0 is proved by PARI `ellrank` on the curve or on a curve isogenous to it. No torsion point gave a non-degenerate square.

| Set | Pairs | Closed | Open |
|---|---|---|---|
| trivial point, `m ≤ 200` | 60 | 44 | 16 |
| trivial point, `m ≤ 1000` | 306 | 183 | 123 |
| no trivial point, no congruence, `m ≤ 1000` | 119 | 85 | 34 |

- Of the 85, 42 are closed without `E_B`, by `f0·g4·f8`, `g4·K` or `w·K`.
- 19 of the closures use `L(E,1) ≠ 0` (Kolyvagin) where 2-descent gives only `[0, 2]`: 18 pairs with a trivial point, among them the couple (112,15)~(127,97) with `L(E,1) ≈ 1.026`, and (700,529). Only curves of conductor below `10^14` were tried.

**Open for `m ≤ 200` (16 pairs, 10 up to the symmetry):** (15,8)~(23,7), (21,20)~(41,1), (72,65)~(137,7), (91,60)~(151,31), (99,20)~(119,79), (117,44)~(161,73), (156,133), (176,57), (187,84), (195,28). Kepler71 reports a computational proof for (15,8) (see the correction above); if it stands, nine remain.

**Open for `m ≤ 1000` without a trivial point (34 pairs):** (320,289), (388,287), (423,332), (427,320), (540,239), (541,476), (560,33), (593,527), (609,31), (664,301), (675,101), (724,325), (724,665), (740,539), (747,107), (755,91), (773,192), (779,301), (784,391), (788,487), (820,799), (828,173), (836,745), (840,83), (864,151), (889,320), (899,724), (920,371), (923,757), (935,416), (964,541), (965,363), (965,581), (988,27).

All open pairs are listed in `scripts/g1_full/open_pairs_1000.txt`; `global.py` reproduces the closing curve for each closed pair.

**In total: of the 304,191 pairs with `m ≤ 1000`, 157 are open.**

For the 119 pairs without a trivial point, the bad primes and 2 were also tested to prime powers up to 3·10⁶; no further obstruction appeared.

## Limits

- "Open" means not closed by these methods. It does not mean a solution is likely.
- "Not excluded by a congruence" is not a proof of local solubility: at primes of bad reduction only small prime powers were tested.
- The rank bounds depend on PARI 2.17.2. The `L`-value closures also depend on its numerical `L`-function code.
- The Weil-bound step is a standard argument, written by an AI model and not checked by a number theorist.
- In Lean: only what the next section lists. The completeness of the list of 34 primes and every rank statement are not in Lean.

## Why this sieve does not transfer to the other fibrations

The same surface can be cut into curves in other ways, and I tried the same sieve on three of them. It finds nothing, for a simple reason: every fibre has a rational point, so no congruence can exclude it.

- **Slopes** (the Kepler71 census: cells `1, 1 ± sz, 1 ± rz, 1 ± (s−r)z, 1 ± (s+r)z`). At `z = 0` all nine cells are equal.
- **G2** (the two corner pairs are squares identically; centre and edges must be squares). At `(a, b) = (m + n, m − n)` the five values are `4(m² + n²)²` and `4(m² ± 2mn − n²)²`, each of the last two twice.
- **G3** (one corner pair and one edge pair). `families.py` finds no obstructed class at any prime below 700 and no excluded pair for `m ≤ 200`, the same as for G2.

So G1, where the two edge pairs are the identical squares, is the one fibration of these four in which degenerate points are rare and congruences decide almost every fibre. For the others only global methods (ranks, Chabauty) can work, which is what Kepler71 used for the slopes.

## What is checked in Lean

| Theorem | Statement |
|---|---|
| `G1.no_solution_of_ratio` (`G1Full.lean`) | For a prime `p`: if the finite check fails for the ratio `(r : 1)`, no rational `t` makes the five cells squares for any `(m, n)` with `m ≡ r·n (mod p)`, `p ∤ n` |
| `G1.fullTable_valid` (`G1FullCertificates.lean`) | The check fails for each of the 484 classes at the 34 primes. Evaluated by the Lean kernel |
| `G1.no_solution_le_200` | For coprime `1 ≤ n < m ≤ 200` outside the 60 pairs of `exceptions200`, no rational `t` makes the five cells squares |
| `G1.exceptions200_spec` | Each of the 60 exceptions has `m² + n² = k²` or `m² + n² = 2k²` |
| `G1.no_solution_le_1000` (`G1Census1000.lean`) | The same for coprime `1 ≤ n < m ≤ 1000` outside the 425 pairs of `exceptions1000` (the 306 with a trivial point and the 119 without). The kernel evaluates the census in 20 blocks |
| `IsMSS.fullSquares` (`G1Bridge.lean`) | A magic square of squares with roots `a b c / d e f / g h i` and `(b + d)·n = (h + f)·m` gives a solution of the five-form system for `(m, n)` at `t = (h + f)/(b − d)` |
| `IsMSS.ratio_ne_of_le_200` | **In a magic square of squares with `b ≠ d` and `h + f ≠ 0`, `(b + d) : (h + f)` is not `m : n` for any coprime `1 ≤ n < m ≤ 200` outside the 60 exceptions** |

| `IsMSS.ratio_ne_of_le_1000` | The same for `m ≤ 1000` outside the 425 exceptions |

The last two lines are statements about the open problem itself: they exclude every magic square of squares whose edge roots have one of 303,766 ratios, whatever the size of the entries. The signs of the roots are arbitrary, so the same holds for `(b − d) : (h + f)`; in terms of the family, the system is symmetric under exchanging `(m : n)` with `(b : a)`, so the parameter `t` of a solution must avoid the excluded ratios as well. The theorems use only `propext`, `Classical.choice` and `Quot.sound`.

`scripts/g1_full/gen_lean.py` writes `G1FullCertificates.lean` from `modp_classes.txt`.

**A consequence for size (Lean: `IsMSS.edge_root_gt_500`, in `G1Size.lean`).** By the symmetry just mentioned, a solution in which `(m : n)` and `t` both have height at most 1000 needs both of them among the 425 exceptions. `G1SizeCheck.lean` shows that no two exceptional pairs fit together: for each of the 361,250 combinations one of the five forms is a non-residue modulo one of 18 small primes (3, 5, 11, 13 and others, listed in that file; the 34 sieve primes are all ±1 mod 8 and are poor at this job). So in a magic square of squares with `b² ≠ d²`, `h² ≠ f²` and `b² ≠ f²`, one of the edge roots `b, d, f, h` is larger than 500 in absolute value. Published computer searches exclude far larger squares; this is only the bound that follows from the theorems here.

## Reproduce

```
cd scripts/g1_full
gcc -O2 -o classify classify.c && ./classify 9700 > modp_classes.txt
python3 local.py 200            # writes local5_200.json   (needs numpy)
python3 verify_bruteforce.py
python3 local.py 1000           # about 4 minutes
python3 weil49.py
python3 global.py FILE.json OUT.json   # FILE has {"open": [[m, n], ...]}; needs cypari2
python3 gen_lean.py             # rewrites MagicSquareOfSquares/G1FullCertificates.lean
python3 gen_census1000.py       # rewrites MagicSquareOfSquares/G1Census1000.lean
python3 gen_sizecheck.py        # rewrites MagicSquareOfSquares/G1SizeCheck.lean
python3 families.py G2 700 200  # the sister families: no obstruction anywhere
cd ../.. && lake build          # the Lean theorems; G1FullCertificates needs about 5 GB of memory
```
