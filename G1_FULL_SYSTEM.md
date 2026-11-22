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
| An edge entry exceeds 250,000 when the nine entries are distinct | Lean (`IsMagicSquareOfSquares.edge_gt`) |
| `m² + n²` is a square modulo 7, 17, 31, 47 | Lean (`G1.isSquare_norm_mod`, `IsMSS.edge_norm`) |
| No ratio is excluded modulo a prime `≡ ±3 (mod 8)` | Lean (`G1.exists_point_mod_of_not_isSquare_two`) |
| The excluded ratios have density 1 | Argued in Result 1b. Not in Lean |
| 44 of 60, 183 of 306, 85 of 119 closed by rank 0; 157 left | Computed with PARI. Not in Lean |
| 112 of those 157 closed by the Demjanenko–Manin argument; 45 open | Computed with PARI and Sage. Not in Lean |
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

There are 484 obstructed classes in total (`modp_classes.txt`), all at ratios of good reduction. All 34 primes are `±1 mod 8`; Result 1b says why. The classes for 17 and 41 agree with the Kepler71 lists.

## Result 1b: the number `m² + n²`

**Why only primes `±1 (mod 8)` obstruct (Lean: `G1.exists_point_mod_of_not_isSquare_two`).** In integer form, at `t = 1` the five forms are `m² + n²`, `m² + n²`, `4(m² + n²)`, `4(m − n)²`, `4(m + n)²`, and at `t = 0` they are `m²`, `n²` and three times `2(m² + n²)`. So a class `m : n` can be obstructed modulo `p` only if both `m² + n²` and `2(m² + n²)` are non-residues. If `2` is a non-residue this is impossible. Hence no congruence modulo a prime `p ≡ ±3 (mod 8)` excludes a ratio, and every obstructed class has `m² + n²` a non-residue. (Prime powers are different: see the primes `≡ 5 (mod 8)` below.)

**At 7, 17, 31 and 47 the converse holds (Lean: `G1.isSquare_norm_mod`).** At these four primes every class with `m² + n²` a non-residue is obstructed: 4 of 4 classes, 8 of 8, 16 of 16, 24 of 24. So a solution has `m² + n²` a square modulo 7, 17, 31 and 47. At the other 30 primes only some non-residue classes are obstructed, for example 8 of 12 at 23, 8 of 20 at 41, 8 of 304 at 607.

**As conditions on a magic square of squares (Lean: `IsMSS.edge_norm`).** Take any magic square of squares with roots `a b c / d e f / g h i`, `b ≠ d`, `h + f ≠ 0`. Write `b + d = ±g·m`, `h + f = ±g·n` with `g = gcd(b + d, h + f)`. Then

- `m² + n²` is a square modulo 7, 17, 31 and 47;
- every prime `p ≡ 5 (mod 8)` divides `m² + n²` to an even power (this is the statement of the Kepler71 notes, proved in `G1Primes.lean`).

Since `m`, `n` are coprime, no prime `≡ 3 (mod 4)` divides `m² + n²`. So `m² + n² = 2^ε·k²·c`, where every prime factor of `c` is `≡ 1 (mod 8)`, and `c` is a square modulo 7, 31, 47, and modulo 17 unless `17 ∣ c`. The cases `c = 1` are exactly the pairs with a trivial point. The signs of the roots are arbitrary, so the same holds with `b − d` or `h − f` in place of `b + d` or `h + f`, and, rotating the square, for the pair `d + h`, `f + b`. I did not find the conditions at 7, 17, 31, 47 in the Kepler71 notes or in the sources I read; I cannot say they are new.

In Gaussian integers: `β = b + h·i` and `δ = d + f·i` both have norm `2e²`, and `β + δ = 2A·(m + n·i)`. So `m² + n²` is the norm of the primitive part of `β + δ`, and the exchange of `(m : n)` with the parameter is the exchange of the two factors in `β = (m + n·i)(A − B·i)`.

**Density (argued, not in Lean).** The theorems hold for all `m`, `n`, with no bound.

- The 484 classes alone exclude all ratios outside a set of density `∏ (1 − c_p/(p + 1)) ≈ 0.00277`, the product over the 34 primes, `c_p` the number of obstructed classes at `p`.
- The condition at primes `≡ 5 (mod 8)` alone excludes a set of density 1. For one such prime the ratios with `p` dividing `m² + n²` to an odd power have density `2p/(p + 1)²`; these conditions are independent for different primes; and `∑ 1/p` over primes `≡ 5 (mod 8)` diverges. The decay is slow: of 20,000 random coprime pairs, 23% pass this condition at height `10³` and 16% at height `10⁶`.
- For `m ≤ 1000`: the classes exclude 303,204 pairs, the primes `≡ 5 (mod 8)` exclude 237,715, together 303,766.

So the set of ratios not excluded by the Lean theorems has density zero. This says nothing about whether a magic square of squares exists: a single surviving ratio would be enough.

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

**Left by rank 0 for `m ≤ 200` (16 pairs, 10 up to the symmetry):** (15,8)~(23,7), (21,20)~(41,1), (72,65)~(137,7), (91,60)~(151,31), (99,20)~(119,79), (117,44)~(161,73), (156,133), (176,57), (187,84), (195,28). Kepler71 reports a computational proof for (15,8) (see the correction above); if it stands, nine remain.

**Left by rank 0 for `m ≤ 1000` without a trivial point (34 pairs):** (320,289), (388,287), (423,332), (427,320), (540,239), (541,476), (560,33), (593,527), (609,31), (664,301), (675,101), (724,325), (724,665), (740,539), (747,107), (755,91), (773,192), (779,301), (784,391), (788,487), (820,799), (828,173), (836,745), (840,83), (864,151), (889,320), (899,724), (920,371), (923,757), (935,416), (964,541), (965,363), (965,581), (988,27).

These 157 pairs are listed in `scripts/g1_full/open_pairs_1000.txt`; `global.py` reproduces the closing curve for each closed pair.

**In total: of the 304,191 pairs with `m ≤ 1000`, 157 are left by congruences and rank-0 quotients.** Result 4 closes 112 of them.

For the 119 pairs without a trivial point, the bad primes and 2 were also tested to prime powers up to 3·10⁶; no further obstruction appeared.

## Result 4: pairs closed by the Demjanenko–Manin argument

The 157 pairs above have positive rank on all twelve quotients. Two elliptic quotients occur twice in the Jacobian of the genus-49 curve, and a curve that occurs twice can be used when its rank is exactly 1. The Kepler71 notes use this argument for their slope family (`FABLE_JOINT_SIEVE_2026-09-15.md`, with the symmetry `z ↦ −z`); the two curves and the computation below are for G1.

**Curve `E04`.** With `S = m² + n²`, `κ = 2Sn²`: `E: y² = x(x + 2Sm²)(x + 2Sn²)`. A solution `t` gives a point `P` with `x(P) = κt²` (from `F0`, `F4`) and a point `P′` with `x(P′) = κ/t²` (from `F8`, `F4`, after `w ↦ 1/w`). So `x(P)·x(P′) = κ²`.

**Curve `D1`.** With `z = t + 1/t`: `F0·F8·F4·L` and `F0·F8·F4·U` are `t⁴` times a quartic in `z`, and the two quartics differ by `z ↦ −z`. This gives an elliptic curve with two points `P`, `P′` such that `x(P′) = −x(P)`.

**The argument.** If `E(ℚ) = ⟨G⟩ + torsion`, write `P = aG + T`, `P′ = a′G + T′`. The naive heights of `x(P)` and `x(P′)` differ by at most `2 log κ` (by 0 for `D1`), and the naive and canonical heights differ by at most Silverman's bound `B`. So `|a² − a′²|·ĥ(G) ≤ 2B + 2 log κ`. Either `|a|, |a′|` are below an explicit bound (between 19 and 32 in all cases here), and these points are listed; or `a′ = ±a`, and then `x(P)·x(P + T″) = κ²` for a torsion point `T″`, a polynomial equation for `x(P)`. Every candidate `t` is tested in the five-form system.

**Result.** 112 of the 157 pairs are closed: 91 by `E04`, 21 by `D1`, each for the pair or for its partner `(m + n, m − n)`. No candidate gave a non-degenerate square. The list, with the curve and the bound for each pair, is `scripts/g1_full/dm_closures.txt`.

**Open after this, `m ≤ 200` (5 pairs, 3 up to the symmetry):** (72,65)~(137,7), (91,60)~(151,31), (195,28).

**Open after this, `m ≤ 1000`: 45 pairs**, 14 without a trivial point and 31 with one. They are listed at the end of `dm_closures.txt`. For them both curves have rank at least 2, or the rank is not determined, or (one pair, (784,391)) no generator was found.

Kepler71 reports a computational proof for (15,8); the closure here is independent of it and uses a different method.

**Nine further elliptic quotients do not help.** `y² = F0·L`, `F8·L`, `F4·L` (in `t`), `f0·f8·K`, `f0·g4·K`, `f8·g4·K` (in `w`), and three quartics in `t ± 1/t` are also quotients. In every pair I tried each has a rational point of infinite order (for example `t = m/n` on `y² = F0·L`, `w = 1` on the three in `w`), so none of them has rank 0.

## Limits

- "Open" means not closed by these methods. It does not mean a solution is likely.
- "Not excluded by a congruence" is not a proof of local solubility: at primes of bad reduction only small prime powers were tested.
- The rank bounds depend on PARI 2.17.2. The `L`-value closures also depend on its numerical `L`-function code.
- The Weil-bound step is a standard argument, written by an AI model and not checked by a number theorist.
- Result 4 depends on: the upper rank bound 1 from PARI `ellrank`; a generator from eclib (`gens(proof=True)` in Sage) or from a PARI point saturated with `ellsaturation` up to a bound that follows from Sage's lower bound for the canonical height (Cremona–Siksek); Silverman's height bound as implemented in Sage; floating-point heights (the bound on `|a|` is rounded up by one); and my derivation of the two relations, which nobody has checked. `passagemath` 10.8 was used for Sage.
- The density statement in Result 1b is a standard argument (the Chinese remainder theorem and Dirichlet's theorem), not formalized.
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

The last two lines are statements about the open problem itself: they exclude every magic square of squares whose edge roots have one of 303,766 ratios, whatever the size of the entries. The signs of the roots are arbitrary, so the same holds for `(b − d) : (h + f)`; in terms of the family, the system is symmetric under exchanging `(m : n)` with `(v : u)`, where `t = u/v` in lowest terms (the letters of `G1Size.lean`), so the parameter `t` of a solution must avoid the excluded ratios as well. The theorems use only `propext`, `Classical.choice` and `Quot.sound`.

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
python3 dm.py - 15,8            # Result 4, curve E04, for the pair and its partner; needs passagemath-schemes, passagemath-eclib
python3 dm2.py 119,79           # Result 4, curve D1; also needs passagemath-symbolics
cd ../.. && lake build          # the Lean theorems; G1FullCertificates needs about 5 GB of memory
```
