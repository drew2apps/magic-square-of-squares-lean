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
| No ratio is excluded by a congruence modulo a prime `p ≡ ±3 (mod 8)` (modulo `p`, not modulo powers of `p`) | Lean (`G1.exists_point_mod_of_not_isSquare_two`) |
| The excluded ratios have density 1 | Argued in Result 1b. Not in Lean |
| 44 of 60, 183 of 306, 85 of 119 closed by rank 0; 157 left | Computed with PARI. Not in Lean |
| 121 of those 157 closed by the Demjanenko–Manin argument; 36 open | Computed with PARI, Sage and, for two pairs, Magma. Not in Lean, except the curve-membership algebra (`G1Relations.lean`) |
| A prime `≡ 1 (mod 8)` divides the centre root unless `m² + n²` is a square or twice a square | Lean (`IsMSS.norm_isSquare_or_prime_dvd_center`) |
| Exactly 34 primes; the Weil-bound argument | Computed and argued. Not in Lean |

## Correction to the earlier audit

The audit was written from two Kepler71 notes dated 12 September 2026. Later notes in the same repository go further, and I had not read them:

- Their later notes report all 127 pairs with `m ≤ 20` as excluded, not 120. For the last one, (15,8), a note of 13 September gives a computational proof and lists its audit as pending; their README still names (15,8) and (11,4) as open. I have not checked either claim.
- The five-form local test below is theirs (`PARALLEL_RESULT_FULL_NINE_126_2026-09-12.md`). They state the forbidden classes modulo 17 and 41 for all `(m : n)`.
- They also use other elliptic quotients of the same curves, and report 1,715 of the 12,231 pairs with `m ≤ 200` closed by a rank-0 criterion.

So the audit's phrase "seven open pairs" describes their state on 12 September only, and class-wide congruence obstructions were not new. What follows was checked against their notes by text search; I did not find it there, but I have not read all 256 files.

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

For `p ∤ 2mn(m² − n²)(m² + n²)(m⁸ + 14m⁴n⁴ + n⁸)` the five conditions define a smooth curve of genus 49 over `F_p` (a `(Z/2)⁵` cover of the line branched at ten points; the Kepler71 notes give the same genus and cite Auel–Singer, arXiv:2609.09351). The Weil bound gives at least `p + 1 − 98√p` points, which is positive for `p ≥ 9613`. `classify.c` tests every ratio `(m : n)` for every odd prime below 9700.

The octic factor was missing from this statement before version 0.3.2. Up to squares it is the resultant of `F0` and `GL` (and of `F0`, `GU`; `F8`, `GL`; `F8`, `GU`): at a prime dividing it two branch points collide, the curve modulo `p` is singular, and the bound for a smooth curve does not apply. For such a ratio and `p ∤ 2mn(m⁴ − n⁴)`, `F0·F8·GL·GU` is `(2Smn)²` times the square of a quartic modulo `p`, and no other product of some of the five forms is a constant times a square. Weil's bound for character sums then gives `32N ≥ 2p − 100 − 270√p` for the number `N` of `t` modulo `p` at which all five forms are nonzero squares, which is positive for `p ≥ 18400`; such a `t` lifts to `ℚ_p`. `octic_check.py` checks the identity, and checks directly the 2112 ratios of this kind at the 264 primes between 9700 and 20000 that have one: each has such a `t`. Below 9700 these ratios were among those tested by `classify.c`. So the list of 34 primes is unchanged.

*Derivation of the bound, for a ratio on the octic.* (For the other ratios the ten roots are distinct and the bound for the smooth curve is used.) Let `χ` be the quadratic character modulo `p` and `f₁, …, f₅` the five forms. Modulo `p` their ten roots are six points: `±t₀` (roots of `F0`), `±1/t₀` (roots of `F8`), where `t₀, 1/t₀` are the roots of `GL` and `−t₀, −1/t₀` those of `GU`, and `±i` (roots of `G4`); they are distinct because `p ∤ 2mn(m⁴ − n⁴)`. For `t` not a root, `∏(1 + χ(fᵢ(t)))` is 32 if all five values are squares and 0 otherwise; for the at most 6 roots it is at most 16. So `32N ≥ ∑_I ∑_t χ(∏_{i∈I} fᵢ(t)) − 96`, the outer sum over the 32 subsets `I` of the forms. `I = ∅` contributes `p`. `I = {F0, F8, GL, GU}` contributes at least `p − 4`: the product is a nonzero square away from its four roots, by the identity `S²·F0·F8·GL·GU = (2mn·B)² + 4·(m⁸ + 14m⁴n⁴ + n⁸)·t²·B`, `B = S²(t⁴ + 1) + (2S² − 16m²n²)t²` (Lean: `G1.octic_square`; the four resultants are `G1.resultant_F0_GL` and `G1.resultants_octic`). In each of the other 30 products some root has odd multiplicity (the roots `±i` belong to `G4` only, and each of `±t₀, ±1/t₀` to exactly two of the other four forms, in the pattern above), so the product is not a constant times a square, and Weil's bound `|∑_t χ(g(t))| ≤ (r − 1)√p` for a polynomial `g` with `r` distinct roots that is not a constant times a square gives at least `−9√p` each. In total `32N ≥ 2p − 100 − 270√p`.

*Which primes are covered by what.* `p < 9700`: every ratio, including `n ≡ 0` and the ratios on the octic, by `classify.c`. `9700 < p < 20000`: ratios off the octic by the Weil bound for the smooth curve (`p ≥ 9613`), ratios on the octic by the direct check in `octic_check.py` (output in `octic_check.txt`). `p > 20000`: the Weil bound, and for ratios on the octic the bound above (`p ≥ 18400`). `octic_check.py` works with `n = 1`: the forms are homogeneous in `(m, n)`, and a ratio on the octic has `n ≢ 0`, since at `n ≡ 0` the octic is `m⁸ ≢ 0`.

**There are exactly 34 primes with an obstructed class, the largest is 607:**

7, 17, 23, 31, 41, 47, 73, 79, 89, 103, 113, 127, 137, 151, 167, 191, 193, 199, 223, 233, 239, 257, 263, 271, 281, 311, 313, 337, 353, 439, 457, 521, 599, 607.

There are 484 obstructed classes in total (`modp_classes.txt`), all at ratios of good reduction (for the octic factor this is checked by `octic_check.py`). All 34 primes are `±1 mod 8`; Result 1b says why. The classes for 17 and 41 agree with the Kepler71 lists.

## Result 1b: the number `m² + n²`

Everything here is about 3×3 squares. "No bound" means no bound on the size of the entries.

**No obstruction modulo a prime `≡ ±3 (mod 8)` (Lean: `G1.exists_point_mod_of_not_isSquare_two`).** In integer form, at `t = 1` the five forms are `m² + n²`, `m² + n²`, `4(m² + n²)`, `4(m − n)²`, `4(m + n)²`, and at `t = 0` they are `m²`, `n²` and three times `2(m² + n²)`. So a class `m : n` can be obstructed modulo `p` only if both `m² + n²` and `2(m² + n²)` are non-residues. If `2` is a non-residue this is impossible. Three remarks on the scope:

- This is a necessary condition for a prime to obstruct. It does not say which primes `±1 (mod 8)` do; that the list of 34 is complete is computed (Result 1).
- It concerns congruences modulo `p` only. Modulo a power of `p` such a prime can obstruct: the condition at primes `≡ 5 (mod 8)` below is an example. Other prime-power obstructions are not studied here.
- The class `n ≡ 0` is never obstructed either: for `(m, n) ≡ (1, 0)` the forms at `(u, v) = (1, 1)` are `1, 4, 1, 4, 4`. So modulo `p` there are `p + 1` classes, and the obstructed ones, `c_p` in number, are among the `p` classes with `n ≢ 0`.

**At 7, 17, 31 and 47 the converse holds (Lean: `G1.isSquare_norm_mod`).** At these four primes every class with `m² + n²` a non-residue is obstructed: 4 of 4 classes, 8 of 8, 16 of 16, 24 of 24. So a solution has `m² + n²` a square modulo 7, 17, 31 and 47 (if `p` divides `m² + n²` this says nothing). This is a reading of the class table, checked by the kernel; I know no reason why it holds at these four primes and fails at the others. The numbers `c_p` are:

7:4, 17:8, 23:8, 31:16, 41:8, 47:24, 73:8, 79:24, 89:28, 103:20, 113:8, 127:8, 137:8, 151:36, 167:24, 191:24, 193:8, 199:8, 223:16, 233:12, 239:8, 257:16, 263:16, 271:24, 281:16, 311:16, 313:16, 337:16, 353:8, 439:16, 457:8, 521:8, 599:8, 607:8 (total 484).

**As conditions on a 3×3 magic square of squares (Lean: `IsMSS.edge_norm`; for the matrix form `IsMagicSquareOfSquares.edge_norm`).** Take a magic square of squares with roots `a b c / d e f / g h i` such that `b ≠ d` and `h + f ≠ 0`. This holds when the nine entries are distinct and the roots are taken positive. Write `|b + d| = g·m`, `|h + f| = g·n` with `g = gcd(b + d, h + f)`. Then

- `m² + n²` is a square modulo 7, 17, 31 and 47;
- every prime `p ≡ 5 (mod 8)` divides `m² + n²` to an even power (this is the statement of the Kepler71 notes, proved in `G1Primes.lean`).

Since `m`, `n` are coprime, `m² + n²` is odd or twice an odd number, and no prime `≡ 3 (mod 4)` divides it. So `m² + n² = 2^ε·k²·c` with `ε ∈ {0, 1}`, `c` squarefree, every prime factor of `c` `≡ 1 (mod 8)`, and `c` a square modulo 7, 31, 47, and modulo 17 unless 17 divides `m² + n²`. `c = 1` means that `m² + n²` is a square or twice a square; these are exactly the pairs with a trivial point. Replacing a root by its negative, or rotating the square, gives another magic square of squares, so the same Lean theorem applies to the ratios `(b − d) : (h + f)`, `(b + d) : (h − f)` and `(d + h) : (f + b)`; these variants are not stated separately.

**The centre (Lean: `IsMSS.norm_dvd_center`, `IsMSS.norm_isSquare_or_prime_dvd_center`).** With `m`, `n` as above there are integers `u`, `v` with `b + d = m·u`, `h + f = n·u`, `b − d = n·v`, `f − h = m·v`, and `8e² = (m² + n²)(u² + v²)`. In Gaussian integers: `2(b + h·i) = (m + n·i)(u − v·i)` and `2(d + f·i) = (m + n·i)(u + v·i)`. Hence every odd prime dividing `m² + n²` divides the centre root `e`, and:

> either `m² + n²` is a square or twice a square, or some prime `p ≡ 1 (mod 8)` divides the centre root and divides `m² + n²` to an odd power.

So a magic square of squares whose centre root has no prime factor `≡ 1 (mod 8)` has both ratios of trivial type (the second ratio is `u : v`, obtained by negating `h`). The theorems need no primitivity and no distinctness beyond `b ≠ d`, `h + f ≠ 0`; these force `e ≠ 0`. The matrix form is `IsMagicSquareOfSquares.centre`.

**Rabern's paper.** I have now read Rabern (2003). Its theorems are about primes dividing entries; nothing there concerns `m² + n²` or the ratio of edge sums. It ends with a question: whether in a solution all prime divisors of all entries are `≡ 1 (mod 8)`, or whether, on the contrary, some entry must have a prime divisor `≡ 5` or `7 (mod 8)`. The centre statement above is a small step of a related kind: unless both ratios are of trivial type, the centre root does have a prime factor `≡ 1 (mod 8)`.

**Other sources.** C. Woll, *A Partial Residue Categorization of the Magic Square of Squares* (arXiv:1809.03067), studies the square modulo the primes dividing the centre root; its Corollary 4.3 says that a prime dividing both the centre and a mid-edge entry of a primitive square is `≡ 1 (mod 8)`. I went through its list of statements with an AI summary, not line by line, and did not find `m² + n² ∣ 8e²` or the statement on the centre proved here. O. Hill, arXiv:2510.08286 (October 2025), claims a proof that no 3×3 magic square of distinct squares exists; I have not verified it and do not know its status, and nothing here depends on it.

**What I searched.** The Kepler71 repository has 256 notes. I searched them by keyword (non-residue, quadratic, `m² + n²` together with residue or Legendre, `±1 (mod 8)`) and read about ten in part. I found the statement on primes `≡ 5 (mod 8)` and none of the others. I have not read Zimmermann's, Morgenstern's or Robertson's papers. No novelty is claimed for anything in this section.

**Density (argued; not in Lean).** Order the coprime pairs `1 ≤ n < m` by `m`. The density of a set of pairs is the limit, if it exists, of its proportion among the pairs with `m ≤ X`. A condition on `m : n` modulo a fixed integer `N` has a density, equal to the proportion of the points of the projective line over `ℤ/N` that satisfy it; for coprime moduli these proportions multiply (Chinese remainder theorem). Sketch of a proof of the first statement (written with an AI model, not checked by anyone else; below, `N` is a product of primes and squares of primes, since the conditions used are conditions on the ratio modulo `p` or `p²`): fix `(a, b)` modulo `N` with `gcd(a, b, N) = 1`. The pairs `(m, n) ≡ (a, b)` in the triangle `1 ≤ n < m ≤ X` number `X²/(2N²) + O(X)`. A common divisor `d` of such a pair is prime to `N`, and the pairs divisible by `d` number `X²/(2d²N²) + O(X/d)`, so by Möbius inversion the coprime ones number `(X²/(2N²))·∏_{p∤N}(1 − p⁻²) + O(X log X)`. This does not depend on `(a, b)`. A point of the projective line over `ℤ/N` has `φ(N)` representatives `(a, b)`, all with `gcd(a, b, N) = 1`, so every point receives the same share of the coprime pairs; there are `N·∏_{p∣N}(1 + 1/p)` points.

- The ratios in none of the 484 classes have density `∏ (1 − c_p/(p + 1)) = 0.0027661`.
- Let `P` be a finite set of primes `≡ 5 (mod 8)`. The ratios for which no `p ∈ P` divides `m² + n²` exactly once form a set defined modulo `∏ p²`, of density `∏ (1 − 2(p − 1)/(p(p + 1)))`. Every ratio that passes the condition at primes `≡ 5 (mod 8)` lies in this set, for every `P`. As `P` grows the product tends to 0, because `∑ 1/p` over primes `≡ 5 (mod 8)` diverges (the primes `≡ a (mod m)`, `gcd(a, m) = 1`, have Dirichlet density `1/φ(m)`: Serre, *A Course in Arithmetic*, chapter VI, §4, Theorem 2; a set `A` of primes of positive Dirichlet density has a divergent sum of reciprocals, because `∑_{p∈A} p⁻ˢ → ∞` as `s → 1⁺` and `∑_{p∈A} 1/p ≥ ∑_{p∈A} p⁻ˢ` for every `s > 1`). So the ratios that pass have density 0. Only finite sets of primes are used, so no sieve with error terms is needed.
- The decay is slow. `scripts/g1_full/density.py` draws 20,000 pairs uniformly from the coprime pairs `1 ≤ n < m ≤ H` (seed 1) and prints the proportion that the second condition does not exclude for `H = 10³, 10⁴, 10⁵, 10⁶`: 22.05%, 18.70%, 17.36%, 16.08%.
- For `m ≤ 1000` (same script): 303,204 pairs lie in one of the classes, 237,715 fail the condition at primes `≡ 5 (mod 8)`, 237,153 do both, 303,766 do at least one, 425 do neither. So 562 pairs are excluded by the second condition only. The Lean census uses the second condition only in the form "`p` divides `m² + n²` exactly once", for the 19 primes of `fivePrimes1000` (5, 13, 29, 37, 53, 61, 109, 149, 157, 173, 181, 269, 373, 541, 557, 701, 773, 3461, 12277); the generator picks, for each pair not in a class, the first prime that works, and these 19 are the ones it needed. The 425 exceptions are the same.

This says nothing about whether a magic square of squares exists: a single surviving ratio would be enough.

## Result 2: census

Local tests used: the 34 primes above, the 2-adic test modulo 32, 64, 128, and powers of the odd primes dividing `mn(m² − n²)(m² + n²)`.

**Second implementation.** `local_all.py` is a second program for solubility in `ℚ_p`, written separately from the census (the five forms are typed in again; nothing else is shared). It does not stop at a fixed power of `p`. Its answer is soluble, insoluble or undecided: for `p > 5000` the search among the residues that are not roots is by random sampling, which can only answer soluble or undecided, so an answer "insoluble" always comes from the exhaustive and the recursive branches. It was run on all 304,191 pairs (`local_all_census.py`, output in `local_all_census.txt`) with the bound 700 passed explicitly, that is with the primes below 700 and every prime dividing `2mn(m² − n²)(m² + n²)(m⁸ + 14m⁴n⁴ + n⁸)`, it finds a prime of insolubility for exactly the 303,766 pairs excluded by the two conditions of Result 1b, and none for the other 425. No pair was undecided at any prime. With Result 1 for the remaining primes, the 425 pairs are soluble in every `ℚ_p`: no congruence closes any of them.

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

**In total: of the 304,191 pairs with `m ≤ 1000`, 157 are left by congruences and rank-0 quotients.** Result 4 closes 121 of them.

For the 119 pairs without a trivial point, the bad primes and 2 were also tested to prime powers up to 3·10⁶; no further obstruction appeared.

## Result 4: pairs closed by the Demjanenko–Manin argument

The 157 pairs above have positive rank on all twelve quotients. Some elliptic quotients occur more than once in the Jacobian of the genus-49 curve, and such a curve can be used when its rank is exactly 1. The Kepler71 notes use this argument for their slope family (`FABLE_JOINT_SIEVE_2026-09-15.md`, with the symmetry `z ↦ −z`); the curves and the computation below are for G1.

| Curve | Comes from | Times in the Jacobian | Use |
|---|---|---|---|
| `w·f0·f8` (`E_B`), `w·K`, and the other curves of Result 3 | `F0·F8`, `L·U`, … | 1 | rank 0 |
| `E04` | `F0·F4` and `F8·F4` | 2 | this section |
| `F4·L`, `F4·U` | | 2 | it is `E04` of the partner pair |
| `D1` | `F0·F8·F4·L` and `F0·F8·F4·U` | 2 | this section |
| `f0·g4·K`, `f8·g4·K` | | 2 | it is `D1` of the partner pair |
| `E0L` | `F0·L`, `F0·U`, `F8·L`, `F8·U` | 4 | this section, with two of the four points |
| `f0·f8·K` and two quartics in `t ± 1/t` | `F0·F8·L·U` | 1 each | none: each has the rational point `w = 1` or `t = ±1`, and PARI never gave rank 0 in 15 pairs tried |

**Curve `E04`.** With `S = m² + n²`, `κ = 2Sn²`: `E: y² = x(x + 2Sm²)(x + 2Sn²)`. A solution `t` gives a point `P` with `x(P) = κt²` and a point `P′` with `x(P′) = κ/t²`. So `x(P)·x(P′) = κ²`.

**Curve `D1`.** With `z = t + 1/t` and `c0 = −16Smn(m² − n²)²`, a solution gives two points on one cubic with `x(P) = c0/z` and `x(P′) = −c0/z`.

Lean checks the algebra of curve membership only (`G1Relations.lean`): the abscissas above give a square times the stated product of forms, so the points lie on the stated curves, on the same curve and not on a twist (`e04_point`, `e04_point'`, `d1_point`, `d1_point'`, with the same cubic as in `dm2.py`); `F0·GL` at `−t` is `F0·GU` at `t`, and `t = m/n` is a point (`e0l_neg`, `e0l_base`); and the forms of the partner pair are the forms of the pair (`partner_forms`). Ranks, generators and heights are not in Lean.

**Curve `E0L`.** `C: y² = F0(t)·GL(t)` has the rational point `t = m/n`. A solution `t` gives rational points on `C` with `t`-coordinates `t` and `−t`. Here `t` is not an even function on the curve, so its height is compared with canonical heights through the fibre `{Q, R0 − Q}` of `t`: the sum and the product of the two Weierstrass abscissas are quadratic polynomials in `σ = 1/(t − m/n)`, which gives `|2h(σ(Q)) − ĥ(Q) − ĥ(R0 − Q)| ≤ K` with an explicit `K` (`dm3.py`, which also checks this inequality numerically on every point it lists and checks that every listed point is on `C`).

**The argument.** If `E(ℚ) = ⟨G⟩ + torsion`, write `P = aG + T`, `P′ = a′G + T′`. The naive heights of the two abscissas differ by a bounded amount, and naive and canonical heights differ by at most Silverman's bound `B`. For `E04`: `|a² − a′²|·ĥ(G) ≤ 2B + 2 log κ`. Either `|a|, |a′|` are below an explicit bound, and these points are listed; or `a′ = ±a`, and then `P′ = ±P + T″` with `T″` torsion, a polynomial equation for `x(P)`, solved for every torsion point. Every candidate `t` is tested in the five-form system, so a superset of candidates is harmless. The partner of `(m, n)` is `(m + n, m − n)`, divided by 2 when both are even.

**Result.** 121 of the 157 pairs are closed: 91 by `E04`, 24 by `D1`, 6 by `E0L`, each for the pair or for its partner. No candidate gave a non-degenerate square. The bound on `|a|` used is 27 to 34 for `E04`, 21 or 22 for `D1`, 11 or 12 for `E0L`. `scripts/g1_full/dm_evidence.txt` has one line per closed pair: generator, heights, bounds, torsion, number of candidates, how the generator was proved, and the software versions.

**The bound for `E0L` in symbols.** With `τ = t − m/n`, `σ = 1/τ`, `q0 = 2m(m² − n²)/n` and `q(m/n + τ) = aτ⁴ + bτ³ + cτ² + dτ + q0²`, Connell's formulas give a Weierstrass model in which `X(Q) + X(R0 − Q) = 4q0²σ² + 2dσ` and `X(Q)·X(R0 − Q) = (d² − 4q0²c)σ² − 4q0²bσ − 4q0²a`, where `R0` is the second point with `t = m/n`. `C1` is the logarithm of the larger of the row-sum norms of the 3×3 integer matrix of these quadratics (after the change to the minimal model) and of its adjugate; `K = C1 + log 2 + 2B`; `ρ` is the integer with `R0 = ρG + torsion`; and with `W = (2K + 4(log max(m, n) + log 2))/(2ĥ(G))` the proved bound is `|a| ≤ W + |ρ|/2`. The script lists `|a| ≤ ⌊1.02·(W + |ρ|)⌋ + 3`, which is larger; for `E04` and `D1` it lists the proved bound times 1.02, rounded down, plus 2.

**A worked example for `E0L`, (423,332).** Minimal model `y² = x³ − x² − 10810266634575227101441x + 180302397457159147910876363714305`, rank bounds `[1, 1]`, torsion of order 4, generator with `x(G) = −27869819136`, `ĥ(G) = 25.54`, `B = 42.53`, `C1 = 94.27`, `K = 180.02`, `ρ = 1`, `W = 7.58`. The script lists `|a| ≤ 11` with all four torsion translates, 92 points; the largest observed `|2h(σ) − ĥ(Q) − ĥ(R0 − Q)|` is 35.4, below `K`. No listed `t`, and no `t` from the torsion equations, makes the five forms squares.

**A worked example, (15,8).** `S = 289`, `κ = 36992`, `E: y² = x(x + 130050)(x + 36992)`. PARI gives rank bounds `[1, 1]`; eclib proves that `G` with `x(G) = −83521 = −S²` generates `E(ℚ)` modulo its torsion subgroup of order 4. `ĥ(G) = 3.1576`, `B = 29.242`, so `|a| ≤ (2B + 2 log κ)/ĥ(G) = 25.2`; the script lists `a ≤ 27`. Among the 111 candidate abscissas, the only ones of the form `κt²` at which the five forms are squares are `t = ±1`, the degenerate square of the pair.

**The line for (784,391).** PARI's `ellrank` finds no point on this `D1` curve at efforts 0 to 3 and finds `x = 418733994302308977769` (minimal model) at effort 4. Its canonical height is 13.3913, the same as that of the point `x = 0` of the `D1` model, which is the generator in the other 23 `D1` lines as well. The saturation bound 122 in that line is large only because the proved lower bound for the height of a non-torsion point on this curve is small (0.0009): `⌊√(13.39/0.0009)⌋ + 1 = 122`. PARI's `ellsaturation` up to 122 leaves the point unchanged. This is the line I would check first. Magma's `TwoPowerIsogenyDescentRankBound` also gives rank at most 1 for this curve.

**The lines for (631,71) and (351,280).** For the `D1` curve of (631,71) PARI gives the rank bounds `[1, 3]`. Magma's `TwoPowerIsogenyDescentRankBound` (public Magma calculator, version 2.29-10) gives rank at most 1. With that bound, `dm2.py --upper=1 631,71` closes the pair and its partner (351,280): the generator is again `x = 0`, and the only candidate is the degenerate `t = 0`. These two lines rest on Magma alone for the upper bound and are marked `upper=magma`.

**Open after this: 36 pairs.** None can be closed by a congruence: 29 have a rational point, and the seven others are soluble in every `ℚ_p` (`local_all.py`, with all primes below 20000). For `m ≤ 200`: (72,65)~(137,7), (91,60)~(151,31), (195,28). By what blocks them:

| Pairs | Situation | What would be needed |
|---|---|---|
| 7 without a trivial point: (724,665), (740,539), (820,799), (828,173), (840,83), (889,320), (923,757) | `E0L` has rank 2; `E04` has rank 2 or 3; `D1` has rank at least 2, or bounds `[1, 3]`, `[0, 2]`, `[0, 4]` that PARI does not improve with more effort | the argument with all four points on `E0L`, which needs explicit bounds on how far the four maps are from orthogonal; or a second descent on the `D1` curves with bounds `[1, 3]` |
| (72,65)~(137,7), (425,168)~(593,257), (920,129) | `E0L` has rank 3, the others at least 2 | the same four-point argument |
| the other 24, all with a trivial point | `E0L` has rank 4 to 6 (for two pairs only bounds `[3, 5]`) | curves of genus 2 (Chabauty), or another fibration |

Kepler71 reports a computational proof for (15,8). The closure here uses a different method; both computations await outside checking.

## Limits

- "Open" means not closed by these methods. It does not mean a solution is likely.
- That the 425 pairs are soluble in every `ℚ_p` rests on `local_all.py` (primes below 700 and primes of bad reduction) and on Result 1 (the other primes). Both were written with an AI model and checked by nobody else.
- The rank bounds of Result 3 were computed with PARI 2.17.2; Result 4 was computed with PARI 2.17.4 as shipped with passagemath 10.8.14. The `L`-value closures also depend on its numerical `L`-function code.
- The Weil-bound step and the character-sum step for the singular fibres are standard arguments, written by an AI model and not checked by a number theorist. The singular fibres were overlooked until version 0.3.2, which shows the kind of error that may remain.
- Result 4 is not in Lean, except the identities of `G1Relations.lean` (curve membership for `E04` and `D1`; for `E0L` the base point, the symmetry of the quartic and the map from the quartic to a Weierstrass model, `G1.quartic_to_weierstrass` and `G1.e0l_weierstrass`, but not the change to the minimal model that Sage makes; the partner forms, their halving and the coprimality of the partner pair). What an outside reader would have to check: (1) the upper rank bound 1 from PARI `ellrank` (2-descent with the Cassels pairing). A second implementation, eclib's 2-descent without that pairing, gives the same bound for 82 of the 91 `E04` lines, on the curve or on an isogenous curve (`eclib_check.py`). For the other 9 lines, on 7 curves, (320,289), (609,31), (660,259), (700,429), (840,559), (864,151), (919,401), (935,416), (988,27), eclib gives `[1, 3]`, and for the `D1` and `E0L` lines it gives no usable bound, so there eclib does not confirm the bound (the curve of (320,289) serves the lines (320,289) and (609,31), the curve of (660,259) the lines (660,259) and (919,401); the other five curves serve one line each). eclib bounds the rank from above only; the lower bound is the recorded point of positive height. A third program, Magma's `TwoPowerIsogenyDescentRankBound` through the public Magma calculator (`magma_check.py`, answers in `magma_check.txt`), gives rank at most 1 for all 121 lines (87 distinct curves), among them (784,391), the 9 lines above and every `D1` and `E0L` line. So each upper rank bound used in Result 4 now comes from two independent programs, PARI and Magma, except for (631,71) and (351,280), where it comes from Magma alone. This routine bounds the rank from above by repeated descents along the 2-isogeny given by a rational point of order 2 and along its dual; every curve here has such a point, and the rank does not change under isogeny or under the change to the minimal model, to which the routine is applied. The requests and the verbatim responses are in `magma_raw.txt` from 22 November 2026, 06:30 PT, on; the first 80 answers were kept only as the parsed bounds; (2) that `G` is a generator: eclib (`gens(proof=True)`, or its saturation with its own index bound), or PARI `ellsaturation` with a bound `N` such that the index of the starting point is at most `N − 1` (`N − 1 = ⌊√(ĥ(P)/λ)⌋`, `λ` the height lower bound), which follows from Sage's lower bound for the canonical height (`height_function().min(0.01, 20)`; the two arguments are a tolerance and a cap on the number of multiples used, and the returned value is documented as a proved lower bound); (3) Silverman's bound as implemented in Sage, applied to the integral model whose points are listed; (4) the height comparison for `E0L` in `dm3.py`; (5) floating-point heights: the proved bound on `|a|` is multiplied by 1.02 and increased by 2 before listing (see the header of `dm_evidence.txt` for `E0L`), which is a safety margin and not interval arithmetic; the numerical check of the `E0L` inequality on the listed points is a sanity check, not a proof of it; (6) the torsion case, solved for every torsion point; (7) the final test of each candidate (`classify_t`, exact rational arithmetic; unit tests in `check_evidence.py`); (8) the passage to the partner pair (in Lean: `G1.partner_forms`, `G1.partner_forms_half`, `G1.partner_coprime`, `G1.partner_halved_coprime`). Nobody else has checked any of these; only (7) and (8) have the automatic checks named.
- The density statement in Result 1b uses only the counting statement proved there, the Chinese remainder theorem for finitely many primes, and the divergence of `∑ 1/p` over the primes `≡ 5 (mod 8)` (Dirichlet); it is not formalized.
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
python3 density.py              # the counts and the sample of Result 1b; no dependencies
python3 dm.py - 15,8            # Result 4, curve E04, for the pair and its partner; needs passagemath-schemes, passagemath-eclib
python3 dm2.py 119,79           # Result 4, curve D1; also needs passagemath-symbolics
python3 dm3.py 423,332          # Result 4, curve E0L
python3 check_evidence.py       # tests of classify_t and consistency of the tables; no dependencies
python3 octic_check.py          # Result 1, the singular fibres; no dependencies, 5 seconds
python3 eclib_check.py          # Result 4: eclib's rank bounds for the E04 lines, about 80 seconds
python3 magma_check.py          # Result 4: Magma's rank bounds through the public calculator; slow, rate-limited
python3 dm2.py --upper=1 631,71 # Result 4: the two lines that use Magma's bound
python3 local_all_census.py 0 4 # Result 2: the second local program on a quarter of the pairs (needs cypari2)
python3 local_all.py 724,665    # solubility in every Q_p for one pair (needs cypari2)
cd ../.. && lake build          # the Lean theorems; G1FullCertificates needs about 5 GB of memory
```
