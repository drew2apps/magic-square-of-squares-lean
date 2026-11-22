# Magic square of squares: known facts, checked in Lean

[![CI](https://github.com/drew2apps/magic-square-of-squares-lean/actions/workflows/ci.yml/badge.svg)](https://github.com/drew2apps/magic-square-of-squares-lean/actions/workflows/ci.yml) [![DOI](https://zenodo.org/badge/1466891268.svg)](https://doi.org/10.5281/zenodo.24435351)

Archived at Zenodo: [10.5281/zenodo.24435351](https://doi.org/10.5281/zenodo.24435351) resolves to the latest archived version, and each release has its own DOI listed there. A short paper-style summary is in [PAPER.md](PAPER.md).

Is there a 3×3 magic square whose nine entries are distinct perfect squares? The question is open. It appears in [formal-conjectures](https://github.com/google-deepmind/formal-conjectures) as `MagicSquares.exists_magic_square_squares`.

This repository proves, in Lean 4 with Mathlib, the elementary facts that are known about such a square. It does not solve the problem.

Every theorem here is checked by the Lean kernel and depends only on the standard axioms `propext`, `Classical.choice` and `Quot.sound`. There is no `sorry` and no `native_decide`.

## What is proved

The layout is

```
a² b² c²
d² e² f²
g² h² i²
```

with centre `e`, corners `a c g i` and edges `b d f h`. A square is *primitive* if no prime divides all nine entries.

### The nine necessary conditions listed on Wikipedia

The [Wikipedia article](https://en.wikipedia.org/wiki/Magic_square_of_squares) (read on 21 November 2026) lists nine properties of a primitive magic square of squares with sum `S`. All nine are proved here. Each has two versions: one for nine integer roots (`IsMSS`), and one for a function `m : Fin 3 → Fin 3 → ℕ` as in formal-conjectures (`IsMagicSquareOfSquares`).

| # | Property | Source cited by Wikipedia | Theorem |
|---|---|---|---|
| 1 | All entries are odd | Rabern | `IsMSS.odd`, `IsMagicSquareOfSquares.odd` |
| 2 | All entries are ≡ 1 (mod 24) | Zimmermann | `IsMSS.sq_emod_24`, `IsMagicSquareOfSquares.emod_24` |
| 3 | `S` is three times the centre | Gardner | `IsMagic.sum_eq`, `IsMagicSquareOfSquares.sum_eq` |
| 4 | `S` ≡ 3 (mod 72) | Zimmermann | `IsMSS.sum_emod_72`, `IsMagicSquareOfSquares.sum_emod_72` |
| 5 | No entry has a prime factor ≡ 3 (mod 8) | Rabern | `IsMSS.not_dvd_of_three_mod_eight` |
| 6 | Every prime factor of the centre is ≡ 1 (mod 4) | Rabern | `IsMSS.center_prime_factors` |
| 7 | No prime ≡ 5 (mod 8) divides an edge | Rabern | `IsMSS.not_dvd_edge` |
| 8 | A prime ≡ 3 (mod 4) dividing a corner divides the two edges that do not touch it | Rabern | `IsMSS.dvd_edges_of_dvd_corner` |
| 9 | A prime ≡ 5 (mod 8) dividing a corner divides the centre and the opposite corner | Morgenstern | `IsMSS.dvd_center_of_dvd_corner` |

Properties 3, 8 and 9 do not need the square to be primitive. Theorems 7 and 9 are proved for primes ≡ 3 (mod 8) as well; for property 9 this is already in Rabern's paper.

Checked against Rabern (2003), read on 21 November 2026: property 1 is his Theorem 1.1, property 6 is Theorem 1.2, property 5 is Corollary 1.1, property 7 is Theorem 1.4, property 8 is Theorem 1.5, each with the same conclusion, and with the same hypothesis or a weaker one (Rabern assumes entries without a common factor; the Lean proof of property 8 does not use it). Property 9, which Wikipedia cites to Morgenstern, is contained in Rabern's Theorem 1.3: a prime ≡ 3 or 5 (mod 8) dividing a non-centre entry divides the centre and the other entry in that line.

### Other results

| Result | Theorem |
|---|---|
| Every 3×3 magic square is `E+u, E-u-v, E+v / E-u+v, E, E+u-v / E-v, E+u+v, E-u` | `IsMagic.iff_param` |
| If one entry is 0 then all nine are 0 | `IsMSS.eq_zero_of_exists_eq_zero` |
| Distinct entries are automatically positive, so the positivity hypothesis in the formal-conjectures statement is redundant | `IsMagicSquareOfSquares.pos`, `openStatement_iff` |
| If a solution exists then a primitive solution exists | `IsMSS.exists_primitive`, `openStatement_iff_primitive` |
| Robertson's first reformulation: a solution exists if and only if there are three 3-term arithmetic progressions of squares with the same common difference whose middle terms form an arithmetic progression | `openStatement_iff_three_progressions` (integer squares), `openStatement_iff_rat_progressions` (rational squares) |
| Robertson's second reformulation: a solution exists if and only if there are three rational right triangles with the same area whose squared hypotenuses form an arithmetic progression | `openStatement_iff_three_triangles` |
| Robertson's third reformulation: a solution exists if and only if, for some rational `n`, the curve `y² = x³ - n²x` has three rational points, each the double of a rational point, whose `x`-coordinates form an arithmetic progression | `openStatement_iff_three_doubles` |
| For rational `n ≠ 0`: `x - n`, `x`, `x + n` are all rational squares if and only if `x` is the `x`-coordinate of the double of a rational point of `y² = x³ - n²x` | `isSquare_iff_isDoubleX` |
| The formal-conjectures statement is equivalent to the nine-integer-roots form used here | `openStatement_iff_exists_isMSS` |

In the three reformulations the nine squares are required to be distinct. The informal statements leave this implicit. The curve uses Mathlib's `WeierstrassCurve.Affine.Point` and its group law. In the curve version `n` is any rational number: the statement does not say that `n` is a congruent number.

### Checks against published examples

| Example | Theorem |
|---|---|
| Bremner's square is magic and has exactly two non-square entries | `bremner_isMagic`, `bremner_not_isSquare` |
| Sallows' square: rows, columns and one diagonal sum to 147²; the other diagonal does not | `sallows` |
| The Parker square: rows, columns and one diagonal sum to 3051; the other diagonal sums to 4107 | `parker` |
| Euler's 4×4 magic square of distinct squares (1770) | `exists_magic_square_squares_four` |

### Conditions on the edge roots

All of this is about 3×3 squares. Write the roots of a magic square of squares as `a b c / d e f / g h i`, assume `b ≠ d` and `h + f ≠ 0` (true when the nine entries are distinct and the roots are taken positive), and let `m : n` be the ratio `|b + d| : |h + f|` in lowest terms. [G1_FULL_SYSTEM.md](G1_FULL_SYSTEM.md) has the details; the forms, the five-form local test and the statement on primes `≡ 5 (mod 8)` come from the [Kepler71 project](https://github.com/Kepler71/magic-square-of-squares).

| Statement | Lean name |
|---|---|
| `m² + n²` is a square modulo 7, 17, 31 and 47, and every prime `≡ 5 (mod 8)` divides it to an even power | `IsMSS.edge_norm`, `IsMagicSquareOfSquares.edge_norm` |
| `m² + n²` divides `8e²`. Either `m² + n²` is a square or twice a square, or a prime `≡ 1 (mod 8)` divides the centre root `e` | `IsMSS.norm_dvd_center`, `IsMSS.norm_isSquare_or_prime_dvd_center`, `IsMagicSquareOfSquares.centre` |
| `m : n` lies in none of 484 residue classes at 34 primes | `IsMSS.ratio_not_in_class`, `G1.fullTable_valid` |
| In particular `m : n` is not one of the coprime pairs `1 ≤ n < m ≤ 1000` outside 425 listed pairs (60 of them up to 200) | `IsMSS.ratio_ne_of_le_1000` |
| The same for the five forms of the G1 family: no rational `t` makes them squares | `G1.no_solution_le_1000` |
| If `b² ≠ d²`, `h² ≠ f²` and `b² ≠ f²`, one of the edge roots `b, d, f, h` is larger than 500 in absolute value; so one of the four edge entries is larger than 250,000, in particular when the nine entries are distinct | `IsMSS.edge_root_gt_500`, `IsMagicSquareOfSquares.edge_gt_of_ne`, `edge_gt` |
| For a prime `p ≡ ±3 (mod 8)` no ratio is excluded by a congruence modulo `p` (powers of `p` are not covered) | `G1.exists_point_mod_of_not_isSquare_two` |

The first three rows have no bound on the size of the entries. The signs of the roots are arbitrary, so the theorems also apply with `b − d` or `h − f`. The ratios excluded by rows 1 and 3 have density 1 among all ratios (argued in G1_FULL_SYSTEM.md, not in Lean); this does not decide the problem. The bound 500 is far below what computer searches have covered.

Computed, not in Lean: exactly 34 primes can obstruct the system modulo `p`, the largest is 607. Of the 425 listed pairs, 268 are closed by elliptic curves of rank 0, 121 more by the Demjanenko–Manin argument on three curves of rank 1, and 36 are open (5 of them up to 200). These computations depend on PARI and Sage (with eclib), on Magma (it confirms every upper rank bound, and two pairs rest on it alone), and on derivations that nobody has checked; `scripts/g1_full/dm_evidence.txt` records each one.

Also computed, not in Lean (`scripts/search/seven.c`, summary in `scripts/search/seven_search.txt`): among the 3×3 magic squares with nine distinct positive entries, a square centre of root at most 10⁹ and at least seven square entries, the only ones are Bremner's square (1999) and its multiples; none has eight or nine squares. The reduction: a line through the centre with both ends square is a Pythagorean triangle with hypotenuse the centre root; seven squares force two such lines; the two lines can be any two of the four, which gives six cases for the other two deviations, all tested for every pair of triangles (the first version of the program missed two of the six; that was found in review and the whole range was run again). A second program without a table of cases, `seven_brute.c`, agrees up to centre root 1300. Scope: squares with seven square entries and a centre that is not a square are not covered. For nine squares this is much weaker than what is already known: Lee Morgenstern showed in 2006 that all entries must exceed 10¹⁴ and searched triples of arithmetic progressions of squares with common difference up to 2.4·10²¹ (reported on multimagie.com). For seven squares I know of searches by Rabern, Boyer, MacLeod and Morgenstern but not of their exact ranges; no record is claimed.

An earlier audit of the same family, restricted to three of the five forms, is in [archive/](archive/KEPLER71_G1_AUDIT.md). Its numbers are superseded.

## What is not proved

- Whether a 3×3 magic square of distinct squares exists.
- Any lower bound on the size of the entries beyond "an edge root exceeds 500" (an edge entry exceeds 250,000), which is far below what computer searches have covered.
- Anything about the 425 excepted ratios up to 1000, or about the ratios above 1000 that lie in none of the 484 classes and pass the condition at primes `≡ 5 (mod 8)`. For 306 of the 425 the system has a degenerate solution and no congruence can help. This sieve cannot settle the problem by itself, and it does not transfer to the other natural ways of cutting the problem into curves (see [G1_FULL_SYSTEM.md](G1_FULL_SYSTEM.md)).
- That the list of 34 obstructing primes is complete, and anything about obstructions modulo prime powers other than the statement on primes `≡ 5 (mod 8)`. The first is computed, the second is not studied.
- The Lean proofs use congruences and algebraic identities only. The rank, generator and height computations behind [G1_FULL_SYSTEM.md](G1_FULL_SYSTEM.md) (268 pairs by rank 0, 121 by the Demjanenko–Manin argument) are not formalized, and the derivations behind the 121 have not been checked by anyone. The density statement is argued in prose.
- Any search bound, and any result that needs the arithmetic of curves (Bremner 1999 and 2001) or the geometry of the surface of magic squares of squares (Bruin–Thomas–Várilly-Alvarado, Auel–Singer).
- The existence of magic squares of squares for every size from 4×4 up (Rome–Yamagishi). Only Euler's 4×4 example is checked.

## A remark on the formal-conjectures docstring

The docstring of `MagicSquares.exists_magic_square_squares` says that 0 is excluded because a magic square of squares with 0 and eight distinct squares is known. `IsMSS.eq_zero_of_exists_eq_zero` shows that this cannot happen for a 3×3 magic square of integer squares: an entry 0 forces `y² = 2e²` for the opposite entry `y²` and the centre `e²`, so all nine entries are 0. The formal statement itself is correct. Only the reason given in the docstring is not. We could not open the source cited there to see which example it means.

## Files

| File | Contents |
|---|---|
| `MagicSquareOfSquares/Defs.lean` | `IsMagic`, `IsMSS`, `Primitive`; facts true for every 3×3 magic square |
| `MagicSquareOfSquares/Basic.lean` | Identities, zero entries, parity |
| `MagicSquareOfSquares/Primes.lean` | Prime divisors of the entries |
| `MagicSquareOfSquares/Residues.lean` | Residues modulo 24 and 72 |
| `MagicSquareOfSquares/Primitive.lean` | Reduction to primitive squares |
| `MagicSquareOfSquares/Progressions.lean` | Robertson's arithmetic-progression reformulation |
| `MagicSquareOfSquares/Examples.lean` | Bremner, Sallows, Parker, Euler |
| `MagicSquareOfSquares/Statement.lean` | The formal-conjectures statement and the results restated for it |
| `MagicSquareOfSquares/Curve.lean` | Squares in progression and doubles on `y² = x³ - n²x` |
| `MagicSquareOfSquares/Elliptic.lean` | Rational squares; the elliptic-curve reformulation |
| `MagicSquareOfSquares/Triangles.lean` | The right-triangle reformulation |
| `MagicSquareOfSquares/G1.lean` | A finite check modulo `M` that rules out squares on the diagonal of the G1 family |
| `MagicSquareOfSquares/G1Certificates.lean` | The 114 certificates (generated) |
| `MagicSquareOfSquares/G1Classes.lean` | Remainder classes of `(m, n)` with no rational point: modulo 8, 32, 7, 23, 47 |
| `MagicSquareOfSquares/G1Primes.lean` | Primes `p ≡ 5 (mod 8)` dividing `m² + n²` to an odd power |
| `MagicSquareOfSquares/G1Full.lean` | The G1 family with all nine cells: a finite check that rules out whole ratio classes |
| `MagicSquareOfSquares/G1FullCertificates.lean` | The 484 obstructed classes and the census for `m ≤ 200` (generated) |
| `MagicSquareOfSquares/G1Bridge.lean` | From a magic square of squares to the G1 system; the theorem on edge-root ratios |
| `MagicSquareOfSquares/G1Census1000.lean` | The census for `m ≤ 1000`, evaluated by the kernel in 20 blocks (generated) |
| `MagicSquareOfSquares/G1SizeCheck.lean`, `G1Size.lean` | No two exceptional pairs fit together (generated); the bound on the edge roots |
| `MagicSquareOfSquares/G1Norm.lean` | The conditions on `m² + n²`; why only primes `±1 (mod 8)` obstruct |
| `MagicSquareOfSquares/G1Statement.lean` | The size bound and the conditions on `m² + n²` for the matrix form |
| `MagicSquareOfSquares/G1Centre.lean` | `m² + n²` divides `8e²`; a prime `≡ 1 (mod 8)` divides the centre root unless `m² + n²` is a square or twice a square |
| `MagicSquareOfSquares/G1Relations.lean` | Curve-membership algebra for the Demjanenko–Manin computation, and the partner symmetry |
| `scripts/g1_census.py`, `g1_good_primes.py`, `g1_three_kinds.py` | Further scans behind the archived audit |
| `PAPER.md` | A short note: what is proved, what is computed, what is unchecked |
| `scripts/g1_full/` | The scans behind [G1_FULL_SYSTEM.md](G1_FULL_SYSTEM.md); `genus2/` holds the genus-2 experiment |
| `scripts/search/seven.c`, `seven_brute.c`, `seven_search.txt` | The search for squares with seven square entries, its independent check on a small range, and the summary of the full run |
| `scripts/g1_audit.py` | The scan behind the audit; writes the certificates |
| `scripts/Audit.lean` | Lists every theorem and fails if any uses a non-standard axiom |

## Reproduce

Lean `v4.33.1`, Mathlib commit `0df444a360eaa60ab8c11dca51a86af692955474` (pinned in `lake-manifest.json`).

```
lake exe cache get
lake build
lake env lean scripts/Audit.lean
```

On GitHub's runner a CI run takes about 15 minutes from an empty cache. With the cache, the four runs of 21 and 22 November 2026 took between 3 min 47 s and 5 min 9 s in all, of which the audit step took between 1 min 49 s and 3 min 20 s; on a laptop the audit takes about four minutes. `G1FullCertificates.lean` needs about 5 GB of memory. The last line of the audit has the form `AUDIT OK: N theorems (including auxiliary ones), standard axioms only`, and CI copies that line into the summary of each run. At this commit `N` is 494, against 241 theorems declared in the sources (the rest are auxiliary theorems that Lean generates). The audit checks axioms only. CI also checks that the generated certificate files are the ones their scripts produce, and runs `scripts/g1_full/check_evidence.py` and `scripts/g1_full/octic_check.py`.

## How this was made, and its limits

- The Lean code was written by Claude (an AI model made by Anthropic), directed by Drew (`drew2apps`). The proofs do not rest on trust in either: the Lean kernel checks them.
- The kernel checks proofs, not whether a statement says what was intended. The statements are short. Please read them, and open an issue if one does not match its source.
- The statements follow the Wikipedia list. Rabern's paper was read afterwards and its theorems agree with properties 1, 5, 6, 7, 8 and contain 9. The papers by Zimmermann, Morgenstern and Robertson were not consulted. The proofs here were written independently of all of them.
- Related work found on 21 November 2026 and only skimmed: C. Woll, *A Partial Residue Categorization of the Magic Square of Squares* (arXiv:1809.03067), studies the square modulo the primes dividing the centre root; its Corollary 4.3 says that a prime dividing both the centre and a mid-edge entry of a primitive square is `≡ 1 (mod 8)`. I did not find `m² + n² ∣ 8e²`, or the statement on the centre proved here, in that paper, in Rabern's paper or in the Kepler71 notes. I have not made a full literature search and claim no novelty.
- O. Hill, arXiv:2510.08286 (October 2025), claims a proof that no 3×3 magic square of distinct squares exists. I have not verified it and do not know its status. Wikipedia (read on 21 November 2026) lists the problem as open. Nothing in this repository depends on it.
- Read on 22 November 2026, by going through their statements, not line by line: A. Bremner, *On squares of squares* (Acta Arith. 88, 1999), which constructs families with six squares on elliptic curves and the square with seven; D. Weisenberg, *Some Thoughts on the 3×3 Magic Square of Squares Problem* (Rose-Hulman Undergrad. Math. J. 24, 2023), on the primes dividing the entries; and C. Boyer's page on multimagie.com, which reports Morgenstern's bounds and the congruences of Zimmermann, Pierrat and Thiriet (magic sum `≡ 3 (mod 72)`, entries `≡ 1 (mod 24)`). I did not find `m² + n² ∣ 8e²` or the statement on the centre in any of them. Not read: Robertson (1996), Bremner's second paper (2001), Boyer's article (2005), the paper of Zimmermann, Pierrat and Thiriet itself, and Morgenstern's own pages.

## Sources

- L. Rabern, *Properties of Magic Squares of Squares*, Rose-Hulman Undergraduate Mathematics Journal 4(1).
- P. Zimmermann, *Magic Squares of Squares*, Inria/LORIA note, 2015.
- M. Gardner, *Riddles of the Sphinx*, 1987, pp. 136–137.
- Morgenstern, *Properties of a 3x3 Magic Square of Squares*, multimagie.com, 2007.
- J. P. Robertson, *Magic Squares of Squares*, Mathematics Magazine 69(4), 1996, 289–293.
- A. Bremner, *On squares of squares*, Acta Arithmetica 88, 1999; *On squares of squares II*, Acta Arithmetica 99, 2001.
- N. Bruin, J. Thomas, A. Várilly-Alvarado, *Explicit computation of symmetric differentials and its application to quasihyperbolicity*, Algebra & Number Theory 16(6), 2022, 1377–1405.
- A. Auel, B. Singer, *The algebraic geometry of 3-by-3 magic squares of squares*, arXiv:2609.09351, 2026.
- L. Sallows, *The Lost Theorem*, The Mathematical Intelligencer 19(4), 1997.
- C. Boyer, *Some Notes on the Magic Squares of Squares Problem*, The Mathematical Intelligencer 27(2), 2005, 52–64.

## Licence

Apache-2.0.
