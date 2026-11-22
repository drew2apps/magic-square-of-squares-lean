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
| If also `f ≠ h`, let `m′ : n′` be the ratio `(b + d) : (f − h)` in lowest terms. Then `8e² = (m² + n²)·k²·(m′² + n′²)` with `k ≠ 0`; in terms of the entries, `((b + d)² + (h + f)²)·((b + d)² + (f − h)²) = 8e²(b + d)²` | `IsMSS.center_factorization`, `IsMSS.edge_sums_identity` |
| `2(m² + n²)(m′² + n′²)` is a square; exactly one of the two numbers is even; every odd prime divides them to exponents of the same parity | `IsMSS.two_mul_norms_isSquare`, `IsMSS.norms_parity`, `IsMSS.norms_padicVal_parity` |
| Either a prime `≡ 1 (mod 8)` divides `e` and divides both numbers to an odd power, or one of them is a square and the other is twice a square | `IsMSS.centre_two_ratios` |
| `m : n` lies in none of 484 residue classes at 34 primes | `IsMSS.ratio_not_in_class`, `G1.fullTable_valid` |
| In particular `m : n` is not one of the coprime pairs `1 ≤ n < m ≤ 1000` outside 425 listed pairs (60 of them up to 200) | `IsMSS.ratio_ne_of_le_1000` |
| The same for the five forms of the G1 family: no rational `t` makes them squares | `G1.no_solution_le_1000` |
| If `b² ≠ d²`, `h² ≠ f²` and `b² ≠ f²`, one of the edge roots `b, d, f, h` is larger than 500 in absolute value; so one of the four edge entries is larger than 250,000, in particular when the nine entries are distinct | `IsMSS.edge_root_gt_500`, `IsMagicSquareOfSquares.edge_gt_of_ne`, `edge_gt` |
| For a prime `p ≡ ±3 (mod 8)` no ratio is excluded by a congruence modulo `p` (powers of `p` are not covered) | `G1.exists_point_mod_of_not_isSquare_two`, `G1.exists_point_mod_of_three_five_mod_eight` |

The first three rows have no bound on the size of the entries. The signs of the roots are arbitrary, so the theorems also apply with `b − d` or `h − f`. The ratios excluded by rows 1 and 3 have density 1 among all ratios (argued in G1_FULL_SYSTEM.md, not in Lean); this does not decide the problem. The bound 500 is far below what computer searches have covered.

Computed, not in Lean: exactly 34 primes can obstruct the system modulo `p`, the largest is 607. Of the 425 listed pairs, 268 are closed by elliptic curves of rank 0, 121 more by the Demjanenko–Manin argument on three curves of rank 1, and 36 are open (5 of them up to 200). These computations depend on PARI and Sage (with eclib), on Magma (it confirms every upper rank bound, and two pairs rest on it alone), and on derivations that nobody has checked; `scripts/g1_full/dm_evidence.txt` records each one.

Also computed, not in Lean (`scripts/search/seven.c`, summary in `scripts/search/seven_search.txt`): among the 3×3 magic squares with nine distinct positive entries, a square centre of root at most 10⁹ and at least seven square entries, the only ones are Bremner's square (1999) and its multiples; none has eight or nine squares. The reduction: a line through the centre with both ends square is a Pythagorean triangle with hypotenuse the centre root; seven squares force two such lines; the two lines can be any two of the four, which gives six cases for the other two deviations, all tested for every pair of triangles (the first version of the program missed two of the six; that was found in review and the whole range was run again). A second program without a table of cases, `seven_brute.c`, agrees up to centre root 1300 (output in `seven_brute_1300.txt`; CI runs it to 430). Scope: squares with seven square entries and a centre that is not a square are not covered. None of this is new. Ben Chaffin's page (benchaffin.com, published in 2019, work of about 2008) reports a search of the same family (a square centre, and pairs of decompositions of twice the centre into two squares) up to centre entry 10²⁸, that is centre root 10¹⁴, and of non-square centres up to 10¹⁴, with no square with seven or more square entries other than Bremner's (read from a PDF of the page with an AI model; for a square centre the page describes the same six placements of two complete lines as the program here). Pierrat, Thiriet and Zimmermann (2015) report searches of particular seven-square configurations and cite Buell for a centre above 25·10²⁴ in the nine-square case; I do not quote their ranges, because I have them only from an AI summary of the preprint and two figures in it disagree. Morgenstern showed in 2006 that all nine entries must exceed 10¹⁴ and searched triples of arithmetic progressions of squares with common difference up to 2.4·10²¹ (as reported on multimagie.com). The search here re-verifies a small part of the range reported by Chaffin with an independent program. In Lean (`SevenSquares.lean`): a complete line is a Pythagorean triangle with hypotenuse the centre root and its deviation is twice the product of the legs (`full_line_legs`), and the six cases are all the cases (`two_lines_exhaust`, `six_cases`; the second case has two alternatives for the sign, so there are seven pairs of formulas); Bremner's square is in `Examples.lean`. The program and its result are not in Lean. `seven_search.txt` records the commands and the SHA-256 of the source file that produced it, and CI checks that `seven.c` still has that hash. The run of `seven_brute.c` to 1300 is not repeated in CI, which runs it to 430.

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
| `MagicSquareOfSquares/G1Centre.lean` | `m² + n²` divides `8e²`; with the second ratio, `2(m² + n²)(m′² + n′²)` is a square, and a prime `≡ 1 (mod 8)` divides the centre root and both numbers to an odd power unless one of them is a square and the other twice a square |
| `MagicSquareOfSquares/G1Relations.lean` | Curve-membership algebra for the Demjanenko–Manin computation, and the partner symmetry |
| `MagicSquareOfSquares/SevenSquares.lean` | Algebra behind the search for seven squares: the triangle of a complete line; the six cases are all the cases |
| `scripts/g1_census.py`, `g1_good_primes.py`, `g1_three_kinds.py` | Further scans behind the archived audit |
| `AUDIT_REQUEST.md` | What an outside reader is asked to check, with three closures worked out |
| `PAPER.md` | A short note: what is proved, what is computed, what is unchecked |
| `scripts/g1_full/` | The scans behind [G1_FULL_SYSTEM.md](G1_FULL_SYSTEM.md); `genus2/` holds the genus-2 experiment |
| `scripts/search/seven.c`, `seven_brute.c`, `seven_search.txt`, `seven_brute_1300.txt` | The search for squares with seven square entries, its independent check on a small range, and the summary of the full run |
| `scripts/g1_audit.py` | The scan behind the audit; writes the certificates |
| `scripts/Audit.lean`, `scripts/audit_output.txt`, `scripts/audit_line.txt` | Lists every theorem and fails if any uses a non-standard axiom; its output for the current sources, and the last line of it |

## Reproduce

Lean `v4.33.1`, Mathlib commit `0df444a360eaa60ab8c11dca51a86af692955474` (pinned in `lake-manifest.json`).

```
lake exe cache get
lake build
lake env lean scripts/Audit.lean
```

On GitHub's runner a CI run takes about 15 minutes from an empty cache. With the cache, the runs of 21 and 22 November 2026 took between about 4 and 15 minutes (it depends on the state of the cache and on the runner); on a laptop the audit alone takes about four minutes. `G1FullCertificates.lean` needs about 5 GB of memory. The last line of the audit has the form `AUDIT OK: N theorems (including auxiliary ones), standard axioms only`, and CI copies that line into the summary of each run. The whole output for the current sources is committed as `scripts/audit_output.txt` and its last line as `scripts/audit_line.txt`; CI fails if the audit prints anything else, and shows the line as a notice on the page of each run. `N` was 494 at release 0.3.3 (commit 5707ec8), 500 at commit 53d2c89 and 578 at commit 9ce6fe5; it is 592 at release 0.3.4, against 257 theorems declared in the sources (the audit counts the theorems of the Lean environment, which include auxiliary theorems that Lean generates). The audit checks axioms only. CI also checks that the generated certificate files are the ones their scripts produce, and runs `scripts/g1_full/check_evidence.py`, `scripts/g1_full/octic_check.py` and `scripts/g1_full/height_sides.py` (exact rational arithmetic), and compares their output with the committed `octic_check.txt` and `height_sides.txt`. `scripts/g1_full/height_precision.py` needs Sage and is not run in CI; its output is committed as `height_precision.txt`.

## How this was made, and its limits

- The Lean code was written by Claude (an AI model made by Anthropic), directed by Drew (`drew2apps`). The proofs do not rest on trust in either: the Lean kernel checks them.
- The kernel checks proofs, not whether a statement says what was intended. The statements are short. Please read them, and open an issue if one does not match its source.
- The statements follow the Wikipedia list. Rabern's paper was read afterwards and its theorems agree with properties 1, 5, 6, 7, 8 and contain 9. The papers by Zimmermann and Robertson and Morgenstern's own pages were not consulted. The proofs here were written independently of all of them.
- Literature, read with an AI model from the PDFs on 21 and 22 November 2026. In full: Rabern, Woll (arXiv:1809.03067), Hill (arXiv:2510.08286). Statements and the relevant sections: Weisenberg (2023), Bremner (1999, 2001), Auel and Singer (arXiv:2609.09351), Kominers, Bruin, Thomas and Várilly-Alvarado, Várilly-Alvarado's survey, and C. Boyer's page on multimagie.com, which reports Morgenstern's bounds and the congruences of Zimmermann, Pierrat and Thiriet. Also read: Ben Chaffin's page, and the statements used from Cremona, Prickett and Siksek (2006). Not read: Pierrat, Thiriet and Zimmermann (2015; known to me only through an AI summary), Silverman (1990), Robertson (1996), Boyer's article (2005), Buell's preprint, Morgenstern's own pages, Labruna. The comparison is in [G1_FULL_SYSTEM.md](G1_FULL_SYSTEM.md), Result 1b.
- On the statement about the centre. A prime is `≡ 1 (mod 8)` exactly when `−1` and `2` are both squares modulo it, that is, when it splits completely in `ℚ(i, √2) = ℚ(ζ₈)`; Rabern's proofs work in `ℤ[i]` and `ℤ[√2]`. Published statements go from a prime to its residue class: every prime divisor of the centre is `≡ 1 (mod 4)` (Rabern, Woll); a prime dividing the centre and a mid-edge entry is `≡ 1 (mod 8)` (Woll); a prime dividing only the centre has three consecutive nonzero quadratic residues (Weisenberg, Woll). The statements here say where such primes come from (all in Lean, see the table above): the odd primes dividing `m² + n²` and `m′² + n′²` divide the centre root, those dividing to an odd power are `≡ 1 (mod 8)` and are the same for the two numbers, and if there are none, one number is a square and the other twice a square. Behind this is that `2(m² + n²)(m′² + n′²)` is a square, which is the condition that the centre is a square, written multiplicatively; searches already impose it, so this is not a new search method. I did not find `m² + n² ∣ 8e²` or this statement in the texts listed above. I do not know whether the exceptional case can occur, and I claim no novelty while the texts listed as not read remain unread.
- The family of curves used here is a pencil of hyperplane sections of the surface studied by Auel and Singer: the hyperplane `n(b + d) = m(h + f)` cuts their surface in the curve of the ratio `m : n` and 16 conics of squares with repeated entries (my remark, checked by hand; see G1_FULL_SYSTEM.md).
- O. Hill, arXiv:2510.08286 (version 3, April 2026), claims a proof that no 3×3 magic square of distinct squares exists. It has not been verified here, and this repository takes no position on it. Wikipedia (read on 21 November 2026) lists the problem as open. Nothing in this repository depends on it.

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
- C. Woll, *A Partial Residue Categorization of the Magic Square of Squares*, arXiv:1809.03067, 2018.
- D. Weisenberg, *Some Thoughts on the 3×3 Magic Square of Squares Problem*, Rose-Hulman Undergraduate Mathematics Journal 24(1), 2023.
- S. D. Kominers, *Center-zero magic squares of squares over number fields*, preprint.
- A. Várilly-Alvarado, *The geometric disposition of Diophantine equations*, survey article.
- J. E. Cremona, M. Prickett, S. Siksek, *Height difference bounds for elliptic curves over number fields*, Journal of Number Theory 116, 2006, 42–68.
- J. H. Silverman, *The difference between the Weil height and the canonical height on elliptic curves*, Mathematics of Computation 55, 1990, 723–743.

## Licence

Apache-2.0.
