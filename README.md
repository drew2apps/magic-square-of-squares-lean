# Magic square of squares: known facts, checked in Lean

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

Properties 3, 8 and 9 do not need the square to be primitive. Theorems 7 and 9 are proved for primes ≡ 3 (mod 8) as well.

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

In the three reformulations the nine squares are required to be distinct. The informal statements leave this implicit. The curve uses Mathlib's `WeierstrassCurve.Affine.Point` and its group law.

### Checks against published examples

| Example | Theorem |
|---|---|
| Bremner's square is magic and has exactly two non-square entries | `bremner_isMagic`, `bremner_not_isSquare` |
| Sallows' square: rows, columns and one diagonal sum to 147²; the other diagonal does not | `sallows` |
| The Parker square: rows, columns and one diagonal sum to 3051; the other diagonal sums to 4107 | `parker` |
| Euler's 4×4 magic square of distinct squares (1770) | `exists_magic_square_squares_four` |

## What is not proved

- Whether a 3×3 magic square of distinct squares exists.
- Any search bound, and any result that needs the arithmetic of curves (Bremner 1999 and 2001, Bruin–Thomas–Várilly-Alvarado, Auel–Singer).
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
| `scripts/Audit.lean` | Lists every theorem and fails if any uses a non-standard axiom |

## Reproduce

Lean `v4.33.1`, Mathlib commit `0df444a360eaa60ab8c11dca51a86af692955474` (pinned in `lake-manifest.json`).

```
lake exe cache get
lake build
lake env lean scripts/Audit.lean
```

The audit takes about three minutes and ends with `AUDIT OK`.

## How this was made, and its limits

- The Lean code was written by Claude (an AI model made by Anthropic), directed by Drew (`drew2apps`). The proofs do not rest on trust in either: the Lean kernel checks them.
- The kernel checks proofs, not whether a statement says what was intended. The statements are short. Please read them, and open an issue if one does not match its source.
- The statements follow the Wikipedia list. The original papers by Rabern, Zimmermann, Morgenstern and Robertson were not consulted, and the proofs here were written independently of them.

## Sources

- L. Rabern, *Properties of Magic Squares of Squares*, Rose-Hulman Undergraduate Mathematics Journal 4(1).
- P. Zimmermann, *Magic Squares of Squares*, Inria/LORIA note, 2015.
- M. Gardner, *Riddles of the Sphinx*, 1987, pp. 136–137.
- Morgenstern, *Properties of a 3x3 Magic Square of Squares*, multimagie.com, 2007.
- J. P. Robertson, *Magic Squares of Squares*, Mathematics Magazine 69(4), 1996, 289–293.
- A. Bremner, *On squares of squares*, Acta Arithmetica 88, 1999.
- L. Sallows, *The Lost Theorem*, The Mathematical Intelligencer 19(4), 1997.
- C. Boyer, *Some Notes on the Magic Squares of Squares Problem*, The Mathematical Intelligencer 27(2), 2005, 52–64.

## Licence

Apache-2.0.
