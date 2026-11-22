# The 3×3 magic square of squares: a machine-checked reduction and a computed census

Drew Carlson. Version 0.3.3, 22 November 2026.
Repository: https://github.com/drew2apps/magic-square-of-squares-lean. Archive: https://doi.org/10.5281/zenodo.24435351.

The Lean code and this text were written with Claude, an AI model made by Anthropic. I am not a number theorist. Nobody else has checked the mathematics outside Lean. This note says what is proved, what is only computed, and what is unchecked, so that someone who knows the subject can check it.

## 1. The problem

Is there a 3×3 magic square whose nine entries are distinct perfect squares? The question is open. Write the square of roots as

```
a b c
d e f
g h i
```

so that the rows, columns and two diagonals of `a², …, i²` have the same sum `3e²`. Lee Morgenstern showed by computation (2006) that every entry of a solution exceeds 10¹⁴, and searched much further in other parametrisations; Pierrat, Thiriet and Zimmermann (2015) cite Buell for a centre above 25·10²⁴. The only known square with seven square entries is Bremner's (1999). A 2025 preprint (O. Hill, arXiv:2510.08286) claims a proof of nonexistence; I have not verified it, and nothing here depends on it.

This note does not solve the problem. It contributes three things: a Lean library of the known facts, a reduction of the problem to a one-parameter family of curves with a kernel-checked list of excluded cases, and a computed census of the remaining cases.

## 2. Proved in Lean

Lean 4 with Mathlib, 241 declared theorems, no `sorry`, no `native_decide`, standard axioms only (checked in CI for every theorem). Statements exist both for nine integers and for the matrix form used in google-deepmind/formal-conjectures.

**Known facts.** The nine necessary conditions listed on Wikipedia (most of them are in Rabern 2003), and the classical reformulations: three arithmetic progressions of squares with a common difference, rational right triangles of equal area, and points on the curves `y² = x³ − n²x`.

**The ratio of the edges.** Let `m : n` be the ratio `(b + d) : (h + f)` in lowest terms. Then there are integers `u`, `v` with

`b + d = mu`, `h + f = nu`, `b − d = nv`, `f − h = mv`, `8e² = (m² + n²)(u² + v²)`,

and the square is determined by `m`, `n` and a rational parameter `t`, for which five quadratic forms must be squares: `F0 = m² + n²t²`, `F8 = n² + m²t²`, `G4 = 2(m² + n²)(t² + 1)`, `GL = G4 − 8mnt`, `GU = G4 + 8mnt`. The pair `(m + n, m − n)` (halved if both are even) gives the same system; this is the *partner*.

**Excluded ratios.** 484 residue classes of `m : n` at 34 primes (the largest is 607) are impossible; each class is certified by a finite computation that the Lean kernel checks.

**The norm.** `m² + n²` is a square modulo 7, 17, 31 and 47; every prime `≡ 5 (mod 8)` divides it to an even power; and `m² + n²` divides `8e²`. Hence: either `m² + n²` is a square or twice a square, or a prime `≡ 1 (mod 8)` divides the centre root `e` (and divides `m² + n²` to an odd power). I have not found this statement in Rabern (2003), Bremner (1999), Woll (2018), Weisenberg (2023) or the multimagie pages, but I have not made a full literature search and claim no novelty.

**Size.** One of the four edge roots exceeds 500. This is far below Morgenstern's bound; it is included because it is machine-checked.

**Census.** For coprime `1 ≤ n < m ≤ 1000` (304,191 ratios), no solution exists outside an explicit list of 425 ratios. The kernel evaluates the census.

**Algebra for Section 3.** Identities used by the computations below are checked in Lean: membership of the points on the three elliptic curves, the map from the quartic `v² = F0·GL` to a Weierstrass model, the partner identities, and the identities behind the singular fibres.

## 3. Computed, not formalized

Details, limits and commands are in `G1_FULL_SYSTEM.md`. Programs: PARI 2.17, Sage (passagemath) with eclib, Magma 2.29 through its public calculator, and C and Python programs in the repository.

**Which primes can obstruct.** For a ratio of good reduction the five conditions define a smooth curve of genus 49 over `F_p`. The Weil bound, a separate character-sum bound for the fibres where two branch points collide (`p` dividing `m⁸ + 14m⁴n⁴ + n⁸`), and a direct computation below those bounds show that exactly the 34 primes above have an obstructed class. All 34 are `±1 (mod 8)`; there is a short proof that no other residue is possible, and that part is in Lean.

**Local solubility.** The 425 ratios are exactly those not excluded by the classes and by the condition on primes `≡ 5 (mod 8)`. A second program, written separately, decides solubility in every `ℚ_p` and agrees on all 304,191 ratios. So no congruence closes any of the 425. The excluded ratios have density 1 (argued, with a proof sketch).

**Rank 0.** 268 of the 425 ratios are closed by a quotient elliptic curve of rank 0.

**Rank 1 with multiplicity.** Three quotient elliptic curves occur more than once in the Jacobian: `E04: y² = x(x + 2Sm²)(x + 2Sn²)` with `S = m² + n²` (twice), a second curve `D1` (twice), and `E0L: y² = F0·GL` (four times). When such a curve has rank exactly 1, a solution gives two points whose heights are tied to each other, and the Demjanenko–Manin argument leaves finitely many candidates, which are tested exactly. This closes 121 more ratios (91, 24 and 6 on the three curves). For every one of them the bound "rank at most 1" is given by Magma (87 distinct curves); for all but two it is also given by PARI, and for most of the first curve by eclib. For two ratios, (631,71) and its partner, PARI gives only `[1, 3]` and the bound rests on Magma alone.

**Open.** 36 ratios, 5 of them with `m ≤ 200`: (72,65)~(137,7), (91,60)~(151,31), (195,28). Seven have no trivial point. For these the curves above have rank at least 2. A quotient of genus 2 was tried through the Magma calculator; for 20 of its 23 instances the calculator's limits were reached.

**A side search.** Among magic squares with nine distinct positive entries, a square centre of root at most 10⁹ and at least seven square entries, the only ones are Bremner's square and its multiples. This re-verifies, on a smaller range and with a different program, searches reported by Pierrat, Thiriet and Zimmermann (2015) for several configurations; it is not a new bound.

## 4. What has not been checked

In order of importance:

1. **Heights and generators.** The Demjanenko–Manin closures need the generator of each curve (saturation with a proved height lower bound) and explicit comparisons between canonical and naive heights. The derivations were written with an AI model. The heights are floating-point with a safety margin, not interval arithmetic. An automatic checker verifies the recorded points, models, candidates and bounds for consistency, but not the heights or that the points generate.
2. **The character-sum bound** for the singular fibres and the Weil-bound step.
3. **The density argument**, which is a sketch.
4. **Rank bounds** are computed, by two or three programs, not proved here.
5. **The literature.** Robertson (1996), Bremner's second paper (2001) and Boyer's article (2005) have not been read; the preprint of Pierrat, Thiriet and Zimmermann (2015) only through an AI summary.

Two errors were found and corrected during the work, both recorded in `CHANGES.md`: the singular fibres were overlooked in the first statement of the 34 primes, and the first version of the side search missed two of six cases. Others may remain.

## 5. Reproducing

`lake build` and `lake env lean scripts/Audit.lean` check the Lean part (about 15 minutes from an empty cache). CI also regenerates the certificate files, runs the evidence checker, the check of the singular fibres and the side search on a small range. The elliptic-curve computations are rerun with the scripts in `scripts/g1_full/`; their outputs are committed.

## References

- A. Bremner, On squares of squares, Acta Arith. 88 (1999), 289–297.
- C. Boyer, 3×3 magic square of squares search, multimagie.com (read 22 November 2026); reports Morgenstern's results.
- P. Pierrat, F. Thiriet and P. Zimmermann, Magic squares of squares, preprint, 2015.
- L. Rabern, Properties of magic squares of squares, Rose-Hulman Undergrad. Math. J. 4 (2003).
- J. P. Robertson, Magic squares of squares, Math. Mag. 69 (1996), 289–293 (not read).
- D. Weisenberg, Some thoughts on the 3×3 magic square of squares problem, Rose-Hulman Undergrad. Math. J. 24 (2023).
- C. Woll, A partial residue categorization of the magic square of squares, arXiv:1809.03067.
- A. Auel and B. Singer, The algebraic geometry of 3-by-3 magic squares of squares, arXiv:2609.09351.
- Kepler71, notes on the full nine-cell system (cited in `G1_FULL_SYSTEM.md`).
