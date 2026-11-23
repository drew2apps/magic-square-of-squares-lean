# A reading of the last step of Hill's preprint

O. Hill, *On arithmetic progressions and a proof of the nonexistence of magic squares of squares*, arXiv:2510.08286, claims that no 3×3 magic square of nine distinct squares exists. If that proof is right, it settles the problem and everything in this repository is a footnote to it. So the preprint has to be read. This file records one reading, of the last step only, with two small programs that anyone can run.

**Status.** Written by Claude (an AI model made by Anthropic) on 23 November 2026, from the HTML of **version 2** on arXiv and the PDF of **version 3** (stamped 7 April 2026, the current one), both read through an automated reader that returned the formulas; no person has compared them with the PDF. In version 3 the reader returned the same equations (11), (12), (30), the same right side of (29) and the same inference; the left side of (29) was legible only in the HTML of version 2. The author has not been contacted. This note does not say that the theorem is false, and it does not say anything about the steps before equation (29). It says that one inference, as transcribed, does not hold, and gives a way to check that in a minute.

## 1. The setting, in the notation used here

Every 3×3 magic square of positive numbers has entries

    v, v+D₁, v+2D₁,   v+2D₁+D₂, v+3D₁+D₂, v+4D₁+D₂,   v+4D₁+2D₂, v+5D₁+2D₂, v+6D₁+2D₂

(equation (15) of the preprint): three progressions of three terms with the same difference `D₁`, whose middle terms are again in arithmetic progression. For a magic square of squares each progression is `x² < y² < z²`. The preprint attaches to a progression the numbers

    n₂ = z − y,    α = D₁ / n₂²,

writes `α = α_n / α_d`, and introduces `N` by `α_n = 4α_d²/N + 3α_d + N/2` (its equations (11), (12)). With indices 1, 2, 3 for the three progressions, `β_{1n}² = α_{2n}/α_{1n}`, `β_{1d}² = α_{2d}/α_{1d}`, `β_{2n}² = α_{3n}/α_{2n}`, `β_{2d}² = α_{3d}/α_{2d}` ((23) to (26)), and `β₁ = β_{1n}/β_{1d}`, so `β₁ = 1` means `α₁ = α₂`.

The proof ends with an equation (29) between `N₂, N₃`, the `β` and `α_{1n}, α_{1d}`, whose right side is rewritten as (30): a factor `(β_{1d}² − β_{1n}²)` times a polynomial of degree 4 in `α_{1d}` with coefficients built from `N₁` and the `β`. The left side of (29) contains only even powers of `α_{1d}`. The preprint then compares coefficients of the odd powers of `α_{1d}`, concludes `β_{1d}² = β_{1n}²`, hence `β₁ = 1`, and from that the three progressions coincide, against distinctness. (This description of the inference comes from two automated readings that agree; it is the part a person should confirm first.)

## 2. What equation (29) is

`scripts/hill/hill_identity.py` (sympy, exact) takes `α_{1d}, α_{2d}, α_{3d}, N₁, N₂, N₃` as six independent indeterminates, defines `α_{in}` by (12) and the `β²` by (23) to (26), and finds, as identities of rational functions:

    RHS(29) − RHS(30) = 0,

    LHS(29) − RHS(29) = 4·N₂²·N₃²·α_{3n}·α_{3d} · ( 2·h(α₂) − h(α₁) − h(α₃) ),        h(a) = a + 1/a,

    bracket of (30) = β_{1n}²·α_{1n}² − β_{1d}²·α_{1d}².

The first line says that the two transcribed formulas are consistent with each other, which a mistranscription of either would almost certainly break. The second says what (29) is: since the factor in front is not zero,

> (29) holds exactly when `h(α₁), h(α₂), h(α₃)` are in arithmetic progression.

And `h(α) = 4y²/D₁ + 2` (from `y = n₂(α − 1)/2`, `D₁ = α·n₂²`). So (29) is the statement that the middle terms `y₁², y₂², y₃²` of the three progressions are in arithmetic progression: the magic condition itself, written in the variables `α`. It is one equation in three unknowns. Nothing in it forces `α₁ = α₂`.

The third line explains the odd powers: the bracket of (30) is `β_{1n}²α_{1n}² − β_{1d}²α_{1d}²` with `α_{1n}` expanded by (12). Its coefficients contain `N₁`, which is tied to `α_{1d}` and `α_{1n}` by (12), and the `β`, which are ratios involving `α_{1n}` and `α_{1d}`. The equation is a relation between numbers that depend on each other, not an identity in a free variable `α_{1d}`, so coefficients cannot be compared.

## 3. A test on real magic squares, and what it is not

The preprint works with integers: the progressions have integer lengths, and `α = α_n/α_d` with `α_n, α_d` positive integers. The formulas of section 1 still make sense for positive real numbers. `scripts/hill/hill_check.py` evaluates them, with 40 digits, on magic squares of positive real numbers, read as squares of their square roots, with an arbitrary choice of each `α_d`:

- the square with entries 1, …, 9 (progressions 1, 2, 3; 4, 5, 6; 7, 8, 9): `α = 5 + 2√6, 11 + 2√30, 17 + 12√2`, `h(α) = 10, 22, 34`;
- five random ones.

In all six, (29) and (30) hold to 30 digits, and `β₁` is 1.489, 1.669, 1.391, 1.275, 1.072, 1.436 (`hill_check.txt`). **What this is not.** The square of 1, …, 9 is not a counterexample to the theorem, and it does not satisfy the hypotheses of the preprint: its `α` are irrational, and the program uses real `α_d` and `N` where the preprint has integers. An example that satisfies all the hypotheses would be a magic square of squares.

**What it shows.** The relations that the last step uses, (11), (12), (23) to (26) and (29), (30), read as relations between positive real numbers, do not imply `β₁ = 1`. So the step cannot be a consequence of those relations alone; it would have to use, at that point, that the quantities are integers. The readings found no such use after Lemma 3.2 of the preprint: positivity and the comparison of coefficients only.

**And the comparison fails for integers too, for a reason that needs no example.** For the fixed integers of a hypothetical square, (29) is one equality between numbers. `α_{1d}` is one of those numbers, and `N₁`, `β_{1n}`, `β_{1d}`, which appear in the coefficients, are tied to it by (12) and (23), (24). An equality that holds at one value of a quantity says nothing about the coefficients of a polynomial in that quantity. This is the objection of section 2, and it does not depend on section 3. We looked for integer data that would make it concrete, three progressions of rational squares with (29) but not the same difference (section 5): none exists with shapes up to 2000, so we have no such example to offer, and say so.

## 4. What this does and does not show

- It shows that, in versions 2 and 3 as transcribed, the last inference is not valid as written, and that equation (29) carries no more information than the magic condition from which it was derived.
- It does not show that the theorem is false, or that the approach cannot be completed.
- It rests on formulas returned by an automated reader. The identities of section 2 would be very unlikely to hold for a wrong transcription, but the description of *how the preprint uses* (29) and (30) is the reader's, and should be confirmed against the PDF by a person. Three automated readings (two of version 2, one of version 3) describe it the same way: coefficients of odd powers of `α_{1d}` are compared, and integrality is not used at that step.
- There is also an argument that does not depend on the exact form of (29): every step from the definitions to the conclusion `β₁ = 1` is algebra that is valid for positive real numbers, and for the square of 1, …, 9 the conclusion is false. So some step must fail, whatever the formulas are, unless one of them uses integrality; and the readings found none that does.
- Until a person has confirmed the reading, the README keeps its sentence that the preprint has not been verified here. What changes is that there is now a specific, checkable reason for not relying on it.

## 5. A weaker problem that the same variables suggest

Dropping the condition that the three progressions have the same difference leaves: three arithmetic progressions of three rational squares, with differences `D⁽ⁱ⁾` and middle terms `yᵢ²`, such that `y₁²/D⁽¹⁾, y₂²/D⁽²⁾, y₃²/D⁽³⁾` are distinct and in arithmetic progression. Every magic square of squares gives one (in two ways, one for each direction of its broken diagonals), and the condition does not change when a progression is multiplied by a square, so it depends only on the three shapes `(p, q)` of the progressions, `y²/D = (p² + q²)² / (4pq(p² − q²))`. Said with right triangles: `y` is the hypotenuse and `D/4` the area of a rational right triangle with legs `a, b`, and `4y²/D = 2(a/b + b/a)`; so the weaker problem asks for three Pythagorean triples whose values of `a/b + b/a` are distinct and in arithmetic progression, and a magic square of squares needs, in addition, that the three areas agree up to squares (as far as we know this is the formulation of Robertson, 1996, which was not read here). Two programs look for solutions among all shapes with `q < p ≤ P`, each with a planted progression as a self-test. `scripts/hill/ap_relaxed.py` was the first (its run to 1000 is `ap_relaxed_1000.txt`): it finds candidates with floating-point numbers and confirms them with exact fractions, so it cannot report a false solution, but nothing proves that it cannot miss one. `scripts/hill/ap_exact.c` has no floating point: every value is reduced modulo the prime `2⁶¹ − 1`, where `w₁ + w₃ = 2w₂` implies the same congruence between residues, the residues of distinct values are checked to be distinct, and the sum of every pair is looked up in a table of the doubled residues; a solution cannot be missed, and a candidate would be printed for exact confirmation. **There is no solution, and not even a candidate, for `P = 2000`**: 811,155 shapes, 328,985,811,435 pairs (`ap_exact_2000.txt`, about 40 minutes on two cores, not rerun in CI; `ap_exact_1000.txt`, 202,861 shapes, is rerun in the full CI and `ap_exact_250.txt` on every push). So in a magic square of squares, the three progressions of either direction cannot all have shape `p ≤ 2000`. This is a bounded search, exhaustive within its bound as far as the program is right; it says nothing beyond the bound. This is a search, with the usual caveats, in a region that the searches by size of the entries do not cover: the scales of the three progressions are free.

## 6. Reproduce

```
cd scripts/hill
python3 hill_identity.py     # the three identities, exact (sympy); hill_identity.txt
python3 hill_check.py        # the six real squares (mpmath); hill_check.txt
gcc -O2 -fopenmp -o ap_exact ap_exact.c && ./ap_exact 1000    # the weaker problem, exact; about a minute on two cores; ap_exact_1000.txt
python3 ap_relaxed.py 250    # the first version (floating-point candidates, exact confirmation); ap_relaxed_250.txt
```
