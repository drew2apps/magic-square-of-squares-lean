# The rational distance problem for the unit square, line by line

**Status.** An experiment, started on 24 November 2026 with the descent of `LAST_PAIR.md`, sections 11 to 14. Written by Claude (an AI model made by Anthropic), directed by Drew (`drew2apps`). Nobody has read it. It is not part of the census of magic squares and nothing here is in Lean. **Most of the structure used here is not new: see "Prior work" below, which was written after the first version of this note.**

**The problem.** Is there a point of the plane at rational distance from the four corners of the unit square? It is open (Guy, *Unsolved Problems in Number Theory*, D19). Known: there is no such point on the edges, the midlines or the diagonals of the square, and points at rational distance from three corners are plentiful (T. G. Berry, *Points at rational distance from the corners of a unit square*, Ann. Scuola Norm. Sup. Pisa (4) 17 (1990), 505–529; not read by us).

## Prior work

The first version of this note (commit `9d59c67`) said that the literature had not been searched. A reviewer then pointed to a preprint, and we read the parts of it that overlap:

**Y. Si, *Arithmetic constraints and covering curves for rational distances in a square*** (University of Waterloo; version 1.0.0 of 25 October 2026, version 1.1.0 of 14 November 2026, doi:10.5281/zenodo.24212887, sources at `github.com/SamSi0322/fourcorner-arithmetic-geometry`; read from the LaTeX sources of version 1.1.0, sections 1 and 12, by an AI model; the rest of its 83 statements was not read). For a point strictly inside the square, with `L` the common denominator of its coordinates:

- *Theorem "Canonical corner classes"* (his section 1): `z − c_j` is `L·m_j` times a square of `ℚ(i)`, where `m_0, …, m_3` are positive, odd, squarefree, **pairwise coprime**, prime to `L`, with all prime factors `≡ ±1 (mod 8)`. These are the four classes of section 1 below, in a sharper form than ours. His proof is the same Gaussian factorisation.
- The *kernel* `n = m_0m_1m_2m_3` is the squarefree part of `t² − 1`, `t = r_0r_3 + r_1r_2`, and `n > 1` (used in the proof of his Corollary "A uniform rank requirement"). In his terms **our Proposition 1 is the statement `n ≠ 1`**: four equal classes are four equal `m_j`, hence all equal to 1.
- *Corollary "A uniform rank requirement"* (his section 12, new in version 1.1.0; a written proof, not among the statements he has formalised in Lean): every such point gives two independent points of `E_n: Y² = X³ − n²X`, so `rank E_n(ℚ) ≥ 2`. This is stronger than the remark on congruent number curves in section 1 below.

Also relevant and not read by us: Berry's paper above (a Kummer surface and curves of genus five, according to Si's introduction), and A. Bremner and M. Ulas, *Points at rational distances from the vertices of certain geometric objects* (arXiv:1502.07312).

**So what is left of this note as possibly additional:** (a) the census of section 3, line by line: we did not find a statement of this kind in the parts of Si's paper that we read; (b) a second proof of `n ≠ 1`, over `ℚ(i)` (Proposition 1); (c) a test of Si's criterion on curves of rank 1 and the remark that, with his corollary, the one line our census leaves open is closed (section 4). No priority is claimed for any of it, and the comparison rests on a partial reading of one preprint.

## 1. The four classes

A point with four rational distances has rational coordinates (the differences of the squared distances give `x` and `y`). Put `z = x + iy`. The four distances are rational exactly when

    z,   z − 1,   z − i,   z − 1 − i

have square norms. A nonzero element of `ℚ(i)` of square norm is a positive squarefree integer, its *class*, times a square of `ℚ(i)`, and the class is unique (`LAST_PAIR.md`, Lemma S1). So a solution has four classes `D₁, D₂, D₃, D₄`, up to a common rational factor (the denominator of `z`). In Si's notation `D_j = L·m_j` up to squares.

**Proposition 1.** There is no `z ∈ ℚ(i)`, other than the corners, for which the four numbers `z, z − 1, z − i, z − 1 − i` are in the same class of `ℚ(i)*` modulo squares. In particular, at a point with four rational distances the four classes are not all equal.

*Proof.* Write `K = ℚ(i)`. Suppose `z − a_j = D·w_j²` with one `D ∈ K*` and `w_j ∈ K`, where `a_1, …, a_4 = 0, 1, i, 1 + i`. Then `w₁² − w₂² = 1/D`, `w₁² − w₃² = i/D`, `w₁² − w₄² = (1 + i)/D`, and eliminating `1/D`:

    C:   w₃² = a·w₁² + b·w₂²,     w₄² = c·w₁² + d·w₂²,        a = 1 − i,  b = i,  c = −i,  d = 1 + i.

*Step 1, a map to an elliptic curve.* If `w₂ ≠ 0`, put `X = (w₁/w₂)²` and `Y = w₁w₃w₄/w₂³`. Multiplying the two equations of `C` and then by `w₁²/w₂⁶` gives

    E′:   Y² = X(aX + b)(cX + d).

*Step 2, the curve.* With `κ = −(1 + i)`, `β = 1 − i`, `x = κ(X − β/2)`, `y = κY`, the equation `E′` becomes `E: y² = x³ − x` (an identity of polynomials: `κ²·X(aX + b)(cX + d) = x³ − x`).

*Step 3, the points of `E` over `K`.* The rank of `E` over `K` is the rank of `E` over `ℚ` plus the rank of its twist by `−1`; the twist `−y² = x³ − x` is `E` again (`x ↦ −x`), and `E(ℚ)` has rank 0 (Fermat: 1 is not a congruent number). So `E(K)` is finite. Its torsion: `K` is unramified at 5 and at 13, with residue fields `F₅` and `F₁₃` at the primes above them, and `E` has good reduction there; reduction is injective on torsion of order prime to the residue characteristic, and `E` has 8 points over `F₅` and 8 over `F₁₃`. A torsion point of order prime to 5 thus has order dividing 8, and the 5-part of the torsion embeds in a group of order 8 and is trivial. So `E(K)` has at most 8 points, and it has these 8: `O`, `(0, 0)`, `(±1, 0)`, `(i, ±(1 − i))`, `(−i, ±(1 + i))`.

*Step 4, back to `C`.* `X = (β/2)(1 − x)`, so the `x`-coordinates `∞, 0, 1, −1, i, −i` give `X = ∞, β/2, 0, β, −i, 1`. But `X = (w₁/w₂)²` is a square of `K`, and `β/2`, `β`, `−i` are not squares of `K` (the norms of the first two are `1/2` and `2`, and `−i = u²` has no solution in `K`). `X = 0` means `w₁ = 0` and `w₃² = i·w₂²`, and `w₂ = 0` (the case left out in step 1, `X = ∞`) means `w₃² = (1 − i)·w₁²`: impossible, since `i` and `1 − i` are not squares of `K`. So `X = 1`: `w₁² = w₂²`, which contradicts `w₁² − w₂² = 1/D`. ∎

Every computation of this proof is checked in exact arithmetic by `scripts/rational_distance/rd_prop1.py` (output `rd_prop1.txt`); no Jacobian and no elliptic-curve package is used (sympy expands the two polynomial identities; the rest is arithmetic with fractions). What the proof takes from the literature: the rank of `y² = x³ − x` over `ℚ` is 0; the rank over a quadratic field is the sum of the ranks of the curve and of its twist; reduction at a prime of good reduction is injective on torsion of order prime to the residue characteristic (Silverman, *The Arithmetic of Elliptic Curves*, VII.3.1). The first version of this note identified the Jacobian of `C` with Sage and compared point counts, which was evidence and not a proof; the proof above replaces it. (`rd.py curve` still compares the numbers of points of `C` and `E` modulo 11 primes; the proof no longer uses it.)

**The general case, as a remark.** For classes `D₁, …, D₄` that are not equal the same computation gives a curve `C_D` whose Jacobian is the twist of `E` by `δ = D₁D₂D₃D₄`; for a point with four rational distances the twist is the congruent number curve `y² = x³ − d²x` of the squarefree part `d` of `δ` (`d = n` in Si's notation). We did not work this out further. Si's corollary quoted above does: it needs rank at least 2 over `ℚ`.

## 2. The first step, on a line

Fix `y = r/s` in lowest terms and write `x = u/v`. Up to the common factor `sv`, the four numbers are the linear forms

    su + i·rv,    s(u − v) + i·rv,    su + i·(r − s)v,    s(u − v) + i·(r − s)v

over `ℤ[i]`, whose classes are supported on the primes of `2rs` (first two) and of `2s(r − s)` (last two), by Lemmas S1 to S3 of `LAST_PAIR.md` applied to each form. This is a system of the kind that `scripts/cover_descent/tuple_descent.py` treats: it computes the tuples of four classes that pass every local test, at the primes of the supports and all primes up to 150. The tuple of the point at infinity of the line always remains.

**A line is closed** when every tuple that remains has four equal classes: by Proposition 1 there is then no point with four rational distances on it.

## 3. Results

`scripts/rational_distance/rd.py`, outputs in the same directory.

- **Lines `y = r/s` with `0 < r/s < 1/2` and `s ≤ 300`: 13,698 lines; 13,697 are closed** (`rd_scan_2_300.txt`, one line of output per line `y = r/s`). The exception is `y = 26/225`: the twelve tuples left are `(g, g, g, g)`, `(g, g, g, 199g)` and `(g, g, 199g, g)` for `g = 1, 3, 5, 15` (`225 − 26 = 199`, a prime `≡ 7 mod 8`).
- The midline `y = 1/2` is closed too (one tuple, four equal classes); this is a known case.
- By the symmetries of the square (`y ↦ 1 − y`, `x ↔ y`) the same holds for the lines `1/2 < y < 1` and for the vertical lines. **So, if the experiment is right, a point strictly inside the unit square at rational distance from the four corners has both coordinates with denominator above 300, or a coordinate equal to 26/225 or 199/225.**
- Outside the square: the 2,204 lines `y = r/s` with `1 < r/s ≤ 3` and `s ≤ 60` are all closed (`rd_scan_outside_60.txt`).

**Control with real points** (`rd.py control`, `rd_control.txt`). Points at rational distance from *three* corners exist; a search finds 14 with denominator below 400. With the third form in place of the four, the first step on the line through each of them keeps the tuple of the point, 14 times of 14. For none of the 14 are the three classes equal. So equal classes are not a feature of the method: it is the fourth distance that forces them, on all but one of the lines tried.

## 4. Against Si's results

- **The open line has the shape his theorem predicts.** On `y = 26/225` the tuples that are not constant are a common factor `g` (a divisor of the denominator) times `(1, 1, 1, 199)` or `(1, 1, 199, 1)`: pairwise coprime classes, a prime `≡ 7 (mod 8)` prime to the denominator, at a corner with `y`-coordinate 1 (199 is the numerator of `1 − y`). The kernel is `n = 199`. Our first step was written without knowledge of his theorem, so this is an independent agreement, on one case; on the other 13,697 lines the first step leaves `n = 1` only.
- **With his corollary the open line is closed, inside the square.** His Corollary "A uniform rank requirement" needs `rank E_199(ℚ) ≥ 2`. PARI's `ellrank` gives lower and upper bound 1 for `E_199: y² = x³ − 199²x` (`si_criterion_test.txt`, last line). So, **if his corollary and our first step are both right, no point strictly inside the square has a coordinate with denominator at most 300.** This rests on a statement of a preprint that is not formalised and that we have not checked line by line.
- **A test of his criterion** (`si_criterion_test.py`, `si_criterion_test.txt`). His corollary follows from a criterion for two points of `E_n` to be independent. On 72 curves `E_n` of rank exactly 1 every two points are dependent, so the hypotheses of the criterion must never hold together: of 17,280 ordered pairs of points tried, 608 satisfy its first hypothesis and none satisfies all three. This is a test that could have refuted the criterion and did not; it is not a proof of it.

## 5. What it rests on, and what could be wrong

- The first step is the program of `LAST_PAIR.md`, section 12, with its tests (section 14 there). **Exhaustion on this system** (`rd_exhaust.py`, `rd_exhaust.txt`): on six lines, `1/3, 2/5, 3/7, 1/8, 7/16` and the open line `26/225`, the exhaustion program (no code in common with the descent) keeps exactly as many tuples as the descent leaves, and the descent's tuples are among them: two on each of the first five lines, twelve on `26/225`. None of this has been read by anyone else.
- The support of the classes (primes of `2rs` and `2s(r − s)`) is Lemma S2 there, applied to these forms; the common-divisor step is three lines (a common divisor of `su + irv` and its conjugate divides `2su` and `2rv`, hence `2rs` when `u, v` are coprime) and is not in Lean for these forms. Si's theorem implies a sharper support inside the square (apart from the common factor, primes `≡ ±1 mod 8` of `r` and of `s − r`).
- Proposition 1: the proof above is complete as far as we can see, with the three quoted facts; it has not been read by anyone else.
- "Closed" is a statement about one line. It says nothing about lines with larger denominators. This is a census, not an approach to the problem itself.

## 6. Reproduce

```
cd scripts/rational_distance
python3 rd_prop1.py                                  # Proposition 1, exact arithmetic; ends with PROP1 OK
python3 rd.py line 1,3
python3 rd.py scan 2 300 > rd_scan_2_300.txt        # about 9 minutes
python3 rd.py scan 1 60 3 > rd_scan_outside_60.txt  # about 2 minutes
python3 rd.py curve
python3 rd.py control
python3 rd_exhaust.py 1,3 2,5 3,7 1,8 7,16 26,225 > rd_exhaust.txt      # about 5 minutes
python3 si_criterion_test.py 4 2>/dev/null | tail -2 > si_criterion_test.txt   # needs cypari2; about 1 minute
```

CI reruns `rd_prop1.py`, the lines with `s ≤ 40`, the curve check and the control on every push, and both scans, the exhaustion and the test of the criterion in the long jobs.
