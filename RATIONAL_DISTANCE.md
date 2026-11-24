# The rational distance problem for the unit square, line by line

**Status.** An experiment, started on 24 November 2026 with the descent of `LAST_PAIR.md`, sections 11 to 14. Written by Claude (an AI model made by Anthropic), directed by Drew (`drew2apps`). Nobody has read it. It is not part of the census of magic squares and nothing here is in Lean.

**The problem.** Is there a point of the plane at rational distance from the four corners of the unit square? It is open (Guy, *Unsolved Problems in Number Theory*, D19). Known: there is no such point on the edges, the midlines or the diagonals of the square, and points at rational distance from three corners are plentiful (Berry, 1990). We have not read Berry's paper and have not searched the literature for results line by line; what follows may be known.

## 1. The four classes

A point with four rational distances has rational coordinates (the differences of the squared distances give `x` and `y`). Put `z = x + iy`. The four distances are rational exactly when

    z,   z − 1,   z − i,   z − 1 − i

have square norms. A nonzero element of `ℚ(i)` of square norm is a positive squarefree integer, its *class*, times a square of `ℚ(i)`, and the class is unique (`LAST_PAIR.md`, Lemma S1). So a solution has four classes `D₁, D₂, D₃, D₄`, up to a common rational factor (the denominator of `z`).

**Proposition 1.** There is no `z ∈ ℚ(i)`, other than the corners, for which the four numbers `z, z − 1, z − i, z − 1 − i` are in the same class of `ℚ(i)*` modulo squares. In particular, at a point with four rational distances the four classes are not all equal.

*Proof.* Suppose `z − a_j = D·w_j²` with one `D ∈ ℚ(i)*` and `w_j ∈ ℚ(i)`. Put `W_j = w_j²`. Then `W₁ − W₂ = 1/D`, `W₁ − W₃ = i/D`, `W₁ − W₄ = (1 + i)/D`, so

    C:   w₃² = (1 − i)·w₁² + i·w₂²,     w₄² = −i·w₁² + (1 + i)·w₂².

`C` is a smooth intersection of two quadrics in `P³`, a curve of genus 1 over `K = ℚ(i)`, with the eight points `(1 : ±1 : ±1 : ±1)`. Its Jacobian is the Jacobian of `Y² = det(x·Q₁ + Q₂) = x((1 − i)x − i)(ix + 1 + i)`, which has `j = 1728` and is isomorphic over `K` to `E: y² = x³ − x` (computed with Sage; `rd.py curve` compares the numbers of points of `C` and of `E` modulo 11 primes, in both embeddings: equal every time). `E` has rank 0 over `ℚ`, and its twist by `−1` is `E` again, so `E` has rank 0 over `K`; its torsion over `K` has 8 points (it injects into `E(F₅)`, which has 8). `C` has a point over `K`, so it is isomorphic to `E`, and `C(K)` is the eight points above. At each of them `w₁² = w₂²`, which contradicts `W₁ − W₂ = 1/D`. ∎

The rank 0 of `y² = x³ − x` over `ℚ(i)` is classical; it is the statement that `x⁴ − y⁴ = z²` has only trivial solutions in `ℤ[i]`.

**The general case, as a remark.** For classes `D₁, …, D₄` that are not equal the same computation gives a curve `C_D` whose Jacobian is the twist of `E` by `δ = D₁D₂D₃D₄`. For a point with four rational distances `δ` is rational, and the twist is `y² = x³ − d²x`, the curve of the congruent number problem for the squarefree part `d` of `δ`. Its rank over `ℚ(i)` is twice its rank over `ℚ`. So a point with four rational distances gives a point over `ℚ(i)` on a torsor of the congruent number curve of `d`, and if `d` is not a congruent number that torsor has finitely many points. We have not worked out what those finitely many points give; only the case `d = 1` with equal classes (Proposition 1) is used below.

## 2. The first step, on a line

Fix `y = r/s` in lowest terms and write `x = u/v`. Up to the common factor `sv`, the four numbers are the linear forms

    su + i·rv,    s(u − v) + i·rv,    su + i·(r − s)v,    s(u − v) + i·(r − s)v

over `ℤ[i]`, whose classes are supported on the primes of `2rs` (first two) and of `2s(r − s)` (last two), by Lemmas S1 to S3 of `LAST_PAIR.md` applied to each form. This is a system of the kind that `scripts/cover_descent/tuple_descent.py` treats: it computes the tuples of four classes that pass every local test, at the primes of the supports and all primes up to 150. The tuple of the point at infinity of the line always remains.

**A line is closed** when every tuple that remains has four equal classes: by Proposition 1 there is then no point with four rational distances on it.

## 3. Results

`scripts/rational_distance/rd.py`, outputs in the same directory.

- **Lines `y = r/s` with `0 < r/s < 1/2` and `s ≤ 300`: 13,698 lines; 13,697 are closed** (`rd_scan_2_300.txt`, one line of output per line `y = r/s`). The exception is `y = 26/225`: tuples with one class multiplied by 199 remain (`225 − 26 = 199`, a prime `≡ 7 mod 8`), and the congruent number curve of 199 has rank 1 (PARI), so the remark of section 1 gives nothing there.
- The midline `y = 1/2` is closed too (one tuple, four equal classes); this is a known case.
- By the symmetries of the square (`y ↦ 1 − y`, `x ↔ y`) the same holds for the lines `1/2 < y < 1` and for the vertical lines. **So, if the experiment is right, a point strictly inside the unit square at rational distance from the four corners has both coordinates with denominator above 300, or a coordinate equal to 26/225 or 199/225.**
- Outside the square: the 2,204 lines `y = r/s` with `1 < r/s ≤ 3` and `s ≤ 60` are all closed (`rd_scan_outside_60.txt`).

**Control with real points** (`rd.py control`, `rd_control.txt`). Points at rational distance from *three* corners exist; a search finds 14 with denominator below 400. With the third form in place of the four, the first step on the line through each of them keeps the tuple of the point, 14 times of 14. For none of the 14 are the three classes equal. So equal classes are not a feature of the method: it is the fourth distance that forces them, on all but one of the lines tried.

## 4. What it rests on, and what could be wrong

- The first step is the program of `LAST_PAIR.md`, section 12, with its tests (section 14 there): random local points, points near the zeros of the forms, the class rules against lists of squares, and an exhaustion program with no code in common that agrees on the cases compared. None of this has been read by anyone else, and the exhaustion has not been run on this system.
- The support of the classes (primes of `2rs` and `2s(r − s)`) is Lemma S2 there, applied to these forms; the common-divisor step is three lines (a common divisor of `su + irv` and its conjugate divides `2su` and `2rv`, hence `2rs` when `u, v` are coprime) and is not in Lean for these forms.
- Proposition 1 uses the formula for the Jacobian of an intersection of two quadrics and Sage's identification of the curve; the comparison of point counts is evidence, not a proof. The rest of its proof is classical.
- "Closed" is a statement about one line. It says nothing about lines with larger denominators, and one line with `s ≤ 300` is already open. This is a census, not an approach to the problem itself.

## 5. Reproduce

```
cd scripts/rational_distance
python3 rd.py line 1,3
python3 rd.py scan 2 300 > rd_scan_2_300.txt        # about 9 minutes
python3 rd.py scan 1 60 3 > rd_scan_outside_60.txt  # about 2 minutes
python3 rd.py curve
python3 rd.py control
```

CI reruns the lines with `s ≤ 40`, the curve check and the control, and checks the counts of the committed files.
