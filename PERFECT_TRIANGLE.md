# The perfect triangle, line by line

**Frozen on 24 November 2026.** This note is a stress test of the descent programs on another problem. It is kept as it is and no further work is planned here: several reviewers pointed out that side problems add to what has to be checked without helping the main one, and they are right.

**Status.** An experiment of 24 November 2026: the descent of `LAST_PAIR.md`, sections 11 to 14, on a third open problem. Found by a reviewer instance of Claude (an AI model made by Anthropic) that was asked where the method transfers; checked and written up by another instance, directed by Drew (`drew2apps`). Nobody has read it. It is not part of the census of magic squares and nothing here is in Lean.

**The problem.** Is there a triangle with rational sides, rational medians and rational area (a Heron triangle with three rational medians, a "perfect triangle")? It is open (Guy, *Unsolved Problems in Number Theory*, D21). Triangles with rational sides, area and **two** rational medians exist (Buchholz and Rathbun).

**Prior work, as far as we looked.** A. Hone, *Heron triangles and the hunt for unicorns* (arXiv:2401.05581, 2024; read in part) calls the non-existence a longstanding conjecture and says that the proofs in the literature have turned out to be incorrect. M. Makhul (arXiv:1910.06888; abstract) proves finiteness for a fixed angle. A preprint of L. Shihaliev claims non-existence by an elementary identity and a descent (arXiv:2111.05118, first version 2021, thirteenth version July 2026; we read the abstract and a few lines, and have not checked it; the arXiv page gives no journal). S. Ismail, *Perfect triangles on the curve `C₄`* (J. Aust. Math. Soc., 2019; abstract), searches one of the eight curves of Buchholz and Rathbun that carry the Heron triangles with two rational medians. **This note proves nothing about the problem as a whole and does not depend on any of these.** We made about ten searches and did not find the formulation below or a census of this kind; that is weak evidence, and the reformulation with five collinear points may well be folklore. No priority is claimed.

## 1. Five norms

Scale the triangle so that `A = 0`, `B = 1`, and let `C = z = x + iy`, `y > 0`. Then the sides are `1, |z|, |z − 1|`, the medians are `|z + 1|/2`, `|z − 2|/2`, `|2z − 1|/2`, and the area is `y/2`. So the triangle is perfect exactly when `x, y` are rational and

    z,   z − 1,   z + 1,   z − 2,   2z − 1

have square norms: **`C` is at rational distance from the five collinear points `−1, 0, 1/2, 1, 2`.** As in `RATIONAL_DISTANCE.md`, each of the five numbers is then a positive squarefree integer, its class, times a square of `ℚ(i)`.

**Proposition T.** There is no `z ∈ ℚ(i)` outside `ℚ` for which the four numbers `z, z − 1, z + 1, z − 2` are in the same class of `ℚ(i)*` modulo squares.

*Proof.* If they are, their product is a square of `K = ℚ(i)`: `z` is the first coordinate of a point over `K` of `Y² = z(z − 1)(z + 1)(z − 2)`. With `X = 2/z` this is the elliptic curve `E: W² = (X − 1)(X − 2)(X + 2)` (conductor 24). Its rank over `K` is the sum of its rank over `ℚ` and that of its twist by `−1` (conductor 48); PARI's `ellrank` gives the bounds 0 and 0 for both (`pt_second.txt`, first line). Its torsion over `K`: `E` has 8 points over `F₅` and 16 over `F₁₃`, the residue fields of `K` at primes above 5 and 13, where it has good reduction, so as in Proposition 1 of `RATIONAL_DISTANCE.md` the torsion has order dividing 8, and the eight points with `X ∈ {∞, 0, 1, 2, −2, 4}` are all of it (`pt_tors.txt`). So `z = 2/X ∈ {0, ∞, 2, 1, −1, 1/2}`, all rational. ∎

Here the rank 0 is PARI's (two curves of conductor 24 and 48, which are in every table), not proved by hand.

## 2. The first step on a line, and the result

On the line `y = r/s`, with `x = u/v`, the five numbers are, up to the factor `sv`, linear forms in `(u, v)` over `ℤ[i]` whose classes are supported on the primes of `2rs`; `scripts/perfect_triangle/pt.py` gives them to `tuple_descent.py`. The point at infinity of the line always leaves the tuples `g·(1, 1, 1, 1, 2)`. `y = r/s` is twice the area divided by the square of the side taken as base, and with the longest side as base `y ≤ √3/2`.

- **Every line `y = r/s` with `0 < r/s ≤ 0.867` and `s ≤ 200`: 10,604 lines, and on every one of them the only tuples left are of the shape `g·(1, 1, 1, 1, 2)`** (`pt_scan_200.txt`). By Proposition T there is then no perfect triangle on the line. Unlike the unit square, where one line in 13,698 stayed open, no line is left here. The same run to `s ≤ 300`: 23,754 lines, no exception (about 16 minutes; the output is not committed, only its last line and its SHA-256, `pt_scan_300_summary.txt`; CI reruns the part to 200).
- **So, if the experiment is right: in a perfect triangle, for the longest side `c`, the number `2·Area/c²` in lowest terms has denominator above 300.** This is a census of lines, of the same kind as the census of ratios for the magic square, and not a step towards a proof: nothing says that the next line behaves in the same way.
- *Read as three separate things:* a theorem with a short proof that rests on two ranks from PARI (Proposition T); a computation, line by line, that is only as good as the program of `LAST_PAIR.md`, section 12; and the consequence above, which needs both. This is a restriction of a different kind from a bound on the sides: each line contains infinitely many rational triangles.
- **Control with real triangles** (`pt_control.py`, `pt_control.txt`): for five Heron triangles with two rational medians, placed so that the first four numbers have square norms, the first step with four forms keeps the tuple of the real point, five times of five; its four classes are never equal; and the fifth form removes it.

## 3. What it rests on, and what could be wrong

- The first step is the program of `LAST_PAIR.md`, section 12, with the tests listed there and in section 16. **No exhaustion has been run on the system of five forms**: the exhaustion programs and the two reviewer implementations handle four forms. The four-form subsystems are of the kind compared there, but the five-form run itself has no independent check.
- `pt_tors.txt` lists each value of `X` once for each point of the group that has it (a point and its negative have the same `X`), which is why values repeat; the six values of `z` are `0, ∞, 2, 1, −1, 1/2`, all rational, so none is a triangle.
- The numbers of points over `F₅` and `F₁₃` used for the torsion (8 and 16) are counted by plain enumeration in `pt_count.py` (`pt_count.txt`).
- The support of the classes: a common divisor of a form `s(u − av) + irv` and its conjugate divides `2s(u − av)` and `2rv`, hence `2rs` (for `2z − 1`, `4rs`), when `u, v` are coprime. Three lines, not in Lean.
- Proposition T rests on PARI for two ranks.
- The second step uses only the equality of four of the five classes; the other four quadruples of the five numbers give curves of rank at least 1 over `ℚ(i)` (`pt_second.txt`) and are not used.
- A uniform result to denominator 200 is not evidence of a theorem: on the unit square the first exception appeared at denominator 225.

## 4. Reproduce

```
cd scripts/perfect_triangle
python3 pt.py line 1,3
python3 pt.py scan 1 200 0.867 > pt_scan_200.txt      # about 15 minutes
python3 pt_second.py 2>/dev/null | grep -v "PARI stack" > pt_second.txt     # needs cypari2
python3 pt_tors.py 2>/dev/null | grep -v "PARI stack" > pt_tors.txt
python3 pt_control.py > pt_control.txt
python3 pt_count.py > pt_count.txt
```
