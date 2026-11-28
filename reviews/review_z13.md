# Review z13: PAIRING_FORMULA_2.md (second note on the Cassels–Tate pairing on E_08, E_LU)

Written by an AI model (reviewer z13), not read by a person. 27 November 2026 (UTC).
Scratch: /home/claude/conductor/z13/experiments/referee_z13/ (branch run11_z13).

> Committed as written apart from this note. In this repository the scratch folder is `experiments/referee_z13/`; the note reviewed is `experiments/research-pairing/note2/` (revision 2, and revision 3 in the section "Revision 3").

(Findings appended as they are made; bottom line at the end once known.)

## Part 1: the proofs, line by line (written first; checked against computation in Part 2)

Notation as in the note: E: y^2 = x(x - A^2)(x - C^2), (A, B, C) primitive, C odd, A, C > 0; a = (1, d, d) in Sel^2,
b = (b_r, b_m, b_*) in Sel^2; the covering u^2 = x, d u_m^2 = u^2 - A^2, d u_*^2 = u^2 - C^2.
Facts used throughout: gcd(A, B) = gcd(B, C) = gcd(A, C) = 1, so an odd p divides at most one of A, B, C, and an odd
p | B divides exactly one of C - A, C + A (their gcd divides 2); 2 | B iff A is odd; A even => 4 | A
(A^2 = C^2 - B^2 = 0 mod 8). Every Pythagorean triple has 3 | AB and 5 | ABC (squares mod 5 are 0, +-1).

- **Lemma 1 (p | B odd; p = 2, A odd, H2): holds.** (x - C^2)/(x - A^2) = 1 - B^2/(x - A^2). Branch
  v(x - A^2) <= v(B^2) - 1 (odd p), resp. <= v(B^2) - 3 (p = 2): the ratio is 1 mod p, resp. 1 mod 8, a square, so
  b_m = b_* and b_r = square. Otherwise v(x - A^2) >= v(B^2) >= 2 (odd p), so x = A^2 (mod p), a unit square; for
  p = 2, v(x - A^2) >= v(B^2) - 2 >= 4 because H2 gives v2(B^2) >= 5 + 1, so x = A^2 = 1 (mod 8). Thresholds: the two
  branches are complementary at the boundary values (v = v(B^2), resp. v(B^2) - 2 fall in the second branch and are
  handled there). v(x) < 0 falls in the first branch. The 2-torsion images (1,-1,-1), (1,-1,-1), (1,1,1) are right
  (T_r: (A^2C^2, -A^2, -C^2); T_m: (A^2, -A^2B^2, -B^2); T_*: (C^2, B^2, C^2B^2)).
- **Lemma 2: holds.** Odd p | A: v(x) < 0 gives all three = x times squares, so x and x - C^2 are squares; x a unit,
  x != 0, C^2 mod p: x - A^2 = x(1 - A^2/x) = x mod squares, so x - C^2 is a (unit) square; x = C^2 mod p: x and
  x - A^2 unit squares, so x - C^2 is a square; v(x) >= 1: x - C^2 = -C^2 mod p. So b_* is a unit, in {1, -1} mod
  squares. p | C: the same with A, C exchanged. p = 2, A even (so 4 | A, B and C odd): I redid the five cases
  (v(u) <= -2: both classes 1; v(u) = -1: classes 1 and 5; u odd: u^2 - A^2 = 1 mod 8; v(u) = 1: u^2 - C^2 = 3 mod 8 and
  u^2 - A^2 in 4{1,5}; v(u) >= 2: u^2 - C^2 = -1 mod 8); the points over 2-torsion (u = +-A: class -1; u = +-C: class 1)
  fit. All as stated.
- **(i) d a square in Q_v: holds (given Cassels' theorem).** For |u|_v large, u - A, u - C are u times squares (for
  v = 2 need |A/u|_2, |C/u|_2 <= 2^-3), u_m - u_* = (B^2/(2su))(1 + eps) with 1 + eps a square, and
  (u, b_m b_*) (2su, b_r) = (u, b_r)(2su, b_r) = (2s, b_r). The other sign of u_* gives 2u/s instead of B^2/(2su),
  the same class mod squares times s^2, so the same value. Real place: d > 0, s > 0, (2s, b_r)_inf = 1.
- **(ii) -d a square in Q_v: holds (given Cassels' theorem).** u = 0 is itself a point of the covering (over T_r):
  u_m = A/r, u_* = C/r, and the three entries are -A, -C, (A - C)/r = (A - C) r mod squares. The other branch
  (u_* = -C/r) changes the term by ((A - C)(A + C), b_r) = (-1, b_r)_v; by Cassels' theorem the two must agree, so this
  is a consequence, not a gap. Real place: d < 0; with A > 0 all three entries are negative and the term is
  (-1, b_r b_m b_*)_inf = 1. Needs A > 0 and C > A (true: C is the hypotenuse); both are in the note's setting.
- **(iii) p | B, C = A: holds at every point of the covering where the entries are non-zero (no Cassels' theorem).**
  Lemma 1 for b gives t_p = ((u - A)(u - C), beta). Branch v(u - A) <= v(C - A) - 1 (resp. - 3): ratio 1 mod p (resp.
  mod 8). Otherwise v(u - A) >= v(C - A) >= 1 (odd p), resp. >= v(C - A) - 2 >= 3 (p = 2), so v(u + A) = 0, resp. 1,
  and (u + C)/(u + A) = 1 + (C - A)/(u + A) = 1 mod p, resp. mod 2^(v(C-A)-1) with v(C - A) - 1 >= 4. With
  (u^2 - A^2)(u^2 - C^2) = d^2 u_m^2 u_*^2 a square, (u - A)(u - C) is a square. At p = 2 the argument needs
  v2(C - A) >= 4, i.e. >= 5: exactly H2.
- **(iv) p | B, C = -A: holds (no Cassels' theorem).** The note's "small valuation" means v(u + A) <= v(C + A) - 1
  (odd p), resp. - 3 (p = 2); then u - C = (u + A)(1 - (C + A)/(u + A)) and (u - A)(u - C) = (u^2 - A^2) = d mod
  squares. Otherwise v(u + A) >= v(C + A) (resp. - 2 >= 3), v(u - A) = 0 (resp. 1), u + C = (u - A)(1 + (A + C)/(u - A))
  with the bracket 1 mod p (resp. mod 16), and (u - A)(u - C) = (u^2 - C^2) = d mod squares. So t_p = (d, b_*)_p.
  [fix] Write the threshold out ("small" is undefined in the text).
- **"An odd p not dividing ABC gives 1": holds, but the justification in the note is incomplete [fix].** d, b are
  units there (Sel^2 is supported on -1, 2 and the primes of ABC). If p = 3 mod 4, d or -d is a square and (i)/(ii)
  give Hilbert symbols of units: 1. If p = 1 mod 4 and d is a non-square unit, neither (i) nor (ii) applies and the
  note's "all entries are units for a suitable point" needs an argument: every Q_p-point has u integral (a
  non-integral u forces d = square), the reduction of the covering is a smooth genus-one curve with
  #D(F_p) >= p + 1 - 2 sqrt(p) affine points, at most 4 of which have u = A or u = C, and u_m - u_* is then
  automatically a unit (d(u_m^2 - u_*^2) = B^2). This leaves p >= 13 fine and p = 5 open, but 5 | ABC for every
  Pythagorean triple, so p = 5 is never a good prime. (Or quote the standard good-reduction lemma of Cassels 1998.)
- **The real place: holds** (see (i), (ii); needs A > 0, stated).
- **[fix] The representative of d.** (2 sqrt(d), b_r)_v depends on the representative of d mod squares: d -> k^2 d
  multiplies T_v by (k, b_r)_v, and the good-prime terms by (k, b_r)_p, which is not 1 when p | k is a good prime
  with b_r a non-residue. The formula is right for d squarefree (or square part supported on the primes of ABC); the
  scripts use squarefree d, the text does not say so. (The sign of sqrt(d) does not matter, by Cassels' theorem.)
- **H2 and its derivation (b): holds.** E_08, m, n of opposite parity: C - A = 2n^2, C + A = 2m^2; 32 | 2n^2 iff 4 | n;
  so H2 iff 4 | mn; m = +-n mod 8 is impossible. E_08 both odd: A = (m^2 - n^2)/2 is even, H2 is empty. E_LU opposite
  parity: A = 2mn even, empty. E_LU both odd: A = mn, C -+ A = (m -+ n)^2/2, and 32 | (m -+ n)^2/2 iff 8 | m -+ n. So H2
  is equivalent to (4 | mn or m = +-n mod 8) on both curves. v2(C -+ A) is odd (v2(B^2) is even and the other
  difference has v2 = 1), so "5 or at least 7" is automatic.
- **(c) S' = {b in Sel^2 : b_r = 1}: holds (given H2 at 2 when A is odd).** "Subset": b_r = 1 gives b_m = b_* = beta;
  Lemma 2 at odd p | A, p | C and at 2 (A even) makes beta a unit there and +- a square; so beta is supported on -1
  and the primes of B. "Superset": at p | B (and 2 | B) the local image has |E(Q_p)/2E(Q_p)| = 4 (resp. 8) elements,
  Lemma 1 puts all of them among the 4 (resp. 8) distinct classes (1, beta, beta), so it is all of them; at v | AC
  (and 2 | A) beta = +-1 gives (1,1,1) = image of O or (1,-1,-1) = image of T_r; at good odd p beta is a unit,
  i.e. unramified, and the local image at a good odd prime is the unramified subgroup; at R the image of E(R) is
  {(1,1,1), (1,-1,-1)} (x >= C^2, resp. 0 <= x <= A^2). Every step holds.
- **"S' is the image of the Selmer group of the 2-isogeny E -> E/<(0,0)>": holds** (my check, not in the note):
  for phi with kernel <T_r>, iota: H^1(E[phi]) -> H^1(E[2]) is beta -> (1, beta, beta), injective because E[2] is
  rational, and iota^{-1}(local image of [2]) = local image of phi by the diagram E'/phi E -> E/2E -> E/phi-hat E'. So
  {b_r = 1} part of Sel^2 = iota(Sel^phi(E)).
- **Theta (3) from (2) with b_r = 1: holds** (the T_v with d a square is (2s, 1) = 1; the other is (-A, beta)(-C, beta)
  = (AC, beta)). It is bilinear: the factor at p = 3 mod 4 is included iff (d/p) = -1, a character of d; at p = 1 mod 4
  "not a square but -d a square" never happens; at 2 (A even) d is in {1, -1} mod squares and the indicator of -1 is
  a character.

**What has to be trusted (not proved in the note):** (T1) Cassels 1998's formula: for a, b in Sel^2 of a curve with
full 2-torsion, the Cassels-Tate pairing is prod_v prod_i (l_i(P_v), b_i)_v with l_i the tangent at a rational point of
the conic that omits u_i, and the component of b that goes with each l_i (z12 derived this assignment from the Kummer
map and found it is the one rank_second uses; neither z12 nor I compared it with Cassels' or FSS's own indexing);
(T2) Cassels' local constancy (t_v does not depend on P_v in D_a(Q_v) when b_v is in the local image); (T3) FSS 2010
for "Cassels' pairing = Cassels-Tate pairing"; (T4) standard facts: |E(Q_p)/2E(Q_p)| = 4 (odd p), 8 (p = 2), 2 (R);
the local image at a good odd prime is the unramified subgroup. Given these, every step of section 2 and of the
"superset" direction holds; the gaps are of exposition (the good-prime argument, the threshold in (iv), squarefree d).
No step is wrong.

## Part 3 (first half): the note's numbers, recomputed

Reruns of the note's scripts (copies in a scratch folder outside the repository, with `ct08.py`, `parse.py` of the
first note, z12's `full0/1.jsonl` and the pilot's `t10_*_[ab].txt`; /tmp/t10 does not exist, so the hard-coded path is
inert). Times on this machine: closed_formula 33 s, isogeny_part 11 s, shape 2 s, two 1 s, two2 2 s.

- 69,304 (E_08, agree), 35,096 (E_LU, agree), 316 (E_LU rows skipped), 1,780 + 784 rows without the table:
  reproduce exactly (`closed_formula.txt` identical). [fix] "agree with the committed values" means: with ct08's
  NUMERICAL values stored by z12 (`ct_all`), which z12 found equal to the linear extension of rank_second; the note's
  test is therefore against rank_second only through z12's check. Also, the 316 skipped rows are every tuple class
  on E_LU, so on E_LU the test covers only the (at most 8 per curve) non-tuple (1, d, d) rows.
- 720 and 1,434 (two.py), 1,434 = 320 + 184 + 472 + 458 (two2.py): reproduce (identical output). [fix] These are
  pilot kill lines, i.e. repetitions: 720 lines = 24 (E_08) + 66 (E_LU) = 90 distinct (ratio, curve, class, b); 1,434
  lines = 109 + 45 = 154 distinct (`claims_z13.txt`). The note should give distinct counts (z12's correction 2,
  repeated).
- 146 + 133 curves with S' = {b_r = 1}, 29,484 + 11,928 Theta pairs, dimensions 1..7 (E_08) and 1..6 (E_LU):
  reproduce (identical). Same caveat: "the pairing" is ct08's stored value.
- 172 of 202, 86 of 134 (shape.py): reproduce. [fix] `shape.txt` is not the script's full output: the rerun prints a
  fourth example line, `(1323, 803, '08', [1, -2409], ...)`, missing from the committed file. "b_r = +-1" in the
  text: b_r = -1 cannot occur in Sel^2 (the image of E(R) has b_r > 0), so this is b_r = 1, i.e. Theta. Correct.
- 232 of 1,762 (section 6): reproduces (1,762 tuples killed on some pair curve, 232 on mixed curves only).
- (5623,1017): the two surviving tuples have d0 d8 = 2,541,596 = 4 * 635,399 and 635,399; [note] "the two tuples
  have d0 d8 = 4 * 635399" is true only up to squares for the second. The squarefree part of mn = 5623 * 3^2 * 113 is
  635,399, so d = 1 and the class is trivial: right. [fix] "This is the shape of the pilot's exceptions" describes
  one of the two exceptional ratios only: at (3320,2303) the same tuples are "kernel" (not "image") on E_08, with
  d = sqf(mn * 2) = 19,505, not a torsion class; there it is E_LU that shows the torsion class. ((3320,2303) and
  (5623,1017) are related by (m, n) -> (m + n, m - n), which exchanges E_08 and E_LU, so the two ratios are one
  example seen twice.)

## Part 2: the independent test, outside the pilot

Driver `experiments/referee_z13/z13_indep.py` (job list `gen_jobs.py`, seed 20261127; my own implementation of (2),
S', Theta in `z13_formula.py`, written from the note's text, the note's scripts not imported). Curve built directly
as y^2 = x(x - A^2)(x - C^2), shifted as PairData does; Sel^2 by `rank_second.full_selmer2` (seed 1); the pairing by
`t10_pairct.DFullX` for EVERY a = (1, d, d) in Sel^2 (all a evaluated, none extrapolated) against a basis of Sel^2,
extended to all of Sel^2 by multiplicativity in b; formula (2) on every (a, b), b over all of Sel^2.

- **378 curves of 197 ratios (181 with both curves, 16 controls with one), none in the pilot**: random m <= 20000 (80 curves) and m <= 2000 (60), both odd (30),
  A odd with v2(C -+ A) = 5 exactly (v5_08: 24 curves, v5_LU: 24), square type m^2 + n^2 = square (28), twice-square
  m^2 + n^2 = 2 square (28), every ratio m <= 18 with a 2-adic point (88), and 16 controls OUTSIDE H2 (no 2-adic point).
- **Within H2: 362 curves, 336 distinct curves (A, C)** (the map (m, n) -> (m + n, m - n) exchanges E_08 and E_LU, so
  some curves occur twice). Sel^2 of dimension 2..8, S' of dimension 1..7; A odd on 168 distinct curves, of them
  v2(C -+ A) = 5 exactly on 86; 3 | B on 168.
  - **formula (2) = rank_second on 194,504 of 194,504 pairs (a, b)** (193,576 on distinct curves), 16,296 of them with b
    a basis element (rank_second's own values, no extension). 0 disagreements, 0 undefined.
  - **S' = {b in Sel^2 : b_r = 1} on 362 of 362 curves.**
  - **Theta = the pairing on S' x S' on 96,632 of 96,632 pairs** (96,216 on distinct curves).
  - rank_second is bilinear in a on all 11,585 checks (values on non-basis a against products over a basis).
- **Controls outside H2 (16 curves, A odd, v2(C - A), v2(C + A) = {1, 3})**: formula (2) is WRONG on exactly half of
  the pairs on every curve (3,188 of 6,376); smallest: E_LU at (725, 497), A, B, C = 360325, 139308, 386317,
  a = (1, 1, 1) (trivial), b = (5365, 435): rank_second 1, formula (2) -1 (b_r = 5365 = 5 mod 8 is not a 2-adic square,
  so Lemma 1 fails at 2). S' = {b_r = 1} and Theta still agreed on all 16. So H2 is a real hypothesis of (2), not an
  artefact of the proof; this is not a counterexample to the note, which assumes H2.
- Side checks (`analyse_z13.py`): (e1) where both branches (i), (ii) apply (p = 1 mod 4 dividing AC) they agree on
  388,216 of 388,216 evaluations, as Cassels' theorem requires; (e2) (-1, b_r)_p = 1 at every p | AC (and 2 | A) where
  +-d is a square, on 862,880 evaluations, so the sign of sqrt(d) is immaterial; (e3) on E_08 at all 14 square-type
  ratios, (1, 2mn, 2mn) is in S' and in the radical of Theta, and 2mn is not +-square; (e4) the representative of d
  matters: at E_08 (535, 23), d = 1, b = (41, 1), formula (2) with d replaced by 7^2 gives -1, with d squarefree 1
  (rank_second: 1); (e5) on E_08 at the other ratios, (1, 2mn, 2mn) is in S' on 32 of 93 opposite-parity curves and in
  the radical on 16 of those 32 (never in S' for m, n both odd, where 2 does not divide B).
- Reproducibility: `rerun_z13.py` reruns 25 curves (every 25th job and the first of each category) from scratch:
  identical records (89 s). `analyse_z13.py` 30 s; `claims_z13.py`, `gen_jobs.py` 1 s each. The full run: 1,778 s of
  per-curve wall time on two shared cores.

## Part 3 (second half): status words and citations

- **Section 5, "candidate answer".** The displayed statement ("a tuple that survives the pairing ... has a torsion
  class there, or the radical of Theta is larger than +-1") follows in two lines from (c) and (3): surviving the pairing
  implies surviving Theta, i.e. d in the radical; if the radical is {+-1}, the class is (1, +-1, +-1). So "Not
  proved" is too weak for the implication itself (it is as proved as (2), given Cassels), and "the candidate answer to
  the question" is too strong for its content: it is a restatement, and by section 7 the second disjunct holds on
  about two thirds of the non-degenerate ratios per curve, where it decides nothing. [fix] Say "consequence of (c)
  and (3)", drop "answer". [note] "which only a rank statement removes": a torsion class is the class of a rational
  point, so rank 0 alone does not remove it; one needs E(Q) and the fibre over its points.
- **Section 7, "The reason is exact".** Exact for the square-type statement, with one missing step: (1, 2mn, 2mn) is
  the class of the point of order 4 (x = C^2 + BC gives (C(C + B), B(B + C), BC) = (S, 2mn, 2mnS) on E_08), it pairs
  trivially with everything, and it is not +-1 because 2mn = square with m^2 + n^2 = square would give
  a^4 + 4b^4 = square (Fermat); the note does not say why 2mn is not +-square [fix]. The non-degenerate part is an
  explanation, not a proof, and it skips a step: before formula (2) applies to (1, 2mn, 2mn) it must be in S' (a
  Legendre condition at the primes of m^2 +- n^2; true on 32 of 93 of my opposite-parity E_08 curves), and only then
  is "killed or not" the quartic-symbol question (radical on 16 of those 32). "(2 sqrt(b), S)_p ... the quartic
  residue symbol of b there" is (2/p) times it [note]. "It is the class of the degenerate solutions at t = +-1 (there
  d0 d8 = 2)" I did not check. "These are counts on a sample, not a density statement": right.
  The 400 "square type" ratios of `radical.py` are the first 400 of a fixed enumeration (a < 60), not random [note].
- **Radical table**: reproduces exactly (`radical.py`, 10 s); 482/1500 and 494/1500, "about a third": right.
- **Fisher 2017**: T. Fisher, Higher descents on an elliptic curve with a rational 2-torsion point, Math. Comp. 86
  (2017) 2493-2518 (arXiv:1509.03234). Its Theorem 2.1 defines, for phi: E -> E' of degree 2, alternating pairings
  Theta_1 on S^(phi)(E/K) x S^(phi)(E/K) from the Cassels-Tate pairing, with kernel S_2 (the next descent). The note's
  "Theta is the pairing on it whose kernel is the next step of a descent by 2-isogenies (Fisher 2017)" is a correct use,
  given that S' = iota(S^(phi)(E/Q)) for phi with kernel <(0,0)> (I checked that identification; injectivity of iota
  holds because E[2] is rational). [fix] Give the full reference (it is not in the first note's list).
- **"Monsky-style"**: Monsky's appendix to Heath-Brown, Invent. Math. 118 (1994) 331-370, computes the 2-Selmer rank
  of the congruent number curve as the corank of a matrix of Legendre symbols and n mod 8. The analogy is fair for
  S' (a kernel of Legendre conditions; it is a 2-isogeny Selmer group, not all of Sel^2 as in Monsky) and for Theta
  (a Gram matrix). [note] "matrices of Legendre symbols" omits the 2-adic entries (b = +-1 mod 8 when A is even; the
  2-adic Hilbert symbol in Theta when 2 | B and C = -A mod 8).
- **Fisher 2022** (section 6): Res. Number Theory 8 (2022) 74: three locally soluble quartics with the same
  invariants and z(g1) z(g2) z(g3) a square in L = K[phi]; Remark 3.2(i): finding g3 from g1, g2 needs one conic. The
  note's description is right.
- **Header**: "Tier D; no count changes", "not reviewed": fine. Title "every local term derived": right given
  Cassels' theorem.

## Part 4: tier

The formula now merits the written-proof tier (C) as an experiment, conditional on: (T1) Cassels 1998's formula in
the form (1), including which component of b goes with which tangent (this is the one input pinned only by z12's own
derivation and by agreement with rank_second; for tier C it should be matched to a numbered statement of Cassels 1998
or FSS 2010 / Shukla-Stoll); (T2) Cassels' local constancy for b in the local image; (T3) FSS 2010 (Cassels' pairing
= Cassels-Tate); (T4) the standard local image sizes (4, 8, 2) and the unramified image at good odd primes. Given
these, every step of section 2 and of (c) holds; the fixes are expository: (F1) the good-prime argument
(p = 3 mod 4 by (i)/(ii); p = 1 mod 4 by the Hasse bound and 5 | ABC); (F2) "d squarefree" (counterexample (e4) for
other representatives); (F3) the threshold in (iv); (F4) distinct counts in section 4 and "against ct08's stored
values". Sections 5 and 7 stay tier D (sections 5's implication is a corollary; 7 is counts plus an explanation).

## Bottom line

- **Proofs: hold.** Every step of section 2 (Lemmas 1, 2, cases (i)-(iv), the real place, the good primes), H2 and its
  derivation, and both directions of S' = {b in Sel^2 : b_r = 1} hold, given Cassels' theorem in the form (1) and
  standard local facts. No step is wrong. Three expository gaps: the good-prime claim has no argument for p = 1 mod 4
  (one exists: Hasse bound plus 5 | ABC); d must be taken squarefree (otherwise (2) is false, example (e4)); the
  threshold of (iv) is not written.
- **Independent test: 0 disagreements.** On 362 curves (336 distinct) of 181 ratios outside the pilot, within H2:
  formula (2) = rank_second on 194,504 of 194,504 pairs (a over all of {a_r = 1}, b over all of Sel^2), S' = {b_r = 1}
  on 362 of 362, Theta = the pairing on 96,632 of 96,632. Outside H2 (16 control curves) formula (2) fails on exactly
  half the pairs (smallest: E_LU (725, 497), a trivial, b = (5365, 435)); S' and Theta still agree there.
- **Numbers:** all reproduce exactly; several are pilot kill lines (repetitions) and the "agreement" of section 4 is
  with ct08's stored values, not rank_second directly; one example line is missing from `shape.txt`.

## Findings

1. [fix] Good primes: "all entries are units for a suitable point" needs an argument when p = 1 mod 4 and d is a
   non-square unit; Hasse bound (p >= 13) plus 5 | ABC for every Pythagorean triple closes it.
2. [fix] State d squarefree in (2) (T_v depends on the representative; at E_08 (535, 23) with d = 1 written as 7^2 the
   formula gives -1 against the true 1).
3. [fix] (iv): write the threshold (v(u + A) <= v(C + A) - 1, resp. - 3).
4. [note] H2 is necessary for (2) (controls); S' and Theta survived its failure on 16 of 16 controls (unexplained,
   not claimed by the note).
5. [fix] Section 4: give distinct counts (720 lines = 90 distinct, 1,434 = 154) and say the comparison is with ct08's
   stored values (equal to rank_second by z12); on E_LU the tuple rows were all skipped.
6. [fix] `shape.txt` lacks the fourth example line the script prints.
7. [fix] Section 5: the displayed statement is a corollary of (c) and (3), not a "candidate answer"; it says nothing on
   the two thirds of ratios where the radical is larger; "only a rank statement removes" is inexact.
8. [fix] Section 5: "the pilot's exceptions" fits (5623,1017) only; at (3320,2303) the E_08 class is in the kernel and
   not torsion (d = 19,505); the two ratios are related by (m, n) -> (m + n, m - n). "d0 d8 = 4 * 635399" holds up to
   squares only for the second tuple.
9. [fix] Section 7: add why 2mn is not +-square at square type (a^4 + 4b^4 = square has no solution); the non-degenerate
   explanation omits the S'-membership (Legendre) step; "(2 sqrt(b), S)_p" is (2/p) times the quartic symbol.
10. [fix] Give the full reference for Fisher 2017 (Math. Comp. 86, 2493-2518; Theorem 2.1, Theta_1 on S^(phi) x S^(phi));
    the use is correct. [note] "matrices of Legendre symbols" omits the 2-adic entries.
11. [note] Tier: formula (2) and (c) can go to tier C conditional on (T1)-(T4); (T1)'s component assignment should be
    matched to a numbered statement of Cassels 1998 or FSS 2010 before that.

## Revision 3 (PAIRING_FORMULA_3_2026-11-27.md, h2_check.py), follow-up review, 27 November 2026

Written by an AI model (reviewer z13), not read by a person.

**Bottom line.** Revision items (1), (2), (3), (4), (6) and (8) are applied correctly. (7) is applied correctly as far as it goes.
(5) has the right references and the right conclusion, but one sentence is false as written. Two problems are new.
First, the Prior art paragraph's "correction" of the arXiv number is wrong: the first note's 2303.05059 was right.
Second, h2_check.py cannot be run as supplied. Not applied from my first review: finding 5 (the §4 counts and
wording), finding 6 (`shape.txt`), and the S′-membership step of finding 9. The header still says "not reviewed".

Item by item, numbered as in the revision:

- **(1) Good primes: applied, correct.** I checked each step. Classes of Sel² and squarefree δ are units at p ∤ 2ABC.
  If δ is a square, (i) gives a symbol of two units; if −δ is, (ii) does, with A − C a unit because
  (A − C)(A + C) = −B². The remaining case is p ≡ 1 mod 4 with δ a non-square unit.
  - The curve u² − A²T² = δu_m², u² − C²T² = δu_*² is smooth: the pencil's singular members are at four distinct
    values because A² ≢ C² and A, C ≢ 0 mod p. So it has genus 1, and Hasse–Weil gives p + 1 − 2√p > 0 points.
  - At every F_p-point the four exclusions hold. T = 0 would make δ a square. u ≡ ±A gives u_m ≡ 0 and −δ a square.
    u ≡ ±C gives δ a square. u_m ≡ u_* contradicts δ(u_m² − u_*²) = B², a unit.
  - The point lifts by Hensel, and all entries are units. Because every F_p-point is good, no lower bound on p is
    needed. This is better than my version (which subtracted 4 points and needed 5 | ABC).
- **(2) Squarefree: applied, correct.** The (k, b_r)_v remark and my example (δ = 7² at E_08 (535, 23)) are both right.
- **(3) Threshold in (iv): applied, correct.** Case 1 is v(u + A) ≤ v(C + A) − 1, resp. − 3. Case 2: u ≡ −A, so
  v(u − A) = 0, resp. 1, and (u + C)/(u − A) ≡ 1 mod p, resp. mod 16. Case 2 now argues directly through (u − A)(u + C),
  which is fine. "C ≡ −A mod 8, and then H2 gives 32 | C + A" is right, because the other difference is ≡ 2 mod 4.
- **(4) H2 necessary: applied, numbers correct.** I could not rerun h2_check.py: it imports `pipeline` (and
  `closed_formula`) from /tmp/rn, /tmp/z12, /tmp/mix, which do not exist, and fails with `ModuleNotFoundError`.
  [fix] Supply that module or say where it comes from. My recomputation of its five curves (`h2_rerun_z13.py`,
  rank_second against my own formula (2), 54 s) gives exactly the numbers in `h2_check.txt`:
  - (725,497) E_LU fails H2: 32 agree, 32 disagree.
  - (725,497) E_08: 64 agree. (535,23) E_08: 64. (465,416) E_LU: 256. (465,416) E_08: 512.
  - The (A, C) values and Sel² dimensions also agree.
- **(5) Tangent-form assignment: references right, conclusion right, one sentence false as written [fix].**
  - I fetched arXiv:0710.2079 (Fisher–Schaefer–Stoll, LMS J. Comput. Math. 13 (2010) 451–460). Its sections are
    1 Introduction, 2 Two definitions, 3 The Cassels pairing, 4 The main diagram, 5 The main theorem.
  - §3 defines ⟨a, a′⟩_Cas = ∏_v (f(P_v), M)_{A_v}, with f: T_j ↦ f_j, f_j defined over K(θ_j),
    div f_j = 2𝔇_j, [𝔇_j] ↦ T_j, and M the Kummer image of a′. §5 builds the covering obtained by setting each f_j
    equal to a square. Theorem 5.2 (with Lemma 5.1, the only numbered statements of §5) says
    Cassels = ⟨,⟩_2 = ⟨,⟩_1 = Cassels–Tate. So §3, §5 and Theorem 5.2 are cited correctly.
  - The false sentence: "the conic that omits u_j … its tangent form (over a fixed linear form) is f_j". Under FSS's
    normalisation it does not hold. With F_j the fibre over the tangency point, F_u ~ (T_1) + (O) and so on, so
    l_u/l_w has divisor class T_2, not T_1. FSS's f is (l_v l_w, l_u l_w, l_u l_v) mod squares. No choice of l_0 makes
    (l_u, l_v, l_w)/l_0 have classes (T_1, T_2, T_3).
  - The conclusion still holds: (l_v l_w, l_u l_w, l_u l_v) is (l_u, l_v, l_w) times the common factor l_u l_v l_w,
    and a common factor g changes the pairing by (g, b_1 b_2 b_3) = 1. So ∏_j (l_j, b_j) is FSS's pairing, which is
    formula (1). This agrees with z12's derivation and with rank_second.
  - [note] "the only one of the three defined over K(θ_j)" fails for a cyclic cubic, where all three conics are
    defined over the same field. In the split case it is a heuristic. The intrinsic reason is that the conic omitting
    u_j is the quotient of the covering by R ↦ T_j − R.
- **(6) §5 rewritten: applied, correct.**
  - "Exactly when its δ lies in the radical" holds.
  - (1, ±1, ±1) is the class of a point of order 2 (T_r, resp. T_*).
  - "a rank statement … does not change that" is right about the class. [note] The tuple itself can still be removed
    when E(Q) is known (for example rank 0), through the finitely many fibres over E(Q). That was the first note's
    "rank-0 statement" at (5623,1017).
  - "decides nothing" where the radical is larger: right. The partners sentence: right.
  - New, unverified claim: the EXPLORATION_4 paragraph (6,904 tuples, 937 ratios, 10 survivors). That file was not
    supplied, so I did not check it. It is consistent with the pilot's four survivors, which I checked: δ_08 = 1,
    δ_LU = sqf(2(m² − n²)) at (5623,1017), and δ_08 = sqf(2mn) = 19,505, δ_LU = 1 at (3320,2303).
- **(7) §7: Fermat line and symbol applied, correct.** mn = 2k² with gcd 1 forces {m, n} = {2a², b²}, and
  4a⁴ + b⁴ = □ has no solution with ab ≠ 0. −□ is excluded because 2mn > 0. "(2/p) times the quartic residue symbol"
  is right (p | S odd, so p ≡ 1 mod 4, and β is a unit at p). [fix, not applied] Before (2) can be used on
  (1, 2mn, 2mn) at a non-square ratio, the text should say that class must be in S′, a Legendre condition (in S′ on 32
  of 93 of my opposite-parity E_08 curves).
- **(8) Fisher 2017: applied, correct** (Math. Comp. 86 (2017) 2493–2518, Theorem 2.1).
- **Prior art: one new error [fix].**
  - arXiv:2303.05059 is S. Zhang, "On a comparison of Cassels pairings of different elliptic curves". That is exactly
    the paper the first note described (e1 + e2 + e3 = 0; E and 𝓔 with e1a² + e2b² + e3c² = 0). So the first note's
    number was right.
  - 2303.05058 is a different paper: Z. Wang and S. Zhang, "On the quadratic twist of elliptic curves with full
    2-torsion", y² = x(x − a²)(x + b²) with 2(a² + b²) a square. That description, and "sibling family", are right.
  - The sentence "In my first note the Wang–Zhang number was given as 2303.05059; 2303.05058 is the correct one"
    should go. List both papers instead.
  - Arxiv:1806.06250 is Stevenhagen, "Rédei reciprocity, governing fields, and negative Pell": right. I did not find
    the Mokrani thesis in one search; the note's "reported … I have not opened it" is the right status.
  - [note] The heading credits "the reviewer's separate search". That search was not mine (z13 searched only for
    Fisher and Monsky).
- **Not applied from my first review:** finding 5 (§4: 720 lines = 90 distinct and 1,434 = 154; "committed values"
  means ct08's stored values; on E_LU every tuple row was skipped); finding 6 (`shape.txt`, not in this delivery);
  the S′ step of finding 9 (above). [note] "Matrices of Legendre symbols" still omits the 2-adic entries. Header,
  line 3: "not reviewed" is now out of date.

Script: `experiments/referee_z13/h2_rerun_z13.py` → `h2_rerun_z13.txt` (54 s; deterministic, can run in CI).
