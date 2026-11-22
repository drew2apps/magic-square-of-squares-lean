# The 3×3 magic square of squares: a statement on the centre, a machine-checked reduction, and a computed census

Drew Carlson. 22 November 2026, release 0.3.4.
Repository: https://github.com/drew2apps/magic-square-of-squares-lean. Archive: https://doi.org/10.5281/zenodo.24435351.

The Lean code and this text were written with Claude, an AI model made by Anthropic. I am not a number theorist. Nobody else has checked the mathematics outside Lean. This note says what is proved, what is only computed, and what is unchecked, so that someone who knows the subject can check it.

## 1. The problem

Is there a 3×3 magic square whose nine entries are distinct perfect squares? The question is open. Write the square of roots as

```
a b c
d e f
g h i
```

so that the rows, columns and two diagonals of `a², …, i²` have the same sum `3e²`. Lee Morgenstern showed by computation (2006) that every entry of a solution exceeds 10¹⁴. The only known square with seven square entries is Bremner's (1999). A 2025 preprint (O. Hill, arXiv:2510.08286) claims a proof of nonexistence; it has not been verified here, this note takes no position on it, and nothing here depends on it.

This note does not solve the problem. It contains:

- Theorems 1 to 3, on the centre of a solution. They are short, their proofs are below, and Lean checks them.
- Theorems 4 to 6, a reduction of the problem to a one-parameter family of curves with a kernel-checked list of excluded cases.
- Computations A to E, a census of the remaining cases, not formalized.

## 2. Proved in Lean

Lean 4 with Mathlib, 257 declared theorems, no `sorry`, no `native_decide`, standard axioms only (checked in CI for every theorem). Statements exist both for nine integers and for the matrix form used in google-deepmind/formal-conjectures. The library also contains the nine necessary conditions listed on Wikipedia (most of them are in Rabern 2003) and the classical reformulations: three arithmetic progressions of squares with a common difference, rational right triangles of equal area, and points on the curves `y² = x³ − n²x`.

### 2.1 The centre

Throughout, `a, …, i` are integers whose squares form a magic square; nothing is assumed about distinctness or common factors unless stated. For a square with nine distinct entries one can take the roots positive and, after transposing the square if necessary, `f > h`; then all hypotheses below hold.

**Theorem 1** (`IsMSS.edge_params`, `IsMSS.center_factorization`, `IsMSS.edge_sums_identity`). Let `h + f ≠ 0`, and let `m ≥ 1`, `n ≥ 0` be coprime with `(b + d)·n = (h + f)·m`. Then there are integers `u`, `v` with

`b + d = mu`, `h + f = nu`, `b − d = nv`, `f − h = mv`, `8e² = (m² + n²)(u² + v²)`.

If also `f ≠ h`, and `m′ ≥ 0`, `n′ ≥ 0` are coprime with `(b + d)·n′ = (f − h)·m′`, then there is an integer `k ≠ 0` with

`8e² = (m² + n²)·k²·(m′² + n′²)`;

in fact `u = k·m′` and `v = k·n′` (a step of the Lean proof, not part of the Lean statement).

*Proof (a sketch in Gaussian integers; the Lean proof is organised differently).* Opposite entries satisfy `b² + h² = d² + f² = 2e²`. Put `β = b + hi`, `δ = d + fi`, `μ = m + ni` in `ℤ[i]`. Then `β + δ = uμ` with `u` rational, and `u` is an integer because `m`, `n` are coprime. Since `β` and `δ` have the same norm, `(β − δ)·conj(β + δ)` is purely imaginary, so `β − δ = −iv·μ` with `v` rational, and again an integer. Hence `2β = μ(u − iv)`, `2δ = μ(u + iv)`; real and imaginary parts give the four linear relations, and norms give `8e² = (m² + n²)(u² + v²)`. The second ratio is `(b + d) : (f − h) = u : v`. ∎

In terms of the entries the factorisation is the identity `((b + d)² + (h + f)²)·((b + d)² + (f − h)²) = 8e²·(b + d)²`, valid whenever `b² + h² = d² + f² = 2e²` (`IsMSS.edge_sums_identity`).

**Theorem 2** (`IsMSS.norm_dvd_center`, `IsMSS.norm_isSquare_or_prime_dvd_center`). Let `h + f ≠ 0`, `b ≠ d`, and `m ≥ 0`, `n ≥ 1` coprime with `(b + d)·n = (h + f)·m`. Then `m² + n²` divides `8e²`, every odd prime dividing `m² + n²` divides `e`, and either `m² + n²` is a square or twice a square, or some prime `p ≡ 1 (mod 8)` divides `e` and divides `m² + n²` to an odd power.

*Proof.* The first two statements follow from Theorem 1 (they are trivial if `m = 0`). Let `p` be an odd prime dividing `S = m² + n²` to an odd power. It is `≡ 1 (mod 4)`, since `m`, `n` are coprime. From `S(u² + v²) = 8e²`, `p` divides `u² + v²`. The corner `a` satisfies `2a² = f² + h²`, that is, `4a² = n²u² + m²v²`. Modulo `p`, `n² ≡ −m²` and `v² ≡ −u²`, so `4a² ≡ −2m²u²`. If `p ∤ u`, then `−2` is a square modulo `p`, hence so is `2`, and `p ≡ 1 (mod 8)`. If `p ∣ u`, then `p ∣ v`, and `p` divides all nine roots; divide the square by `p` and repeat, which stops because `u ≠ 0`. So every prime `≡ 5 (mod 8)` divides `S` to an even power (a statement of the Kepler71 project), and if `S` is not a square or twice a square, some prime divides it to an odd power. This is a sketch; the Lean proof is organised differently. ∎

**Theorem 3** (`IsMSS.centre_two_ratios`). Under the hypotheses of Theorem 2, let also `f ≠ h`, and `m′ ≥ 0`, `n′ ≥ 1` coprime with `(b + d)·n′ = (f − h)·m′`. Then either some prime `p ≡ 1 (mod 8)` divides the centre root `e` and divides both `m² + n²` and `m′² + n′²` to an odd power, or one of these two numbers is a square and the other is twice a square.

*Proof.* The square with `h` replaced by `−h` is again a magic square of squares, and its first ratio is `m′ : n′`. So Theorem 2 applies to both numbers: if no prime `≡ 1 (mod 8)` divides `e`, each is a square or twice a square. By Theorem 1, `8e² = (m² + n²)·k²·(m′² + n′²)` with `k ≠ 0` and both numbers positive. If both were squares, `8e²` would be a nonzero square; if both were twice a square, `32e²` would be. Neither is possible, since `x² = 2y²` has only the zero solution. ∎

**Theorem 3′** (`IsMSS.two_mul_norms_isSquare`, `IsMSS.norms_parity`, `IsMSS.norms_padicVal_parity`). Under the same hypotheses, `2(m² + n²)(m′² + n′²)` is a square; exactly one of the two numbers is even; and for every odd prime `p` the exponents of `p` in the two numbers have the same parity.

*Proof.* By Theorem 1, `2(m² + n²)(m′² + n′²) = (4e/k)²`. A sum of two coprime squares is odd or twice an odd number, and neither twice a product of two odd numbers nor eight times one is a square. The last statement compares the exponents of `p` on the two sides. ∎

In other words the two numbers have the same odd squarefree part `c`. By Theorem 2 every prime factor of `c` is `≡ 1 (mod 8)` and divides `e`, and the two alternatives of Theorem 3 are `c ≠ 1` and `c = 1`. The sentence about `c` is a restatement in words; Lean has the three statements of Theorem 3′.

**Remarks.**

1. *The exceptional case (this remark is not in Lean).* If no prime `≡ 1 (mod 8)` divides `e`, write the two numbers as `r²` and `2r′²`. Then `2e = ±r·k·r′`, `r ≥ 5`, and every prime divisor of `r` and of `r′` is `≡ 5 (mod 8)` and divides `e`. The proof does not exclude this case, and I do not know whether it can occur.
2. *Relation to known results.* Every prime divisor of the centre of a primitive solution is `≡ 1 (mod 4)` (Rabern, Theorem 1.2; Woll, Theorem 3.1). A prime dividing the centre and a mid-edge entry is `≡ 1 (mod 8)` (Woll, Corollary 4.3). A prime dividing only the centre has three consecutive nonzero quadratic residues (Weisenberg, Theorem 2.10; Woll). These statements start from a prime dividing the centre. Theorems 2 and 3 say where such primes come from: the odd primes dividing `m² + n²` and `m′² + n′²` divide the centre root, those dividing to an odd power are `≡ 1 (mod 8)` and are the same for the two numbers, and if there are none the solution is of the exceptional type. As a statement about the centre root alone, "it has a prime divisor `≡ 1 (mod 8)`, or all its prime divisors are `≡ 5 (mod 8)`" already follows from Rabern's theorem; the content here is the link with the two ratios. Rabern closes his paper with the remark that one might conjecture that all prime divisors of all entries of a solution are `≡ 1 (mod 8)`, and that it would be interesting to prove the opposite, that some entry has a prime divisor `≡ 5` or `7 (mod 8)`. A solution of the exceptional type would be of the second kind.
3. *Novelty.* I did not find Theorem 1, the divisibility `m² + n² ∣ 8e²` or Theorem 3 in Rabern (2003), Woll (2018), Weisenberg (2023), Bremner (1999, 2001), Auel and Singer (2026), Kominers, or the multimagie pages. Pierrat, Thiriet and Zimmermann (2015), Robertson (1996) and Boyer (2005) have not been read. I claim no novelty; the proofs are elementary and the statements may be known in another form.
4. *What the factorisation is.* `2(m² + n²)(m′² + n′²)` being a square is the condition that the form `G4` of Section 2.2 is a square at `t = m′/n′`, that is, that the centre entry is a square, written multiplicatively. Every search already imposes it: a search by centre starts from a square centre and pairs of ways of writing twice the centre as a sum of two squares, and the five-form test contains `G4`. So these statements give an arithmetic reading of the centre condition, not a new search method.

### 2.2 The reduction

By Theorem 1 a solution is determined by the ratio `m : n` and the rational number `t = u/v`, and five quadratic forms in `t` must be squares: `F0 = m² + n²t²`, `F8 = n² + m²t²`, `G4 = 2(m² + n²)(t² + 1)`, `GL = G4 − 8mnt`, `GU = G4 + 8mnt` (`IsMSS.fullSquares`). For fixed `m : n` this is a curve of genus 49. In the coordinates `a, …, i` it is a component of the hyperplane section `n(b + d) = m(h + f)` of the surface studied by Auel and Singer (my remark, checked by hand; see `G1_FULL_SYSTEM.md`). The pair `(m + n, m − n)`, halved if both are even, gives the same system; it is the *partner*. The forms and the five-form local test come from the Kepler71 project.

**Theorem 4** (`G1.fullTable_valid`, `G1.no_solution_of_ratio`). For each of 484 residue classes of `m : n`, at 34 primes of which the largest is 607, the five forms are not simultaneously squares for any rational `t`. Each class is certified by a finite computation that the Lean kernel evaluates.

**Theorem 5** (`IsMSS.ratio_ne_of_le_1000`). In a magic square of squares with `b ≠ d` and `h + f ≠ 0`, the ratio `(b + d) : (h + f)` is not `m : n` for any coprime `1 ≤ n < m ≤ 1000` outside an explicit list of 425 pairs. There are 304,191 such pairs in all. The size of the entries plays no role.

**Theorem 6** (`IsMSS.edge_root_gt_500`). If `b² ≠ d²`, `h² ≠ f²` and `b² ≠ f²`, one of the edge roots `b, d, f, h` exceeds 500 in absolute value. This is far below Morgenstern's bound; it is included because it is machine-checked.

Also in Lean: `m² + n²` is a square modulo 7, 17, 31 and 47 (`IsMSS.edge_norm`); for a prime `p ≡ ±3 (mod 8)` every ratio has a point modulo `p`, so no ratio is excluded by a congruence modulo such a prime (`G1.exists_point_mod_of_not_isSquare_two`, and `G1.exists_point_mod_of_three_five_mod_eight` for the statement in terms of `p` modulo 8; powers of `p` are not covered); and the algebra used in Section 3 (membership of the points on the three elliptic curves, the map from the quartic `v² = F0·GL` to a Weierstrass model, the partner identities, the identities behind the singular fibres, and the case analysis of the side search).

**What Lean does not prove about that algebra.** These are statements about points and polynomials. `e0l_weierstrass` sends a point of the quartic with `t ≠ m/n` to a point of a Weierstrass curve; it does not treat the base point, does not prove that the two curves are isomorphic, and does not cover the change to the minimal model on which the computations run (a script checks that the recorded models are isomorphic). The resultant and octic statements are polynomial identities; what they imply modulo a prime is argued in the text. No Lean theorem transfers a closure from a ratio to its partner. Ranks, generators and heights are not in Lean at all.

## 3. Computed, not formalized

Details, limits and commands are in `G1_FULL_SYSTEM.md`. Programs: PARI 2.17, Sage (passagemath) with eclib, Magma 2.29 through its public calculator, and C and Python programs in the repository.

In this section a number is *computed*: a program returned it, and anyone can rerun the program. Whether the algorithm behind the program proves the number unconditionally is a separate question, and it is the subject of Section 4.

**A. Which primes can obstruct.** For a ratio of good reduction the five conditions define a smooth curve of genus 49 over `F_p`. The Weil bound, a separate character-sum bound for the fibres where two branch points collide (`p` dividing `m⁸ + 14m⁴n⁴ + n⁸`), and a direct computation below those bounds show that exactly the 34 primes of Theorem 4 have an obstructed class. All 34 are `±1 (mod 8)`; that no other residue is possible is in Lean.

**B. Local solubility.** The 425 ratios of Theorem 5 are exactly those not excluded by the classes and by the condition on primes `≡ 5 (mod 8)`. A separately written local test agrees with the classification on all 304,191 ratios. With A, the 425 ratios are soluble in every `ℚ_p`, so no congruence closes any of them. The excluded ratios have density 1 (argued, with a proof sketch).

**C. Rank 0.** 268 of the 425 ratios are closed by a quotient elliptic curve of rank 0.

**D. Rank 1 with multiplicity.** Three quotient elliptic curves occur more than once in the Jacobian: `E04: y² = x(x + 2Sm²)(x + 2Sn²)` with `S = m² + n²` (twice), a second curve `D1` (twice), and `E0L: y² = F0·GL` (four times). When such a curve has rank exactly 1, a solution gives two points whose heights are tied to each other, and the Demjanenko–Manin argument leaves finitely many candidates, which are tested exactly. This closes 121 more ratios (91, 24 and 6 on the three curves). For every one of them the bound "rank at most 1" is given by Magma (87 distinct curves); for all but two it is also given by PARI, and for most of the first curve by eclib. For two ratios, (631,71) and its partner, PARI gives only `[1, 3]` and the bound rests on Magma alone.

**E. Open.** 36 ratios, 5 of them with `m ≤ 200`: (72,65)~(137,7), (91,60)~(151,31), (195,28). Seven have no trivial point. For these the curves above have rank at least 2. An experiment with a quotient of genus 2 through the Magma calculator closed nothing: for 20 of its 23 instances the calculator's limits were reached.

**A side search.** Among magic squares with nine distinct positive entries, a square centre of root at most 10⁹ and at least seven square entries, the only ones are Bremner's square and its multiples. This is not new: B. Chaffin's page reports a search of the same family to centre root 10¹⁴ (about 2008, published 2019) with nothing else found (read from his page with an AI model). The search here re-verifies a small part of that range with an independent program. That its six cases are all the cases is proved in Lean; the program is not verified, and its output is evidence of a computation, not a theorem.

## 4. What has not been checked

In order of importance:

1. **Heights and generators.** The closures in D need the generator of each curve (saturation with a proved height lower bound) and explicit comparisons between canonical and naive heights. The derivations were written with an AI model. Scripts check what can be checked without a specialist: the recorded points, models, candidates and bounds are consistent; Silverman's constant, which the software documents for one side of the height difference, bounds the other side as well for all 121 models (in exact rational arithmetic, from statements of Cremona, Prickett and Siksek, read from their paper with an AI model, and a bound of Silverman as quoted there); the listed range of multiples exceeds the bound computed from the recorded heights by at least 9% in every case; and the heights of the generators are stable from 53 to 200 bits. The heights themselves are floating-point values without a proved error bound. `AUDIT_REQUEST.md` lists what a reader is asked to check, with three closures worked out.
2. **The character-sum bound** for the singular fibres and the Weil-bound step in A.
3. **The density argument** in B, which is a sketch.
4. **Rank bounds** are computed, by two or three programs, not proved here. That a program returns a bound is not a proof that its algorithm gives that bound unconditionally; the audit request asks about this.
5. **The literature.** Robertson (1996), Boyer's article (2005) and the preprint of Pierrat, Thiriet and Zimmermann (2015) have not been read.

Three errors were found and corrected during the work, all recorded in `CHANGES.md`: the singular fibres were overlooked in the first statement of the 34 primes; the first version of the side search missed two of six cases; and the curve of genus 2 in the experiment lacked a factor `m² + n²` (no recorded number was affected). Others may remain.

## 5. Reproducing

`lake build` and `lake env lean scripts/Audit.lean` check the Lean part (about 15 minutes from an empty cache). CI also regenerates the certificate files and runs the evidence checker, the check of the singular fibres, the height check and the side search on a small range; these are regression tests, not a rerun of the computations. The elliptic-curve computations are rerun with the scripts in `scripts/g1_full/`; they need PARI and Sage, and Magma for the rank bounds; their outputs are committed, with verbatim Magma responses for 41 of the 87 curves.

## References

- A. Auel and B. Singer, The algebraic geometry of 3-by-3 magic squares of squares, arXiv:2609.09351.
- C. Boyer, 3×3 magic square of squares search, multimagie.com (read 22 November 2026); reports Morgenstern's results.
- A. Bremner, On squares of squares, Acta Arith. 88 (1999), 289–297; On squares of squares II, Acta Arith. 99 (2001), no. 3.
- N. Bruin, J. Thomas and A. Várilly-Alvarado, Explicit computation of symmetric differentials and its application to quasihyperbolicity, Algebra Number Theory 16 (2022), 1377–1405.
- B. Chaffin, Magic square of squares, https://benchaffin.com/magic-squares/magic-squares.html (2019).
- J. E. Cremona, M. Prickett and S. Siksek, Height difference bounds for elliptic curves over number fields, J. Number Theory 116 (2006), 42–68.
- O. Hill, On arithmetic progressions and a proof of the nonexistence of magic squares of squares, arXiv:2510.08286.
- Kepler71, notes on the full nine-cell system, https://github.com/Kepler71/magic-square-of-squares (cited in `G1_FULL_SYSTEM.md`).
- S. D. Kominers, Center-zero magic squares of squares over number fields, preprint.
- P. Pierrat, F. Thiriet and P. Zimmermann, Magic squares of squares, preprint, 2015 (not read).
- L. Rabern, Properties of magic squares of squares, Rose-Hulman Undergrad. Math. J. 4 (2003).
- J. P. Robertson, Magic squares of squares, Math. Mag. 69 (1996), 289–293 (not read).
- J. H. Silverman, The difference between the Weil height and the canonical height on elliptic curves, Math. Comp. 55 (1990), 723–743 (not seen; the bound used is quoted in Cremona, Prickett and Siksek, Theorem 15).
- D. Weisenberg, Some thoughts on the 3×3 magic square of squares problem, Rose-Hulman Undergrad. Math. J. 24 (2023).
- C. Woll, A partial residue categorization of the magic square of squares, arXiv:1809.03067.
