# The class of a solution: what is proved for every ratio

Two pages, without the history. Written by Claude (an AI model made by Anthropic), directed by Drew (`drew2apps`), 24 November 2026. Every statement below is a theorem of `MagicSquareOfSquares/CoverDescent.lean` or `CoverDescentTable.lean` (Lean 4 with Mathlib, no `sorry`, no `native_decide`, the three standard axioms). Theorems 2 to 5 also have a proof by hand in `LAST_PAIR.md`, section 17; Theorem 1 in sections 11 and 15; Theorems 6 and 7 are finite checks by the kernel on top of them. **No mathematician has read any of it.** The problem of the 3×3 magic square of squares is open and is not solved here. The methods are classical (a 2-descent written with norms; quadratic residues); we do not know whether the statements are new, and we claim only that they are proved and checked.

## 1. Setting

`m > n > 0` are coprime integers (the ratio `m : n`), and `u, v` are coprime integers. Four pairs of integers, read as Gaussian integers `x + iy`:

    P0 = (mv, nu),   P8 = (nv, mu),   PL = ((m + n)(u − v), (m − n)(u + v)),   PU = ((m − n)(u − v), (m + n)(u + v)).

Their norms `x² + y²` are the forms `F0, F8, GL, GU`. In a 3×3 magic square of squares with `(b + d) : (h + f) = m : n` there are coprime `u, v` at which all four are squares (in Lean: `G1.exists_coprime_squares_full`, from `IsMSS`). Call such `(u, v)` a *solution on the four forms*. The **eight products** are `P0 · P8^± · PL^± · PU^±`, each factor or its conjugate.

For an integer `d`, "`d` times `x + iy` is a Gaussian square" means `d·x = a² − b²` and `d·y = 2ab` for some integers `a, b` (`GSq`). Nothing else about Gaussian integers is used: no unique factorisation, no ideal.

## 2. The theorems

**Theorem 1 (the class exists, and where its primes are).** At a solution on the four forms there is an integer `d > 0` without square factors such that `d` times each of the eight products is a Gaussian square, and every prime factor of `d` divides `2mn(m − n)(m + n)`. — `exists_class_prod`. (In Lean the hypotheses are: `u, v` coprime, the four norms squares, and none of the four pairs is `(0, 0)`, which holds for `m > n > 0` as soon as `u, v` are not both 0; `m, n` need not be coprime for this one.) Such a `d` is unique (two of them have a quotient that is a rational square in `ℚ(i)`, hence `±` a rational square); Lean does not state the uniqueness and does not need it, because Theorems 3, 4, 5 are proved for every `d` with the property of the eight products.

This is the step that a 2-descent usually takes from unique factorisation in `ℤ[i]`. Here it is one identity (if `x² + y² = w²` then `2(x + w)(x + iy) = (x + w + iy)²`) and two divisibility lemmas, and it is in Lean. The list of candidates for `d` is therefore produced inside Lean: a checker receives a list of primes, and the kernel verifies that it contains every prime of `2mn(m − n)(m + n)`.

**Theorem 2 (E).** The `d` of Theorem 1 can be taken with every odd prime factor `≡ ±1 (mod 8)` (and it is the same `d`, by the uniqueness). — `exists_class_prod8`; the step is `mod_eight_of_dv`: an odd prime dividing both coordinates of one of the four pairs is `≡ ±1 (mod 8)`, which uses `F0` and `GL` squares and both coprimalities.

**Theorem 3 (F).** `d` is even. — `two_dvd_class`. Uses: the four forms squares; the eight products; both coprimalities.

**Theorem 4 (G).** `d` is a square modulo every odd prime dividing `m² + n²`. — `isSquare_class_mod`. Uses: both coprimalities and the eight products only (in fact one of them, `P0·conj(P8)·PL·conj(PU)`; Lean's proof uses two).

**Theorem 5 (O).** `d` is a square modulo every odd prime dividing `m⁸ + 14m⁴n⁴ + n⁸`. — `isSquare_class_mod_octic`. Uses: both coprimalities and one of the eight products only (`P0·P8·PL·PU` for a prime of one quartic factor of the octic, `P0·P8·conj(PL)·conj(PU)` for the other).

**Theorem 5b (the law at 5).** If 5 divides `mn(m − n)(m + n)` and 25 does not, then `d ≢ ±1 (mod 5)`; for the `d` of Theorem 1, which 5 does not divide (Theorem 2), `d ≡ ±2 (mod 5)`: a non-residue. — `class_mod_five`. Uses: both coprimalities and the eight products only. The proof is a complete check modulo 25 by the kernel (`law5Check`): multiplying `(m, n)` and `(u, v)` by numbers prime to 5 multiplies each product by a fourth power, so it is enough to take both pairs in the normal form `(1, t)` or `(5s, 1)` modulo 25: 16 pairs `(m, n)`, 30 pairs `(u, v)`, 10 residues of `d`. When 25 divides one of `m, n, m ± n` there is no such law (both kinds of class pass modulo 625; tested, not proved), and at the other primes `q ≡ 1 (mod 4)` of `mn(m² − n²)` the analogous statement restricts nothing (tested for `q` to 20000).

Together: `exists_class_struct` (with the hypotheses of Theorem 1 and `m, n` coprime). Theorems 3, 4 and 5 hold for any `d` with the property of the eight products, not only for the one of Theorem 1. In Lean "a square modulo `p`" is `IsSquare (d : ZMod p)`, which is also true when `p` divides `d`; for the `d` of Theorem 1 that does not happen (an odd prime of `m² + n²` or of the octic divides none of `2, m, n, m ± n`), so it means that `d` is a quadratic residue. Each of the four resultants of `F0` or `F8` with `GL` or `GU` is 4 times the octic `W` (and the resultants of `F0` with `F8` and of `GL` with `GU` are `(m² − n²)²(m² + n²)²` and `1024·m²n²(m² + n²)²`, which is where `m² + n²` comes from); the odd primes of `W` do not divide `mn(m² − n²)(m² + n²)`. It is Klein's octahedral invariant: the Hessian of `f = mn(m⁴ − n⁴)`, with `W³ − 108f⁴` a square, and the discriminant of `F0·F8·GL·GU` is `2⁴⁸(fW)⁸`. So its primes are primes of bad reduction of the curves of the problem, where a condition is to be expected; Theorems 4 and 5 say what the condition is.

**The criterion.** So a solution on the four forms gives a set `S` of primes `≡ ±1 (mod 8)` dividing `mn(m² − n²)` with

    ∏_{q ∈ S} (q/p) = (2/p)    for every odd prime p dividing (m² + n²)(m⁸ + 14m⁴n⁴ + n⁸),

a system of linear equations over `𝔽₂`; the law at 5 adds the equation `∏_{q ∈ S} (q/5) = 1` when 5 divides `mn(m² − n²)` exactly once. If it has no solution, the ratio `m : n` does not occur. If it has exactly one, the class of every solution is known.

**Theorem 6 (the class, for 185 ratios with a degenerate solution; the last pair).** For `m : n = 627 : 364` and for `991 : 263`, at every solution on the four forms, the degenerate ones included (any coprime `u, v` with the four norms squares), `521266` times each of the eight products is a Gaussian square. — `class_of_mem_structTable`, `structTable_last_pair`. The proof is Theorems 1 to 5 and three primes: 5 and 29 divide `m² + n²`, 37 divides the octic, and of the eight even candidates only `521266 = 2·263·991` is a square modulo all three. The same theorem covers 185 of the 306 ratios to 1000 for which `m² + n²` is a square or twice a square (`structTable`, with the class of each; 150 without the law at 5): for these the criterion leaves one candidate, with primes below 10⁵, and the kernel checks it (`structOk`: the primes are primes, they divide `(m² + n²)(m⁸ + 14m⁴n⁴ + n⁸)`, and every other candidate fails the law at 5 or Euler's criterion at one of the primes).

**Theorem 6b (no solution, from the class alone; 67 ratios).** For 67 of the 119 ratios of Theorem 7 below, the four forms are squares at no coprime `u, v`, and the proof uses no computation at the points: the criterion leaves no candidate. — `no_point_of_structOk`, `structZero`, `no_point_of_mem_structZero`. The kernel checks the same things as for Theorem 6, with no class left. For these 67 the certificates of Theorem 7 are not needed; for the other 52 they are.

**Theorem 7 (ratios to 1000).** In a 3×3 magic square of squares (rows, columns and both diagonals of the nine squares with one sum; entries not assumed distinct) with `h + f ≠ 0`, `b ≠ d` and `(b + d)·n = (h + f)·m` for coprime `1 ≤ n < m ≤ 1000`, the number `m² + n²` is a square or twice a square. — `IsMSS.ratio_norm_sq_or_two_sq_of_le_1000`.

The theorem is about the one ratio `(b + d) : (h + f)` with the signs and the order as written (`0 < n < m`); the roots may have any signs, so it applies after a change of signs of `b, d, h, f`, but that reduction is not part of the statement. It is not vacuous: the squares built from three squares in arithmetic progression (entries not distinct) satisfy the hypotheses (a reviewer's example: rows 25 1 49 / 49 25 1 / 1 49 25).

Theorem 7 is a finite statement. Its proof: congruences leave 425 ratios (`G1.exceptions1000`, in release 0.5.0); for the 119 of them with `m² + n²` neither a square nor twice a square the kernel checks a certificate that no class is possible: for 115 with one class `d` (Theorems 1 and 2, then finitely many points modulo small prime powers, or modulo a prime by Euler's criterion), for 4 with one class for each of the four pairs (`exists_classes4`). For the 306 other ratios a degenerate solution exists (for instance `u = v` when `m² + n²` is a square) and no such certificate can.

## 3. What the criterion decides (in Lean: Theorems 6 and 6b; the rest computed)

`scripts/cover_descent/structure/class_struct.py`, with the primes of the octic below 10⁵:

- Of the 119 ratios of Theorem 7, **67 are excluded by the criterion alone** (Theorem 6b; 54 without the law at 5). For the other 52 the certificate needs a small prime outside the criterion, which we do not explain.
- Of the 306 ratios with a degenerate solution, the system is solvable for all (as it must be), and **for 185 it has exactly one solution**: the class of any solution is that of the degenerate one (these 185 are Theorem 6, in Lean; 150 without the law at 5).
- With the octic factored completely (not in Lean, where the primes are below 10⁵) the 54 become 58 and the 67 become 71.
- **Theorem C (the five laws are all that the bad primes say). A proof by hand, not read by a person: [`COMPLETE.md`](COMPLETE.md). Partly in Lean (`ClassLocal.lean`: the class formula `class_iff`, and two of the four families of points of the lower bounds, modulo every power of `p`); `COMPLETE.md`, section 5, says which part.** Call a prime *bad* if it divides `2mn(m² − n²)(m² + n²)(m⁸ + 14m⁴n⁴ + n⁸)`. A `p`-adic point at which the four norms are nonzero squares has a class in `ℚ_p^×` modulo squares and `−1`, and `COMPLETE.md` determines the set of classes of such points at every bad prime, for every ratio: one class at 2 and at the odd primes of `(m² + n²)(m⁸ + 14m⁴n⁴ + n⁸)`, and at an odd prime `p` of `mn(m² − n²)` a set that depends only on `p` modulo 8 (and, for `p = 5`, on whether 25 divides). Consequence: if every bad prime has such a point, a squarefree `d` supported on the primes of `2mn(m² − n²)` is the class of a point at every bad prime if and only if it satisfies Theorems 2 to 5b. "Only if" is those theorems read `p`-adically; "if" was a conjecture in commit `8cc39fd` and is proved by explicit points (a conic and the bound of Hasse for `y² = v³ + v`; points with `p | u`; for `p = 5` points with `u² + v² ≡ 0`). Checks: `class_local.py` computes exact `p`-adic classes at 211,894 random points and at the points of the proof; `class_complete.py` and `class_complete_random.py` are the earlier tests modulo small prime powers (425 ratios and 1,952 random ratios, about a million candidates, no exception). So what the certificates of Theorem 7 use beyond the criterion is at primes of good reduction.
- **Theorem F′ (Lean: `class_mod_sixteen`).** Under the hypotheses of Theorem 3, any `d` with the property of the eight products is `≡ 0, 2, 8` or `14 (mod 16)`. For the class of Theorem 1 this follows from Theorems 2 and 3; it is here because Theorem C needs it at the prime 2 for every `d`. A check modulo 64 by the kernel.
- On random ratios with complete factorisations the criterion excludes about half of those that survive every congruence condition of the repository, with no visible trend from 10³ to 10⁶ (`LAST_PAIR.md`, section 17). We have no theorem about the proportion.

## 4. What this is not

- Not a proof that there is no magic square of squares, and not a step that reduces the problem to something known. For the 306 ratios to 1000 with a degenerate solution, and for every ratio above 1000, nothing is proved here beyond the criterion.
- Theorem 6 gives the class for the last pair. The argument that then excludes that pair (`LAST_PAIR.md`, Theorem A: a curve of genus 9 with two maps to an elliptic curve of rank 1, an explicit height bound from the literature read by an automated reader, a Mordell–Weil sieve) is not in Lean and has not been read by anyone.
- The census of the repository (release 0.5.0) is not changed by any of this.

## 5. For a reader with one hour

Theorem 5 has a short proof: `LAST_PAIR.md`, section 17, "Theorem O", fifteen lines and two polynomial identities. Theorem 4 is of the same kind and shorter. In the language of the subject these are local conditions of a two-cover descent at primes of bad reduction, made explicit for this family; the conditions that remain, at small primes of good reduction, are the usual ones and have no such law. If either theorem is wrong, or is in the literature, we would like to know.
