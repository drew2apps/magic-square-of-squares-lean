# The class of a solution: what is proved for every ratio

Two pages, without the history. Written by Claude (an AI model made by Anthropic), directed by Drew (`drew2apps`), 24 November 2026. Every statement below is a theorem of `MagicSquareOfSquares/CoverDescent.lean` or `CoverDescentTable.lean` (Lean 4 with Mathlib, no `sorry`, no `native_decide`, the three standard axioms), and each has a proof by hand in `LAST_PAIR.md`, section 17. **No mathematician has read any of it.** The problem of the 3×3 magic square of squares is open and is not solved here. The methods are classical (a 2-descent written with norms; quadratic residues); we do not know whether the statements are new, and we claim only that they are proved and checked.

## 1. Setting

`m > n > 0` are coprime integers (the ratio `m : n`), and `u, v` are coprime integers. Four pairs of integers, read as Gaussian integers `x + iy`:

    P0 = (mv, nu),   P8 = (nv, mu),   PL = ((m + n)(u − v), (m − n)(u + v)),   PU = ((m − n)(u − v), (m + n)(u + v)).

Their norms `x² + y²` are the forms `F0, F8, GL, GU`. In a 3×3 magic square of squares with `(b + d) : (h + f) = m : n` there are coprime `u, v` at which all four are squares (in Lean: `G1.exists_coprime_squares_full`, from `IsMSS`). Call such `(u, v)` a *solution on the four forms*. The **eight products** are `P0 · P8^± · PL^± · PU^±`, each factor or its conjugate.

For an integer `d`, "`d` times `x + iy` is a Gaussian square" means `d·x = a² − b²` and `d·y = 2ab` for some integers `a, b` (`GSq`). Nothing else about Gaussian integers is used: no unique factorisation, no ideal.

## 2. The theorems

**Theorem 1 (the class exists, and where its primes are).** At a solution on the four forms there is an integer `d > 0` without square factors such that `d` times each of the eight products is a Gaussian square, and every prime factor of `d` divides `2mn(m − n)(m + n)`. — `exists_class_prod`, for all integers `m, n`.

This is the step that a 2-descent usually takes from unique factorisation in `ℤ[i]`. Here it is one identity (if `x² + y² = w²` then `2(x + w)(x + iy) = (x + w + iy)²`) and two divisibility lemmas, and it is in Lean. The list of candidates for `d` is therefore produced inside Lean: a checker receives a list of primes, and the kernel verifies that it contains every prime of `2mn(m − n)(m + n)`.

**Theorem 2 (E).** Every odd prime factor of `d` is `≡ ±1 (mod 8)`. — `mod_eight_of_dv`. Uses: `F0` and `GL` squares; `m, n` coprime; `u, v` coprime.

**Theorem 3 (F).** `d` is even. — `two_dvd_class`. Uses: the four forms squares; the eight products; both coprimalities.

**Theorem 4 (G).** `d` is a square modulo every odd prime dividing `m² + n²`. — `isSquare_class_mod`. Uses: the eight products and both coprimalities only.

**Theorem 5 (O).** `d` is a square modulo every odd prime dividing `m⁸ + 14m⁴n⁴ + n⁸`. — `isSquare_class_mod_octic`. Uses: the eight products and both coprimalities only.

Together: `exists_class_struct`. Theorems 3, 4 and 5 hold for any `d` with the property of the eight products, not only for the one of Theorem 1. The octic is, up to constants, the product of the resultants of `F0, F8` with `GL, GU`; its primes do not divide `2mn(m² − n²)(m² + n²)`.

**The criterion.** So a solution on the four forms gives a set `S` of primes `≡ ±1 (mod 8)` dividing `mn(m² − n²)` with

    ∏_{q ∈ S} (q/p) = (2/p)    for every odd prime p dividing (m² + n²)(m⁸ + 14m⁴n⁴ + n⁸),

a system of linear equations over `𝔽₂`. If it has no solution, the ratio `m : n` does not occur. If it has exactly one, the class of every solution is known.

**Theorem 6 (the last pair).** For `m : n = 627 : 364` and for `991 : 263`, at every solution on the four forms `521266` times each of the eight products is a Gaussian square. — `class_of_mem_structTable`. The proof is Theorems 1 to 5 and three primes: 5 and 29 divide `m² + n²`, 37 divides the octic, and of the eight even candidates only `521266 = 2·263·991` is a square modulo all three.

**Theorem 7 (ratios to 1000).** In a 3×3 magic square of squares (rows, columns and both diagonals of the nine squares with one sum; entries not assumed distinct) with `h + f ≠ 0`, `b ≠ d` and `(b + d)·n = (h + f)·m` for coprime `1 ≤ n < m ≤ 1000`, the number `m² + n²` is a square or twice a square. — `IsMSS.ratio_norm_sq_or_two_sq_of_le_1000`.

Theorem 7 is a finite statement. Its proof: congruences leave 425 ratios (`G1.exceptions1000`, in release 0.5.0); for the 119 of them with `m² + n²` neither a square nor twice a square the kernel checks a certificate that no class is possible: for 115 with one class `d` (Theorems 1 and 2, then finitely many points modulo small prime powers, or modulo a prime by Euler's criterion), for 4 with one class for each of the four pairs (`exists_classes4`). For the 306 other ratios a degenerate solution exists (for instance `u = v` when `m² + n²` is a square) and no such certificate can.

## 3. What the criterion decides, computed (not in Lean except Theorem 6)

`scripts/cover_descent/structure/class_struct.py`, with the primes of the octic below 10⁵:

- Of the 119 ratios of Theorem 7, **54 are excluded by the criterion alone**. For the other 65 the certificate needs a small prime outside the criterion, which we do not explain.
- Of the 306 ratios with a degenerate solution, the system is solvable for all (as it must be), and **for 150 it has exactly one solution**: the class of any solution is that of the degenerate one.
- For random coprime ratios with `m` between 10² and 10⁶ the criterion excludes about four in five, most of them by Theorem 4 at a prime `≡ 5 (mod 8)`, which is a condition already known from the fifth form. We have no theorem and no conjecture about the proportion.

## 4. What this is not

- Not a proof that there is no magic square of squares, and not a step that reduces the problem to something known. For the 306 ratios to 1000 with a degenerate solution, and for every ratio above 1000, nothing is proved here beyond the criterion.
- Theorem 6 gives the class for the last pair. The argument that then excludes that pair (`LAST_PAIR.md`, Theorem A: a curve of genus 9 with two maps to an elliptic curve of rank 1, an explicit height bound from the literature read by an automated reader, a Mordell–Weil sieve) is not in Lean and has not been read by anyone.
- The census of the repository (release 0.5.0) is not changed by any of this.

## 5. For a reader with one hour

Theorem 5 is the least expected and has a short proof: `LAST_PAIR.md`, section 17, "Theorem O", fifteen lines and two polynomial identities. Theorem 4 is of the same kind and shorter. If either is wrong, or is in the literature, we would like to know.
