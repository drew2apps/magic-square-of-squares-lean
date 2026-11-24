# The class at the bad primes: a complete description (Theorem C)

**Status.** A proof by hand, written on 24 November 2026 by an AI model (Claude) in the session that also wrote `CLASS.md`. **No person has read it.** Part of it is in Lean (section 5 says exactly which part); the theorem as stated is not. Every step that is a computation is repeated by `scripts/cover_descent/structure/class_local.py` (exact `p`-adic classes of points, not a test modulo small powers), and section 6 says what that program does and does not check. Experiments: nothing here is counted in the census, and nothing here shows that a magic square of squares does or does not exist.

It answers the question left open in `CLASS.md`, section 3 (there a conjecture, in commit `8cc39fd`): **at the primes of bad reduction, the five laws on the class are the whole local information.**

## 1. Statement

`m, n` are coprime integers with `mn(m² − n²) ≠ 0`. The four pairs, as Gaussian integers `x + iy`:

    α0 = mv + i·nu,   α8 = nv + i·mu,   αL = (m+n)(u−v) + i·(m−n)(u+v),   αU = (m−n)(u−v) + i·(m+n)(u+v),

with norms `F0, F8, GL, GU`, and the eight products `Q_f = α0 · α8^± · αL^± · αU^±` (`^−` is the conjugate). `W = m⁸ + 14m⁴n⁴ + n⁸`. A prime is *bad* if it divides `2mn(m² − n²)(m² + n²)W`.

For a prime `p` let `X_p` be the set of pairs `(u, v)` of `p`-adic integers, not both divisible by `p`, at which the four norms are nonzero squares in `ℤ_p`. An integer solution on the four forms is a point of every `X_p`.

**Lemma 1 (the local class).** For `(u, v) ∈ X_p` there is `D ∈ ℚ_p^×` such that `D·Q_f` is a square in `ℤ_p[i] = ℤ_p ⊗ ℤ[i]` for all eight `f`, and the `d ∈ ℚ_p^×` with that property are exactly the `±D·t²`, `t ∈ ℚ_p^×`. One can take `D = ∏_j (w_j + x_j)`, where `x_j` is the real part of the `j`-th pair and `w_j² = x_j² + y_j²`, the sign of `w_j` chosen with `w_j + x_j ≠ 0`.

So a point has a class `cl_p(u, v)` in `Γ_p = ℚ_p^× / ⟨−1, squares⟩`. For odd `p` write `ε` for a non-residue unit. `Γ_p = {1, ε, p, εp}` if `p ≡ 1 (mod 4)`, `Γ_p = {1, p}` if `p ≡ 3 (mod 4)` (there `−1` is a non-residue, so all units are trivial), and `Γ_2 = {1, 5, 2, 10}`. The sign is not a loss: `−1 = i²`, so over `ℤ` too `d` and `−d` have the property together, and `CLASS.md` takes `d > 0`.

**Theorem C.** Let `p` be a bad prime with `X_p ≠ ∅`. Then the image `cl_p(X_p)` is:

| `p` | image |
|---|---|
| odd, divides `(m² + n²)·W` | `{1}` |
| `2` | `{2}` |
| odd, divides `mn(m² − n²)`, `p ≡ 3 (mod 8)` | `{1}` |
| …, `p ≡ 7 (mod 8)` | `{1, p}` |
| …, `p ≡ 5 (mod 8)` | `{1, ε}`; but `{ε}` if `p = 5` and 25 does not divide `mn(m² − n²)` |
| …, `p ≡ 1 (mod 8)` | all of `Γ_p` |

and `X_p ≠ ∅` for every odd `p` dividing `mn(m² − n²)`.

**Corollary.** Suppose `X_p ≠ ∅` at every bad prime (otherwise the ratio `m : n` has no solution, by a congruence on the ratio in which no class appears). Let `d > 0` be squarefree with all prime factors dividing `2mn(m² − n²)`. Then `d` is the class of a point of `X_p` at every bad prime `p` **if and only if** `d` satisfies Theorems 2 to 5b of `CLASS.md`:

- E: every odd prime factor of `d` is `≡ ±1 (mod 8)`;
- F: `d` is even;
- G, O: `d` is a square modulo every odd prime dividing `(m² + n²)·W`;
- the law at 5: `d ≡ ±2 (mod 5)` if 5 divides `mn(m² − n²)` and 25 does not.

*Proof of the corollary from the theorem.* At an odd `p` dividing `(m² + n²)W`, `p` does not divide `d`, and the class of `d` is `1` exactly when `d` is a residue. At 2: if E and F hold, `d = 2d′` with `d′ ≡ ±1 (mod 8)`, which is the class `2`; conversely the class `2` means `d` even. At an odd `p` dividing `mn(m² − n²)`: the class of `d` has odd valuation exactly when `p` divides `d`, and the table allows that exactly for `p ≡ ±1 (mod 8)`, which is E; for `p ≡ 3 (mod 4)` units are trivial; for `p ≡ 5 (mod 8)` both unit classes are allowed except in the case of the law at 5, where the class must be `ε`, that is `d ≡ ±2 (mod 5)`; for `p ≡ 1 (mod 8)` everything is allowed. ∎

So what excludes a ratio beyond the criterion of `CLASS.md` is at primes of good reduction, or global. The theorem does not say that a class allowed at all bad primes is the class of a solution.

## 2. Proof of Lemma 1

In any commutative ring, if `w² = x² + y²` then

    2(w + x)·(x + iy) = (w + x + iy)²,        2(w + x)·(x − iy) = (w + x − iy)².

At a point of `X_p` each norm is a nonzero square `w_j²`, and `(w_j + x_j)(w_j − x_j) = y_j²`, so one of the signs of `w_j` gives `w_j + x_j ≠ 0` (if both vanish, `x_j = y_j = 0`, and the norm is zero). With `D = ∏ (w_j + x_j)`: `16·D·Q_f = (∏ (w_j + x_j ± i·y_j))²`. The square root of `D·Q_f` lies in `ℚ_p(i)` (or `ℚ_p × ℚ_p`) and is integral over `ℤ_p[i]`, which is integrally closed, for `p = 2` too. If also `d·Q_f` is a square, then `d/D` is a square in `ℚ_p ⊗ ℚ(i)`, because `Q_f` is invertible there (its norm is a product of nonzero norms). An element `c` of `ℚ_p^×` is a square in `ℚ_p ⊗ ℚ(i)` exactly when `c` or `−c` is a square in `ℚ_p`: in the field case `(a + bi)² ∈ ℚ_p` forces `ab = 0`; in the split case `−1` is a square. Changing the sign of one `w_j` multiplies `D` by `(x_j − w_j)/(x_j + w_j) = −y_j²/(w_j + x_j)²`, a square times `−1`. ∎

## 3. The upper bounds

That the image is contained in the set of the table is, for integer points, Theorems 1 to 5b of `CLASS.md`, which are in Lean, with one addition:

**Theorem F′ (Lean: `class_mod_sixteen`).** For coprime `m, n` and coprime `u, v` at which the four forms are squares, a `d` such that `d` times each of the eight products is a Gaussian square is `≡ 0, 2, 8` or `14 (mod 16)`; a squarefree one is `≡ ±2 (mod 16)`. The proof is a check modulo 64 by the kernel on pairs in normal form (`law16Check`: 1,024 points with the four norms squares). For the class of a solution this follows from E and F; what is new is that it holds for every `d` with the property, which is what the prime 2 needs.

For `p`-adic points the same statements hold, for these reasons; this transfer is by inspection and is not in Lean.

- *F′ and the law at 5* are finite checks on residues (modulo 64 and modulo 25). Take the representative `D₀` of the class among `1, 3, 5, 7, 2, 6, 10, 14` (for `p = 2`), reduce `(m, n, u, v)` and the square roots modulo 64: the check says `D₀ ∈ {2, 14}`, the class `2`. The same modulo 25 for the unit classes at 5.
- *The valuation at a prime `p` not dividing `2mn(m² − n²)`* (the class is a unit there). For `p ≡ 1 (mod 4)` write `α_j = (α_j′, α_j″)` in `ℤ_p × ℤ_p`. If `p` divided both components, it would divide `2x_j` and `2y_j`, hence `u` and `v` (for `α0`: `p | mv` and `p | nu`; for `αL`: `p | (m+n)(u−v)` and `p | (m−n)(u+v)`). So one component is a unit, the sum of the two valuations is even, both are even, and the class of `α_j`, which is that of either component, is a unit class. For `p ≡ 3 (mod 4)` the norm of a primitive element is a unit.
- *G, at an odd `p` dividing `m² + n²`.* Let `ι ∈ ℤ_p` be the square root of `−1` with `ι ≡ m/n (mod p)`, and use `i ↦ ι` for the first component, `i ↦ −ι` for the second. Modulo `p`:

      α0′ ≡ ιn(u+v),  α0″ ≡ ιn(v−u),  α8′ ≡ n(v−u),  α8″ ≡ n(u+v),
      αL′ ≡ −2n(1+ι)v,  αL″ ≡ 2n(1+ι)u,  αU′ ≡ 2n(ι−1)u,  αU″ ≡ −2n(ι−1)v.

  The class is that of any product of one unit component from each pair. `α0′·α8″·αL″·αU′ ≡ −8ι·n⁴(u+v)²u²` and `α0″·α8′·αL′·αU″ ≡ −8ι·n⁴(u−v)²v²`, and `−8ι = (2(1−ι))²`. One of `u(u+v)`, `v(u−v)` is a unit (if `p | u` then `v` and `u − v` are units; if `p | u + v` then `v` and `u − v = −2v + (u+v)` are). So the class is `1`. These are the identities of the Lean proof of Theorem 4.
- *E and O.* The proofs in `LAST_PAIR.md`, section 17, are congruences modulo `p` for any `d, u, v` and any square root of `d·Q_f`; we have not rewritten them over `ℤ_p`.

`class_local.py`, part 1, computes the exact class at 211,894 random points of the sets `X_p` (odd `p < 100`, random ratios with `m ≤ 400`): none is outside the table, and every class of the table occurs.

## 4. The lower bounds: every class of the table occurs

This is the part that was a conjecture. Three rows need nothing: for `p` dividing `(m² + n²)W`, for `p = 2`, and for `p ≡ 3 (mod 8)` the set has one element and `X_p ≠ ∅`. For the others, `p` is odd and divides one of `m, n, m − n, m + n`.

**Reduction to `p | m`.** Exchanging `m` and `n` exchanges `α0` and `α8` and replaces `αL`, `αU` by `ᾱL`, `−ᾱU`. The substitution `(m, n, u, v) ↦ (m + n, m − n, u + v, u − v)` sends `(α0, α8, αL, αU)` to `(αL, αU, 4α0, 4α8)`. Multiplying `(m, n)` by `λ` multiplies every product by `λ⁴`. None of this changes classes at an odd prime, so we may assume `m = p^e·m₁` with `e ≥ 1` and `p ∤ m₁n`. Below `χ` is the quadratic character modulo `p`.

**4.1 Points of first order.** Let `u, v` be units with `2(u² + v²) ≡ s² ≢ 0 (mod p)`. Modulo `p` the norms are `n²u², n²v², n²s², n²s²`: nonzero squares, so `(u, v) ∈ X_p`. The four factors `w_j + x_j` are, modulo `p`: `±nu` (as `x₀ = mv ≡ 0`); `2nv` (taking `w₈ ≡ nv`); `n(s + u − v)` and `n(s − u + v)`, whose product is `n²(u + v)²`. If `p | u + v`, take `s ≡ 2u` in the third and `−s` in the fourth: `4nu` and `−4nu`. In both cases

    cl_p(u, v) = the unit class of 2uv.

With `u = v = 1`: `X_p ≠ ∅`. For `p ≡ 3 (mod 4)` this is the class `1`.

**4.2 Both unit classes, `p ≡ 1 (mod 4)`, `p ≠ 5`.** Take `u = 1`. For `c = ±1` the number `N_c` of `v ≠ 0`, `v² ≠ −1`, with `χ(2(1 + v²)) = 1` and `χ(2v) = c` is

    N_c = ¼ · Σ_v [1 + χ(2(1+v²))]·[1 + c·χ(2v)],

the sum over those `p − 3` values. The sum of `χ(2(1 + v²))` over them is `−2χ(2)`; the sum of `χ(2v)` over them is at most 2 in absolute value; the sum of `χ(4v(1 + v²))` is the sum for the elliptic curve `y² = v³ + v`, at most `2√p` in absolute value (Hasse). So `N_c ≥ (p − 7 − 2√p)/4 > 0` for `p ≥ 17`. For `p = 13`: `(u, v) = (1, 1)` gives `χ(2) = −1`, and `(1, 2)` gives `2(1 + 4) = 10 ≡ 6²` and `χ(4) = 1`. For `p = 5` the points of first order are `v = ±u` and give only `χ(±2) = −1`: this is the law at 5. (`class_local.py`, part 3: for the 329 primes `≡ 1 (mod 4)` below 5000, only 5 has one value.)

**4.3 Classes of odd valuation, `p ≡ ±1 (mod 8)`.** Let `r² ≡ 2 (mod p)`, and take `v = 1`, `u = p·u₁` with `u₁` a unit. Then `F8 ≡ n²`, `GL ≡ GU ≡ 2n² = (nr)²`, with factors `2n`, `n(r − 1)`, `n(r + 1)`: product `2n³`, a unit of class `n` (as 2 is a square). This is where `p ≡ ±1 (mod 8)` is used; it is Theorem E read backwards.

- If `e ≥ 2`: `F0 = p²(n²u₁² + p^{2e−2}m₁²)`, a square, with `w₀ + x₀ = p·(nu₁ + …)`. So `cl_p = p·u₁`: both `p` and `εp` occur if `p ≡ 1 (mod 8)`, and `p` if `p ≡ 7 (mod 8)`.
- If `e = 1`: `F0 = p²(m₁² + n²u₁²)`. For `t` a unit with `t² ≢ ±m₁²` put `nu₁ ≡ (t − m₁²/t)/2` and `ω ≡ (t + m₁²/t)/2`; then `ω² ≡ m₁² + n²u₁²`, `ω ≢ 0`, `u₁ ≢ 0`, so `F0` is a square with `w₀ ≡ pω (mod p²)`, and `w₀ + x₀ = p(ω + m₁ + …)` with `ω + m₁ ≡ (t + m₁)²/(2t)`. So `cl_p = p·nt`. At most four values of `t` are excluded, all in one class of squares when `p ≡ 1 (mod 8)`, and `(p − 1)/2 ≥ 8` then; for `p ≡ 7 (mod 8)` any admissible `t` will do.

**4.4 The class `1` at `p = 5` when 25 divides `m`.** Here `χ(2) = χ(−2) = −1`. Take `u = 1` and `v` with `v² ≡ −1 (mod 25)`, so `v ≡ ±2 (mod 5)`; write `1 + v² = 5^b·c` with `b ≥ 2`, `5 ∤ c`. Then `F0 ≡ n²`, `F8 ≡ n²v²`, and

    GL = A − B,   GU = A + B,   A = 2(m² + n²)(1 + v²),   B = 8mnv,

with valuations `b` and `e`. If `GL` and `GU` are squares, the point is in `X_5`, and its factors are, modulo 5: `±n`, `2nv`, and `n(1 − v)`, `−n(1 − v)` (as `w_L ≡ w_U ≡ 0`, and `v ≢ 1`), so `cl_5` is the class of `2v`, and `χ(2v) = χ(2)χ(±2) = 1`: **the class `1`**. It remains to choose `v`. Replacing `v` by `v + 25s` changes `c` by `2vs` modulo 5, so `c` modulo 5 is free when `b = 2`.

- `e ≥ 3`: take `b = 2` and `c ≡ ±2`. Then `A = 25·2(m² + n²)c` is a square, because `χ(2c) = 1`, and `B/A ∈ 5ℤ_5`, so `A ∓ B = A(1 ∓ B/A)` are squares.
- `e = 2`, `χ(m₁n) = 1`: take `b ≥ 3`. Then `A/B ∈ 5ℤ_5`, `GU = B(1 + A/B)`, `GL = −B(1 − A/B)`, and `B = 25·8m₁nv` is a square, because `χ(8m₁nv) = χ(2)χ(m₁n)χ(v) = 1`; `−1` is a square.
- `e = 2`, `χ(m₁n) = −1`: take `b = 2` and `c ≡ 4m₁v/n (mod 5)`, a residue because `χ(c) = χ(m₁n)χ(v) = 1`. With `A = 25A′`, `B = 25B′`: `A′ ≡ 2n²c ≡ B′`, so `GU = 25(A′ + B′)` with `A′ + B′ ≡ 4n²c`, a square unit. For `GL = 25(A′ − B′)`: with `v = v₀ + 125σ` the number `A′ − B′` is `κ + 5λσ + 25(…)`, where `5 | κ` and `λ ≡ 4n²v₀` is a unit, so by Hensel's lemma there is `σ ∈ ℤ_5` with `A′ − B′ = 25`, and `GL = 625`.

For `e = 1` none of this is possible (`GL` or `GU` gets odd valuation), in accordance with the law at 5.

**4.5 Conclusion.** `p ≡ 7 (mod 8)`: `1` by 4.1, `p` by 4.3. `p ≡ 5 (mod 8)`, `p ≥ 13`: `1` and `ε` by 4.2. `p = 5`: `ε` by 4.1; `1` by 4.4 when 25 divides. `p ≡ 1 (mod 8)`: `1, ε` by 4.2 and `p, εp` by 4.3. With section 3 this proves Theorem C. ∎

## 5. What is in Lean

`MagicSquareOfSquares/ClassLocal.lean` (since 24 November 2026), besides the theorems for integer points (1 to 5b of `CLASS.md`, and F′):

- **Lemma 1 over the integers, both directions** (`class_iff`; `gsqP_class`, `class_unique`): at a solution on the four forms with roots `w_j` of the norms, none zero, a `d ≠ 0` has the property of the eight products if and only if `d·D` or `−d·D` is a square, `D = ∏ (w_j + x_j)`. So the class of a solution is given by a formula.
- **The same identity modulo any `M`** (`gsqM_class`): roots of the four norms modulo `M` and `d ≡ σ²·16D (mod M)` give `d·Q_f ≡ (a + ib)² (mod M)` for the eight products. This is the form in which the lower bounds are proved.
- **Hensel's lemma for squares** modulo powers of an odd prime, by induction (`exists_sq_modEq_pow`).
- **Section 4.1 for every power of `p`** (`pass_first_order`): `p` odd, `p | m`, `p ∤ n`, `u, v, u + v` units, `2(u² + v²)` a nonzero square modulo `p`. Then every `d` with `d·2uv` a nonzero square modulo `p` passes at `(u, v)` modulo `p^k` for all `k`: the four norms have square roots and `d` times each product is a Gaussian square. (The case `p | u + v` of 4.1 is not in Lean.)
- **Section 4.3, the case `p² | m`, for every power of `p`** (`pass_odd_valuation`): if 2 is a square modulo `p`, every `d = p·d₁` with `d₁·2u₁` a nonzero square modulo `p` passes at `(p·u₁, 1)`.
- **Two statements with no hypothesis on a point.** `unit_passes_of_three_mod_four`: at a prime `p ≡ 3 (mod 4)` dividing `m`, every `d` prime to `p` passes at `(1, 1)` modulo every power of `p`. `prime_mul_unit_passes_mod_eight`: at a prime `p ≡ ±1 (mod 8)` with `p² | m`, every `d = p·d₁`, `p ∤ d₁`, passes at `(2p·d₁, 1)`. Together with Theorem E for integer points these are the rows `p ≡ 3 (mod 8)` and (when `p² | m`) `p ≡ 7 (mod 8)` of the table, in the sense below.

**In what sense.** Lean's statements say "passes modulo `p^k` for every `k`" at a fixed integer point. A `p`-adic point of `X_p` with that class gives this, and conversely by compactness; the converse step is not in Lean, and Lean has no statement about `ℤ_p` here. The upper bounds in Lean are for integer solutions only.

**Not in Lean:** the conic and the bound of Hasse (4.2); the case `p ∥ m` of 4.3; the prime 5 (4.4); the reduction from `p | n`, `p | m ± n` to `p | m`; the rows for `p = 2` and for the primes of `(m² + n²)W` as statements about local points; the transfer of the upper bounds to `p`-adic points. An estimate, not a promise: 4.3 for `p ∥ m` and 4.4 are of the same kind as what is done and a few hundred lines each; 4.2 needs a bound for the character sum of `y² = v³ + v` (Jacobsthal's elementary proof would do) and is the largest piece.

## 6. What the program checks

`scripts/cover_descent/structure/class_local.py` (output `class_local.txt`; about three minutes; rerun by CI in the job `scan-5000`). Its function `cls` computes the class of Lemma 1 exactly for odd `p`: the norms at integers `u, v` are integers, a `p`-adic integer is a square when its valuation is even and its unit part a residue, and the square roots are taken to 30 digits.

1. 211,894 random points with a class (ratios with `m ≤ 400`, odd bad primes below 100): none outside the table; every class of every row is seen.
2. The points of section 4 have the classes that section 4 says: 1,728 of first order, 990 with `u = pu₁` and `e ≥ 2`, 11,208 with `e = 1`, and for nine ratios with 25 dividing `m` a point of class `1` among `v² ≡ −1 (mod 25)` (`p` up to 97 and 193, 1009, 2017; `e = 1, …, 4`).
3. The conic of 4.2 for the primes below 5000.
4. `p = 2`: the residues of `d` modulo 16 that pass somewhere modulo 64 are `0, 2, 8, 14` (the content of F′).
5. The table against the weaker test modulo `p^K` of `lean_cert.py`, for `p` dividing each of `m, n, m − n, m + n`: 180 cases, no difference.

What it does not check: the general statements. Parts 1, 2 and 5 are samples; the bound of Hasse in 4.2, Hensel's lemma in 4.4, and the transfer in section 3 are used as stated. An error in the argument that happens to agree with the samples would not be seen.
