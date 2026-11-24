# The class at the bad primes: a complete description (Theorem C)

**Status.** A proof by hand, written on 24 November 2026 by an AI model (Claude) in the session that also wrote `CLASS.md`. **No person has read it.** Part of it is in Lean (section 5 says exactly which part); the theorem as stated is not. The points used in the proof are computed again by `scripts/cover_descent/structure/class_local.py` (exact `p`-adic classes of points, not a test modulo small powers); section 6 says what that program does and does not check (it does not check the general statements). Section 7, added later the same day, treats the primes of good reduction. Experiments: nothing here is counted in the census, and nothing here shows that a magic square of squares does or does not exist.

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

**Reduction to `p | m`.** Exchanging `m` and `n` exchanges `α0` and `α8` and replaces `αL`, `αU` by `ᾱL`, `−ᾱU`. The substitution `(m, n, u, v) ↦ (m + n, m − n, u + v, u − v)` sends `(α0, α8, αL, αU)` to `(αL, αU, 4α0, 4α8)`. Multiplying `(m, n)` by `λ` multiplies every product by `λ⁴`. None of this changes classes at an odd prime, so we may assume `m = p^e·m₁` with `e ≥ 1` and `p ∤ m₁n`. (Since `fdaf0b5` the first two steps are in Lean, at the level of the eight products modulo any `M`: `gsqM_swap`, `gsqM_rot`; `Passes.rot` also checks that an odd prime does not divide both `u + v` and `u − v` when it does not divide both `u` and `v`. The factor `λ` is not needed there, since Lean's statements are for coprime `m, n`. Whether 25 divides is read off the factor that 5 divides: `Row 5 (m + n) d` and so on.) Below `χ` is the quadratic character modulo `p`.

**4.1 Points of first order.** Let `u, v` be units with `2(u² + v²) ≡ s² ≢ 0 (mod p)`. Modulo `p` the norms are `n²u², n²v², n²s², n²s²`: nonzero squares, so `(u, v) ∈ X_p`. The four factors `w_j + x_j` are, modulo `p`: `±nu` (as `x₀ = mv ≡ 0`); `2nv` (taking `w₈ ≡ nv`); `n(s + u − v)` and `n(s − u + v)`, whose product is `n²(u + v)²`. If `p | u + v`, take `s ≡ 2u` in the third and `−s` in the fourth: `4nu` and `−4nu`. In both cases

    cl_p(u, v) = the unit class of 2uv.

With `u = v = 1`: `X_p ≠ ∅`. For `p ≡ 3 (mod 4)` this is the class `1`.

**4.2 Both unit classes, `p ≡ 1 (mod 4)`, `p ≠ 5`.** Take `u = 1`. For `c = ±1` the number `N_c` of `v ≠ 0`, `v² ≠ −1`, with `χ(2(1 + v²)) = 1` and `χ(2v) = c` is

    N_c = ¼ · Σ_v [1 + χ(2(1+v²))]·[1 + c·χ(2v)],

the sum over those `p − 3` values. The sum of `χ(2(1 + v²))` over them is `−2χ(2)` (over all `v` it is `−χ(2)`; the term of `v = 0` is `χ(2)` and the two terms with `v² = −1` are 0); the sum of `χ(2v)` over them is at most 2 in absolute value; the sum of `χ(4v(1 + v²))` is the sum for the elliptic curve `y² = v³ + v`, at most `2√p` in absolute value (Hasse). So `N_c ≥ (p − 7 − 2√p)/4 > 0` for `p ≥ 17`. For `p = 13`: `(u, v) = (1, 1)` gives `χ(2) = −1`, and `(1, 2)` gives `2(1 + 4) = 10 ≡ 6²` and `χ(4) = 1`. For `p = 5` the points of first order are `v = ±u` and give only `χ(±2) = −1`: this is the law at 5. (`class_local.py`, part 3: for the 329 primes `≡ 1 (mod 4)` below 5000, only 5 has one value.) The Lean proof of this step (`ClassConic.lean`) does not use the bound of Hasse: an identity of Jacobsthal gives a weaker bound that is enough for `p ≥ 37`, and 13, 17, 29 are done by explicit points (section 5).

**4.3 Classes of odd valuation, `p ≡ ±1 (mod 8)`.** Let `r² ≡ 2 (mod p)`, and take `v = 1`, `u = p·u₁` with `u₁` a unit. Then `F8 ≡ n²`, `GL ≡ GU ≡ 2n² = (nr)²`, with factors `2n`, `n(r − 1)`, `n(r + 1)`: product `2n³`, a unit of class `n` (as 2 is a square). This is where `p ≡ ±1 (mod 8)` is used; it is Theorem E read backwards.

- If `e ≥ 2`: `F0 = p²(n²u₁² + p^{2e−2}m₁²)`, a square, with `w₀ + x₀ = p·(nu₁ + …)`. So `cl_p = p·u₁`: both `p` and `εp` occur if `p ≡ 1 (mod 8)`, and `p` if `p ≡ 7 (mod 8)`.
- If `e = 1`: `F0 = p²(m₁² + n²u₁²)`. For `t` a unit with `t² ≢ ±m₁²` put `nu₁ ≡ (t − m₁²/t)/2` and `ω ≡ (t + m₁²/t)/2`; then `ω² ≡ m₁² + n²u₁²`, `ω ≢ 0`, `u₁ ≢ 0`, so `F0` is a square with `w₀ ≡ pω (mod p²)`, and `w₀ + x₀ = p(ω + m₁ + …)` with `ω + m₁ ≡ (t + m₁)²/(2t)`. So `cl_p = p·nt`. At most four values of `t` are excluded, all in one class of squares when `p ≡ 1 (mod 8)`, and `(p − 1)/2 ≥ 8` then; for `p ≡ 7 (mod 8)` any admissible `t` will do.

**4.4 The class `1` at `p = 5` when 25 divides `m`.** Here `χ(2) = χ(−2) = −1`. Take `u = 1` and `v` with `v² ≡ −1 (mod 25)`, so `v ≡ ±2 (mod 5)`; write `1 + v² = 5^b·c` with `b ≥ 2`, `5 ∤ c`. Then `F0 ≡ n²`, `F8 ≡ n²v²`, and

    GL = A − B,   GU = A + B,   A = 2(m² + n²)(1 + v²),   B = 8mnv,

with valuations `b` and `e`. If `GL` and `GU` are squares, the point is in `X_5`, and its factors are, modulo 5: `±n`, `2nv`, and `n(1 − v)`, `−n(1 − v)` (as `w_L ≡ w_U ≡ 0`, and `v ≢ 1`), so `cl_5` is the class of `2v`, and `χ(2v) = χ(2)χ(±2) = 1`: **the class `1`**. It remains to choose `v`. Replacing `v` by `v + 25s` changes `c` by `2vs` modulo 5, so `c` modulo 5 is free when `b = 2`.

- `e ≥ 3`: take `b = 2` and `c ≡ ±2`. Then `A = 25·2(m² + n²)c` is a square, because `χ(2c) = 1`, and `B/A ∈ 5ℤ_5`, so `A ∓ B = A(1 ∓ B/A)` are squares.
- `e = 2`, `χ(m₁n) = 1`: take `b ≥ 3`. Then `A/B ∈ 5ℤ_5`, `GU = B(1 + A/B)`, `GL = −B(1 − A/B)`, and `B = 25·8m₁nv` is a square, because `χ(8m₁nv) = χ(2)χ(m₁n)χ(v) = 1`; `−1` is a square.
- `e = 2`, `χ(m₁n) = −1`: take `b = 2` and `c ≡ 4m₁v/n (mod 5)`, a residue because `χ(c) = χ(m₁n)χ(v) = 1`. With `A = 25A′`, `B = 25B′`: `A′ ≡ 2n²c ≡ B′`, so `GU = 25(A′ + B′)` with `A′ + B′ ≡ 4n²c`, a square unit. For `GL = 25(A′ − B′)`: with `v = v₀ + 125σ` the number `A′ − B′` is `κ + 5λσ + 25(…)`, where `5 | κ` and `λ ≡ 4n²v₀` is a unit, so by Hensel's lemma there is `σ ∈ ℤ_5` with `A′ − B′ = 25`, and `GL = 625`. In detail: `5 | κ` because `A′ ≡ B′ (mod 5)` by the choice of `c`; the equation is `g(σ) = (κ − 25)/5 + λσ + 5(…) = 0`, a power series in `σ` with coefficients in `ℤ_5` whose derivative is `≡ λ`, a unit, so it has a root; `v₀ + 125σ ≡ v₀ (mod 125)` leaves `b = 2` and `c` modulo 5 as they were, so `GU = 25(A′ + B′)` with `A′ + B′ ≡ 4n²c` is still a square; and no norm is zero (`F0 ≡ n²`, `F8 ≡ n²v²`, `GL = 625`, `GU = 25·unit`).

For `e = 1` none of this is possible (`GL` or `GU` gets odd valuation), in accordance with the law at 5.

Since `fdaf0b5` the three cases are in Lean with explicit points (`ClassFive.lean`: `(1, 7)`; `(1, 57)`; and in the third case a root `v ≡ 3 (mod 5)` of `GL(1, v) ≡ 625` modulo each power of 5, by the quadratic formula), as statements modulo every power of 5; `7² + 1 = 50` and `57² + 1 = 3250 = 2·5³·13` are the choices `b = 2` and `b = 3` above.

**4.5 Conclusion.** `p ≡ 7 (mod 8)`: `1` by 4.1, `p` by 4.3. `p ≡ 5 (mod 8)`, `p ≥ 13`: `1` and `ε` by 4.2. `p = 5`: `ε` by 4.1; `1` by 4.4 when 25 divides. `p ≡ 1 (mod 8)`: `1, ε` by 4.2 and `p, εp` by 4.3. With section 3 this proves Theorem C. ∎

## 5. What is in Lean

`MagicSquareOfSquares/ClassLocal.lean` (since 24 November 2026), besides the theorems for integer points (1 to 5b of `CLASS.md`, and F′):

- **Lemma 1 over the integers, both directions** (`class_iff`; `gsqP_class`, `class_unique`): at a solution on the four forms with roots `w_j` of the norms, none zero, a `d ≠ 0` has the property of the eight products if and only if `d·D` or `−d·D` is a square, `D = ∏ (w_j + x_j)`. So the class of a solution is given by a formula.
- **The same identity modulo any `M`** (`gsqM_class`): roots of the four norms modulo `M` and `d ≡ σ²·16D (mod M)` give `d·Q_f ≡ (a + ib)² (mod M)` for the eight products. This is the form in which the lower bounds are proved.
- **Hensel's lemma for squares** modulo powers of an odd prime, by induction (`exists_sq_modEq_pow`).
- **Section 4.1 for every power of `p`** (`pass_first_order`): `p` odd, `p | m`, `p ∤ n`, `u, v, u + v` units, `2(u² + v²)` a nonzero square modulo `p`. Then every `d` with `d·2uv` a nonzero square modulo `p` passes at `(u, v)` modulo `p^k` for all `k`: the four norms have square roots and `d` times each product is a Gaussian square. (The case `p | u + v` of 4.1 is not in Lean.)
- **Section 4.3, the case `p² | m`, for every power of `p`** (`pass_odd_valuation`): if 2 is a square modulo `p`, every `d = p·d₁` with `d₁·2u₁` a nonzero square modulo `p` passes at `(p·u₁, 1)`.
- **Section 4.2 without the bound of Hasse** (`MagicSquareOfSquares/ClassConic.lean`). An elementary identity of Jacobsthal replaces it: with `S(k) = ∑_v χ(v(k + v²))`, `∑_k S(k)² = 2p(p − 1)` when `−1` is a square (`jac_sum_sq`), and `S(t²)² = S(1)²` (`jac_sq`), so `3·S(1)² ≤ 2p(p − 1)` from the squares 1, 4, 9; that is enough for the count of 4.2 when `p ≥ 37` (`exists_conic_point`), and 13, 17, 29 are done by explicit points. Result, `unit_passes_of_one_mod_four_ne_five`: **at every prime `p ≡ 1 (mod 4)` dividing `m`, except 5, every `d` prime to `p` passes at some point `(1, v)` modulo every power of `p`.** With the next item: the unit classes pass at every odd prime of `m` except 5, where the law at 5 forbids one of them.
- **Section 4.3 for `m = p·m₁`, and the prime 5** (`MagicSquareOfSquares/ClassOdd.lean`). `pass_odd_valuation_gen`: at `(p·u₁, 1)` with `ω₀² ≡ m₁² + n²u₁²`, the classes `p·d₁` with `d₁·2n(ω₀ + m₁)` a nonzero square modulo `p` pass (this is the computation of 4.3 for `e = 1`, and it holds for any `m₁`). `prime_mul_passes_at_zero_one`: the point `(0, 1)`. `prime_mul_unit_passes_seven_mod_eight`: at `p ≡ 7 (mod 8)` with `p ∥ m`, every `p·d₁` passes. `passes_at_five`: at `p = 5` dividing `m`, every `d ≡ ±2 (mod 5)` passes at `(1, 1)`. `prime_mul_unit_passes_one_mod_eight`: at `p ≡ 1 (mod 8)` with `p ∥ m`, every `p·d₁` passes: with `λ = d₁nm₁`, at `(0, 1)` if `λ` is a square, and otherwise at the point of 4.3 with `t = λ` (for a non-residue `λ` no value of the parameter has to be avoided, because `±√−1` are squares when `p ≡ 1 (mod 8)`).
- **Two statements with no hypothesis on a point.** `unit_passes_of_three_mod_four`: at a prime `p ≡ 3 (mod 4)` dividing `m`, every `d` prime to `p` passes at `(1, 1)` modulo every power of `p`. `prime_mul_unit_passes_mod_eight`: at a prime `p ≡ ±1 (mod 8)` with `p² | m`, every `d = p·d₁`, `p ∤ d₁`, passes at `(2p·d₁, 1)`. Together with Theorem E for integer points these are the rows `p ≡ 3 (mod 8)` and (when `p² | m`) `p ≡ 7 (mod 8)` of the table, in the sense below.
- **Section 4.4, the class `1` at `p = 5` when 25 divides `m`** (`MagicSquareOfSquares/ClassFive.lean`), in three cases, each modulo every power of 5: `class_one_passes_at_five_of_125` (125 divides `m`: the point `(1, 7)`), `class_one_passes_at_five_of_25` (`m = 25m₁` with `m₁n` a square modulo 5: the point `(1, 57)`), `class_one_passes_at_five_of_25_nonres` (`m₁n` not a square modulo 5: a point `(1, v)` with `v ≡ 3 (mod 5)`, where `v` is a root modulo `5^k` of `GL(1, v) ≡ 625`, found by the quadratic formula and Hensel's lemma). The text of 4.4 takes `v² ≡ −1 (mod 25)` and leaves the choice to Hensel's lemma; Lean's points are explicit, and they are points of that kind.
- **The whole table for an odd prime of `mn(m² − n²)`, in one statement** (`MagicSquareOfSquares/ClassTable.lean`). `Passes p k d m n`: some `(u, v)`, not both divisible by `p`, at which `d` times each of the eight products is a Gaussian square modulo `p^(k+1)`. `Row p t d`: the condition of the table on `d` at a prime `p` dividing the factor `t` (`d` prime to `p`, and if `p = 5` then `d ≡ ±2 (mod 5)` or 25 divides `t`; or `d = p·d₁` with `p ∤ d₁` and `p ≡ ±1 (mod 8)`). `passes_of_row`: for `p | m`, every class of the row passes modulo every power of `p`. **The reduction of section 2 is in Lean too**: `gsqM_swap` (exchanging `m` and `n` permutes the eight products up to sign and conjugation), `gsqM_rot` (the substitution `(m + n, m − n, u + v, u − v)` permutes them up to the factor 16 and conjugation), `Passes.swap`, `Passes.rot`. `passes_of_table`: **for coprime `m, n` and an odd prime `p` dividing one of `m`, `n`, `m + n`, `m − n`, every class of the row of that factor passes modulo every power of `p`.**

**In what sense.** Lean's statements say "passes modulo `p^k` for every `k`" at an integer point; in some of them the point is the same for all `k`, in others (the conic of 4.2, the third case at 5, and so `Passes`) it may depend on `k`. A `p`-adic point of `X_p` with that class gives this. Conversely a limit of such points, by compactness, is a `p`-adic point at which `d` times each product is a Gaussian square, and it has the class `d` if none of the four forms vanishes there. `Passes` does not record that last condition (in the proofs the four norms have roots of fixed valuation, so it holds for the points that are used, but that is not what the Lean statement says). The converse step is not in Lean, and Lean has no statement about `ℤ_p` here. The upper bounds in Lean are for integer solutions only.

**The table, row by row, as lower bounds in Lean** (modulo every power of `p`, no hypothesis on a point, for `p` dividing any of `m`, `n`, `m ± n`): `p ≡ 3 (mod 8)`: complete. `p ≡ 7 (mod 8)`: complete. `p ≡ 5 (mod 8)`: complete, the prime 5 included. `p ≡ 1 (mod 8)`: complete. So **every lower bound of section 4 at an odd prime of `mn(m² − n²)` is in Lean**, in the sense above (`passes_of_table`). One remark on `Row`: it is stated for `d` prime to `p` or exactly divisible by `p`; a `d` divisible by `p²` has the class of `d/p²`, and that step is not a Lean statement.

**Not in Lean:** every statement about `ℤ_p` or `ℚ_p` (that a class which passes modulo every power of `p` is the class of a `p`-adic point; Lemma 1 over `ℤ_p`; the upper bounds E, F′, G, O and the law at 5 for `p`-adic points, which Lean has for integer solutions only); the rows for `p = 2` and for the primes of `(m² + n²)W` as statements about local points (their lower bound is that some local point exists, which is a hypothesis of the theorem, not a claim); the case `p | u + v` of 4.1, which the proof does not need. So Theorem C as stated, a statement about `p`-adic points, is not a Lean theorem; what Lean has is its lower-bound half at the odd primes of `mn(m² − n²)` in the form "modulo every power of `p`", and its upper-bound half for integer solutions. The `p`-adic statements are a different kind of work (Mathlib's `ℤ_p`), not started.

## 6. What the program checks

`scripts/cover_descent/structure/class_local.py` (output `class_local.txt`; about three minutes; rerun by CI in the job `scan-5000`). Its function `cls` computes the class of Lemma 1 exactly for odd `p`: the norms at integers `u, v` are integers, a `p`-adic integer is a square when its valuation is even and its unit part a residue, and the square roots are taken to 30 digits.

1. 211,894 random points with a class (ratios with `m ≤ 400`, odd bad primes below 100): none outside the table; every class of every row is seen.
2. The points of section 4 have the classes that section 4 says: 1,728 of first order, 990 with `u = pu₁` and `e ≥ 2`, 11,208 with `e = 1`, and for nine ratios with 25 dividing `m` a point of class `1` among `v² ≡ −1 (mod 25)` (`p` up to 97 and 193, 1009, 2017; `e = 1, …, 4`).
3. The conic of 4.2 for the primes below 5000.
4. `p = 2`: the residues of `d` modulo 16 that pass somewhere modulo 64 are `0, 2, 8, 14` (the content of F′).
5. The table against the weaker test modulo `p^K` of `lean_cert.py`, for `p` dividing each of `m, n, m − n, m + n`: 180 cases, no difference.

What it does not check: the general statements. Parts 1, 2 and 5 are samples; the bound of Hasse in 4.2 (not used by the Lean proof), Hensel's lemma in 4.4 (not used by the Lean proof either), and the transfer in section 3 are used as stated. An error in the argument that happens to agree with the samples would not be seen.

## 7. The primes of good reduction, and a criterion for every ratio

Added on 24 November 2026, after `fdaf0b5`. **A finite computation with a short argument, written by an AI model; no person has read it; nothing of it is in Lean.** It is not counted in the census.

Theorem C is about the bad primes. At a prime `p` that does not divide `2mn(m² − n²)(m² + n²)W` the class of a point is a unit, so for `p ≡ 3 (mod 4)` there is one class, and for `p ≡ 1 (mod 4)` two, `1` and `ε`. Section 12 of `LAST_PAIR.md` observed that classes "die at small primes that change from ratio to ratio". Those primes are of this kind, and what happens there depends only on `m/n` modulo `p`, so it can be tabulated once.

**Theorem C′ (the class at the good primes).** Let `p` be a prime of good reduction for `m : n`.

1. If `p ≡ 1 (mod 4)`: both classes `1` and `ε` occur in `X_p`, unless `m/n` modulo `p` is in a row of `scripts/cover_descent/structure/class_good_table.txt`. The table has rows at 40 primes `≡ 1 (mod 4)`, from 17 to 617: in a row `R` every point has the class `1` (484 ratios modulo `p` in all), in a row `N` every point has the class `ε` (332), in a row `E` there is no point (36, at 17, 41, 89, 137).
2. If `p ≡ 3 (mod 4)`: `X_p` is not empty, unless `m/n` modulo `p` is in a row `E` of the table at one of 7, 23, 31, 47, 79, 103, 127, 151 (84 ratios modulo `p` in all).
3. In particular a good prime above 617 gives no condition on the class, for any ratio.

*Proof.* At a good prime the eight zeros of the four forms are distinct modulo `p` (`LAST_PAIR.md`, section 9), so the curve `w_j² = N_j(u, v)`, `j = 0, 8, L, U`, has good reduction; it is a cover of the line with group `(ℤ/2)⁴` branched at 8 points, of genus 17. Its points over `𝔽_p` lift to `X_p` and every point of `X_p` reduces to one. By the identity of Lemma 1, `D = ∏ (w_j + x_j)` is, modulo squares, `∏ (x_j + i·y_j)` up to a constant: a product of four linear forms, one zero of each `N_j`; this is not a square times a product of the `N_j`, so `z² = D` is an unramified double cover that is geometrically irreducible, of genus `2·17 − 1 = 33`. A point of the curve has the class `1` when it lifts to this cover over `𝔽_p` and the class `ε` when it lifts to its twist. By the bound of Weil each of the two has a point over `𝔽_p` when `p + 1 > 66√p`, that is for `p > 4354`, and the curve of genus 17 has one when `p + 1 > 34√p`, that is for `p > 1154`. Below these bounds `class_good.py` goes through every ratio modulo `p` and every point: the class of a point modulo `p` is `χ(∏ f_j)` with `f_j = w_j + x_j`, or `2x_j` where `y_j = 0`, or `x_j` where `N_j = 0` (at a `p`-adic point near a zero of `N_j` the root `w_j` is divisible by `p`, and `x_j` is a unit there). ∎

What was checked (`class_good.txt`): the table against the test modulo `p²` of `lean_cert.py` for the five primes up to 41 (132 cases, points near the zeros included, no difference), and against exact `p`-adic classes at integer points, computed by the function `cls` of `class_local.py`, for integer ratios at all 40 primes (196 cases, no difference). What was not: nobody has read the argument; the program is the only implementation of the enumeration above 41.

**Where there is no local point at a bad prime.** Theorem C gives the classes at a bad prime if there is a point. Computed (`class_selmer.py`, and the two short programs quoted in `class_selmer_notes.txt`):

- `p = 2`: there is a point modulo 64 exactly when `4 | mn` or `m ≡ ±n (mod 8)` (all classes modulo 64; a sample modulo 128 agrees). That there is no point otherwise is certain; that a point modulo 64 lifts is not proved here.
- An odd prime `q` of `m² + n²`: modulo `q` the forms are `m²(v² − u²)`, its negative, and `∓8mn·uv`, and `2mn` is a square, so a point of first order is a `t = u/v` with `t` and `1 − t²` nonzero squares; their number is at least `(q − 7 − 2√q)/4`, positive for `q ≥ 17` (the count of 4.2). For `q = 5` and `q = 13` there is none, and **if 5 or 13 divides `m² + n²` exactly once there is no point at all** (every class of ratios modulo 125, modulo 169 and beyond). If 25 or 169 divides, the tests modulo `5⁵` and `13³` have points; that case is not decided here.
- A prime `q` of `W`: a point of first order exists for every root, for the 792 pairs (prime below 3000, root). Not proved for larger primes.
- An odd prime of `mn(m² − n²)`: there is always a point (section 4).

**The criterion.** For coprime `m > n > 0` let `S(m, n)` be the set of squarefree `d > 0` that are the class of a point of `X_p` at every prime `p`. Putting Theorem C, Theorem C′ and the list above together: `S(m, n)` is empty if there is no point at 2, at 5 or 13, or at a good prime (rows `E`); otherwise it is the set of

    d = 2·∏_{p ∈ T} p,   T a set of primes ≡ ±1 (mod 8) dividing mn(m² − n²),

such that `d` is a square modulo every odd prime of `(m² + n²)W`, `d ≡ ±2 (mod 5)` if 5 divides `mn(m² − n²)` exactly once, and `d` is a residue (a non-residue) modulo `p` when `m/n` is in a row `R` (`N`) at a good prime `p`. These are linear equations over `𝔽₂` in the exponents of the primes of `T`, so `S(m, n)` is empty or a coset of a subgroup, and it is computed by factoring three numbers and looking up 48 residues: no local point is searched for. A solution on the four forms has its class in `S(m, n)`; the converse is not claimed (a class in `S` need not come from a rational point).

*Status.* This rests on Theorem C (by hand, unread), on Theorem C′ (above), and on the list of bad primes without a point, where two cases are open: `q² | m² + n²` for `q = 5, 13`, where the program falls back on the test modulo `5³`, `13²`; and primes of `W` above 3000, where it assumes a point. So it is a criterion with stated gaps, and the tests below are what supports it as a whole.

*Tests* (`class_selmer.py`; outputs `class_selmer.txt`, `class_selmer_lists.txt`, `class_selmer_all_5000.txt`).

1. The 425 ratios of `exceptions1000`, with every squarefree product of the primes of `2mn(m² − n²)` as a candidate: the criterion and the local test of `lean_cert.py` at every prime below 620 leave the same set **for all 425**. Of the 119 ratios with `m² + n²` neither a square nor twice a square, 115 have no class and 4 have one (these are the four that needed four classes in section 15 of `LAST_PAIR.md`); of the 306 others, 300 have one class and 6 have two.
2. The 21,104 ratios with `m ≤ 5000` of section 12 of `LAST_PAIR.md`. Of the 19,581 that are not of the two types, 19,478 have no class (13,820 of them because there is no point at 2, 5 or 13: the 484 classes of that list are about another curve), and 103 have a class; the local test below 620 gives the same for these 103 and for 1,500 drawn from the others. Of the 1,523 of the two types, 1,482 have exactly the class of their degenerate solutions, 40 have two classes and one has four; a class of a degenerate solution is never missing, and the local test below 620 gives the same set for all 1,523.
3. **Every coprime ratio with `m ≤ 5000`** (7,600,457; no list is used; about four minutes): 7,530,318 have no local point at some prime; 67,619 have points everywhere but no class; 1,523 are of the two types, and each of them has a class; **997 are not of the two types and have a class**. In the four quarters of the range the last number is 60, 191, 322, 424, that is 1.5, 1.5, 1.5 and 1.4 per cent of the ratios that are not of the two types and have points everywhere: the proportion does not go to zero in this range, which agrees with section 17 of `LAST_PAIR.md` (no statement of density 1 for one class). The first five are (291,224), (292,133), (308,185), (340,53), (393,287).
4. For these 997 ratios (`class_selmer_997.txt`) the descent with four classes of `LAST_PAIR.md`, section 12 (`tuple_descent.py`, good primes up to 150) leaves no tuple, for every one (`scripts/cover_descent/tuple_selmer_997.txt`). With test 3 this says: **for every coprime ratio with `m ≤ 5000`, if `m² + n²` is neither a square nor twice a square, then `F0, F8, GL, GU` are not all squares at any coprime `(u, v)`**, provided these programs are right. Until now the repository had this for the 19,581 ratios that avoid the 484 classes, and those classes come from the fifth form of the family; here the fifth form is not used, and every ratio is covered. What the exclusion of a ratio uses: either that there is no point at 2, at 5, at 13 or at a good prime (finite computations); or the laws E, F, G, O and at 5, which are Lean theorems for integer solutions, with the rows `R`, `N` of the table (a computation); or, for 997 ratios, the program with four classes. It does not use the lower bounds of Theorem C, which matter only for the claim that the criterion loses nothing.

*What this is good for.* It turns the first descent into arithmetic of residue symbols for every ratio at once, so the question "which ratios have a class" can be asked on millions of ratios, and questions about the group of which `S` is a coset (its rank, how it varies in the two families of degenerate type) become questions about a matrix of Legendre symbols. *What it is not:* a result on the magic square. The census is not changed, and for the ratios with `m ≤ 1000` it says again what the certificates in Lean already prove.
