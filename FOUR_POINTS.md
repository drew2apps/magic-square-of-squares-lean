# Four points on one elliptic curve: ratios of rank 2 and 3

A self-contained account of the argument behind `scripts/g1_full/dm4.py`, written so that a reader who knows elliptic curves over ℚ can check it without the rest of the repository. It was written with Claude (an AI model made by Anthropic). **Nobody else has read it.** The closures it produces are listed separately from the others everywhere in this repository until somebody has.

What is machine-checked and what is not is stated in each section and collected in Section 9.

## 1. Setting

Fix coprime integers `m > n ≥ 1` and put `S = m² + n²`. A 3×3 magic square of squares whose edge roots satisfy `(b + d) : (h + f) = m : n` gives a rational number `t0` at which five quadratic forms are squares (`G1_FULL_SYSTEM.md`, "The system"):

`F0 = m² + n²t²`, `F8 = n² + m²t²`, `G4 = 2S(t² + 1)`, `GL = G4 − 8mnt`, `GU = G4 + 8mnt`.

Degenerate squares correspond to `t0 ∈ {0, ±1}`. The claim to be proved for a given pair `(m, n)` is: **no rational `t0` other than these makes the five forms squares.**

Let `C: y² = q(t) = F0(t)·GL(t)`. It has the rational point `t1 = m/n`, `y1 = 2m(m² − n²)/n`. Three substitutions map the other products to the same curve:

`q(−t) = F0·GU`, `t⁴q(1/t) = F8·GL`, `t⁴q(−1/t) = F8·GU`.

So a solution `t0 ≠ 0` gives **four rational points of `C`**, with `t`-coordinates

`Q1: t0`, `Q2: −t0`, `Q3: 1/t0`, `Q4: −1/t0`.

(Each `Q_i` is one of the two points above its `t`; the choice does not matter below.)

`E` denotes `C` as an elliptic curve with origin `O = (t1, y1)`. The involution `y ↦ −y` is `Q ↦ R0 − Q`, where `R0 = (t1, −y1)`: the fibres of `t` are the divisors of a pencil of degree 2 that contains `O + R0`. Canonical heights are normalised as in Sage, `ĥ(P) = lim h(x(2ᵏP))/4ᵏ` with `h(u/v) = log max(|u|, |v|)`; this is twice Silverman's.

Notation: `a_i = Q_i − R0/2` in `E(ℚ) ⊗ ℚ`; `N(v) = ĥ(v)`, extended to a positive definite quadratic form on `E(ℚ) ⊗ ℝ`, and `⟨ , ⟩` its bilinear form, so `N(a_i) = ĥ(2Q_i − R0)/4`. `H = h(t0)`.

## 2. The model and three identities

With Connell's formulas for the quartic and its point `(t1, y1)`, `E` is the curve `E0` with `a1 = 6m² − 2n²`, `a2 = 5m⁴ − 2m²n² + n⁴`, `a3 = 32m⁴(m² − n²)`, `a4 = −32m²(m² − n²)²(m² + n²)`, `a6 = a2·a4`. Its three points of order 2 are rational. With `x = X + 6m⁴ + 2n⁴` the curve is

`E_T: Y² = x(x − 2A)(x − 2B)`, `A = m⁴ − 2m³n + 2m²n² + 2mn³ + n⁴`, `B = m⁴ + 2m³n + 2m²n² − 2mn³ + n⁴`, `A·B = m⁸ + 14m⁴n⁴ + n⁸`.

(The octic is the one whose prime divisors give the singular fibres in Result 1 of `G1_FULL_SYSTEM.md`.)

Let `y1² = F0·GL`, `y2² = F0·GU`, `y3² = F8·GL`, and `Q1 = (t, y1)`, `Q2 = (−t, y2)`, `Q3 = (1/t, y3/t²)`. Then, with `Π(t) = 2mn(n²t² − m²) + (m⁴ − n⁴)t`:

- **(I0)** `x(2Q − R0) = 4Π(t)²/(F0·GL)` for `Q = (t, y)` on `C`;
- **(I±)** `x(Q1 + Q2 − R0) = ((y1 − y2)/(2t))²` and `x(Q1 − Q2) = ((y1 + y2)/(2t))²`;
- **(Iinv)** `x(Q1 + Q3 − R0) = ((y1 − y3)/(t² − 1))²` and `x(Q1 − Q3) = ((y1 + y3)/(t² − 1))²`.

*How checked.* `scripts/g1_full/dm4_identities.py` verifies the model and the five equalities by exact computer algebra in the fields `ℚ(m, n)(t)(y1, y2)` and `ℚ(m, n)(t)(y1, y3)`, with `m`, `n` as indeterminates (output in `dm4_identities.txt`). `dm4.py` checks them again for each pair it is run on. They are not in Lean. As a test on rational points: for a degenerate solution `t0 = 1` (possible when `S` is a square) the consequence `x(Q1 + Q2 − R0)·x(Q1 − Q2) = (4mnS)²` holds exactly (`dm4_degenerate.py`).

*Remarks, not used below.* For a magic square with roots `a b c / d e f / g h i` and `t0 = u/v` the right sides are `(4a(c ∓ g)/(uv))²` and `(8c(a ∓ i)/(u² − v²))²`, up to the signs of the roots. And since `x` is a square exactly on the image of the 2-isogeny whose dual has kernel `{O, (0, 0)}`, the four points of a solution lie in one class of `E(ℚ)` modulo that image.

## 3. Heights of the three fractions

Let `t0 = p/q` in lowest terms, `M = max(p², q²)`, so `H = ½ log M`. Put

`Ψ0 = m²q² + n²p²`, `Ψ8 = n²q² + m²p²`, `Ψ = S(p² + q²) − 4mnpq`, `Π = 2mn(n²p² − m²q²) + (m⁴ − n⁴)pq`.

These are `q²F0`, `q²F8`, `q²GL/2`, `q²Π` at `t0`. For a rational number `ξ = u/v` in lowest terms write `H(ξ) = max(|u|, |v|)`.

**Lemma A (from (I±)).** Let `pq ≠ 0`, `G = 2S(p² + q²)`, and let `r` be an integer with `r² = G² − 64m²n²p²q²`. Then `ξ = 32m²n²Ψ0/(G + r)` satisfies `M/(2m²) ≤ H(ξ) ≤ 32m²n²S·M`.

*Proof.* `|r| < G`, so `0 < G + r < 2G ≤ 8SM`, and `n²M ≤ Ψ0 ≤ SM`. Hence numerator and denominator are positive and at most `32m²n²S·M`. Since `(G + r)(G − r) = 64m²n²p²q²`, the denominator divides `64m²n²p²q²`, and a common divisor of `32m²n²Ψ0` and `64m²n²p²q²` divides `64m⁴n⁴` (Lean: `G1.fourPoint_pm_dvd`). So `H(ξ) ≥ 32m²n²·n²M/(64m⁴n⁴) = M/(2m²)` (Lean: `G1.fourPoint_pm_lower`). ∎

That the two squares of (I±) are the fractions `32m²n²Ψ0/(G ± r)`, with `r = Y1Y2/Ψ0` for `Y_i = q²y_i`, is `G1.fourPoint_pm_fraction`.

**Lemma B (from (Iinv)).** Let `p² ≠ q²`, `G′ = S(p² + q²)`, and let `r′` be an integer with `r′² = G′² − (m² − n²)²(p² − q²)²`. Then `ξ = 2(m² − n²)²Ψ/(G′ + r′)` satisfies `M/(8(m + n)²) ≤ H(ξ) ≤ 4(m² − n²)²(m + n)²·M`.

*Proof.* `|4mnpq| ≤ 2mn(p² + q²)` gives `(m − n)²M ≤ Ψ ≤ 2(m + n)²M`. `0 < G′ + r′ < 2G′ ≤ 4SM`. So numerator and denominator are positive and at most `4(m² − n²)²(m + n)²M` (for the denominator because `S ≤ (m² − n²)²(m + n)²`). The denominator divides `(m² − n²)²(p² − q²)²`, and a common divisor of that and of `2(m² − n²)²Ψ` divides `16(m² − n²)⁴` (Lean: `G1.fourPoint_inv_dvd`). So `H(ξ) ≥ 2(m² − n²)²(m − n)²M/(16(m² − n²)⁴) = M/(8(m + n)²)` (Lean: `G1.fourPoint_inv_lower`). ∎

The corresponding statement about the squares is `G1.fourPoint_inv_fraction`.

**Lemma C (from (I0)).** `ξ = 2Π²/(Ψ0Ψ)` satisfies `n²(m − n)²M²/(2R0²R1²) ≤ H(ξ) ≤ 2(2m³n + m⁴ − n⁴)²M²`, where `R0 = m²n²(m⁸ + 14m⁴n⁴ + n⁸)` and `R1 = (m² − n²)²(m⁸ + 14m⁴n⁴ + n⁸)`.

*Proof.* `|Π| ≤ (2m³n + m⁴ − n⁴)M` and `n²(m − n)²M² ≤ Ψ0Ψ ≤ 2S(m + n)²M²`, which is at most `2(m⁴ − n⁴)²M²` because `m⁴ − n⁴ ≥ S(m + n)`; this gives the upper bound. For every prime `ℓ`, `v_ℓ(gcd(2Π², Ψ0Ψ)) ≤ v_ℓ(2) + 2 min(v_ℓ(Π), v_ℓ(Ψ0)) + 2 min(v_ℓ(Π), v_ℓ(Ψ))`, so the common factor divides `2·gcd(Π, Ψ0)²·gcd(Π, Ψ)²`. For coprime `p`, `q`, `gcd(Π, Ψ0)` divides `R0` and `gcd(Π, Ψ)` divides `R1` (Lean: `G1.fourPoint_res_psi0_dvd`, `G1.fourPoint_res_psi_dvd`). So a common divisor of `2Π²` and `Ψ0Ψ` divides `2R0²R1²` (Lean: `G1.fourPoint_lemmaC_dvd`, with the step on common factors as `G1.gcd_two_sq_mul_dvd`), and `H(ξ) ≥ Ψ0Ψ/(2R0²R1²)`. ∎

In Lean, then: every divisibility statement of this section and the lower bounds of Lemmas A and B. Not in Lean: the upper bounds (they only use the sizes of numerators and denominators stated in the proofs), and the last line of Lemma C, which divides the lower bound for `Ψ0Ψ` by `2R0²R1²`.

## 4. From naive to canonical heights

For every rational point `P ≠ O` of an integral Weierstrass model, `−U ≤ h(x(P)) − ĥ(P) ≤ BUP`, with

`BUP = Σ_p (α_p + ord_p(Δ/Δ_min)/6) log p + ¼ log⁺|j| + log⁺|b2/12| + log 2 + 1.946`, `U = ⅙ log|Δ| + ⅙ log⁺|j| + log⁺|b2/12| + log 2 + 2.14`

(Cremona, Prickett and Siksek 2006: equation (1), Proposition 8 and Table 1 for the finite places, `α_p` depending on the Kodaira type and the Tamagawa number; their Theorem 15, quoting Silverman 1990, for the real place; the `log 2` is present when `b2 ≠ 0`). These statements were read from the paper with an AI model, and Silverman's paper has not been seen. `hbounds.py` evaluates both constants for `E_T` in rational arithmetic; the reduction data come from Sage.

With Lemmas A to C, for a solution `t0 = p/q` with `pq ≠ 0` and `p² ≠ q²`:

**Proposition 1.** Put `lo4 = log(n²(m − n)²) − log(2R0²R1²)`, `up4 = log(2(2m³n + m⁴ − n⁴)²)`, `w = max(log(64m⁴n²S), log(32(m² − n²)²(m + n)⁴))`, and

`eL = (BUP − lo4)/4`, `eU = (up4 + U)/4`, `κ = (BUP + U + w)/4`.

Then (i) `H − eL ≤ N(a_i) ≤ H + eU` for `i = 1, …, 4`; (ii) `|⟨a_1, a_2⟩|`, `|⟨a_3, a_4⟩|`, `|⟨a_1, a_3⟩|`, `|⟨a_2, a_4⟩|` are at most `κ`.

*Proof.* (i) `N(a_i) = ĥ(2Q_i − R0)/4`, and by (I0) and Lemma C `4H + lo4 ≤ h(x(2Q_i − R0)) ≤ 4H + up4`, for each of the four values `±t0`, `±1/t0`, which have the same `M`. (ii) `4⟨a_1, a_2⟩ = ĥ(a_1 + a_2) − ĥ(a_1 − a_2) = ĥ(Q1 + Q2 − R0) − ĥ(Q1 − Q2)`. By (I±) and Lemma A both naive heights lie in `[2H − log(2m²), 2H + log(32m²n²S)]`, so the canonical heights differ by at most `BUP + U + log(64m⁴n²S)`. The pair `(a_3, a_4)` is the same with `1/t0` in place of `t0`. For `(a_1, a_3)` use (Iinv) and Lemma B at `t0`, and for `(a_2, a_4)` at `−t0`. None of the points whose `x` is used is `O`, because `t0 ∉ {0, ±1}`. ∎

So the four points have almost the same height `h(t0)`, and four of the six pairs are almost orthogonal. The pairs `(Q1, Q4)` and `(Q2, Q3)` are not controlled: their sums and differences are defined over the function field of the curve `w² = F0·F8·GL·GU`, of genus 3.

## 5. The lattice step

Let `L = H − eL`, `e = eL + eU`, and write `N_i = N(a_i)`, `x = ⟨a_1, a_4⟩`, `y = ⟨a_2, a_3⟩`. The number `N_1 + N_4 − 2|x|` is `N(a_1 − a_4)` or `N(a_1 + a_4)`, that is, the canonical height of the rational point `z = Q1 − Q4` or `z = Q1 + Q4 − R0`. Likewise for `N_2 + N_3 − 2|y|` and `Q2 ∓ Q3`.

*These combinations are rational points, not half-points.* Each `a_i` lies in `E(ℚ) ⊗ ℚ` (in general not in `E(ℚ)`, because of the term `R0/2`), but `a_1 − a_4 = Q1 − Q4` and `a_1 + a_4 = Q1 + Q4 − R0` with `R0 ∈ E(ℚ)`: both signs give elements of `E(ℚ)`, with coefficients ±1. Only such combinations are used, in rank 2 and in rank 3, so the bound on `ĥ(z)` is a bound for a point of `E(ℚ)`, and the list of all points of `E(ℚ)` up to that height (torsion included, and `z = O`) is the right list. No vector is divided by anything: `N_i ≥ L > 0`, so no `a_i` is zero, and the propositions make no other assumption on the four vectors than the rank.

**Proposition 2 (rank at most 2).** If `E(ℚ)` has rank at most 2 and `L > 0`, then `ĥ(z) ≤ (e²/2 + 4κ²)/L` for `z = Q1 − Q4` or `z = Q1 + Q4 − R0`.

*Proof.* `a_1, a_2, a_4` are linearly dependent, so their Gram matrix has a nonzero kernel vector and its determinant vanishes: `N_2(N_1N_4 − x²) = N_1⟨a_2, a_4⟩² + N_4⟨a_1, a_2⟩² − 2⟨a_1, a_2⟩⟨a_2, a_4⟩x ≤ κ²(N_1 + N_4 + 2|x|)`. Now `(N_1 + N_4 − 2|x|)(N_1 + N_4 + 2|x|) = (N_1 − N_4)² + 4(N_1N_4 − x²) ≤ e² + 4κ²(N_1 + N_4 + 2|x|)/N_2`, and `N_1 + N_4 + 2|x| ≥ 2L`, `N_2 ≥ L`. ∎ (Lean: `G1.fourPoint_rank2_bound`, as a statement about the entries of a 3×3 symmetric matrix with a kernel vector.)

**Proposition 3 (rank 3).** If `E(ℚ)` has rank at most 3 and `L > 2κ`, then `ĥ(z) ≤ 4κ + (e² + 16κe)/(4(L − κ))` for one of `z = Q1 − Q4`, `Q1 + Q4 − R0`, `Q2 − Q3`, `Q2 + Q3 − R0`.

*Proof.* The Gram matrix of `a_1, a_4, a_2, a_3` is singular; write it in blocks `[[A, B], [Bᵀ, D]]`, `A = [[N_1, x], [x, N_4]]`, `D = [[N_2, y], [y, N_3]]`. The entries of `B` are the four products of Proposition 1. Let `(u, v) ≠ 0` be in the kernel. Then `uᵀAu = −uᵀBv = vᵀDv =: s`, and `|uᵀBv| ≤ κ(|u_1| + |u_4|)(|v_2| + |v_3|)`, so `s² ≤ 4κ²|u|²|v|²`. If `A − 2κ` and `D − 2κ` were both positive definite, then `s ≥ 2κ|u|²` and `s ≥ 2κ|v|²`, so `s² ≥ 4κ²|u|²|v|² ≥ s²` and all three are equal: `uᵀ(A − 2κ)u = 0` and `vᵀ(D − 2κ)v = 0`, hence `u = 0` and `v = 0`. So one of them is not positive definite: since `N_i − 2κ > 0`, `x² ≥ (N_1 − 2κ)(N_4 − 2κ)` or `y² ≥ (N_2 − 2κ)(N_3 − 2κ)` (Lean: `G1.fourPoint_blocks_not_both_definite`). In the first case `|x| ≥ L − 2κ`, and `(N_1 + N_4 − 2|x|)(N_1 + N_4 + 2|x|) = (N_1 − N_4)² + 4(N_1N_4 − x²) ≤ e² + 8κ(N_1 + N_4) − 16κ² ≤ e² + 16κ(L + e) − 16κ²`, while `N_1 + N_4 + 2|x| ≥ 4(L − κ)` (Lean: `G1.fourPoint_block_bound`). This is the stated bound, since `16κ(L + e) − 16κ² = 16κ(L − κ) + 16κe`. The second case is the same. ∎ (Lean: `G1.fourPoint_rank3_bound`.)

Both bounds decrease when `H` grows, visibly. The first tends to 0, the second to `4κ`.

In rank 4 four vectors can be orthogonal, and the curve occurs only four times in the Jacobian of the curve of the five forms, so there is no fifth point: the method stops at rank 3.

What Lean does not contain here: that dependent vectors have a singular Gram matrix, and the identification of `N_1 + N_4 − 2|x|` with a canonical height. Both are standard.

## 6. The two lists

Fix a rational `H*` with `H* > eL` (rank at most 2) or `H* − eL > 2κ` (rank 3), and let `c*` be the bound of Proposition 2 or 3 at `H = H*`. For a solution `t0 ∉ {0, ±1}`:

- **(A)** `H < H*`. Then `ĥ(2Q1 − R0) = 4N(a_1) ≤ 4(H* + eU)`. All `Q ∈ E(ℚ)` with this property are listed, and for each `t(Q)` is tested.
- **(B)** `H ≥ H*`. Then some `z ∈ E(ℚ)` with `ĥ(z) ≤ c*` satisfies `Q4 = Q1 − z` or `Q4 = R0 − (Q1 − z)`, hence `t(Q1 − z)·t(Q1) = −1`; or the same with `Q2`, `Q3` in place of `Q1`, `Q4`, and then `t(Q2 − z)·t(Q2) = −1` with `t(Q2) = −t0`. All such `z` are listed. For each, `t(Q − z)·t(Q) + 1` is a nonconstant function on `E` (the program checks that it is not identically zero); writing `t(Q − z) = α(t) + β(t)·y` on the generic point `Q = (t, y)`, a rational solution has `t` a root of the numerator of `(αt + 1)² − β²t²q(t)` or of a denominator of `α`, `β`. All rational roots are taken, and also `t(w)` for `w = ±z`, `±z ± R0`, `z + T` with `T` of order 2 (the positions where the generic formula is not valid).

Every candidate, with both signs, is put through the exact test that the five forms are squares. `t0 ∈ {0, ±1, ±m/n, ±n/m}` are tested directly. If no candidate passes other than degenerate ones, the pair is closed.

*Why no rational solution of the equation in (B) is missed.* Fix `z ≠ O` and put `f = t(Q − z)·t(Q) + 1`, a function on `E` of degree at most 4, and `ι: Q ↦ R0 − Q`. Let `Q` be a rational point with `f(Q) = 0` and `t(Q) = τ` finite. With `g(t) = (αt + 1)² − β²t²q(t)`, which is the product `f·(f∘ι)` as a function of `t`: if `f∘ι` has no pole at `Q`, then `g(τ) = 0`, so `τ` is a root of the numerator of `g`. If `f∘ι` has a pole at `Q`, then `t(Q′ − z)` has a pole at `Q′ = ιQ`, a point above the finite value `τ`, so `α + βy` is not regular above `τ`, and then `α` or `β` has a pole at `τ`: `τ` is a root of one of their denominators. The values of `t` at the points where the addition formula for the generic point degenerates are added as well, although the identity `t(Q − z) = α + βy` is one of functions and does not need them. `z = O` gives `t(Q)² = −1`, which has no rational solution.

As a test of the solver: for about 40 pairs of known points `(Q, z)` per run, with `Q` from list (A) and `z` from list (B), `t(Q)` must be among the candidates for the equation `t(Q − z)·t(Q) = c` with `c` their actual product; this is checked in every run.

## 7. The group

Needed: the rank `r ≤ 3` and a basis of `E(ℚ)` modulo torsion.

- Upper bound for the rank: PARI's `ellrank` (2-descent and the Cassels pairing, on the curve and its 2-isogenous curves). Magma's `TwoPowerIsogenyDescentRankBound` gives the same bounds for the curves used (`magma_e0l.txt`).
- Lower bound: `r` independent points are exhibited; independence is certified by a positive lower bound for the regulator in exact arithmetic (Section 8). PARI's own lower bound is **not** used, because it can exceed the number of points found, by parity.
- Bound on the index: if `λ` is a lower bound for `ĥ` on non-torsion points, the lattice `Λ = E(ℚ)/torsion` satisfies `Reg(Λ) ≥ λʳ/γ_rʳ` with `γ_1 = 1`, `γ_2² = 4/3`, `γ_3³ = 2` (Hermite), so the index of the lattice of the exhibited points is at most `K = √(γ_rʳ·Reg′/λʳ)`. `λ` is the value returned by Sage's `height_function().min(0.01, 20)`, which its documentation describes as a number for which it has been established rigorously that every point of infinite order has greater canonical height (the algorithm of Cremona and Siksek 2006). The program lowers it by 1% before use, to be safe against any rounding of the number as returned.
- Saturation, twice. PARI's `ellsaturation` produces a basis saturated at all primes up to `K`. Independently of PARI, the program then writes a certificate that this basis is saturated at each prime `p ≤ K`, and a checker in pure Python verifies it (`satcert.py`, run in CI by `check_evidence.py`): for auxiliary primes `q` of good reduction, a homomorphism `E(ℚ) → E(F_q) → ℤ/p` that kills `pE(ℚ)` is evaluated on the basis and on the torsion points; if the matrix of these values has rank `r` plus the rank of the torsion rows, no point `c_1P_1 + … + c_rP_r + T` with some `c_i ≢ 0` is `p` times a rational point. For odd `p` the map is `P ↦ (N/p)·P̄` for primes `q` with `p` dividing `N = #E(F_q)` exactly once (`N` is counted naively), and for `p = 2` it is the 2-descent map `P ↦ (x − e)` modulo squares in `F_q`. The checker also confirms that the listed torsion points are all of them (the torsion subgroup injects into `E(F_q)`, and the greatest common divisor of a few `#E(F_q)` equals the number of listed points). In the 52 records the index bound is between 3 and 116 and the auxiliary primes are at most 6173.
- For (740,539), PARI finds one point. The second was found with Magma (`TwoDescent`, `FourDescent` with class groups under GRH, `PointsQI`) through the public calculator. GRH enters only the search: the program checks that the point is on the curve.

Taken on trust from software: the upper bound for the rank, and `λ`. The saturation no longer is.

## 8. Exact arithmetic

The decisions of the program are made in rational arithmetic.

- `BUP`, `U`, `lo4`, `up4`, `w`, hence `eL`, `eU`, `κ`, `c*`, are fractions; logarithms of integers are enclosed in rational intervals (a series for `atanh`).
- For a point `P`, `2ᵏP` is computed exactly with the duplication formula on `E_T` until `4ᵏĥ(P)` is about 2·10⁶, and `(h(x(2ᵏP)) − BUP)/4ᵏ ≤ ĥ(P) ≤ (h(x(2ᵏP)) + U)/4ᵏ`. This encloses every entry of the height pairing matrix `G` of the basis (through `ĥ(G_i)`, `ĥ(G_i + G_j)`).
- Let `G0` be a rational matrix inside the enclosure, `δ` the largest distance from an entry of `G0` to an end of its interval, and `μ > 0` rational with `G0 − μ` positive definite (checked exactly by Sylvester's criterion). Then for every real vector `v`: `|vᵀ(G − G0)v| ≤ δ(Σ|v_i|)² ≤ rδ|v|² ≤ (rδ/μ)·vᵀG0v`. With `ε = rδ/μ`: `vᵀGv ≥ (1 − ε)·vᵀG0v`.
- So `{v : vᵀGv ≤ b}` is contained in `{v : vᵀG0v ≤ b/(1 − ε)}`. The second set is listed exactly: the coordinates are bounded through the diagonal of the exact inverse of `G0`, and membership is tested in integer arithmetic. Lists (A) and (B) are therefore supersets of the sets in Section 6. In the runs `ε` is below 5·10⁻⁴.
- The index bound of Section 7 uses an upper bound for `Reg′` from the same enclosures and the floating-point number `λ` as an exact rational.

- Each record carries `G0` (as integers, scaled by 2⁶⁰), the two thresholds and the number of integer vectors in each list. `check_evidence.py` counts the vectors again with a different algorithm (the form is diagonalised over ℚ and the coordinates are bounded one after another, instead of filtering a box) and must find the same numbers.

Floating point is used to choose the basis, `k`, `μ` and `H*`, and for a diagnostic (the observed `N(a) − H` on the listed points against `[−eL, eU]`). No decision depends on it.

## 9. Status of each ingredient

| Ingredient | Status |
|---|---|
| Identities (I0), (I±), (Iinv), the model | Exact computer algebra with `m`, `n` indeterminate (`dm4_identities.py`); again per pair in `dm4.py`. Not in Lean |
| The fractions; common-factor bounds, also for Lemma C; the two resultants; the lower bounds of Lemmas A and B | Lean (`G1FourPoints.lean`) |
| The upper bounds in Lemmas A to C; the last line of Lemma C | Written out above. Not in Lean |
| `BUP`, `U` | Statements read from Cremona–Prickett–Siksek with an AI model; evaluated exactly; reduction data from Sage |
| Lattice step (Propositions 2 and 3) | Lean, as statements about the entries of a singular Gram matrix |
| Lists (A) and (B) | Exact arithmetic; superset property argued in Section 8; the vectors counted a second time by a different enumeration in CI |
| The equation in (B) | Computer algebra per `z`; completeness argued in Section 6; self-test on about 40 pairs of known points in every run |
| Proposition 1 on actual points | Checked on 72 configurations of two rational points on quadratic twists, `h(t0)` up to 230.8 (Section 11) |
| Reproduction | CI regenerates all 52 records with PARI and Sage and compares them with the committed ones |
| Rank upper bound | PARI; Magma agrees |
| Basis | Points exhibited; saturated up to a bound from Hermite's constant and Sage's `λ` (lowered by 1%); saturation by PARI and, independently, by a certificate checked in pure Python in CI |
| Torsion | Sage; completeness checked in pure Python in CI |
| The whole argument | Read by nobody but its authors (a person who is not a number theorist, and an AI model) |

## 10. Results

Nine ratios with `m ≤ 1000` that the rank-1 method leaves open are closed (`dm4_log.jsonl`): no candidate gave a non-degenerate solution, and for the ratios with a trivial point the only values found are the degenerate ones.

| Pair | rank | heights of the basis | saturated up to | `eL`, `eU` | `κ` | `H*` | `c*` | points in (A) | points `z` in (B) |
|---|---|---|---|---|---|---|---|---|---|
| (740,539) | 2 | 27.7, 86.1 | 33 | 92.5, 28.7 | 47.6 | 753 | 24.8 | 192 | 4 |
| (828,173) | 2 | 27.3, 111.7 | 10 | 91.6, 28.6 | 47.5 | 752 | 24.6 | 176 | 4 |
| (840,83) | 2 | 27.3, 40.3 | 6 | 91.2, 28.5 | 47.3 | 748 | 24.6 | 328 | 4 |
| (923,757) | 2 | 27.3, 40.3 | 6 | 95.7, 29.7 | 49.0 | 808 | 24.5 | 336 | 4 |
| (72,65) | 3 | 10.8, 11.6, 11.7 | 5 | 61.7, 19.6 | 32.8 | 722 | 151.0 | 11208 | 1052 |
| (137,7) | 3 | 10.8, 11.6, 11.7 | 9 | 66.2, 21.2 | 35.3 | 776 | 162.3 | 12456 | 1156 |
| (425,168) | 3 | 14.7, 15.3, 15.7 | 5 | 83.3, 26.2 | 43.5 | 973 | 200.1 | 11384 | 1036 |
| (593,257) | 3 | 14.7, 15.3, 15.7 | 5 | 87.8, 27.6 | 45.6 | 1026 | 209.6 | 12456 | 1132 |
| (920,129) | 3 | 16.9, 17.3, 17.3 | 5 | 92.7, 28.9 | 48.0 | 1081 | 220.8 | 11376 | 1036 |

`H*` is the threshold between lists (A) and (B), `c*` the bound of Proposition 2 or 3 at `H*`. In rank 2, `c*` is below the smallest height of a non-torsion point, so list (B) consists of the 4 torsion points. The largest of `eL`, `eU`, `κ` is 95.7; the points in (A) have `ĥ(2Q − R0)` up to about 4400. A run takes 15 to 60 seconds.

The method does not apply to the other 26 open ratios: for 22 of them the upper bound for the rank of `E` is 4 or more, for (609,200) and its partner (809,409) the bounds are `[3, 5]`, and for (820,799) and (889,320) the bound is 2 but only one point of infinite order is known.

## 11. Tests

- On every point of list (A), `N(a) − H` is compared with `[−eL, eU]`: always inside, and far from the ends (for (840,83): `[−26.9, 13.5]` against `[−91.2, 28.5]`).
- Degenerate solutions are the only rational configurations available. For `t0 = 1` and six pairs with `S` a square, `⟨a_1, a_2⟩` is between −0.16 and −0.01, where `κ` is between 30 and 50 (`dm4_degenerate.py`).
- Over quadratic and biquadratic fields, for `t0` of height 0.7 to 9.9, `N(a_1) − H` and `ĥ(Q1 ± Q2 …) − 2H` stay in a band of width 2.5 (`dm4_validate.py`). There the two points of a pair are conjugate, so this tests boundedness only.
- Cross-check: the argument for rank at most 2 covers rank 1. Run on the other 148 pairs of the census, the program closes 43 of the 122 pairs that had been closed by the rank-1 method, with the same outcome each time (`dm4_cross.jsonl`).
- **Proposition 1 on actual rational points** (`dm4_twist.py`, output in `dm4_twist.txt`). If `GL·GU` is a square at `t0`, then `Q1` and `Q2` are rational points of one quadratic twist of `C`, and `2Q1 − R0`, `2Q2 − R0`, `Q1 + Q2 − R0`, `Q1 − Q2` are rational points of the corresponding twist of `E`, with the abscissas given by (I0) and (I±). Such `t0` are the abscissas of the rational points of the genus-1 curve `w² = GL·GU`, which has positive rank for most of the pairs tried, so they exist with any height. The same holds for `Q1`, `Q3` and (Iinv) when `F0·F8` is a square. For 72 such configurations, on 8 pairs, with `h(t0)` from 2.8 to 230.8 (numerators and denominators of up to 100 digits), the script (a) checks the identities against the group law, exactly: the abscissas of `(Q1 + Q2 − R0) ± (Q1 − Q2)` computed from the addition formula are those of `2Q1 − R0` and `2Q2 − R0`; (b) encloses the canonical heights by exact doubling. Result: `N(a_i) − H` lies between 2.4 and 13.5 in all 144 cases, where the bounds are `[−61.7, 19.6]` to `[−92.7, 28.9]`; `|⟨a_1, a_2⟩|` and `|⟨a_1, a_3⟩|` are at most 10.2, where `κ` is 32.8 to 48.0. The constants do not grow with the height. (For points that are rational only on a twist the height-difference bounds have to be valid over number fields; the script uses Silverman's, in the form that Sage's documentation gives as valid over all number fields. For these curves they agree with the constants of Section 4 to one decimal.)
- All four points are never rational on a twist in the range searched: for `m, p, q ≤ 400` no `t0 = p/q ≠ 1` makes `F0·F8` and `GL·GU` both squares (`scripts/search/twist_find.c`, `twist_search.txt`). So Propositions 2 and 3, which are in Lean, cannot be tested on an actual four-point configuration; an earlier search without twists found none either (`m ≤ 250`).

## 12. Relation to known work

The principle is not new. The method of Demjanenko and Manin applies to a curve with `k` independent morphisms to an elliptic curve whose rank over ℚ is less than `k`, and bounds the height of its rational points. Here the curve is the curve of the five forms, the morphisms are `t0 ↦ Q1, …, Q4`, and `k = 4`: the condition "rank at most 3" is the classical one, and the method stops where the classical method stops. Explicit applications that we found use two morphisms and rank 1 (Kulesz 1999; Kulesz, Matera and Schost 2004) or three morphisms and rank at most 2 (Girard and Kulesz 2005). We have not read these papers; this description is from their abstracts and from the introduction of Viada 2018.

What this document adds, as far as we can tell: the three identities, which give the heights and four of the six pairings up to bounded errors with constants that are explicit in `m`, `n`; the lattice step in the form of Propositions 2 and 3, which needs only those four pairings; exact arithmetic and certificates for the two lists; and applications with four morphisms to curves of rank 3. Whether any of that is new in a sense that matters to specialists is their call, not ours.

## 13. What a reader is asked to check

1. Section 2: are (I0), (I±), (Iinv) plausible as stated, and is exact computer algebra in the function field an acceptable proof of them?
2. Section 3: the three lemmas.
3. Section 4: the use of the height-difference bounds for a model that is integral and not minimal, at the real place and at the primes of bad reduction; Proposition 1.
4. Section 5: Propositions 2 and 3 (their inequalities are in Lean; the question is whether they are the right statements).
5. Section 6: the completeness of the case distinction and of the solution of the equation in (B).
6. Section 7: the index bound; whether Sage's `height_function().min` returns a proved lower bound; whether PARI's upper bound for the rank is unconditional.
7. Section 8: the superset argument.
8. Section 12: is this explicit form of the method, with four morphisms and rank 3, in the literature?

Answers of the form "this step is wrong" are as welcome as confirmations: https://github.com/drew2apps/magic-square-of-squares-lean/issues.

## References

- J. E. Cremona, M. Prickett, S. Siksek, Height difference bounds for elliptic curves over number fields, J. Number Theory 116 (2006).
- V. A. Demjanenko (1966) and Yu. I. Manin (1969), on rational points of curves with several maps to an elliptic curve; see J.-P. Serre, *Lectures on the Mordell–Weil theorem*, Section 5.2.
- M. Girard, L. Kulesz, Computation of sets of rational points of genus-3 curves via the Dem'janenko–Manin method, LMS J. Comput. Math. 8 (2005) (three independent morphisms to an elliptic curve of rank at most two; not read).
- L. Kulesz, Application de la méthode de Dem'janenko–Manin à certaines familles de courbes de genre 2 et 3, J. Number Theory 76 (1999) (not read).
- E. Viada, An explicit Manin–Dem'janenko theorem in elliptic curves, Canad. J. Math. 70 (2018) (introduction read).
- J. E. Cremona, S. Siksek, Computing a lower bound for the canonical height on elliptic curves over ℚ, ANTS VII (2006) (the algorithm behind Sage's `height_function().min`; not read).
- T. Fisher, Higher descents on an elliptic curve with a rational 2-torsion point, Math. Comp. 86 (2017) (the routine behind Magma's `TwoPowerIsogenyDescentRankBound`).
- J. H. Silverman, The difference between the Weil height and the canonical height on elliptic curves, Math. Comp. 55 (1990) (as quoted by Cremona, Prickett and Siksek).
