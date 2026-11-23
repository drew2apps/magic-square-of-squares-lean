# The last pair: (627,364) and (991,263)

23 November 2026, after release 0.4.0. Written by Claude (an AI model made by Anthropic), directed by Drew (`drew2apps`).

**Status. Nobody has read this.** It is an argument written out step by step. The numbers that depend on the curves are computed by programs in `scripts/second_descent/`; the constants of the height bound (`C₁ = 72.251` and those inside `c₁`, `c₂`) are copied from one published paper, as read by an automated reader, and the argument depends on that reading. There is a positive control on real points, and two programs for each arithmetic input where that was possible. Until a number theorist has checked it, it is a proposed proof and **the census is not changed**: (627,364) and (991,263) stay open there.

The file has two parts that should be judged separately.

- **Theorem A** is a statement about two explicit curves of genus 1 over `ℚ`. Its proof does not use the Gaussian descent or the second descent of this repository.
- **The Corollary** (no magic square of squares has ratio (627,364) or (991,263)) follows from Theorem A *and* from the two descents (`GAUSS_DESCENT.md`, `SECOND_DESCENT.md`), which are experiments that nobody has reviewed.

## 1. Statement

Let `m = 627`, `n = 364`, so `m² + n² = 725²`, and `d = 521266`. With `s` a variable,

    q₊ = (m² + n²)s² − 2(m² − n²)s + (m² + n²),
    Q₁ = q₊·(m − ns)(n + ms),        Q₂ = q₊·(m + ns)(n − ms).

Let `K = ℚ(i)`, with conjugation `σ`.

**Theorem A.** The rational numbers `t` for which `d·Q₁(−it)` and `d·Q₂(−it)` are both squares in `K` are `t = 1` and `t = −1`. (At `t = ∞` and `t = 0` neither is a square.)

**Corollary (conditional on the two descents).** There is no 3×3 magic square of nine distinct squares with `(b + d) : (h + f)` equal to `627 : 364` or to `991 : 263`.

*Why the Corollary follows.* By the first descent (if it is right), at a solution `t` the classes of the ten linear factors in `K*/K*²` are those of a degenerate solution; for `Q₁` and `Q₂` this class is the rational number `d` (`SECOND_DESCENT.md`, sections 1 and 2). So a solution is a `t` as in Theorem A, hence `t = ±1`, and these are the degenerate solutions (`classify_t`). The ratio (991,263) is the partner of (627,364), and a solution for one gives a solution for the other (`partner_forms` in Lean).

## 2. The objects

- `C₁: Y² = d·Q₁(s)` and `C₂: Y² = d·Q₂(s)`, curves of genus 1 over `ℚ`, each with a rational point (`s = m/n`, resp. `s = −m/n`, `Y = 0`).
- `ψⱼ: Cⱼ → Eⱼ` the isomorphism to the elliptic curve given by `weierstrass()`: with `r` the rational root and `s = r + A₀/x`. In minimal form

      E₁: y² = x³ + x² + 15223208213202973031750280786182975·x + 372517508636065850025943055293341065999858640587423,

  and `E₂` is the twist of `E₁` by `−1` (checked by the programs). Fix an isomorphism `θ: E₂ → E₁` over `K`; then `σθ = −θ`, so `θ` carries the points `R` of `E₂(K)` with `σR = −R` onto `E₁(ℚ)`.
- For a `t` as in Theorem A put `s = −it`, so `σ(s) = −s`. Choose `Y₁, Y₂ ∈ K` with `Yⱼ² = d·Qⱼ(s)` and let

      P = ψ₁(s, Y₁) ∈ E₁(K),    R = ψ₂(s, Y₂) ∈ E₂(K),
      A = P + σP ∈ E₁(ℚ),       B = θ(R − σR) ∈ E₁(ℚ).

  Changing the sign of `Y₁` or `Y₂` changes the sign of `A` or `B`.

## 3. The seven lemmas

**Lemma 1 (the group).** `E₁(ℚ) = ℤ·G ⊕ {O, T}`, where `T` is the point of order 2 and `G = (83214074135015341, 47069367709023802879687500)` on the minimal model. `G` is the point `A` of the degenerate solution `t = 1`. Its canonical height is `ĥ(G) = 20.00665…` in the normalisation of PARI, Sage and Magma (twice Silverman's).

*Proof, by computation.* Rank at most 1: PARI's `ellrank` gives the bounds 1 and 1, and Magma's `RankBounds` gives 1 and 1 and its `TwoPowerIsogenyDescentRankBound` gives 1 (the 2-Selmer group has order 64, so a 2-descent alone would not do; `magma_last_pair_raw.txt`). `G` has infinite order, so the rank is 1. Torsion: `{O, T}` (Sage; Magma: invariants `[2]`). Generator: Sage's `height_function().min` gives `0.045162` as a lower bound for the height of a point of infinite order; if `G = k·G₀ + torsion` then `k² ≤ 20.007/0.045162`, so `k ≤ 21`, and `G` is not `ℓ` times a point plus torsion for `ℓ = 2, 3, 5, 7, 11, 13, 17, 19`, by reduction modulo a prime in each case (`saturation.py`). Independently, Magma's `Saturation([G], 200)` returns the two points `T` and `G`: `G` is not divisible by any prime up to 200. (The last line of that first request prints `false`: it compared the first point returned, which is `T`, with `±G`. The request was sent again with a correct test, also with the primes up to 1000: one point of infinite order is returned, it is `±G` plus torsion, and the other is a torsion point; the last request in the transcript.) About the lower bound: in the source of Sage (10.8, `height.py`) `min` is `min_gr`, the bound of Cremona and Siksek for points of good reduction, divided by the square of the least common multiple of the Tamagawa exponents and the number of real components, and its documentation states the result as a proved bound for all points of infinite order. The height of `G` itself is enclosed in exact arithmetic (`height_enclosure.py`: the point `2⁷·G` computed exactly, and the bounds of Cremona, Prickett and Siksek for the difference between the naive and the canonical height): `20.002899 ≤ ĥ(G) ≤ 20.009091`. The index bound uses the upper end, `k² ≤ 20.0091/0.045161`, and Lemma 5 uses the lower end. Since Magma finds `G` saturated at all primes up to 1000, the conclusion would survive a lower bound for the height two thousand times smaller than Sage's (any bound above `2.1·10⁻⁵` gives `k < 1000`). ∎

So for a `t` as in Theorem A there are integers `a, b ≥ 0` with `A = ±(aG + εT)`, `B = ±(bG + ε′T)`.

**Lemma 2 (the curve of solutions).** Let `f₁ = d·Q₁(−it)`, `f₂ = d·Q₁(it)`, `f₃ = d·Q₂(−it)`, `f₄ = d·Q₂(it)` in `K[t]`, and let `V` be the curve of the `(t, y₁, y₁′, y₂, y₂′)` with `y₁² = f₁`, `y₁′² = f₂`, `y₂² = f₃`, `y₂′² = f₄`.

1. `f₁f₂f₃f₄ = c·g²` with `c` a square in `K`, and no other product of the `fᵢ` is a constant times a square. So `V` has two components `W`, each irreducible, a cover of the line of `t` with group `(ℤ/2)³` and degree 8.
2. The `fᵢ` have 8 distinct roots in all, each a root of exactly two of them. So the cover has 8 branch points with inertia of order 2, and `2g − 2 = 8·(−2) + 8·4`: **`W` has genus 9.** (The first version of this argument said 16 branch points and genus 25. That was wrong; a reviewer found it.)
3. On `W` define `α = ψ₁(−it, y₁) + ψ₁(it, y₁′)` and `β = θ(ψ₂(−it, y₂) − ψ₂(it, y₂′))`, morphisms `W → E₁`. For integers `u, v`, the map `uα + vβ` is constant only if `u = v = 0`.

*Proof.* 1 and 2 are checked by factoring (`closure_estimate.py` prints both). (Which of the two components a solution lies on can be decided, since for rational `t` the product of the four square roots is the positive number `N(Y₁)·N(Y₂)`; this is not needed, because everything below is proved for both.) For 3: let `ω` be the invariant differential of `E₁`. For morphisms `φ, φ′` to an elliptic curve, `(φ + φ′)*ω = φ*ω + φ′*ω`. The pull-back of the invariant differential of a curve `y² = f(t)` of genus 1 is a nonzero multiple of `dt/y`. So

    α*ω = c₁·dt/y₁ + c₂·dt/y₁′,        β*ω = c₃·dt/y₂ + c₄·dt/y₂′,        all cᵢ ≠ 0.

The group of the cover acts on the four differentials `dt/y₁, dt/y₁′, dt/y₂, dt/y₂′` by four characters: those of `y₁`, `y₁′`, `y₂`, and, since `y₂′ = ±√c·g/(y₁y₁′y₂)`, their product. These are four different characters, so the four differentials, which are not zero, are linearly independent. Hence `(uα + vβ)*ω = 0` forces `uc₁ = uc₂ = vc₃ = vc₄ = 0`, that is `u = v = 0`; and in characteristic 0 a morphism to an elliptic curve is constant exactly when it pulls `ω` back to 0. ∎

Everything in Lemma 2 holds for each of the two components (they differ by the sign of `y₂′`, and the argument does not see it). Let `W` be either of them and `C ⊂ E₁ × E₁` its image under `(α, β)`; Lemmas 3 to 5 are proved for both images, with the same numbers. `C` is an irreducible curve. `E₁` has no complex multiplication (its `j`-invariant is not an integer), so the proper algebraic subgroups of `E₁²` of positive dimension are given by one equation `ux + vy = 0`, and by Lemma 2 **`C` is not contained in a translate of a proper algebraic subgroup** (it is "transverse"). Then `C` has geometric genus at least 2: an abelian variety contains no curve of genus 0, and a curve of genus 1 in it is a translate of an elliptic curve.

For a `t` as in Theorem A, the point `(t, Y₁, σY₁, Y₂, σY₂)` lies on one of the two components, and its image is `(A, B)`: a point of `C` with both coordinates in `ℤG ⊕ {O, T}`.

**Lemma 3 (degree).** `deg C ≤ 96` in the embedding of `E₁²` in `ℙ⁸` (Weierstrass coordinates, then Segre).

*Proof.* The map `W → C₁`, `(t, …) ↦ (−it, y₁)`, has degree `8/2 = 4`, and `ψ₁` is an isomorphism; the same for the three other maps. The degree is a positive definite quadratic form on the group of morphisms from the Jacobian of `W` to `E₁`, so `deg(φ + φ′) ≤ 2 deg φ + 2 deg φ′`, and `deg α, deg β ≤ 16`. The two projections of `C` have degree at most 16, the hyperplane section of `E₁ ⊂ ℙ²` has degree 3, and the degree of `C` under the Segre embedding is `3·deg(pr₁) + 3·deg(pr₂) ≤ 96`. ∎ (The four maps go through different characters, so the true degrees of `α` and `β` should be 8 and the bound 48. Nothing below needs it; the numbers are also given for 192.)

**Lemma 4 (height of the curve).** `h₂(C) ≤ 7.78·10⁵`.

*Proof.* Notation and facts from [CVV] (section 5 below): `h₂` is the Weil height with the Euclidean norm at the infinite places; on `E₁²` heights are sums over the two factors; `ĥ` is their Néron–Tate height, 1.5 times PARI's; (P) for a point `x` of `E₁`: `−c₂ ≤ h₂(x) − ĥ(x) ≤ c₁`, with `c₁ = 236.45`, `c₂ = 216.64` for algebraic points; (Z) `h₂(C) ≤ 2·deg C·μ` whenever the points of `C` with `h₂ ≤ μ` are Zariski dense, which on an irreducible curve means infinitely many.

Take `t = ζ` a root of unity and a point `w` of `W` above it. Then `s = ∓iζ` is a root of unity, and the point `x = ψ₁(s, y)` has, on the model `y² = x³ + Ax + B` with integers `A, B`, the coordinate `X = a·A₀/(s − r) + b` with rational `a, b, r, A₀`. So `h(X) ≤ h(aA₀) + h(r) + h(b) + 2 log 2 = 110.21`; `h(X : Y : 1) ≤ 1.5·h(X) + 0.5·log(1 + |A| + |B|)`, because `|Y|² ≤ (1 + |A| + |B|)·max(1, |X|)³` at the infinite places and `|Y|² ≤ max(1, |X|)³` at the finite ones; `h₂ ≤ h + 0.5·log 3`; and by (P) `ĥ(x) ≤ 446.3`. On `E₂`: `448.4`. By the parallelogram law `ĥ(α(w)) ≤ 2·446.3 + 2·446.3` and `ĥ(β(w)) ≤ 4·448.4` (`θ` is an isomorphism and does not change `ĥ`). By (P) again, `h₂(α(w), β(w)) ≤ 4·446.3 + 4·448.4 + 2·236.45 = 4051.4`. There are infinitely many such points, so by (Z) `h₂(C) ≤ 2·96·4051.4`. ∎ (In the program every quantity that enters a bound is a rational number, computed exactly: a logarithm of an integer is replaced by its number of bits times `0.6931471806 > log 2`, `log 3` by `1.0987`, and the constants of the paper are decimal fractions. The comparison of the bound with the range of the sieve is made on exact rationals. The height of `G` enters through the proved lower bound 20.002899 of Lemma 1; no floating-point number enters a bound.)

**Lemma 5 (the bound).** For a `t` as in Theorem A, `a² + b² ≤ 6.03·10⁹/30.004 ≈ 2.0·10⁸`; in particular `a, b ≤ 14,175`.

*Proof.* [CVV, Theorem 1.2] for `N = 2`: if `E: y² = x³ + Ax + B` has no complex multiplication, `C ⊂ E²` is an irreducible curve of genus at least 2, and `P ∈ C` has rank at most 1, that is, lies in an algebraic subgroup of dimension at most 1 (which is the case for every point with coordinates in a group of rank 1, as here: `2b·x ∓ 2a·y = 0`), then

    ĥ(P) ≤ C₁·h₂(C)·deg C + C₂(E)·(deg C)² + C₃(E),    C₁ = 72.251,  C₂(E) = C₁(6.019 + 4c₁(E)),  C₃(E) = 4c₂(E).

With Lemmas 3 and 4 and the larger of the two sets of constants, the right side is at most `6.03·10⁹` (`2.42·10¹⁰` if one only assumes `deg C ≤ 192`). And `ĥ(A, B) = 1.5·ĥ(G)·(a² + b²) ≥ 30.004·(a² + b²)`. ∎

**Lemma 6 (the sieve).** For a `t` as in Theorem A, if `a, b < 59,631,000` then `(a, b) = (1, 1)`.

*Proof, by computation.* Let `p ≡ 1 (mod 4)` be a prime of good reduction for `E₁`, not dividing `d` nor the numerators or denominators of the `r` and `A₀`, and with `j(E₁) ≢ 0, 1728 (mod p)`; let `ι² = −1` in `𝔽_p` and `𝔭` the prime of `K` with `i ↦ ι`. Reduction at `𝔭` is a homomorphism `E₁(K) → E₁(𝔽_p)` that commutes with the map `x ↦ s = r + A₀/x` to the projective line (an automorphism of the line over `ℤ_(p)` composed with the coordinate `x`). So the reduction of `P` lies above `s̄ = −ι·t̄` and that of `σP` above `+ι·t̄`, where `t̄ ∈ ℙ¹(𝔽_p)` is the reduction of `t`. Hence the reduction of `A` is `X + X′` for some point `X` of `E₁(𝔽_p)` above `−ι·t̄` and some `X′` above `ι·t̄`; and the reduction of `B` is `±τ(Z − Z′)` for points `Z, Z′` of `E₂(𝔽_p)` above the same two values, with `τ` an isomorphism `E₂ → E₁` over `𝔽_p` (there are two, `±τ`). The program lists, for every `t̄ ∈ ℙ¹(𝔽_p)` for which such points exist, all these sums and differences, and for each the residues `a, b` modulo the order of `G` with `±(aG + εT)` equal to it. A pair `(a, b)` that is in none of these lists for some `p` cannot occur. With `p = 10009, 20089, 30013, 40009, 50021, 60013` the only pair with `0 ≤ a, b < 59,631,000` that is left is `(1, 1)` (`mw_sieve_big.py`); with the first four, the only pair up to `596,309` (`mw_sieve.py`), which is already more than Lemma 5 needs: the proof of Theorem A uses only this smaller run; the larger run is a margin, and a second implementation of the lifting (arrays instead of a loop) that gives the same answer on the common range. Each of the two programs tests its own lifting: pairs spread over the whole range (five up to `(596,309, 3)`, six up to `(59,630,999, 3)`) are planted by adding their residues to every local set, and the program checks that all of them come back and nothing else does. ∎

**Lemma 7 (the fibre).** A `t` as in Theorem A with `a = 1` is `1` or `−1`.

*Proof, by computation.* `A = ±G` or `±G + T` means `σP = c − P` with `c` one of four points. Since `σ(s) = −s`, the coordinate `x` of `P` satisfies `x(c − P) = −A₀x/(2rx + A₀)`: a polynomial equation over `K`, not identically satisfied, whose roots give four values of `s`; the rational `t = is` among them are `±1` (`fibre.py`; the equation is squared, so `−c` is covered with `c`). The same computation for `B` gives the same answer. ∎

*Proof of Theorem A.* For a `t` as in the theorem, Lemma 5 gives `a, b ≤ 14,175`, Lemma 6 then gives `(a, b) = (1, 1)`, and Lemma 7 gives `t = ±1`. At `t = ±1` both numbers are squares (the degenerate solution). `t = ∞` and `t = 0` are not points: `d` times the leading coefficient, and `d·Qⱼ(0)`, are not squares in `K` (`fibre.py`). ∎

## 4. The numbers

| Quantity | Value | Program |
|---|---|---|
| `ĥ(G)` (PARI, Sage, Magma agree) | 20.006653 | `saturation.py`, `magma_last_pair_raw.txt` |
| `ĥ(G)`, proved enclosure | between 20.002899 and 20.009091 | `height_enclosure.py` |
| lower bound for heights on `E₁` | 0.045162 | `saturation.py` |
| genus of `W`; branch points | 9; 8 | `closure_estimate.py` |
| `deg C` | ≤ 96 | Lemma 3 |
| height of points above roots of unity | ≤ 4051.4 | `closure_estimate.py` |
| `h₂(C)` | ≤ 7.78·10⁵ | `closure_estimate.py` |
| bound of [CVV] for `ĥ(A, B)` | ≤ 6.03·10⁹ | `closure_estimate.py` |
| sieve: `ĥ(A, B)` if not `(1, 1)` | ≥ 1.07·10¹⁷ | `mw_sieve_big.py` |
| so `a, b` are at most | 14,175 (28,350 with `deg C ≤ 192`) | `closure_estimate.py` |
| the sieve covers | 59,630,999 (596,309 with four primes) | `mw_sieve_big.py`, `mw_sieve.py` |

## 5. What is used from the literature

[CVV] S. Checcoli, F. Veneziano, E. Viada, *The explicit Mordell conjecture for families of curves* (with an appendix by M. Stoll), Forum of Mathematics, Sigma 7 (2019), e31; arXiv:1602.04097. Used: Theorem 1.2 (the bound of Lemma 5), Proposition 3.2 (the comparison (P), with the constants of the paper for algebraic points and for rational points), Zhang's inequality in the form of their Theorem 2.2 (Z), their definition of `h₂` and of `ĥ` (which they note is three times Silverman's), and the remark that for a curve in `E²` transverse is the same as genus at least 2. **These statements were read from the arXiv version (v5) through an automated reader, in two separate readings with the same result; they have not been compared with the published paper by a person.** That is the first thing to check. The paper itself notes that PARI's canonical height is two thirds of its `ĥ`, which is the factor 1.5 used here.

## 6. Tests on real values

- **The degenerate solution.** `(1, 1)` survives every prime of the sieve, as it must; the fibre equations return `t = ±1`; `G` satisfies (P): `h₂(G) − ĥ(G) = 34.5`, between `−172.9` and `192.7`.
- **A planted solution that is not degenerate** (`planted_control.py`, `planted_control.txt`). For any ratio, at `t₀ = (m + n)/(m − n)` the norms of `Q₁(−it₀)` and `Q₂(−it₀)` are squares, so with other twists `d₁, d₂` the same system has the solution `t₀`. (1) For six ratios and 14 primes each, the reductions of the true points `A` and `B` lie in the sets that the sieve computes from `t₀` modulo `p`: 84 of 84. (2) For (2,1), (3,2), (10,3) both curves have rank 1; the planted solution has `(a, b) = (2, 1)`, computed exactly, and the whole sieve leaves, among the pairs up to 30, exactly `(2, 1)` and `(0, 1)`. So the sieve neither removes a true point nor leaves anything unexplained. (The twists differ there, so the second point lives on another curve with its own generator; the code for the reductions, the signs and the twist is the same.)
- **The whole route on every ratio where it applies** (`route_all.py`, `route_all.txt`). `E₁` has rank exactly 1 for 24 of the 74 ratios of square type. On each of them: the degenerate point generates `E₁(ℚ)` modulo torsion; the bound of Lemma 5 is between 13,262 and 14,725; the sieve covers it and leaves, within the bound, only `(1, 1)`; the fibre is `t = ±1`. So the analogue of Theorem A holds for all 24 by this route alone. For the 23 other than (627,364) this agrees with what is known: each of them is closed either in the census (ranks of other curves) or by the second descent on another quartic (for (440,279) with an `L`-value), and neither uses this curve or this theorem. Details: 18 of the 24 sieve runs lift from the first prime only (the orders of `G` modulo the first two primes are not coprime, and the output says so), a simpler path than the one used for (627,364); for (156,133) the prime 10009 is not admissible and 10037 is used; for (99,20), (493,276) and (684,37) four primes leave `(0, 0)` as well, the pair of `t = 0` and `∞`, which are not points over `ℚ(i)` but are points modulo each of the four primes (the value there has a squarefree part that is a square modulo all four); one or two more primes remove it. The files `*_440_279.txt` and `mw_sieve_99_20.txt` are the single runs for two of these ratios.
- **The count of points.** If `f₁f₂f₃f₄` is a square, the `t` modulo `p` at which all four `fᵢ` are squares are 1/8 of all, not 1/16: 1,208 of 10,009 and 2,420 of 20,089.
- **Independence, numerically** (`two_maps_independent.py`): modulo two primes, a relation `A = ±B + torsion` holds for 8 and 4 values of `t` among 1,208 and 2,420. This only illustrates Lemma 2; the proof is the argument with differentials.

## 7. What a reader has to check, in the order of risk

1. **[CVV]**: the three statements and their normalisations, against the published paper (section 5).
2. **Lemma 6**, the reduction argument: that the points of the sieve's lists are all the possible reductions (signs, the two embeddings, `t̄ = ∞`, the points with `s = r`), and the identification of `E₂` with the twist. The planted control tests this on real points; it is not a proof.
3. **Lemma 2**: short, and the genus was wrong once.
4. **Lemma 4**: the chain of inequalities, in particular the use of (P) for algebraic points and of (Z) with the points above roots of unity.
5. **Lemma 1**: the lower bound for the height is Sage's (its source and documentation say it is a proved bound for all points of infinite order; the computation behind it is Sage's and was not repeated by hand); the rank bound is PARI's and Magma's; the saturation is checked by two programs.
6. **For the Corollary only**: the first descent, which gives the class `d`, and the second descent's set-up. These are the experiments of `GAUSS_DESCENT.md` and `SECOND_DESCENT.md`, with their own lists of what could be wrong.

The route has been run on all 24 ratios of square type for which `E₁` has rank exactly 1 (section 6); it gives the known answer on the 23 that are closed by other means. That tests the programs and the constants; it does not test whether the sieve could wrongly remove a true point, which is what the planted control is for, and it does not replace a reader.

## 8. Reproduce

```
cd scripts/second_descent              # Python with passagemath and cypari2, see scripts/g1_full/requirements.txt
python3 saturation.py 627,364 521266                                              # Lemma 1
python3 closure_estimate.py 627,364 521266 59630999                               # Lemmas 2 to 5
python3 mw_sieve.py 627,364 521266 1000000 10009 20089 30013 40009                # Lemma 6, to 596,309 (under a minute)
python3 mw_sieve_big.py 627,364 521266 100 10009 20089 30013 40009 50021 60013    # Lemma 6, to 59,631,000 (2 minutes)
python3 fibre.py 627,364 521266                                                   # Lemma 7
for r in 2,1 3,2 10,3 4,1 4,3 5,2; do python3 planted_control.py $r 30; done      # the planted control
python3 route_all.py 4                                                            # the route on all 24 ratios (10 to 20 minutes)
```

The full CI run executes all of these in a job of its own (`last-pair`) and compares the outputs with the committed files; the job `quick`, on every push, runs all but the two sieves. Magma is not run in CI; its answers are in `magma_last_pair_raw.txt`.

## 9. The construction for every ratio

Nothing in Lemma 2 is special to (627,364). `lemma2_general.py` checks the following as identities in the variables `m, n, d, t` over `ℚ(i)`. Put `A = q₊(−it)`, `A′ = q₊(it)`, `B = (m + int)(n − imt)`, `B′ = (m − int)(n + imt)`. Then

    f₁ = d·A·B,    f₂ = d·A′·B′,    f₃ = d·A·B′,    f₄ = d·A′·B,
    A·A′ = GL·GU/4,        B·B′ = F0·F8,        f₁f₂f₃f₄ = d⁴·(F0·F8·GL·GU/4)²,
    Q₂(s) = −s⁴·Q₁(1/s),

where `F0, F8, GL, GU` are four of the five forms of the family, in `t`. The discriminant of `F0·F8·GL·GU` in `t` is a constant times a power of

    m·n·(m − n)·(m + n)·(m² + n²)·((m² − mn − n²)² + 3m²n²)·((m² + mn − n²)² + 3m²n²),

which is not zero for integers `m > n > 0`. So for every ratio and every `d`:

**Theorem B (structure).** (1) The four polynomials `A, A′, B, B′` are squarefree and pairwise coprime, each root of `f₁f₂f₃f₄` is a root of exactly two of the `fᵢ`, and a product of some of the `fᵢ` is a constant times a square only for none or all of them (the pairs of `fᵢ` that share roots form a cycle `1–3–2–4–1`). (2) The curve of solutions `V` has two components, each a cover of the line of `t` with group `(ℤ/2)³`, branched exactly over the eight roots of `F0·F8·GL·GU`, of genus 9. (3) Each component is an **unramified double cover of the curve `X: u² = F0·F8, v² = GL·GU`**, which has genus 5 (`16 = 2·8 + 0` in the formula of Riemann and Hurwitz). (4) `Y² = d·Q₂` is the twist by `−1` of `Y² = d·Q₁`, so both maps `α, β` of section 3 go to the same curve `E₁`, and by the argument with four characters no `uα + vβ` is constant. (5) For a rational `t` the product of the four square roots is `N(Y₁)·N(Y₂) > 0`, which says on which component the point lies.

**Theorem C (the degenerate point, for every ratio).** `Q₁(∓i) = 2(m² − n²)·(m ± in)²`, so `t = ±1` satisfies the two conditions with `d = 2(m² − n²)` for every ratio, and `G = P + σP` is a rational point of `E₁`. With `τ = m/n` (and `n = 1` after scaling, which changes `d·Q₁` by a square) it is, on the model of `weierstrass()`,

    G = ( 2τ(τ² − 1)(τ² + 1)²,  4τ(τ² − 1)²(τ² + 1)² ).

**`G` has infinite order for every rational `τ` other than `0, ±1`**, that is, for every ratio. Proof by computation over `ℚ(τ)` (`degenerate_section.py`): a rational torsion point has order at most 10 or 12 (Mazur); the discriminant of the model has no rational root other than `0, ±1`, so specialisation at any other `τ₀` is a homomorphism; and if `G(τ₀)` had order `k`, a polynomial in `τ` (the numerator of `y(jG)` for `k = 2j`, of `x((j+1)G) − x(jG)` for `k = 2j + 1`) would vanish at `τ₀`. These ten polynomials, of degree 9 to 162, have no rational root other than `0, ±1`. For (627,364) the formula gives the point `G` of Lemma 1 (same height, 20.006653). So `E₁` has rank at least 1 for every ratio. That explains an observation of `SECOND_DESCENT.md` (this curve never has rank 0: 74 ratios of 74 below 1000 and 40 of 40 between 1000 and 3000), and it says where the effective statement below applies: exactly on the fibres of rank 1 of this family, with a point of infinite order that costs nothing to find.

*What this says about the problem.* A magic square of squares with ratio `m : n` gives a rational point of `X` (two of the conditions "product of two forms is a square"). At a rational point of `X` the numbers `N(A(t)) = GL·GU/4` and `N(B(t)) = F0·F8` are squares, so `A(t)·B(t)` has square norm, and an element of `ℚ(i)` with square norm is a rational number times a square: **the class `d` of `A(t)·B(t)` is a rational number for every rational point of `X`**, whatever the type of the ratio, and the point lifts to the double cover with that `d`. So the curves of this file are the twists of one unramified double cover of the classical curve `X`, and the two descents of this repository are used, for the Corollary, only to say which twists `d` can occur. That is a much smaller statement to check than the whole Gaussian descent on five forms: a 2-descent on one double cover. It has not been written out separately, and it is the natural next step.

**Corollary of Theorem B and of [CVV] (effective, rank 1).** For every `m, n, d` such that `E₁ = Jac(Y² = d·Q₁)` has no complex multiplication and rank at most 1 over `ℚ`, the set of rational `t` with `d·Q₁(−it)` and `d·Q₂(−it)` both squares in `ℚ(i)` is finite and can be determined: the points `(A, B)` have height at most the bound of Lemma 5, computed by `closure_estimate.py` for the given numbers, and a sieve covers that range. The statement and its proof are those of Theorem A, with the same lemmas; what depends on the numbers is computed by the same programs (they take the degenerate point as generator, so as written they need a ratio of square type and that point of infinite order). Nobody has read this either.

## 10. Twelve more ratios, between 1000 and 3000

`SECOND_DESCENT.md` leaves 67 ratios with `1000 < m ≤ 3000` (two descents, then rank 0 from descent bounds or `L`-values). Of these, 40 are of square type, and for 12 of them `E₁` has rank exactly 1 (`two_map_scope_wide.txt`; PARI's bounds, not checked with Magma). On each of the 12 the whole route gives the same answer as for (627,364): the degenerate point generates, the bound of Lemma 5 (13,651 to 14,132) is covered by the sieve, only `(1, 1)` is left within it, and the fibre is `t = ±1` (`route_wide.txt`):

    (1064,423) (1332,1045) (1395,532) (1785,688) (2015,1632) (2107,276) (2108,195) (2135,1248) (2331,1300) (2752,825) (2765,2508) (2867,756).

Six of them have their partner in the list of 67: (1487,641), (2377,287), (1927,863), (2473,1097), (2383,1831), (2303,1913). So, if the two descents and this route are right, 18 of the 67 are closed and 49 are left with `1000 < m ≤ 3000`; with (627,364) and (991,263), none is left with `m ≤ 1000`. For 15 other ratios of square type PARI's first bounds for the rank of `E₁` are 1 and 3; with more effort one of them, (2257,1176), is found to have rank 3, and 14 stay undecided. These are results of the experiment, not of the census.

