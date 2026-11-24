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

*Proof, by computation.* Rank at most 1: PARI's `ellrank` gives the bounds 1 and 1, and Magma's `RankBounds` gives 1 and 1 and its `TwoPowerIsogenyDescentRankBound` gives 1 (the 2-Selmer group has order 64, so a 2-descent alone would not do; `magma_last_pair_raw.txt`). `G` has infinite order, so the rank is 1. Torsion: `{O, T}` (Sage; Magma: invariants `[2]`). Generator: Sage's `height_function().min` gives `0.045162` as a lower bound for the height of a point of infinite order; if `G = k·G₀ + torsion` then `k² ≤ 20.007/0.045162`, so `k ≤ 21`, and `G` is not `ℓ` times a point plus torsion for `ℓ = 2, 3, 5, 7, 11, 13, 17, 19`, by reduction modulo a prime in each case (`saturation.py`). Independently, Magma's `Saturation([G], 200)` returns the two points `T` and `G`: `G` is not divisible by any prime up to 200. (The last line of that first request prints `false`: it compared the first point returned, which is `T`, with `±G`. The request was sent again with a correct test, also with the primes up to 1000: one point of infinite order is returned, it is `±G` plus torsion, and the other is a torsion point; the last request in the transcript.) About the lower bound: in the source of Sage (10.8, `height.py`) `min` is `min_gr`, the bound of Cremona and Siksek for points of good reduction, divided by the square of the least common multiple of the Tamagawa exponents and the number of real components, and its documentation states the result as a proved bound for all points of infinite order. The height of `G` itself is enclosed in exact arithmetic (`height_enclosure.py`: the point `2⁷·G` computed exactly, and the bounds of Cremona, Prickett and Siksek for the difference between the naive and the canonical height): `20.002899 ≤ ĥ(G) ≤ 20.009091`. Lemma 5 uses the lower end. For the index, `saturation.py` does not rely on the exact value of Sage's bound: it takes half of it (room for the conversion of Sage's number to a float), which gives `k ≤ 29`, and then shows by reduction that `G` is not `ℓ` times a point plus torsion for every prime `ℓ ≤ 199`. So `G` generates as soon as the points of infinite order have height above `20.0091/211² < 0.00046`, a hundred times less than Sage's bound; with Magma's saturation at the primes up to 1000, above `2·10⁻⁵`. What remains trusted here is that Sage's bound is not wrong by a factor of a hundred. ∎

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

Take `t = ζ` a root of unity and a point `w` of `W` above it. Then `s = ∓iζ` is a root of unity, and the point `x = ψ₁(s, y)` has, on the model `y² = x³ + Ax + B` with integers `A, B`, the coordinate `X = a·A₀/(s − r) + b` with rational `a, b, r, A₀`. So `h(X) ≤ h(aA₀) + h(r) + h(b) + 2 log 2 = 110.21`; `h(X : Y : 1) ≤ 1.5·h(X) + 0.5·log(1 + |A| + |B|)`, because `|Y|² ≤ (1 + |A| + |B|)·max(1, |X|)³` at the infinite places and `|Y|² ≤ max(1, |X|)³` at the finite ones; `h₂ ≤ h + 0.5·log 3`; and by (P) `ĥ(x) ≤ 446.3`. On `E₂`: `448.4`. By the parallelogram law `ĥ(α(w)) ≤ 2·446.3 + 2·446.3` and `ĥ(β(w)) ≤ 4·448.4` (`θ` is an isomorphism and does not change `ĥ`). By (P) again, `h₂(α(w), β(w)) ≤ 4·446.3 + 4·448.4 + 2·236.45 = 4051.4`. There are infinitely many such points, so by (Z) `h₂(C) ≤ 2·96·4051.4`. ∎ (In the program every quantity that enters a bound is a rational number, computed exactly: a logarithm of an integer is replaced by its number of bits times `0.6931471806 > log 2`, `log 3` by `1.0987`, and the constants of the paper are taken exactly as the decimal numbers printed there (whether they are the right numbers is the question of section 5, not a question of rounding). The comparison of the bound with the range of the sieve is made on exact rationals. The height of `G` enters through the lower bound 20.002899 of Lemma 1. That enclosure itself uses two constants (1.946 and 2.14) and a table from a paper of Cremona, Prickett and Siksek, also read by an automated reader.)

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
6. **For the Corollary only**: that a solution has class `d`. Section 11 gives this as Theorem D. To check there: **Lemma S2**, the support step, with its three cases (it uses unique factorisation in `ℤ[i]` and is on paper only; Lean has the common divisors of Lemma S3 and nothing else of that section); the identity for `A·B` in the setting; that the local programs test what part (c) says; and the thresholds of section 14. Before section 11 the class came from the two descents of `GAUSS_DESCENT.md` and `SECOND_DESCENT.md`, which say the same.

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

**`G` has infinite order for every rational `τ` other than `0, ±1`**, that is, for every ratio. Proof, with the computations of `degenerate_section.py` over `ℚ(τ)`. (i) *The model.* The coefficients of the model are polynomials in `τ` with rational coefficients, and so are the coordinates of `G`; the program checks symbolically that `P` lies on the curve and computes `G` by the chord through `P` and `σP`. For a ratio `m : n` the curve with `d = 2(m² − n²)` is the one for `τ = m/n` up to the factor `n⁶`, a square, in `d·Q₁`. (ii) *Specialisation.* The discriminant of the model is a constant times `τ²(τ − 1)⁶(τ + 1)⁶(τ² + 1)²` times the squares of two quartics without rational roots, so for a rational `τ₀ ≠ 0, ±1` the fibre is an elliptic curve, `G(τ₀)` is a point of it with finite coordinates, and specialisation at `τ₀` is a homomorphism from the points over `ℚ(τ)` whose multiples used below are defined at `τ₀`; a multiple `jG` specialises to `O` exactly when its `x` has a pole at `τ₀`. (iii) *Torsion.* A rational torsion point has order 1 to 10 or 12 (Mazur). If `G(τ₀)` had order `2j`, then `jG(τ₀)` would be a point of order 2, finite, with `y = 0`; if it had order `2j + 1`, then `jG(τ₀)` and `(j+1)G(τ₀)` would be finite with the same `x`. In each case the numerator of a rational function of `τ` (`y(jG)`, resp. `x((j+1)G) − x(jG)`) vanishes at `τ₀`. (iv) *The ten polynomials.* These numerators, of degree 9 to 162, are written out in factored form in `degenerate_section_polys.txt` (6 kB; Sage and PARI give the same factor degrees): apart from `τ`, `τ ± 1` and `τ² + 1` their factors are irreducible of degree 6 or more, so they have no rational root other than `0, ±1`. This proves rank at least 1; it says nothing about rank exactly 1, which Theorem A needs and which is computed curve by curve. For (627,364) the formula gives the point `G` of Lemma 1 (same height, 20.006653). So `E₁` has rank at least 1 for every ratio. That explains an observation of `SECOND_DESCENT.md` (this curve never has rank 0: 74 ratios of 74 below 1000 and 40 of 40 between 1000 and 3000), and it says where the effective statement below applies: exactly on the fibres of rank 1 of this family, with a point of infinite order that costs nothing to find.

*What this says about the problem.* A magic square of squares with ratio `m : n` gives a rational point of `X` (two of the conditions "product of two forms is a square"). At a rational point of `X` the numbers `N(A(t)) = GL·GU/4` and `N(B(t)) = F0·F8` are squares, so `A(t)·B(t)` has square norm, and an element of `ℚ(i)` with square norm is a rational number times a square: **the class `d` of `A(t)·B(t)` is a rational number for every rational point of `X`**, whatever the type of the ratio, and the point lifts to the double cover with that `d`. So the curves of this file are the twists of one unramified double cover of the classical curve `X`, and the two descents of this repository are used, for the Corollary, only to say which twists `d` can occur. That is a much smaller statement to check than the whole Gaussian descent on five forms: a 2-descent on one double cover. It has not been written out separately, and it is the natural next step.

**Corollary of Theorem B and of [CVV] (effective, rank 1).** For every `m, n, d` such that `E₁ = Jac(Y² = d·Q₁)` has no complex multiplication and rank at most 1 over `ℚ`, the set of rational `t` with `d·Q₁(−it)` and `d·Q₂(−it)` both squares in `ℚ(i)` is finite and can be determined: the points `(A, B)` have height at most the bound of Lemma 5, computed by `closure_estimate.py` for the given numbers, and a sieve covers that range. The statement and its proof are those of Theorem A, with the same lemmas; what depends on the numbers is computed by the same programs (they take the degenerate point as generator, so as written they need a ratio of square type and that point of infinite order). Nobody has read this either.

## 10. Twelve more ratios, between 1000 and 3000

`SECOND_DESCENT.md` leaves 67 ratios with `1000 < m ≤ 3000` (two descents, then rank 0 from descent bounds or `L`-values). Of these, 40 are of square type, and for 12 of them `E₁` has rank exactly 1 (`two_map_scope_wide.txt`; PARI's bounds, not checked with Magma; the list of 67 is regenerated from the tallies by `wide_open.py`, which also computes the partners and the counts below, `route_wide_summary.txt`). On each of the 12 the whole route gives the same answer as for (627,364): the degenerate point generates, the bound of Lemma 5 (13,651 to 14,132) is covered by the sieve, only `(1, 1)` is left within it, and the fibre is `t = ±1` (`route_wide.txt`):

    (1064,423) (1332,1045) (1395,532) (1785,688) (2015,1632) (2107,276) (2108,195) (2135,1248) (2331,1300) (2752,825) (2765,2508) (2867,756).

Six of them have their partner in the list of 67: (1487,641), (2377,287), (1927,863), (2473,1097), (2383,1831), (2303,1913).

*The twist.* The two conditions do not change when `d` is replaced by `−d`, because `−1` is a square in `ℚ(i)`; but the curves over `ℚ` do: `Y² = −d·Q₁` has Jacobian `E₂`, and the two curves of section 2 change places. So the whole argument, and every program, can be run with `−d`, and it then needs `E₂` of rank 1 instead of `E₁`; the point of the degenerate solution is again a point of infinite order, now on `E₂`. For three of the 40 ratios `E₂` has rank exactly 1 while `E₁` does not (`rank_pairs_wide.txt`): (2365,588) and (2405,1428), where `E₁` has rank 2, and (2640,511). On each the route with `−d` leaves only `(1, 1)` and `t = ±1` (`route_twist.txt`). The partner of (2365,588), (2953,1777), is in the list.

So, if the two descents and this route are right, 22 of the 67 are closed and 45 are left with `1000 < m ≤ 3000` (`route_wide_summary.txt`); with (627,364) and (991,263), none is left with `m ≤ 1000`. Where `E₁` and `E₂` both have rank 2 or more (10 of the 40 for certain) nothing here applies. For 15 other ratios of square type PARI's first bounds for the rank of `E₁` are 1 and 3; with more effort one of them, (2257,1176), is found to have rank 3, and 14 stay undecided (`rank_effort.py`, `rank_effort_wide.txt`; not rerun in CI). These are results of the experiment, not of the census.

## 11. The class `d` without the two descents

The Corollary of section 1 used the two descents of this repository for one thing: that a solution has class `d = 521266`. Section 9 said that this is a 2-descent on one double cover. Here it is, as one theorem and two short programs. The programs (`scripts/cover_descent/cover_descent.py` and, by a different method, `cover_bruteforce.py`) compute with Python integers; the first also uses the library sympy to factor integers and to list primes. No polynomial algebra, number field or elliptic curve is in them.

**Setting.** Let `m > n > 0` be coprime integers and `u, v` coprime integers with `F0, F8, GL, GU` all squares of integers (a magic square of squares with ratio `m : n` gives such a pair, with `t = u/v`; `G4` is not used). In `ℤ[i]` put

    γ0 = mv + inu,   γ8 = nv − imu,   γL = (1 + i)((mu − nv) + i(mv − nu)),   γU = (1 + i)((mu + nv) + i(mv + nu)),
    N(γ0) = F0,   N(γ8) = F8,   N(γL) = GL,   N(γU) = GU,
    A = −Su² + 2Diuv + Sv²,   B = γ0·γ8,   B′ = conj(γ0)·conj(γ8),   A·B = γ0·γ8·conj(γL)·conj(γU) / (1 + i)²,

with `A, B` as in section 9, made homogeneous (`S = m² + n²`, `D = m² − n²`). The four forms are positive, because `m ≠ n` and `(u, v) ≠ (0, 0)`.

**Theorem D (the class of a solution).** In this setting:

- (a) there is exactly one squarefree positive integer `d` with `A·B = d·z²`, `z ∈ ℚ(i)`; and `A·B′` is `d` times a square as well;
- (b) every prime that divides `d` divides `2mn(m² − n²)`;
- (c) for every prime `p`, `d` is in the set `Loc_p` of the squarefree positive integers `d′` for which some point `(u′ : v′)` of the projective line over `ℚ_p` has `F0, F8, GL, GU` squares in `ℚ_p` and `d′·A·B`, `d′·A·B′` squares in `ℚ_p ⊗ ℚ(i)`;
- (d) `t = u/v` satisfies the two conditions of Theorem A with this `d`.

So: **a rational solution gives one of finitely many `d`, and that `d` passes every local test.** For (627,364) there are 256 candidates by (b), and the programs say that only `521266` is in `Loc_p` for all of 2, 3, 11, 13, 17, 37, 53.

*Lemma S1 (a number of square norm is a rational number times a square).* Let `g ∈ ℚ(i)`, `g ≠ 0`, with `N(g) = w²`, `w` rational. Then `g = d·z²` with `z ∈ ℚ(i)` for exactly one squarefree positive integer `d`. *Proof.* If `g = −w`, `g` is rational. If not, put `z₀ = g + w ≠ 0`; then `z₀² = g² + 2gw + g·conj(g) = g·r` with `r = g + conj(g) + 2w` rational and not zero, so `g = r·(z₀/r)²`. A rational number is `±d·q²` with `d` squarefree and positive, and `−1 = i²`. If `d·z² = d′·z′²` then `d/d′ = (a + bi)²` is rational, so `ab = 0` and `d/d′` is `a²` or `−b²`; it is positive, so it is the square of a rational number, and `d = d′`. ∎

*Lemma S2 (support; unique factorisation in `ℤ[i]`).* Let `γ ∈ ℤ[i]`, `γ ≠ 0`, with `γ = d·z²`, `d` a squarefree positive integer, `z ∈ ℚ(i)`. If an odd prime `p` divides `d`, then some prime `π` of `ℤ[i]` above `p` divides both `γ` and `conj(γ)`. *Proof.* `ℤ[i]` is a unique factorisation domain; write `v_π` for the exponent of a prime `π`. Then `v_π(γ) = v_π(d) + 2·v_π(z)`.
- *`p ≡ 3 mod 4` (inert).* `p` is a prime of `ℤ[i]` and `v_p(d) = 1`, so `v_p(γ)` is odd, hence at least 1: `p` divides `γ`, and `p = conj(p)` divides `conj(γ)`.
- *`p ≡ 1 mod 4` (split).* `p = π·conj(π)` with `π`, `conj(π)` not associated, and `v_π(d) = v_conj(π)(d) = 1`. So `v_π(γ)` and `v_conj(π)(γ)` are both odd, hence both at least 1: `π` divides `γ`, and `conj(π)` divides `γ`, that is, `π` divides `conj(γ)`.
- *`p = 2` (ramified)* is excluded by the hypothesis; `2` is always among the primes allowed in (b). ∎

*Lemma S3 (common divisors; Lean).* For `u, v` with `au + bv = 1`: a common divisor in `ℤ[i]` of `mv + inu` and its conjugate divides `2mn`; the same for `nv + imu`; a common divisor of `(mu − nv) + i(nu − mv)` and its conjugate divides `2(m² − n²)`; the same for `(mu + nv) + i(mv + nu)` (`support_alpha0`, `support_alpha8`, `support_gammaL`, `support_gammaU` in `GaussDescent.lean`). **These four statements are the only part of this section that is checked in Lean.** Lemmas S1 and S2, the valuations, and everything local are not; the header of the Lean file says the same.

*Proof of Theorem D.* (a) `N(A·B) = F0·F8·GL·GU/4` is the square of a rational number and is not zero, so S1 gives `d`. And `(A·B)·(A·B′) = A²·F0·F8` is a square, so `A·B′` has the same `d`. (b) Each of `γ0, γ8, γL, γU` has square norm, so by S1 it is `d_j` times a square, with squarefree positive integers `d0, d8, dL, dU`. Let `p` be an odd prime dividing `d0`. By S2 a prime `π` above `p` divides `γ0` and its conjugate, by S3 it divides `2mn`, so `N(π)`, which is `p` or `p²`, divides `(2mn)²`, and `p` divides `2mn`. The same for `d8`; for `dL` and `dU` (`π` is odd, so it divides the form without the factor `1 + i`) `p` divides `2(m² − n²)`. By the last identity of the setting, and because conjugation does not change a rational class and `(1 + i)²` is a square, `d` is the squarefree part of `d0·d8·dL·dU`. (c) The solution is such a point, for every `p`. (d) `Q₁(−it) = A·B/v⁴` and `Q₂(−it) = A·B′/v⁴`, so both are `d` times a square by (a). ∎

**The local test.** For a prime `p`, the program `cover_descent.py` covers the projective line over `ℤ_p` by discs `t₀ + p^k ℤ_p` in the two charts `(t, 1)` and `(1, t)`. On a disc a value is congruent to its value at `t₀` modulo `p^k`, so its square class is constant once its valuation is small against `k` (the exact thresholds, and why they are enough, are in section 14); a disc is dropped when a constant class is wrong, kept when all are constant and right, and cut into `p` discs otherwise, down to a fixed depth, where everything not yet constant is counted as possible. The set it returns contains `Loc_p`; it can be too large but not too small.

**The same by exhaustion.** `cover_bruteforce.py` shares no code with the first program and uses no valuation, residue symbol or Hensel lemma. For `M = p^N` it goes through every point of the projective line over `ℤ/M` and keeps a candidate `d` when at some point the four forms are in a list of all squares modulo `M` and `d·A·B`, `d·A·B′` are in a list of all squares of `ℤ[i]` modulo `M`. A candidate that is not kept is not in `Loc_p`: if it were, the six values at a primitive `(u, v)` would be squares of integral elements, hence squares modulo `M`.

**Result for (627,364)** (`cover_627_364.txt`, `cover_certificate_627_364.txt`). Of the 256 classes, 64 pass at 2, then 32 (at 3), 16 (at 11), 8 (at 13), 4 (at 17), 2 (at 37) and 1 (at 53): **`d = 521266`, the class of the degenerate solution, and no other.** The certificate file lists, for each prime up to 61, the candidates removed there, from both programs. They remove the same candidates at 17 of the 18 primes; at 5 the exhaustion modulo `5⁵` removes half of the candidates and the first program reaches its depth and removes none. Both end with `521266` alone.

**Corollary, in its new form.** There is no magic square of squares with ratio `627 : 364` or `991 : 263`, given:

1. the five forms and the partner map (Lean);
2. Theorem D, proved on this page. Its part (b) uses unique factorisation in `ℤ[i]` (Lemma S2), which is not in Lean, and the four Lean lemmas of S3;
3. the output of the local programs for part (c): seven primes, two programs that agree;
4. Theorem A, with everything it rests on: the rank bound and the generator (Lemma 1: PARI, Sage, Magma, the saturation), the genus and degree of the curve of solutions (Lemmas 2 and 3), the height bound from the theorem of Checcoli, Veneziano and Viada (Lemmas 4 and 5), the sieve (Lemma 6) and the fibre (Lemma 7).

**The new descent replaces the route to `d` and nothing else.** `GAUSS_DESCENT.md` and `SECOND_DESCENT.md` are no longer used for this pair; Theorem A is used exactly as before. What is still unread is listed in section 7.

**Checks of the programs.** The class of every known point passes at every prime, in every run (the first program stops otherwise). Random points: for (627,364) and 14 primes, 43,228 random `(u, v)` at which the four forms are squares in `ℚ_p`, with the class computed from the exact integers: none outside the allowed set, and `A·B` and `A·B′` always in the same class (`cover_selftest_627_364.txt`). The two programs give the same final set for (627,364), (820,799) and (889,320) (`cover_certificate_*.txt`). And the first agrees with the Gaussian descent wherever both say something, on the 157 ratios of the census list (`cover_34.txt`, `cover_123.txt`, `cover_summary.txt`):

- **The 34 ratios without a trivial point: no class is possible, for every one of them.** Among them are (820,799) and (889,320), which the census counts on Magma's two-cover descent alone (`cover_820_799.txt`, `cover_889_320.txt`: empty after the primes up to 61; the exhaustion gives the same, `cover_certificate_820_799.txt`, `cover_certificate_889_320.txt`). For these 34 ratios this is a complete argument of its own: Theorem D and the local programs, and nothing else.
- The 123 ratios with a trivial point: for 111 the only class left is that of the known points (65 of square type, 46 of the other type); for 12 one or three further classes are left, also with good primes up to 1500. This descent follows one class where the Gaussian descent follows five, so it can leave more; those 12 would need the further classes treated by ranks.
- Of the 24 ratios of section 6, 19 have exactly one class here, (627,364) among them. The finer descent of section 12 leaves exactly the known points for all 24.

Nobody has read this section or the programs. It is the smallest complete argument in the repository for a ratio that the census counts as open: for (820,799) or (889,320), a reader has to check Theorem D, one page, and either 150 lines of Python or the 60 lines of the exhaustion.

## 12. One class for each form, and every ratio of the list

The descent of section 11 follows one class, that of the product `A·B`. A finer one follows four (`scripts/cover_descent/tuple_descent.py`; Python integers, and sympy to factor integers and list primes): each of the four linear forms of section 11, with `γ_L, γ_U` multiplied by `1 + i` so that their norms are `GL, GU`, has square norm at a solution, hence a rational class `d0, d8, dL, dU`, supported on the primes of `2mn` for the first two and of `2(m² − n²)` for the others (Lemmas S1 to S3 of section 11, applied to each form; only S3 is in Lean). The local test is the one of section 11 for each form separately, with one more point: the forms are linear and the coprime integers `(u, v)` are `(t, 1)` or `(1, t)` times a `p`-adic unit, so at a prime the set of allowed tuples is closed under multiplying the four classes by the class of a unit. (The first version of the program forgot this. Its test with random local points failed 7,706 times in 43,228, which is how it was found; with the unit, 0 failures on four ratios, `tuple_selftest.txt`. The descent of section 11 uses forms of even degree and is not affected.)

Results (`tuple_*.txt`; under a second per ratio; the lists up to 1000 are rerun in the job `quick`, and the first lines of the two files up to 5000):

- **(627,364): of 262,144 tuples exactly two are left, those of `t = 1` and `t = −1`**, with product `d = 521266`.
- The 123 ratios with a trivial point: **exactly the tuples of the known points, for all 123.** The further classes that section 11 left for 12 of them are gone.
- The 34 ratios without a trivial point, and all 681 ratios with `m ≤ 1000` that avoid the 484 classes and are not of the two exceptional types: **no tuple is left, for every one.**
- The 15 ratios between 1000 and 3000 of section 10: exactly the known points.
- **Every ratio with `m ≤ 5000`** (`tuple_scan_5000.txt`, `tuple_trivial_5000.txt`, `pairs_trivial_5000.txt`). 21,104 coprime ratios with `m ≤ 5000` avoid the 484 residue classes. For the 19,581 with `m² + n²` neither a square nor twice a square (the 681 above and 18,900 with `1000 < m ≤ 5000`) **no tuple is left**, each time with good primes up to 150. For the 1,523 others (306 with `m ≤ 1000`, 1,217 above) **exactly the tuples of the degenerate solutions are left**; 129 of them need good primes above 150, up to 1500. About half an hour in all. So, if this program and the 484 classes are right, a magic square of squares with `m ≤ 5000` has one of 1,523 ratios, and there it has the class `d` of a degenerate solution: what remains for each of these ratios is the curve of section 9 at that one class.
- **The same for (627,364) by exhaustion** (`tuple_bruteforce.py`, 70 lines, no code shared, no valuations or residue symbols: every point modulo `p^N`, a unit, and a list of all squares of `ℤ[i]` modulo `p^N`): of the 262,144 tuples the same two are left, after the primes up to 41 (`tuple_bruteforce_627_364.txt`).

This is the conclusion of the Gaussian descent of `GAUSS_DESCENT.md` on these lists, obtained again by a different and much smaller program that follows rational classes of four forms (the Gaussian descent follows classes in `ℚ(i)` of five). The two agree everywhere. With it, **for all 24 ratios of section 6 and all 15 of section 10 the argument of this file does not use the two large descents**: it rests on Theorem D applied to each form (on paper, with four Lean lemmas inside it), this program, and the route of sections 2 to 8, which is used exactly as before (rank and generator by PARI, Sage and for (627,364) Magma, the theorem of Checcoli, Veneziano and Viada, the sieve). Nobody has read it.

*What holds for every ratio, and what does not.* For every ratio: the classes exist and are supported on the primes of `2mn(m² − n²)` (Theorem D); the curve of solutions with a given class has genus 9 and two maps to one elliptic curve (Theorem B); that curve has a point of infinite order (Theorem C). Ratio by ratio: which classes pass the local tests (a computation of under a second), and the rank of the elliptic curve. No uniform reason was found for the first: for ratios that are not of the exceptional types the classes die at small primes that change from ratio to ratio.

## 13. Two side remarks

**Why the bound for `a, b` is always near 14,000.** In the 39 ratios where it was computed the bound of Lemma 5 lies between 13,262 and 14,725. Both sides of the inequality grow like the logarithm of the ratio: the canonical height of the degenerate point is close to `3·log m` (fitted on 36 ratios: `2.995·log m + 0.67`; this slope is the height of the section `G` of Theorem C on the family, presumably exactly 3), and the bound of the theorem is close to `8.4·10⁸·log m`, because every height that enters it is a number of digits of a polynomial in `m, n`. The logarithms cancel, and `a² + b²` is bounded by about `8.4·10⁸/4.5`, that is `a, b ≤ 13,700` in the limit. So the number is not a property of the problem: it is 72.251 (the constant of the theorem) times the square of the degree 96 times crude bounds for the coefficients, divided by the height of one section. What it is good for: one range of the sieve serves every ratio, so the cost of a ratio does not grow with its size; and, if the bound were proved with one constant for all ratios (it is not: that needs an explicit lower bound for the height of the section in the whole family), every solution on a curve of rank 1 would be one of finitely many multiples of `G`, the same for all ratios.

**The perfect cuboid.** The same program has a second system (`tuple_descent.py cuboid m,n`): for a Pythagorean triple `(a, b, c)`, the three forms `t² + a², t² + b², t² + c²`, which are norms of `t + ia` and so on; a perfect cuboid with that face gives a point. The first descent works there too (the known points always pass, and the tests of section 14 have no failure), but it leaves more than for the magic square: 4 tuples for (3,4,5), 9 for (161,240,289), between 3 and 36 on the bases tried.

*A correction.* An earlier version of this paragraph described a second step on the curves `Y² = δ·(s − a)(s − b)(s − c)`, with `s = it`. That argument was not valid: the class of a product of three of the linear factors depends on the denominator of `t`. It had closed no base, and it is withdrawn.

*The second step that works* is in the repository [perfect-cuboid-census](https://github.com/drew2apps/perfect-cuboid-census) (`scripts/second_step.py`). With four factors the denominator goes away: `(s² − x²)(s − ey)(s − fz)` is `δ_y·δ_z` times a square of `ℚ(i)`, so the solution is a point over `ℚ(i)` on a curve of genus 1 defined over `ℚ`, and if its Jacobian or the twist by `−1` has rank 0 over `ℚ`, the conjugate of the point is `±P + T` with `T` of finite order, which leaves finitely many `t`. This is the argument of `SECOND_DESCENT.md`, in the form it takes there. **Of the 411 bases with `m ≤ 100` that the census by rank 0 leaves open, it closes 274**; the only rational `t` on all curves used are `0` and `∞`. A control with planted points (480 curves through real solutions of similar systems: none has a rank 0 with the right class, more than 220 with a wrong one) is there as well. It is an experiment, counted apart, and nobody has read it.

## 14. Appendix: when a square class is constant on a disc

The two descent programs of sections 11 and 12 decide square classes from congruences. Here are the rules they use, each with its proof. `p` is a prime, `k ≥ 1`, the disc is the set of `t ∈ ℤ_p` with `t ≡ t₀ mod p^k`, and `g` is one of the forms at `(t, 1)` or `(1, t)`: a polynomial in `t` with coefficients in `ℤ` or `ℤ[i]`. So `g(t) ≡ g(t₀)` modulo `p^k` on the disc, in `ℤ_p` or in `ℤ_p[i] = ℤ_p ⊗ ℤ[i]`. Write `g₀ = g(t₀)` and suppose `g₀ ≠ 0`.

- **R1 (`p` odd, values in `ℚ_p`; the programs ask `v < k`).** Let `v = v_p(g₀)`. If `v < k` then `g(t)/g₀ ∈ 1 + pℤ_p`, and every element `a` of `1 + pℤ_p` is a square (Hensel: `x² − a` has the simple root 1 modulo `p`). So `g(t)` has the class of `g₀ = p^v·ε`, which is given by `v mod 2` and by whether `ε` is a square modulo `p`.
- **R2 (`p = 2`, values in `ℚ_2`; `v + 3 ≤ k`).** If `v + 3 ≤ k` then `g(t)/g₀ ∈ 1 + 8ℤ_2`, and every such `a` is a square (Hensel in the form `|a − 1| < |2|²`: here `|a − 1| ≤ 1/8 < 1/4`). The class of `g₀ = 2^v·ε` is given by `v mod 2` and `ε mod 8`.
- **R3 (`p ≡ 3 mod 4`, values in `ℚ_p(i)`; `v < k`).** `ℚ_p(i)` is a field, unramified over `ℚ_p`, with integers `ℤ_p[i]` and residue field `F_p[i]` of `p²` elements; `v` is the smaller of the valuations of the two coordinates of `g₀`. If `v < k` then `g(t)/g₀ ∈ 1 + pℤ_p[i]`, a square as in R1. The class of `g₀ = p^v·ε` is given by `v mod 2` and by whether `ε` is a square in `F_p[i]`, that is, whether `ε^((p²−1)/2) = 1` there (`cover_descent.py`), or, which is the same because `ε^((p²−1)/2) = N(ε)^((p−1)/2)`, whether the norm of `ε` is a square in `F_p` (`tuple_descent.py`). A rational unit is always a square in `F_p[i]`, so a rational `d` has the label `(v_p(d) mod 2, square)`.
- **R4 (`p ≡ 1 mod 4`; `v < k` in each component).** Let `ι ∈ ℤ_p` with `ι² = −1`. Then `x + iy ↦ (x + ιy, x − ιy)` identifies `ℚ_p ⊗ ℚ(i)` with `ℚ_p × ℚ_p`, and an element is a square exactly when both components are. Each component of `g(t)` is congruent to that of `g₀` modulo `p^k`, and R1 applies to it. A rational `d` has both components equal to `d`. `tuple_descent.py` computes both components of each linear form; if only one of them is constant on the disc (the other one is close to its zero), a rational class on the disc can only be the class of that one, and at the largest depth the program takes that class in place of "every class is possible". `cover_descent.py` computes the first component of `A·B` and of `A·B′` and asks that both have the class of `d`; this is the same condition, because the second component of `A·B` is the first component of `conj(A)·B′ = N(A)·(A·B′)/A²`, and `N(A) = GL·GU/4` is a square in `ℚ_p` at the points considered. `ι` is computed by Newton's iteration modulo `p^(K+6)` (or `p^(K+4)`), with `K` the largest depth, and a value whose component has valuation too large for that precision is treated as not constant.
- **R5 (`p = 2`, values in `ℚ_2(i)`; `w + 6 ≤ 2k`).** Let `π = 1 + i`, so `2 = −iπ²`, and `w = v_π(g₀)`. Modulo `2^k` is modulo `π^(2k)`. If `w + 6 ≤ 2k` then `g(t)/g₀ ∈ 1 + π⁶ℤ_2[i] = 1 + 8ℤ_2[i]`, and every such `a` is a square (Hensel: `|a − 1| ≤ |8| < |2|²`). The class of `g₀ = π^w·ε` is given by `w mod 2` and the class of the unit `ε` modulo squares of units. *Why modulo 8 is enough for that:* two units `ε, ε′` have the same class exactly when `ε ≡ ε′·σ² mod 8` for some unit `σ` (if so, `ε/(ε′σ²) ∈ 1 + 8ℤ_2[i]` is a square; the converse is a reduction modulo 8). So the class of `ε` is the orbit of `ε mod 8` under the squares of the units modulo 8, and the function `key2` returns the smallest element of that orbit. (`w + 5 ≤ 2k` would do; the programs ask for more.)

**Why the allowed set cannot be too small.** A disc is dropped only when a class that is constant on it by R1 to R5 is not the required one; then no point of the disc is a local solution. At the largest depth, a value that is not yet constant is counted as possible. For `tuple_descent.py` the set is then closed under the classes of units (section 12). So each program returns a set that contains every class that occurs at a point over `ℚ_p`.

**Tests of these rules.** (1) `tuple_descent.py … oracle`: for 2⁶, 3⁴, 5³, 7³, 13², a list of all squares of `ℤ[i]` modulo `p^N` is made, and for about 34,000 values of the forms the rational classes `d` with `d·γ` in that list are compared with the label of the program: no disagreement on four ratios (`tuple_neartest.txt`). With a rule deliberately made wrong (units of `ℤ_2[i]` compared modulo 4; the residue symbol reversed) the same test reports 3,010 and 14,309 disagreements (`tuple_mutation.py`, `tuple_mutation.txt`), so it can fail. (2) `… neartest`: points `t* + p^j·r`, `j = 1, …, 6`, around the zeros `t*` of each form in `ℚ_p` and around infinity, where the classes change fastest: between 2,154 and 5,641 points for each of four ratios, none outside the allowed set. (3) The programs `cover_bruteforce.py` and `tuple_bruteforce.py` reach the same final sets by exhaustion over all points modulo `p^N`, without any of R1 to R5 (sections 11 and 12). *What the exhaustion found.* On the system of the perfect cuboid (section 13) it kept fewer tuples than `tuple_descent.py` at the primes 5, 13 and 29 for the base (161, 240, 289): 12 tuples against 33. Nothing was wrong, but the program was weaker than it had to be: near the zero of a form in one of the two components it counted every class as possible, where the other component already fixes the class (the last sentence on `tuple_descent.py` in R4). With that rule the allowed sets can only shrink (checked prime by prime on three ratios), the two programs give the same set at every prime up to 61 on the three bases compared, and the descent leaves 9 tuples on that base. For the magic square no result changes: every file was made again; the only differences are lines where good primes up to 150 are now enough (2 lines of `tuple_123.txt`, and 18 of `tuple_trivial_5000.txt` between its first and its last version).
