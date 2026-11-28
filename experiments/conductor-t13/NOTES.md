# (627,364): a point search for rank 2 of the corner curve over ℚ(i) (build agent t13)

**Written by an AI model (build agent t13), not read by a person.** Tier D (experiment). Not reviewed by a
separate agent. No census count changes. 27 November 2026 (UTC).

## Bottom line

**Rank E(K) ≥ 2 is not proved: no point of infinite order was found.** The 2-coverings of every non-torsion class
of `Sel²(E/K)` were searched on `E = L1·L3·M1·M2` and on `E3 = E/⟨T3⟩`, the curve of the class with the smallest
expected regulator, with no hit:

- about 1.7·10¹⁰ pairs `(t1, t2) ∈ ℤ[i]²` with `h_p(t) ≤ 5·10⁴` on the six coset coverings;
- on `E3` also the shell `5·10⁴ < h_p(t) ≤ 10⁵`, about 2.8·10¹⁰ more pairs;
- on the other minimal models at the large primes, 46 models with `h_p ≤ 10⁴`, about 5.4·10⁹ pairs;
- 5.04·10¹⁰ pairs in all (see the table).

The same pipeline certifies rank ≥ 2 on 7 of 7 control curves of rank 2 over `K`, from points it finds itself. On
the torsion classes of the target it finds the torsion points. So the elliptic-Chabauty route over `K` is **not
retired** by this work, and (627,364) is not closed.

## The question: 0, 1 or 2?

From the repository's numbers (`KCTP_627.md`, `kctp_recheck.txt`):

- `dim Sel²(E/K) = 4`;
- `E(K)_tors = ℤ/4 × ℤ/2`, so its image in `E(K)/2E(K)` has dimension 2.

Hence `r + dim Ш(E/K)[2] = 2`. The zero Cassels–Tate pairing on `Sel²` gives no further bound.

- `r = 1` would need `dim Ш[2] = 1`. That is impossible if `Ш[2^∞]` is finite: the pairing is then alternating and
  non-degenerate, so `Ш[2^∞] ≅ A × A` and its 2-rank is even.
- With `Ш[2^∞]` finite, `r ∈ {0, 2}`, and PARI's root number +1 (a guide) agrees.
- **Unconditionally `r ∈ {0, 1, 2}`: the bound does not exclude 1.** Rank 1 would leave elliptic Chabauty over `K`
  open in principle; rank 2 closes it; rank 0 would close the ratio.
- Only exhibited points can give `r = 2` without a hypothesis. Two independent points are needed: one point gives only
  `r ≥ 1`.

## What was computed

### `covers.py`: the 2-coverings as binary quartics over ℤ[i]

Exact, Sage over `K`; about 40 s.

- For a class `(d1, d2, d3)`, take a point of the conic `d1 u² − d2 v² + (e1 − e2) z² = 0` (the lattice finder
  `kctp_recheck.conic_point`, checked exactly). It parametrises the conic.
- The covering is then `W² = G(t) = d3 (d1 u(t)² + (e1 − e3) z(t)²)`.
- `G` is made integral and minimised at every prime of `2·Π(e_i − e_j)·d1 d2 d3`, and at every prime of the ideal
  `(I, J)`. The moves are invariant-preserving `π²`-moves at multiple roots mod `π`, and division by `π²`.
- `G` is then reduced: its covariant point in hyperbolic 3-space is moved into the Picard fundamental domain. *(Run 12, after review referee_run12.)* Two of the ten models are not reduced in this sense: for `E` "coset B" (`T` about 2.2·10⁴) and `E3` "torsion t1" (`T` about 3.5·10⁵) the floating-point minimiser's error in `z` is about 10⁻³ to 5·10⁻², the loop alternates translations by ±1 and stops at its cap of 500 moves, so the committed models and covariants are the state at the cap. The searches use the stored models and covariants, so the regions searched are the ones stated; only the word "reduced" is wrong for these two. The tolerance of the moves is `10⁻⁶` when `covers.py` runs (since run 12; one covariant lies on the corner `1/2 + i/2`, to `10⁻⁷`) and `10⁻⁹` when `altmodels.py` and `controls.py` import it, as when their outputs were made.
- The change of variable is kept, and `G(t) = s² G0(M t)` is asserted exactly.
- After minimisation, the odd part of `N(4I³ − J²)` is the same for every class of a curve:
  - on `E`: `3⁶·5¹⁶·29⁸·37⁴·181⁴·263⁴·991⁴·23490241⁴·395249137273⁴`;
  - on `E3`: `3⁶·5³²·29¹⁶·37²·181²·263⁸·991⁸·23490241²·395249137273²`.
- The 2-part varies between `2⁴` and `2²⁸`: the minimisation at `1 + i` is incomplete (open item below).

### `search.py`: the sieve

numpy only.

- **Region.** All `(t1, t2)` with `h_p(t) = |t1 − z t2|²/T + T|t2|² ≤ H`, where `(z, T)` is the covariant. This is
  the ellipsoid adapted to `G`, with one `t2` per unit class.
- **Sieve.** 44 split prime ideals of norm below 1000.
- **Exact test.** On the survivors, `G(t1, t2)` must be a square in `ℤ[i]` (Python integers).

### `certify.py`: back to `E`

Exact, Sage.

- Each hit is mapped back to `E`, or to `E3` and then to `E` by the dual isogeny (checked on the points of order 4).
- The point is checked on the curve, and its Kummer class is checked against the covering's class.
- Torsion is decided exactly. Independence goes through `kctp_points.certify`.
- `hints_t13.json` holds the points found. It is empty.

### `controls.py`

- 4 base changes of rank 2 over `K` (the controls of `kctp_points.py`).
- 3 curves with `j ∉ ℚ` and two independent points from a naive search.
- On each, the pipeline (Sel², cosets, quartics, search at `H = 2000`, map back, certify) finds points of infinite
  order on all three cosets and **certifies rank ≥ 2: 7 of 7**.
- On the controls, the smallest points sit at `h_p ≈ 1` to `10`, with canonical heights 0.7 to 4.5 (Sage's
  normalisation).

### Results on the target (no exact hit anywhere)

| curve | coset | T (covariant) | region | pairs | hits |
|---|---|---|---|---|---|
| E | A (991i, 991i, 1) | 1.18e4 | h_p ≤ 5e4 | 2,500,560,299 | 0 |
| E | B (145, 145, 1) | 2.17e4 | h_p ≤ 5e4 | 2,375,919,282 | 0 |
| E | A+B (143695i, 143695i, 1) | 976 | h_p ≤ 5e4 | 3,038,009,691 | 0 |
| E3 | A | 1.07 | h_p ≤ 5e4, then 5e4 < h_p ≤ 1e5 | 3,084,208,492 + 9,252,710,421 | 0 |
| E3 | B | 5.53 | h_p ≤ 5e4, then 5e4 < h_p ≤ 1e5 | 3,084,026,639 + 9,252,542,040 | 0 |
| E3 | A+B | 1.00 | h_p ≤ 5e4, then 5e4 < h_p ≤ 1e5 | 3,084,209,985 + 9,252,709,706 | 0 |

**Other minimal models** (`altmodels.py`). These are the invariant-preserving `π²`-moves at the odd primes of norm
≥ 10⁴ (the two split primes of norm 23490241 and 395249137273, and the inert primes 263 and 991), taken one prime at a
time and at the two largest primes together. Each model is reduced again, and models with equal coefficients up to
sign are dropped. Some survivors are still equivalent under diagonal unit changes, so the counts overstate the
distinct models.

- **Control.** On `E3`'s torsion covering `t1`, the base model (`T = 3.5·10⁵`) shows one of the four torsion
  `t`-values at `H = 2000`. The alternative model `m01` (`T = 0.56`) shows all four.
- **Shape.** The coset coverings of `E` had `T ≈ 10³` to `2·10⁴` on the base model. Their alternative models have
  `T ≈ 1` to `3`.
- **Search.** With the base models excluded, there are 30 models on `E` and 16 on `E3`. All were searched with
  `h_p ≤ 10⁴` (about 1.23·10⁸ pairs each), and on `E3` also the shell `10⁴ < h_p ≤ 3·10⁴`, for two models only. **No exact hit.**
  - The shell was started on all 16 models of `E3` and stopped by the time box after two: `coset A m01`, which is the
    base model turned by `t1 → i t1` and so already covered to `10⁵`, and `coset B m01`. These are 9.9·10⁸ pairs
    each, no hit, recorded in `logs/` only (no JSON written).

**Total** in the committed outputs: 50,351,347,885 pairs `(t1, t2)` over the coset coverings (`search_*coset*.json`),
with no exact hit. All torsion-class hits, on the base and the alternative models, map back exactly to torsion
points of `E` with the right Kummer class (`certify.txt`).

### What the regions mean (a guide, not a bound)

- **Reach.** A covering parameter `t` with `h_p(t) ≈ H` has height `h(t) ≈ ½ log H`. A 2-covering in its quartic
  model gives `ĥ(P) ≈ 4 h(t)` (PARI normalisation, absolute), so `ĥ ≈ 2 log H`. This is theory up to `O(1)`, not
  measured; on the controls the smallest points sit at `h_p` = 1 to 10. The searches reach roughly 22–23 on the base
  models (about 45 relative to `K`) and about 18 on the alternative models, up to the local corrections in the
  next item.
- **The BSD scale** of `kctp_isog.txt` (a guide) puts `Reg·|Ш|` at about `4.6·10⁵·L''(E/K,1)/2` on `E`, and a
  quarter of that on `E3`. With `Ш = 1` and balanced generators, each generator would have height of order a few
  hundred relative to `K`. The search could only have succeeded if `L''` were small or the generators very
  unbalanced. The negative result is what that scale predicts.
- **The naive height is not uniform on these coverings.** `E` has multiplicative reduction with large Tamagawa
  numbers: I4 at the primes of norm 23490241 and 395249137273, and I8 and I4 above 5 and 29.
- **A point whose lift lies on a non-identity component there has an inflated naive height** (`calib_torsion.py`,
  `calib_torsion.txt`; on the `t1` coverings it lists nothing, because it skips the points with `x = e_k`). On the torsion
  covering `t2` of `E`, the torsion preimages have `h_p` = 0.84, 1.95 and 2.1·10¹¹. So a point of small canonical
  height can lie far outside the region on a given model.
- `altmodels.py` searches the other minimal models at the large primes (above) for this reason. It is still not a
  complete search.

## Relation to the separate session's search (`untrusted/explore5/explore5/kcurve/`, read, not modified)

- **The cosets agree.** That session's Sel² cosets on `E` have the same representatives `(d, d, 1)` as here,
  `d = 991i, 145, 143695i` up to torsion (its `selmer.log`).
- **Its searches.**
  - On `E' = E3`: `x' = d (U/W)²` for the six non-torsion `d`, box `|re|, |im| ≤ 280`, no hit.
  - On `E`: `X − e_j = b_j (U/W)²`, box 100, only the points of order 4.
  - In height, these reach `ĥ ≈ 2 log(280√2) ≈ 12` on `E'` and about 10 on `E` (PARI normalisation, absolute).
    Ours reach about 21–23. Under the uniform heuristic their regions lie inside ours.
  - The models differ (one coordinate of the degree-4 embedding against the quartic parameter). Points with large
    local corrections are therefore placed differently, and the two searches are complementary rather than
    redundant.
- **Its sieve passes a planted test** (its own caveat was that it had never been tested on a non-torsion point).
  - `search2.c` was compiled outside the repository (`gcc -O2 -o /tmp/t13_ext/search2 …/search2.c`), source
    unchanged.
  - It was run on our control `y² = x(x − 6i)(x − 4 + 5i)` (`j ∉ ℚ`, rank ≥ 2 certified by `controls.py`) with
    `H = 10`: `./search2 10 tag 0 1 -4 -1 0 -6 -4 5`, and with `0 -1`.
  - It finds the planted points `x = −8` (`b = i`, `U = 2 + 2i`, `W = 1`) and `x = 4 + 3i` (`b = −i`, `U = 1 + 2i`).
  - Every hit of these runs, checked exactly with Sage, is a point of infinite order (2 distinct `x` per class, no
    false positive).
  - This one-off check is not committed as a script.

## What is not done

- **A 4-descent.** It would double the reach: canonical height about `4 log H`. In the form of Merriman–Siksek–Smart,
  one writes `t1 − θ t2 = ξ·z²` in `L = K[θ]/G(θ)`; here `L = K(√((e3 − e1)/d1), √((e3 − e2)/d2))`, a field of
  degree 8. That needs `L(S, 2)`: class group and S-units of a degree-8 field with large discriminant (PARI `bnfinit`
  under GRH would do as a finder).
  - A route without `L` is Cassels' conic method: the 4-coverings are `D_a ∩ {l_u l_w = λ1 □, l_v l_w = λ2 □}`, with
    the `λ` of `kctp_recheck.lift_data`. It needs a Riemann–Roch model of degree 4 in `ℙ³`.
  - Neither was started.
- **The minimisation at `1 + i`** (2-parts `2⁴` to `2²⁸` above) and a height-difference bound that would turn
  "`h_p ≤ H`" into "`ĥ ≤ X`" for every point.
- The four curves of the class without full 2-torsion over `K`.

## What a proof of rank 2 would mean (if it is ever found)

With two independent points, `rank E(K) = 2` exactly (bound of `Sel²`), unconditionally. Elliptic Chabauty over `K`
needs `rank E(K) < [K:ℚ] = 2`, so the corner-curve route over `ℚ(i)` could never close (627,364) and could be retired.
The same holds for the other three corner curves, since `KCTP_627.md` shows the four have one rank. Rank 0, which
would close the ratio by `kclose.py`'s torsion check, would be excluded too. This run proves neither.

## Files and reproduction

Times are wall-clock on a shared 2-core machine (load 4 to 8).

    python3 covers.py --write              # covers.txt, covers.json            (about 40 s)
    python3 controls.py --write            # controls.txt                       (about 15 min wall, 5 min CPU)
    python3 search.py --curve=E3 --H=2000 --classes=torsion --write             (seconds; also --curve=E)
    python3 search.py --curve=E  --H=50000 --classes=coset --write             (about 13 min)
    python3 search.py --curve=E3 --H=50000 --classes=coset --write             (about 15 min)
    python3 search.py --curve=E3 --H=100000 --H0=50000 --classes="coset A" --write   (coset A and A+B)
    python3 search.py --curve=E3 --H=100000 --H0=50000 --classes="coset B" --write
    python3 altmodels.py --write           # altmodels.txt, altmodels.json   (3 min CPU, 9 min wall)
    python3 search.py --models=altmodels.json --curve=E  --H=10000 --classes=coset --skipbase --write   (about 8 min)
    python3 search.py --models=altmodels.json --curve=E3 --H=10000 --classes=coset --skipbase --write   (about 5 min)
    python3 search.py --models=altmodels.json --curve=E  --H=2000 --classes=torsion --skipbase --write  (seconds; also E3)
    python3 certify.py search_E_H50000_coset.json search_E3_H50000_coset.json search_E3_H100000_from50000_cosetA.json \
        search_E3_H100000_from50000_cosetB.json search_E_alt_H10000_coset.json search_E3_alt_H10000_coset.json \
        search_E_H2000_torsion.json search_E3_H2000_torsion.json search_E_alt_H2000_torsion.json \
        search_E3_alt_H2000_torsion.json --write     # certify.txt, hints_t13.json (no points)
    python3 calib_torsion.py > calib_torsion.txt   # h_p of torsion preimages (seconds; a guide)
    python3 reach.py --H=50000 --write     # reach_H50000.txt (naive heights of x on the region's boundary; a guide)

All seeds are fixed (`random.Random(20261127)`). The outputs carry no times (they go to stderr).
Re-checked here: `covers.py` reproduces `covers.json` byte for byte, and the torsion searches reproduce their `.txt`.
`controls.py` and `altmodels.py` were run once (deterministic by construction, not re-run). `logs/` holds the run
logs with times.
