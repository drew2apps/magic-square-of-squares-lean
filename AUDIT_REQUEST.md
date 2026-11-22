# Request for an audit

For a reader who knows elliptic curves over ℚ. The Lean part of this repository is checked by the Lean kernel and needs no audit beyond reading the statements. The part described here is computed, was written with Claude (an AI model made by Anthropic), and has been checked by nobody else. I am not a number theorist. Each question below can be answered separately; an answer to one of them is useful.

## What is claimed

For coprime `1 ≤ n < m ≤ 1000` there are 425 ratios `m : n` for which congruences do not exclude a 3×3 magic square of squares with `(b + d) : (h + f) = m : n`. A solution is a rational `t` at which five quadratic forms in `t` are squares (`G1_FULL_SYSTEM.md`, "The system"). 268 ratios are closed by a quotient elliptic curve of rank 0. 121 more are closed by the Demjanenko–Manin argument on a quotient elliptic curve that occurs at least twice in the Jacobian and has rank exactly 1. This request is about those 121. `scripts/g1_full/dm_evidence.txt` has one line for each.

## Three closures, worked out

Heights are Sage's (twice Silverman's `ĥ`). `B` is Sage's `silverman_height_bound()`. `S = m² + n²`.

**1. Curve `E04`, pair (15,8)** (`python3 dm.py - 15,8`). `κ = 2Sn² = 36992`, `E: y² = x(x + 130050)(x + 36992)`. A solution `t` gives `P, P′ ∈ E(ℚ)` with `x(P) = κt²`, `x(P′) = κ/t²`, so `x(P)·x(P′) = κ²`. Rank 1 (PARI `ellrank` gives `[1, 1]`; eclib and Magma agree). Torsion of order 4. Generator `G` with `x(G) = −83521 = −289²` (eclib, `gens(proof=True)`), `ĥ(G) = 3.1576`, `B = 29.242`. Write `P = aG + T`, `P′ = a′G + T′`. From `|h(x(Q)) − ĥ(Q)| ≤ B` and `|h(x(P)) − h(x(P′))| ≤ 2 log κ`: `|a² − a′²|·ĥ(G) ≤ 2B + 2 log κ = 79.52`. So either `|a|, |a′| ≤ 25`, or `|a| = |a′|`, and then `P′ = ±P + T″` is an equation for `x(P)` for each torsion point `T″`. The script lists `|a| ≤ 27`, solves the equations, and gets 111 candidate abscissas. The only rational `t` at which the five forms are squares are `t = ±1`, which give degenerate squares.

**2. Curve `D1`, pair (631,71)** (`python3 dm2.py --upper=1 631,71`). `E: Y² = x³ + c1x² + c0c2x + c0²c3` with `c3 = 4S²m²n²`, `c2 = −16Sm³n³`, `c1 = 4S²(m² − n²)²`, `c0 = −16Smn(m² − n²)²`. With `z = t + 1/t`, a solution gives `P, P′` with `x(P) = c0/z`, `x(P′) = −c0/z`, so `|a² − a′²|·ĥ(G) ≤ 2B`. `G` is the point with `x = 0`; `ĥ(G) = 12.5482`, `B = 121.374`, so `|a|, |a′| ≤ 19` or `|a| = |a′|`; the script lists `|a| ≤ 21`; torsion of order 2; 43 candidates; only the degenerate `t = 0`. That `G` generates: PARI `ellsaturation` at the primes up to 9; the index is at most 8 because `ĥ ≥ 0.1816` on non-torsion points (Sage, `height_function().min(0.01, 20)`). Rank: PARI gives only `[1, 3]`; Magma's `TwoPowerIsogenyDescentRankBound` gives 1. This line and its partner (351,280) are the only two that rest on Magma alone.

**3. Curve `E0L`, pair (423,332)** (`python3 dm3.py 423,332`). `C: y² = F0(t)·GL(t)` with the rational point `t = m/n`; a solution `t` gives points of `C` at `t` and at `−t`. `t` is not an even function on the curve, so with `σ = 1/(t − m/n)` the script uses `|2h(σ(Q)) − ĥ(Q) − ĥ(R0 − Q)| ≤ K = C1 + log 2 + 2B`, where `R0` is the second point over `t = m/n` (derivation in `G1_FULL_SYSTEM.md`, "The bound for `E0L` in symbols"). Here `ĥ(G) = 25.5397`, `B = 42.530`, `C1 = 94.27`, `K = 180.02`, `R0 = G + torsion`, and the proved bound is `|a| ≤ 8.08`; the script lists `|a| ≤ 11`, 92 points, none of which gives a square. Saturation at the primes up to 3 (`ĥ ≥ 6.3039`).

## Questions

1. **Rank.** Are PARI's `ellrank` (2-descent with the Cassels pairing, run on the curve and on the isogenous curves) and Magma's `TwoPowerIsogenyDescentRankBound` unconditional upper bounds, as used here? For (631,71) and (351,280), is Magma's answer alone enough, and how would you confirm it?
2. **Generators.** Is the value returned by Sage's `height_function().min(0.01, 20)` a proved lower bound for `ĥ` on all non-torsion points of `E(ℚ)`, so that saturation at the primes up to `⌊√(ĥ(P)/λ)⌋` proves that the saturated point generates?
3. **Silverman's bound on both sides.** The software documents `B` as a bound for `h − ĥ`. `scripts/g1_full/height_sides.py` argues from Silverman's Theorem 1.1 that `ĥ − h ≤ B − h(j)/12 + 0.218`, and checks this for all 121 models. Is the theorem quoted correctly there, and does it apply to the integral models of `E04` and `D1`, which are not minimal?
4. **The comparison for `E0L`** in `dm3.py`: sum and product of the two abscissas as quadratics in `σ`, the constant `C1`, and the resulting bound on `|a|`.
5. **Floating point.** Heights are in double precision, and the proved bound on `|a|` is multiplied by 1.02 and increased by 2 (by 3 for `E0L`) before listing. Is that margin acceptable, or is interval arithmetic needed?
6. **The singular fibres** (`G1_FULL_SYSTEM.md`, Result 1): the character-sum bound `32N ≥ 2p − 100 − 270√p` and the count of roots behind it.
7. **The density sketch** (Result 1b).

## What already has an automatic check

`check_evidence.py` re-tests every recorded `t`, checks that every recorded generator lies on its curve and that the `E0L` models are the right curves, and compares the table with the JSON output of the scripts. `height_sides.py` is described in question 3. `eclib_check.py` and `magma_check.py` are the second and third programs for the rank bounds. None of these checks a height, a rank or a saturation.

## How to answer

Please open an issue at https://github.com/drew2apps/magic-square-of-squares-lean/issues. "This step is wrong" and "this step is fine" are both welcome, with or without a reference.
