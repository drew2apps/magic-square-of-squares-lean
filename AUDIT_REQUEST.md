# Request for an audit

For a reader who knows elliptic curves over ℚ. The Lean part of this repository is checked by the Lean kernel and needs no audit beyond reading the statements. The part described here is computed, was written with Claude (an AI model made by Anthropic), and has been checked by nobody else. I am not a number theorist. Each question below can be answered separately; an answer to one of them is useful. This file is a working checklist for the repository, not part of a paper.

**Two kinds of answer.** For every number below there are two different questions: did the program return this value (anyone can rerun it), and does the algorithm behind it prove the bound unconditionally. The second is what I cannot judge. Please say which one an answer is about.

## What is claimed

For coprime `1 ≤ n < m ≤ 1000` there are 425 ratios `m : n` for which congruences do not exclude a 3×3 magic square of squares with `(b + d) : (h + f) = m : n`. A solution is a rational `t` at which five quadratic forms in `t` are squares (`G1_FULL_SYSTEM.md`, "The system"). 268 ratios are closed by a quotient elliptic curve of rank 0. 122 more are closed by the Demjanenko–Manin argument on a quotient elliptic curve that occurs at least twice in the Jacobian and has rank exactly 1. This request is about those 122, and about 9 more that are closed by a version of the argument with four points, for rank 2 and 3 (closure 4 below). `scripts/g1_full/dm_evidence.txt` has one line for each, and `dm_log.jsonl` the raw output of the scripts.

## Software

Python 3.13.18; passagemath 10.8.14 (Sage), installed with `pip install passagemath-schemes passagemath-eclib passagemath-symbolics cypari2`; PARI 2.17.4 through cypari2; eclib as shipped with passagemath-eclib; Magma 2.29-10 through the public calculator at magma.maths.usyd.edu.au/calc. Commands are run in `scripts/g1_full/`. Each script prints one JSON record per pair; the records below are the ones that were printed, copied from `dm_log.jsonl` without its bookkeeping fields (`of`, `run`, `method`, `sec`).

## Four closures, worked out

Heights are Sage's (twice Silverman's `ĥ`). `B` is Sage's `silverman_height_bound()`. `S = m² + n²`.

**1. Curve `E04`, pair (15,8).** Command: `python3 dm.py - 15,8`.

```
{"pair": [15, 8], "rank": [1, 1], "gens": "mwrank", "G": "-83521", "hG": 3.1576167237813184, "B": 29.241874456487658, "amax": 27, "ntors": 4, "ncand": 111, "t": {"1": "degenerate", "-1": "degenerate"}, "nt": 57, "status": "closed"}
```

`κ = 2Sn² = 36992`, `E: y² = x(x + 130050)(x + 36992)`, in Sage `EllipticCurve([0, 167042, 0, 4810809600, 0])`. A solution `t` gives `P, P′ ∈ E(ℚ)` with `x(P) = κt²`, `x(P′) = κ/t²`, so `x(P)·x(P′) = κ²`. Rank 1 (PARI `ellrank` gives `[1, 1]`; eclib and Magma agree). Torsion of order 4. Generator `G` with `x(G) = −83521 = −289²` (eclib, `gens(proof=True)`), `ĥ(G) = 3.1576`, `B = 29.242`. Write `P = aG + T`, `P′ = a′G + T′`. From `|h(x(Q)) − ĥ(Q)| ≤ B` and `|h(x(P)) − h(x(P′))| ≤ 2 log κ`: `|a² − a′²|·ĥ(G) ≤ 2B + 2 log κ = 79.52`. So either `|a|, |a′| ≤ 25`, or `|a| = |a′|`, and then `P′ = ±P + T″` is an equation for `x(P)` for each torsion point `T″`. The script lists `|a| ≤ 27`, solves the equations, and gets 111 candidate abscissas. The only rational `t` at which the five forms are squares are `t = ±1`, which give degenerate squares.

**2. Curve `D1`, pair (631,71).** Command: `python3 dm2.py --upper=1 631,71` (the option supplies the upper rank bound 1 from Magma, because PARI gives `[1, 3]`).

```
{"pair": [631, 71], "curve": "D1", "rank": [1, 1], "upper": "external: --upper=1 (Magma, see magma_check.txt)", "gens": "pari point, saturated at primes <= 9 (height lower bound 0.1816)", "G": "0", "hG": 12.548181437530385, "B": 121.37383378653512, "amax": 21, "ntors": 2, "ncand": 43, "nt": 1, "t": {"0": "degenerate"}, "status": "closed"}
```

`E: Y² = x³ + c1x² + c0c2x + c0²c3` with `c3 = 4S²m²n²`, `c2 = −16Sm³n³`, `c1 = 4S²(m² − n²)²`, `c0 = −16Smn(m² − n²)²`. With `z = t + 1/t`, a solution gives `P, P′` with `x(P) = c0/z`, `x(P′) = −c0/z`, so `|a² − a′²|·ĥ(G) ≤ 2B`. `G` is the point with `x = 0`; `ĥ(G) = 12.5482`, `B = 121.374`, so `|a|, |a′| ≤ 19` or `|a| = |a′|`; the script lists `|a| ≤ 21`; torsion of order 2; 43 candidates; only the degenerate `t = 0`. That `G` generates: PARI `ellsaturation` at the primes up to 9; the index is at most 8 because `ĥ ≥ 0.1816` on non-torsion points (Sage, `height_function().min(0.01, 20)`). Rank: Magma's `TwoPowerIsogenyDescentRankBound` gives 1; the request, as sent, is

```
m:=631; n:=71; S:=m^2+n^2; D2:=(m^2-n^2)^2; c3:=4*S^2*m^2*n^2; c2:=-16*S*m^3*n^3; c1:=4*S^2*D2; c0:=-16*S*m*n*D2; E:=EllipticCurve([0,c1,0,c0*c2,c0^2*c3]); E:=MinimalModel(E); G,mp:=TorsionSubgroup(E); T:=mp((Order(G.1) div 2)*G.1); hi:=TwoPowerIsogenyDescentRankBound(E,T); print "BOUNDS", 0, hi;
```

and it is in `magma_raw.txt` with the response. This line, its partner (351,280) and the line for (724,665) (run on its partner (1389,59)) are the only three that rest on Magma alone; the last one is also closed by the method of closure 4, with PARI's bounds.

**3. Curve `E0L`, pair (423,332).** Command: `python3 dm3.py 423,332`.

```
{"pair": [423, 332], "curve": "E0L", "rank": [1, 1], "model": ["0", "-1", "0", "-10810266634575227101441", "180302397457159147910876363714305"], "ntors": 4, "gens": "pari point, saturated at primes <= 3 (height lower bound 6.3039)", "G": ["-27869819136", "21446103082364927"], "hG": 25.539739392571548, "B": 42.52991460522845, "C1": 94.27115544209109, "K": 180.02413183310793, "rho": 1, "W": 7.576630582558059, "amax": 11, "nalg": 0, "npts": 92, "worst_observed": 35.40671163564093, "t": {}, "status": "closed"}
```

`C: y² = F0(t)·GL(t)` with the rational point `t = m/n`; a solution `t` gives points of `C` at `t` and at `−t`. `t` is not an even function on the curve, so with `σ = 1/(t − m/n)` the script uses `|2h(σ(Q)) − ĥ(Q) − ĥ(R0 − Q)| ≤ K = C1 + log 2 + 2B`, where `R0` is the second point over `t = m/n` (derivation in `G1_FULL_SYSTEM.md`, "The bound for `E0L` in symbols"). `model` is the minimal model on which `G` is given. Here `R0 = G + torsion`, and the bound computed from these values is `|a| ≤ W + 1/2 = 8.08`; the script lists `|a| ≤ 11`, 92 points, none of which gives a square. `worst_observed` is the largest value of the left side on the listed points; it must be below `K`.

**4. Curve `E0L` with four points, pair (840,83), rank 2.** Command: `python3 dm4.py 840,83`. The record, without the coordinates of the points:

```
{"pair": [840, 83], "rank_bounds": [2, 2], "rank": 2, "gens": "points saturated at primes <= 6 (height lower bound 6.7272, regulator of the points before saturation 918.282)", "gram": [[27.33573160470327, 13.55818659312613], [13.55818659312613, 40.31743168527588]], "rho": [1, 0], "BUP": 72.46748012016593, "U": 59.15875684694732, "lo4": -292.5078744823475, "up4": 54.92087403996096, "eL": 91.24383865062836, "eU": 28.51990772172707, "kappa": 47.329118695951244, "Hstar": 1053.0644658966871, "cstar": 24.464379015342047, "nA": 472, "observed N(a) - H": [-26.914238473697424, 13.52452656518119], "nz": 4, "nroots": 8, "t": {}, "status": "closed"}
```

A solution `t0` gives four rational points `Q1, …, Q4` of `C: y² = F0·GL`, above `t0, −t0, 1/t0, −1/t0`. With `a_i = Q_i − R0/2`: (i) `N(a_i) = ĥ(2Q_i − R0)/4` is within `[−eL, eU]` of `H = h(t0)`; (ii) `|⟨a_1, a_2⟩|` and `|⟨a_2, a_4⟩|` are at most `κ`. Both come from identities of the form `x(…) = (…)²` on the model `Y² = x(x − α)(x − β)`, from elementary bounds for the heights of the resulting fractions, and from the bounds `−U ≤ h(x) − ĥ ≤ BUP` (`G1_FULL_SYSTEM.md`, "Rank 2 and 3 on `E0L`"; the header of `dm4.py`). In rank 2 the Gram determinant of `a_1, a_2, a_4` vanishes, so `z = a_1 ∓ a_4` has `ĥ(z) ≤ c*(H)`. With `H* = 1053`: either `H < H*` and `Q1` is one of the 472 listed points, or `ĥ(z) ≤ 24.5`, which is below the smallest height of a non-torsion point (27.3), so `z` is one of the 4 torsion points, and `t(Q − z)·t(Q) = −1` is solved for each. No candidate passes the five-form test.

## Questions

1. **Rank.** Are PARI's `ellrank` (2-descent with the Cassels pairing, run on the curve and on the isogenous curves) and Magma's `TwoPowerIsogenyDescentRankBound` unconditional upper bounds, as used here? For (631,71), (351,280) and (724,665), is Magma's answer alone enough, and how would you confirm it?
2. **Generators.** Is the value returned by Sage's `height_function().min(0.01, 20)` a proved lower bound for `ĥ` on all non-torsion points of `E(ℚ)`, so that saturation at the primes up to `⌊√(ĥ(P)/λ)⌋` proves that the saturated point generates?
3. **Silverman's bound on both sides.** The software documents `B` as a bound for `h − ĥ`. `scripts/g1_full/height_sides.py` derives `ĥ − h ≤ B − h(j)/12 + 0.218` from two statements: at the finite places the naive local height is at least the canonical one, for any integral model (Cremona, Prickett and Siksek 2006, Proposition 8), and Silverman's bound at the real place (1990, p. 737, as quoted in their Theorem 15). It checks the inequality for all 122 models in exact arithmetic. Are the two statements used correctly in that script (normalisation of the heights, the real place of ℚ treated with the bound for complex points)? The first paper was read from its text with an AI model; Silverman's paper has not been seen, only the quotation in their Theorem 15.
4. **The comparison for `E0L`** in `dm3.py`: sum and product of the two abscissas as quadratics in `σ`, the constant `C1`, and the resulting bound on `|a|`.
5. **Heights of the generators.** The heights `ĥ(G)` used by the scripts are double-precision values from Sage. `scripts/g1_full/height_certify.py` replaces them by enclosures in exact rational arithmetic: `ĥ(G)` lies between `(h(x(2ᵏG)) − BUP)/4ᵏ` and `(h(x(2ᵏG)) + U)/4ᵏ`, where `2ᵏG` is computed exactly with the duplication formula (`k = 6` or `7`), `U` is as in question 3, and `BUP` bounds `h − ĥ` by the local bounds of Cremona, Prickett and Siksek (Proposition 8 and Table 1: `α_p` by Kodaira type and Tamagawa number, plus `ord_p(D/D_min)/6` for a non-minimal model) and Silverman's bound at the real place. With the lower ends of the enclosures the bound on `|a|` stays below the listed range in every line, by at least 9%. Is the use of Table 1 correct (in particular for a model that is integral and not minimal), and is the doubling enclosure the right way to certify a height?
6. **The singular fibres** (`G1_FULL_SYSTEM.md`, Result 1): the character-sum bound `32N ≥ 2p − 100 − 270√p` and the count of roots behind it.
7. **The density sketch** (Result 1b).
8. **The argument with four points** (closure 4): the three identities are checked by the program; the questions are the elementary height bounds in the header of `dm4.py` (the divisibility statements for the common factors are in Lean, `G1FourPoints.lean`; the size estimates and their combination are not), the use of the Gram determinant in rank 2 and of the block argument in rank 3, the completeness of the two listings, and the saturation bound from Hermite's constant. For (740,539) the second independent point came from Magma's four-descent with class groups under GRH; the point is checked to lie on the curve, so GRH should not matter. Is that right?

A question for a number theorist of any kind, not only for this audit: are Theorem 1, Theorem 1′, the Lemma and Theorem 3′ of `PAPER.md`, on the two ratios of a magic square of squares, known? (The congruences in Theorems 2 and 3 follow from Rabern's Theorem 1.4 and Woll's Corollary 4.3 once the Lemma is known; `PAPER.md` says so.)

## What already has an automatic check

`check_evidence.py` re-tests every recorded `t`, checks that every recorded generator lies on its curve and that the `E0L` models are the right curves, and compares the table with the JSON output of the scripts. `height_sides.py`, `height_precision.py` (the same heights at 200 bits) and `height_certify.py` are described in questions 3 and 5. `eclib_check.py` and `magma_check.py` are the second and third programs for the rank bounds. None of these proves a height, a rank or a saturation.

## How to answer

Please open an issue at https://github.com/drew2apps/magic-square-of-squares-lean/issues. "This step is wrong" and "this step is fine" are both welcome, with or without a reference.
