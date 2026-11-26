# Ratios with 5000 < m ≤ 7000 (experiment, 26 Nov 2026)

An AI model (Claude) wrote this, and no person has read it. It changes nothing in the repository, whose census stops at m ≤ 1000; any change there is the owner's call. Nothing has been committed or pushed. It continues `CENSUS5000.md` with the same programs, unchanged.

## Result
**Every coprime ratio m : n with 5000 < m ≤ 7000 is excluded, subject to the conditions of the label table below: PARI's rank bounds throughout, Demjanenko–Manin generators saturated only at primes ≤ 1000, and for the rank-1 route (a) and (b), four of them only conditionally on a rank being 1, and (6557,324), (6881,6233) on a numerical L-value.** The four are (6121,959), (6671,2209), (6480,4961) and (6903,5104).
- **What the four depend on:** an elliptic curve over ℚ having rank exactly 1, where PARI only bounds the rank between 1 and 3.
- **The last hard pair, (6557,324) and its partner (6881,6233), is closed.** A quotient curve has L(E,1) ≈ 1.44409 ≠ 0, so its rank is 0 by Gross–Zagier–Kolyvagin, and the repository's rank-0 step closes both. The theorem is proved; the L-value is computed numerically (see below).
- **What every exclusion depends on:** the programs below being right, and the repository's tier C/D steps.
- **Tier:** all of this is tier D, computed rather than proved. Nothing above m = 1000 is in Lean.

## How the 4116 ratios split
| Step | Ratios left | Label |
|---|---|---|
| Congruences (`../census2000/gen_pairs.py 5000 7000`): **4116** pass. Of these, 3508 are of other type, 357 of square type and 251 of twice-square type. | 4116 | computed (exact) |
| **Other type (3508):** the repository's class criterion with four-class descent, run to m ≤ 15000 | 608 | repository, tier D |
| **Second descent** on the 608 (`wide_pari.py`): closes 439 directly (by '1356' 286, '+14' 91, '+23' 62) and 68 through their partner | 101 | computed; PARI rank bounds |
| **Rank-0 quotients** (the repository's `global.py`, unchanged) on the 101 and their partners (176 ratios): closes 54 | 47 | computed |
| **Demjanenko–Manin** (`../census2000/dm_pari.py`) on the 47 and their partners: closes 38, no alerts | 9 | computed; generators saturated at primes ≤ 1000 |
| **Route, partner closed in `CENSUS5000.md`:** (5761,881), partner (3321,2440) | 8 | computed, conditional as there |
| **Route (`route/route_pari.py`), new:** (5073,2336) with d = 40556866, and (6588,805) with d = −85507438. Each: rank 1 over the whole 2-isogeny class, G saturated at primes ≤ 2000, sieve leaves [(1,1)] under the bound (14063 and 14222), fibre gives only t = ±1 | 6 | computed; same conditions (a), (b) as in `CENSUS5000.md` |
| **Conditional on a rank being 1** (see below): (6121,959), (6671,2209), (6480,4961), (6903,5104) | 2 | conditional |
| **(6557,324) and its partner (6881,6233)** | 0 | rank 0 from a numerical L-value; see below |

Files:
- listing: `pairs_5000_7000.json`, `degen_5000_7000.json`, `gen7k.log`;
- second descent: `s7_*.jsonl`;
- rank 0: `g7k*.txt`, `after_rank0_7k.json`;
- Demjanenko–Manin: `dm7k*.jsonl`, `after_dm_7k.json`;
- route: `route/scope_7k.jsonl`, `route/log_5073.txt`, `route/log_6588.txt`, `route/log_6480.txt`;
- the hard cases: `hard7k/`.

## The conditional four
All four need an elliptic curve over ℚ to have rank exactly 1, where PARI only bounds it in [1,3]. In each case the root number is −1 (computed, exact), so with PARI's bounds [1,3] the rank is 1, 2 or 3, and 1 or 3 if Ш(E)[2^∞] is finite (2^∞-parity); it is 1 if moreover Ш(E)[2] is as large as the Selmer group allows.
- **(6121,959) and (6671,2209):** partners of (3540,2581) and (4440,2231), the two conditional ratios of `CENSUS5000.md`. Rank 0 and Demjanenko–Manin fail on them too.
- **(6480,4961), d = 709342:** the same picture as (3540,2581). The sieve to 10⁶ leaves [(1,1)], the bound is 13684, the fibre gives t = ±1; every curve in the 2-isogeny class has PARI bounds [1,3], also at effort 8 (`route/rh_6480_e8.txt`). With −d the rank is 2. Its partner (11441,1519) is of twice-square type, and its rank-0 quotients have ranks ≥ 1 even at effort 4 (`hard7k/gh_6480.jsonl`).
- **(6903,5104):** every route curve has rank ≥ 2. What works is Demjanenko–Manin on the D1 curve of its partner (12007,1799): PARI bounds [1,3] at effort 6, root number −1. **Assuming rank 1**, the method closes it: ĥ(G) = 18.44, height bound 174.25, multiples a ≤ 21, 43 candidate x-values, no rational t (`hard7k/dm_12007_D1.json`, made by `hard7k/dm_assume1.py`).

The conductors of these curves have 49 to 55 digits, too large for PARI to compute L′(E,1) and prove rank 1 through Gross–Zagier–Kolyvagin.

## (6557,324) and (6881,6233)
The square-type ratio (6557,324) and its partner (6881,6233) resisted every other step: every route curve and every Demjanenko–Manin curve has rank ≥ 2 or an uncertain rank. Two rank-0 candidates were left, the quotient w·f0·f8 of (6557,324) and w·K of (6881,6233), each with PARI bounds [0,2] over its whole 2-isogeny class, also at effort 8 (`hard7k/gh_6557.jsonl`, `hard7k/gh_6881.jsonl`).
- **The two curves are isogenous.** Both have conductor 11077451294972790 = 2·3·5·7·13·23·79·83·101·271·983, the same a_p, and w·K's j-invariant is in the 8-curve isogeny class of w·f0·f8 (computed, exact). So one L-value decides both.
- **Root number +1** for both (computed, exact), so the rank is 0 or 2 by parity if Ш[2^∞] is finite (not needed below: L(E,1) ≠ 0 gives rank 0 directly).
- **PARI's `ellL1` could not do it:** it ran past 2 hours on a 12 GB stack with no result (`hard7k/l1.py`, `l1_6557.txt`).
- **So I summed the series myself** (`hard7k/l1_sieve.py`): L(E,1) = 2 Σ a_n/n · e^{−2πn/√N}, with a_p from PARI's `ellap` for all primes up to 3.5·10⁸ and a_n by a multiplicative sieve. The tail factor at the cutoff is e^{−2πX/√N} ≈ 8·10⁻¹⁰; an elementary bound (|aₙ| ≤ d(n)√n, d(n) ≤ 6720 for n ≤ 10¹²) gives a tail below 0.01 at X = 3.5·10⁸ (a reviewer's bound; the sum itself is in floating point with aₚ from PARI).
  - **Validated** against `ellL1` on 11a1 and on three curves of conductor 2·10⁸, 5·10⁸ and 1.4·10¹²: 15–16 matching digits each (the session's runs; their logs are not in this folder; a reviewer reproduced 11a1 and a zero for 389a).
  - **Result: L(E,1) ≈ 1.444085956** (`hard7k/l1s_6557.txt`). The partial sums from 3·10⁸ to 3.5·10⁸ agree to 8 digits.
  - **Consistency check (not used for the rank):** PARI's `ellbsd` gives Ω·∏c_p/|T|² = 0.0902553722…, and L(E,1) divided by it is 15.9999999997, close to the square 16 that BSD would predict for the order of Ш.
- **So rank E(ℚ) = 0** (Gross–Zagier–Kolyvagin: L(E,1) ≠ 0 implies rank 0). The repository's torsion check on this curve (`hard7k/tors_check.py`) then finds only degenerate t, which closes (6557,324), and its partner (6881,6233) with it.

## Labels
| Claim | Label |
|---|---|
| 4116 ratios pass the congruences; 608 are of degenerate type | computed (exact) |
| Second descent, rank 0 and Demjanenko–Manin leave 9 | computed; PARI rank upper bounds; Demjanenko–Manin conditional on generators |
| (5073,2336), (6588,805) closed by the route; (5761,881) through its partner | computed, conditional on (a) G generating E1(ℚ) modulo torsion (saturation checked at primes ≤ 2000), (b) the reading of Checcoli–Veneziano–Viada's constants |
| (6121,959), (6671,2209), (6480,4961) closed | **conditional on a rank being 1** (root number −1, PARI bounds [1,3]), plus (a) and (b) |
| (6903,5104) closed | **conditional on a rank being 1** (root number −1, PARI bounds [1,3]) for the curve D1 of its partner (12007,1799), plus the Demjanenko–Manin generator condition (PARI generators, saturation at primes ≤ 1000) |
| (6557,324) and (6881,6233) closed | rank 0 is a theorem once L(E,1) ≠ 0; L(E,1) ≈ 1.44409 is **computed (numerical, floating point; aₚ from PARI)**, with an elementary tail bound below 0.01; the BSD ratio 15.9999999997 is a consistency check only |

## What would firm it up
1. **Rank 1** for the four curves above: Magma's `FourDescent` or Sage's `mwrank` on curves over ℚ. With (3540,2581) and (4440,2231) that makes six such curves, the same kind of task.
2. **A rigorous L(E,1) bound** for (6557,324): the same sum with interval arithmetic and an explicit tail bound, or Magma/Sage's `AnalyticRank` as a second implementation.
3. The items listed in `CENSUS5000.md`: a height lower bound for full saturation, proven Demjanenko–Manin generators, a second implementation of the rank bounds.
