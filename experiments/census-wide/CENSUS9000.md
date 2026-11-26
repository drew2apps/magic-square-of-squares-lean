# Ratios with 7000 < m ≤ 9000 (experiment, 26 Nov 2026)

An AI model (Claude) wrote this, and no person has read it. It changes nothing in the repository, whose census stops at m ≤ 1000; any change there is the owner's call. Nothing has been committed or pushed. It continues `CENSUS7000.md` with the same programs, unchanged.

## Result
**Every coprime ratio m : n with 7000 < m ≤ 9000 is excluded except one, (7735,4968), subject to the conditions of the label table below: PARI's rank bounds throughout, Demjanenko–Manin generators saturated only at primes ≤ 1000, and for the rank-1 route (a) and (b), and one more, (8249,3720), only conditionally on a rank being 1.**
- **(8249,3720)** is closed by the rank-1 route **if** its curve E1 has rank exactly 1. PARI bounds it in [1,3] on both curves of its 2-isogeny class, also at effort 8, and the root number is −1. This is the same situation as (6480,4961) in `CENSUS7000.md`.
- **(7735,4968) is open.** Every curve the pipeline can use has rank ≥ 2, or a 51-digit conductor that rules out an L-value (details below).
- **What every exclusion depends on:** the programs below being right, and the repository's tier C/D steps.
- **Tier:** all of this is tier D, computed rather than proved. Nothing above m = 1000 is in Lean.

## How the 5366 ratios split
| Step | Ratios left | Label |
|---|---|---|
| Congruences (`../census2000/gen_pairs.py 7000 9000`): **5366** pass. Of these, 4757 are of other type, 357 of square type and 252 of twice-square type. | 5366 | computed (exact) |
| **Other type (4757):** the repository's class criterion with four-class descent, run to m ≤ 15000 | 609 | repository, tier D |
| **Second descent** on the 609 (`wide_pari.py`): closes 447 directly (by '1356' 295, '+14' 88, '+23' 64) and 73 through their partner | 89 | computed; PARI rank bounds |
| **Rank-0 quotients** (the repository's `global.py`, unchanged) on the 89 and their partners: closes 49 | 40 | computed |
| **Demjanenko–Manin** (`../census2000/dm_pari.py`) on the 40 and their partners: closes 33, no alerts | 7 | computed; generators saturated at primes ≤ 1000 |
| **Route, partner closed in `CENSUS7000.md`:** (7393,5783), partner (6588,805); (7409,2737), partner (5073,2336) | 5 | computed, conditional as there |
| **Route (`route/route_pari.py`), new:** (7811,3060) with d = 103296242, (8364,4123) with d = −105914734, (8687,4416) with d = 111925826. Each: E1 rank 1 over the whole 2-isogeny class, G saturated at primes ≤ 2000, sieve leaves [(1,1)] under the bound (13993, 14174, 14006), fibre gives only t = ±1 | 2 | computed; conditions (a), (b) of `CENSUS5000.md` |
| **Conditional on a rank being 1:** (8249,3720), d = −108415202. Same route: sieve to 10⁶ leaves [(1,1)] under the bound 14079, fibre gives t = ±1 | 1 | conditional |
| **(7735,4968)** | 1 | open |

Files:
- listing: `pairs_7000_9000.json`, `degen_7000_9000.json`, `gen9k.log`;
- second descent: `s9_*.jsonl`, `s9r_*.jsonl` (the second batch finished the ratios the first run's time limit cut off);
- rank 0 and Demjanenko–Manin: `pipe9k.py`, `g9k*.txt`, `after_rank0_9k.json`, `dm9k*.jsonl`, `after_dm_9k.json`;
- route: `route/scope_9k.jsonl`, `route/log_7811.txt`, `route/log_8364.txt`, `route/log_8687.txt`, `route/log_8249.txt`;
- the hard cases: `hard9k/`.

The route was also run on the twice-square partners (11969,4529) and (12703,2767) (`route/log_11969.txt`, `route/log_12703.txt`). Those runs are not counted: the route is built for square type, and its fibre step reports "no" rather than "degenerate" there.

## (8249,3720)
- Square type; with −d, E1 has PARI bounds [1,3] on both curves of its 2-isogeny class at effort 8 (`hard9k/rh_8249_e8.txt`); with +d the bounds are [0,4]. Root number −1 (computed, exact), so the rank is 1 or 3.
- Assuming rank 1, the route closes it (`route/log_8249.txt`).
- Its partner (11969,4529) is of twice-square type. Its rank-0 candidates (f0·g4·f8 and w·f0·g4·f8, bounds [0,2], root number +1) have 27-digit conductors, far beyond the series method used for (6557,324).

## (7735,4968)
- **Route:** E1 has rank [2,2] with both d and −d.
- **Demjanenko–Manin:** every curve (E, D1, Σ±, for the ratio and for its partner (12703,2767)) has rank ≥ 2. One run gave D1 of the partner [0,2], another [2,2]; PARI's lower bound 2 may come from parity (one point plus the root number), so it shows rank ≥ 2 only if Ш is finite; no exclusion depends on it.
- **Rank-0 quotients:** the partner's w·K looked promising, with bounds [0,2], root number +1 and an 18-digit conductor, small enough to sum its L-series. But two curves in its isogeny class have PARI bounds [2,2] at effort 8 (`hard9k/rn8_12703.txt`), so its rank is 2. I stopped the L-series run. f0·K has PARI bounds [2,2] with one point (`hard9k/rn8_more.txt`), so rank 2 there needs Ш finite. The remaining [0,2] quotients have 51-digit conductors.
- **What might work:** the rank-1 quotients g4·K and w·g4·K of (7735,4968) (bounds [1,1]). `global.py` uses only rank-0 quotients; a Chabauty or Demjanenko–Manin step on those would be new.

## Labels
| Claim | Label |
|---|---|
| 5366 ratios pass the congruences; 609 are of degenerate type | computed (exact) |
| Second descent, rank 0 and Demjanenko–Manin leave 7 | computed; PARI rank upper bounds; Demjanenko–Manin conditional on generators |
| (7811,3060), (8364,4123), (8687,4416) closed by the route; (7393,5783), (7409,2737) through their partners | computed, conditional on (a) G generating E1(ℚ) modulo torsion (saturation checked at primes ≤ 2000), (b) the reading of Checcoli–Veneziano–Viada's constants |
| (8249,3720) closed | **conditional on a rank being 1** (root number −1, PARI bounds [1,3]), plus (a) and (b) |
| (7735,4968) | **open** |

## Running total, 1000 < m ≤ 9000
All but one ratio excluded: (7735,4968) is open. Conditions: seven exclusions are conditional on a rank being 1 ((3540,2581), (4440,2231), (6121,959), (6671,2209), (6480,4961), (6903,5104), (8249,3720); they rest on five curves, since partners share one); nine on the rank-1 route's (a) and (b) ((3321,2440), (5073,2336), (6588,805), (7811,3060), (8364,4123), (8687,4416), and through partners (5761,881), (7393,5783), (7409,2737)); every Demjanenko–Manin exclusion on PARI generators saturated only at primes ≤ 1000; (6557,324) and (6881,6233) on a numerical L-value; and all on PARI's rank bounds. (627,364), below 1000, is still open in the repository.
