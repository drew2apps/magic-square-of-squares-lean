# Ratios with 3000 < m ≤ 5000 (experiment, 26 Nov 2026)

An AI model (Claude) wrote this, and no person has read it. It changes nothing in the repository, whose census stops at m ≤ 1000; any change there is the owner's call. Nothing has been committed or pushed. It continues `../census2000/CENSUS3000.md`.

## Result
**Every coprime ratio m : n with 3000 < m ≤ 5000 is excluded, subject to the conditions of the label table below: PARI's rank bounds throughout, Demjanenko–Manin generators saturated only at primes ≤ 1000, and for the rank-1 route (a) and (b), and two more only conditionally on a rank being 1.** Those two are (3540,2581) and (4440,2231).
- **What the two depend on:** E1(ℚ) has rank exactly 1. PARI only bounds that rank between 1 and 3.
- **What every exclusion depends on:**
  - the programs below being right;
  - the repository's tier C/D steps.
- **Tier:** all of this is tier D, computed rather than proved. Nothing above m = 1000 is in Lean.

## How the 2918 ratios split
| Step | Ratios left | Label |
|---|---|---|
| Congruences (`gen_pairs.py`, Result 1b): **2918** pass. Of these, 2305 are of other type, 356 of square type and 257 of twice-square type. | 2918 | computed (exact) |
| **Other type (2305):** the repository's class criterion with four-class descent, run to m ≤ 15000 | 613 | repository, tier D |
| **Second descent on the 613 degenerate-type ratios** (`sd_pari.py`, my PARI port of `second_descent.py`): it closes 433 directly and 92 through their partner. | 88 | computed; PARI rank bounds |
| **Rank-0 quotients** (the repository's `global.py`, unchanged), on each ratio or its partner: closes 36 | 52 | computed |
| **Demjanenko–Manin** (`../census2000/dm_pari.py`), on each ratio or its partner: closes 47, with no alerts | 5 | computed; conditional on PARI's generators (saturation ≤ 1000) |
| **The repository's own route closes the partners of two of the five.** (3833,977)'s partner is (2405,1428) (`route_twist.txt`). (3631,1031)'s partner is (2331,1300) (`route_wide.txt`). | 3 | repository, tier C/D |
| **The rank-1 route ported to PARI** (`route/route_pari.py`, below) closes (3321,2440) | 2 | computed; see the labels |
| **(3540,2581) and (4440,2231):** the same route closes them **if E1(ℚ) has rank 1** | 0 (conditional) | conditional |

Files:
- second-descent sweep: `sw_*.jsonl`;
- rank 0: `g5k*.txt`;
- Demjanenko–Manin: `dm5k*.jsonl`;
- route: `route/out_*.json`.

## The route port (`route/route_pari.py`)
The repository closes square-type ratios with E1(ℚ) of rank 1 by a Mordell–Weil sieve plus a height bound (Checcoli–Veneziano–Viada). It does this with four Sage scripts: `saturation.py`, `mw_sieve.py`, `closure_estimate.py` and `fibre.py`. Sage isn't installed here, so I ported them to PARI and plain Python.
- **What matches the repository:** the local sets, the lifting, the exact-rational height bound, Cremona–Prickett–Siksek's enclosure of ĥ(G) (its helpers imported from the repository), and the fibre equations over ℚ(i).
- **Validation:** the port was run on the repository's two route closures in this range's partners.
  - **(2331,1300), d = 7487122.** Every number matches `route_wide.txt`: the primes 10009, 20089, 30013 and 40009; the height bound a, b ≤ 14111; the sieve to 30000 from one prime, leaving [(1,1)]; and the fibre, which gives t = ±1, degenerate.
  - **(2405,1428), d = −7489682:** matches `route_twist.txt`, with bound 14101 and [(1,1)] left.
  - **The lifting** also passes the repository's planted-pairs test on every run.
- **(3321,2440):** run with d = 10150882 and again with −d.
  - Both runs: E1 has rank 1. PARI bounds E1 itself at [1,3], but the 2-isogenous curve gives [1,1].
  - Both runs: G is saturated at every prime ≤ 2000, and the sieve leaves only [(1,1)].
  - The height bounds are 14042 and 14062, both within the sieve range of 30000.
  - The fibre gives only t = ±1, which is degenerate.
- **(3540,2581) with −d, and (4440,2231) with −d:** the sieve, to 10⁶, leaves only [(1,1)], under the bounds 14116 and 13394, and the fibre gives t = ±1.
  - What fails is the rank. Every curve in the 2-isogeny class gets PARI bounds [1,3], and [1,5] on one curve for 4440, even at effort 8.
  - With d (not −d), E1 has rank [2,2], so the route can't be used.
  - The partners, (6121,959) and (6671,2209), are of twice-square type, where the route doesn't apply. Rank 0 and Demjanenko–Manin fail on them too.

## Labels
| Claim | Label |
|---|---|
| 2918 ratios pass the congruences; 613 are of degenerate type | computed (exact) |
| The second-descent port agrees with the repository on 48 of 48 test ratios | computed (exact) |
| Second descent, rank 0 and Demjanenko–Manin leave 5 | computed; PARI rank upper bounds; Demjanenko–Manin conditional on generators |
| (3321,2440) is closed by the route | computed, conditional on (a) G generating E1(ℚ) modulo torsion, (b) the reading of Checcoli–Veneziano–Viada's constants |
| (3540,2581) and (4440,2231) are closed | **conditional on rank E1(ℚ) = 1**, plus (a) and (b) |

About (a): PARI checks saturation only at primes ≤ 2000. The repository bounds the index with Sage's Cremona–Siksek lower bound for heights, which PARI lacks. Without it, the index could still be a product of primes above 2000. That would mean a generator of height below about 25/2003² ≈ 6·10⁻⁶.

## What would firm it up
1. **Rank 1 for E1 at (3540,2581) and (4440,2231), with −d.** This needs a 4-descent or Cassels–Tate pairing, from Magma or Sage's `mwrank`, which is not available here. These two then join (627,364) as the cases where only a better rank bound is missing. Here, though, the curves are over ℚ, not ℚ(i), so Magma's `FourDescent` should be routine.
2. **A height lower bound** (Sage `height_function().min` or Magma) to turn saturation ≤ 2000 into a full proof for the three route closures.
3. **The usual items from CENSUS3000.md:** proven generators for Demjanenko–Manin and a second implementation of the rank bounds.
