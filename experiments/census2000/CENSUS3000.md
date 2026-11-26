# Ratios with 1000 < m ≤ 3000: what is left once every method is combined (experiment, 26 Nov 2026)

Written by an AI model (Claude), not read by a person. Nothing here changes the repository's census, which stops at m ≤ 1000 and changes only by the owner's decision. Nothing has been committed or pushed.

## The result
**If every program involved is right, every coprime ratio m : n with 1000 < m ≤ 3000 is excluded** (subject to the conditions of the label table: PARI's rank bounds throughout, Demjanenko–Manin generators saturated only at primes ≤ 1000, and for the rank-1 route (a) and (b)) by at least one of the repository's methods or by the two sweeps run here.

That is conditional in tier D (computed, not proved):
- above 1000 nothing is in Lean;
- the Demjanenko–Manin step here takes its generators from PARI, saturated only at primes ≤ 1000.

Until this run, the repository left **45 degenerate-type ratios open with 1000 < m ≤ 3000**: 26 below 2000 and 19 above, the "left" list of `scripts/second_descent/route_wide_summary.txt`. The repository's census methods (rank 0 and Demjanenko–Manin) had not been run above 1000. Run here, they close all 45.

## How the ratios split
1. **Congruences** (Result 1b of `G1_FULL_SYSTEM.md`: the 484 classes at 34 primes, plus the 5 mod 8 condition). `gen_pairs.py` reproduces the 425 ratios for m ≤ 1000 and gives **731 for 1000 < m ≤ 2000**: 177 of square type, 125 twice-square, 429 other.
2. **Rank-0 quotients** (the repository's `global.py`, unchanged): close **470 of the 731**, leaving 261. *(computed; PARI rank bounds, as in the census)*
3. **Demjanenko–Manin on four rank-1 curves**, for the ratio and for its partner: close **226 of the 261**, leaving 35. *(computed; conditional on the generator; see below)*
   - The curves are E (the repo's dm.py), D1 (dm2.py), and Σ⁺ and Σ⁻ (dm_sigma.py).
   - Code: `dm_pari.py`, a port of those scripts to PARI because Sage isn't installed here.
   - No alert came up: no candidate t was a non-degenerate solution.
4. **The 35 left:**
   - **16 are not of degenerate type.** The class criterion with the four-class descent (`COMPLETE.md`, section 7, test 7, run to m ≤ 15000) excludes every such ratio. *(repository, tier D)*
   - **19 are of degenerate type.** Of these, 17 are closed by the repository's second descent or its L-values (`wide_1000_2000.jsonl`, `wide_tally.txt`).
   - **The other two, (1064,423) and (1487,641), are partners.** They play the role that (627,364) and (991,263) play below 1000. The repository's Theorem A route closes them: class d = 1906334, sieve and height bound (`route_wide.txt`). *(repository, tier C/D)*
   - Conversely, every one of the 30 ratios with 1000 < m ≤ 2000 that the repository's descents leave open is closed by step 2 or step 3. Some closures overlap.
5. **2000 < m ≤ 3000:** I did not redo the congruences or the second descent here; the repository's `wide.py` runs cover that range. The **19 ratios it leaves** were run through steps 2 and 3, together with their partners.
   - Rank 0, on the ratio or its partner, closes 9 of them.
   - Demjanenko–Manin, on the ratio or its partner, closes the other 10.
   - That leaves none. Outputs: `w3kg*.txt`, `w3kdm*.jsonl`.

## Checks on the port (`dm_pari.py`)
- **Agreement with the repository** on its 157 ratios left after rank 0 (m ≤ 1000):
  - Every curve with rank bounds [1,1] closes. Where the repository records the same curve, the generator heights, the height bounds and the candidate counts agree.
    - For example, (15,8) on E: ĥ(G) = 3.1576, B = 29.27 (the repository has 29.24), a ≤ 27, 111 candidates.
    - (15,8) on Σ⁺: ĥ(G) = 1.94633, as in `dm_sigma_306.jsonl`.
  - The closures the repository has and the port lacks all need Magma upper rank bounds, which aren't available here: (760,39), (920,129), (945,248), (987,884) and (836,123).
  - `compare_val.py` prints every other difference. They are bookkeeping: the repository's logs record only the first curve that closed a ratio.
- **Height bounds:** the Cremona–Prickett–Siksek / Silverman bounds of the repository's `hbounds.py`, in floating point with a margin of 0.01%.
- **What isn't proved:** that the PARI point generates E(ℚ) modulo torsion. `ellsaturation` checks primes ≤ 1000 only. The repository's own dm.py uses mwrank's proven generators, or saturation where those fail.

## Labels
| Claim | Label |
|---|---|
| 731 ratios pass the congruences in 1000 < m ≤ 2000 | computed (exact) |
| Rank 0 closes 470 of them | computed: PARI `ellrank` upper bounds, as in the census |
| Demjanenko–Manin closes 226 more, and all 19 open ratios in 2000–3000 with their partners | computed, **conditional on the generators** (saturation at primes ≤ 1000) |
| No ratio with 1000 < m ≤ 3000 is left | conditional on all of the above and on the repository's tier C/D steps (class criterion to 15000, second descent, L-values, the Theorem A route) |

## What would firm it up
1. Proven generators for the rank-1 curves used: mwrank or Magma, or an index bound from a lower bound on ĥ.
2. Rank upper bounds from a second program: Magma, or the repository's exact `indep_rank.py`.
3. The congruences and the second descent rerun for 2000 < m ≤ 3000 by a second implementation.
