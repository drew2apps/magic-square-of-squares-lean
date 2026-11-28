# Roadmap I: hand-off, part 2 (28 Nov 2026, about 00:10 UTC)

> **Read `NOTES.md` first.** This is the session's report, committed unchanged apart from this line; where they differ, `NOTES.md` (after review z17) holds.


Written by an AI model (Claude); not read by a person. Experiments, not counted. Nothing here changes the census. Nothing was committed, pushed, posted or sent; Drew publishes. Θ₄ was not started (it is the conductor's), and neither were I6 or I7.

## "Neither" first (I8, task 8 of the addendum)

On lines 1 to 360 of `i2/i2_needgood_list.txt`, 6 residual tuples in 3 ratios are removed neither by the Cassels–Tate pairing on any of the six pair curves nor by a pair curve of rank 0. E± was tried on each.

| Ratio | Tuples (good prime that kills them) | Pair curves | E± |
|---|---|---|---|
| (18313,120) | [1,2,1,335351569], [1,2,670703138,2] (61) | kernel on 08, 0L, 0U, 8L, 8U; on LU, the class of a point found | E+ rank 1 by PARI; Demjanenko–Manin (z ↦ −z): 119 candidates, no four-corner point |
| (18433,18193), the partner | [2,670703138,2,1], [335351569,1,2,1] (61) | 08: class of a point found; kernel on the other five | no curve of rank 1 (E+ [1,3], E− [2,2]); closed through its partner's E+ |
| (20056,19825) | [2,1,1,1], [2,1,2,2] (17) | kernel on all six | E− rank 1 by PARI; DM: 113 candidates, no point. Its partner (39881,231), not in the list, closes the same way on its curve E |

All three closures are **conditional**:
- the rank 1 is PARI's `ellrank` alone, and the repository's exact descent (`rank_indep.rank_bound`) gives only rank ≤ 3 on both curves;
- the generator is saturated only at primes ≤ 1000;
- the height bounds are in floating point.

The same ranks come out under PARI 2.17.2. So the pairing does not reproduce every good-prime kill: these 6 tuples fall to a good prime and to a rank-1 E± curve, but to no pairing. *(computed; conditional as stated)*

## For the conductor

1. **I8 checkpoint:**
   - Coverage is lines 1 to 360 of the list, all from H2's exhaustive list, in `i8/i8_first360.jsonl`.
   - 3,160 of 3,166 tuples are pairing-removed, the 6 above are "neither", and 11 ratios are also excluded by a pair curve of rank 0.
   - 462 tuples are killed only on the mixed curves 0L, 0U, 8L, 8U.
   - The run continues in list order on all 3,876 ratios; G2's 937 are lines 2,940 to 3,876.
   - Each tuple's killing curves are recorded (`killed_by`).
2. **I2 is complete:**
   - `i2/i2_residual.jsonl` covers all 3,876 ratios, with the same fields in the same order as the interim file: 32,474 tuples, every |S| even, every tuple killed by a good prime ≤ 229.
   - The shape run did 3,838 (38 over its memory cap), and REPORT_I2 has the full counts: 70% cosets, killing sets of 3 to 5 primes for 2,557, and the separating quotient at P never 1 coordinate.
3. **PARI versions:**
   - The diffed outputs (`i8_negative.txt`, `i8_compare50.txt`, `i2_residual_negative.txt`) print only facts that don't depend on PARI's choice of Sel² basis or points.
   - They print exactly the committed text under PARI 2.17.2 (cypari2 2.2.3 in a venv) as well as 2.19.1.
   - 2.17.4 was not available here, so nothing was run under 2.17.4 itself.
   - In the jsonl, a per-curve label "image" versus "kernel" depends on the points PARI finds; the verdicts do not.
4. **I4:**
   - Chunk 16000 is done, with 2 of 305 open: (16065,2032) and (16873,11161).
   - (16873,11161) is the partner of chunk 14000's open (14017,2856), so that pair stays open.
   - The exchange test finds 0 failures over 1,927 reported points.
   - Chunk 18000 is paused part way so that I8 gets the cores.
5. **I3:**
   - The 50,014 new ratios started with the two largest ranges (k = 29, 28), and neither has finished. No new result yet.
   - One chain is running again; the other is paused for I8.

## Per task

- **I2 (done):** i2/REPORT_I2.md (updated), i2_residual.jsonl, the shape outputs `i2_i2_needgood_list_0of2.txt` and `_1of2.txt`.
- **I8 (checkpoint at 360):** i8/REPORT_I8.md, i8_pairct.py, i8_checkpoint.py, i8_first360.jsonl and its summary, the controls, and neither/ (E± scripts and outputs).
- **I4 (interim 2):** i4/REPORT_I4.md (updated), stage_16000.json, i4_exchange.txt (snapshot).
- **I3:** unchanged since part 1 (REPORT_I3.md).

## Times and proposed CI lines

| Command | Time |
|---|---|
| `python3 experiments/roadmap-i/i2/i2_residual.py --negative` | 3 s |
| `python3 experiments/roadmap-i/i8/i8_pairct.py --negative` | 13 s |
| `python3 experiments/roadmap-i/i8/i8_pairct.py --compare50` | 5 s |
| `python3 experiments/roadmap-i/i4/i4_exchange.py --negative` | 3 s |

The run times are in each folder's times.txt. I8 needs `experiments/conductor-t10` (main at f0d7edc) and `roadmap-i/i2` beside it. PARI 2.19.1 with cypari2 produced every output; real sympy 1.14.0.
