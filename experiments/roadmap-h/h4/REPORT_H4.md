# H4: the three ratios of G1 that ran out of memory

> **Read `NOTES.md` first.** This is the session's report, committed unchanged apart from this line; where they differ, `NOTES.md` (after review) holds.

Written by an AI model (Claude); not read by a person. Experiment, not counted. Nothing here changes the census, which stops at m ≤ 1000.

Machine: PARI 2.19.1 with cypari2, not CI's PARI 2.17.4. Real sympy 1.14.0. The programs are G1's exact_local.py and G2's g2_killsets.py from the roadmap G bundle, unchanged.

## Bottom line

- **All three are emptied by their bad primes alone.** *(computed, exact: every bad prime tested exactly; at 2 the repository's test at depth 12)*

| ratio | log₂ of the tuple space | bad primes alone | irredundant kill set |
|---|---|---|---|
| (88373348,70002375) | 32 | empty | 2, 13, 1777 |
| (120756965,14608109) | 34 | empty | 2, 5, 17 |
| (711343101,552021980) | 36 | empty | 2, 3, 13, 17, 53 |

- **G1's 63 kept ratios are now all done.** With G2's 60, bad primes alone empty 60 of the 63 and 3 need good primes, the largest 61.

## How

- **What it computes:** G2's computation exactly (g2_killsets.run):
  - every bad prime, exact;
  - then good primes in increasing order;
  - then an irredundant kill set.
- **One change, to save memory:** tuple_descent's join is run on the tuple space split by the class d0, one chunk for each divisor d0 of the first form's support, and the chunks are joined one at a time.
  - The union of the chunks is the join of the whole space, so the sets are the same.
  - An emptiness test stops at the first chunk that keeps a tuple.
- **Memory:** the run fits in 8 GB (`ulimit -v 8000000`).

## Checks

- **Negative mode (h4_negative.txt):**
  - **Split against unsplit:** on four ratios small enough for G2's unsplit join, the split gives the same line.
    - (291,224) and (308,185) are emptied by bad primes.
    - (465,416) and (511,260) need the good primes 41 and 17.
  - **Degenerate control:** on (12,5) with the primes 2, 3, 5, 13, the split join keeps 8 tuples, among them the 2 known tuples of its degenerate points, as it must.
- **Known points:** G2's run asserts at every prime that each known point's tuple passes.

## What failed

Nothing in this part. The three ratios failed in G2 only because the unsplit join held the whole space in memory.

## Proposed LEDGER / CHANGES text

> **Experiment, not counted.** The three kept ratios of G1 that ran out of memory in G2 (experiments/roadmap-h/h4): with the join split by d0, the bad primes alone empty the four-class set of (88373348,70002375), (120756965,14608109) and (711343101,552021980). So all 63 of G1's kept ratios are done: 60 by bad primes alone, 3 with good primes up to 61. Computed (exact). Written by an AI model; not read by a person.

## Proposed CI lines

```
cd experiments/roadmap-h/h4        # F = grep -v '^PARI stack size'; needs ../../roadmap-g/g1 and g2 (or G1_DIR, G2_DIR)
python3 h4_split.py --negative 2>&1 | $F | tail -1 | grep -q '^# matches'    # about 6 s
# not for CI: python3 h4_split.py   (411 s under ulimit -v 8000000; times.txt)
```
