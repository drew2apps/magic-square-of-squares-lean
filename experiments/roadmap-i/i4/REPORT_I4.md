# REPORT_I4 (interim 2, 28 Nov 2026 about 00:10 UTC): G3 for 14000 < m ≤ 20000

> **Read `../NOTES.md` first.** This is the session's report, committed unchanged apart from this line; where they differ, `../NOTES.md` (after review z17) holds.


Written by an AI model; not read by a person. Experiment, not counted. Nothing here changes the census. **Interim:** chunks 14000 to 15000 and 16000 to 17000 are done. Chunk 18000 is part way through its second descent (111 of its ratios) and is paused to give I8 the cores; 15000, 17000 and 19000 follow.

## Bottom line so far

- **Chunk 14000 to 15000: 301 ratios, 2 open.** 45 are left after the second descent, 11 after rank 0 and 2 after Demjanenko–Manin. The two open rows are (14017,2856) and (14640,10679). No DM alert. *(computed; the rank-0 stage is "consistent with PARI" where it uses ellrank, as in G3)*
- **Chunk 16000 to 17000: 305 ratios, 2 open.** 47 are left after the second descent, 22 after rank 0 and 2 after Demjanenko–Manin: (16065,2032) and (16873,11161). No DM alert. *(computed, as above)*
- **(14017,2856) and its partner (16873,11161) are both open,** so that pair stays open. On both, every curve that DM tries has rank bounds other than exactly 1 (for (16873,11161): E [2,2], D1 [1,3], Sig+ [2,4], Sig− [2,2]). The partner of (16065,2032) is (18097,14033), in chunk 18000.
- **Partners of the two open rows of chunk 14000:** (14017,2856) ↔ (16873,11161) and (14640,10679) ↔ (25319,3961), with the partner map (m + n, m − n), halved when both are odd. (16873,11161) is in chunk 16000 and is open too (above). (25319,3961) is outside I4's range.
- **Their curves:** for (14640,10679), several curves have rank bounds [0, 2] (f0·K, g4·K, f8·K, w·K among them). For (14017,2856), w·f0·f8 and w·K have [0, 2] and g4·K has [1, 3]. Both need H7's step or the g4·K route next.
- **The exchange test finds no failure.** Every point the point searches reported in chunks 14000 and 16000, and in chunk 18000's second descent so far, passes: 828 'no' and 1,099 trivial points (t ∈ {0, ±1, ∞}), 0 failures. The same test on G3's own outputs (roadmap-g/g3) gives 767 'no', 1,098 trivial and 0 failures. Every reported point so far is trivial, so the ratio check can only bite on a SOLUTION, and none has come up. *(computed, exact)*

## The exchange test (from the addendum)

`i4_exchange.py` reads the second descent's `rational_t` and Demjanenko–Manin's `t` from every sd_*.jsonl and dm_*.jsonl, partners included. For each point t = v/u it:

1. recomputes the verdict at (m, n, t) in exact rationals;
2. recomputes it at the exchanged point (u, v, m, n), ratio u : v with t′ = n/m, and requires the same verdict;
3. sorts the point: trivial when uv(u² − v²) = 0; otherwise u : v reduced to u > v > 0 coprime, and a SOLUTION's u : v must be a valid ratio (not degenerate, outside modp_classes.txt, keeping a class by class_selmer's criterion), or it prints EXCHANGE FAILS.

Negative controls (`--negative`): the same t = 1 point at (14744,3633) claimed as a SOLUTION is refused; t = 2/3 claimed degenerate is refused; u : v = 7 : 1 (degenerate) and 3 : 2 (in a mod-p class) are refused; 308 : 185 (I1, design B) is accepted. `i4_exchange.txt` is a snapshot taken with chunk 18000 part way through, and will be rewritten at the end.

## What failed

- Nothing in the results. One bookkeeping slip: a restart script of mine wrote duplicate records into sd_14000_1.jsonl while the worker was still alive. I removed the script, and the worker's own writes left 150 unique records, the same 301 ratios in all as the input.

## Files and checks

- run_i4.sh, pipe_g3.py, gen_family.py, fam_11000_30000.json, chunk_*.json (inputs), sd_*, r0_*, dm_* (stage outputs), stage_14000.json, stage_16000.json, i4_exchange.py with i4_exchange.txt and i4_exchange_negative.txt, times.txt.
- It runs from the repository's experiments/roadmap-i/i4, so G3's relative paths (../../census-wide, ../../census2000) resolve.
- PARI 2.19.1 here; CI has 2.17.4. PARI finds points and rank bounds; every reported point is reclassified exactly.

```
python3 experiments/roadmap-i/i4/i4_exchange.py --negative   # about 3 s
python3 experiments/roadmap-i/i4/i4_exchange.py              # about 2 s, prints i4_exchange.txt
```

Times: chunk 14000 to 15000, 13,055 s on 2 low-priority cores, and its later stages 615 s. Chunk 16000 to 17000 took 10,093 s, and its later stages 870 s.

## Proposed LEDGER/CHANGES text (draft)

> **Experiment, not counted.** I4 (interim): G3's pipeline leaves 2 of 301 ratios open on 14000 < m ≤ 15000, (14017,2856) and (14640,10679), and 2 of 305 on 16000 < m ≤ 17000, (16065,2032) and (16873,11161), the partner of (14017,2856). Every point the point searches report passes the exchange test (a point (m, n, u, v) gives the ratio u : v). Computed. Written by an AI model; not read by a person.
