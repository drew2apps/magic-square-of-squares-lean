# REPORT_I1: a designed survivor hunt, and the constraints that bite

> **Read `../NOTES.md` first.** This is the session's report, committed unchanged apart from this line; where they differ, `../NOTES.md` (after review z17) holds.


Written by an AI model; not read by a person. Experiment, not counted. Nothing here changes the census.

## Bottom line

- **No survivor.** Three designed hunts took 1,569 ratios. 50 keep a class (29 distinct ratios). Bad primes alone empty the four-class set for 49 of the 50 (28 of the 29 distinct ratios). The last one, (651432,64265), is emptied at the good prime 29. None reached the check of good primes between 5000 and 67,583. *(computed, exact)*
- **What bites, with counts.**
  1. **The class criterion.** Designs that make the tuple space big lose the class. Forcing 17, 41 and 73 into m² + n² (design A) keeps a class for 4 of 900 ratios (0.44%). A smaller m² + n², with every prime ≡ 1 mod 8 (design B), keeps it for 36 of 519 (6.9%). H2 had 4.4% (excess heuristic) and 1.3% (biased random).
  2. **Bad primes kill the rest.** Of the 19 distinct ratios that keep a class, 18 die at bad primes. The irredundant killing sets are short (KILL_SUMMARY). 2 is in almost every set.
  3. **Good primes rule out constant label vectors.** From H3's table, at every good residue the vectors 0000 and 1111 are never allowed at 41 and 73. At 17 only the eight odd-weight vectors are ever allowed, and 8 of 17's 12 good residues allow nothing (i1_laws.txt, part a). So a tuple whose labels are constant at 17, 41 and 73 needs all three bad. Making them bad is what design A did, and it cost the class.
- **The laws at large bad primes are rigid** (part b, exact local test, 24 ratios per row):

  | q divides | q ≡ 1 mod 8 | q ≡ 5 mod 8 | q ≡ 3, 7 mod 8 |
  |---|---|---|---|
  | m or n | L = U (8 unit vectors), 8 with odd valuation | L = U, none with odd valuation | 0000 only |
  | m − n or m + n | bit 0 = bit 8, 8 with odd valuation | bit 0 = bit 8, none odd | 0000 only |
  | m² + n² | {0000, 0011, 1100, 1111} | {0101, 0110, 1001, 1010} | (no such q) |
  | Q+ | {0000, 1111} | {0110, 1001} | (no such q) |
  | Q− | {0000, 1111} | {0101, 1010} | (no such q) |

  Bits are in class order (d0, d8, dL, dU). Each prime of m² + n² or of the octic halves the unit vectors, and fixes how they pair. With these laws, at good 17, 41 and 73, a survivor needs labels that match the octic's law at every prime of Q±, the class law at every prime of S, and H3's tables at every good prime. Constant labels are ruled out at the good primes, and alternating labels are ruled out at the bad primes ≡ 1 mod 8. The hunt shows the two together almost always leave nothing.
- **Correction to REPORT_H2 (for the HANDOFF):** the complete good-prime bound is 67,583 (genus-129 cover), not about 66,562. REPORT_H2's line "By H3, a cutoff above 1153 is never needed" is wrong: H3's method-B tables are right for p ≤ 1153, but a survivor must be checked to 67,583. The hunt does that (extend in i1_hunt.py), but no ratio got that far.

## The hunts

Each round draws 1500 ratios m + i·n = (unit) × (product of Gaussian primes), skips degenerate ratios and ratios that fail modp_classes.txt, and keeps the one with the highest excess (h2_design: log2 of the tuple space plus H3's score at good primes ≤ 1153). Then h2_hunt's pipeline: the class criterion, then the four-class descent with every bad prime and good primes to 5000 (P = 0: bad primes alone), then every good p ≡ 1 mod 4 in (5000, 67583].

| Design | Ratios | No local point | No class | Class kept (distinct) | Emptied by bad primes alone | Survivors |
|---|---|---|---|---|---|---|
| A: 17·41·73 in m² + n², plus 1 to 6 Gaussian primes below 400 | 900 (3 seeds) | 486 | 409 | 4 (3) | 4 | 0 |
| B: m² + n² = 1 to 3 primes ≡ 1 mod 8 below 10⁶ | 519 (3 seeds, 173 each) | 0 | 463 | 36 (16) | 35 | 0 |
| C: 17, 41, 73 in m, n or m ± n, residues at 29, 37, 53, 61 that allow a target pattern (CRT backwards from H3's tables) | 150 (1 seed) | 103 | 36 | 10 (10) | 10 | 0 |

- Design B by K (primes in m² + n²): K = 1: 178 taken, 22 keep a class, but only 2 distinct ratios, (308,185) and (340,53), small and both emptied by bad primes. K = 2: 183 taken, 11 keep. K = 3: 158 taken, 3 keep. More primes in m² + n² means fewer classes, as the class criterion predicts.
- 6 rounds in B passed nothing through the congruences (lines "- K none").
- Design C (`I1_MODE=C`) works backwards from H3's tables, as the updated prompt asks: it picks a target label pattern, puts each of 17, 41 and 73 into m, n, m − n or m + n (so that each is bad, which a constant pattern needs), chooses residues of m/n at 29, 37, 53 and 61 whose H3 mask allows the pattern, and solves by CRT; of 1500 candidates per round it keeps the one with the highest excess. It keeps a class as often as design B (10 of 150, 6.7%), but 103 of the 150 have no local point at all: with 17, 41 and 73 dividing mn(m² − n²) the class criterion's local conditions fail far more often. All 10 kept ratios are emptied by their bad primes alone. One ratio, (1720428819828,5633), hit the class criterion's time limit and is not counted as kept.
- Design B was stopped at 173 per seed and replayed from its .part files. The ratios are the same as the first 173 of the stopped run.

## Killing sets of the class-keeping ratios (I2's script on i1_classkept_list.txt)

25 of the 29 distinct class-keeping ratios are done; 6 of them are over the 5 GB memory cap per process (the four-class set of those was emptied in the hunt itself, with H4's split join, so they are not open), and 4 are still running ((23746525140,545237), (45364602280,56287), (559041464,23486071), (69243112,5127009)). For the 19 done, P = 0 means the bad primes alone empty the set; the killing set is irredundant (each prime in it is needed).

| Ratio | P | Irredundant killing set |
|---|---|---|
| (1097705,550604) | 0 | 2, 5, 13, 73 |
| (156887577,143562092) | 0 | 2, 5, 13, 17, 43, 73 |
| (172983,148520) | 0 | 2, 13, 1597 |
| (183124,539) | 0 | 2, 3, 5, 17 |
| (225924,69649) | 0 | 3, 5, 13, 73 |
| (234913,61148) | 0 | 2, 29 |
| (238989,117824) | 0 | 5, 11, 13, 29, 41 |
| (2674377,1885892) | 0 | 2, 5, 41 |
| (308,185) | 0 | 2, 5, 29 |
| (340,53) | 0 | 2, 13, 17 |
| (358148,250087) | 0 | 2, 5, 13, 37 |
| (392932,125777) | 0 | 2, 5, 17 |
| (449909,413264) | 0 | 2, 5, 181 |
| (471597,172960) | 0 | 2, 13 |
| (621437309,63174836) | 0 | 2, 5, 17, 37, 277 |
| (651432,64265) | 29 | 2, 3, 29 |
| (84697,78540) | 0 | 2, 3, 11, 13, 17, 37, 47, 131 |
| (1247577551500,7023) | 0 | 2, 3, 5, 7, 11, 13, 17, 41, 73, 937, 1321, 2341, 54367 |
| (1306936960,7596161) | 0 | 2, 3, 5, 11, 13, 17, 541 |

Over the memory cap: (335143916,279093591), (726100,480963), (1037411335213,9085), (1419446173223,463), (1576718793964,3471), (233570277,677075). 2 is in 17 of the 19 sets, and 5 or 13 in 17. Primes in these sets with P = 0 are all bad primes of their ratio (29, 37, 41 and 73 among them). At every prime of every set, the smallest label quotient that separates uses 2 coordinates or more, except at 3 and 5 for (1247577551500,7023), where 1 coordinate does.

## What holds and how it was checked

- The class criterion and the four-class descent are the repository's (class_selmer, tuple_descent, exact_local) through h2_hunt.py, unchanged. Every bad prime is tested exactly (`assert P.isprime`). Every known tuple of a ratio is checked to be allowed at every prime tested.
- i1_negative.txt (`--negative` checks): a forged survivor claim, (308,185) from design B, is refused on recomputing (P = 0). A forged local table for the degenerate ratio (7,1), with its known tuple removed at 5, is refused by the pipeline's own check. Positive control: H1's phantom (204724,97893) is not emptied by bad primes and good primes to 5000, as its phantom tuples require.
- i1_laws_negative.txt: a forged claim that 0000 is allowed at good 17 is refused, and a forged unit vector 0110 at a prime ≡ 1 mod 8 of Q+ is refused. 0001 at 17 and 1111 at that prime pass.
- Rests on: H3's h3_weil.json (method B, exact, good primes ≤ 1153), the repository's local test, and the genus-129 bound for the complete check.
- PARI: 2.19.1 here (CI has 2.17.4). PARI is used only for factoring and prime lists, and the outputs don't depend on its order.

## What failed

- No ratio got past the bad primes and the good primes to 5000, so the 67,583 check never ran on a hunt ratio. The extend step is exercised only by H2's ratios, not here.
- A design that keeps both the class and a big tuple space was not found. The two pull against each other: each prime of m² + n² or of the octic adds tuples but halves the chance of keeping a class, and its law then cuts the tuples back.

## Files

- i1_hunt.py, run_i1_hunt.sh; i1_hunt_s{1,2,3}.txt (design A), i1_hunt_B_s{1,2,3}.txt (design B, `I1_MODE=B`), i1_hunt_C_s1.txt (design C, `I1_MODE=C`).
- i1_laws.py, i1_laws.txt, i1_laws_negative.txt.
- i1_negative.py, i1_negative.txt.
- The killing sets of the 29 class-keeping ratios are in ../i2/i2_i1_classkept_list_0of1.txt (I2's script on the list i1_classkept_list.txt).
- times.txt. i1_explore.py is scratch, not for the repository.

## Proposed LEDGER/CHANGES text

> **Experiment, not counted.** I1: three designed survivor hunts (1,569 ratios, good primes checked to the complete bound 67,583 when needed) found no survivor. 50 ratios keep a class, and bad primes alone empty the four-class set for 49 of them (the other at the good prime 29). H3's tables show 0000 and 1111 are never allowed at good 41 and 73, and only odd-weight vectors at good 17, and the laws at large bad primes fix the unit labels by type and class mod 8 (i1_laws.txt). Computed (exact). Written by an AI model; not read by a person.

## Proposed CI lines

```
python3 experiments/roadmap-i/i1/i1_laws.py              # 1 s
python3 experiments/roadmap-i/i1/i1_laws.py --negative   # 1 s
python3 experiments/roadmap-i/i1/i1_negative.py          # 2 s
```

The hunts take 270 to 420 s per seed (times.txt), so they are not CI lines. Replaying one seed checks against its file with `# matches`, at the same cost.
