# H2: survivor hunt, m from 10⁴ to 10¹², and every ratio with 15000 < m ≤ 30000

> **Read `NOTES.md` first.** This is the session's report, committed unchanged apart from this line; where they differ, `NOTES.md` (after review) holds.

Written by an AI model (Claude); not read by a person. Experiment, not counted. Nothing here changes the census, which stops at m ≤ 1000.

Machine: PARI 2.19.1 with cypari2, not CI's PARI 2.17.4. Real sympy 1.14.0. PARI is a finder only: it factors, and every factor is proved prime (G1's g1_stress sets this).

A **survivor** is a ratio that keeps a class and whose four-class set is not emptied by its bad primes and the good primes up to 5000. **P** is the smallest cutoff such that the bad primes and the good primes ≤ P empty the set (P = 0: the bad primes alone).

## Bottom line

- **No survivor anywhere.** *(computed, exact)*
  - **Exhaustive, 15000 < m ≤ 30000:** every coprime ratio not of the two degenerate types, 205,172,879 in all (15 blocks of 1000 in m). 30,434 keep a class, and the four-class descent empties every one. 2,807 of them need a good prime; the largest P is 229.
  - **Biased hunt, m from 10⁴ to 10¹²:** 174,580 ratios sampled over 30 ranges, 2,251 keeping a class, every one emptied. The largest P is 137.
  - **Excess-guided hunt (h2_design), m from 2¹⁶ to 2²⁸:** 3,849 ratios, each the best of 200 by H3's heuristic excess, 171 keeping a class. 170 are emptied (h2_summary.txt counts these 170); the descent of the other did not finish ("What failed"). The largest P is 41.
  - Together 2,421 ratios keep a class outside the exhaustive range, above the 2,000 the roadmap asks for.
- **P grows slowly with m, and is bounded.** In the exhaustive run the largest P per block rises from 101 to 229 (G2: 181 for all m ≤ 15000). By H3, good primes above 1153 never cut anything, so P ≤ 1153 for every ratio whose set good primes can empty at all. A survivor would have to keep a tuple through all 96 good primes p ≡ 1 mod 4 up to 1153. *(H3 rests on Hasse–Weil, quoted)*
- **Larger tuple spaces need smaller P.** Ratios with a large tuple space have many bad primes, and the bad primes then do the work: in the hunt, P falls from 137 at log₂ 26 to 29 at log₂ 34 and 0 above. The excess-guided ratios, chosen for large tuple spaces and weak good primes, are emptied by their bad primes alone in 169 of the 170 whose descent finished; one needs P = 41.
- **The excess heuristic finds classes 3.4 times as often** as the biased hunt (4.4% against 1.3% of ratios taken keep a class), but its ratios are easier, not harder, to empty.

## Part 1: the biased hunt (h2_hunt.py, h2_hunt_k*.txt)

- Per dyadic range k (2ᵏ ≤ m < 2ᵏ⁺¹, cut to [10⁴, 10¹²]), random coprime m : n passing the congruences of modp_classes.txt, kept only if mn(m² − n²) has at least 6 + ⌊k/4⌋ distinct primes and the odd part of m² + n² at most 2. The bias was chosen on a pilot at m ≈ 10⁶ (it keeps a class about three times as often as plain random ratios).
- Each ratio goes through the repository's class criterion, then the four-class descent with every bad prime exact and good primes by bisection to 1500, then 5000. Time limits 600 s (class) and 1800 s (descent); no ratio hit either.
- Round 0: seed 4000 + k, up to 4000 sampled per range, k = 13 to 39 (not every range ran to the end, see "What failed"). Round 1: seed 5000 + k, up to 12,000 sampled, k = 13 to 22.

```
k round  sampled  kept  [(part) = unfinished]  largest P  P values (count)
13 0    4000    29   37  0:26 29:2 37:1
13 1   12000   116   101  0:109 29:3 37:1 41:1 61:1 101:1
14 0    4000    41   61  0:38 17:1 29:1 61:1
14 1   12000   123   97  0:111 17:1 29:4 37:2 41:2 61:2 97:1
15 0    4000    41   37  0:38 29:1 37:2
15 1   12000   143   37  0:135 17:1 29:4 37:3
16 0    4000    45   89  0:41 17:1 29:2 89:1
16 1   12000   152   97  0:140 17:2 29:5 37:2 41:1 61:1 97:1
17 0    4000    39   41  0:36 37:2 41:1
17 1   12000   142   41  0:133 17:1 29:7 41:1
18 0    4000    53   29  0:48 17:3 29:2
18 1   12000   157   97  0:141 17:1 29:7 37:2 53:2 61:1 73:1 97:2
19 0    4000    64   29  0:61 17:2 29:1
19 1   12000   160   37  0:155 17:1 29:3 37:1
20 0    4000    48   41  0:46 37:1 41:1
20 1   12000   192   137  0:178 17:1 29:5 41:2 53:2 61:1 73:1 113:1 137:1
21 0    4000    43   41  0:42 41:1
21 1   12000   184   73  0:171 29:6 37:3 41:2 61:1 73:1
22 0    4000    50   109  0:43 17:1 29:2 41:2 101:1 109:1
22 1   12000   174   97  0:161 17:1 29:4 41:3 53:3 73:1 97:1
23 0    4000    53   61  0:51 41:1 61:1
24 0    4000    73   89  0:72 89:1
25 0    1040    13   53  0:11 29:1 53:1
26 0    1364    30   61  0:27 37:2 61:1
27 0    1483    30   37  0:29 37:1
28 0    1559    30   0  0:30
29 0     527    15   0  0:15
30 0     265     2   0  0:2
33 0     248     9   0  0:9
36 0      94     0   -  
```

## Part 2: every ratio with 15000 < m ≤ 30000 (h2_exhaust.py, h2_ex_15000_*of3.txt)

- Every coprime m : n, 0 < n < m, not of the two degenerate types (m² + n² a square or twice a square): class criterion (class_selmer.empty, class_selmer_linalg), then the four-class descent as in part 1. Three workers, m split by residue mod 3.

```
block       m done  coprime  degenerate  no point  no class  class  need a good prime  largest P
15000-16000 1000   9422032    304    9334393    85950   1385    149    149
16000-17000 1000  10034106    305    9940644    91634   1523    161    181
17000-18000 1000  10636692    312   10538142    96706   1532    150    137
18000-19000 1000  11244016    306   11138900   103197   1613    131    149
19000-20000 1000  11859234    304   11748291   108857   1782    160    137
20000-21000 1000  12457456    301   12341857   113497   1801    157    101
21000-22000 1000  13073086    303   12951445   119410   1928    152    109
22000-23000 1000  13682912    309   13555633   124909   2061    177    137
23000-24000 1000  14282290    302   14149319   130529   2140    225    173
24000-25000 1000  14893950    305   14755428   136029   2188    175    181
25000-26000 1000  15505848    310   15361098   142132   2308    241    229
26000-27000 1000  16107112    302   15957634   146808   2368    212    193
27000-28000 1000  16716230    310   16560418   152964   2538    223    197
28000-29000 1000  17333698    300   17172467   158292   2639    243    193
29000-30000 1000  17928796    306   17761690   164172   2628    251    181
```

- **Control (h2_control.py, h2_control.txt):** h2_exhaust.py rerun on 14900 < m ≤ 15000, where G2 has the answer. It keeps exactly G2's 128 class-keeping ratios, with the same tuple-space sizes; P = 0 exactly where G2 says the bad primes alone empty the set; and on the 13 that need a good prime, P equals the largest good prime of G2's irredundant kill set (it may not exceed it). No survivor.
- **Negative control (h2_control_negative.txt):** a forged H2 output with one ratio dropped, one with a P changed from 0 to 29, and one with a P raised above G2's kill set are all refused.

## Part 3: the excess-guided hunt (h2_design.py, h2_design_k*.txt)

- H3 gives, for each good p ≡ 1 mod 4 up to 1153 and residue r = m/n mod p, the set A_p(r) of allowed label vectors. The **excess** of a ratio is log₂ of its tuple space plus Σ_p log₂(|A_p(r)|/16): the expected log₂ of the number of tuples left by the good primes if labels were uniform and independent. It is a heuristic, used only to choose ratios; nothing rests on it.
- Per range k = 16 to 28 (seed 6000 + k), 300 ratios, each the best of 200 random ratios by the excess.

```
  k=16 taken 300 kept  14 largest P 0
  k=17 taken 300 kept   8 largest P 0
  k=18 taken 300 kept  20 largest P 0
  k=19 taken 300 kept   9 largest P 41
  k=20 taken 300 kept   7 largest P 0
  k=21 taken 300 kept  13 largest P 0
  k=22 taken 300 kept   6 largest P 0
  k=23 taken 300 kept  15 largest P 0
  k=24 taken 300 kept  17 largest P 0
  k=25 taken 300 kept  18 largest P 0
  k=26 taken 300 kept  13 largest P 0
  k=27 taken 300 kept  20 largest P 0
  k=28 taken 249 kept  10 largest P 0
```

## Part 4: P against the size of the tuple space (h2_summary.txt)

```
largest P needed by log2 of the tuple space (hunt):
  20: ratios 85, need a good prime 9, largest P 97
  22: ratios 312, need a good prime 28, largest P 101
  24: ratios 550, need a good prime 41, largest P 113
  26: ratios 577, need a good prime 35, largest P 137
  28: ratios 419, need a good prime 24, largest P 109
  30: ratios 208, need a good prime 9, largest P 73
  32: ratios 76, need a good prime 5, largest P 97
  34: ratios 18, need a good prime 1, largest P 29
  36: ratios 5, need a good prime 0, largest P 0
  38: ratios 1, need a good prime 0, largest P 0
```

## Checks

- The descent and the class criterion are G1's and the repository's, with their own negative controls (G1's g1_negative.txt, test_exact_local.txt). H2 adds no certifying step: its outputs are counts and cutoffs.
- h2_control.py checks h2_exhaust.py against G2 on a block where G2 has the answer, with its own `--negative`.
- h2_summary.py rereads every committed output; without `--write` it prints exactly h2_summary.txt.

## What it rests on

- **Exact:** the class criterion, the local tests at every bad prime and at the good primes up to P, the primality of every factor.
- **Quoted:** Hasse–Weil, through H3, for "P ≤ 1153".
- **Not a proof of anything beyond the ratios run.** "No survivor" is a statement about the ratios listed. It is computed (exact) for each of them.

## What failed

- The exhaustive run's first start was killed after 12 s by a container restart (times.txt, exit 143); its second start finished.
- Above 2³⁰ the hunt's ratios have large tuple spaces and slow descents, so those ranges stopped at few ratios. Ranges k = 31, 32, 34, 35 and 37 to 39 were never started, and k = 25, 30, 33 and 36 were stopped and finalized from their partial outputs (times.txt), to free cores for round 1, H1 and H7. k = 26 ran last, on a core freed late, and stopped at its cap of 30 class-keeping ratios.
- One excess-guided ratio, (117702451,107508765) at k = 26, keeps a class but its descent hit the 1800 s limit. h2_retry.py reruns it with no limit but it ran 7905 s (about 2 h 12 min) without finishing and was stopped at 12:00 UTC (times.txt, exit 143). No h2_retry.txt was written. That ratio stays unfinished: no survivor claim either way, and it is not among the 170.

## Proposed LEDGER / CHANGES text

> **Experiment, not counted.** Survivor hunt for the four-class descent (experiments/roadmap-h/h2). Every coprime ratio with 15000 < m ≤ 30000 not of the degenerate types (205,172,879 ratios; 30,434 keep a class) is emptied by its bad primes and good primes ≤ 229. A biased random hunt over m from 10⁴ to 10¹² and an excess-guided hunt built on H3's tables add 2,421 class-keeping ratios, all emptied, the largest cutoff 137. No survivor. By H3, a cutoff above 1153 is never needed. Computed (exact). Written by an AI model; not read by a person.

## Proposed CI lines

```
cd experiments/roadmap-h/h2        # F = grep -v '^PARI stack size'; G2_DIR as in the script
python3 h2_summary.py 2>&1 | $F | tail -1 | grep -q '^# matches'                # about 1 s
python3 h2_control.py 2>&1 | $F | tail -1 | grep -q '^# matches'                # about 1 s
python3 h2_control.py --negative 2>&1 | $F | tail -1 | grep -q '^# matches'     # about 1 s
# not for CI (times.txt): h2_exhaust.py (8,300 to 11,600 s per worker, 3 workers), h2_hunt.py (100 to 5,000 s per range),
#   h2_design.py (300 to 3,300 s per range), h2_exhaust.py 14900 15000 0 1 (89 s)
```
