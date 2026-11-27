# G2: which primes empty the four-class set?

> **Read `NOTES.md` first.** This is the session's report, committed unchanged; where they differ, `NOTES.md` (after review x6) holds.

Written by an AI model (Claude); not read by a person. Experiment, not counted. Nothing here changes the census, which stops at m ≤ 1000.

Run at 062c54f (0.6.0 release candidate); the programs used are identical at 24d6ef6. Real sympy 1.14.0 is installed; every factorisation here is made by PARI, and each factor is proved prime by `isprime`. Machine: PARI 2.19.1 and cypari2, not CI's PARI 2.17.4.

## Bottom line

- **Test 6 of COMPLETE.md §7 is decided: no.** The bad primes alone do not always empty the four-class set. They do for 896 of the 997 and 7,903 of the 8,739. For the other 101 and 836 (937 of 9,736), they leave tuples that no bad prime excludes. *(computed, exact: every bad prime is tested exactly, the large ones of m² + n² and of the octic included, and each tuple left has an explicit 2-adic point)*
- **A few good primes always finish the job.** Every one of the 9,736 is emptied by its bad primes plus good primes up to 101 (the 997) or 181 (the 8,739). No ratio is left non-empty. *(computed, exact)*
- **No law found.** Every good prime ever needed is 1 mod 4, but that is forced: at a good prime that is 3 mod 4, every divisor has the same label. No residue-symbol feature tried separates the ratios that bad primes empty from those they don't (table below). The share that bad primes fail on falls as the tuple space grows, from about 23% at 2¹⁶ to 6% at 2²⁶, and stays well above zero. *(computed; no conjecture stated)*

## Kill sets (g2_summary.txt)

For each ratio, g2_killsets.py records an irredundant kill set: no prime of it can be dropped. It is not proved to be of least size. The set is built from all bad primes, then the shortest run of good primes in increasing order. Then each prime, largest first, is dropped if the rest still empties the set.

| | the 997 | the 8,739 |
|---|---|---|
| emptied by bad primes alone | 896 | 7,903 |
| need a good prime too | 101 | 836 |
| not emptied | 0 | 0 |
| largest good prime needed | 101 | 181 |
| kill-set sizes | 1 to 10, mostly 3 (452) | 1 to 12, mostly 3 (3,277) |

Largest good prime needed, with the number of ratios:
- **997:** 17 (5), 29 (42), 37 (24), 41 (10), 53 (6), 61 (6), 73 (2), 89 (1), 101 (5).
- **8,739:** 17 (36), 29 (359), 37 (158), 41 (98), 53 (47), 61 (62), 73 (14), 89 (17), 97 (14), 101 (16), 109 (5), 113 (1), 157 (4), 173 (3), 181 (2).

Bad primes not enough, by log₂ of the tuple space (ratios there, failures):
- **997:** 14: 2, 0; 16: 39, 9; 18: 230, 30; 20: 435, 38; 22: 246, 18; 24: 45, 6.
- **8,739:** 14: 2, 0; 16: 78, 21; 18: 807, 163; 20: 2,409, 252; 22: 3,453, 270; 24: 1,749, 116; 26: 235, 14; 28: 6, 0.

Bad primes used most: 2, 5, 13, 17 and 3, in that order for both lists.

**G1's 63 kept ratios** (m from 10⁴ to 10⁹; g2_g1_kept_0.txt): 60 were run and 3 ran out of memory at 8 GB (listed as SKIPPED). Of the 60, bad primes alone empty 57, and 3 need one or two good primes (largest 61). The pattern is the same as for the two lists above. (A first run at 4 GB lost more ratios to memory. The 8 GB run was stopped by a 30-minute limit on my side after 60 ratios and resumed for the last 3. times.txt has the 4 GB run and the resumed part; the stopped part has no time.)

## What holds and how it was checked

- **Bad primes:** every prime of 2mn(m² − n²), of m² + n² and of the octic m⁸ + 14m⁴n⁴ + n⁸, factorised by PARI and each proved prime.
- **Local test:** exact at every odd prime (exact_local.py from G1, no prime left out, primes of any size handled by the node argument). At 2, the repository's test at depth 12. Every known point's tuple is asserted to pass at every prime used.
- **Test 6 made exact (g2_badalone_witness.py):** a tuple left by the odd bad primes is in each exact local set by construction. At 2, the repository's test is one-sided, so each tuple left is also given an explicit 2-adic point, (u, v) = (t, 1) or (1, 2t) with t < 2¹⁴, at which the four forms' classes are the tuple's (or 3 times them). All 101 + 836 ratios: every tuple left is witnessed (g2_badalone_witness_g2_997_all.txt, g2_badalone_witness_g2_15000_all.txt).
- **Cross-check with the repository's tuple_badprimes_997.txt** (bad primes, with those of m² + n² and the octic only below 2000, tuple space ≤ 2¹⁸):
  - the same 271 ratios;
  - it left 57; with the large bad primes too, 39 are left, so the large bad primes empty 18 more;
  - no ratio has more tuples left here than there, as it must be.
- **Feature scan (g2_features.py, g2_features_997.txt and g2_features_15000.txt):**
  - the features are over the primes q ≡ 1 mod 8 dividing m² + n² to an odd power: whether there is one; q mod 16; whether 2, m, n, m ± n or mn is a fourth power mod q; whether every prime of mn(m² − n²) is a square, or a fourth power, mod q;
  - none, nor its negation, forces either outcome on 5 ratios or more;
  - in the 8,739, "every prime of mn(m² − n²) is a fourth power mod some q" holds for 6 ratios: 5 emptied and 1 left.
- **Negative mode (g2_check.py --negative, g2_negative.txt):**
  - 3 degenerate ratios come out NOT EMPTIED, as they must;
  - for 8 ratios of the 997 with kill sets of 2 or more primes, each kill set is recomputed from scratch and empties the set, and dropping any one of its primes leaves a tuple.
  So a forged, smaller kill set is refused.

## What failed

- An O(p) search for √−1 in the repository's `Local.__init__` hung at p = 22,684,758,289. exact_local.py now sets those fields itself and takes √−1 from PARI (see G1).
- No law. I have no conjecture to state.

## Next step

A uniform statement would need a table at good primes ≡ 1 mod 4 up to about 200, like Theorem C′, for the (ℤ/2)⁴-cover, as COMPLETE.md §7 says. The data here say only that such primes were enough for these 9,736 ratios.

## Proposed COMPLETE.md §7, test 6 (replacing "The test is inconclusive for the 57")

> Run with every bad prime, the large primes of m² + n² and of the octic included (exact local test at a node; experiments/roadmap-g/g2), the bad primes alone empty 896 of the 997 and 7,903 of the 8,739 of class_selmer_15000_list.txt. For the other 101 and 836, every tuple left has a point at 2 and passes the exact test at every odd bad prime, so **the bad primes alone do not suffice**. Bad primes and good primes up to 181, all 1 mod 4, empty all 9,736. No residue-symbol feature tried separates the two cases. Computed (exact); experiment, not counted. Written by an AI model; not read by a person.

## Proposed LEDGER / CHANGES text

> **Experiment, not counted.** Kill sets for the four-class descent (experiments/roadmap-g/g2): for the 997 and the 8,739, an irredundant set of primes that empties the four-class set, with every bad prime tested exactly. Test 6 of COMPLETE.md §7 is decided negatively: the bad primes alone leave tuples, each locally possible at every bad prime, for 101 of the 997 and 836 of the 8,739. Good primes up to 181 finish every case. No law was found. Computed (exact). Written by an AI model; not read by a person.

## Proposed CI lines

```
cd experiments/roadmap-g/g2        # F = grep -v '^PARI stack size'
python3 g2_check.py 2>&1 | $F | tail -1 | grep -q '^# matches'                          # under 1 s
python3 g2_features.py g2_997_all.txt 2>&1 | $F | tail -1 | grep -q '^# matches'        # under 1 s
python3 g2_features.py g2_15000_all.txt 2>&1 | $F | tail -1 | grep -q '^# matches'      # about 1 s
python3 g2_check.py --negative 2>&1 | $F | tail -1 | grep -q '^# matches'               # about 6 s
python3 g2_badalone_witness.py g2_997_all.txt 2>&1 | $F | tail -1 | grep -q '^# matches'  # about 40 s
# not for CI: bash run_g2.sh i N   (the 997 take about 3 min per worker on two workers; the 8,739 about 57 min each; times.txt)
# not for CI: bash run_g1kept.sh   (G1's 63 kept ratios, one process each under 8 GB; times.txt)
```

g2_check.py finds the repository through MSS_REPO, or else `../../..`.
