# REPORT_I3 (part 1, 27 Nov 2026 about 20:15 UTC): phantoms with every bad prime

> **Read `../NOTES.md` first.** This is the session's report, committed unchanged apart from this line; where they differ, `../NOTES.md` (after review z17) holds.


Written by an AI model; not read by a person. Experiment, not counted. **Part 1:** the pairs and the rescan of H1's 20,009 ratios are done. The 50,014 new ratios are running and follow in the addendum.

## Bottom line so far

- **7 of H1's 12 pairs are true phantoms, and 5 are not.** With every bad prime tested (the primes of 2mn(m² − n²), of m² + n² and of the octic, each factor proved prime), 14 of the 24 ratios still keep extra tuples. The other 10 (5 pairs) lose every extra tuple at a prime of the octic that H1 never tested. Both sides of each pair lose them at the same prime. *(computed, exact)* This agrees with review z8, and was found independently.

  | Pair | Killed at (octic prime) |
  |---|---|
  | (11875080,6497849) / (18372929,5377231), H1's "open" pair | 5462461 |
  | (32396425,3694248) / (36090673,28702177) | 218077 |
  | (175775580,15972419) / (191747999,159803161) | 10057297 |
  | (234319635,115640732) / (349960367,118678903) | 21062406847889269 |
  | (518804652,31743595) / (550548247,487061057) | 1650997632721 |

  The 7 true pairs are (204724,97893)/(302617,106831), (413004,148603)/(561607,264401), (756492,255635)/(1012127,500857), (853587,255916)/(1109503,597671), (1992873,999664)/(2992537,993209), (189089400,156063271)/(345152671,33026129) and (231641256,72371617)/(304012873,159269639). H1 settles all 7 (rank-0 quotient or pair curves), so no pair is open.
- **The rescan of H1's 20,009 ratios is complete and agrees with H1 except for the 5 false pairs.** Per range k = 13 to 29 the phantom counts are 0, 0, 0, 0, 1, 0, 3, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0 (H1: the same up to k = 22, then 1 at k = 23, 1 at 25, 4 at 27 and 2 at 28). The 6 phantoms that disappear are H1's (11875080,6497849), (36090673,28702177), (175775580,15972419), (191747999,159803161), (349960367,118678903) and (518804652,31743595), all ratios of the 5 false pairs (their other ratios were not drawn by H1's scan). At k = 27 the 2 phantoms left are (189089400,156063271) and (231641256,72371617), both true pairs. No new phantom appears: testing more primes can only remove them. *(computed, exact)*

## How it was checked

- i3_scan.py runs the four-class descent (exact_local.descend_exact) on the supports' primes, every prime of m² + n² and of the octic (PARI factors, each factor checked with a primality proof), and good primes to 150, 1500, 5000 and 20000, while more than the known tuples is left.
- i3_negative.txt: with the new primes left out, the phantoms come back as in H1 ((11875080,6497849) among them). A forged extra tuple at (1109503,597671) is refused.
- PARI 2.19.1 here; CI has 2.17.4. Only factorisations come from PARI, and the printed sets don't depend on its order.

## Still to do

- 50,014 new degenerate ratios (seed 5000 + k, 2942 per range k = 13 to 29), as the updated prompt asks: i3_scan.py k --new, run_i3_new.sh. That covers 50,014, not exactly 50,000, so that each range gets the same count.
- Settling any new pair, and explicit local points for the extra tuples (z8's witness scripts would help; they aren't in the repository yet).

## Files

i3_scan.py, run_i3.sh, run_i3_new.sh, i3_pairs.txt, i3_negative.txt (rewritten; `--negative` prints it exactly), i3_scan_k13 to k29.txt, times.txt (pairs 449 s; ranges 88 s at k = 13 up to 4,828 s at k = 29, 16,595 s in all).

## Proposed CI lines (draft)

```
python3 experiments/roadmap-i/i3/i3_scan.py --negative   # 27 s
python3 experiments/roadmap-i/i3/i3_scan.py 13           # 88 s (one range as a sample; the others take up to about 5 minutes)
```
