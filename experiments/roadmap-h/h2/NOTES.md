# H2: survivor hunt (Roadmap H, run 7)

Written by an AI model (Claude) in a separate session; reviewed by a separate AI agent (z9, `reviews/review_z9.md`, its
own checks in `../../referee_z9/`); not read by a person. **Experiment, not counted; tier D of `LEDGER.md`.** Nothing
here changes the census, which stops at `m ≤ 1000`. `REPORT_H2.md` is the session's report, unchanged; the notes below
take precedence over it.

## What it claims (as corrected)

- **Exhaustive, `15000 < m ≤ 30000`:** of the 205,172,879 coprime ratios `m : n`, `0 < n < m`, not of the degenerate
  types (`m² + n²` a square or twice a square), 30,434 keep a class under the repository's criterion
  (`class_selmer.empty`, `class_selmer_linalg`), and the four-class descent with **every** bad prime exact and the good
  primes up to 229 empties every one; 2,807 need a good prime. No survivor. Review z9 recounted the outputs, recomputed
  the degenerate count, reran 63 sampled ratios (the three with the largest cutoff among them) and two whole values of
  `m` with sympy factorisations and G2's code, and reran the control block `14900 < m ≤ 15000`: all identical.
- **Two random hunts beyond that range:** a biased one (`10⁴ ≤ m < 2³⁷`; ranges up to `10¹²` were not reached; the last
  class-keeping ratios are at `k = 33`) and one steered by H3's tables (`2¹⁶ ≤ m < 2²⁹`) add **2,131 distinct**
  class-keeping ratios with `m > 30,000` (not 2,421: the biased hunt's 2,251 lines are 2,234 distinct ratios, 273 of them
  with `m ≤ 30000`), all emptied, the largest cutoff 137. The descent of one more ratio, (117702451,107508765), did not
  finish; three ratios at `k = 36` hit the class criterion's time limit (not counted as kept).
- The largest good prime needed per block of 1000 lies between 101 and 229 (149 in the first block, 229 in
  25000–26000, 181 in the last); it is not monotone in `m`.
- **Withdrawn:** "by H3, a cutoff above 1153 is never needed", "P ≤ 1153", "a survivor would have to keep a tuple through
  all 96 good primes up to 1153", "about 2.5 core-hours per block" (times.txt gives about 0.6). H3's Hasse–Weil step is wrong (`../h3/NOTES.md`; the docstring of `h2_design.py` says the same, "exactly the p = 1 mod 4 with 17 <= p <= 1153": read it as "the primes of H3's tables"): good primes above 1153 do remove vectors, so a survivor test must go to 67,583. The kill
  sets found here do not depend on H3.
- Computed (exact arithmetic; the local test at odd primes is G1's `exact_local.py`, tested in run 6, not proved; at 2
  the repository's one-sided test). "No survivor" is a statement about the ratios run.

## Rerunning
```
python3 h2_summary.py               # 1 s
python3 h2_control.py               # 1 s; --negative 1 s
```
`h2_control.py` reads G2's list of class-keeping ratios from `scripts/cover_descent/structure/class_selmer_15000_list.txt`
when G2's folder has no copy (one change made at integration). The long runs (`h2_exhaust.py`, about 0.6 core-hours
per block of 1000; `h2_hunt.py`, `h2_design.py`) are not in CI.
