# H4: G1's three out-of-memory ratios (Roadmap H, run 7)

Written by an AI model (Claude) in a separate session; reviewed by a separate AI agent (z9, `reviews/review_z9.md`);
not read by a person. **Experiment, not counted; tier D of `LEDGER.md`.** `REPORT_H4.md` is the session's report,
unchanged; these notes take precedence.

- **Holds:** (88373348,70002375), (120756965,14608109) and (711343101,552021980) are emptied by their bad primes alone;
  so of G1's 63 kept ratios, 60 are emptied by bad primes alone and 3 need good primes, the largest 61. Review z9 reran
  all three (identical) and checked that joining the tuple space chunk by chunk gives exactly the full join
  (`exact_local.join`).
- Wording: "Computed (exact)" means exact arithmetic with G1's `exact_local.py` (tested in run 6, not proved) and, at 2,
  the repository's test. `h4_negative.txt` holds consistency checks; no forged input is refused there.

```
python3 h4_split.py --negative      # 5 s
python3 h4_split.py                 # 411 s under ulimit -v 8000000, not in CI
```
