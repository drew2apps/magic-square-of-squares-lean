# table: every ratio 1000 < m ≤ 11000 that passes the congruences, with its closing step (run 5, F0)

**Experiment, not counted; tier D of `LEDGER.md`.** Script and table written by an AI model (Claude, the "Fable" session
of run 5); reviewed by a separate AI agent (`reviews/review_x4.md`: the bookkeeping holds, with the gaps below). Not
read by a person. Nothing here changes the census, which stops at m ≤ 1000.

`census_table.py` builds `CENSUS_TABLE_1000_11000.tsv` (one row per ratio: m, n, type, the first step that closes it, the
record, the label of that step) and `census_table.txt` (the counts) from records already committed. It runs no curve
computation and imports no PARI. Default: prints the summary and asserts that both files are unchanged; `--write`
rewrites them; `--negative`: for one ratio of each machine-checked step kind, every record of the ratio and of its
partner is removed and the ratio must come out UNCLOSED (7 of 7 refused, `census_table_negative.txt`).

**What it shows.** 20,682 ratios pass the congruences (other type 17,636; m² + n² a square 1,780; twice a square 1,266).
First closing step: class 17,636; second descent 2,582; repository L-value 23; rank 0 231; Demjanenko–Manin 175;
repository route 9; special 26. Every row has a step and a label; that is all "UNCLOSED 0" means. 20,681 rows have a
closing record of the tier their label states; (7735,4968) has none (label OPEN here; `../g4k/` excludes it only
conditionally on a rank being 1). The labels, not the counts, carry the tiers:
- **conditional on a rank being 1** (9): (3540,2581), (4440,2231), (6121,959), (6480,4961), (6671,2209), (6903,5104),
  (8249,3720), (10512,145), (10657,10367);
- **numerical L-value** (25): (6557,324), (6881,6233) (`../hard7k/`; `../lvalue/` now gives an exact enclosure
  conditional on PARI's `ellap`, not yet used by this table), and the 23 rows "repository (wide L-values)" for
  m ≤ 3000, PARI's value with a tail bound (`scripts/second_descent/wide*_lvalue*.jsonl`);
- **Demjanenko–Manin** (175): conditional on PARI's generators (saturation at primes ≤ 1000);
- **rank-1 route** (10): conditional on (a) the generator generating E1(ℚ) modulo torsion, (b) an AI reading of the
  constants of Checcoli–Veneziano–Viada; rank 1 there is PARI's upper bound over the isogeny class;
- every rank bound in every row is PARI's (unverified program output); the class step is the repository's criterion
  and four-class descent (tier D).

**Corrections made when this folder was committed (run 6).**
- The label of (6903,5104) said "plus (a), (b)"; its conditions are rank 1 of D1 of the partner (12007,1799) and the
  Demjanenko–Manin generator condition (as `../CENSUS7000.md` already says). Fixed in the script; the table and the
  summary were rewritten with `--write`; that row and one summary line are the only changes.
- The class step is now read from the repository's committed outputs (`scripts/cover_descent/tuple_scan_5000.txt`,
  `structure/class_selmer_15000_list.txt`, `tuple_selmer_15000.txt`). The session's copies of its rerun
  (`class_step/`, 2.1 MB) were byte-identical to these files on this range (its log) and are not committed; the
  table is the same either way (checked). `MSS_CLASS_DIR=DIR` builds the table from a rerun made by `class_step.sh`.
- `class_step.sh` ran the repository's programs through a PARI stand-in for sympy; it now calls them with sympy.
  `class_step/log.txt` is the session's log as it ran (PARI 2.19.1, the stand-in); its "summary diff exit 1" in step
  (2) is a bug of the first driver (it compared 3 of the 4 summary lines; `linalg_5000_10000.diff` shows the missing
  line); the 4 lines are identical. The review reran steps (3) and (4) for 10000 < m ≤ 11000 with sympy: identical
  (`../../referee_x4/samples_11k/`). Steps (1), (2) (about 60 min) were not rerun in review.
- The session's report gave "second descent 2,580"; the table has 2,582 (the counts sum to 20,682 only so).

**Gaps (stated, not fixed).**
- The 26 "special" rows (route, conditional, numerical, open) are typed into `SPECIAL` from the write-ups
  `CENSUS3000/5000/7000/9000/11000.md`; the script checks no record for them and `--negative` cannot refuse them. The
  review checked each against its route record and partner relation (`../../referee_x4/spot.py`).
- For other-type ratios with m > 5000, "closed by the criterion" means "absent from the criterion's list". That is
  sound because `class_selmer_linalg.py` runs over every coprime ratio of its range and lists each that keeps a class.
- `scanned` (`scan_certificates.txt`) is loaded and not used.

Files: `census_table.py`, `CENSUS_TABLE_1000_11000.tsv`, `census_table.txt`, `census_table_negative.txt`;
`pairs_2000_3000.json` and `gen_2000_3000.log` (the 2000–3000 congruence list, `gen_pairs.py 2000 3000` at 24d6ef6; the
review regenerated all six lists with separate code, `../../referee_x4/congr.py`: the same 20,682 ratios);
`class_step.sh` and `class_step/` (the session's log, commit, the two criterion summaries and the diffs).
