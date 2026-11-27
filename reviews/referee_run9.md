# Referee run 9 (hostile claims editor)

Start: Fri Nov 27 20:42:06 UTC 2026. Written by an AI model; not read by a person.

## ~20:45-21:00 UTC: setup, counts, cleanup

Scratch: /tmp/ref9/ (clean `git archive run9` in /tmp/ref9/tree). I edit no repository file. Repository at d073fd7
(run9), base run8 = abe7d7b (content-equal to public main 97c8564 except the mode of two referee_x3 .sh files).

### (1) Counts, own parser (/tmp/ref9/recount.py, written before reading z11's recount.py)
- 201 ratios (distinct), 1,766 tuples (= sum of header counts), 1,762 PAIRING-KILLED, 4 not: (3320,2303) x2,
  (5623,1017) x2, as stated. Kill lines 8,930 (= sum of per-curve 1280 + 4x1694 + 874); local -1 terms 18,222
  (1x5182, 3x2886, 5x826, 7x36); distinct (ratio, curve, class) kills 4,976; first killing primes as in NOTES;
  dim Sel^2 = 2 curves exactly (3320,2303) LU, (5623,1017) 08, (9191,7688) LU, (9532,4445) LU, (14885,6076) LU.
  Labels on (tuple, curve): only image / kernel / killed (no "outside Sel^2"). HOLDS.
- The 101: m from 465 to 4997, among the 997 of tuple_selmer_997.txt (291 <= m <= 5000): HOLDS as worded in LEDGER.
- The 836: of the 8,739 of class_selmer_15000_list.txt, whose m runs 5005..15000; the 836 have 5029 <= m <= 14977.
  LEDGER's "the 836 with m <= 15000" is wrong in scope (with m <= 15000 there are 101 + 836 = 937). SHOULD-FIX.
- CHANGES "201 of about 3,700 such ratios": not in NOTES or review z11, which say "201 of the 937 ratios that need a
  good prime" (G2's 101 + 836, m <= 15000). 3,700 appears to be 937 + H2's 2,807 (15000 < m <= 30000) = 3,744,
  a population the sample was not drawn from. MUST-FIX (two denominators for one sample in one run's docs).
- "(R)" in the LEDGER title and CHANGES: defined in no repository document (only used, undefined, in G1's
  REPORT_G1.md / g1_stress.py docstring as "the reduction (R)"); NOTES and the review say "a candidate global
  mechanism, not a proof of one" for the good-prime removals, not "for (R)". MUST-FIX (undefined, and stronger).

### (2) Cleanup
Method: applied /home/claude/conductor/patch_r7_referee.py ONCE to a clean archive of 060b1f1 (run 7 before the fixes)
-> I7; 3-way merged (60e01ad -> I7) into each run8 file; compared with run9. Result: g1, h1, h7, roadmap-h NOTES
identical; LEDGER, experiments/README, CHANGES, ci.yml differ only by run 9's own additions (and, for CHANGES, by
60e01ad's hand-written "After the pre-publication review" bullet, which my method strips and run9 rightly keeps); g3
NOTES differs by one blank line (run9 has one blank line after "(Run 6's range ...)", the single application would
leave two; renders the same). Whole-tree diff I7 vs 60e01ad: only these nine files plus 60e01ad's hand edits
(referee_z7/z9 README, reviews/README, CHANGES bullet). No sentence missing, none repeated; nothing else changed.
Copies in run8: LEDGER S1 4, ER g4k 4, CI line 4; all others 3 (CHANGES x2, ER block, g1, g3, h1 x3 bullets, h7,
roadmap-h). So CHANGES' "two to four times" should be "three or four times", and "one CI line ... three times" should
be "four times". NIT (but a stated count).

### (3) OUTPUTS totals
Own parser over all eight tables: 176 rows: A 71 + 5 part-A rows = 76, B 57, C 27, D 16. Run 9 adds 5 rows (A 2,
B 2, D 1). HOLDS. Row 1 time "about 25 s together" contradicts NOTES ("--negative 1 min, --h1control 5 s").
To check against the CI dry run.

### (4) CI, first findings
- YAML parses (PyYAML); the new lines are in the only step of `descents`, whose run starts `set -e -o pipefail`.
- **The `--summarize` line fails**: `t10_pairct.py` calls `PARI.allocatemem`, which prints "PARI stack size set to
  500000000 bytes, maximum size set to 4000002048" on STDOUT; the line has no `grep -v 'PARI stack'`, so diff
  exits 1 and the step fails (reproduced from my clean archive; with the grep it matches). MUST-FIX.
- Duplicate `cd "$GITHUB_WORKSPACE"` after the new block (harmless). NIT.
- z11 recount, wash_check, rank0, pmat rerun from the clean archive: all "# matches" (0.1 s, 0.4 s, 0.5 s, 4.2 s).

## ~20:50-20:57 UTC: CI dry runs, z11 spot checks, HEAD 6f1d3ee
- Conductor's first dry run (/tmp/r9ci/run.log): run 7 part passed (6 min); failed on the run-9 `--summarize`
  line (PARI stack notice), as I found independently. Fixed in 6f1d3ee (only change d073fd7..6f1d3ee: that one line
  gains `| grep -v 'PARI stack'`). HEAD YAML parses; the line is in the pipefail step. /tmp/r9ci/run9.log: BLOCK9_OK.
- My own run of the run-9 block (with the fix) from my clean archive, nice'd, alongside the conductor's:
  every line matches; 236 s total (t10 --negative 16 s, --h1control 4 s, --summarize 1 s: OUTPUTS "about 25 s"
  holds, NOTES "--negative 1 min" is pessimistic; z11 checks 215 s vs OUTPUTS "about 4.5 minutes": fine).
  Only file created: __pycache__ (untracked; `git diff --exit-code` unaffected).
- (5) referee_z11/README vs outputs: rank0, recount, wash_check, pmat, h1_mine, forge_t10, z11_compare (+mid),
  z11_pair, forge_mine all rerun and match; README claims (108 verdicts small set, 108 mid, 1,296 big, 960 batch30,
  gcd 8 over good p <= 199, matrices dims 2/6/5 with ranks 0/0/2, H1 kills on 0U/0L as h1_pairs.txt, nothing of a
  rational point killed) equal the committed .txt. HOLDS. (Independence: z11's Sel^2 and pairing are z7's code.)
- Controls (t10_negative.txt) and t10_deep.txt agree with NOTES. Root identification "alone unique" is
  (True, True) on all 1,206 curves of the 201; (True, False) on 8 control curves only.
- review_z11.md says "11 ratios ... (first 3 by t10, 6 middle and 5 last by me)": 3 + 6 + 5 = 14; z11 reran 11
  (oom_check.txt), so OUTPUTS' "11 of its ratios rerun identical by review z11" is right. Review is committed
  unedited; no change needed.

## Final (written ~20:57 UTC)
**Verdict: ready after fixes** (no number wrong; CI fixed in 6f1d3ee and passes from a clean archive).

MUST-FIX
- M1 (done in 6f1d3ee, verified): CI `--summarize` line needed `| grep -v 'PARI stack'`.
- M2 CHANGES, Ninth run: "201 of about 3,700 such ratios" is in neither NOTES nor review z11 (937 = 101 + 836,
  m <= 15000; 3,700 seems to add H2's 2,807 with m > 15000, not sampled). M3 "(R)" is defined in no repository
  document and is stronger than NOTES/z11 ("a candidate global mechanism, not a proof of one"). Replace
  "a candidate global mechanism for (R), not a proof (it does not reproduce the good-prime removals at (627,364),
  (991,263); 201 of about 3,700 such ratios)." by
  "a candidate global mechanism for the good-prime removals, not a proof (it does not reproduce the good-prime
  removals at the degenerate controls (627,364), (991,263); 201 of the 937 ratios with m ≤ 15000 that need a good
  prime)."
- M3 LEDGER title: replace "**Good-prime removals shadowed by the Cassels–Tate pairing (a candidate global mechanism
  for (R), not a proof)**" by "**Good-prime removals shadowed on the pair curves, by the Cassels–Tate pairing or by
  rank 0 (a candidate global mechanism for the good-prime removals, not a proof)**".

SHOULD-FIX
- S1 LEDGER: "and 100 of the 836 with m ≤ 15000)" -> "and 100 drawn with seed 1 from the 836 of the 8,739 ratios
  of `class_selmer_15000_list.txt`, 5000 < m ≤ 15000; 201 of the 937 with m ≤ 15000 that need a good prime)".
- S2 LEDGER: "with its own Kummer map and pairing on 38 ratios" -> "with its own forms, root roles and Kummer map
  and review z7's separate descent and pairing (`referee_z7`) on 38 ratios, 2,472 (tuple, curve) verdicts, all
  equal".
- S3 experiments/README row: "On 201 ratios, every tuple is removed by the pairing or by a pair curve of rank 0 (a
  candidate, not a proof)" -> "On 201 ratios that need a good prime, each of the 1,766 tuples left by the bad primes
  is also removed by the pairing (1,762) or by a pair curve of rank 0 (4); not at the degenerate controls
  (627,364), (991,263) (a candidate mechanism, not a proof; review z11: holds with gaps)".

NIT
- N1 CHANGES cleanup bullet: run 8 had each repeated passage 3 or 4 times and the CI line 4 times: "two to four
  times" -> "three or four times"; "one CI line of `descents` three times" -> "one CI line of `descents` four times";
  "eleven sentences" -> "eleven passages" (h1's S4 inserted two bullets).
- N2 experiments/README: "(z7–z11)" -> "(z7–z9, z11)" (there is no referee_z10).
- N3 ci.yml: `cd "$GITHUB_WORKSPACE"` now appears twice in a row after the run-9 block (harmless).
- N4 NOTES Times: "`--negative` 1 min" measured 16 s (OUTPUTS' 25 s is right).
- N5 t10_pairct.py docstring says mixed-pair base points come from hyperellratpoints; the code tries ±m/n, ±n/m
  first (as NOTES says) and the outputs use them.

Holds: every count (201; 1,766 = 598 + 1,168; 1,762; 4 by rank 0 on five rank-0 curves; 8,930 kill lines / 18,222
local terms; 4,976 distinct; per curve; first primes); "candidate"; no "first"; tier D; "not released"; the
degenerate caveat in LEDGER/CHANGES/NOTES; "101 of the 997 ..., 291 ≤ m ≤ 5000"; all seven z11 replacements
(a)-(g) applied in NOTES; cleanup exact (single application of patch_r7_referee.py, 3-way-merged, equals run9 up to
run 9's additions and one blank line in g3); OUTPUTS 176 = A 76 / B 57 / C 27 / D 16; CI YAML, pipefail, files.
Not checked: the manifest commit (not yet made); the 163 ratios no reviewer recomputed; G2's residual sets.
