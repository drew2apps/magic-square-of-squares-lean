# Referee, run 6 (pre-publication; hostile claims editor)

Started 05:21 UTC 27 Nov 2026. Hard stop 07:30 UTC. Repository /home/claude/conductor/int6, branch run6 = 19f58a5;
base dabf7b4 (= public 062c54f, release 0.6.0). I edit no repository file. Written by an AI model; not read by a person.

## Log

### 05:40 UTC — first pass
- OUTPUTS recount (script over every table with a Class column): 5 tables, 147 rows; A 57 + 5 mixed "A for part" = 62,
  B 44, C 28, D 13. New table: 30 rows = 15 A, 12 B, 2 C, 1 D; 117 + 30 = 147. HOLDS.
- CHANGES: no '-' line in the diff; the 0.6.0 entry is byte-unchanged; new entry headed "Unreleased, after 0.6.0". HOLDS.
- Θ₄ outputs: RESULT lines give Θ₄ = 1/2 (product −1) on 920,129 / 760,39 / 631,71 D1 (E′), bounds 1+2−2, 1+2−2,
  2+1−2 = 1. HOLDS.
- LEDGER tier C new row has FOUR cells in a 3-column table: the cell "not released" is dropped by GitHub (GFM ignores
  excess cells). Tier C and D headings still say "in the release tree as experiments". -> MUST-FIX (release metadata).
- Stale "Magma alone" counts found: README:119 "So four of those six still rest on Magma alone"; README:100,
  G1_FULL_SYSTEM:21, LEDGER:52 "three pairs rest on it alone" (122 tier; (631,71),(351,280) now have Θ₄);
  G1_FULL_SYSTEM:252 "These two lines rest on Magma alone"; G1:254 "Like the two lines above"; G1:343 "so three do";
  LEDGER:53 "so for four counted ratios with the partners, from Magma alone"; AUDIT_REQUEST:39 "the only three".
- CENSUS11000.md:42 cites REPORT_F0.md, not committed.

### 05:55 UTC — numbers traced
- Θ₄: dims and bounds from the RESULT lines; (945,248) fields from rank_fisher_theta4_gen.py header and
  selmer_square_5000.json (S contains 1356488820560645889749762; Sp contains ±5820647). "class group (2,2,2,2) under
  GRH" (LEDGER:134) is only in the conductor's reports/y5.md, not in a committed file (NIT).
- Fable part 2: census_table.txt 20,682; steps 17,636/2,582/23/231/175/9/26; 9 conditional; open (7735,4968);
  lvalue record.json L_y1 = [1.444085955623, 1.444085956081], Mmax 1,209,869,493, apcheck 20,000/20,000; review_x4
  line 107/262: one label fixed ((6903,5104)). experiments/README "50,600 primes" = 20,000 + 30,600. HOLDS.
- G1: g1_summary.txt: 17 ranges x 20,000 = 340,000; 63 kept; 3,400 controls, 3,399; but 3 of the 340,000 drawn are of
  the degenerate types ("degen" column 1+2) and are not "excluded" -> LEDGER/experiments README "each excluded" SHOULD-FIX.
  control_1109503.txt: 6 tuples, explicit points at 2 and every odd bad prime; 23,397 good primes tested. HOLDS.
- G2: g2_summary.txt 896+101 = 997, 7,903+836 = 8,739, largest good prime 181, all 1 mod 4; cross-check 57 -> 39/18.
  badalone_witness_all.txt: 937 ratios, 6,904 tuples, all witnessed. HOLDS.
- G3: g3_table.txt 611 = 364 square + 247 twice; 4 open, 3 conditional, route (12416,5313), (12487,4241). HOLDS.
- G4: g4_summary.txt 3,573; full 1,741 / partial 1,473 / none 359 / 23 no curve; 12 quotient curves; 64x. HOLDS
  (LEDGER row omits the coverage figures other than 1,741: SHOULD-FIX).
- MANIFEST: verify_manifest.py on a clean `git archive run6`: OK, 1308 files.
- YAML: parses (yaml.safe_load), descents timeout-minutes 150; f5 defined in the same step, set -e -o pipefail.
- CI lines rerun in a clean archive copy (/tmp/ref6/tree): block 1 (Fable part 2) passed in 92 s.
(Correction: the two section times above were estimates; `date -u` at this point reads 05:32 UTC. Another agent's
`quick` rerun (/tmp/qrun6) shares the 2 cores, so my CI timings are upper bounds.)

### CI rerun (item 5)
All 42 added commands of `descents` (ci.yml 580–634), extracted verbatim and run with `set -e -o pipefail`, `f5` as in
the job, GITHUB_WORKSPACE = a clean `git archive run6` copy, PARI 2.17.4, sympy 1.14, numpy 2.5.5, Python 3.13:
**ALL_OK in 448 s** (Fable block 92 s, Θ₄ block 199 s, roadmap-G block 157 s; 2 cores shared with another agent's
run). No file of the copy changed. So every referenced file exists and every diff is empty. YAML parses; the run-6 lines
are inside the step that defines f5 and sets pipefail. Timeout: 150 min vs ~104 min + 7.5 to 12 min = ~112 to 116 min:
enough (~34 min margin). HOLDS. NIT: the block comments say 2.5 / 4.5 / 2.5 minutes; measured 1.5 / 3.3 / 2.6.

### reviews/README (item 6)
Rows for review_x4, review_y5, review_x6 present and correct (y5's file calls itself "Review z5"; its scripts are in
/tmp/z5, not committed, and the row claims none). HOLDS. NIT: preamble says reports "committed on 26 November 2026".

## Findings

(In the replacements, the outer backticks delimit the old and new text; backticks inside are part of the text.)

Every number in the new CHANGES entry, the new LEDGER rows, the changed README/LEDGER/G1_FULL_SYSTEM/COMPLETE
sentences and the OUTPUTS totals traces to a committed output or to the review that recomputed it (sections above). No
number is wrong. The defects are stale counts, release labelling and three uncorrected session reports.

### MUST-FIX 1 (item 2): sentences that still give six, four or three as the current Magma-alone count
After run 6 the counted ratios resting on Magma alone for an upper rank bound: **one, (945,248)** (six before: (760,39),
(799,721), (945,248), (920,129) on the F0·F8·GL curves, (631,71), (351,280) on one D1 curve; five now also have Θ₄).
Within the 122, (724,665) (D1 of (1389,59)) still rests on Magma alone but is also closed by the four-point argument,
so it is not counted, as before. Exact replacements:
- a. README.md:119 `for each of the four curves on which six counted ratios rest.` ->
  `for each of the four curves on which six counted ratios rested on Magma alone until 26 November 2026.`
- b. README.md:119 `So four of those six still rest on Magma alone for it.` ->
  `So from 26 to 27 November 2026 four of those six still rested on Magma alone for it.`
- c. README.md:100 `within the 122 three pairs rest on it alone; one of the three, (724,665), is also closed by the four-point argument;` ->
  `within the 122 three pairs rested on it alone until 27 November 2026; since then the `D1` curve of two of them, (631,71) and (351,280), also has Fisher's `Θ₄` (tier C, not released), and the third, (724,665), which still rests on it alone within the 122, is also closed by the four-point argument;`
- d. G1_FULL_SYSTEM.md:21 `and three pairs rest on Magma alone for that bound within this method (one of the three, (724,665), is also closed in the four-point tier)` ->
  `and three pairs rested on Magma alone for that bound within this method until 27 November 2026 (one of the three, (724,665), is also closed in the four-point tier and still rests on Magma alone here; for the other two, (631,71) and (351,280), the `D1` curve also has Fisher's `Θ₄` since then, `rank_fisher_theta4_gen.py`, tier C, not released)`
- e. G1_FULL_SYSTEM.md:24 `and four of the curves rest on Magma alone for it, so six counted ratios do:` ->
  `and four of the curves rested on Magma alone for it until 26 November 2026, so six counted ratios did:`; and
  `so four counted ratios of this row rest on Magma alone;` -> `so four counted ratios of this row rested on Magma alone;`
- f. G1_FULL_SYSTEM.md:252 `These two lines rest on Magma alone for the upper bound and are marked `upper=magma`.` ->
  `These two lines rested on Magma alone for the upper bound and are marked `upper=magma`; since 27 November 2026 (not released) the bound 1 for this `D1` curve also comes from Fisher's `Θ₄` on its `E′` side (`rank_fisher_theta4_gen.py`, tier C), so they no longer rest on Magma alone.`
- g. G1_FULL_SYSTEM.md:343 `so those four rest on Magma for it, as three pairs of the 122 already do (since 26 November 2026 the curve of (72,65) also has bound 1 from Fisher's `Θ₄`, `rank_fisher_theta4.py`, tier C, so three do).` ->
  `so those four rested on Magma alone for it, as three pairs of the 122 did (since 26 November 2026 the curve of (72,65) also has bound 1 from Fisher's `Θ₄`, `rank_fisher_theta4.py`, tier C, so three did; since 27 November 2026 only `E⁺` of (945,248) does, see the end of this item).`
- h. G1_FULL_SYSTEM.md:365 `except for (631,71), (351,280) and (724,665), where it comes from Magma alone (for the `D1` curve of (631,71), (351,280) also from Fisher's `Θ₄` since 27 November 2026, tier C).` ->
  `except for (724,665), where it comes from Magma alone, and for (631,71), (351,280), where it came from Magma alone until 27 November 2026 and since then also comes from Fisher's `Θ₄` (`rank_fisher_theta4_gen.py`, tier C, not released).`
- i. LEDGER.md:52 `(it confirms every upper rank bound; three pairs rest on it alone)` ->
  `(it confirms every upper rank bound; three pairs rested on it alone until 27 November 2026; since then (631,71), (351,280) also have Fisher's `Θ₄` (tier C, not released), and (724,665), also closed in the four-point tier, is the one left)`
- j. LEDGER.md:53 `so for four counted ratios with the partners, from Magma alone;` ->
  `so for four counted ratios with the partners, from Magma alone until 27 November 2026;`; and
  `from Magma alone for the curves of (760,39), (945,248) and (920,129); for the curve of (72,65) from Magma and from Fisher's `Θ₄`, tier C; since 27 November 2026 `Θ₄` for (760,39) and (920,129) too)` ->
  `from Magma alone for the curve of (945,248); for the curves of (72,65), (760,39) and (920,129) from Magma and from Fisher's `Θ₄`, tier C (the last two since 27 November 2026, not released))` (review_y5 FIX 3, not fully done).

### MUST-FIX 2 (item 4): the run-6 LEDGER rows are not visibly marked unreleased
- LEDGER.md:65 has four cells in a three-column table; GitHub drops excess cells, so `not released` is invisible.
  Delete the trailing ` not released |` (end the row after the third cell) and change
  `(since 27 November 2026; the curves of the counted closures` -> `(since 27 November 2026, after release 0.6.0: not released; the curves of the counted closures`.
- LEDGER.md:60 and :82 headings say "in the release tree as experiments", false for the run-6 rows. Replace
  `(not counted; in the release tree as experiments)` (both lines) with
  `(not counted; in the release tree as experiments, except the rows marked not released)`.
- The five new tier D rows (LEDGER.md:103–107): replace (both occurrences) `(run 6; a separate AI session;` ->
  `(run 6, after release 0.6.0, not released; a separate AI session;`; `no law found (reviewed: holds) |` ->
  `no law found (run 6, after release 0.6.0, not released; reviewed: holds) |`;
  `other types not covered (reviewed: holds with gaps) |` ->
  `other types not covered (run 6, after release 0.6.0, not released; reviewed: holds with gaps) |`;
  `| Search (numerical, not a certificate): on the 3,573` ->
  `| Search (numerical, not a certificate; run 6, after release 0.6.0, not released): on the 3,573`.
- Also LEDGER.md:53 `for the curves of (760,39) and (920,129) also from `Θ₄` since 27 November 2026,` ->
  `for the curves of (760,39) and (920,129) also from `Θ₄` since 27 November 2026 (after release 0.6.0, not released),`
  (that row sits under the heading "release 0.6.0").
CHANGES (0.6.0 entry byte-unchanged; new entry "Unreleased, after 0.6.0") and CITATION.cff (0.6.0) are right.

### MUST-FIX 3 (item 3): three committed session reports keep claims stronger than the tiers, with no warning in the file
experiments/roadmap-g/g1/REPORT_G1.md:15 "It does show the descent can leave classes that no point explains" (and :60
"2·(5·7129)²", :80 "no prime left out"); g4/REPORT_G4.md:9 "No solution found. On every coprime square-type ratio";
g2/REPORT_G2.md:9 "is decided". NOTES.md corrects them, but each REPORT is new in run 6 and reads alone. Keep the text,
insert after line 1 of each:
- REPORT_G1.md: `> **Read `NOTES.md` first.** This is the session's report, committed unchanged; where they differ, `NOTES.md` (after review x6) holds: `m² + n² = 2·(5³·7129)²`; "no prime left out" is false for G1; the four extra tuples of (1109503,597671) are locally soluble everywhere, and whether a rational point has one of these classes is open ("classes that no point explains" is not shown).`
- REPORT_G4.md: `> **Read `NOTES.md` first.** This is the session's report, committed unchanged; where they differ, `NOTES.md` (after review x6) holds: "no solution" means none of the points listed gives one; 23 ratios had no curve searched; a search, not a certificate.`
- REPORT_G2.md and REPORT_G3.md: `> **Read `NOTES.md` first.** This is the session's report, committed unchanged; where they differ, `NOTES.md` (after review x6) holds.`

### SHOULD-FIX
- S1 (release dates). 0.6.0 is dated 27 November 2026, so "since 27 November 2026" does not tell a reader that the
  change is after the release. README.md:100 `(eight until 26 November 2026, six until 27 November)` ->
  `(eight until 26 November 2026, six until 27 November 2026 and in release 0.6.0; the change is not released)`;
  README.md:102 `and since 27 November 2026 the curves of (920,129)` -> `and since 27 November 2026 (after release 0.6.0; not released) the curves of (920,129)`;
  README.md:119 `since 27 November 2026 `rank_fisher_theta4_gen.py` computes` -> `since 27 November 2026 (after release 0.6.0; not released) `rank_fisher_theta4_gen.py` computes`;
  G1_FULL_SYSTEM.md:24 and :343 `since 27 November 2026 the curves of (920,129), of (760,39)` -> `since 27 November 2026 (not released) the curves of (920,129), of (760,39)`;
  COMPLETE.md:197 `**No** (computed; run 6,` -> `**No** (computed; run 6, after release 0.6.0, not released;`.
- S2. LEDGER.md:134 "The one curve whose rank 1 rests on Magma alone" leaves out D1 of (1389,59) (the (724,665) line;
  not among rank_fisher_curves.json's five). Replace
  `The one curve whose rank 1 rests on Magma alone, `E⁺` of (945,248) (the curves of (72,65), (760,39), (920,129) and the `D1` curve of (631,71) have `Θ₄`):` ->
  `The one curve whose rank 1 rests on Magma alone for a counted ratio that no other method closes, `E⁺` of (945,248) (the curves of (72,65), (760,39), (920,129) and the `D1` curve of (631,71) have `Θ₄`; the `D1` curve of (1389,59), which closes (724,665) in the 122, also rests on Magma alone, but (724,665) is closed again by the four-point argument):`
- S3. G1_FULL_SYSTEM.md:254 `Like the two lines above, this one rests on Magma alone for the upper bound.` ->
  `Like the two lines above until 27 November 2026, when their curve got Fisher's `Θ₄`, this one rests on Magma alone for the upper bound.`
- S4. AUDIT_REQUEST.md:39 (stale since 25 November, more so now) `are the only three that rest on Magma alone; the last one is also closed by the method of closure 4, with PARI's bounds.` ->
  `are the only three lines of the 122 that rested on Magma alone; the last one is also closed by the method of closure 4, with PARI's bounds, and since 27 November 2026 (not released) the `D1` curve of the first two also has Fisher's `Θ₄` (`scripts/g1_full/rank_fisher_theta4_gen.py`, tier C). For the whole census see `README.md`: one counted ratio, (945,248), rests on Magma alone.`
- S5. 3 of the 340,000 G1 draws are of the degenerate types (g1_summary.txt "degen" 1, 2) and are not "excluded".
  LEDGER.md:104 `(20,000 per dyadic range, `10⁴ ≤ m ≤ 10⁹`): each excluded by` ->
  `(20,000 per dyadic range, `10⁴ ≤ m ≤ 10⁹`; 3 of them of the degenerate types): each of the other 339,997 excluded by`;
  experiments/README.md:127 `` `10⁴ ≤ m ≤ 10⁹`): each is excluded by`` -> `` `10⁴ ≤ m ≤ 10⁹`; 3 of them of the degenerate types): each of the other 339,997 is excluded by``.
- S6. LEDGER.md:103 (census table row) says what the 9 conditional and 2 L-value rows rest on but not the rest.
  `(700,529)'s big `L`-value reruns to its committed record |` ->
  `(700,529)'s big `L`-value reruns to its committed record; every rank bound in the table is PARI's, the 175 Demjanenko–Manin rows rest on PARI's generators, the 10 route rows on (a), (b), 23 rows on PARI's numerical `L`-value with a tail bound, and the table keeps the numerical label of (6557,324), (6881,6233) |`
- S7. LEDGER.md:107 `(1,741 ratios with a full search on E+ or E−)` ->
  `(wherever PARI's rank bounds meet, rank ≤ 3: 1,741 ratios with a full search on E+ or E−, 1,473 partial, 359 neither, 23 no curve searched)`.
- S8. LEDGER.md:65, third cell, the main gap review_y5 named is missing:
  `the four outside controls through the general code (BSD-conditional);` ->
  `the four outside controls through the general code (BSD-conditional; like the (72,65) control, all in ℚ(√−2), ℚ(√2) with a paired class ≡ 1 mod 8, so no control with a known answer exercises another field or the 2-adic symbol);`
- S9. scripts/g1_full/rank_fisher_theta4_gen.py:15–16 and :47 (review_y5 FIX 2 not done): lines 15–16
  `Only side E' runs in CI and in the committed outputs for (631,71); for the other curves the side is the one of rank_fisher.py.` ->
  `Every committed output and CI line reduces on the E' side (the lifted class on E', as for (72,65)): for (631,71) by --side=E', for the other curves because that is rank_fisher.py's side. The branch for a reduction on the E side (the else of prepare()) is run by no output and no CI line: untested, and its K2 list holds only the class of T'.`;
  line 47 `for a reduction on the E side the roles of E and E' are` -> `for a reduction on the E side (untested, see above) the roles of E and E' are`.
- S10. experiments/census-wide/CENSUS11000.md:42 `(see REPORT_F0.md)` -> `(see `table/NOTES.md`)` (REPORT_F0.md is not committed).

### NIT
- N1. OUTPUTS.md:207 lists `reviews/review_x4.md` (class C) but not review_y5/x6 or any earlier review: delete the row
  (then the totals read 146 = 62 A / 44 B / 27 C / 13 D) or add rows for the other two (149 = 62/44/30/13).
- N2. README.md:102 `so four do,` -> `so four did,`.
- N3. REFERENCES.md:131 `| The route for the ratios that rest on Magma alone |` -> `| The route for the ratios that rested on Magma alone (Θ₄; one, (945,248), is left) |`.
- N4. CHANGES.md:8 `(the review found explicit local points at every bad prime and gave the argument at the good primes)` ->
  `(the review found explicit local points at every bad prime; at the good primes, `exact_local.py` to 600,000 and a short argument above)`.
- N5. LEDGER.md:134 "class group (2,2,2,2) under GRH" is only in the conductor's reports/y5.md; no committed output.
- N6. reviews/README.md preamble "committed on 26 November 2026": add "(run 6's on 27 November)".
- N7. CI block comments' times (2.5 / 4.5 / 2.5 min) vs measured 1.5 / 3.3 / 2.6; CHANGES' "about 12 minutes" is safe.

### Checked and holding
OUTPUTS totals; CHANGES 0.6.0 unchanged; every Θ₄ value and bound; (945,248) fields; Fable table and L-value numbers;
G1–G4 numbers; COMPLETE §7.6 numbers and tier; tiers: Θ₄ tier C everywhere, nothing "proved"; all experiments tier D;
(6557,324) "conditional on PARI's ellap" in CHANGES, LEDGER, experiments/README, lvalue/NOTES; (7735,4968) open /
conditional on rank 1; G4 "a search, not a certificate"; (1109503,597671) "locally soluble everywhere; whether a rational
point has one is open" in CHANGES, LEDGER, COMPLETE, NOTES, experiments/README (only the REPORT_G1 text is stronger, MF3);
no "Sha" for it anywhere outside REPORT_G1:76 ("I haven't checked whether they come from a Selmer element ... (open)",
which is not a claim). CI: all 42 lines pass, YAML valid, timeout enough; MANIFEST verifies (1308 files); reviews/README
rows right. After the fixes the MANIFEST must be regenerated.

## Final (05:39 UTC)
Bottom line: **not ready as is; ready after doc-text fixes only.** No number is wrong, CI passes, the tiers are right in
CHANGES, LEDGER and the NOTES. MUST-FIX: (1) ten stale "Magma alone" counts (six/four/three as current) in README,
G1_FULL_SYSTEM, LEDGER; (2) the run-6 LEDGER rows are not visibly unreleased (the tier C row's fourth cell is dropped by
GitHub; headings say "in the release tree"); (3) REPORT_G1/G2/G4 keep overstatements with no banner. 10 SHOULD-FIX,
7 NIT. I edited no repository file.
