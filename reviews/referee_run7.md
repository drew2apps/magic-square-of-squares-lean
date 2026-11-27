# Referee report, run 7 (hostile claims editor)

Start: Fri Nov 27 15:00:58 UTC 2026. Reviewer of `git diff run6..run7` in /home/claude/conductor/int7 (branch run7, HEAD 060b1f1).
Written by an AI model; not read by a person. I edit no repository file; scratch work in /tmp/referee7/.

## Log

15:10 Read BRIEF_COMMON7, the three reviews (z7, z8, z9), the diffs of LEDGER/CHANGES/OUTPUTS/SECOND_DESCENT_CHAIN/
experiments README/ci.yml/g3_sample.py/reviews README, roadmap-h NOTES + HANDOFF, h1/h2/h3/h4/h6/h7 NOTES, g3 NOTES.
15:15 Started the run-7 CI block (lines 635-656 of ci.yml) under `set -e -o pipefail` from a clean `git archive run7`
in /tmp/referee7/ws (git-initialised to detect rewrites), nice 10, with per-line timestamps.

Early findings (to be confirmed):
- REPORT_H5.md and REPORT_H8.md have no "Read NOTES.md first" banner (all others do). REPORT_H8 still says
  "numerical L(E,1)" and omits the CT pairing from what H7 rests on (z7 asked for both changes).
- Stale older text (run 7 settles or changes it): LEDGER row "Random ratios to m = 10^9" and g1/NOTES.md and
  experiments/README.md (G1 bullet) still say for (1109503,597671) "whether a rational point has one of these classes is
  open" (H1, reviewed by z8, settles it twice); LEDGER row "Beyond the census, 1000 < m <= 9000" ("(7735,4968), open";
  "(6557,324), (6881,6233) on a numerical L(E,1)"), LEDGER row "Census table 1000 < m <= 11000" ("(7735,4968) ...
  conditional on rank 1 by the g4.K route"), experiments/README.md run-5 F4 bullet ("Rank 1 is PARI's ellrank bound
  only"), census-wide/g4k/NOTES.md ("Rank 1 is not proved"), census-wide/table/NOTES.md, README.md line 115,
  experiments/README.md roadmap-g row ("11000 < m <= 13000 (4 open, 3 conditional)"), g3/NOTES.md head.
- (13348,3795)'s rank 1 (in LEDGER) has no CI line: h7_b13000_rank1.txt is class B though it takes 6 s.
- h1/NOTES and LEDGER "explicit points to 67,583": z8's own good-prime test (phantom_mid) covers p < 1200 and
  20000 < p < 67583 only; 1200 < p < 20000 rests on G1's exact_local (the scan, full_descent), "tested, not proved".

15:40 CI (item 4). The run-7 block of `descents` (ci.yml lines 635-656, verbatim), under `set -e -o pipefail`, from a
clean `git archive run7` (PARI 2.17.4, sympy 1.14.0, numpy 2.5.5): **exit 0 in 7 min 10 s** (nice 10, 2 shared cores);
`git status` of the copy clean afterwards (no file rewritten). Longest lines: h7_e04_rank 89 s, h1_iso 67 s,
h1_pairs --negative 44 s, run_ct_g4k 23 s, h3_tables 20 s. YAML parses (yaml.safe_load); the lines are in the one
`run:` step of `descents` that starts `set -e -o pipefail` and defines f5; numpy is installed earlier in that step.
Timeout 150 min against ~116 min after run 6 (referee_run6) + 7 = ~123 min: enough, margin ~27 min.
FAIL (harness): `hm` checks only the last line. `h3_tables.py` prints "# table matches h3_tables.json" / "# TABLE
DIFFERS" before its last line, so a forged `h3_tables.json` (one entry at p = 41 changed) passes `hm` (tested: prints
"# TABLE DIFFERS", then "# matches", hm returns 0), while OUTPUTS calls h3_tables.json class A in descents. A stricter
hm (fail on any "DIFFERS" line) fails on the forged json; on the clean archive all 24 hm calls pass with it (4 min).
(13348,3795): `h7_g4k_rank1.py 13348,3795 13899,1180` reproduces h7_b13000_rank1.txt under 2.17.4 (6 s); no CI line.

Own recomputations: Σ φ(m), 15000 < m ≤ 30000 = 205,177,458; degenerate coprime pairs (by Pythagorean parameters)
4,579 (2,678 square, 1,901 twice); non-degenerate 205,172,879 = H2's number. Primes p ≡ 1 (mod 4) in [17, 1153]: 91
(93 with 5, 13). Least p with (p − 511)² > 66,564 p: 67,583 (not prime); genus of a (ℤ/2)⁷-cover of ℙ¹ branched at 8
points: 129; 8·2⁶ = 512 points over the branch points: z8's threshold arithmetic holds.
H6 conductor's parity remark (unreviewed text in h6/NOTES): PARI 2.17.4 ellrootno: c3 (F0·L class) root number +1,
c5 (L·U) root number −1 (conductor 34,500,574,027,920, = z9). So with 2-parity (quoted) rank 3 for c3 and rank 0 for
c5 each force Ш[2^∞] of odd corank: the remark holds.
Reviews: reviews/review_z7.md, z8, z9 are byte-identical to conductor/reports/z7.md, z8.md, z9.md; experiments/referee_z7,
z8, z9 identical to conductor/r7/referee_z*. OUTPUTS totals: six tables, 168 rows with a Class column: 67 A + 5 mixed
"A for ..." = 72, B 54, C 27, D 15 — as stated; new table 22 rows (10 A, 10 B, 2 D); every file named exists.

(Correction to the clock above: the CI block ran 15:06-15:13 UTC, not "15:40"; the container's local clock is UTC−8.)

15:25 Own code (/tmp/referee7/own/, written from the prose before reading z8's scripts):
- `h3_own.py` (labels χ(l_j⁺) at points of ℙ¹(F_p) with every F_j a non-zero square, plus the 8 root discs, where the
  vanishing form's label is the class of the other factor): (5,1) at p = 1181 misses exactly vectors 5, 10; (17,1) at
  1181 misses 4, 7, 8, 11 — z8's counterexamples to H3's "good primes above 1153 remove nothing" reproduce.
- `p3_own.py`: at p ≡ 3 (mod 4) no (u:v) with all four forms non-zero squares for (2,1) and (9,1) at 7, (3,1) at 23,
  (5,1) at 47, (25,1) at 127, (14,1) and (151000467,1) at 151 (so the set is empty there); points exist for (2,1) at 11,
  (3,2) at 167. z8's claim 2 reproduces.
- `octic_own.py`: (11875080,6497849), the session's "open pair", at the octic prime 5462461 (≡ 5 mod 8; full numpy
  enumeration of ℙ¹(F_p), single-root discs, and the four collision points of linear factors, each with v_p(Res) = 1
  so that disc carries no allowed point): allowed set {5, 10}; both known tuples allowed (label 5); both extra tuples
  (2,1,1,98795483379599), (2,1,197590966759198,2) NOT allowed (label 9). z8's finding reproduces: not a phantom.
H3 tables and the corrections hold; H1's "one open pair" is withdrawn rightly.
- `gap_own.py` (numpy, same local test): the 14 ratios of H1's 7 true pairs at every good p ≡ 1 (mod 4) with
  1200 < p < 20000 (14,404 ratio-prime pairs, the range z8's phantom_mid.py skips and that rests only on G1's
  exact_local.py): known tuples allowed at all, extra tuples cut at none. Positive controls (`gap_control.py`): forged
  tuples (1109503,597671,6,1), (3,1,2,1) and one for (204724,97893) are cut (at 41, 29, 101); the vectorised test
  reproduces (5,1), (17,1) at 1181. So "locally soluble everywhere" for the 7 pairs holds; but the committed record
  for 1200 < p < 20000 is exact_local only, so "explicit points to 67,583" overstates what is committed.
15:30 Bundle provenance: unpacked /mnt/user-data/uploads/40b9ad47-roadmap-h.tar.bin in /tmp/referee7/bundle and compared
  with int7: the only changed bundle files are the banners of HANDOFF and REPORT_H1/H2/H3/H4/H6/H7, h2_control.py's
  fallback, and the docstrings of h3_weil.py, h3_sizes.py; the G3 block records (sd/dm/r0/stage 13000, times.txt) are
  byte-identical to the bundle's; CENSUS_TABLE_G3.tsv differs from the bundle's only in the (12007,1799) row.
MANIFEST: MANIFEST.sha256 is run 6's (1309 lines); verify_manifest.py on the clean archive: 288 files changed or
  unlisted -> MANIFEST FAILED; job `quick` runs it on every push. (Presumably the conductor's last step; listed as a
  MUST-FIX so it is not forgotten.)

## Checks by item (what holds)

1. NOTES / LEDGER / CHANGES / OUTPUTS against z7, z8, z9 and the outputs. Numbers checked and holding: H2 205,172,879 /
   30,434 / 2,807 / 229 / 2,131 (hunts) / 137 / (117702451,107508765) (h2_summary.txt, own recount of the totals); H3 91
   primes, 49,764 / 17,224 / 15,984 / 36 empty, last at 137 (h3_sizes.txt, z8's h3_counts.txt, own prime count),
   losses to 3853 (scan_losses.txt: 91 primes with a loss from 1153, the largest 3853), empties at p ≡ 7 (mod 8) to 151
   (p3mod4.txt, own p3_own.py); H1 7 true pairs and 5 removed pairs (full_descent.txt, badprimes_left.txt; own
   octic_own.py on the former open pair), 5 rank-0 + 2 pair-curve settlements (z8); H7 3/0, 3/0, 7/4, 5/2 and 5/0
   (h7_g4k_rank1.txt, h7_b13000_rank1.txt, rerun), E04 91/66/65/90 with 63 on E04 itself (46 + 16 + 1) and 2 on an
   isogenous curve ((320,289), (1015,961)), dm_evidence 91 E04 / 25 D1 / 6 E0L; (724,665), the one Magma-only tier-B
   pair left, is D1, so "the one left" stays true; rank0_ct 6/4 (rank0_ct.txt = z7); H6 bounds 6 to 9, A's bound 8
   (h6_richelot_627_364.txt); G3 916 = 540 + 376, block 305 = 270 + 23 + 9 + 1 + 2 (= 176 square + 129 twice, from
   the g3_table.txt differences), open/conditional/route lists. Tier D everywhere for roadmap H and rank0_ct: yes (all
   new LEDGER rows sit in the tier-D table; the tier-B and second-descent additions say tier D); nothing "proved";
   every new LEDGER row says "run 7, after release 0.6.0, not released" (tier-B row: "run 7, not released").
   Gaps: see MUST-FIX 4, 5 and SHOULD-FIX 5, 6, 9.
2. Withdrawn claims: corrected in the NOTES — H3's ">1153 removes nothing" and "96" (h3/NOTES), H1's "12 pairs / one
   open / every local test" (h1/NOTES), H2's "2,421", "P ≤ 1153", "96", "2.5 core-hours" (h2/NOTES), H7's
   "unconditional", "three", "(1015,713) all 6", "rank-1-route", "numerical L-value" (h7/NOTES). Banner present on
   HANDOFF and REPORT_H1/H2/H3/H4/H6/H7; missing on REPORT_H5, REPORT_H8 (MUST-FIX 3). Still readable as standing:
   docstrings of h1_partners.py, h1_scan.py, h2_design.py and line 2 of the committed h3_weil.txt (SHOULD-FIX 4); older
   text made stale (MUST-FIX 4, SHOULD-FIX 1-3). REPORT_G3.md is run 6's (banner and paths kept): z9's regression fixed.
3. OUTPUTS totals hold (168 = 72 + 54 + 27 + 15; new table 22 rows). CHANGES: the seventh-run entry sits under
   "Unreleased, after 0.6.0" before the sixth, same shape (bold run header with provenance, italic items, CI item).
4. CI: holds with one harness gap (MUST-FIX 2); timings and file checks above; MANIFEST (MUST-FIX 1).
5. rank0_ct.py: header's PARI sentence is z7's replacement verbatim in substance; outputs equal z7's (6, 4; controls
   (861,620) 4/2 → 0, (440,279) 4/0 → 2); --negative as z7 describes. z7's nit (torsion images in the kernel not
   asserted) not addressed (NIT).

## Findings, with exact replacement text

### MUST-FIX

M1. `MANIFEST.sha256` is run 6's (1309 lines): `python3 scripts/verify_manifest.py` on `git archive run7` reports 288
changed or unlisted files ("MANIFEST FAILED"); job `quick` runs it on every push. Fix: after every other edit,
`python3 scripts/verify_manifest.py --write` and commit it last (as run 6's "run 6: manifest").

M2. `.github/workflows/ci.yml`, the `hm` line: it checks only the last line, and `h3_tables.py` prints its JSON check
("# table matches h3_tables.json" / "# TABLE DIFFERS") before the last line. Tested: with one entry of h3_tables.json
changed (p = 41), `hm $H/h3 h3_tables.py` returns 0, while OUTPUTS lists h3_tables.json as class A in `descents`.
Replace
`hm() { (cd "$1" && shift && python3 "$@" 2>/dev/null | grep -v 'PARI stack' | tail -1 | grep -q '^# matches'); }`
by
`hm() { (cd "$1" && shift && python3 "$@" 2>/dev/null | grep -v 'PARI stack' > /tmp/hm.out && ! grep -q 'DIFFERS' /tmp/hm.out && tail -1 /tmp/hm.out | grep -q '^# matches'); }`
(tested: every one of the 24 hm calls passes with it on a clean archive, 4 min; the forged json fails).

M3. `experiments/roadmap-h/h5/REPORT_H5.md` and `h8/REPORT_H8.md` have no banner (the brief requires it on every
session report). h7/NOTES.md states the L-value's status correctly, but nothing points REPORT_H8's reader there, and
REPORT_H8 also says "numerical L(E,1)" (since run 6 an enclosure conditional on PARI's `ellap`), "resting on Vélu,
Mazur, Silverman and Cremona–Prickett–Siksek" (omits the Cassels–Tate pairing) and "would not have fitted anyway"
(run 6's streamed series of `census-wide/lvalue/` did fit). Insert after the title of each (line 2):
`> **Read \`../NOTES.md\` first.** This is the session's report, committed unchanged apart from this line; where they differ, \`../NOTES.md\` (after review) holds.`
and in `experiments/roadmap-h/NOTES.md` replace the H8 row's text
"not run; its purpose is met through (6557,324)'s partner (`h7/NOTES.md`)" by
"not run; its purpose is met through (6557,324)'s partner (`h7/NOTES.md`). In `REPORT_H8.md` read \"the L-value (conditional on PARI's `ellap`, `../census-wide/lvalue/`)\" for \"the numerical L(E,1)\", add the Cassels–Tate pairing as computed by `rank_second.py` to what H7 rests on, and note that run 6 already summed the series for this curve (1.21·10⁹ terms, `../census-wide/lvalue/`), so only the PARI-free `aₚ` of H8 is missing"

M4. LEDGER contradicts itself: the row "Random ratios to `m = 10⁹`" says of (1109503,597671) "whether a rational point
has one of these classes is open", while the run-7 row "Phantoms" (and h1/NOTES, z8) settle it twice. In LEDGER.md
line 104 replace
"whether a rational point has one of these classes is open (run 6, after release 0.6.0, not released;"
by
"whether a rational point has one of these classes was left open in run 6; since run 7 (row \"Phantoms\" below, `experiments/roadmap-h/h1`, tier D, not released) none does: `y² = w·f0·f8` has rank 0 and the pair curves `E_0L`, `E_0U` rank exactly 3 (run 6, after release 0.6.0, not released;"

M5. LEDGER row "Rank exactly 1 for F4's quotient `y² = g4·K` …" states "every solution of these ratios is degenerate"
without what it rests on (z7, section 8 item 4 and its replacement text, required it). In LEDGER.md line 108 replace
"(run 7, after release 0.6.0, not released; a separate AI session; reviewed by a separate AI agent: holds, with its own descent and pairing)"
by
"(computed, exact; resting on Vélu's formula, isogeny invariance, Mazur, the Cassels–Tate pairing in its Brauer–Manin form (Cassels 1998; Skorobogatov 6.2) as computed by `rank_second.py` (`reviews/review_t2.md`), Silverman and Cremona–Prickett–Siksek Theorems 2 and 15 (quoted, tier E); PARI only finds the root, points, conic points and factorisations, each checked; run 7, after release 0.6.0, not released; a separate AI session; reviewed by a separate AI agent: holds, with its own descent and pairing)"

### SHOULD-FIX

S1. Older text that run 7's H7 makes stale ((7735,4968) "open", the three ratios "conditional on rank 1" or "on a
numerical L-value"). Append, keeping the old text as history:
- LEDGER.md line 86, after "(6557,324), (6881,6233) on a numerical `L(E,1)`**": " (since 27 November 2026, run 7, after release 0.6.0, not released: (7735,4968) and (6881,6233) excluded by F4's Demjanenko–Manin chain with rank 1 computed exactly, (6557,324) through its partner; row \"Rank exactly 1 for F4's quotient\" below)".
- LEDGER.md line 103: "except (7735,4968) (open in the table; conditional on rank 1 by the `g4·K` route)" → "except (7735,4968) (open in the table; conditional on rank 1 by the `g4·K` route, a rank computed exactly since run 7, not released: row \"Rank exactly 1 for F4's quotient\" below)".
- README.md line 115: "and two on a numerical `L`-value." → "and two on a numerical `L`-value (since 27 November 2026, after release 0.6.0 and not released, (7735,4968) is excluded too and those two no longer need the `L`-value: F4's Demjanenko–Manin chain with rank 1 computed exactly, `experiments/roadmap-h/h7`, tier D)."
- experiments/README.md line 90: "descent gives 2, 3, 5. Conditional." → "descent gives 2, 3, 5. Conditional. Since run 7 (`roadmap-h/h7/`): rank 1 exactly for (7735,4968) and (6881,6233) by the complete 2-descent with the Cassels–Tate pairing on the 2-isogenous curve with full 2-torsion (review z7: holds; tier D); (6557,324) through its partner."
- experiments/census-wide/g4k/NOTES.md, append a paragraph: "**Run 7.** `../../roadmap-h/h7/h7_g4k_rank1.py` gives rank C(ℚ) = 1 exactly for (7735,4968) and (6881,6233): the 2-isogenous curve of C has full rational 2-torsion, and the complete 2-descent with the Cassels–Tate pairing of `scripts/g1_full/rank_second.py` gives dim Sel² = 3 with a pairing of rank 0, and 7 with a pairing of rank 4 (review z7, `reviews/review_z7.md`, with its own descent and pairing: holds; tier D, not released). With the chain here every solution of these two is degenerate; (6557,324), whose own curve gets only 3, follows through its partner (6881,6233)."
- experiments/census-wide/table/NOTES.md lines 16-17: "(`../g4k/` excludes it only conditionally on a rank being 1)" → "(`../g4k/` excludes it conditionally on a rank being 1; since run 7 that rank is computed exactly, `../../roadmap-h/h7/`; the table's labels are unchanged)".

S2. G3 text that still describes 611 ratios as the table:
- experiments/README.md line 17: "the census of the two degenerate types for `11000 < m ≤ 13000` (4 open, 3 conditional; `g3/`)" → "the census of the two degenerate types for `11000 < m ≤ 13000` (4 open, 3 conditional; `g3/`; in run 7 extended to `14000`: 916 ratios, 6 open, 3 conditional)".
- experiments/README.md line 141 (G3 bullet), after "One row's label regenerated from the corrected table.": " Run 7 added the block `13000 < m ≤ 14000` (305 ratios; open (13348,3795), (13899,1180); `g3/NOTES.md`, last section), so the table now has 916 rows."
- experiments/roadmap-g/g3/NOTES.md line 1: "# G3: census along the degenerate-type family, 11000 < m ≤ 13000 (Roadmap G, run 6)" → "# G3: census along the degenerate-type family, 11000 < m ≤ 14000 (Roadmap G, run 6; the block 13000 < m ≤ 14000 in run 7)"; at the start of "## What it claims" add "(Run 6's range; for the block added in run 7 see the last section.)"; line 47 "later blocks are not here." → "blocks from 14000 on are not here (block 14000 was started and stopped)."

S3. G1 text that still calls (1109503,597671) open (the LEDGER row is M4):
- experiments/roadmap-g/g1/NOTES.md line 25-26: after "Whether a rational point has one of these classes is **open**." insert " (Run 7: none does — `y² = w·f0·f8` has rank 0 and the pair curves `E_0L`, `E_0U` rank exactly 3; `../../roadmap-h/h1/NOTES.md`, review z8.)"
- experiments/README.md line 136: "whether a rational point has one of these classes is open." → "whether a rational point has one of these classes is open (run 7, `roadmap-h/h1/`: none does)."
- (optional) REPORT_G1.md's banner says the same; the NOTES fix suffices.

S4. Withdrawn claims still standing, unannotated, in committed bundle files: h1_partners.py docstring ("Since H3, good
primes to 1153 are all that can cut, so 'good primes to 20000' is the full local test"), h1_scan.py and h1_partners.py
("the exact local test at every odd prime", i.e. G1's one(), which skips the primes of m² + n² and of the octic above
20,000), h2_design.py docstring ("H3 shows that the good primes that can cut the four-class set are exactly the p = 1
mod 4 with 17 <= p <= 1153"), and line 2 of the committed output h3_weil.txt ("Hasse-Weil threshold (genus 17): every
good prime p >= 1154 keeps all 16 vectors"). Without touching the bundle, add:
- to h1/NOTES.md, "Corrections" list: "- The docstrings of `h1_scan.py` and `h1_partners.py` ('the exact local test at every odd prime'; 'Since H3, good primes to 1153 are all that can cut, so good primes to 20000 is the full local test') are wrong in the same way: G1's `one()` skips the primes of `m² + n²` and of the octic above 20,000, and H3's bound is withdrawn."
- to h2/NOTES.md, "Withdrawn" bullet, after "H3's Hasse–Weil step is wrong (`../h3/NOTES.md`)": "; the docstring of `h2_design.py` says the same (\"exactly the p = 1 mod 4 with 17 <= p <= 1153\"): read it as \"the primes of H3's tables\""
- to h3/NOTES.md, after "the STATEMENT / PROOF in the docstring of `h3_weil.py`": ", and line 2 of `h3_weil.txt` (\"every good prime p >= 1154 keeps all 16 vectors\")"

S5. "Explicit points to 67,583" overstates what is committed: z8's own good-prime test (phantom_mid.py) covers
p < 1200 and 20000 < p < 67,583; 1200 < p < 20000 rests on G1's exact_local.py (the scan, full_descent.py), "tested,
not proved". (My own vectorised test, /tmp/referee7/own/gap_own.py, finds no extra tuple cut there either, so the claim
is true; it is the wording that is too strong.)
- LEDGER.md line 109: "(explicit points to 67,583; Hasse–Weil above, quoted)" → "(explicit points at 2 and at every odd bad prime; at the good primes `p ≡ 1 (mod 4)` below 67,583 the review's own local test, except `1200 < p < 20000`, where G1's `exact_local.py` (tested, not proved) decides; Hasse–Weil above, quoted)".
- h1/NOTES.md line 28: "at every good prime `p ≡ 1 (mod 4)` below 67,583 (its own enumeration)" → "at every good prime `p ≡ 1 (mod 4)` below 1200 and between 20,000 and 67,583 (its own enumeration; between 1200 and 20,000 G1's `exact_local.py`, as in the scan, and an unpublished check by the pre-publication referee agrees)".

S6. "Computed (exact)" without its quoted inputs:
- h3/NOTES.md line 16: "Computed (exact); tier D." → "Computed (exact arithmetic; method A is G1's `exact_local.py`, tested in run 6, not proved; method B is an enumeration of ℙ¹(F_p), lifted by Hensel's lemma, quoted); tier D."
- h1/NOTES.md, end of "What it claims": add "- Computed (exact), resting on the soundness of the complete 2-descent and of the Cassels–Tate pairing as `rank_second.py` computes them (Cassels 1998; Skorobogatov 6.2), Vélu's formula and isogeny invariance of the rank, the injectivity of torsion into `E(F_p)` at good odd `p` (quoted, tier E); the local tests on G1's `exact_local.py` (tested, not proved) and the repository's one-sided test at 2. Tier D."
- LEDGER.md line 100: "so that rank no longer rests on PARI or Magma" → "so that rank no longer rests on PARI or Magma (resting on Vélu's formula, isogeny invariance and the Cassels–Tate pairing as computed by `rank_second.py`, quoted)".

S7. (13348,3795)'s rank 1 is a LEDGER claim with no CI line (`h7_b13000_rank1.txt` is class B though it takes 6 s; I
reran it under PARI 2.17.4: identical). Add after the `hm $H/h7 …` line of ci.yml:
`(cd $H/h7 && python3 h7_g4k_rank1.py 13348,3795 13899,1180 2>/dev/null | grep -v 'PARI stack' | diff - h7_b13000_rank1.txt)`
and in OUTPUTS move `h7_b13000_rank1.txt` from the class-B h7 row to the class-A h7 row (command
`python3 h7_g4k_rank1.py 13348,3795 13899,1180`, compared with `diff`); the counts line then reads 168 rows, 72 A, 54 B
(a file moves within rows; no row count changes).

S8. experiments/roadmap-h/NOTES.md lines 21-22: "One change at integration: `h2/h2_control.py` falls back to
`scripts/cover_descent/structure/class_selmer_15000_list.txt`." → "Changes at integration (checked against the
bundle): `h2/h2_control.py` falls back to `scripts/cover_descent/structure/class_selmer_15000_list.txt`; a note heads
the docstrings of `h3/h3_weil.py` and `h3/h3_sizes.py` (96 → 91); a 'Read NOTES.md first' line heads the hand-off and
each report. Nothing else in the bundle's `h1`–`h8` was changed." (Same correction for experiments/README.md's
"committed with the session's text unchanged and a `NOTES.md` in each folder": "… unchanged apart from those lines, and
a `NOTES.md` in each of `h1`–`h4`, `h6`, `h7`".)

S9. CHANGES.md line 9 (H3 item): "tables at the 91 primes `p ≡ 1 (mod 4)` up to 1153 hold" → "tables at the 91 primes
`p ≡ 1 (mod 4)` with `17 ≤ p ≤ 1153` hold" (there are 93 up to 1153); "and gives a complete local test to 67,583." →
"and argues (Hasse–Weil quoted) that the good primes `p ≡ 1 (mod 4)` below 67,583 and `p ≡ 7 (mod 8)` below 1154 make
the local test complete." CHANGES.md line 7 (phantoms), after "Corrected in `h1/NOTES.md`.": " Among the 7 is G1's
control (1109503,597671), whose extra classes run 6 left open."

S10. LEDGER tier E has no row for what run 7's tier-D rows newly lean on. Add:
"| The Cassels–Tate pairing on Sel² and its Brauer–Manin formula by local Hilbert symbols at points of the conics (used by `rank_second.py` and, in run 7, by `roadmap-h/h1`, `h7` and `rank0_ct.py`); the Hasse–Weil bound (used for the genus-129 threshold 67,583 in `roadmap-h/h3/NOTES.md`); the 2-parity theorem for elliptic curves over ℚ (Monsky; Dokchitser–Dokchitser), used in `roadmap-h/h6/NOTES.md` | Cassels 1998; Skorobogatov, *Torsors and rational points*, 6.2; Weil; Monsky 1996, Dokchitser–Dokchitser 2010 | Read by AI models; the pairing formula re-derived in `reviews/review_t2.md` and `review_z7.md`; not read by a person |"

### NIT

N1. OUTPUTS run-7 table: the h6 class-B row's `h6_ranks_*.txt`, `h6_richelot_*.txt` also match the class-A files
`h6_ranks_negative.txt`, `h6_richelot_negative.txt`, `h6_richelot_controls.txt` → write `h6_ranks_627_364.txt`,
`h6_ranks_991_263.txt`, `h6_richelot_627_364.txt`, `h6_richelot_991_263.txt`. The h7 class-B row says "(F4's scripts,
print-only …)" but `h7_b13000_rank1.txt` comes from `h7_g4k_rank1.py` (moot after S7). h7 class-A time "about 80 s
together": measured 111 s (h7_e04_rank 89 s).
N2. experiments/referee_z7/README.md ("`reports/z7.md` of run 7"), referee_z8/README.md ("`reviews/review_z8.md` if
committed (made as `conductor/reports/z8.md`)"), referee_z9/README.md ("`reviews/` (z9) or the conductor's
`reports/z9.md`") → "`reviews/review_z7.md`" etc. (the conductor's paths are not in the repository).
N3. SECOND_DESCENT_CHAIN.md row (609,200): "on the curve of the class with full 2-torsion" → "on a curve of the class
with full 2-torsion" (there are three; z7 got 6/4 on all three). Step-7 bullet: add the tag "the Cassels–Tate pairing
and its Brauer–Manin formula **[quoted: Cassels 1998; Skorobogatov 6.2]**". SECOND_DESCENT.md line 96 ("no 2-descent
gives rank 0") is literally true; a pointer "(the pairing does: `rank0_ct.py`, run 7)" would help.
N4. CHANGES seventh-run items: H7 "Reviewed: holds." → "Reviewed: holds (wording corrected in `h7/NOTES.md`)."; add
"Reviewed: holds with gaps." to the H6 item and "Reviewed: holds." to the G3 item; H2 item "205,172,879
non-degenerate ratios" → "205,172,879 coprime non-degenerate ratios".
N5. The banners say "committed unchanged" — they are themselves an added line: "committed unchanged apart from this
line".
N6. rank0_ct.py: z7's nit (bound() does not assert that the torsion images lie in the pairing's kernel) not done;
header "the pairing matrix symmetric with 1 on the diagonal" → "(values ±1) symmetric, with +1 on the diagonal".
N7. roadmap-h/NOTES.md H3 row and experiments/README.md roadmap-H paragraph: "a correct complete local test goes to
67,583" / "needs good primes ≡ 1 (mod 4) below 67,583" → add "and good primes ≡ 7 (mod 8) below 1154 (z8's argument,
Hasse–Weil quoted)".
N8. h6/NOTES.md: the conductor's parity remark is now checked: PARI 2.17.4 `ellrootno` gives +1 for the class of
y² = F0·L and −1 for y² = L·U (conductor 34,500,574,027,920); you may replace "(this remark is the conductor's, not the
review's)" by "(this remark is the conductor's; the root numbers +1 and −1 were checked with PARI's `ellrootno` by the
pre-publication referee)".
N9. LEDGER tier-B row: "since 27 November 2026, run 7, not released" → "since 27 November 2026, run 7, after release
0.6.0, not released".
N10. h1/NOTES.md: "A proof on either side settles a pair" → "An exclusion on either side settles a pair"; and REPORT_H1's
"Phantoms … occur only at degenerate ratios" is a statement about a scan of degenerate ratios only: add to the
corrections "the scan sampled degenerate ratios only; for the others see H2 and G1 (no survivor)".
N11. COMPLETE.md §7 test 6 ("this is not attempted") and test 7 could point to `experiments/roadmap-h/h3` (the
review's argument) and `h2` (no survivor to 30,000), tier D, not released. Optional.
N12. HANDOFF "Ideas tried" bullets on Hasse–Weil ("It worked, and it bounds H2's P") and the CT kernel test ("kills
the extras at … 550548247", a pair removed by a local test) are superseded only implicitly; one line in roadmap-h/NOTES.md
would make it explicit.

S11 (new; a consequence nobody stated). (12703,2767), open in G3's table, is the partner of (7735,4968)
((7735 + 4968, 7735 − 4968) = (12703, 2767), both odd, not halved; g3/NOTES.md itself says so). H7's step with F4's
chain makes every solution of (7735,4968) degenerate, and the partner map carries solutions and degenerate solutions
both ways (z7's partner_cells.py checks the nine cells symbolically for all m, n, halved and unhalved; in CI). So
(12703,2767) is excluded through its partner exactly as (6557,324) is — the session looked only at (12703,2767)'s own
g4·K (PARI [2,4]) and the reviews did not notice. Either state it or say why not. Text if stated (tier D, same
status as (6557,324); not separately reviewed):
- LEDGER.md G3 row: "(12201,7360) and (13348,3795) are excluded separately by H7's g4·K step, `experiments/roadmap-h/h7`)" → "(12201,7360) and (13348,3795) are excluded separately by H7's g4·K step, and (12703,2767) through its partner (7735,4968), which that step excludes, `experiments/roadmap-h/h7`)".
- LEDGER.md H7 row: after "(6557,324) excluded through its partner (6881,6233)," add " and (12703,2767), open in G3's table, through its partner (7735,4968),".
- CHANGES.md line 11: "6 open in the table, two of them excluded by H7's step" → "6 open in the table, three of them excluded by H7's step ((12703,2767) through its partner (7735,4968))".
- roadmap-g/g3/NOTES.md line 56: "(12201,7360) and (13348,3795) are excluded separately by H7's g4·K step (`../../roadmap-h/h7/NOTES.md`), not in this table" → "(12201,7360) and (13348,3795) are excluded separately by H7's g4·K step (`../../roadmap-h/h7/NOTES.md`), and (12703,2767) through its partner (7735,4968), which that step excludes; not in this table".
- roadmap-h/h7/NOTES.md, the (6557,324) bullet: add "The same move excludes (12703,2767), open in G3's table, through its partner (7735,4968) (the session's 'G3's other four open ratios … need another route' overlooks it)."

## Could not check
The long runs (H2's exhaustive blocks and hunts, H1's scan and pair-curve runs, G3's block, h3_weil.py's 320 s run,
h4_split.py, h6 class-B runs), rank_second.py's pairing code beyond the reviews' independent reimplementations and the
reruns, Richelot internals, anything needing Magma; CI jobs other than the run-7 lines of `descents` (plus g3_table,
g3_sample and verify_manifest run by hand).

## Final
(15:37 UTC.) **Verdict: ready after fixes.** No number in the new LEDGER rows, CHANGES lines, OUTPUTS rows or NOTES
was found wrong against z7, z8, z9 or the committed outputs; tier D throughout, nothing "proved", release labels
present; OUTPUTS totals (168 = 72 + 54 + 27 + 15) hold; the run-7 CI block passes from a clean `git archive run7` in
7 min 10 s with no file changed (YAML valid, inside the pipefail/f5 step, ~27 min margin to 150); rank0_ct.py's header
and outputs match z7. Own recomputations agree: H2's totals (205,172,879; 4,579 degenerate), the 91 primes, z8's
67,583/genus 129 arithmetic, z8's H3 counterexamples ((5,1), (17,1) at 1181; empty sets at p ≡ 7 mod 8), the former
"open pair" cut at 5462461, the 7 true pairs' extras surviving at every good p ≡ 1 (mod 4) in (1200, 20000), the H6
root numbers (+1, −1).
MUST-FIX: M1 MANIFEST.sha256 not regenerated (quick fails); M2 `hm` ignores "# TABLE DIFFERS" (forged
h3_tables.json passes; stricter hm given and tested); M3 no banner on REPORT_H5/H8, H8's stale wording; M4 LEDGER's G1
row still calls (1109503,597671)'s classes open; M5 LEDGER's H7 row omits what it rests on. SHOULD-FIX S1–S11 (stale
(7735,4968)/G3/G1 text; unannotated withdrawn docstrings and h3_weil.txt; "explicit points to 67,583"; quoted inputs;
a CI line for (13348,3795); "one change at integration"; CHANGES H3 wording; tier E row; and S11, new: (12703,2767) is
excluded through its partner (7735,4968), which nobody stated). NITs N1–N12.
