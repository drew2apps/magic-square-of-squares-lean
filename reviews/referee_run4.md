# Referee run 4 (hostile pre-publish) — diff 0a52347..HEAD in /home/claude/conductor/int4

Started 01:49 UTC. Findings written as found; verdict at top updated at end.

## Verdict: PUBLISH AFTER FIXES (final, 02:00 UTC)

The Lean and computed results hold as far as checked here (UniformRank.lean says exactly what it should; counts 147,
13/37/58 and 53 match; tables and paths fine; CI hygiene fine). What blocks publication is prose: the headline Magma
count is stale in four places, Theorem 2 drops `m ≠ ±n` in four places, UNIFORM_RANK.md still says Lemma D's descent
is not in Lean, and CHANGES.md carries a placeholder over an unresolved census contradiction.

FIX-FIRST: F1 (+addendum), F2, F3, F7.   FIX: F4, F5, F6, F8, F9.   NIT: N1-N9.

## Findings

### F1 (FIX) README.md:102 and G1_FULL_SYSTEM.md:24 still say six counted ratios rest on Magma alone *including (72,65), (137,7)*
README.md:118 (new) says Θ₄ removes (72,65)/(137,7), "four of those six still rest on Magma alone", total six (was eight).
README.md:102 (unchanged) still lists "(760,39), (799,721), (945,248), (72,65), (137,7), (920,129)" as resting on Magma alone,
and G1_FULL_SYSTEM.md:24 the same ("six counted ratios do: ... (72,65) and its partner (137,7) ..."). A reader sees "six" twice
with two different lists. Replace in README.md:102
  "so six counted ratios rest on Magma alone for it ((760,39), (799,721), (945,248), (72,65), (137,7), (920,129));"
with
  "so six counted ratios rested on Magma alone for it ((760,39), (799,721), (945,248), (72,65), (137,7), (920,129)); since 26 November 2026 the bound for the curve of (72,65), (137,7) also comes from Fisher's `Θ₄` (exact, tier C; below), so four do;"
and in G1_FULL_SYSTEM.md:24 append after "and (920,129)." : " (Since 26 November 2026 the bound for the curve of (72,65), (137,7) also comes from Fisher's `Θ₄`, `rank_fisher_theta4.py`, tier C; so four counted ratios of this row rest on Magma alone; README.md.)"
LEDGER.md:53 is already correct.

### F2 (FIX-FIRST) Theorem 2 prose drops the hypothesis `m² ≠ n²`
`UniformRank.theorem2` assumes `IsCoprime m n`, `m² + n² = 2k²` AND `m² ≠ n²`. (1,1) is a coprime pair with 1 + 1 = 2·1², so
"every coprime pair with m² + n² twice a square" is false as written (curve singular there). Affected:
- README.md:171 "Theorem 2 at every coprime pair with `m² + n²` twice a square" -> "Theorem 2 at every coprime pair with `m² + n²` twice a square and `m ≠ ±n`"
- UNIFORM_RANK.md:3 "Theorem 2 with `m² + n²` twice a square" -> "Theorem 2 with `m² + n²` twice a square and `m ≠ ±n`"
- LEDGER.md:37 "at every coprime pair with `m² + n²` twice a square (`theorem2`)" -> "at every coprime pair with `m² + n²` twice a square and `m ≠ ±n` (`theorem2`)"
- CHANGES.md:4 "Theorem 2 (rank `≥ 2` at every coprime pair with `m² + n²` twice a square)" -> "... twice a square and `m ≠ ±n`)"
(UNIFORM_RANK.md:87 item (3) already says `m² ≠ n²`; README/LEDGER row 35 of run 3 "every coprime twice-square pair" has the same gap — Theorem2Family's own signature should be checked; see F5.)

### F3 (FIX-FIRST) UNIFORM_RANK.md still says, twice, that the descent over Q(√2) / Lemma D are not in Lean
- UNIFORM_RANK.md:3 (status line, same paragraph that announces the 26 Nov result): "the descents, the torsion bounds, the local tests and the identities over `ℚ(i)` are not (section 7)."
  -> "the descents over `ℚ` and `ℚ(i)`, the torsion bounds over `ℚ`, the local tests and the identities over `ℚ(i)` are not (section 7; the descent and the torsion over `ℚ(√2)` of Lemma D are, since 26 November 2026, below)."
- UNIFORM_RANK.md:85 (section 7, *Lean.*): "Lemma D, the descent over `ℚ(√2)` and the local tests are not in Lean (Lemma B is, since 25 November 2026: `lean_rank/LeanRank/LemmaB.lean`)."
  -> "The local tests are not in Lean; Lemma B is (since 25 November 2026: `lean_rank/LeanRank/LemmaB.lean`), and Lemma D with the descent and the torsion over `ℚ(√2)` are since 26 November 2026 (next paragraph)."
- Same paragraph :85, run-3 heading "**Theorems 1 and 2 for every ratio in Lean, with Lemma D as the only hypothesis**" -> "**Theorems 1 and 2 at every coprime pair in Lean, with Lemma D as the only hypothesis (discharged since 26 November 2026, below)**".
This is exactly the stale-contradiction class the brief names.

### F4 (FIX) LEDGER.md stale rows
- LEDGER.md:77 (F0·F8·GL row, "Where"): "the 2-descent step itself in Lean for the 254 curves of `m ≤ 1000` (`TwoDescentCertificates.lean`, tier A), not Lemmas B and D"
  -> "the 2-descent step itself in Lean for the 254 curves of `m ≤ 1000` (`TwoDescentCertificates.lean`, tier A); since 26 November 2026 Lemmas B and D and the statements at every coprime pair as well (`lean_rank/LeanRank/UniformRank.lean`; the scaling to non-coprime `m, n` on paper)".
- LEDGER.md:35 (run-3 row) "`LemmaD` (D2 has only `(±1, ±4)`) is a hypothesis, proved in tier C" -> "... is a hypothesis of these theorems, proved in tier C and, since 26 November 2026, in Lean (`UniformRank.lemmaD`, two rows below)". (NIT-level but it is the row a reader lands on first.)

### F5 (FIX) CHANGES.md:8, (627,364) bullet (w5, unreviewed) overstates `kctp_isog.txt`
(a) "still rank 0 or 2 over `ℚ(i)`" drops the finiteness of Ш that h4 restored everywhere else (SECOND_DESCENT.md §6,
KCTP_627.md "Unconditionally the result is only rank ≤ 2") -> "still rank ≤ 2 over `ℚ(i)` (0 or 2 if Ш is finite)".
(b) "so no 2-descent with its pairing in the isogeny class can prove rank 0": kctp_isog.txt:52-58 shows this only for E,
E/⟨T1⟩, E/⟨T2⟩, E3 and says for E3/⟨S1⟩, E3/⟨S2⟩ "it could only if n3 = 2" (n3 = n − 1 ≥ 2, so n3 = 2 is possible).
Replace with: "so the 2-descent with its Cassels–Tate pairing cannot prove rank 0 on four of the six curves of the
isogeny class (`E`, `E/⟨T1⟩`, `E/⟨T2⟩`, `E/⟨T3⟩`), and on the other two only if `n = 3` (`kctp_isog.txt`, item 5; the
step rests on Cassels' isogeny invariance of the BSD quotient, quoted)".
KCTP_627.md:7 (run 3 text, now published) "so the rank over `K` is 0 or 2" -> "so the rank over `K` is at most 2, and 0 or 2 if Ш is finite (parity)"; also NIT: it names internal branch `sess_t5` twice (":7", ":11") — say "an earlier AI slice (25 November 2026)".

### F6 (FIX) "the end-to-end statement [was] reviewed by separate AI agents" (UNIFORM_RANK.md:3) is stronger than the record
review_w2.md:109-118,153: UniformRank.lean was untracked then; the w2 reviewer checked the composition at type level and
read it, but did not compile it ("not compiled by me"); no reviewer built it with w1's files. Replace
"Each of the two halves and the end-to-end statement were reviewed by separate AI agents; no person has read them."
with "Each of the two halves was reviewed by a separate AI agent; the composition in `UniformRank.lean` was read by
the torsion review and by the pre-publication referee and built by CI's `rank` job; no person has read them."
(This referee read UniformRank.lean: the three statements are verbatim `rank_one_family`, `P_infinite_order`,
`rank_family_all` with `hD := lemmaD`; no hypothesis added or dropped; `LemmaD` is Lemma D itself, Theorem2Family.lean:35.)

### F7 (FIX-FIRST) CHANGES.md:9 leaves a placeholder and exposes an unresolved census contradiction
CHANGES.md:9 ends "...as that document says it does (see the entry below if it was settled in this run)." No entry below
settles it, and G1_FULL_SYSTEM.md:211 still says "`global.py` reproduces the closing curve for each closed pair" while
CHANGES now says (review_w6 FIX-FIRST 1) it does not for (112,15), (127,97), (700,529) — three COUNTED census closures.
Publishing both sentences is a self-contradiction on a census claim. Minimum fix before publishing:
- CHANGES.md:9: replace "(see the entry below if it was settled in this run)" with "; not settled in this run: the three
  stay counted as before, and `G1_FULL_SYSTEM.md` now says that their closing curve is not reproduced by `global.py`".
- G1_FULL_SYSTEM.md:211: "`global.py` reproduces the closing curve for each closed pair." -> "`global.py` reproduces the
  closing curve for each closed pair except (112,15), (127,97) and (700,529), whose `L(E,1) ≠ 0` closures (above) it
  does not find again (found by a review on 26 November 2026; open)."
- Add the item to the run-4 "Not done" bullet of CHANGES.md and to LEDGER's census row caveats.
(Better: settle it — recompute those three curves and certify `L(E,1)` with `lball_exact53`'s machinery — but that is new work.)
lball_exact53 wording itself (CHANGES:9, LEDGER new tier-D row, SECOND_DESCENT.md:112, OUTPUTS rows) correctly says
"closures of the experiment of SECOND_DESCENT.md, not of the census"; I found no sentence calling the 53 counted.

### F8 (FIX) SELMER_SQUARE.md:13 and CHANGES.md:13 (h4): "14 of 120" undercounts, and PARI's bound is not "the 2-Selmer bound"
selmer_square_pari.txt: PARI's upper bound is below `s` on 100 of the 120 curves (s=3:[1,1] 14; s=4:[2,2] 42; s=5: 25;
s=6: 16; s=7: 3). "14" is only the s = 3 row. And `ellrank`'s upper bound may use more than Sel².
- SELMER_SQUARE.md:13 "and can be smaller: on 14 of the 120 random curves of `selmer_square_pari.txt`, `s = 3` while PARI's rank bounds are `[1, 1]` (guide only)."
  -> "and can be smaller: on 100 of the 120 random curves of `selmer_square_pari.txt` PARI's upper rank bound (`ellrank`, a guide only; it is not exactly the 2-Selmer bound) is below `s`, for example `s = 3` with PARI's `[1, 1]` on 14."
- CHANGES.md:13 "which can exceed the 2-Selmer bound (14 of 120 curves of `selmer_square_pari.txt`)" -> "which can exceed the full 2-descent bound (PARI's upper bound, a guide, is below `s` on 100 of the 120 curves of `selmer_square_pari.txt`)".

### F9 (FIX) REPRODUCE.md CI table: times and contents stale after run 4
- REPRODUCE.md:55 `quick` "a few minutes" — the job limit is now 45 minutes (README.md:198 says 13 at `78659d1`, before
  run 4 added selmer_square --max 300, review_t8, CensusClosure regeneration, e1/torsion checks, lball_exact53).
  -> time cell "about 15 to 20 minutes (13 at `78659d1`; limit 45)", and add to "What it reruns": "; since run 4 the
  exact checks of `E1(ℚ(√2))` (`e1_rank0_check.py`, `torsion_e1_check.py`), a regeneration of `CensusClosure.lean`
  with `census_closure_check.py`, `selmer_square.py --max 300`, `review_t8_rank1.py`, and `lball_exact53.py` on three curves".
- REPRODUCE.md:56 `descents` "up to about an hour each": descents was 51 min at `78659d1` and run 4 adds ~25 (lball_exact53
  51 curves) + 4 (aptest) + ~4 (Θ₄, Fisher) minutes; limit raised to 150. -> "up to about an hour each (`descents` about
  90 minutes since run 4; limit 150)".
- NIT: table has no `scan-15000` row; REPRODUCE.md:61 list of jobs ending with `git diff --exit-code` omits `square-section`.

### F1 addendum (FIX-FIRST, upgraded): README.md:100 still says "eight"
README.md:100: "for the whole census, eight counted ratios rest on Magma alone for an upper rank bound and two on its
two-cover descent, see below" — contradicts README.md:118 ("**six** ... (eight until 26 November 2026)") in the same file.
-> "for the whole census, six counted ratios rest on Magma alone for an upper rank bound (eight until 26 November 2026) and two on its two-cover descent, see below".
Also G1_FULL_SYSTEM.md:343 "for `E⁺` of (760,39), (945,248), (72,65) and (920,129) PARI's bound is 3 and Magma's is 1, so those four rest on Magma for it" -> append "(since 26 November 2026 the curve of (72,65) also has bound 1 from Fisher's `Θ₄`, `rank_fisher_theta4.py`, tier C, so three do)".
With README.md:102 and G1_FULL_SYSTEM.md:24 (F1) that is four stale copies of the headline count; the headline count
"six" appears with three different meanings. Upgraded to FIX-FIRST.

## Checked and found consistent (no action)
- 147 `#print axioms` lines: `grep -c '^#print axioms' lean_rank/LeanRank/*.lean` sums to 147 (UniformRank 5, CensusClosure 26, E1Descent 13, TorsionE1 13, ...); all at column 0; lakefile globs `LeanRank.+`, so UniformRank is built.
- UniformRank.lean: `theorem1`, `theorem1_P`, `theorem2` are verbatim `rank_one_family`, `P_infinite_order`, `rank_family_all` (Theorem1Family.lean:191, Theorem1Order.lean:103, LemmaDEquiv.lean:183) with `hD := lemmaD`; same binders (`IsCoprime`, `m ≠ 0`, `n ≠ 0`, `m² ≠ n²`; resp. `m² + n² = 2k²`, `m² ≠ n²`). `LemmaD` (Theorem2Family.lean:35) is Lemma D itself. `lemmaD := lemmaD_of_torsion torsion_abscissa`; no `sorry`/`axiom`/`native_decide` in LeanRank/*.lean or the Checks files.
- CensusClosure: six `ratio_ne_*` with the eight line equations and `b² ≠ d²`, `b² ≠ f²`; reduction at 11; partners 12:5→17:7, 24:7→31:17. Prose "with distinct entries" is implied-by, so OK.
- Square section: square_section.tsv statuses 13 point / 37 sha / 58 open; txt tally "(point: 3 in the box, 10 from hints)"; the by-s table rows sum to 13/37/58 and the column sums (8/26/44/22/7/1) are unchanged from 391e693.
- Markdown tables: every row of LEDGER, OUTPUTS, REPRODUCE, README, UNIFORM_RANK, SQUARE_SECTION, SELMER_SQUARE, KCTP_627, SECOND_DESCENT has the header's cell count (code spans masked).
- Paths: every backticked file name added in the diffed .md exists in the tree, except Stoll's `WeierstrassFormalGroup/Reduction.lean` and Kepler71's `bridge/qc40/...` (external, fine).
- CI: no new top-level `a && b`, no `head`/`grep -q` on a live pipe (ccc.txt, r2neg.txt read from files), every new step under `set -e -o pipefail`; temp files `cc.lean ccgen.txt ccc.txt` and `$T` removed; new axiom regex reads wrapped lists (`[^\]]*` spans newlines) and primed names (`'.*?' ` extends past an inner `'`); a list wrapped before `[` would fail closed, not open. Forge pins compare exact line lists incl. the duplicate `30`. Timeouts: quick 45, descents 150 (est. ~85 min), square-section 30 (5–7 min) are adequate. Writes in the merged Python (rank_fisher*, rank_second, kctp*, theta4*, selmer_square dump modes, lball_exact53) are all under `--write` or a user-given path.
- lball_exact53 wording (CHANGES:9, LEDGER tier-D row, SECOND_DESCENT.md:112, OUTPUTS) says "closures of the experiment, not of the census" everywhere; nothing calls the 53 counted.

## NIT
- N1 ci.yml (descents, Θ₄/Fisher loops): `2>/dev/null` discards tracebacks; write stderr to `$RUNNER_TEMP/$s.err` and `cat` it on failure.
- N2 ci.yml rank_fisher5_local `--negative | grep -v '^rank_fisher5_local.py: local points'`: the committed file is not what the command prints; drop that header line from the script's `--negative` mode and remove the filter.
- N3 rank_fisher_theta4.py:531 (`qfsolve`), `_k4.py` (`bnfinit`, `bnfisprincipal`): diffed outputs depend on PARI's choices; pinned by cypari2==2.2.9 — say "A (with cypari2 2.2.9 / PARI 2.17.4)" in OUTPUTS.md.
- N4 ci.yml rank: `E1RankChecks`, `CensusClosureChecks`, `RankChecks` pass even with a `sorry` warning; capture output and fail on `declaration uses 'sorry'`.
- N5 LEDGER.md:38 "generated by `scripts/gen_census_closure.py`" -> "`lean_rank/scripts/gen_census_closure.py`".
- N6 CHANGES.md:4 "The scaling to coprime pairs" -> "The scaling from any `m ≠ ±n` to coprime pairs".
- N7 README.md:171 "six census ratios excluded end to end" -> add "(for squares with distinct entries; the Lean hypotheses are `b² ≠ d²`, `b² ≠ f²`)" (review_w3 FIX 3).
- N8 SECOND_DESCENT.md:112 "every enclosure lies inside the one of `lball.py`": the program only checks that they meet (lball_exact53.py:447); containment does hold at the printed digits for all 53 (checked here) — say "meets, and at the printed digits lies inside".
- N9 OUTPUTS.md run-4 table lists no class-D hint files (kctp_hints.json, rank_fisher_hints.json, rank_second_hints.json, rank_fisher_curves.json, rank_fisher5_hints.json, rank_fisher5_pushout.json, rank_fisher_theta4_{k4,k6,control}_hints.json).
