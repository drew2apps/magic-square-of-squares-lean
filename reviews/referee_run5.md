# Referee report, run 5 (RC 0.6.0), branch run5 HEAD 419c91a vs base 01e35dd

Written by an AI model (pre-publication referee); not read by a person.
Started 21:48 UTC.

## 1. Headline numbers (checked so far, 21:55 UTC)

- 177 `#print axioms` lines: `grep -c '^#print axioms' lean_rank/LeanRank/*.lean` sums to 177 at HEAD, 147 at 01e35dd
  (new: CensusClosureMSS 9, EminusModel 10, UniformRankAll 11). CI compares the set, not the number. HOLDS.
- ten forged-claim files, TwoDescentCertForge 14 pinned lines (8 + 6), UniformRankAllForge 9: ci.yml lines 120-124. HOLDS.
- 22 forged torsion certificates: torsion_exact.py `negative()` = 3 curves x 6 + 11a2 + 2 (empty, composite) + x^3-2x = 22;
  CI greps `negative tests: 22 forged certificates, 22 refused`. HOLDS.
- torsion 8/8/3: census_lvalue_exact_torsion.txt: Z/4 x8, Z/2xZ/4 x8, Z/2xZ/8 x3. HOLDS.
- census 268+122+18+6+2 = 416, 425-416 = 9; at 0.5.0 268+122+9+2 = 401, 24 open. HOLDS.
- "six counted ratios rest on Magma alone (eight before)": 3 curves (760,39),(945,248),(920,129) -> 4 ratios with
  partners, + (631,71),(351,280) of the 122 = 6; before, + (72,65),(137,7) = 8. HOLDS (README/LEDGER/G1 consistent).
- OUTPUTS totals: rows per table 80 + 14 + 12 + 11 = 117; pure A 42, B 32, C 26, D 12, plus 5 rows "A for ..., B for
  the rest", which the totals count as A (47). Before: 106 = 41/30/24/11 on the same rule. HOLDS, but the rule is
  unstated (NIT N1).
- Action pins: `git ls-remote` upstream: actions/checkout v5.1.0 = fbc6f3992d24... (also = v5), actions/setup-python
  v6.3.0 = ece7cb06caef... ; all 10 checkout and 7 setup-python uses carry the same hash and comment; no old hash left.
  HOLDS.

## 4. CI lines rerun from a clean clone (/tmp/ref5 = git clone of int5 at 419c91a), 21:52-22:00 UTC

quick (new lines): torsion_exact.py self-test (0.7 s), --negative "22 forged certificates, 22 refused" (0.5 s),
census_lvalue_exact.py --torsion (1.8 s) and --torsion-negative (1.4 s), rank_fisher_theta4_order.py and --negative
(<1 s each), kctp_isog_exact.py and --negative (Sage, 4 s, 1 s): all PASS (diff empty).
quick (lines whose committed outputs changed): lball_exact.py --all vs lball_exact.txt + _negative.txt (57 s),
lball_exact53.py --negative (11 s), gen_census_closure.py regeneration cmp CensusClosure.lean: PASS.
"PARI nowhere": with sys.modules['cypari2'] = None, --torsion, --torsion-negative, torsion_exact --negative,
rank_fisher_theta4_order.py and census_lvalue_exact.py --pairs '112,15' give byte-identical output. HOLDS at code
level (the CLOSURES table pair -> curve was found by the PARI-guided census_lvalue_identify.py; a finder only).
descents (Fable lines, run with the CI's f5 filter and pipefail): see the end of this section when done.

descents (Fable, 14 new lines, exactly as in ci.yml with `2>/dev/null | f5 | diff`): all 14 PASS, 4.3 min in all
(g4k_exact --negative 71 s the longest); `git status` of the clone clean afterwards. With bash -o pipefail (GitHub's
`shell: bash`), a crash after full output would still fail the line, so `2>/dev/null` (not x3's `2>&1`) is safe.

## 5. reviews/ (22:00 UTC)

No secrets: no tokens/keys/passwords/emails (other than none), no IPs; no personal data beyond the owner's name
(public in CITATION.cff) and machine paths. README.md lists exactly the 31 report files present (diff of names: empty)
and each one-line description matches the report's own title/scope (spot-read every header). Gap: see S6.

## 6. Fable F4 conditional everywhere: HOLDS
CHANGES.md:12 ("excluded **if** ... rank 1"), experiments/README.md:13, :81-83, g4k/NOTES.md, g4k_exact.py docstring
("given rank C(Q) = 1"), g4k_exact.txt ("IF rank C(Q) = 1" on every conclusion line). LEDGER.md:85, README.md:115,
CENSUS9000.md:8/23/52 keep (7735,4968) "open", which is consistent. Only dm_g4k.jsonl has a bare `"status": "closed"`
(a program field; its `gens` field says conditional) - acceptable.

## Findings (22:02 UTC)

### MUST-FIX
M1. CITATION.cff:24-30 (abstract), now `version: 0.6.0`, still gives the 0.5.0 census and "nobody has reviewed"
  ("268 ... 122 by rank 1; 9 more by an argument with four maps ...; 2 more ...; 24 are open"). Zenodo shows this.
  Replace the lines from "Computed and not formalized:" to "...last pair of ratios left by them." with:
  ```
  published search bounds). Computed and not formalized: a census of the ratios with m <= 1000
  (of 425 ratios left by congruences, 268 are closed by rank 0, 122 and 18 by the Demjanenko-Manin
  argument on curves of rank 1; 6 more by an argument with four maps derived here and 2 more by a
  two-cover descent run in Magma, both counted apart; 9 are open). A Lean sub-project on Lean v4.36.0
  with M. Stoll's Mordell-Weil library proves rank >= 1 on the two elliptic quotients of
  y^2 = F0 F8 GL for every ratio, and rank >= 2 when m^2 + n^2 is twice a square. The repository also
  contains experiments that are not counted in the census and that no person has reviewed (reports of
  AI reviews are in reviews/).
  ```
  (keep the first sentence part "One edge root exceeds 500 (far below" as is; only the tail changes).
M2. CHANGES.md: inside the 0.6.0 entry, 18 sub-entries say "nothing released" / "nothing is released"
  (lines 16, 25, 31, 38, 47, 76, 82, 102, 104, 111, 118, 129, 133, 137, 140, 145, 150, 159); a reader of the v0.6.0
  tag reads "not released" for things in it. Minimal fix, CHANGES.md:3, replace
  `Changes since 0.5.0 (2026-11-23), newest first;` with
  `Changes since 0.5.0 (2026-11-23), newest first; each entry below was written when its work was committed, so its "nothing (is) released" means "not released at that date": all of it is in 0.6.0;`
  Also CHANGES.md:102 "`PAPER.md` and `CITATION.cff` describe release 0.5.0 (24 open)" becomes false once M1 is
  applied: append " (`CITATION.cff` was updated for 0.6.0; `PAPER.md` was not)".

### SHOULD-FIX
S1. README.md:102 (end): "kernel-checked with no hypothesis since 26 November 2026, for coprime pairs:
  `lean_rank/LeanRank/UniformRank.lean`). Written by an AI model; nobody else has read the bound." contradicts
  README.md:172, LEDGER.md:37, UNIFORM_RANK.md:3. Replace with: "kernel-checked with no hypothesis since 26 November
  2026, for every rational pair: `lean_rank/LeanRank/UniformRank.lean` for coprime pairs, `UniformRankAll.lean` for all
  (run 5)). Written by an AI model; no person has read the bound."
S2. README.md:172: "The sub-project compiles `G1Sigma`, `G1SigmaGeneric`, `TwoDescentMap` and
  `TwoDescentCertificates` from this tree unchanged." -> "... `TwoDescentCertificates` and, since run 5, `Defs` from
  this tree unchanged." (lakefile.toml:27-28 has five roots; the CHANGES claim "the stale 'four modules' wording" fixed
  is not true here.)
S3. LEDGER.md:128 still a "four curves ... Magma alone" tail: "The four curves whose rank 1 rests on Magma alone: a
  deeper 2-power descent or a Cassels–Tate pairing (the exact second descent of `rank_indep.py` stops at 3)." ->
  "The three curves whose rank 1 rests on Magma alone (`E⁺` of (760,39), (945,248), (920,129); the curve of (72,65)
  has `Θ₄` since 26 November 2026): `Θ₄` as for (72,65) (it needs points on conics over `ℚ(√±799)`, `ℚ(√±23)` and
  larger fields), or a deeper 2-power descent (the exact second descent of `rank_indep.py` stops at 3)."
S4. LEDGER.md:33-35 still say `Eplus m (−n) = E⁻` is sympy-only, against row 37 ("no longer on paper").
  Row 33: "that the list of 127 is complete, that `Eplus m (−n)` is `E⁻` (sympy, and a cross-check on the 254
  curves), and Lemmas B and D (uniformity) are not Lean statements" -> "that the list of 127 is complete is not a Lean
  statement; `Eplus m (−n)` = `E⁻` was sympy (and a cross-check on the 254 curves) until run 5 (`EminusModel.lean`),
  and Lemmas B and D are in the rows below". Row 35: "`Eplus m (−n) = E⁻` by sympy" -> "`Eplus m (−n) = E⁻` by sympy
  then, in Lean since run 5 (`EminusModel.lean`)".
S5. "read by nobody (else)" where a committed AI review exists, against "reviewed by a separate AI agent" (the same
  issue as "not checked by anyone", fixed only in UNIFORM_RANK): LEDGER.md:12 "checked by nobody" -> "checked by no
  person (the 18 of 25 November 2026 by a separate AI review)"; LEDGER.md:13 "read by nobody" -> "read by no person";
  LEDGER.md:53 and G1_FULL_SYSTEM.md:24 "read by nobody else" -> "read by no person (one AI review,
  `reviews/review_s5.md`)"; SELMER_SQUARE.md:75 "have been read by nobody" -> "have been read by no person (one AI
  review, `reviews/review_t4.md`)"; SQUARE_SECTION.md:52 "has been read by nobody" -> "has been read by no person
  (one AI review, `reviews/review_t8.md`)".
S6. CHANGES.md:3 "The AI reviews behind the words "reviewed by a separate AI agent" are committed in `reviews/`."
  is too strong: LEDGER.md:32 / CHANGES.md:47 (three AI reviews of `ClassPadic.lean`) and LEDGER.md:69-70 (the
  review of `0b34613`, Theorems C, C′) have no report there (grep of reviews/ for ClassPadic, 0b34613: none; none in
  conductor/reports either). Replace with "The reports of the AI reviews of 25 and 26 November 2026 are committed in
  `reviews/` (a few earlier ones, among them those of `ClassPadic.lean` and of `0b34613`, were not saved)." and add
  the same sentence to reviews/README.md after "has no report here."
S7. REPRODUCE.md:55 (`quick` row) does not list the run-5 lines: append "; since run 5 the torsion step of the 19
  census `L`-value closures (`torsion_exact.py`, `census_lvalue_exact.py --torsion`), the order of `G` on `E⁺` of
  (72,65) (`rank_fisher_theta4_order.py`) and the (627,364) period ratio (`kctp_isog_exact.py`, Sage)". REPRODUCE.md:56
  (`descents`): after "the `Θ₄` and Fisher scripts" add "; run 5 adds about 5 minutes for `experiments/census-wide/g4k/`
  and `rank1/`".
S8. OUTPUTS.md:165 (run-4 table) says every `.json/.jsonl/.txt/.log` under `census-wide/`, `open-ratio/` is class B
  with CI job none; the run-5 table puts g4k/ and rank1/ outputs in A (descents), the hints in D, val_Q_run* in C.
  Add to that row's first cell: "(except the files of the run-5 table below)". CHANGES.md:13 "(the page lists 106
  files before this release)" -> "(the page lists 106 rows before this release)" (OUTPUTS.md:188 counts rows).
S9. CHANGES.md:12: "A second implementation of the Cassels–Tate pairing on `Sel²` over `ℚ(i)` agrees that it vanishes
  for (627,364) (`open-ratio/ctp/`)." sits under "reviewed: ..." but experiments/README.md:13 calls it unreviewed ->
  "... (`open-ratio/ctp/`; not reviewed beyond its wording and the real-place fix)".

### NIT
N1. OUTPUTS.md:187-188: say how mixed rows count: after "(a row may name several files)" add "; the five rows that
  are A for part of a file and B for the rest count as A".
N2. LEDGER.md:51: "...PARI's `elltors` agrees on all 19, a check; reviewed by a separate AI agent: holds); reviewed by
  a separate AI agent: holds with gaps)" -> "...a check; this step reviewed by a separate AI agent: holds); the
  `L`-value step reviewed by a separate AI agent: holds with gaps)".
N3. G1_FULL_SYSTEM.md:24: the cell still leads with "four of the curves rest on Magma alone for it, so six counted
  ratios do" and corrects it in parentheses; lead with the present: "three of the curves rest on Magma alone for it,
  so four counted ratios do: (760,39) and its partner (799,721), (945,248) and (920,129) (until 26 November 2026 also
  the curve of (72,65), (137,7), which now has bound 1 from Fisher's `Θ₄` as well, tier C)".
N4. CLASS.md:65 "The census of the repository (release 0.5.0) is not changed by any of this." -> "The census of the
  repository is not changed by any of this."
N5. "PARI nowhere" (CHANGES.md:3, :8; LEDGER.md:51; census_lvalue_exact.py:26) is true of the certificate (checked
  with cypari2 blocked); the pair -> curve table was found by the PARI-guided census_lvalue_identify.py. Optional:
  "use PARI nowhere in the certificate (the curves were found with PARI, a finder)".
N6. New scripts torsion_exact.py, rank_fisher_theta4_order.py, kctp_isog_exact.py: headers say "not read by a
  person" but not that they were reviewed (LEDGER says reviewed: holds). Not a contradiction; optional.
N7. ci.yml:5-6 header: long-job list omits `square-section`; "`rank` ... about 15 minutes" vs REPRODUCE.md:54
  ("6 minutes at run 3; the run-4 modules add about 15").

### Checked and holding
177 / 147 axiom lines; ten forge files, 14 + 9 pinned lines; 22 forged torsion certificates; torsion 8/8/3;
416 = 268+122+18+6+2, 9 open, 24 at 0.5.0; six (eight) ratios on Magma alone, consistent in README:100/119,
LEDGER:53, G1:343; OUTPUTS 117 = 47/32/26/12 and 106 = 41/30/24/11 (rule N1); action pins = upstream v5.1.0, v6.3.0,
consistent in all 17 uses, comment correct; all new quick lines and all 14 Fable descents lines pass from a clean
clone; PARI-free torsion step confirmed with cypari2 blocked; LEDGER release column 0.6.0 throughout tier A; README
"last release is 0.6.0", REPRODUCE "last tag v0.6.0", CITATION version/date match CHANGES; Θ₄ headers now "reviewed
with gaps"; reviews/ has no secrets and an accurate index; F4 conditional everywhere.

## Final (22:04 UTC)
Verdict: PUBLISH AFTER FIXES. Two MUST-FIX items are release metadata (CITATION abstract; "nothing released" inside
the 0.6.0 entry). No numerical claim found wrong; no CI line fails.

Addenda (22:03 UTC)
- M1: there is no .zenodo.json, so Zenodo's GitHub integration takes the record's metadata (abstract included) from
  CITATION.cff; that is why M1 is MUST.
- N8. reviews/README.md: this report (referee_run5.md) is not yet in reviews/ or its index; if it is committed, add the
  row "| `referee_run5.md` | the integrated diff of run 5 (release candidate 0.6.0) |".
- Forge headers state the expected error counts (TwoDescentCertForge: 14; UniformRankAllForge: 9, with lines) as the
  brief requires; the pinned lines were confirmed by review_x1 on the Mac (not rerun here: no Lean in this container).
