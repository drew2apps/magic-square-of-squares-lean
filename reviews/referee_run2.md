# Referee run2 (04deeaa..b4da2f2)

Started Wed Nov 25 20:04:04 UTC 2026

## Numbers / hygiene (checked 20:06 UTC)
- declared_theorems.txt 961 = recount 961 (922 `theorem` + 39 `private theorem`). OK
- audit_line.txt "AUDIT OK: 1659" = last line of audit_output.txt; audit_output = 949 names + "axioms used" + AUDIT line (951 lines); audit_internal = 710 + summary "949 printed, 710 auxiliary, 1659 in all". 949 = 922+1+5+3+13+5; 5 `curve_a*` names present. OK
- verify_manifest.py: MANIFEST OK 507/507. No __pycache__/.pyc/.lake tracked; `.gitignore` has `.lake/` (unanchored, covers lean_rank/.lake). OK
- MagicSquareOfSquares.lean imports TwoDescentMap, TwoDescentCertificates; no `#print`/sorry/native_decide in MagicSquareOfSquares/TwoDescent*.lean. OK
- ci.yml parses (yaml.safe_load), 9 jobs; rank/quick/build have no `if:`.

## Prose
- FIX-FIRST (wording) LEDGER.md:34: "given that `T` is a nonsingular point (a hypothesis of that statement, proved elsewhere in Lean)". No Lean theorem states T nonsingular at every ratio with m,n,m±n != 0 (nonsingular_T needs qe != 0; qe != 0 appears only inside the proof of Eplus_isElliptic, and theorem2_of only at twice-square ratios). review_t7 says "not yet discharged in Lean". Say "true at every such ratio (qe != 0 is shown inside Eplus_isElliptic) but not discharged in Lean".
- NIT README.md:198 stale after this diff: "`SECOND_DESCENT.md`, section 8, lists the committed outputs that no job reruns" (now OUTPUTS.md, as ci.yml:7 and REPRODUCE.md:60 say), and "`height_precision.py` and `height_certify.py` need Sage and are not run in CI" (now quick ci.yml:133 and four-points). Same stale sentence in G1_FULL_SYSTEM.md:240.
- FIX-FIRST OUTPUTS.md:13: "whose name did not occur anywhere in `ci.yml` at commit `8214be3`": 8214be3 is not an ancestor of b4da2f2 (an earlier version of the lock-in; 04deeaa is the one on this branch), so the published reference dangles. Cite the published 9d0e7e3 (check ci.yml is the same there).
- FIX-FIRST (wording) CHANGES.md:11 "Not integrated" bullet states results of work not in the repo as facts ("it gives 3, like PARI, on all five curves", "rank exactly 3 for `E0L` of (920,129)") while saying the review is still open. Say "an unintegrated, not yet fully reviewed run reports ..." or drop the numbers.
- NIT OUTPUTS.md:3 "Item B4 of the task list": internal reference, meaningless to a reader.

## Dry runs (scratch copy /tmp/ref_run2/wt = git archive b4da2f2, GITHUB_WORKSPACE set, bash -c 'set -e -o pipefail; ...')
PASS: two_descent_certs; selmer_square --negative; selmer_square_rule --max 5000 (62 s); rank_indep (+ --negative, --all) (16 s); lball_exact --census, --tatetest 4000, --rowtest 3000; awk sample diff (cover_descent); wide.py --summary 2000_2500, 2500_3000; genus2/frobenius; indep_rank --ainvs 440,279.
- FIX-FIRST (wording) UNIFORM_RANK.md (section 7 Lean paragraph): "the hypothesis that `T` is a nonsingular point of the curve, which `TwoDescentMap.lean` proves" — same issue as LEDGER.md:34: no theorem of TwoDescentMap.lean states it at every ratio (nonsingular_T needs qe != 0, shown only inside Eplus_isElliptic's proof). Reword ("follows from nonsingular_T, qe_P, genH_pos; not stated as a theorem").
- NIT REPRODUCE.md:60 "each long job ends with `git diff --exit-code`": twist-search has none and scan-5000 has one but not at its end (OUTPUTS.md:7 lists the jobs correctly). Say "the jobs four-points, descents, last-pair, scan-15000 end with ...".
- NIT SELMER_SQUARE.md:3 "item C1 of the task list" (internal reference).
- NIT SECOND_DESCENT_CHAIN.md (440,279) bullet: the interval quoted with the tag [exact script, one] is lball.py's [..., 1.104924976]; lball_exact's is [..., 1.104924975].
- OK: no census count changed; all new dates 25 November 2026 (UTC); other dates are of external sources; rank >= 2 in Lean stated for the 254 curves of m <= 1000 with what is not in Lean (completeness of the 127, Eplus m (-n) = E-, Lemmas B, D); review claims in CHANGES/LEDGER match review_t1/t3/t4/t7; t4's tier-C FIX-FIRST resolved (LEDGER rows in tier C and D); t7's NITs 2 (dead `!`), 4 (lean-action, header) applied; t1's "not in CI" resolved.
- height_precision.py (quick, "no eclib"): run here with eclib/mwrank imports blocked: PASS (1 s).
- All files named by new CI lines exist; cwd tracked for each new line in quick, four-points, descents, last-pair, scan-15000: correct. No new top-level `a && b`, no `head` on a live pipe, no `! cmd`. Background jobs in descents are waited with `wait $PID` (fails under -e).
- FIX-FIRST (wording) README.md:198 audit sentence: "the 5 lemmas `curve_a₁`, …, `curve_a₆` that Lean generates for the definition `curve`" is wrong: they are written in the source, `@[simp] theorem curve_a₁ … curve_a₆` (TwoDescentMap.lean:79-83), and are missed by the line count for the same reason as `alpha_zero` (:218). The numbers hold (949 = 922 + 6 declared after `@[simp]` + 3 + 13 + 5). Say "`alpha_zero` and the 5 lemmas `curve_a₁`, …, `curve_a₆` (declared after an attribute, so not counted by the line count)".
- PASS lproof curve (440,279) and pt_tors lines (quick). All writing scripts in quick gate writes behind --write (or write to /tmp); after the dry runs `git status`/`git diff --exit-code` clean.
- lean-action at the pinned ba3fdfe: has `lake-package-directory` (action.yml:169), `lake build` with no args (no --wfail). Lean not installed here: the build/rank Lean steps (TwoDescentChecks, RankForge lines 19/23/27/33, 9 `depends on axioms`) not rerun; rely on review_t3/t7.
- Timeouts: added foreground work small (quick +~2-3 min vs 30-min limit; descents +~7 min foreground, 16 min background, 90-min limit; four-points +~4.5 min). OK.
- Pre-existing, out of range (not a FIX-FIRST for this diff): ci.yml:135 `sigma_exact.py --negative | grep -q` on a live pipe.
- FULL quick step ("transcript checkers ...") extracted from ci.yml and run on a fresh `git archive b4da2f2` copy with GITHUB_WORKSPACE set, `bash -e -o pipefail`, pip lines skipped (local: cypari2 2.2.9, sympy 1.14.0, passagemath 10.8.14, numpy 2.5.3 instead of the pinned 2.5.5): EXIT 0 in 17 min, final `git diff --exit-code` clean.

## Verdict (20:26 UTC)
FAIL as it stands, on wording only: README.md:198 (curve_a lemmas "generated by Lean"), LEDGER.md:34 + UNIFORM_RANK.md:85 (T nonsingular "proved in Lean"), OUTPUTS.md:13 (dangling 8214be3), CHANGES.md:11 (unintegrated results stated as facts). PASS once these are fixed; no number, CI or Lean change needed. NITs above.
