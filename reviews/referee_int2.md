# Referee report 2: integrate branch at dc56fe9 (after FIX-FIRST on 060ece0)

Repo /home/claude/review/repo2, HEAD dc56fe9, worktree clean before and after. All runs in clones (/tmp/r2c, /tmp/mut/*);
the repo itself was not modified.

## 1. Issues of referee_int1

1. (blocking) Self-writing scripts: FIXED (b1847ae). Every write in exact_last_pair.py, rank_indep.py, rank_indep_heights.py,
   rank_indep_control.py, exact_last_pair_route.py (outputs and hints JSON) is now under `if WRITE` (`--write` in argv);
   save_hints() is gated; route script prints and only appends with --write; sigma_exact.py (imported) writes nothing.
   - The three quick CI lines (exact_last_pair +--negative, rank_indep +--negative, rank_indep_heights +--negative), run
     with pipefail from the CI directories in a clean clone: all exit 0; `git status --untracked-files=all` clean afterwards.
   - Mutation test, one fresh clone per file, last line forged (past 4 KB for exact_last_pair.txt and rank_indep_heights.txt;
     rank_indep.txt summary forged to "<= 1"): all 6 CI diffs now FAIL (exit 1), and the forged file stays forged
     (script no longer overwrites it). The diffs are no longer vacuous.
2. (gating) Lean build: no Lean toolchain here, so not rebuilt by me. The regenerated audit is consistent with the current
   sources (section 2): the 74 added printed names are exactly the 74 public theorems of the three new files, and no name
   from the pre-deletion SecondDescentIdentities appears. Lean sources unchanged since 071fc05 (which removes only the
   39 `#print axioms` lines). CI's Lean job on push is the remaining confirmation.
3. Stale sentences: FIXED. LEDGER.md:45 now "the one floating-point input left in this route; an exact alternative ...",
   consistent with :93; G1_FULL_SYSTEM.md:343-344; SECOND_DESCENT.md "Not done: a third program" gone, :93 PARI wording
   corrected (also indep_rank.py:35); AUDIT_REQUEST.md (e) now asks about CPS Thm 15/Thm 2, minimality, m = 24, 4.444, 110;
   LAST_PAIR.md:100, :128, :231 point to section 18; BARRIER.md:60-61 reworded.
4. LAST_PAIR.md:451: FIXED ("300 curves (100 with PARI rank at least 2, 200 of any rank; ... PARI's lower bound").
5. Uncovered outputs: FIXED, all six listed with their commands in SECOND_DESCENT.md:205.
6. Minor: README Lean rows for G1Sigma/G1SigmaGeneric and SecondDescentIdentities added; tier E rows added (LMFDB
   576.a.576.1; CPS 2006 Thms 2 and 15, Silverman 1990, Kodaira-Neron, AEC X.4 / VII.3.1); "PARI and Magma" for the
   eight fixed in README and G1_FULL_SYSTEM; indep_class.py:570-575 now counts only the descent's own
   ("planted tuple removed at p =", p) assertion (matches line 408) and re-raises others; `#print axioms` removed.
   Residuals (non-blocking):
   - scripts/g1_full/rank_indep_heights.py:135 and rank_indep_heights.txt:14 summary still "conditional on rank 1 from
     PARI/Magma/eclib"; eclib gave no bound on any of the 8 (the per-curve lines say so honestly). Fix needs a
     regenerated .txt and manifest.
   - SECOND_DESCENT_CHAIN.md:35 "the review of the slice", :42 "the review of slice s3", :47 "The review of the slice":
     internal slice names remain (lines 3, 29, 71 were fixed).
   - CHANGES.md has no entry with 642 declared / AUDIT OK 1315 (house style; optional).

## 2. Audit files

- audit_line.txt = last line of audit_output.txt = "AUDIT OK: 1315 theorems (including auxiliary ones), standard axioms
  only"; axioms line [Classical.choice, Quot.sound, propext].
- audit_output.txt: 624 printed names, unique, sorted; old 550 names all still present; 74 added = the 74 public
  theorems of ClassHensel/G1SigmaGeneric/SecondDescentIdentities (0 private there); none missing.
- 624 = 603 public declared + 21 non-theorem entries (3 fact instances, 5 eq_def, 13 structure fields), the same 21 as
  before.
- audit_internal.txt: 691 unique entries, disjoint from the printed ones, 624 + 691 = 1315, trailer matches; old 627 all
  present, 64 added (all ClassHensel); all 39 declared private theorems present as _private names.
- README.md:195 matches: 1315, against 642, 691, 39 private, 603 and 39, 624 names = 603 + 3 + 13 + 5.

## 3. Manifest, counts, forbidden tokens

- `python3 scripts/verify_manifest.py`: MANIFEST OK, 472 listed, 472 tracked (also OK in the clone after all runs).
- grep count of `^theorem |^private theorem ` = 642 = scripts/declared_theorems.txt.
- No `#print`, sorry, admit, native_decide, axiom, implemented_by, extern, unsafe, opaque in the three new files (no
  `#print` anywhere in the library).

## 4. CI

- ci.yml parses (8 jobs). No plain `cd` in the quick step before the new lines; new quick lines use root paths or
  `(cd ... )` subshells; the descents step is in scripts/second_descent, where indep_*.py and magma_second_ranks.json live.
- Dry runs from the CI directories, clean clone: sigma_generic PASS; the three slice-4 lines PASS; indep_class (9 ratios)
  PASS; indep_class --selftest PASS; indep_fourclass, indep_rank, indep_rank --selftest, indep_second,
  indep_second --selftest PASS. git status clean afterwards.

## Verdict: PASS

Nothing blocking. Three minor residuals, which can wait for the next commit: rank_indep_heights.py:135 (and .txt:14)
"PARI/Magma/eclib"; SECOND_DESCENT_CHAIN.md:35, :42, :47 slice names; no CHANGES entry for 642 / 1315. The Lean build
itself is confirmed only by the consistency of the audit files; CI's Lean job is the confirmation.
