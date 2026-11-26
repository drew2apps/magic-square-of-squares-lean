# Referee report: integrate branch (claims editor / reproducibility)

Repo /home/claude/review/repo2, main d2fe4c6 vs integrate 060ece0. Findings appended as found.

## Progress notes (appended as found)

- Reviews s1/s3/s4/s6 read. Diff main..HEAD: 51 files. YAML parses (8 jobs).
- Lean: root imports ClassHensel, G1SigmaGeneric, SecondDescentIdentities (39 imports = 39 files in MagicSquareOfSquares/).
  No sorry/admit/native_decide/axiom/implemented_by/extern/unsafe/opaque in the three files; `set_option maxHeartbeats`
  twice in G1SigmaGeneric (harmless); `decide +kernel` only on the ClassHensel controls. Count
  `cat MagicSquareOfSquares/*.lean | grep -c '^theorem \|^private theorem '` = 642 = declared_theorems.txt
  (568 + 29 + 14 + 31; public 603, private 39).
- NEW: the three new library files contain 39 `#print axioms` commands (6 + 14 + 19); no library file on main has any.
  Harmless for soundness (Audit.lean covers them) but noisy build output and against house style.
- MANIFEST.sha256 not regenerated: verify_manifest.py on HEAD: 12 hash diffs + 39 unlisted files, so the first step of
  `quick` fails as committed. (Treat as "regenerate last", with the audit files.)
- Quick-job lines run in a clean `git archive` copy: sigma_generic 37 s PASS; exact_last_pair (+ --negative) 145 s PASS;
  rank_indep (+ --negative) 5 s PASS. After the run the committed .txt and hints .json in the copy are byte-identical
  to HEAD, so the outputs genuinely reproduce.
- NEW, CI design defect: exact_last_pair.py, rank_indep.py, rank_indep_heights.py (and their --negative modes) WRITE
  the very file they are diffed against (exact_last_pair.py:945,948; rank_indep.py:726,743; rank_indep_heights.py:113,135)
  and rewrite the hints JSON (rank_indep.py:727,741; rank_indep_heights.py:112,134; exact_last_pair.py:721).
  `python3 X.py | diff - X.txt` then races the rewrite (mutation test below). sigma_generic.py only writes with --write (fine).
- rank_indep_heights (+ --negative) 66 s PASS. All four new quick lines pass on a clean copy and the rewritten files are
  byte-identical to HEAD (so the results do reproduce).
- **MUTATION TEST (blocking): the diff in those CI lines is vacuous.** In a fresh copy I changed the committed
  `scripts/g1_full/rank_indep.txt` summary to `760,39 E+ <= 1` (a false claim that would remove the "Magma alone"
  caveat) and ran the exact CI line `python3 rank_indep.py | diff - rank_indep.txt`: exit 0, three times out of three.
  The script prints with flush, diff blocks on stdin, the script then rewrites rank_indep.txt before exit, and diff
  compares the output with itself. So these CI lines cannot fail on a changed result (only on a crash/assertion).
  Fix: write outputs only under an explicit `--write` flag (as sigma_generic.py already does), or have CI copy the
  reference to /tmp first, or run `git diff --exit-code scripts/` after the scripts. Same for the hints JSON.
- Mutation results (forged committed file, then the exact CI line): rank_indep.txt exit 0; rank_indep --negative exit 0;
  rank_indep_heights --negative exit 0; exact_last_pair --negative exit 0 (all VACUOUS); exact_last_pair.txt and
  rank_indep_heights.txt caught (exit 1) only because the forged line was in the first 4 KB. Synthetic reproduction
  (/tmp/ri1/syn/gen.py, GNU diff 3.10): a script that prints with flush and then rewrites its reference file passes
  `| diff - ref` with ANY change when the output is < 4 KB, and with any change after the first ~4 KB block otherwise.
  So 6 of the 8 new quick-job diffs check nothing, and the other 2 check their first 4 KB only. BLOCKING.
- descents-job lines (clean copy, cwd scripts/second_descent): indep_class 91 s PASS; indep_fourclass 0 s; indep_rank 5 s;
  indep_second 2 s; indep_second --selftest 1 s; indep_rank --selftest 0 s (all PASS; these scripts only print, no
  self-writing). indep_rank*.txt contain PARI's ellrank/elltors/conductor (library output, diffed verbatim; deterministic
  with the pinned cypari2 2.2.9, but it is the kind of output the house rule says to grep, not diff).
- indep_class --selftest 360 s PASS. All 7 new descents lines and all 4 new quick lines pass now (on the real outputs).
- class_hensel_witness.py: exit 0, "all checks passed"; its 40 big integers = the 40 in ClassHensel.lean. Not in CI
  (fine: the kernel checks the witnesses). py_compile of every new script OK.

## Headline numbers against their source files (all agree unless marked)
- 20 forged inputs (LAST_PAIR §18, CHANGES, LEDGER): exact_last_pair_negative.txt "20 of 20 refused". OK.
- 13,406 / 26,813 / c1 197.62 / c2 177.86 / mu 3651.3 / h2(C) 7.01e5 / hhat(A,B) 5.39e9 / hhat(G) [20.002899, 20.009091] /
  4.444 / 110 / 210^2 / 0.443 / Sel 2 and 64 / 15 cosets / sieve sizes 1632, 204, 909, 35858 / 596,309: exact_last_pair.txt. OK.
- 4 of 12 curves (G generates unconditionally), 8 conditional: rank_indep_heights.txt summary. OK. rank 3 for the four
  Magma-alone curves: rank_indep.txt; 1 for (91,60) E+/E-, (572,555) E-, (741,580) E+, 3 for the other eight:
  rank_indep_all.txt. OK. (187,84) ok, (440,279) and (420,341) rank <= 3: exact_last_pair_route.txt. OK.
- 36 of 36 (indep_fourclass.txt), 12 planted, 4 mutated tuples (indep_class_selftest.txt), 131,072/524,288 and the nine
  "primes of last change" 37,37,41,41,29,29,37,37,37 and the d values 1, 37241, 661762, 713842, 1 (indep_class.txt):
  all OK. 50 (=R+s) + 40 (>= lower bound) curves: indep_rank_selftest.txt OK. Torsion 4,4,2,2 / 8,4,4,4: indep_rank.txt OK.
  Sel 4 and 16, bound 4 for (609,200): indep_second.txt OK; bound 2 for (440,279): indep_rank_440_279.txt OK.
- UNIFORM_RANK: 127/254, 200 random (largest 983351), 35 square-type, 1085 coprime m < 60, 84 forged, 186 quartics
  (167/19): sigma_generic.txt OK. "under a minute": 37 s OK. LMFDB isomorphism identities rechecked by me in sympy: true.
- 642 = grep count = declared_theorems.txt (568 + 29 + 14 + 31). No prose states 642 (README:193 still says 568, see below).
- "31/35": occurs in no changed file; nothing to check.
- Reviewer counts quoted in prose: ">320 curves" (s3: 160+63+9+40+60) OK; "868 classes" OK; BUT LAST_PAIR.md:451 "300 curves
  of rank 2 and 3 (the bound never fell below the true rank" is wrong: s4 ran 40 + 60 curves with PARI rank >= 2 and
  200 of ANY rank, and checked "not below PARI's lower bound", not the true rank.

## Stale or contradictory sentences
- LEDGER.md:45 (tier B, the 18) "Sage's lower bound for hhat (Cremona-Siksek) is the one floating-point input left" vs
  LEDGER.md:93 (tier E) "No longer used ... for the 12 generators of the 18 closures". Same in G1_FULL_SYSTEM.md:344.
  Also note sigma_exact.py (the CI route of the census) still takes lam from the evidence; say "an exact alternative
  (rank_indep_heights.py) removes the need for it".
- SECOND_DESCENT.md:77 "Not done: a third program for the class step ..." now contradicted by the paragraph right after it.
- SECOND_DESCENT.md:93 "exclusions are certified and re-verified, and PARI is used only to factor" (the chain's step 7
  was corrected: certificates only for 2-isogeny descents; PARI also factors polynomials and runs polrootsmod).
  Same "only to factor integers" in indep_rank.py:35.
- AUDIT_REQUEST.md:67 (e) still asks outside readers to check Sage's 0.045162 for Theorem A, which section 18 no longer
  uses; the question should now be CPS Thm 15/Thm 2, minimality, m = 24 and the search to 110.
- LAST_PAIR.md:100 (table "0.045162 | saturation.py"), :128 (item 5, "the lower bound for the height is Sage's"), :231
  ("Lemma 1: PARI, Sage, Magma") contradict section 18 without a pointer to it.
- BARRIER.md:60-61 "This covers the 127 ratios to 1000 only" now sits two lines above the UNIFORM_RANK sentence that
  says the opposite; reword.
- README.md:118 / G1_FULL_SYSTEM.md:344 "for the other eight ... rests on PARI and Magma": for four of the eight PARI
  says [1,3], so Magma alone; rank_indep_heights.py:131 summary "PARI/Magma/eclib": eclib gave no bound on any of the 8.
- README Lean table (README.md:130-173) has no rows for G1SigmaGeneric.lean and SecondDescentIdentities.lean (nor for
  the older G1Sigma.lean).
- SECOND_DESCENT_CHAIN.md:3, :29, :71 name internal "slice s3"/"slice s2" and promise s2's Lean class step; a reader
  cannot see either. Drop the slice names and the promise.
- LEDGER tier E (LEDGER.md:87-98) lacks the newly quoted inputs: LMFDB 576.a.576.1, CPS 2006 Thms 2 and 15 (only the
  constants are listed), Silverman 1990 Thm 1, Kodaira-Neron, Silverman AEC X.4 / VII.3.1.

## Reproducibility gaps (no CI job regenerates these new committed outputs; not listed in SECOND_DESCENT.md section 8)
- scripts/g1_full/rank_indep_all.txt, rank_indep_control.txt, rank_indep_control_census.txt, sigma_generic_pari.txt;
  scripts/second_descent/exact_last_pair_route.txt (its script APPENDS, exact_last_pair_route.py:74, so a rerun
  duplicates it), indep_rank_440_279.txt. Add them to the list of SECOND_DESCENT.md:205 or to CI.
- indep_class.py:570-573 counts ANY AssertionError from descent() as "removed"; check the message
  ("planted tuple removed at p =").
- The three new library files carry 39 `#print axioms` (house style: none in the library; the audit covers them).

## Lean status (gating, not a defect)
- ClassHensel.lean: built at ad737f3 (sha 16599247..), committed d06f523 differs by one comment line only. OK.
- G1SigmaGeneric.lean: compiled by the conductor at c73a9ea = committed sha 9e55153f. OK (record not in the repo).
- SecondDescentIdentities.lean: the committed file (sha 5403637a) was NEVER compiled; s3's compile was of 4d3643e7,
  before a theorem was deleted and docstrings changed ("not recompiled -- Mac offline", reports/s3.md:151). Its
  tier-A row (LEDGER.md:33) and every [Lean] tag in SECOND_DESCENT_CHAIN.md rest on the pending build.
- To regenerate on the build machine after `lake build` passes: scripts/audit_output.txt, scripts/audit_line.txt,
  scripts/audit_internal.txt; README.md:193 ("is 1177 with G1Sigma.lean, against 568"; "529 and 39" -> 603 and 39;
  "The audit prints 550 names: the 529 public theorems ..."; "lists by name the 627"); optionally the CHANGES entry
  (house style gives declared count and audit line); MANIFEST.sha256 last (now: 12 hash diffs, 39 unlisted, so `quick`
  fails at its first step).
- Line-number corrections: the "(PARI, Magma)" remark is G1_FULL_SYSTEM.md:343 and the floating-point sentence :344;
  the heights summary is rank_indep_heights.py:132-133.

## Verdict: FIX-FIRST
1. (blocking) Self-writing scripts make 6 of the 8 new quick diffs vacuous and the other 2 partial (mutation test above).
2. (gating) SecondDescentIdentities.lean as committed never compiled; build + audit regeneration + manifest pending.
3. Stale/contradictory: LEDGER.md:45 vs :93; G1_FULL_SYSTEM.md:344; SECOND_DESCENT.md:77, :93; AUDIT_REQUEST.md:67(e);
   LAST_PAIR.md:100, :128, :231; BARRIER.md:60-61.
4. LAST_PAIR.md:451 misreports the review ("300 curves of rank 2 and 3", "true rank").
5. Six new committed outputs regenerated by no job and not listed (SECOND_DESCENT.md:205).
6. Minor: README Lean table rows, LEDGER tier E entries, "PARI and Magma" for the eight, slice names in the chain,
   indep_class.py:570-573 catch-all, 39 #print axioms in library files.
Finished.
