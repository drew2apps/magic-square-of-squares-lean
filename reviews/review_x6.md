# x6 (run 6) — hostile review, then packaging, of the Roadmap G bundle (experiments/roadmap-g)

Written by an AI model (review agent x6); not read by a person. Committed as written (paths are those of the machine it ran on).

Started 04:20 UTC 27 Nov 2026. Worktree /home/claude/conductor/x6, branch run6_x6 on run6_x4 (f54bfac).
Checkpoint 06:30 UTC, hard stop 08:30 UTC.

## Review

Inputs: the bundle /home/claude/untrusted/roadmapg/experiments/roadmap-g (copied to /home/claude/x6work/roadmap-g; every
script I ran was read first: they compute, import the repository's modules through MSS_REPO, and write only into their
own folder and only with --write). Repository = worktree x6 at f54bfac (tree of 0.6.0 + run6_x4). This container: PARI
2.17.4, cypari2, sympy 1.14.0, numpy 2.5.5, 2 cores. My own test scripts: /home/claude/x6work/mytests/ (outputs copied
to /home/claude/x6work/evidence/). They are committed, with relative paths, in experiments/referee_x6/ under these
names: cmp30.py + cmp_more.py -> local_compare.py; weilcheck.py -> weil_check.py; recheck5.py, recheck63.py ->
kept63_repo_only.py; g2_recheck.py -> killsets_witness.py; g2_witness_all.py -> badalone_witness_all.py; control.py ->
control_1109503.py; fam_check.py -> g3_family.py; g3_sample.py; g4_enum.py; the stage reruns (cw/, route/) ->
samples/; new after the first draft: local_mutation.py.

### 1. G1 exact_local.py (the exact local test at an odd prime) -- HOLDS (with gaps of wording)

Read line by line against tuple_descent.Local and LAST_PAIR.md section 14 (R1, R3, R4).
- Charts (t,1), t in Z_p, and (1,px), x in Z_p, cover P^1(Q_p); sub-discs are exact substitutions A+Bx0, pB in Z[i]
  (no truncation). Labels and the residue-symbol convention are those of Local.dlabel (0 = square), and the closure
  under a non-square unit of Z_p is applied for p = 1 mod 4 only, as in Local (units are squares in Q_p(i) for p = 3 mod 4).
- p = 1 mod 4: each component a + b x; v(a) < v(b) gives a constant class (R1); otherwise p^v(b)(a' + b'x) and on a
  residue x != r the label is (v(b) mod 2, chi(b') xor chi(x - r)) -- checked. A tuple is emitted only when the two
  components of every form agree (the class of a rational number). Root residue x0: if all functions vanishing there
  have exactly the same root (A1 B2 = A2 B1, or A1 conj(B2) = conj(A2) B1 for opposite components: I checked this
  algebra; comp_s is injective on Z[i] since iota is irrational), x - x* = p^e w takes every e >= 1 and every unit w,
  so the four (parity, symbol) pairs are all realised and shared by the functions with that root -- correct; otherwise
  the disc is refined (roots that differ separate at the depth of their difference; depth < 400 asserted).
- p = 3 mod 4: label (v mod 2, 0), alive iff the norm of the unit part is a square (R3); a form with root rho = r0 + i r1,
  r1 != 0, is a unit on all of F_p and needs chi((x - r0)^2 + r1^2) = chi(N(b')); forms with the same {rho, conj rho}
  are grouped and conflicting wants kill the residues -- correct. No form has a root in Z_p (it would be a rational zero;
  F0, F8, GL, GU have none for m != n), so the refinement ends.
- p | coefficients: valuations are pulled out of a and b separately (va >= vb >= 1 handled), so p^k | m, n, m +- n,
  m^2 + n^2 or the octic needs no special case; precision: values are mod p^N, a valuation >= N - 2 raises and the
  whole prime is recomputed with N doubled (allowed_at).
- **The Weil step.** For g distinct roots the number of x outside the roots with a given sign pattern is
  2^-g [ (p - g) + sum over 2^g - 1 nonempty subsets S of +- sum_x chi(prod_S (x - r)) ]; each product is squarefree of
  degree |S|, not c*square, so |sum over F_p| <= (|S| - 1) sqrt p, and removing the g excluded points costs at most g.
  weil_all checks p - g > sum_S ((|S| - 1)(isqrt(p) + 1) + g): this is that inequality with sqrt p rounded up --
  correct and conservative. Same for p = 3 mod 4 with irreducible quadratics (2|S| roots). The number of characters is
  2^g - 1, g <= 8 (four forms, two components). Thresholds computed from weil_all itself (weilcheck.out): 8 roots from
  p >= 595,717; 7 from 105,222; 6 from 17,542; quadratics (p = 3 mod 4) at most 2,564. So the assertion "p <= 5*10^6
  when listing is needed" can never fire. **Direction:** the Weil step only ever adds patterns, so an error there
  would make a set too large, never cause a false exclusion; it matters in the other direction for G2's claim that
  bad primes do not suffice, and for that it is right.
- **Wording gap:** the docstring and REPORT_G1 say patterns are "settled by Weil's bound for p > 5*10^6 and by listing
  below that". In the code Weil is tried first at every p and decides whenever its inequality holds (from p = 3 for one
  root); listing (numpy, O(p)) is used only below the threshold, and 5*10^6 is only the cap asserted for listing.
  MixedLocal is dead code. At p = 2 the repository's one-sided test (depth 12) is kept, as said.
- **Can it answer "no local point" wrongly?** Only through a wrong label, a missed residue or a precision slip; the
  structure above visits every residue of every disc, and the following found no case:
  - test_exact_local.py rerun on PARI 2.17.4: `# matches test_exact_local.txt` (2.4 s; 118 (ratio, prime) cases equal
    to the repository's set, 240,000 random points).
  - my cmp30.py / cmp_more.py (own label code: own sqrt(-1) mod p^60 from PARI and Newton; own valuations and symbols;
    random points in both charts and points t* + p^e w, e = 1..11, around every p-adic root of every component, or
    of the real part of the Gaussian root for p = 3 mod 4): **107 (ratio, prime) cases**, ratios (627,364), (291,224),
    (1109503,597671), (993260025,766111057), (10125,4459) (3^4 | m, 5^3 | m, 7^3 | n), (626,1) (5^4 | m - n), (7,4)
    (5 and 13 exactly once in m^2 + n^2: both sets empty, as Proposition A says), (1189,820), (4999,1802); primes from 3
    to 2^61 - 1, every odd bad prime of the last six ratios (octic primes up to 1.9*10^15 included). Wherever the
    repository's test could run (p up to about 2*10^5, depth 2 or 3) **the exact set equals it**; and in **all 107 cases
    the labels of explicit points are exactly the exact set** (nothing outside it, every element of it realised). So on
    these cases it is neither too small nor too large.
  - weilcheck.py: 593 (ratio, prime) cases with weil_all forced False (listing only) against the default: 0 differences.
  - local_mutation.py: five mutations of exact_local.py that would make a set too small (tail with one parity of e; no
    unit closure; listing keeps one pattern; dead test at p = 3 mod 4 reversed; symbol of b' forgotten): each is caught
    by explicit points (on 1, 12, 8, 9, 7 of 24 cases); the unmutated code on none. So the check has teeth (the tail
    mutation only weakly).
  Verdict: HOLDS as a computed (exact) local test, tested, not proved; written by an AI model and read by no person.

### 2. G1 stress test -- HOLDS WITH GAPS (one false number, one overstatement)

- Byte-for-byte reruns from the seeds on PARI 2.17.4: **k = 13** (17 s) and **k = 29** (804 s, the largest range):
  both `# matches g1_k<k>.txt`. (Also --negative and the summary, see the CI section.)
- **All 63 kept ratios** (not only 5) through the descent twice (recheck5.out, recheck63.out): (a) exact_local
  descent with the supports and good primes <= 150: nothing left for 63/63; (b) **the repository's own
  tuple_descent.Local only**, at the same primes up to 2*10^5 with its own depths (12 at 2, 6 below 50, 3 above), the
  larger support primes simply left out (30 of the 63 lose one to three of them): **nothing left for 63/63**. So G1's
  "no survivor" does not depend on exact_local at all; it rests on the repository's criterion (class_selmer,
  class_selmer_linalg, unchanged) and the repository's four-class test, as in COMPLETE.md section 7.
- Scope as stated is right: 340,000 random coprime ratios (20,000 per dyadic range, seed 1000 + k), the criterion, and
  for the 63 kept ones the four-class descent at the primes of the supports and the good primes <= 150, with the 2-adic
  test at depth 12. **But REPORT_G1's proposed LEDGER text says "exact_local.py, no prime left out": false for G1** --
  g1_stress.py tests only the support primes and primes <= pmax, so primes of m^2 + n^2 and of the octic above 150 are
  left out (harmless for "nothing is left", since leaving a prime out can only keep more tuples). "No prime left out"
  is true of G2, not G1.
- The control (1109503,597671): gcd 1, m and n odd, m^2 + n^2 = 1,588,207,531,250 = 2 * 891125^2 -- so twice-square
  type, coprime, with the rational points t = 0 and t = infinity, hence points everywhere for its two degenerate
  tuples. **REPORT_G1 says m^2 + n^2 = 2*(5*7129)^2: wrong; 891125 = 5^3 * 7129** (2*(5*7129)^2 = 2,541,132,050).
  The six tuples are reproduced (control.out); each has an explicit 2-adic point (G2's witnessed_at_2 and my own
  search) and explicit p-adic points at every odd bad prime (my own label code). See the control paragraph below for
  the good primes.
- **Overstatement:** REPORT_G1 and the hand-off say "the descent can leave classes that no point explains" and
  "survive every good prime up to 20,000". What is computed is: the four extra tuples pass every bad prime and every
  good prime <= 20,000. That they come from no rational point is not shown (it would need the rational points of the
  curve), and "survive" is bounded by the primes tested. Proposed wording: "the four-class descent, run to good primes
  20,000, keeps four tuples besides those of the degenerate points; whether a rational point has one of these classes
  is open".
- **Control, further (control_1109503.txt):** explicit points for all six tuples at 2 and at all 18 bad primes (up to
  1.19*10^15), and exact_local at the 23,397 good primes = 1 mod 4 in (20000, 600000]: no tuple fails (807 s). Above
  595,717 every tuple passes at a good prime (8 distinct root residues, Weil gives every pattern; checked that no root is
  at infinity: the alpha components are units at good p), and at good p = 3 mod 4 all candidate tuples have the labels of
  the degenerate ones. So the four extra tuples are **locally soluble everywhere** (computed + this short argument): they
  are in the set that the four-class descent with every local condition leaves. Whether a rational point has one of
  these classes is open -- so the true statement is stronger than G1's "good primes to 20,000" and weaker than its "no
  point explains".

### 3. G2 -- HOLDS

- g2_check.py, g2_features.py (both lists), g2_check.py --negative: all `# matches` on PARI 2.17.4 (0.1 s, 0.1 s,
  0.6 s, 3.9 s).
- **20 kill sets recomputed** (g2_recheck.out, my own driver; 10 "bad only", 10 with good primes; random, seed 66):
  emptiness with **the repository's Local** at every kill prime <= 3000 (19 of the 20 sets need nothing else; 4241
  of (4241,631) by exact_local); **irredundancy proved, not just tested:** after dropping each prime, one tuple left
  has explicit points (my code) at every remaining prime of the set -- 20/20 sets, every drop.
- **10 "bad primes alone fail" ratios** (random, from both lists; 13 to 17 bad primes each, the largest up to
  1.2*10^17): every tuple the bad primes leave (2 to 10 per ratio, 34 in all) has an explicit 2-adic point and an
  explicit p-adic point at every odd bad prime, found by my own search and label code. So the negative answer to
  test 6 holds for these independently of exact_local's Weil step.
- **Then all 937** (badalone_witness_all.txt, about 30 min): the tuples left by the bad primes number exactly as G2
  says for 937 of 937 ratios, 6,904 tuples in all, and **every one of the 6,904 has explicit local points at 2 and at
  every odd bad prime**. So "the bad primes alone do not suffice" for 101 of the 997 and 836 of the 8,739 is computed
  exactly, independently of exact_local (which was used only to list the tuples).
- The feature scan is reported honestly: counts per feature, the rule (holds for >= 5 ratios and forces one outcome)
  stated, nothing passes it, no conjecture. "Every good prime needed is 1 mod 4 -- forced" is right: at a good prime
  = 3 mod 4 every candidate d has the label (0, 0).
- Gap (wording): "kill sets ... not proved of least size" is said; good. "Every bad prime is tested exactly" rests on
  exact_local (item 1) for primes the repository cannot reach.

### 4. G3 -- HOLDS WITH GAPS (one stale label, one misleading "next step")

- Family regenerated by brute force over all pairs (fam_check.out, numpy, my code): 611 ratios in 11000 < m <= 13000,
  identical to fam_11000_30000.json there.
- g3_table.py --negative: `# matches` (0.3 s). **g3_table.py (default) FAILS against the repository**: the bundle's
  TSV was built from F0's table before run6_x4 fixed the label of (6903,5104); with the repository's table, row
  (12007,1799) (closed through (6903,5104)) gets "CONDITIONAL on a rank being 1 (D1 of the partner: root number -1,
  PARI bounds [1,3]); plus the DM generator condition" instead of "...; plus (a), (b)". The new label is the right
  one; packaged with --write (a label propagated from the fixed table, no mathematical change). All other 610 rows and
  g3_table.txt are unchanged.
- 10 random rows (seed 2026) and all 11 open/special rows against the raw records (g3_sample.out): type right, the
  step's record present and "closed" (or open in every stage for the 4 open ones), partner arithmetic right.
- Reruns on PARI 2.17.4 with the repository's programs (paths patched in a scratch copy only): second descent
  (wide_pari.py) for (11383,7769), (12719,10001), (12927,5336), (11165,7548): 4/4 identical apart from "sec";
  global.py for (11009,840), (11207,6983): identical; dm_pari.py for (11300,10269), (11177,553) (16 records with
  partners): identical apart from "sec"; route_pari.py (12416,5313), d = 251858174: identical apart from "sec", closed.
- The 4 open and 3 conditional are as stated; the route is counted only for square type (g3_table.py drops route logs
  of twice-type pairs; log_11441.txt is kept as a record and not counted); (12487,4241) is closed through its
  square-type partner (8364,4123), closed by the route in F0's table -- consistent.
- **Misleading "next step":** "its partner (18713,3617) ... is reached in block 18000" (also 19000, 22000). The
  partners of the open ratios were already run through second descent, rank 0 and DM in G3 (as partners: sd, r0 and
  dm records exist for (18713,3617), (19561,4841), (22057,2351)), all "open"/"rank not 1". Reaching their blocks reruns
  the same programs, so nothing new is expected from that alone.
- Cosmetic: a detail "curve E-partner" without "partner (...)" means the curve E of the pair itself, computed when it
  ran as another pair's partner; and the detail of (12007,1799) quotes F0's detail of (6903,5104), "D1 of the partner
  (12007,1799)", i.e. its own D1.
- Bookkeeping only, PARI rank bounds (not proved), DM on PARI generators: labels say so per row.

### 5. G4 -- HOLDS WITH GAPS (two sentences overstate)

- g4_summary.py, g4_search.py --negative, gen_square.py 20000: all `# matches` on PARI 2.17.4 (0.5 s, 4.0 s, 8.6 s).
- REPORT_G4 labels it "computed (numerical): a search, not a certificate" and lists the reasons (PARI ranks, saturation
  at primes <= 100, floating-point heights and Fincke-Pohst bound, coverage). Every t is tested exactly (Fraction,
  global.classify_t). All 35,310 recorded hits are t with five squares and verdict "degenerate" (my tally); no 4-square
  hit, no solution.
- **Overstatements:** the hand-off's bottom line "no solution, and no non-degenerate t with four of the five forms
  squares, on the 3,573 coprime square-type ratios" and REPORT_G4's proposed LEDGER sentence "No solution, and no
  non-degenerate t with four of the five forms squares." read as statements about the ratios. They are about the
  points listed. Corrected text: "None of the points listed gives a solution, or a non-degenerate t with four of the
  five forms squares." (g4_summary.txt's line "SOLUTIONS: none" is a header of the bookkeeping and is explained in
  NOTES.md.)
- The planted test does exercise the search code (curve models, to_t, the Fincke-Pohst list, the exact tests and the
  assertion), but weakly: the point of t = +-1 is in PARI's list or is appended to the generators before saturation, so
  its coefficient vector is (nearly) a unit vector. My own check (g4_enum.out): on the curves searched for 7 ratios,
  every combination sum a_i g_i with |a_i| <= 4 and height <= 0.999 H (9,495 points) is in the list `points()`
  returns: none missing. So the enumeration is not just "the generators back".
- Not reproducible bit for bit on PARI 2.17.4 (as REPORT_G4 says): rerun of 5 ratios (g4_rerun.jsonl): (12,5) and
  (40,9) identical; (195,28), (440,279), (7735,4968) differ on one quotient curve each (PARI's points / effort; on
  (195,28) w*f8*K is searched under 2.19.1 and not under 2.17.4, so 6 hits instead of 8, the two missing ones the
  degenerate t = +-1 on that curve). No solution and no non-degenerate hit in either. The committed search
  outputs are 2.19.1 data; CI can check only the bookkeeping and --negative.

### 6. Wording against LEDGER.md tiers

- G1, G2: tier D ("experiment, computed"); the local test at odd primes is exact arithmetic with a written argument
  (exact_local.py docstring and the Weil step), tested, read by no person -- tier D, not C (no written proof in the
  repository's sense of a document of proof; I would not call exact_local "proved"). G1's "no survivor" is tier D and
  does not need exact_local (item 2). G2's "test 6: no" is tier D; for 10 ratios independently witnessed (item 3; and
  see the full witness run below).
- G3: tier D, like the 1000 < m <= 11000 table: per-row labels; conditional rows say on what; "open" rows stay open.
- G4: tier D as a **search** (like the LEDGER row "Searches"), "computed (numerical)", not a certificate, not
  reproducible across PARI versions except its bookkeeping.
- Words to avoid in the packaged text: "decided" for test 6 is acceptable only as "answered (computed)"; "the descent
  can leave classes that no point explains" (G1) is an inference beyond the computation (item 2); "no prime left out"
  (G1); "no solution" without "among the points listed" (G4). No "first" claims, no census change in any of the four.

## Package (commit on run6_x6)

- `experiments/roadmap-g/g1` ... `g4`: the bundle's files unchanged, except: `g3/CENSUS_TABLE_G3.tsv` regenerated with
  `--write` against the repository's F0 table (one row, (12007,1799), its label); `g2/run_g2.sh` reads the repository's
  two lists; left out: `g2/tuple_selmer_997.txt`, `g2/class_selmer_15000_list.txt` (byte-identical copies of repository
  files) and the four per-worker files `g2/g2_tuple_selmer_997_{0,1}.txt`, `g2/g2_class_selmer_15000_list_{0,1}.txt`
  (sorted, their lines are `g2_997_all.txt`, `g2_15000_all.txt`; checked). No absolute path occurs in the bundle (grep);
  g3's pipeline calls census-wide scripts that still name the session's folders (pre-existing; recipe in
  experiments/README.md). No binary files; the largest are g4's three search records (2.5 MB each, needed by
  g4_summary.py). A NOTES.md per folder with the corrected text; REPORT_G*.md kept as written.
- `experiments/README.md`: a row in the table and a section "Roadmap G (committed in run 6)", tier D.
- `experiments/referee_x6/`: my scripts (print only), outputs, NOTES.md.
- `reviews/review_x6.md` (this review section as written) and a row in `reviews/README.md`.

## Proposed CI lines (job `descents`, after the run-5 block, before `git diff --exit-code`)

Run from a clean clone of run6_x6 (/tmp/x6clone), PARI 2.17.4 (cypari2), sympy 1.14.0, numpy 2.5.5 (both already
installed in that job); every line exit 0; then `git diff --exit-code` clean. Total **157 s** on a 2-core container.
A mutated g1_summary.txt, g4_summary.txt or g1_k13.txt makes its line fail.

```
          # run 6 (experiments/roadmap-g, tier D): the exact local test and its controls; two ranges of the stress test
          # (k = 20 holds the control (1109503,597671)); G2's tables, negative test and 2-adic witnesses; G3's table;
          # G4's bookkeeping and negative test; the review's own checks (about 2.5 minutes)
          cd "$GITHUB_WORKSPACE"
          R=experiments/roadmap-g; Q=experiments/referee_x6
          python3 $R/g1/test_exact_local.py 2>/dev/null | f5 | diff - $R/g1/test_exact_local.txt             # 3 s
          for k in 13 20; do python3 $R/g1/g1_stress.py $k 20000 200 2>/dev/null | f5 | diff - $R/g1/g1_k$k.txt; done   # 55 s
          python3 $R/g1/g1_stress.py --negative 2>/dev/null | f5 | diff - $R/g1/g1_negative.txt              # 1 s
          python3 $R/g1/g1_summary.py | diff - $R/g1/g1_summary.txt
          python3 $R/g2/g2_check.py 2>/dev/null | f5 | diff - $R/g2/g2_summary.txt
          python3 $R/g2/g2_check.py --negative 2>/dev/null | f5 | diff - $R/g2/g2_negative.txt               # 4 s
          for L in 997 15000; do python3 $R/g2/g2_features.py g2_${L}_all.txt 2>/dev/null | f5 | diff - $R/g2/g2_features_$L.txt; done
          python3 $R/g2/g2_badalone_witness.py $R/g2/g2_997_all.txt 2>/dev/null | f5 | diff - $R/g2/g2_badalone_witness_g2_997_all.txt   # 30 s
          python3 $R/g3/g3_table.py | f5 | diff - $R/g3/g3_table.txt
          python3 $R/g3/g3_table.py --negative | f5 | diff - $R/g3/g3_table_negative.txt
          python3 $R/g4/g4_summary.py 2>/dev/null | f5 | diff - $R/g4/g4_summary.txt
          python3 $R/g4/g4_search.py --negative 2>/dev/null | f5 | diff - $R/g4/g4_negative.txt            # 4 s
          python3 $R/g4/gen_square.py 20000 | tail -1 | grep -qx '# matches square_20000.json'               # 9 s
          python3 $Q/local_compare.py 2>/dev/null | f5 | diff - $Q/local_compare.txt                         # 40 s
          python3 $Q/weil_check.py 2>/dev/null | f5 | diff - $Q/weil_check.txt                               # 9 s
          python3 $Q/local_mutation.py 2>/dev/null | f5 | diff - $Q/local_mutation.txt                       # 7 s
          python3 $Q/g3_family.py | diff - $Q/g3_family.txt
          python3 $Q/g3_sample.py | diff - $Q/g3_sample.txt
          python3 $Q/samples/compare.py | diff - $Q/samples/compare.txt
```

(`f5` is the job's existing filter `grep -v -e 'PARI stack' -e '^# matches'`. The scripts' default modes also assert
equality with their committed files themselves.) Not for CI: `g1/run_g1.sh` (all ranges, ~1.5 h on two cores; k = 29
alone 13.5 min here), `g2/run_g2.sh`, `g2/run_g1kept.sh`, `g3/run_g3*.sh`, `g4/run_g4*.sh`, and the referee's
`kept63_repo_only.py` (~10 min), `control_1109503.py` (~14 min), `badalone_witness_all.py` (~30 min),
`killsets_witness.py` and `g4_enum.py` (outputs carry timings / PARI-dependent), `samples/samples.sh` (5 min).

**PARI 2.17.4 against the committed 2.19.1 outputs:** every line above reproduces the committed output byte for byte
(also k = 29). The only output that differs is G3's TSV, and not because of PARI (the F0 table was corrected after the
bundle was made; regenerated and committed). G4's search records (not in CI) differ on 3 of 5 rerun ratios in
PARI-chosen data (points, which curves get a full point set); they are left as the 2.19.1 records, since a rerun would
change coverage counts, and the bookkeeping CI checks is over the committed records.

## Proposed OUTPUTS.md rows (a new block "Added on 27 November 2026 (UTC), run 6")

| File | Class | Command | CI job | Time |
|---|---|---|---|---|
| `experiments/roadmap-g/g1/test_exact_local.txt` | A | `python3 experiments/roadmap-g/g1/test_exact_local.py` | descents | 3 s |
| `experiments/roadmap-g/g1/g1_k13.txt`, `g1_k20.txt` | A | `python3 experiments/roadmap-g/g1/g1_stress.py 13 20000 200` (and `20`) | descents | 55 s together |
| `experiments/roadmap-g/g1/g1_k14.txt` ... `g1_k19.txt`, `g1_k21.txt` ... `g1_k29.txt` | B | `bash experiments/roadmap-g/g1/run_g1.sh <k>` (k = 29 rerun by review x6: identical) | none | 15 s to 17 min each |
| `experiments/roadmap-g/g1/g1_negative.txt`, `g1_summary.txt` | A | `g1_stress.py --negative`, `g1_summary.py` | descents | 1 s |
| `experiments/roadmap-g/g1/err_k*.txt`, `times.txt` | C | logs of the session's runs | none | |
| `experiments/roadmap-g/g2/g2_summary.txt`, `g2_negative.txt`, `g2_features_997.txt`, `g2_features_15000.txt`, `g2_badalone_witness_g2_997_all.txt` | A | `g2_check.py` (and `--negative`), `g2_features.py g2_<L>_all.txt`, `g2_badalone_witness.py g2_997_all.txt` | descents | 35 s together |
| `experiments/roadmap-g/g2/g2_997_all.txt`, `g2_15000_all.txt`, `g2_15000_need_good.txt`, `g2_g1_kept_0.txt`, `g1_kept.txt`, `g2_badalone_witness_g2_15000_all.txt` | B | `bash run_g2.sh i N` (the 8,739: ~57 min per worker on two), `bash run_g1kept.sh` (8 GB per process), `g2_badalone_witness.py g2_15000_all.txt` | none | 3 min to 1 h |
| `experiments/roadmap-g/g3/g3_table.txt`, `CENSUS_TABLE_G3.tsv`, `g3_table_negative.txt` | A | `python3 experiments/roadmap-g/g3/g3_table.py` (and `--negative`) | descents | under 1 s |
| `experiments/roadmap-g/g3/` records (`sd_*`, `r0_*`, `dm_*`, `stage_*.json`, `fam_11000_30000.json`, `route/*`, logs) | B | `run_g3.sh`, `run_g3_stages.sh`, `route/run_route.sh` (census-wide paths as in `experiments/README.md`); PARI 2.19.1; samples rerun under 2.17.4 by review x6: identical | none | ~3 h per block |
| `experiments/roadmap-g/g4/g4_summary.txt`, `g4_negative.txt`, `square_20000.json` | A | `g4_summary.py`, `g4_search.py --negative`, `gen_square.py 20000` | descents | 14 s together |
| `experiments/roadmap-g/g4/g4_{0,1,2}_3.json`, `pass2_{0,1}_2.jsonl` | B | `bash run_g4.sh i 3`, `bash run_g4_pass2.sh i 2`; PARI 2.19.1, **not reproduced bit for bit under 2.17.4** (PARI's points) | none | ~4.6 h per worker |
| `experiments/referee_x6/local_compare.txt`, `weil_check.txt`, `local_mutation.txt`, `g3_family.txt`, `g3_sample.txt`, `samples/compare.txt` | A | `python3 experiments/referee_x6/<script>.py` | descents | about 1 minute together |
| `experiments/referee_x6/kept63_repo_only.txt`, `control_1109503.txt`, `badalone_witness_all.txt`, `killsets_witness.txt`, `g4_enum.txt`, `samples/*` (reruns) | B | `python3 experiments/referee_x6/<script>.py`; `bash samples/samples.sh` (outputs with timings, or PARI-dependent) | none | 1 to 30 min |

## Proposed COMPLETE.md section 7, test 6 (replacing the whole item 6)

> 6. *Is the four-class descent a law at the bad primes?* **No** (computed; run 6, `experiments/roadmap-g/g2/`, reviewed
> by a separate AI agent). With every bad prime -- 2, the primes of `mn(m² − n²)`, and every prime of `m² + n²` and of
> the octic `m⁸ + 14m⁴n⁴ + n⁸`, of any size -- the four-class descent leaves tuples for 101 of the 997 of test 4 and
> 836 of the 8,739 of test 7. The local test at the odd primes is `experiments/roadmap-g/g1/exact_local.py`, the test of
> `tuple_descent.py` computed without a depth cap (residue symbols on the residues of a disc, Weil's bound for their
> sign patterns, roots shared exactly in `ℤ[i]` closed in form); it is tested against the program's own test and
> against explicit local points, not proved. The answer does not rest on it: every tuple left has explicit local points
> at 2 and at every odd bad prime (`experiments/referee_x6/badalone_witness_all.txt`). The bad primes and the good
> primes up to 181, all `≡ 1 (mod 4)`, empty all 9,736, with irredundant kill sets (`g2_*_all.txt`; not proved of least
> size); no residue-symbol feature tried separates the two cases, and no conjecture is made. So a uniform statement
> would need the good primes as well, like Theorem C′, for the unramified `(ℤ/2)⁴`-cover of genus 257 of the curve of
> genus 17; this is not attempted. The earlier run (`tuple_badprimes_997.txt`, primes of `m² + n²` and of the octic only
> below 2000) agrees: of its 57, 39 keep tuples with every bad prime and 18 are emptied by the large ones. Experiment,
> not counted. Written by an AI model; not read by a person.

## Proposed LEDGER.md rows (tier D)

| Claim | Where | CI |
|---|---|---|
| **Random ratios to `m = 10⁹`: no survivor** among 340,000 coprime ratios (20,000 per dyadic range, `10⁴ ≤ m ≤ 10⁹`): each excluded by the class criterion (`COMPLETE.md` §7) or, for 63, by the four-class descent at the support primes and good primes `≤ 150` (also with the repository's local test alone). Of 3,400 degenerate controls, 3,399 keep exactly their degenerate tuples; (1109503,597671) keeps four more, **locally soluble everywhere** (explicit points at 2 and at every bad prime; `exact_local.py` and a short argument at the good primes); whether a rational point has one of these classes is open (run 6; a separate AI session; reviewed by a separate AI agent: holds with gaps, corrected in `NOTES.md`) | `experiments/roadmap-g/g1/` | `descents` (k = 13, 20; the exact local test; `referee_x6`) |
| **Test 6 of `COMPLETE.md` §7 answered no**: the bad primes alone leave tuples, each with explicit local points at every bad prime, for 101 of the 997 and 836 of the 8,739; bad primes and good primes `≤ 181` (all `≡ 1 mod 4`) empty all 9,736; no law found (reviewed: holds) | `experiments/roadmap-g/g2/` | `descents` (tables, negative test, 2-adic witnesses for the 997) |
| **Degenerate types, `11000 < m ≤ 13000`**: of 611 ratios, 4 open ((11165,7548), (12201,7360), (12204,9853), (12703,2767)), 3 conditional on a rank being 1 ((11441,1519), (11969,4529), (12007,1799)), the rest closed with the labels of `CENSUS_TABLE_G3.tsv` (PARI rank bounds; DM on PARI generators; (12416,5313) by the route, conditional on (a), (b)); other types not covered (reviewed: holds with gaps) | `experiments/roadmap-g/g3/` | `descents` (the table, bookkeeping) |
| Search (numerical, not a certificate): on the 3,573 square-type ratios `m ≤ 20,000`, the points of E± and 12 quotient curves up to 64 times the largest generator's height, on PARI's subgroup, give no solution and no non-degenerate `t` with four of the five forms squares (1,741 ratios with a full search on E+ or E−) | `experiments/roadmap-g/g4/` | `descents` (bookkeeping, negative test) |

## Proposed CHANGES.md entry (newest first, under the next version)

> - **Sixth run, Roadmap G** (27 November 2026, UTC; a separate AI session's bundle, reviewed by a separate hostile AI
>   agent, `reviews/review_x6.md`; written by an AI model, not read by a person; tier D; no census count changes).
>   - *Random ratios to `m = 10⁹`* (`experiments/roadmap-g/g1/`): no survivor among 340,000 random coprime ratios; a new
>     local test at odd primes without a depth cap (`exact_local.py`; tested, not proved); the degenerate control
>     (1109503,597671) keeps four classes besides its degenerate ones, locally soluble everywhere; whether a rational
>     point has them is open. Corrected after the review: `m² + n² = 2·(5³·7129)²` (the report had `5·7129`); "no prime
>     left out" holds for G2, not G1; "classes that no point explains" overstated.
>   - *Test 6 of `COMPLETE.md` §7 answered no* (`g2/`): the bad primes alone do not empty the four-class set for 101 of
>     the 997 and 836 of the 8,739 (explicit local points at every bad prime, found by the review); good primes `≤ 181`
>     finish every case; no law found. §7, item 6, rewritten.
>   - *Degenerate types, `11000 < m ≤ 13000`* (`g3/`): 4 open, 3 conditional on a rank being 1; one label regenerated
>     from the corrected run-5 table.
>   - *A point search on the square-type fibres to `m = 20,000`* (`g4/`): no solution among the points listed; a search,
>     not a certificate.
>   - *CI*: `descents` adds the exact local test, two ranges of the stress test, G2's, G3's and G4's checks, and the
>     review's own checks (about 2.5 minutes).

## Checkpoint (05:03 UTC)

Review items 1-6 written; package committed locally as a WIP commit (f4eda9b, to be amended into one commit); proposed
CI lines run from a clean clone (157 s, all pass, mutations refused); waiting for badalone_witness_all (937 ratios).
05:13: badalone_witness_all done: 937/937, 6,904 tuples, all witnessed at every bad prime.

## Final

Verdicts: (1) exact_local.py HOLDS (wording gaps: the "Weil above 5*10^6" description; tested, not proved; 107 cases
equal to the repository's test where it runs and to the labels of explicit points everywhere; Weil step and five
mutations checked). (2) G1 stress test HOLDS WITH GAPS: k = 13 and 29 byte for byte on PARI 2.17.4; all 63 kept ratios
emptied again with the repository's local test alone; false numbers/wording corrected in g1/NOTES.md (m^2 + n^2 =
2*(5^3*7129)^2, "no prime left out" false for G1, "classes that no point explains" overstated); the control's four
extra tuples are locally soluble everywhere (stronger than G1 said), their rational origin open. (3) G2 HOLDS: 20 kill
sets recomputed, irredundancy proved by witnesses; all 937 "bad primes fail" ratios (6,904 tuples) witnessed at every
bad prime; feature scan honest. (4) G3 HOLDS WITH GAPS: family and 21 rows checked, stages rerun identical; one stale
label (12007,1799) regenerated; "next step" misleading. (5) G4 HOLDS WITH GAPS: labelled a search; two sentences
overstated ("no solution" -> "none of the points listed"); planted test weak but enumeration checked; not reproducible
across PARI versions (bookkeeping only in CI). (6) Tier D throughout.

Commit: one commit on run6_x6 (hash in the final message). Not done: nothing of the brief is left out; not run: the
full G1 ranges other than 13, 20, 29; G2's kill-set runs (only checked by recomputation of 20 and the witnesses); G3's
second-descent blocks beyond the sampled rows; G4's searches beyond 5 ratios. COMPLETE.md, LEDGER.md, CHANGES.md,
OUTPUTS.md, ci.yml and MANIFEST.sha256 not edited (texts proposed above).
