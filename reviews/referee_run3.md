# Referee report, run 3 (diff f03f8af..80d2b9f: dcd6b56, 32ebb20, 66c590c, 80d2b9f)

VERDICT: PUBLISH AFTER FIXES

No mathematical error found. `E1Points` (hence `LemmaD`) is TRUE: PARI over Q(sqrt2) gives torsion exactly Z/10 (generator (-12-8sqrt2, -32-24sqrt2)), conductor (3) (I1), L(E1/K,1) = 0.22396 != 0; E1 is the twist by d of E0: y^2 = x(x^2+8x-8)/Q, so rank 0 follows from Kato (chi-part of E0 over Q(zeta16)^+) or Zhang GZK (E1/K modular, conductor (3) not a square); LMFDB 576.a.576.1 lists J(Q) = Z/10, MW rank 0. The reduction `E1Points -> LemmaD` is correct, denominators handled, hypothesis neither vacuous nor over-strong. Defects are prose overstatement/staleness, two dangling paths, and two CI checks that are weaker than described.

## Ranked findings

### FIX-FIRST
1. README.md:171 (`lean_rank/` row) and LEDGER.md:35 (row title): "Theorems 1 and 2 for every coprime pair with Lemma D as the only hypothesis" overstates Theorem 2 (Lean `rank_family_all` is for coprime pairs with m^2+n^2 = 2k^2 only) and drops mn(m^2-n^2) != 0 (`rank_one_of_lemmaD`, `R_infinite_order` need m != 0, n != 0, m^2 != n^2). Replace with: "Lemma B, and, with Lemma D as the only hypothesis, Theorem 1 at every coprime pair with `mn(m² − n²) ≠ 0` and Theorem 2 at every coprime pair with `m² + n²` twice a square (on the models `Eplus m (±n)`)".
2. UNIFORM_RANK.md:85 (*Lean.* paragraph): three old sentences now contradict the new text of the same paragraph (Lemma B is in Lean). Replace "Lemmas B and D, the descent over `ℚ(√2)` and the local tests are not in Lean." by "Lemma D, the descent over `ℚ(√2)` and the local tests are not in Lean (Lemma B is: `lean_rank/LeanRank/LemmaB.lean`)."; "Not in Lean: that the list of 127 is complete, that `Eplus m (−n)` is `E⁻` (sympy), and Lemmas B and D, which make the statement uniform." by "Not in Lean: that the list of 127 is complete, that `Eplus m (−n)` is `E⁻` (sympy), and Lemma D."; "(Lemmas B and D in Lean)" by "(Lemma D, i.e. `E1Points`, in Lean)".

### FIX
3. Dangling paths (no `reports/` in the repo): lean_rank/LeanRank/TwoIsogenyMaps.lean:8 "(`reports/v8.md`, "How far, what next")" -> "(a plan kept outside this repository)"; scripts/g1_full/review_t8_rank1.py:14 "(re-derived; see reports/review_t8.md)" -> "(re-derived; the review's notes are not committed)". Update MANIFEST.
4. .github/workflows/ci.yml:99: `len(lines) < 76` accepts >= 76, so k missing lines are masked by k new ones. Require equality with a count computed from the sources (`sum(open(f).read().count('\n#print axioms') for f in glob.glob('LeanRank/*.lean'))`), better: equality of the sets of names.
5. ci.yml:104-108: FamilyForge/LemmaDForge/LemmaDBridgeForge are checked by error count only (RankForge alone has line pins). Pin the error lines for each file as for RankForge, so a forgery accepted plus a broken control cannot keep the count.
6. lean_rank/LeanRank/RankZeroCensus.lean:10 and gen_rank_zero_certs.py:135 (template): "`Y = m³n³y`" is wrong; X(X+m^4)(X+n^4) = m^4n^4 * w(m^2+n^2w)(n^2+m^2w), so "`Y = m²n²y`" (change both, or the CI regeneration diff fails).

### NIT
- ci.yml:92 `lake build > build.txt 2>&1 || { cat build.txt; exit 1; }`; also fail on "declaration uses 'sorry'".
- Headline rank-0 theorems (`RankZero.finrank_eq_zero`, `finrank_4_3`...) state only `finrank = 0`; f.g. is proved (instance / inside `finrank_eq_zero_of`), but put `Module.Finite`/`IsAddFG` in the statement so the type alone means rank 0.
- `R_infinite_order` takes `hR` as hypothesis; restate with `∃ hR` like `P_infinite_order`.
- CHANGES.md:10 "all 44 modules": this tree has 42 (+ root); say which.
- CHANGES.md (third run, positioning) "the forms ... are not in Bremner (1999)" -> "were not found in Bremner (1999) by an AI reading of a text extraction" (as G1_FULL_SYSTEM says).
- CHANGES.md:5 "the statements are closed" -> "the statements have no extra hypotheses".
- SQUARE_SECTION.md §1 "the sixth was not run": PARI ellrank gives [1,1] on E+(180,19) too (mine, and review_t8_rank1.txt already has it).
- SECOND_DESCENT.md:137 "(root number +1 by PARI)": not in any committed output.
- SELMER_SQUARE.md:3 add "(over ℚ there is no constant 2-Selmer parity)" after "take every 2-Selmer rank".
- REPRODUCE.md:55 "the 2-adic classes of the rank-0 certificates" -> "of the certificate of `RankZero.lean`" (rank_zero_check.py covers only that one).
- LEDGER.md:36 lists `TwoIsogenyMaps` as a source of the reduction; it is not used by it.
- LemmaDReduction.lean:23 vs :64 "order dividing 10" / "order 10": harmonise.
- ci.yml:75 rank-job comment describes only rank >= 2.
- quick temp files: use `mktemp -d` instead of the repo root (no false green possible now).

## Raw log (appended during the run)

### P1. LemmaDReduction / LemmaDBridge (plain maths and truth of the hypothesis)

- `LemmaD` (Theorem2Family.lean:35): every (x,w) in Q^2 with w^2 = 2(x^2+1)(x^4-2x^3+2x^2+2x+1) has x = +-1. Affine only (D2 has leading coeff 2, a non-square, so no rational points at infinity; nothing lost).
- `E1Q2 a b c e` (LemmaDReduction.lean:60): (a+b*sqrt2, c+e*sqrt2), a,b,c,e in Q, satisfies Y^2 = X(X^2+8dX-8d^2) in R, d = 2+sqrt2. Since every element of Q(sqrt2) is a+b*sqrt2, this is exactly "affine point of E1(Q(sqrt2))". Not vacuous: `E1_points_on` exhibits 9 points.
- `E1Points` (:66): every affine point of E1(Q(sqrt2)) has X in {0, +-4, +-(12+8sqrt2)}. Given `E1_points_on` this is equivalent to E1(Q(sqrt2)) = the cyclic group of order 10.
- `lemmaD_of_E1Points : E1Points -> LemmaD` (:129). Proof: (x,w) |-> pullback point with explicit rational coordinates (a,b,c,e) (`pullback`, :109); denominators ((x-1)^2-2) and (x-1+sqrt2) are proved nonzero for rational x (`rat_sq_ne_two`, irrationality of sqrt2), so no x/0=0 artefacts. Case split on the 5 abscissas: X=0 impossible (a>0); X=4 forces x=1; X=-4 forces x=1 then a=4, contradiction; X=12+8sqrt2: b=8 gives x(x+1)(x^2-4x+5)=0 (I re-derived: ((x-1)^2-2)^2+(x-1)(x^2+1) = x^4-3x^3+x^2+5x = x(x+1)(x^2-4x+5)), x=0 needs w^2=2; X=-12-8sqrt2: no common solution (linear_combination certificate, kernel-checked). Correct.
- `LemmaDBridge`: `E1PointsK` (Mathlib `E1K.Equation` over `Q2 = Q(sqrt2)` as IntermediateField of R) -> `E1Points` -> `LemmaD`. Direction is the useful one (Mathlib form implies the R form). `E1K_fg` (Stoll fg) is just finite generation, not rank. Correct as stated.
- d2_to_E1 identity: checked by `lemmaD_reduction_check.py` (prints ALL CHECKS PASSED, matches .txt) at the 4 rational points only; the Lean proof is the general identity. Fine.
- **Is E1Points true? Yes (independent check, PARI 2.17.2 via cypari2, not in repo):**
  - `elltors(E1/K)` = [10,[10],[(-12-8t, -32-24t)]] : torsion exactly Z/10, generator matches script's `gen`. (Also forced by Kamienny-Kenku-Momose + `E1_two_torsion`: Z/2xZ/10 excluded by one 2-torsion point; Z/20, Z/30 do not occur over quadratic fields.)
  - Conductor of E1/K = (3) (3 inert), reduction I1 (multiplicative) at 3, good elsewhere; j = 85184/3.
  - `lfun(E1/K,1)` = 0.223962521681858, `lfunorderzero` = 0, L-conductor 576 (= 8^2*9; matches LMFDB 576.a for D2; `genus2red` of D2 gives odd conductor part 9).
  - E1 is the quadratic twist by d of E0: y^2 = x(x^2+8x-8) over Q (X = d x, Y = d^{3/2} y). So E1(K) (x) C = chi + chi-bar parts of E0(L), L = K(sqrt d) = Q(zeta16)^+, chi the even quartic characters of conductor 16, and L(E1/K,s) = L(E0,chi,s)L(E0,chi-bar,s). Rank 0 then follows from Kato (L(E0,chi,1) != 0 => chi-part of E0(L) finite), or from Zhang's GZK over totally real fields (E1/K modular by Freitas-Le Hung-Siksek; conductor (3) is not a square). E0(Q) rank 1, E0^(2)(Q) rank 0 (PARI ellrank), consistent with rank E0(K)=1.
  - So E1Points holds, LemmaD holds; neither is proved in Lean. Hypothesis is neither vacuous nor over-strong: it is exactly the true statement E1(Q(sqrt2)) = Z/10.
  - Also `hyperellratpoints(D2, 10^4)` = only (+-1, +-4).
- Prose in the module docstring (LemmaDReduction.lean:20-23) "rank 0 (C2) and torsion of order dividing 10 (C3)": OK, correctly flagged "not proved here".
- NIT: LemmaDReduction.lean:63-65 docstring says "torsion of order 10" but :23 says "order dividing 10"; both correct in context (with E1_points_on), harmonise.


### P3. Statements of Theorem1Family / Theorem1Order / Theorem2Family / LemmaB / RankZero* / TwoIsogenyMaps

- `rank_one_of_lemmaD`, `rank_one_family` (Theorem1Family.lean:136,191): hypotheses LemmaD, IsCoprime m n, m != 0, n != 0, m^2 != n^2. No sign/order condition. Matches "every coprime pair" only with mn(m^2-n^2) != 0; prose (LEDGER new row, UNIFORM_RANK *Lean.*) drops that. NIT (the excluded pairs are the degenerate (0,1),(1,0),(1,+-1)).
- `R_infinite_order` (Theorem1Order.lean:37) takes the nonsingularity proof `hR` as a hypothesis (satisfiable: `R_add_R` supplies it, and `P_infinite_order` packages it with an exists). Not vacuous. NIT: restate with `∃ hR` like `P_infinite_order`.
- `rank_family_all` (LemmaDEquiv.lean:183): coprime, m^2+n^2 = 2k^2, m^2 != n^2, any signs/order. Good. `rank_family`/`theorem2_family` (Theorem2Family.lean) still need 1 <= n < m; prose that says "every coprime twice-square pair" is backed by `rank_family_all`, fine.
- finrank lower bounds (`1 <=`, `2 <=`) cannot be vacuous: Module.finrank = 0 for non-finite modules, so a lower bound >= 1 is only stronger. OK.
- `RankZero.finrank_eq_zero` (RankZero.lean:604), `finrank_4_3` etc. (RankZeroCensus.lean:74,112,152,191): `finrank = 0` alone would be vacuous for a non-f.g. module, but f.g. is PROVED (`instance : IsAddFG E.Point := fg_point_of_numberField`, RankZero.lean:601; and inside `finrank_eq_zero_of`, which uses Stoll's `card_range_μ`, which needs the FG instance). And `RankZeroPoints.card_E` (#E'(Q) = 8) uses rank 0 + reduction mod 5, so the rank-0 claim is used non-vacuously there. OK. NIT: state the headline as `Module.Finite ℤ E.Point ∧ finrank ℤ E.Point = 0` (or `Finite E.Point`) so the type alone says "rank 0".
- Independent PARI check: y^2 = x(x+m^4)(x+n^4) has rank [0,0,0] at (4,3),(7,1),(12,5),(24,7); rank 1 at (15,8) (so the generator's refusal of (15,8) is correct); y^2 = x(x-12)(x-16) rank 0, torsion Z/4 x Z/2 (8 points), matches `card_E`.
- NIT (math error in docstring): RankZeroCensus.lean:10 "with `X = m²n²w`, `Y = m³n³y`" is wrong: X(X+m^4)(X+n^4) = m^4 n^4 * w(m^2+n^2 w)(n^2+m^2 w) (sympy), so `Y = m²n²y`. (Y = m^3n^3 y lands on the twist by m^2n^2, isomorphic, so no consequence for rank.) Fix: replace `Y = m³n³y` by `Y = m²n²y` here and in the generator's template (see P5).
- LemmaB.lean: `lemmaB` (:306) coprime, mn != 0, m^2 != n^2 -> not square in Z; `lemmaB_rat` in Q. Matches prose. Negative control `control_lemmaB_false` shows m^2 != n^2 is needed.
- TwoIsogenyMaps: field identities only; docstring says so ("that last identification is not made here"). Listed in LEDGER row "Lemma D reduced ..." as a source module though it contributes nothing to the reduction. NIT.

### P4. CI (.github/workflows/ci.yml)

- Axiom count: `grep -c '^#print axioms' lean_rank/LeanRank/*.lean` = 76 exactly (LemmaB 4, LemmaDEquiv 5, Rank 5, RankZeroCensus 4, RankZeroGenCurve 3, Theorem1Family 5, Theorem2Family 9, Index 4, LemmaDBridge 4, LemmaDReduction 6, RankZero 6, RankZeroGen 5, RankZeroPoints 9, Theorem1Order 3, TwoIsogenyMaps 4); the four main-tree roots (G1Sigma, G1SigmaGeneric, TwoDescentMap, TwoDescentCertificates) have 0. REPRODUCE "76 lines" correct.
- The python heredoc fails closed: 0 lines (no log replay / quiet build) -> `len(lines) < 76` -> exit 1; a wrapped/unparsable axiom list -> `re.search(...)` is None -> AttributeError -> non-zero; sorryAx anywhere -> exit 1. It cannot pass vacuously. It relies on Lake replaying cached info logs on the warm build (lean-action has already run `lake build`); if replay ever stops, CI goes red, not green.
- FIX ci.yml:96-97: the test is `len(lines) < 76`, i.e. ">= 76", not "= 76". A future file adding k new `#print axioms` lines while k old ones silently stop being emitted (e.g. a module dropped from the build) still passes. Fix: compute the expected count from the sources and require equality, e.g. add `exp = sum(open(f).read().count('\n#print axioms') for f in glob.glob('LeanRank/*.lean'))` and `len(lines) != exp`. Also check names: the set of names in build.txt equals the set of names in the sources.
- NIT ci.yml:94: `lake build > build.txt 2>&1` under `set -e` hides the build log on failure. Use `lake build > build.txt 2>&1 || { cat build.txt; exit 1; }`.
- NIT: the `sorryAx` check only sees theorems that are #print-axiomed; also fail on `declaration uses 'sorry'` in build.txt.
- `test "$(grep -c ... $f.txt)" = "$n"`: grep -c printing 0 exits 1, but the substitution's status is ignored inside a `test` argument (set -e looks at `test`'s status), so count 0 gives `test 0 = 4` -> false -> step fails. Correct.
- FIX ci.yml:101-106: for FamilyForge/LemmaDForge/LemmaDBridgeForge only the total error count is checked; RankForge alone has line checks. A forged example that is wrongly accepted, compensated by an ACCEPTED control that breaks (or a REFUSED example emitting 2 errors), keeps the count. Fix: pin the error line numbers for each file as for RankForge (e.g. a `lines` map `FamilyForge:"19 23 27 33 ..."`), or split each forge file into one-example files.
- quick temp files (rzc.txt, rz.lean, rzgen.txt, rzbad.lean, rzbad.txt) are written in the repo root after `verify_manifest.py` ran (separate earlier step) and are removed before the final `git diff --exit-code`, which ignores untracked files anyway. No interaction. NIT: write them under `$RUNNER_TEMP`/`mktemp -d` (if a step between creation and `rm` fails the job fails anyway, so no false green).
- NIT: `rank_zero_check.py` is checked only by its last line; no committed transcript to diff.
- `klagsbrun_check.py` default (no --pari) output equals head -18 + tail -1 of the committed .txt (I ran it: CI diff OK, 2.9 s). The PARI lines 19-46 are guide only, not checked in CI; prose says so.

### P2. Prose against Lean/scripts

- FIX UNIFORM_RANK.md:85 (the *Lean.* paragraph, unchanged sentences now contradicted by the new text in the same paragraph): "Lemmas B and D, the descent over `ℚ(√2)` and the local tests are not in Lean." and "Not in Lean: that the list of 127 is complete, that `Eplus m (−n)` is `E⁻` (sympy), and Lemmas B and D, which make the statement uniform." and the last sentence "a Lean proof of Theorem 2 for every ratio (Lemmas B and D in Lean) would be the first family-uniform statement ...". Lemma B is now in Lean (`LemmaB.lemmaB`). Replace with: "Lemma D, the descent over `ℚ(√2)` and the local tests are not in Lean (Lemma B is, since 25 November: `lean_rank/LeanRank/LemmaB.lean`)." / "Not in Lean: that the list of 127 is complete, that `Eplus m (−n)` is `E⁻` (sympy), and Lemma D." / "an unconditional Lean proof of Theorem 2 for every ratio (Lemma D, i.e. `E1Points`, in Lean) would be ...".
- FIX README.md:171 (lean_rank row): "Lemma B and Theorems 1 and 2 for every coprime pair with Lemma D as the only hypothesis" overstates Theorem 2 (Lean: twice-square pairs only, `rank_family_all`) and drops mn(m^2-n^2) != 0. Replace by: "Lemma B, and, with Lemma D as the only hypothesis, Theorem 1 at every coprime pair with `mn(m² − n²) ≠ 0` and Theorem 2 at every coprime pair with `m² + n²` twice a square (on the models `Eplus m (±n)`)". Same fix for the LEDGER.md new row title ("for every coprime pair" -> "(Theorem 1 at every coprime pair with mn(m² − n²) ≠ 0, Theorem 2 at every coprime twice-square pair)").
- FIX dangling references (files not in the repository; `git ls-files | grep ^reports/` is empty): lean_rank/LeanRank/TwoIsogenyMaps.lean:8 "`reports/v8.md`, "How far, what next"" and scripts/g1_full/review_t8_rank1.py:14 "see reports/review_t8.md". Fix: delete the parenthesis, or replace with "(a plan outside this repository)" / "(the review's notes, not committed)". (TwoIsogenyMaps.lean is a source file and is in MANIFEST; changing it changes the hash.)
- NIT RankZeroCensus.lean:10 and gen_rank_zero_certs.py:135: `Y = m³n³y` -> `Y = m²n²y` (see P3). The generator template must change with it, else the CI regeneration diff fails.
- NIT CHANGES.md:9 "the main library on Lean v4.36.0 (all 44 modules build ...)": this tree has 42 modules in MagicSquareOfSquares/ (+ root = 43), and UNIFORM_RANK.md:85 / CHANGES second run say 42. Fix: give the count of this tree ("all 42 modules") or say what the two extra modules are.
- NIT SQUARE_SECTION.md §1 "(PARI's `ellrank`, guide only, also gives `[1, 1]` on the first five; the sixth was not run.)": I ran it: `ellrank` gives [1,1] on all six, including E+(180,19) (models `sigma_exact.model(m, ±n, 2)`; j of (m,-n,2) equals j of (m,n,-2) on all six). Optional update: "on all six".
- NIT SECOND_DESCENT.md:137 "(root number +1 by PARI)": no committed output contains it (review_t5_labels.txt does not). Either commit the PARI line or write "(root number +1 by PARI, not committed)".
- NIT SELMER_SQUARE.md:3 "take every 2-Selmer rank": Klagsbrun's Thm 1.1 says "subject only to constant 2-Selmer parity" (arXiv abstract); over Q constant parity never occurs, so correct, but add "(over ℚ there is no constant 2-Selmer parity)".
- NIT REPRODUCE.md:55 "the 2-adic classes of the rank-0 certificates": `rank_zero_check.py` cross-checks only the certificate of `RankZero.lean` (y^2 = x(x-12)(x-16)); the census certificates are only regenerated. Fix: "the 2-adic classes of the certificate of `RankZero.lean`".
- NIT CHANGES.md:7 "the statements are closed": ambiguous next to census "closures"; write "the statements have no extra hypotheses".
- Counts verified: 76 axiom lines (exact); "four census curves" = (4,3),(7,1),(12,5),(24,7) (RankZeroCensus has 4 `finrank_*`); "254 curves" unchanged claim; square_section.tsv: 108 curves, class2 status sha 37 / point 3 / open 68, alphaE=K on 21, rank exactly 1 on 6 (matches text and table row sums: 0+2+1=3, 0+5+19+11+2=37, 8+19+24+11+5+1=68); klagsbrun "179 square-type ratios, 358 curves; 127 twice-type, 254 curves" printed by the script; review_t5_labels.py output = committed .txt; square_section.py --negative = committed .txt; lemmaD_reduction_check.py = committed .txt; generator regeneration byte-identical to RankZeroCensus.lean and (15,8) refused with exit 1.
- All other repository paths named in the new prose exist.
- LEDGER "Lemma D reduced ..." row cites `TwoIsogenyMaps` as a source; it contributes nothing to the reduction (field identities for a future descent). NIT: move to a parenthesis "(and `TwoIsogenyMaps`, identities for a future descent)".

### Added range: 66c590c, 80d2b9f (review range now f03f8af..80d2b9f)

- 66c590c pt_tors.py: every line now printed by Python (PARI builds strings with `Str`, returned as a vector). I ran the CI command `python3 pt_tors.py 2>/dev/null | grep -v "PARI stack" | diff - pt_tors.txt` 3x buffered and once with PYTHONUNBUFFERED=1: identical each time; pt_tors.txt unchanged by the commit (byte-compare with 66c590c^). The "PARI stack size set to ..." line comes on stdout from PARI's C buffer (its position can still vary) but is removed by `grep -v` before the diff. OK.
- Other CI-diffed scripts: static search of all tracked *.py using cypari for GP-side `print/print1/printf/printsep` inside strings passed to `pari(...)` or triple-quoted GP blocks, and for `pari.print*`: no hit apart from the old pt_tors (one false positive, `second_descent/rank_effort.py`, is a Python print in a subprocess and is not in CI). Remaining mixing risk is only the PARI stack banner, which every CI line filters with `grep -v "PARI stack"` where cypari is used... NIT: `lean_rank`-unrelated, but any new cypari script diffed in CI must keep that filter.
- 80d2b9f: "from source" removed for Mathlib in ci.yml comment, LEDGER row and REPRODUCE row; `grep -rn "from source"` now finds only REPRODUCE.md:33 (the manual alternative), correct. Timeout 60 -> 150. NIT REPRODUCE.md:54 "3 minutes at run 2" is unsourced (no log committed); fine as an estimate. NIT ci.yml:75 job comment still describes the job as "rank E(Q) >= 2 on the 254 twice-square curves" only; it now also builds Lemma B, Theorems 1/2 up to Lemma D, the Lemma D reduction and the rank-0 files.

### Corrections to the raw log above
- CI line numbers in P4 were approximate; at 80d2b9f: `lake build > build.txt` is ci.yml:92, the `len(lines) < 76` test ci.yml:99, the forge loop ci.yml:104-108 (the ranked list at the top uses these).
- P5 spot check: `review_t8_rank1.py` (default, six curves) rerun here: E+(15,8) and E-(15,8) reproduce "every coset excluded ... RANK EXACTLY 1: True"; the run on (56,33) was still going when I stopped reading (not waited for). PARI ellrank [1,1] on all six independently.
