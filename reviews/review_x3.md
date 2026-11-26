# x3 — hostile review of Fable run 5 part 1 (F4, F2, F1, wording)

Written by an AI model; not read by a person. Worktree /home/claude/conductor/x3 (branch run5_x3, tree a15ce41 = 24d6ef6 tree).
Started 20:55 UTC.

## 21:10 UTC — safety read and first CI run (PARI 2.17.4, shim as shipped)

Safety: read every .py in F2/ and F4/ (g4k_exact, dm_g4k, g4k_rank_exact, sympy_pari_shim, run_shim, f2_curves,
f2_general, f2_general_neg, f2_pairing, f2_d1_local2). No network, no subprocess, no deletes; every open(...,'w') is
behind `--write` (f2_curves: behind `--find`). g4k_rank_exact/f2_* do os.chdir into scripts/g1_full (read only).
Fallback REPO path in g4k_exact/g4k_rank_exact points at the author's /tmp scratch dir (harmless; dead in repo layout).
Copied F4 g4k/ and F2 rank1/ into the worktree; ran with `python3 -I` from each dir.

Bundle CI lines, as shipped (shim), PARI 2.17.4, 2 cores idle:
| line | result | wall |
|---|---|---|
| g4k_exact.py | PASS | 4 s |
| g4k_exact.py --negative | (see below) | |
| g4k_rank_exact.py | FAIL: shim `intervals()` AssertionError 'could not isolate a root' | 1.5 s |
| g4k_rank_exact.py --negative | FAIL, same | 1.3 s |
| dm_g4k.py | PASS | 11 s |
| dm_g4k.py --negative | PASS | 31 s |
| run_shim.py rank_indep.py | FAIL, same shim assertion (after 8 lines) | 1.2 s |
| f2_general.py 3 | FAIL, same | 3 s |
| f2_general_neg.py | FAIL, same (part 1 passes, part 2 dies) | 3 s |
| f2_d1_local2.py | FAIL, same | 3 s |
| f2_d1_local2.py --negative | FAIL, same | 3 s |

Cause (checked): under this cypari2/PARI 2.17.4, `pari.set_real_precision(250)` does not raise the working precision
of `polrootsreal` (roots come back with precision 19 digits), so the shim's Sturm-certified isolation loop cannot
converge on the large-coefficient polynomials and asserts. The shim fails CLOSED (assertion, no output), so it did not
produce a wrong result; but "the shim reproduces committed outputs byte for byte" is environment-dependent and false
here. Every shim-dependent line must be re-run with real sympy (below).

## 21:25 UTC — real-sympy re-run (the one-line change, made by a wrapper /tmp/x3/rs.py that maps
`sympy_pari_shim` to real sympy 1.14; scripts unchanged)

| line | shim (as shipped) | real sympy |
|---|---|---|
| g4k_exact.py | PASS 4 s | PASS 4.5 s |
| g4k_exact.py --negative | PASS 70 s | PASS 72 s |
| g4k_rank_exact.py | FAIL (shim) | PASS 24 s |
| g4k_rank_exact.py --negative | FAIL (shim) | PASS 9 s |
| dm_g4k.py | PASS 11 s | (no sympy) |
| dm_g4k.py --negative | PASS 31 s | (no sympy) |
| rank_indep.py vs committed rank_indep.txt | FAIL via run_shim | PASS 4 s (run directly, real sympy) |
| f2_general.py 3 | FAIL (shim) | PASS 14 s |
| f2_general_neg.py | FAIL (shim) | PASS 13 s |
| f2_d1_local2.py | FAIL (shim) | PASS 33 s |
| f2_d1_local2.py --negative | FAIL (shim) | PASS 32 s |

So: all eleven outputs are reproduced byte for byte with real sympy + PARI 2.17.4. With the shim, 7 of 11 fail closed
here. The shim cannot have changed an outcome in these runs: real sympy gives the identical text. Recommendation: do
NOT commit sympy_pari_shim.py or route CI through it; replace each `import sympy_pari_shim; sys.modules['sympy'] = ...`
line by nothing (real sympy is a repository dependency already). REPORT_F2's CI line `run_shim.py
../../scripts/g1_full/rank_indep.py` also has a wrong relative path (HANDOFF has ../../../).

## F4 — independent re-derivation (my code: /tmp/x3/mine/dm_indep.py, index_indep.py; exact rationals)

Algebra (checked as polynomial identities at enough rational points, and by hand):
- global.py: g4 = [2S, 2S], K = [S², 2S²−16m²n², S²] (low degree first), S = m²+n², w = t².
  g4(t²) = 4F4 and K(t²) = 4LU, so g4K(t²) = 16·F4·L·U: a solution (F4, L, U squares) gives (w, Y) = (t², 4√(F4LU)) on
  C: Y² = g4K(w). Only three of the five forms are used; that is fine (weaker hypothesis).
- c0 = c3 = 2S³, c1 = c2 = 2S(3S² − 16m²n²): palindromic, so w³·g4K(1/w) = g4K(w); for w = t² ≠ 0, (1/w, Y/t³) is on C.
  Model x = c3·w, y = c3·Y: y² = x³ + c2x² + c1c3x + c0c3² (= dm_g4k's ainv). x·x′ = c3²; t = ±√(x/c3). HOLDS.
- Rational 2-torsion: K has discriminant −64m²n²(m²−n²)² < 0, so the only rational root is w = −1, x = e = −c3.
  T = (−c3, 0) has w = −1 (not t²). So P = T, P′ = T are impossible; P = O impossible (x finite); P′ = O iff t = 0
  (tested directly, degenerate: L = U). t = ∞ is the image of t = 0 under t → 1/t. HOLDS.
- Dependence: DM here needs only (i) |ĥ(P) − ĥ(P′)| bounded, (ii) the relation not identically true on P′ = ±P + T″.
  (ii) is the explicit quadratic e x² + (β − e² − c3²) x + c3² e = 0 (β = f′(e); I checked x(P+T) = e + β/(x−e) on G
  and 2G) plus x = ±c3 — neither identically zero. Independence of the two maps D → C (bundle's J(D) ~ C×C remark)
  is true but not needed. HOLDS.
- |a| ≠ |b| ⇒ 2·max(|a|,|b|) − 1 ≤ |a²−b²| (B ≤ A−1 ⇒ A²−B² ≥ 2A−1). HOLDS.
- Height relation: x = p/q ⇒ x′ = c3²q/p, so |h(x) − h(x′)| ≤ 2 log|c3|. B, U (Silverman 1990 Thm 1.1, doubled:
  B = h(Δ)/6 + log⁺|j|/6 + log⁺|b2/12| + log 2* + h(j)/12 + 1.946, U = same 2µ + 2.14) are computed on Ed, the
  NON-minimal integral model y² = x³ + c2x² + c1c3x + c0c3² — the model the relation is written in. Correct: Silverman's
  bounds hold for any integral model. HOLDS.

Independent numbers (k = 8 doublings, own ln enclosures; bundle used k = 7):
| | ĥ(G) mine | bundle | DM rhs B+U+2log c3 | amax | #cand x | t not "no" |
|---|---|---|---|---|---|---|
| (7735,4968) | [9.44404, 9.44763] | [9.43923, 9.45355] | 345.592 (=) | 18 (=) | 37 (=) | ±1 degenerate (=) |
| (6557,324) | [8.83130, 8.83480] | [8.82673, 8.84072] | 336.046 (=) | 19 (=) | 39 (=) | ±1 degenerate (=) |
| (6881,6233) | [4.48021, 4.48385] | [4.47447, 4.48900] | 349.143 (=) | 39 (=) | 80 (=) | 0 degenerate (=) |
Robustness: with amax tripled (54/57/117) still only these t. No SOLUTION anywhere.

Index (G generates mod torsion, given rank 1), own code:
- Torsion: odd part 1 (gcd of #C(F_q), q < 600; the 2-part of the gcd is 4, which is why the bundle needs its
  three-step argument), one rational root, T ∉ 2C(Q) (b square? then a ± 2√b not squares). Z/2. HOLDS.
- Saturation at all 46 primes ≤ 199: odd l by (N/l)G ≠ O mod q; l = 2 by a DIFFERENT criterion from the bundle
  (x(P) − e_i all squares in F_q, with the norm reduction when the quadratic is irreducible; cross-checked against the
  order criterion). Found for all three curves (largest q 5779/3229/4391). Note: q with irreducible quadratic never
  work for l = 2 here (b·disc behaves like a square), so the bundle's brute-force over all q is what is needed. HOLDS.
- Index bound: k = [C(Q):⟨G⟩+tors] satisfies k² ≤ ĥ(G)/λ ≤ hi/(hi/210²) ⇒ k ≤ 210 ⇒ every prime factor ≤ 199
  (no prime in 200..210). Complete, with λ from the search. HOLDS.
- CPS search: PARI's minimal model checked (my own Silverman III.1 transcription; Δ = u¹²Δmin; integral); m =
  lcm(12, n_p at multiplicative p) = 24 for all three (Kodaira at 2: I3*, I8* — component groups of order 4, exponent
  | 12). Minimality at 2 is non-trivial for (6557,324) and (6881,6233) (v(Δ)≥12, v(c4)≥4, v(c6)≥6): confirmed by my
  own search over [2, r mod 4, s mod 2, t mod 8] (positive control: the model scaled by 2 is refused) and by PARI
  elllocalred (u = 1). Search bounds X0 = 14, 50, 129 (=), no point found (=). HOLDS (given CPS as quoted, tier E).

F4 DM step: HOLDS. Every step after "rank C(Q) = 1" re-derived and recomputed independently with identical results.

## F4 rank — tier (21:15 UTC)

- Exact bound (rank_indep through real sympy, my re-run): 2, 3, 5 — reproduced (g4k_rank_exact.txt byte for byte).
  The bundle's wrapper around find_conic_point only changes a hint; rank_indep checks every point. OK.
- PARI 2.17.4 ellrank over the isogeny class (effort 0), [r,R,s] per isogenous curve:
  (7735,4968): R = 1 on all 4 curves. (6557,324): R = 1 on ONE of 4 curves, R = 3 on the other three.
  (6881,6233): R = 1 on ONE of 8 curves, R = 3 or 5 on the others. Root number −1 for all three.
- Plain statement: the three closures are CONDITIONAL — on rank C(Q) = 1, which nothing in the repository proves.
  The only support is PARI's ellrank upper bound: an unverified program output whose Selmer computation goes through
  bnfinit of the 2-division algebra (PARI developers, pari-dev Feb 2021), i.e. GRH-dependent class-group/unit data unless
  bnfcertify is run (it was not), plus PARI's own Cassels–Tate implementation; for (6557,324) and (6881,6233) it rests on
  a single isogenous model. This is the same tier as the census's existing "conditional on rank 1" ratios, not a proof,
  and not stronger than the numerical-L-value closure it would sit beside for (6557,324)/(6881,6233) (different
  condition, not a better one).
- Possible upgrade for (7735,4968) only (NOT done, not checked; for the conductor/Drew): exact bound 2, a point of
  infinite order, root number −1. If the second-descent bound of rank_indep dominates the 2^∞-Selmer corank r + δ
  (elements from the divisible part of Ш survive every second descent; on the two sides they contribute
  dim D[φ] + dim D′[φ̂] = δ), then the 2-parity theorem (Monsky 1996; Dokchitser–Dokchitser) gives r + δ odd ≤ 2, so
  r = 1 unconditionally (modulo published theorems and an exact root number). SELMER_SQUARE.md already quotes p-parity
  as "quoted, not used". Bounds 3 and 5 for the other two curves are not decided by parity.
- Seven conditional ratios (sample = all seven, PARI 2.17.4, effort 1 over the isogeny class for the doubtful ones):
  (3540,2581) [2,2]; (6121,959) [1,3] root no. −1; (4440,2231) [2,2]; (6671,2209) [3,3]; (6480,4961) [2,2];
  (8249,3720) [2,2] (effort 0 gives only [0,2]); (6903,5104) [2,2] (effort 0: [0,2]); partner (12007,1799) [2,2].
  So "PARI rank ≥ 2" is right for six of seven; (6121,959) is [1,3] — REPORT_F4 says this correctly ("or bounds
  [1,3]"); any summary that says "all seven have rank ≥ 2" overstates. Note (6121,959) is the partner of (3540,2581):
  if its g4K had rank 1, the same DM would close the pair (not checked).

## F2 — checks (21:20 UTC)

- Reproduction: all F2 lines reproduce their committed logs with real sympy (table above). Shim irrelevant.
- (8249,3720) Θ₂: bound before = log2(2+3·2) + log2(4+3·4) − 2 = 3 + 4 − 2 = 5. Matrix rows (E′ coset reps
  113, 2991753233, 338068115329 = 113·2991753233, so the reps form a group mod known) × 3 open E cosets:
  (+,+,+), (−,+,−), (−,+,−). Bilinearity consistent in both directions (row3 = row1·row2; col3 = col1·col2). F₂-rank 1,
  so 5 − 2·1 = 3. Point classes −1, 108415202 pair +1 with every lift (control). The pairing code itself is the
  repository's rank_second.analyse_general (24d6ef6, unchanged). HOLDS.
- D1 of (12007,1799): before 3 + 3 − 2 = 4 (E: 2 + 3·2 = 8; E′: 4 + 1·4 = 8), after removing the open E′ coset
  3 + 2 − 2 = 3. Independent check with my own 2-adic solver (/tmp/x3/mine/q2_indep.py; P¹ charts, Taylor-expansion
  determination of each disc — a different algorithm from rank_indep.solvable): all four open classes have no
  Q₂-point ("all discs closed"); positive controls: point classes −1, 244673, −244673 all found soluble (rank_indep at
  depth 12 left two undecided). HOLDS.
- F2 bottom line (no rank 1 for any of the seven; bounds 3, 5, 3, 3, 3): HOLDS.

## F1 — checks (21:20 UTC)

Files: bundle ctp.gp vs committed ctp.gp differ only in (i) the real-place skip (`if(ctx[1] == 0, if(abs(L) < ...,
return(-1)); ...)`) and (ii) extra diagnostics in the INCONSISTENT print. ctx[1] == 0 only when isre(pl) (lctx), and
ℚ(i) has no real place, so (i) cannot change E₁/E₄ — claim HOLDS. Note: when local samples disagree the code only
prints INCONSISTENT and uses vals[1] (not fatal) — so "the value the program used" = 1 is correct, but that design
choice should be stated. No INCONSISTENT line in run627.out / run627b.log.
- val_Q.log (bundle) is byte-identical to the committed val_Q.log; val_Q_run1.log = final; bundle val_Q.err is empty,
  committed val_Q.err holds the `INCONSISTENT at inf [1, 1, 0]` line. So "rerun clean, same 9 matrices" HOLDS.
- Controls: over ℚ, 9/9 have pairing rank = dimSel² − 2 − rank and direct == shortcut; over ℚ(i), 8/8 (base changes)
  match, but val_Qi has NO direct-vs-shortcut column: the Hilbert-reciprocity shortcut (the only method used for
  E₁/E₄, "correction 0") is cross-checked against the direct method only over ℚ. Six ℚ and two ℚ(i) controls have
  pairing rank 2, so the code does output non-zero values (a constant-0 bug would be caught). Reference ranks are
  PARI's. Evidence, not proof (as the bundle says).
- run627b.log: bundle copy = committed copy + one hand-written line ("stopped by me ... during <3,4> ...") appended to
  the machine log. Do not replace the committed log with it; put the note in WRITEUP only.
- run627.out (committed): ⟨1,3⟩, ⟨1,4⟩ are 0 twice; ⟨1,2⟩, ⟨2,3⟩, ⟨2,4⟩ 0 twice; run627b: ⟨1,2⟩, ⟨2,4⟩, ⟨2,3⟩ 0 twice.
  Logs support the claims as stated. HOLDS.
- NEW (mine): the torsion image was "heuristic in the is-a-square direction" (torsion_E4.gp tests squares at 300 split
  primes although its own `sq()` is an exact nfroots test). I re-ran torsion_E4.gp's data through cypari2 with an exact
  square test (nfroots of x² − a over ℚ(i)), /tmp/x3/mine/torsion_exact.py: image(T1) = image(T2) = S3·S4, image of the
  4-torsion point over T3 = S1·S4; exact and heuristic agree on every hit. E₄(K)[2] = (ℤ/2)² so the torsion image has
  dimension exactly 2: span{S3·S4, S1·S4} is now exact. This matters: with ⟨3,4⟩ uncomputed and the five other entries 0,
  ⟨3,4⟩ = 1 would force the torsion image to be span{S1, S2} and the Sel² pairing WOULD prove rank 0. With the exact
  torsion image, ⟨3,4⟩ = ⟨S3, S3·(S3S4)⟩ = 0 is forced (alternating, torsion in the kernel), so "check entry" is right.
- Deciding-entry algebra: S1 ≡ S3 ≡ S4 mod T, so Sel²/T = ⟨[S1],[S2]⟩ and ⟨1,2⟩ = ⟨2,3⟩ = ⟨2,4⟩ is the one deciding
  value. HOLDS.
- Hidden hypotheses not stated in REPORT_F1: local terms are computed at sampled local points to finite p-adic
  precision N and mapped back (back()) before nfhilbert; there is no check that N exceeds v(L) + 2v(2) + 1, so
  "all symbols exact" means "exact for the approximation". Agreement across 3 samples and 2 runs is the only guard.
- F1 verdict: HOLDS WITH GAPS (unreviewed implementation; shortcut cross-checked only over ℚ; p-adic precision not
  asserted; ⟨3,4⟩ not computed — now forced to 0 by the exact torsion image).

### WRITEUP_vs_24d6ef6.diff — line by line ("-" = 24d6ef6, "+" = Fable's copy)
1. Hunk @52, "+ Validation (computed, exact) ... four/validation/val4.sh, val4.log ...": SAFE if four/validation/val4.sh
   and val4.log are committed with it; fix val4.sh's hard-coded `cd /mnt/project-files/mss-experiments/open-ratio/four`
   (author's machine) to `cd "$(dirname "$0")/.."`. val4.log matches the claim (bounds 0,1,1,2,2,3 = PARI ranks).
2. Hunk @74, "New shortcut": DO NOT TAKE the "+" line ("mine and unreviewed"); 24d6ef6's "-" line is newer (states the
   hypotheses and the separate AI review). Keep 24d6ef6.
3. Over ℚ(i) line: DO NOT TAKE "+"; it drops 24d6ef6's [0,−22,−33] exception. Keep "-".
4. "+ Caveat on the check" (Ш[2] ∩ 2Ш = 0; evidence not proof; ℚ shortcut check independent of it): SAFE.
5. "+ One inconsistent local value, found and fixed" with Cause / Fix / Effect: SAFE together with the ctp.gp fix,
   val_Q_run1.err, val_Q_run2.err and the empty val_Q.err (replacing the committed one). Add: "when samples disagree
   the program prints INCONSISTENT and uses the first sample".
6. "+ a ninth test curve, [0,240,−283], ran out of time.": DO NOT TAKE (duplicates "Did not finish", misnumbered).
7. Hunk @92 Evidence 2 "+ I stopped it during ⟨3,4⟩ ...": SAFE (replaces the stale "still running"), but change
   "(a pairing between two torsion-related classes, which must be 0)" to "(which must be 0: S3·S4 lies in the torsion
   image; an AI referee confirmed that image with an exact square test, see below)".

## Wording (item 4) — what overstates, and corrected text

Overstatements found:
- HANDOFF table, F4 row: bold "**Closes by Demjanenko–Manin ... The same argument also closes (6557,324) and
  (6881,6233).**" The tier cell is right, the headline is not: nothing closes; it would close IF rank C(ℚ) = 1.
- HANDOFF/REPORT_F4 "Rank 1 is PARI's ellrank upper bound, which uses the Cassels–Tate pairing" — omits that it is
  unverified, GRH-dependent through bnfinit unless certified, and (for 6557/6881) holds on one isogenous model only.
- HANDOFF: "F4 gives (6557,324) a non-numerical route" — true of the DM step only; the route is conditional on a
  program output, so it is not stronger than the numerical L-value closure, just different.
- HANDOFF/REPORT_F2: "The shim reproduces the committed rank_indep and rank_fisher outputs byte for byte" — only in the
  author's PARI 2.19.1 environment; with PARI 2.17.4 it fails closed on 7 of 11 CI lines. REPORT_F2's CI path
  `../../scripts/...` is wrong (`../../../`). REPORT_F4 refers to `../f2/sympy_pari_shim.py` (no such dir).
- HANDOFF F1 tier "computed (exact, given ctp.gp's hypotheses)": the local terms come from sampled points at finite
  p-adic precision; "exact" should go.
- REPORT_F1 "validated on 17 control curves" (in the LEDGER text): "consistent with PARI on 17 controls" — evidence.
- Not overstated: REPORT_F4's "(6121,959) ... bounds [1,3]"; REPORT_F2's "rank 1 is not proved for any".

Proposed LEDGER/CHANGES text (replacing the bundle's):

F4 > **Experiment, not counted.** For a G1 ratio, w = t² gives a point of C: y² = g4(w)K(w), and t → 1/t a second
point with x·x′ = c3² (x = c3·w), so Demjanenko–Manin applies to C. `experiments/census-wide/g4k/g4k_exact.py` does
the DM step in exact arithmetic (torsion ℤ/2; an enclosure of ĥ(G); saturation at every prime ≤ 199; an index bound
from the Cremona–Prickett–Siksek height bound on the minimal model; the DM bound and every candidate). It shows that
**if rank C(ℚ) = 1**, every solution is degenerate for (7735,4968), (6557,324) and (6881,6233). Rank 1 is not proved:
it is PARI's `ellrank` upper bound (an unverified program output; its class groups are GRH-dependent unless certified;
for (6557,324) and (6881,6233) only one curve of the isogeny class gives the bound 1), and the repository's exact
descent gives 2, 3 and 5. Tier: conditional on rank C(ℚ) = 1 and on the published Silverman and CPS statements quoted
(tier E). A separate AI review re-derived the DM step with its own exact code and got the same bounds and candidates.
Written by an AI model; not read by a person.

F2 > **Experiment, not counted.** Exact rank bounds from the repository's rank_indep/rank_second (24d6ef6, unchanged)
for the five curves behind the seven conditional ratios: 3 (E1 of (3540,2581)), 5 (E1 of (4440,2231)), 3 (E1 of
(6480,4961)), 3 (E1 of (8249,3720); was 5: Θ₂ matrix of F₂-rank 1), 3 (D1 of (12007,1799); was 4: no ℚ₂-point on any
class of the open coset). None reaches 1; the seven ratios stay conditional on rank 1. Computed (exact). Written by an
AI model; not read by a person.

F1 > **Experiment, not counted.** A second implementation of the Cassels–Tate pairing on Sel² over ℚ(i)
(`open-ratio/ctp/ctp.gp`, unreviewed; local terms from sampled local points; places of good reduction through a
Hilbert-reciprocity shortcut that is the session's own argument, cross-checked against the direct method over ℚ only)
gives 0 on all six entries of Sel²(E₁) and on the deciding value ⟨S1,S2⟩ of Sel²(E₄) for (627,364), as run 4 found.
On 17 controls (9 over ℚ, 8 base changes to ℚ(i)) its pairing rank agrees with PARI's ranks: evidence, not proof. So
Sel² with its pairing cannot prove rank 0 for E₄; (627,364) stays open. Written by an AI model; not read by a person.

CI: drop every `import sympy_pari_shim; sys.modules['sympy'] = sympy_pari_shim` line (and run_shim.py,
sympy_pari_shim.py) before committing; then the eleven lines pass here as listed in the real-sympy table.

## Final (21:25 UTC)

Verdicts:
1. F4 — DM step: HOLDS (independently re-derived and recomputed: same DM rhs, amax 18/19/39, candidates 37/39/80,
   t = ±1, ±1, 0 degenerate; saturation ≤ 199 and CPS index bound reproduced by own code, l = 2 by a different
   criterion; minimality at 2 checked). Closure claim: HOLDS WITH GAPS — conditional on rank C(ℚ) = 1, which is PARI-only
   (exact 2, 3, 5); headline "closes" FAILS as worded. Seven-ratio claim: holds as REPORT_F4 words it (six ≥ 2,
   (6121,959) [1,3]). Possible upgrade for (7735,4968) via 2-parity (not checked).
2. F2 — HOLDS (all outputs reproduced with real sympy; 8249 matrix arithmetic and bilinearity checked; D1 2-adic
   absence confirmed by an independent solver). Shim claim: FAILS here (fails closed, 7/11 lines), but cannot have
   changed an outcome.
3. F1 — HOLDS WITH GAPS (logs support the claims; torsion image now exact by my re-run, which makes ⟨3,4⟩ = 0 forced;
   shortcut cross-checked over ℚ only; p-adic precision not asserted). Diff lines: take 1 (with val4.sh path fix),
   4, 5, 7 (reworded); do not take 2, 3, 6.
4. Wording — corrected text above.

Files: commit 5289548 on run5_x3: experiments/referee_x3/ (dm_indep.py, index_indep.py, neg_indep.py 4/4 controls,
q2_indep.py, torsion_exact.py, rank_pari.py, run_real_sympy.py, ci_bundle*.sh). The bundle copies
experiments/census-wide/{g4k,rank1}/ are left UNTRACKED in the worktree (not committed). Logs: /tmp/x3/logs/.

## Packaging (21:35 UTC; coordinator's follow-up job)

Committed on run5_x3 (one new commit; hash in the final message):
- `experiments/census-wide/g4k/`: dm_g4k.py, g4k_exact.py, g4k_rank_exact.py and their outputs (dm_g4k.jsonl,
  dm_g4k_negative.txt, g4k_exact.txt, g4k_exact_negative.txt, g4k_rank_exact.txt, g4k_rank_exact_negative.txt), NOTES.md
  (my corrected F4 text). Edits: the shim import lines and the `../f2` path removed (real sympy); the fallback REPO path
  to the session's /tmp folder removed (REPO = repository root, or $MSS_REPO); g4k_rank_exact's docstring now says
  "unchanged, with sympy". Not committed: sympy_pari_shim.py, timings.txt (the session's machine), REPORT_F4.md.
- `experiments/census-wide/rank1/`: f2_curves.py (+ hints json), f2_general.py (+ _3.log, _4.log), f2_general_neg.py
  (+ negative.txt), f2_pairing.py (+ _0, _1, _2 logs), f2_d1_local2.py (+ .txt, _negative.txt), NOTES.md (my corrected
  F2 text). Edits: shim import lines removed; docstrings that described the shim corrected. Not committed:
  sympy_pari_shim.py, run_shim.py, shimcheck.log and the two empty shimcheck diffs, *.time, REPORT_F2.md, and
  f2_pairing_3.log / f2_pairing_4.log: they only show the script's own assertion stop (more than one open coset);
  REPORT_F2 does not rely on them (it uses analyse_general there). f2_pairing_0..2 ARE relied on (3540 and 4440: no open
  coset on one side; 6480: every pairing value +1), so they are kept and rerun. f2_general_4.log is kept (it shows Θ₂
  gives no drop on D1, which is why f2_d1_local2 exists).
- `experiments/open-ratio/ctp/`: ctp.gp with the real-place fix; validation/val_Q_run1.log, val_Q_run2.log,
  val_Q_run1.err, val_Q_run2.err (new) and val_Q.err (now empty). run627b.log and run627b_attempt1.log: kept as
  committed — attempt1 is byte-identical to the bundle's; the bundle's run627b.log is the committed log plus one
  hand-written note line, and that note is now in WRITEUP.md instead (a log should stay program output).
- `experiments/open-ratio/four/validation/`: val4.sh (path → `cd "$(dirname "$0")/.."`), val4.log.
- `experiments/open-ratio/WRITEUP.md`: diff lines 1, 4, 5 (+ "When the samples at a place disagree, the program prints
  INCONSISTENT and uses the first sample; it does not stop."), 7 reworded (⟨3,4⟩ must be 0 because S3·S4 lies in the
  torsion image, confirmed exactly by experiments/referee_x3/torsion_exact.py). Not 2, 3, 6.
- `experiments/README.md`: (7735,4968) "open" → excluded only conditionally on rank C(ℚ) = 1 of g4·K (PARI's bound;
  exact 2); (6557,324)/(6881,6233) get the second, also conditional, route; F2 bounds 3, 5, 3, 3, 3; F1's second
  implementation (deciding value of Sel²(E₄) zero; unreviewed); a "Run 5" paragraph listing what was changed and left out.
- `reviews/review_x3.md`: copy of this report as it stands at commit time.

Pre-run (copy of the working tree, real sympy 1.14, PARI 2.17.4, 2 cores): all 13 lines below PASS; times
5.7, 84.8, 23.7, 7.5, 11.2, 33.6, 11.2, 14.8, 10.4, 26.9, 28.7, 19.2 s (≈ 4.6 min). Clean-checkout run: see below.

Proposed lines for the `descents` job (append before its final `cd "$GITHUB_WORKSPACE"` / `git diff --exit-code`):
```
          # run 5: Demjanenko-Manin on y^2 = g4 K (exact; conditional on rank 1), its exact rank bounds, the PARI cross-check,
          # and the exact rank bounds for the seven conditional ratios (about 5 minutes)
          cd "$GITHUB_WORKSPACE"
          f5() { grep -v -e 'PARI stack' -e '^# matches'; }
          G=experiments/census-wide/g4k; K=experiments/census-wide/rank1
          python3 $G/g4k_exact.py 2>&1 | f5 | diff - $G/g4k_exact.txt
          python3 $G/g4k_exact.py --negative 2>&1 | f5 | diff - $G/g4k_exact_negative.txt
          python3 $G/g4k_rank_exact.py 2>&1 | f5 | diff - $G/g4k_rank_exact.txt
          python3 $G/g4k_rank_exact.py --negative 2>&1 | f5 | diff - $G/g4k_rank_exact_negative.txt
          python3 $G/dm_g4k.py 2>&1 | f5 | diff - $G/dm_g4k.jsonl
          python3 $G/dm_g4k.py --negative 2>&1 | f5 | diff - $G/dm_g4k_negative.txt
          for k in 3 4; do python3 $K/f2_general.py $k 2>&1 | f5 | diff - $K/f2_general_$k.log; done
          python3 $K/f2_general_neg.py 2>&1 | f5 | diff - $K/f2_general_negative.txt
          python3 $K/f2_d1_local2.py 2>&1 | f5 | diff - $K/f2_d1_local2.txt
          python3 $K/f2_d1_local2.py --negative 2>&1 | f5 | diff - $K/f2_d1_local2_negative.txt
          for k in 0 1 2; do python3 $K/f2_pairing.py $k 2>&1 | f5 | diff - $K/f2_pairing_$k.log; done
```
(The scripts also assert equality with their own output files and print `# matches ...`, which is filtered; stderr is
merged so a traceback fails the diff. The job's `shell: bash` runs with pipefail. The quick job's py_compile glob
`experiments/census-wide/*/*.py` already covers the new scripts; they compile.)

Proposed OUTPUTS.md rows (command run from the repository root):
| File | Class | Command | CI job | Time |
|---|---|---|---|---|
| `experiments/census-wide/g4k/g4k_exact.txt` | A | `python3 experiments/census-wide/g4k/g4k_exact.py` | descents | 6 s |
| `experiments/census-wide/g4k/g4k_exact_negative.txt` | A | `... g4k_exact.py --negative` | descents | 85 s |
| `experiments/census-wide/g4k/g4k_rank_exact.txt` | A | `... g4k_rank_exact.py` | descents | 24 s |
| `experiments/census-wide/g4k/g4k_rank_exact_negative.txt` | A | `... g4k_rank_exact.py --negative` | descents | 8 s |
| `experiments/census-wide/g4k/dm_g4k.jsonl` | A | `... dm_g4k.py` | descents | 11 s |
| `experiments/census-wide/g4k/dm_g4k_negative.txt` | A | `... dm_g4k.py --negative` | descents | 34 s |
| `experiments/census-wide/rank1/f2_general_3.log`, `f2_general_4.log` | A | `... rank1/f2_general.py 3` (resp. `4`) | descents | 11 s, 15 s |
| `experiments/census-wide/rank1/f2_general_negative.txt` | A | `... rank1/f2_general_neg.py` | descents | 10 s |
| `experiments/census-wide/rank1/f2_d1_local2.txt` | A | `... rank1/f2_d1_local2.py` | descents | 27 s |
| `experiments/census-wide/rank1/f2_d1_local2_negative.txt` | A | `... rank1/f2_d1_local2.py --negative` | descents | 29 s |
| `experiments/census-wide/rank1/f2_pairing_0.log` … `_2.log` | A | `... rank1/f2_pairing.py k`, k = 0, 1, 2 | descents | 19 s together |
| `experiments/census-wide/rank1/f2_curves_hints.json` | D | read by `f2_curves.py`; written by `python3 f2_curves.py --find` (PARI `ellrank` hints, checked on the curve) | none | |
| `experiments/open-ratio/ctp/validation/val_Q.err` | B | `gp -q val_Q.gp` from `ctp/validation` (writes val_Q.log; stderr → val_Q.err); random local points, so the per-place values vary; about 10 min (session's figure, not rerun here: no `gp` binary) | none | ~10 min |
| `experiments/open-ratio/ctp/validation/val_Q_run1.log`, `val_Q_run2.log`, `val_Q_run1.err`, `val_Q_run2.err` | C | logs of the two runs made with ctp.gp before the real-place fix; not regenerable with the current code | none | |
| `experiments/open-ratio/four/validation/val4.log` | B | `sh experiments/open-ratio/four/validation/val4.sh` (needs `gp`; five run4.py runs with a 3000 s timeout each); not rerun here | none | not measured |
| `reviews/review_x3.md` | C | review note written by an AI agent; no checker | none | |
