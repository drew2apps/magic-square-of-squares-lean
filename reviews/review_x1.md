# Review y1 of run5_x1 (01e35dd..run5_x1, head 7829594)

Hostile referee, written by an AI model; not read by a person. Started 21:24 UTC.

## Progress log

### 21:27 UTC: first reading (source only)
- Diff 01e35dd..run5_x1 = 5 commits (ba44d6a..7829594), 6 files, only new files + lakefile root line. Worktree clean.
- `weierstrassOf q z0` (EminusModel.lean l.38-44): a1=a3=0, a2=q''(z0)/2, a4=q'(z0)q'''(z0)/6, a6=q'(z0)^2 q''''(z0)/24,
  with Mathlib `Polynomial.derivative`/`eval` on an arbitrary `q : ℚ[X]`. This is the standard quartic-with-rational-root
  recipe (u=z-z0, q=Au+Bu^2+Cu^3+Du^4, w=1/u, X=Aw): no answer pasted in. `qPoly s = (X + C s)·gPoly`; `EplusModel = weierstrassOf (qPoly 2) (-2)`,
  `Eminus = weierstrassOf (qPoly (-2)) 2`, matching UNIFORM_RANK.md §1 ("q(z) = (z ± 2)·g(z) and z0 = ∓2").
  `gPoly` = (2Sz − 8mn)(m²n²z² + (m² − n²)²), as §1. Derivatives obtained by simp with product rule + ring (qPoly_d1..d4).
- lakefile: LeanRank lib uses `globs = ["LeanRank.+"]`, so the 3 new LeanRank modules are built by `lake build`; CI's axiom check
  globs `LeanRank/*.lean`, so their `#print axioms` lines are audited automatically. Defs.lean imports only Mathlib, no existing module
  imports it; adding it as a root cannot change any existing module. Forge files are NOT in CI's forge spec list (CI unchanged) —
  needs the conductor to add `'TwoDescentCertForge:…' 'UniformRankAllForge:…'`.
- Stale wording after the lakefile change: CensusClosure.lean header ("lean_rank compiles only four main-library modules"),
  UNIFORM_RANK.md §7 ("compiles the four modules it needs"), ci.yml comment ("compiles four modules of the main tree").

### 21:34 UTC: independent checks run
- Container files identical to the worktree (md5 of the 5 new files, lakefile, Defs.lean).
- Sympy (own script, /tmp/y1/recipe.py, cloud): the §1 recipe applied to (z+2)g at −2 and (z−2)g at 2 gives, coefficient by
  coefficient for symbolic m,n, Eplus(m,n) and Eplus(m,−n) (and NOT Eplus(m,n) for E⁻); A± = ∓4(m±n)²S². At (7,1):
  E⁺ (a2,a4,a6) = (300176, 20572160000, 2007040000000000), E⁻ = (278224, 9596160000, 635040000000000); at (12,5) also matched.
  `sigma_exact.py`'s g (built from the sextic coefficients) equals `gPoly` symbolically. The map from y² = F0·F8·GL through
  z = t + 1/t, Y = (t ± 1)y/t² lands on the model (symbolic identity).
- Lean probe 1 (scripts/y1_probe.lean in /work/mss/x1/lean_rank, `lake env lean`): #check of every headline statement (texts
  recorded below), #print of weierstrassOf/EplusModel/Eminus/qPoly/gPoly/mapA/IsMSS/IsMagic/Theorem2At; axioms of the 11
  headline theorems = [propext, Classical.choice, Quot.sound]. Non-vacuity: theorem1_all at (7,1), (2,3); theorem1_P_all
  and _minus at (4,1); theorem2_all at (7,1,5), (1,−7,−5), (21,3,15); theorem2_all_Eminus at (17/2,31/2,25/2);
  sextic_to_Eplus at (t,y) = (1,600), sextic_to_Eminus at (−1,800) for (7,1): all elaborate. The recipe is general:
  weierstrassOf (X⁴ + X) 0 = y² = x³ + 1 (a2=a4=0, a6=1), and on (X−1)(X³+2) at 1 it gives (3, 9, 9) (my first hand value
  a4 = 3 was wrong and Lean refused it, probe 2 has the corrected value). (Eminus 7 1).a_i and (EplusModel 7 1).a_i equal the
  sympy numbers. Eminus 7 1 ≠ EplusModel 7 1. The 180° rotation and a sign flip of d derive (b+d):(h+f) ≠ 3:4 and
  (b−d):(h+f) ≠ 4:3 from no_mss_4_3 (so the file covers all sign variants, as it should).
- Copies of the three library files compiled from source as scripts/y1_copy_*.lean: no error; 30 `#print axioms` answers,
  all ⊆ {propext, Classical.choice, Quot.sound} (bd_bf_of_nodup: [propext]; control_distinct_not_vacuous: [propext, Quot.sound]).
- Forge files (byte copies as scripts/y1_*Forge.lean): TwoDescentCertForge errors at 35 43 47 51 55 65 69 73 77 81 85 89 97
  102; UniformRankAllForge at 26 30 34 38 41 45 56 60 64 — both = the headers. Columns checked against the source: every error
  is at the `by` of the intended hypothesis (genH value, bracket bound, hk, hne, hz, decide), message "unsolved goals ⊢ False",
  "Type mismatch Theorem2At 7 (-1) 5 / Theorem2At 7 1 5", or "decide proved ... is false". No elaboration/typo errors.
- No `sorry`, `native_decide`, `axiom`, `implemented_by`, `unsafe`, `set_option` (any) in the 5 new files (grep).
- 30 new `#print axioms` lines (10 + 11 + 9); base 147 → 177, matching x1's CI replica (ci_x1.log on the Mac).

### 21:40 UTC: the headline statements in plain mathematics (as Lean has them, `#check` output)
- `theorem1_all`: for all m, n ∈ ℚ with m ≠ 0, n ≠ 0, m² ≠ n²: rank_ℤ E(ℚ) ≥ 1 for E = Eplus m n and for E = Eplus m (−n),
  where Eplus m n is y² = (X + 2S³)(X² + 16m²n²(m+n)²X + 16m²n²(m+n)⁴S²), S = m² + n² (`Module.finrank ℤ` of Mathlib's point
  group; E(ℚ) is f.g. by Stoll, so finrank is the true rank and "≥ 1" is not vacuous).
- `theorem1_P_all`: same hypotheses ⇒ P = (−(m+n)²S², S³(m−n)(m+n)²) is a nonsingular point of Eplus m n and N·P ≠ O for all
  N ≥ 1. `theorem1_P_all_minus`: the same for P(m, −n) on Eplus m (−n) (the image of z = −2 on E⁻, `Eminus_xP`).
- `theorem2_all`: for all m, n, k ∈ ℚ with m² + n² = 2k² and m² ≠ n²: `Theorem2At m n k` and `Theorem2At m (−n) k` (R, Q = (0, yQ),
  T are nonsingular points and c₁R + c₂Q + c₃T ∈ 2E(ℚ) with cᵢ ∈ ℤ forces all cᵢ even), and both ranks ≥ 2. k is an ordinary
  universally quantified variable constrained by hk; instances exist ((7,1,5), (21,3,15), (17/2,31/2,25/2) elaborate).
- `theorem1_all_Eminus`, `theorem2_all_Eminus`: the rank statements on EplusModel m n and Eminus m n (rewritten by the two
  equalities below). `theorem1_int`: the integer special case, no coprimality.
- `EplusModel_eq`: for all m, n ∈ ℚ, weierstrassOf((z+2)g, −2) = Eplus m n as Weierstrass curves (all five aᵢ).
  `Eminus_eq`: weierstrassOf((z−2)g, 2) = Eplus m (−n). Pure polynomial identities, no hypotheses.
- `qPoly_map`: for all s, m, n, z, Y ∈ ℚ with z + s ≠ 0 and Y² = (z+s)g(z): (A/(z+s), A·Y/(z+s)²), A = q′(−s), satisfies the
  equation of weierstrassOf(q, −s). (Equation only; nonsingularity is not stated, automatic when Δ ≠ 0.)
- `sextic_to_Eplus` (`_Eminus`): for t ≠ 0, t ≠ −1 (resp. t ≠ 1) and y² = F0·F8·GL(t) = (m²+n²t²)(n²+m²t²)(2S(t²+1) − 8mnt):
  the point (A/(z − z0), A·Y/(z − z0)²), z = t + 1/t, Y = (t ± 1)y/t², satisfies the equation of Eplus m n (resp. Eminus m n).
- `no_mss_census`: for (m, n) ∈ {(4,3), (7,1), (12,5), (24,7), (17,7), (31,17)} and all integers a…i, T: if a²…i² (rows
  a² b² c² / d² e² f² / g² h² i²) form a magic square with sum T (`IsMSS`, Defs.lean) and the nine squares are pairwise distinct
  (`[a²,…,i²].Nodup`), then (b + d)·n ≠ (h + f)·m. Signs of the roots are free, so all sign variants are covered (probe).

## Verdicts

1. **Statements and vacuity: HOLDS.** Restated above. No unsatisfiable hypothesis: every headline theorem was instantiated at
   concrete pairs (including n > m, negative entries, negative k, non-integral and non-coprime pairs). `weierstrassOf` is not
   degenerate (it gives y² = x³ + 1 on z⁴ + z at 0, (3, 9, 9) on (z−1)(z³+2) at 1, and E⁻ ≠ E⁺ at (7,1)). k is not quantified so
   as to make Theorem 2 vacuous. `Nodup` is on the list of the nine squares, the same expression as `openStatement_iff_exists_isMSS`.
   The only hypotheses of `no_mss_census` that cannot be instantiated are IsMSS ∧ Nodup, which is the open problem itself.
2. **Recipe: HOLDS.** `weierstrassOf q z0` is the §1 recipe for an arbitrary q : ℚ[X] (Mathlib derivative/eval; the answer is not
   pasted in); Eminus = weierstrassOf((z − 2)g, 2), EplusModel = weierstrassOf((z + 2)g, −2), g = (2Sz − 8mn)(m²n²z² + (m²−n²)²),
   exactly §1; `g` equals `sigma_exact.model`'s g and `G1.sigmaG`. Hand/sympy derivation for (7,1) and (12,5) and symbolic m,n
   agrees with Lean. Gap (wording only): the map lemma is proved for qPoly only, not for a general quartic (not claimed); the
   control/forgery "wrong root z0 = 2 for q⁺" is not a root at all (q⁺(2) = 4·g(2) = 16(m − n)²S² ≠ 0); say "at z0 = 2, not a root of q⁺".
3. **CensusClosureMSS: HOLDS WITH GAPS.** Uses `MagicSquareOfSquares.IsMSS`/`IsMagic` compiled from the same Defs.lean (md5 equal),
   fields mapped one to one (`lines`), edge ratio (b + d) : (h + f) as in CensusClosure.lean and G1_FULL_SYSTEM.md (n(b+d) = m(h+f));
   m, n positive coprime literals; distinctness = Statement.lean's. Gap: the link to `OpenStatement` (formal-conjectures form) is
   `openStatement_iff_exists_isMSS` in the main library on v4.33.1, not compiled by lean_rank, so "census closures for the open
   statement" is not one Lean theorem; the doc should say "for IsMSS with distinct entries, the form that the main library proves
   equivalent to the open statement". The reversed ratios n : m and the variants (d+h):(f+b) are not stated (derivable; probe).
4. **Forge files: HOLDS WITH GAPS.** Both files fail exactly at the pinned lines (my run), each at the intended hypothesis, with no
   elaboration error; positive controls compile. Weak tests: in TwoDescentCertForge 7 of 14 claims are TRUE statements with a bad
   certificate (1 Theorem2At 7 1 5, 3, 4, 9 ¬□H(7,1), 6 Theorem2At 7 1 5 again, 10 ¬□5798416, 11 ¬□14); they test that the
   checker refuses wrong witnesses, which is what the header says ("certificates with one thing changed"), not that a false
   theorem is refused. Forgery 6 in particular only shows that `Theorem2At 7 (−1) 5` and `Theorem2At 7 1 5` are different types;
   `local irreducible ... in` is scoped to that one declaration and cannot make a true defeq fail here (−1 ≠ 1 as rationals).
   UniformRankAllForge's 9 claims are all false (3, 4, 5, 6, 7, 8, 9 provably; 1 by PARI's rank 1 for E⁺(12,5); 2 because
   Eplus 3 3 is a nodal curve whose nonsingular points are not finitely generated, so finrank = 0 — my argument, not checked).
5. **lakefile / CI: HOLDS WITH GAPS.** Defs.lean imports only Mathlib and no existing module imports it; adding the root changes no
   existing module and no `#print axioms` line (CI only globs LeanRank/*.lean; Defs has none). `globs = ["LeanRank.+"]` builds the
   three new modules and CI's axiom set comparison covers their 30 lines (177 total). Gaps: (a) the two forge files are not run by CI
   until the conductor adds `'TwoDescentCertForge:35 43 47 51 55 65 69 73 77 81 85 89 97 102' 'UniformRankAllForge:26 30 34 38 41
   45 56 60 64'` (values confirmed by my run); (b) Defs.lean now has to build on two toolchains (v4.33.1 and v4.36.0) — a future
   edit to Defs can break lean_rank; (c) stale "four modules": ci.yml l.78 comment, UNIFORM_RANK.md §7 ("compiles the four modules
   it needs"), and CensusClosure.lean l.38 "(not imported: lean_rank compiles only four main-library modules)" — that file is
   generated from scripts/census_closure_head.lean.in and CI `cmp`s it, so fix the template and regenerate, or leave it.
6. **Forbidden tokens: HOLDS.** None of sorry, native_decide, axiom, implemented_by, unsafe, set_option in the 5 new files; the
   forge files use `decide`/`decide +kernel` only; all 30 axiom answers ⊆ {propext, Classical.choice, Quot.sound}.

## Corrections to the proposed doc text (x1.md), after comparing
- x1's claims (builds, 177 lines, pinned lines, no coprimality needed beyond Lemma B and parity) match what I found.
- "the recipe at the wrong root" (forge list, UNIFORM_RANK §7 paragraph, `control_wrong_root` docstring) → "the recipe at z0 = 2,
  which is not a root of q⁺".
- "Not in Lean: that the list of 127 is complete (`Eplus m (−n) = E⁻` and Lemma D are …)" → "… (`Eplus m (−n)` is, coefficient for
  coefficient, the model of `E⁻` given by the recipe of §1 (`Eminus_eq`), and Lemma D, are …)".
- CHANGES "Theorems 1 and 2 hold in Lean for all rational m, n" → add the hypotheses: "(Theorem 1 when mn(m² − n²) ≠ 0, Theorem 2
  when m² + n² = 2k², m² ≠ n²)".
- LEDGER census row: add "the form `[a², …, i²].Nodup` that `openStatement_iff_exists_isMSS` (main library, v4.33.1, not compiled
  by lean_rank) proves equivalent to the formal-conjectures statement".
- TwoDescentCertForge description: optionally add "in seven of the 14 the statement is true and only the certificate is wrong".
- UNIFORM_RANK.md §7 "compiles the four modules it needs" → "five" (x1's text does not touch this sentence).

## Scratch files left (deleted nothing)
In the container, /work/mss/x1/lean_rank/scripts/: y1_probe.lean, y1_probe2.lean, y1_copy_EminusModel.lean,
y1_copy_UniformRankAll.lean, y1_copy_CensusClosureMSS.lean, y1_TwoDescentCertForge.lean, y1_UniformRankAllForge.lean (the
last two must fail by design). On the Mac: /Users/drew/rh-replay/mss/slices/x1/y1/ (probes and logs). Not in the git worktree;
do not commit them. Cloud: /tmp/y1/recipe.py, /tmp/y1/g_check.py (sympy).

## Final (21:37 UTC)
Bottom line: no mathematical or formal error found. Items 1, 2, 6 HOLD; 3, 4, 5 HOLD WITH GAPS (wording, CI spec lines still
to add, weak-by-design certificate forgeries, stale "four modules", OpenStatement link only in the main library).
