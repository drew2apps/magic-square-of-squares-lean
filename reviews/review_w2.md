# Hostile review of run4_w2 (TorsionE1)

Started 01:14 UTC. Reviewing head aeaca4c (re-check before final).

## Verdict (01:30 UTC, head aeaca4c): HOLDS (mathematics and Lean); integration items open

`TorsionE1.lean` proves, standard axioms only (re-elaborated from source by me, EXIT 0), that the torsion of
`E1(ℚ(√2))`, `E1 : Y² = X(X² + 8dX − 8d²)`, `d = 2 + √2`, over LemmaDBridge's `Q2 = ℚ⟮√2⟯ ⊆ ℝ`, is exactly
`{O, (0,0), (4, ±8√2), (−4, ±(16+8√2)), (12+8√2, ±(80+56√2)), (−12−8√2, ∓(32+24√2))} = ⟨D⟩ ≅ ℤ/10`, and
`torsion_abscissa` is verbatim (and type-checks as) w1's hypothesis of `E1Descent.lemmaD_of_torsion`. The reduction
argument is complete: v = (7, √2 − 3) = ker(ℤ[ω] → 𝔽₇, ω ↦ 3), residue field ≅ ZMod 7 proved (first isomorphism
theorem), model over ℤ[ω] integral with Δ ↦ 3 ≠ 0 in 𝔽₇, #Ẽ(𝔽₇) = 10 and its 9 affine points by `decide` on the
right object (W7, proved equal to the base change of W0), every Stoll hypothesis discharged by a proof
(7 ∉ v⁶ via ℤ[ω] → ℤ/49, ω ↦ 10). 7-torsion needs no separate exclusion: Stoll's lemma is injective on all torsion
when e < p − 1 (here e = 1, p = 7), which is correct mathematics. Headline statements mention only Mathlib's
`E1K.Point`, `IsOfFinAddOrder` and explicit points, so they do not depend on how Stoll defines reduction.

FIX-FIRST: none in the w2 mathematics or Lean.
FIX (integration, before claiming "Theorems 1–2 unconditional in Lean"): (1) the end-to-end theorem is in no
commit (W2EndToEnd.lean uncommitted; int4 UniformRank.lean untracked) — commit it with `#print axioms` and let CI
build it; (2) on run4_w2 alone CI fails twice: `quick` (MANIFEST: 4 files unlisted) and `rank` (the old axiom
parser crashes on TorsionE1's three wrapped `#print axioms` messages; int4 has already replaced the parser);
(3) TorsionE1Forge (7 errors) and torsion_e1_check.py are not in CI on either branch; (4) doc text (w2.md:216-224)
must say "coprime pairs" and keep "on the models Eplus m (±n); E⁻ = Eplus m (−n) by sympy".
NIT: see the end.

## Findings (written as I go)

### 01:17 UTC — statements and the Stoll interface (read, not compiled)

- Head still aeaca4c. Diff 80d2b9f..HEAD: TorsionE1.lean (713), TorsionE1Forge.lean (85), torsion_e1_check.py/.txt.
  scripts/W2EndToEnd.lean is NOT in the branch (Mac/container slice dir only, 01:12, uncommitted; needs w1's files).
- `E1K` is LemmaDBridge's (`⟨0, 8*(2+s2), 0, -8*(2+s2)^2, 0⟩` over `Q2 = ℚ⟮√2⟯ ⊆ ℝ`, `s2 = ⟨√2, _⟩`); TorsionE1
  does not redefine it, and pins the coefficients with `rfl` (TorsionE1.lean:273-277). `s2*s2 = 2` is proved
  (TorsionE1.lean:82, via `Real.mul_self_sqrt`). OK.
- `torsion_abscissa` (TorsionE1.lean:530-532) is token-for-token the hypothesis `h` of w1's
  `E1Descent.lemmaD_of_torsion` (w1 E1RankZero.lean:43-46, head 0da8299): same binders `(x y : Q2) (hxy :
  E1K.Nonsingular x y)`, same `IsOfFinAddOrder (WeierstrassCurve.Affine.Point.some x y hxy)`, same five
  disjuncts in the same order. Both elaborate `Q2/E1K/s2` from `MagicSquareOfSquares.Family`. OK.
- Stoll's `eq_zero_of_isOfFinAddOrder_of_red_eq_zero` (ReductionAtPrime.lean, pkg 2164fb5) needs: Dedekind `R`,
  `IsFractionRing R K`, `hE : W₀.map (algebraMap R K) = E`, `[(redCurve v W₀).IsElliptic]`, `[E.IsElliptic]`,
  `[CharZero K]`, `[DecidableEq (R ⧸ v)]`, `[DecidableEq K]`, `p.Prime`, `(p:R) ∈ v`, `(p:R) ∉ v^(p-1)`, torsion `P`,
  `red P = 0`. No prime-to-p hypothesis: correct mathematics (e = v(p) < p−1 kills all torsion of the formal group,
  Silverman IV.3.2 + IV.6.1 / VII.3.1). Each is discharged in TorsionE1 by a proof: Dedekind via
  `QuadraticAlgebra.Int.isIntegralClosure_iff` + `fund8` (l.62-78); IsFractionRing by `Q2_elem` (l.113-149);
  good reduction by `decide` on Δ(W7) after `map_Δ` (l.247-256); `seven_mem` (l.170), `seven_not_mem_pow` (l.202,
  via ψ₄₉ : ω ↦ 10, 10² = 100 ≡ 2 mod 49, exhaustive 49² `zmod49_key`); `E1K.IsElliptic` from LemmaDBridge:99.
  No `sorry`/`axiom`/`native_decide` in the EllipticCurves package except comparator/Challenge.lean (sorry'd
  restatements, not imported) — grep.

### 01:21 UTC — my one compile (independent re-elaboration of the head file)
- `scripts/ReviewW2.lean` in the w2 slice = byte copy of TorsionE1.lean at aeaca4c (sha256 0548f039…, same as
  worktree, Mac and container) + 39 lines of referee checks; it imports only what the module imports (not the
  TorsionE1 olean), so the whole file was re-elaborated from source. `lean4.sh w2 check`: **EXIT 0, 39 s, no
  error/warning/sorry in last.log**; all 13 `#print axioms` = [propext, Classical.choice, Quot.sound]; Stoll's
  `eq_zero_of_isOfFinAddOrder_of_red_eq_zero` likewise standard.
- Referee checks accepted: (i) w1's hypothesis typed verbatim **inside `namespace MagicSquareOfSquares.E1Descent`
  with `open MagicSquareOfSquares.Family`** (w1's context) is closed by `TorsionE1.torsion_abscissa`;
  (ii) `E1K = ⟨0, 8(2+s2), 0, −8(2+s2)², 0⟩` and `(s2 : ℝ) = √2` by `rfl`; (iii) `⟨−3, 1⟩ = ω − 3 ∈ v7`,
  `⟨3, 1⟩ = ω + 3 ∉ v7`, `Nat.card (O ⧸ v7) = 7`: v7 = (7, √2 − 3), residue field of order 7, proved.
- Elaborated statements (pp): `torsion_E1K : ∀ P : E1K.Point, IsOfFinAddOrder P → P = 0 ∨ P = pT ∨ … ∨ P = pDm`;
  `card_torsion_E1K : Nat.card {P // IsOfFinAddOrder P} = 10`; `torsion_iff_mem_zmultiples : IsOfFinAddOrder P ↔
  P ∈ AddSubgroup.zmultiples pD`; `LemmaD := ∀ x w : ℚ, w² = 2(x²+1)(x⁴−2x³+2x²+2x+1) → x = 1 ∨ x = −1`.

### 01:21 UTC — mathematics re-derived by hand / PARI (independent of the build agent's script)
- Ẽ : y² = x³ + 5x² + 3x over 𝔽₇ (8·5 = 40 ≡ 5, −8·25 = −200 ≡ 3): x=0→1 pt, x=1→2 (y²=2), x=3→2 (y²=4), x=4→2
  (y²=2), x=6→2 (y²=1), x=2,5 non-squares: 9 affine + O = 10. Δ(Ẽ) = 16·3²·(25−12) ≡ 3 ≠ 0.
  Over K: Δ = 2¹⁵·3·d⁶, PARI idealfactor = 𝔭₂³⁶·(3): good at both primes above 7, #Ẽ = 10 at both.
- Reductions under √2 ↦ 3: T→(0,0), ±A→(4,3/4), ±B→(3,5/2), ±C→(1,3/4), ±D→(6,1/6): matches `W7_points` order.
- PARI `ellmul(E, D, k)`, k = 1..9: x = −12−8√2, 12+8√2, −4, 4, 0, 4, −4, 12+8√2, −12−8√2; 10D = O; elltors =
  [10,[10]] with generator D. So ⟨D⟩ = the ten points of `torsion_E1K`, and the abscissa list of
  `torsion_abscissa` is exactly {x(kD) : 1 ≤ k ≤ 9}. Build agent's torsion_e1_check.py re-run here (python3 -I,
  cypari2): output identical to the committed .txt.
- Logic of `torsion_E1K` (l.438-476): for torsion P, case on resPt(red P) ∈ Ẽ(𝔽₇) (10 cases by `decide`); in each
  case an explicit torsion Pᵢ with the same reduction (`resPt_red`, coordinates checked by `decide`), then
  P − Pᵢ is torsion with red = 0 (red is a hom: Stoll's `redHom`), so P = Pᵢ by Stoll. Injectivity holds on the
  whole torsion (7-torsion included) because e(v/7) = 1 < 6 — this is the hypothesis `7 ∉ v⁶`, proved
  (no second prime needed, and no "order prime to 7" argument is used or needed). Sound.

### 01:27 UTC — FIX-FIRST candidate: the CI axiom audit crashes on TorsionE1's wrapped `#print axioms` lines
- Three of TorsionE1's 13 axiom messages are longer than Lean's 120-column format width and are printed on three
  lines: `…TorsionE1.lemmaD_of_finrank_eq_zero' depends on axioms: [propext,⏎ Classical.choice,⏎ Quot.sound]`
  (also `natCard_E1K_of_finrank_eq_zero`, `torsion_iff_mem_zmultiples`; TorsionE1.lean:711-713). Seen in my
  compile output and in the cached lake log `/work/mss/w2/lean_rank/.lake/build/lib/lean/LeanRank/TorsionE1.trace`
  (`"…depends on axioms: [propext,\n Classical.choice,\n Quot.sound]"`), i.e. `lake build` replays them wrapped.
- The `rank` job (.github/workflows/ci.yml:95-99, runs on every push) does, per line,
  `re.search(r'\[(.*)\]', l).group(1)`: on the first line of a wrapped message there is no `]`, `re.search` is
  `None` → AttributeError → job fails. I ran that exact snippet on my compile log: "AttributeError 'NoneType' object
  has no attribute 'group'". No earlier module has a message over 120 columns (longest 113, t7/v8 logs), so this is
  new with w2. The build report's "89 axiom lines, none non-standard" was not produced by the CI parser.
- Also a soundness hazard of the parser once "fixed" naively: a wrapped line that lists a non-standard axiom on a
  continuation line would be invisible to a per-line filter.
- Fix (either): parse build.txt as a whole, `re.findall(r"depends on axioms: \[([^\]]*)\]", text)` with
  whitespace stripped from each item; or put `set_option format.width 1000 in` before each `#print axioms` (or
  shorten nothing and rely on the parser fix). Then re-count (89 on w2 alone).
- 01:25 UTC update: the integration worktree int4 (branch run4, merge 2c05c01 of run4_w2) already replaced the
  parser by a whole-file regex that handles wrapped lists (int4 .github/workflows/ci.yml:95-110, comment "a long
  name makes Lean wrap the axiom list over several lines") and also checks every `#print axioms` source line is
  answered. So: **downgraded to FIX** — true for run4_w2 on its own (merging it onto the old CI breaks `rank`), already
  handled in the integration; make sure the old parser does not come back with any cherry-pick.

### 01:28 UTC — integration hygiene (run4_w2 on its own)
- `python3 scripts/verify_manifest.py` on run4_w2: MANIFEST FAILED, the 4 new files "present but not listed" → the
  `quick` job fails on this branch as committed. The build report (w2.md:155-162) leaves MANIFEST/CI edits to the
  conductor by design; int4's MANIFEST (544 lines) does not list them yet either (int4 is mid-integration: it also
  fails on w5's files). FIX at integration.
- Neither run4_w2 nor int4's ci.yml runs `scripts/TorsionE1Forge.lean` (7 expected errors) or
  `torsion_e1_check.py | diff`; the forge and the Python check are therefore not guarded by CI. FIX at integration.
- `scripts/W2EndToEnd.lean` is not committed anywhere (Mac/container slice dir only, mtime 01:12). The end-to-end
  claim ("`lemmaD_end_to_end : LemmaD`, standard axioms; `rank_family_all lemmaD_end_to_end` etc.") rests on the
  build agent's container run with w1's four files copied in and later deleted; I could not re-run it (my one
  compile was in w2, without w1's files). What I did check: (a) the type of `torsion_abscissa` closes w1's
  hypothesis typed verbatim in w1's namespace/`open` context (compile above); (b) w1's four files add only the
  global instances `Div/Mod/EuclideanDomain (ℤ√2)` (Zsqrt2.lean:33,36,96) and `IsAddFG E1K.Point`
  (E1RankZero.lean:19) — none touches `DecidableEq Q2` or the `AddCommGroup` of `E1K.Point`, so the composition
  cannot fail by an instance mismatch; (c) the conductor's untracked int4 `lean_rank/LeanRank/UniformRank.lean`
  contains the same composition (`lemmaD := E1Descent.lemmaD_of_torsion TorsionE1.torsion_abscissa`) and restates
  `theorem1`, `theorem1_P`, `theorem2` without `hD` — not compiled by me.
- Family statements with `hD` discharged (read from the sources, w2 worktree): `rank_family_all`
  (LemmaDEquiv.lean:183): `IsCoprime m n → m²+n² = 2k² → m² ≠ n² → Theorem2At m n k ∧ Theorem2At m (−n) k ∧
  2 ≤ finrank E⁺ ∧ 2 ≤ finrank E⁻`; `rank_one_family` (Theorem1Family.lean:191), `R_infinite_order`,
  `P_infinite_order` (Theorem1Order.lean:37,103): `IsCoprime m n → m ≠ 0 → n ≠ 0 → m² ≠ n² → …`
  (= mn(m²−n²) ≠ 0); `theorem2_family`, `rank_family` (Theorem2Family.lean:133,149) additionally `1 ≤ n < m`
  (subsumed by `rank_family_all`); `twiceSquareH_of_lemmaD` (LemmaDEquiv.lean:35). No other hypothesis. These match
  the requested hypotheses. Two caveats that the doc text must keep (w2.md:216-224 drops them): Lean's Theorem 1
  is for **coprime** pairs (UNIFORM_RANK Theorem 1 is stated for all non-zero m ≠ ±n and reduces to coprime by
  the scaling isomorphism, which is not in Lean), and the statements are on the models `Eplus m (±n)`; that
  `Eplus m (−n)` is the model of `E⁻` is sympy only (TwoDescentMap.lean:53, UNIFORM_RANK §7).

## NIT
- N1 TorsionE1.lean:59-189: the instances `Algebra O Q2` (l.93), `FaithfulSMul O Q2`, `IsFractionRing O Q2`,
  `Algebra O (ZMod 7)` (l.159), `DecidableEq (O ⧸ v7)`, `Fact (Nat.Prime 7)`, `IsIntegralClosure O ℤ K'` are global
  and reach every importer (UniformRank, later files). `Algebra (QuadraticAlgebra ℤ 2 0) (ZMod 7)` with ω ↦ 3 is
  diamond-prone (the conjugate ω ↦ 4 is equally natural). Make them `local instance` / `scoped`.
- N2 Name clash with w1: `E1PointsK_of_torsion`, `lemmaD_of_torsion` exist in `MagicSquareOfSquares.E1Descent`
  (hypothesis: abscissas of torsion points) and in `MagicSquareOfSquares.Family.TorsionE1` (hypothesis: every point
  torsion) with different statements; with both namespaces open they are ambiguous. Rename w2's
  (e.g. `lemmaD_of_isTorsion`).
- N3 `card_torsion_E1K` (l.616) is on a subtype; also stating `AddCommGroup.torsion E1K.Point =
  AddSubgroup.zmultiples pD` would put the result in Mathlib's vocabulary (one line from
  `torsion_iff_mem_zmultiples`).
- N4 l.342 "the generator of v8's report": an internal report, not a repo file; say "PARI's elltors generator".
- N5 torsion_e1_check.py negative tests are shallow (e.g. "wrong count 11" checks 10 ≠ 11, "wrong lift" checks
  3² − 2 ≢ 0 mod 49): they test constants, not the checking code. The Lean forge is the real negative test.
- N6 w2.md:101-102 proposed §4 sentence "Stoll's injectivity on torsion needs no coprimality": add "since the
  ramification index of v over 7 is 1 < 7 − 1"; in general the condition e < p − 1 is needed.
- N7 w2.md:212-213 "89 axiom lines, none non-standard" was counted by grep, not by CI's parser (which crashes on
  run4_w2's CI; see above).

## Could not check
- Forge file not re-run (my one compile went to the main file); error lines 26/32/35/40/46/50/54 match F1–F7 by
  reading (F5 at l.46 = the y-certificate, as it should be, since 2A and C share x).
- W2EndToEnd.lean / int4 UniformRank.lean not compiled with w1's files (see above for the type-level check done).
- Stoll's proofs (formal group, Filtration, Reduction) taken as compiled with standard axioms; I read only the
  statements of ReductionAtPrime.lean and the adicRed lemma. w1's rank 0 not re-reviewed (out of scope).

### 01:37 UTC — final re-check
- `git log` on run4_w2: head still aeaca4c; working tree clean. Build report w2.md "Checkpoint (01:30 UTC)": no
  further changes to TorsionE1.lean or the forge file. Verdict above stands for aeaca4c.
