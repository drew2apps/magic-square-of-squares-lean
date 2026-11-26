# Review of t7 (`lean_rank/`, branch `sess_t7`, commits 5c9e7b3, d73e585) — hostile referee

Written by an AI model (reviewer); not read by a person. Read-only review; scratch in /tmp/rev_t7/. No Lean build
started. Mac container: read logs only (`/work/mss/t7/checkA*.log`, `checkA_build.txt`, `forge.txt`), sha256 of the
five sub-project files equal to the author's (Rank 571e5f1e…, Index e8a05b91…, RankChecks c60fa673…, RankForge
3a95afb2…, lakefile 4730e143…). Started 19:45 UTC.

## (1) Statements

- `rank_twiceSquare` (Rank.lean:153–155): for each `p ∈ twiceSquareRatios`, `2 ≤ finrank ℤ (Eplus p.1 p.2.1).toAffine.Point`
  and the same for `Eplus p.1 (−p.2.1)`. `Point` is Mathlib's full point group (zero + nonsingular affine points) of a
  `WeierstrassCurve ℚ`, with Mathlib's group law. Curve identity `Eplus m (±n)` = `E±` was settled in review_t3
  (sympy, `Eplus_cubic`); not re-done.
- Meaning of `finrank`: Mathlib at fc4f59e, `Mathlib/LinearAlgebra/Dimension/Finrank.lean:61`:
  `finrank R M := Cardinal.toNat (Module.rank R M)`. So `finrank = 0` exactly when the rank is infinite (or 0);
  a **lower** bound `2 ≤ finrank` therefore forces `Module.rank ℤ E(ℚ)` finite and `≥ 2`, i.e. two ℤ-independent
  rational points. It cannot be vacuous, with or without fg (the danger is only for upper bounds). fg *is* used in
  the proof (`fg_Eplus`, Rank.lean:131–134, feeds `[IsAddFG G]` of `index_range_nsmul_of_fg`), and with fg,
  `finrank ℤ` is the free rank of `E(ℚ) ≅ ℤ^r × torsion` (Stoll's lemma itself goes through
  `AddCommGroup.equiv_free_prod_directSum_zmod`). Good.
- `index_Eplus` (Index.lean:53–56) carries an **extra hypothesis** `hT : Nonsingular (eP m n) 0` that is never
  discharged; see FIX-FIRST 1.
- `index_twiceSquare` (Index.lean:74–76): `8 ≤ index` with `index : ℕ` (0 if infinite), so non-vacuous; it is weaker
  than what `index_Eplus` + `rank_twiceSquare` give (`= 2^(rank+1) ≥ 8`), consistent.

## (2) Group theory and the calls into Stoll

- `two_le_finrank_of_indep` (Rank.lean:78–114): hypotheses `[IsAddFG G]`, `G[2] ⊆ {0, T₀}` (not torsion-freeness),
  and independence mod 2G of `R, Q, T` with integer coefficients. Proof: 8 sums injective into `G/2G`
  (`eightSums_injective`, an honest Fin-2 parity argument), `#G[2] ≤ 2` from the subset, index formula, arithmetic.
  `T` in `hind` and `T₀` in `htors` are independent variables; the lemma is correct for any of them. No hidden
  assumption.
- `index_two_eq` (Index.lean:37–40): needs `T₀ ≠ 0`, `T₀ + T₀ = 0`, `htors`; gives `#G[2] = 2` exactly
  (`Set.ncard_pair`) and then Stoll's formula. Correct.
- Stoll's `AddSubgroup.index_range_nsmul_of_fg` (EllipticCurves/Mathlib/SelmerGroup.lean:317): `[AddCommGroup G]
  [IsAddFG G] (hn : n ≠ 0) : [G : nG] = n ^ finrank ℤ G * #G[n]` — invoked with `n = 2`, `two_ne_zero`.
- Stoll's `fg_point_of_numberField` (MordellWeil.lean:116, variables line 107–108: `[NumberField F] [DecidableEq F]
  {W : Affine F} [W.toAffine.IsElliptic]`): invoked with `F = ℚ`, the `IsElliptic` instance supplied by
  `have := Eplus_isElliptic hm hn hp hd` (Δ ≠ 0 via `Eplus_Δ`, `qe_P`, `genH_pos`, `sP ≠ 0`). Instances unify
  (build is clean).
- `two_torsion_point` (Rank.lean:121–127) lifts t3's `two_torsion_Eplus` to `Point`; correct case split.

## (3) Forbidden constructs and axiom coverage

- grep of `lean_rank/LeanRank/*.lean`, `lean_rank/scripts/*.lean` and the four main-tree modules for
  `sorry|admit|native_decide|decide +native|implemented_by|@[extern|axiom|opaque|unsafe|ofReduceBool|skipKernelTC|csimp`:
  nothing. Only `set_option maxHeartbeats`/`maxRecDepth` (harmless). Certificates use `decide +kernel`
  (`twiceSquareRatios_ok` axioms `[propext]`).
- `#print axioms` covers all headline theorems: Rank.lean:168–172 (5: `two_le_finrank_of_indep`, `fg_Eplus`,
  `two_le_rank_of_theorem2At`, `rank_7_1`, `rank_twiceSquare`), Index.lean:89–92 (4: `index_two_eq`, `index_Eplus`,
  `eight_le_index_of_theorem2At`, `index_twiceSquare`). Container log `checkA_build.txt`: "Build completed
  successfully (18278 jobs)", 0 `error`, axiom lines `[propext, Classical.choice, Quot.sound]` (seen for the
  Index ones; the author's 19:17 run reports all 9). `forge.txt`: exactly 4 errors, at 19:54 (Type mismatch),
  23:43 (unsolved goals: the `n ≠ 0` slot), 27:36 (`decide` proved the membership false), 33:47 (unsolved goals:
  the `43·43 < 1828` slot) — each at the intended argument. (An older `checkA` log shows RankChecks/RankForge
  failing on a missing `LeanRank/Rank.olean` — a pre-`globs` run, superseded.)

## (4) EllipticCurves at 2164fb5

- Clone `/tmp/MichaelStollBayreuth_EllipticCurves` is at 2164fb57… (2026-11-19). README line 7: "The development is
  `sorry`-free." The only `sorry` in any `.lean` file of the repo is `comparator/Challenge.lean` (a deliberate
  comparator challenge, lib `Challenge`, not in `defaultTargets`, not imported).
- Import closure of `EllipticCurves.MordellWeil` + `EllipticCurves.Mathlib.SelmerGroup` (computed): Mathlib.Basic,
  Mathlib.FractionalIdeal, Mathlib.SIntegers, Mathlib.SelmerGroup, MordellWeil, VariableChange, WeakMordellWeil
  (7 modules, 5511 lines) — matches t7.md. grep: no sorry/admit/native_decide/implemented_by/extern/axiom/opaque/
  unsafe in them. Built from source in the same build, so `#print axioms` on `fg_Eplus` certifies the chain.
  Stoll's toolchain v4.36.0 and mathlib rev fc4f59e = lean_rank's manifest; transitive revs identical.

## (5) Proposed CI job (t7.md:212–246)

Works in outline: separate `working-directory: lean_rank`, elan with `--default-toolchain none` then the
`lean-toolchain` file selects v4.36.0 (PATH via `$GITHUB_PATH`, effective from the next step — correct), `lake exe
cache get` then `lake build`, `defaults.run.shell: bash` gives `bash -eo pipefail`. The checkout pin equals the one in
ci.yml. Estimated 10–15 min on a 4-vCPU public-repo runner (16 GB); 45-min timeout ample. (On a private-repo runner,
2 vCPU / 7 GB, slower but `public import Mathlib` is mostly file-backed mmap.)
Issues: NIT 2–5 below.

## FIX-FIRST

1. **`index_Eplus` "at every ratio with m, n, m ± n ≠ 0" is not what Lean states** — Index.lean:53–54 also assumes
   `hT : (Eplus m n).toAffine.Nonsingular (eP m n) 0`, never discharged. It is true at every such ratio (t3 has
   `nonsingular_T (hq : qe e a b ≠ 0)`, TwoDescentMap.lean:193, and `qe ≠ 0` is in `Eplus_isElliptic`'s proof via
   `qe_P`/`genH_pos`), so nothing false is proved and the 254-curve results are unaffected; but the proposed doc
   text (t7.md:261, 269, 276, 297 and t7.md:40) claims more than the theorem. Fix: drop `hT` and use
   `T := Point.some _ _ (nonsingular_T (by rw [qe_P]; …))` inside (a few lines), or reword the docs to "for any
   nonsingular `T = (e, 0)`, i.e. …".

## NIT

2. t7.md:241 `! grep -q sorryAx axioms.txt` is a **dead check**: bash `-e` never exits on a `!`-inverted command
   (tested: `/tmp/rev_t7/errexit_test.sh` continues). Redundant today because the exact-list count on the line
   above (t7.md:240) catches `sorryAx`, but write `if grep -q sorryAx axioms.txt; then exit 1; fi`.
3. t7.md:245 the forge check tests only error **positions**, not reasons; a renamed `rank_7_1` ("unknown
   identifier" at line 19) would still pass. Also grep the messages (`Type mismatch`, `unsolved goals`,
   ``decide` proved that the proposition``).
4. t7.md:227 elan via unpinned `curl | sh` from `master` breaks ci.yml's "actions pinned by commit hash" posture and
   there is no `.lake` cache (full Mathlib cache download + EllipticCurves rebuild every run). The pinned
   `leanprover/lean-action@ba3fdfe…` already used by `build` **does** have `lake-package-directory`
   (action.yml:167 at that commit, checked) and `use-github-cache` (caches `<dir>/.lake`, keyed on toolchain +
   manifest): prefer it. With no `if:`, the job runs on every push; update the ci.yml header ("Every push … runs
   `build` and `quick`") and REPRODUCE.md.
5. t7.md:139–140 "for a non-finitely-generated group `finrank` would be 0" is inaccurate (ℚ is not fg and has
   `finrank ℤ ℚ = 1`); correct reason: `finrank = toNat rank`, 0 only for infinite rank, so a lower bound is never
   vacuous.
6. TwoDescentMap.lean:784–785 docstring ("rank ≥ 2 follows by Mordell–Weil and …, which are not in Lean") becomes
   stale once `lean_rank/` lands; add a pointer to `lean_rank/LeanRank/Rank.lean` (main-tree edit, one line).
7. Optional: `lean-action`'s `leanchecker` input would re-check the oleans of the sub-project in CI.

## What holds

`rank_twiceSquare`, `index_twiceSquare`, `index_two_eq`, `two_le_finrank_of_indep` are stated as claimed, use fg and
`G[2] ⊆ {0, T}` (not torsion-freeness), call Stoll's two results with the right hypotheses, and are free of
`sorry`/`native_decide`/axioms beyond `propext, Classical.choice, Quot.sound`, Stoll's Mordell–Weil included
(sorry-free by README, grep and `#print axioms`). Negative controls are genuine (proved `forged_dependent`,
`no_indep_in_int`; 4 forged claims refused at the intended argument).

## One sentence a doc may claim

"For each of the 127 ratios `(m, n, k)` of `twiceSquareRatios`, Lean v4.36.0 with Mathlib v4.36.0 and Michael Stoll's
EllipticCurves 2164fb5 (Mordell–Weil, compiled from source; axioms `propext`, `Classical.choice`, `Quot.sound` only)
proves `2 ≤ Module.finrank ℤ E(ℚ)` for `E⁺ = Eplus m n` and `E⁻ = Eplus m (−n)` — since `E(ℚ)` is finitely
generated, rank `E(ℚ) ≥ 2` on these 254 curves — and `[E⁺(ℚ) : 2E⁺(ℚ)] = 2^(rank+1)` for every `m, n` with
`m, n, m ± n ≠ 0` once `T = (e, 0)` is given as a nonsingular point (it always is; not yet discharged in Lean)."
