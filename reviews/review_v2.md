# Review of slice v2 (branch r3_v2, diff f03f8af..3fbed44), hostile-referee pass

Written by an AI model (reviewer); not read by a person. 22:42-22:55 UTC 2026-11-25.
Scope: lean_rank/LeanRank/RankZero{,Points,Gen,GenCurve,Census}.lean,
lean_rank/scripts/{RankZeroForge.lean, RankZeroCensusForge.lean, rank_zero_check.py, gen_rank_zero_certs.py}.
Scratch: /tmp/rev_v2/ (pycheck.py, chk_census.py). Worktree not edited.

## Verification limits (read first)
- NO Lean rebuild by the reviewer: this container has no Lean toolchain and the proxy refuses
  GitHub release downloads (elan: CONNECT 403). So "compiles" and `#print axioms` = [propext,
  Classical.choice, Quot.sound] rest on the author's log only (container /work/mss/v2, with a
  lakefile WITHOUT the MagicSquareOfSquares lib). Static review + independent Python/PARI below.

## (1) Statements
- `RankZero.finrank_eq_zero : Module.finrank ℤ E.Point = 0` (RankZero.lean:604), `E = ⟨0,-28,0,192,0⟩ : Affine ℚ`
  (RankZero.lean:372) = y² = x³ − 28x² + 192x = x(x−12)(x−16). Closed theorem, no hypotheses. Point
  group is Mathlib/Stoll `Affine.Point` over ℚ. OK.
- `finrank_eq_zero_of` (RankZeroGenCurve.lean:169): curve `⟨0, −(e₂+e₃), 0, e₂e₃, 0⟩` = x(x−e₂)(x−e₃). Every
  hypothesis is a decidable side condition (primes lists prime and covering e₂e₃ resp. e₂(e₂−e₃),
  `check = true`, `TablesOK`, `#good ≤ 4`, (1,1) ∈ good, three RealizesTriple witnesses); none is
  unsatisfiable, and all four census instances discharge them by decide/norm_num (RankZeroCensus.lean:74-191).
  No hypothesis makes it vacuous.
- Soundness of the certificate does not depend on the certificate being right: `Cert.good_of_check`
  (RankZeroGen.lean:406) is a theorem; `excludedAt_eq_false` proves that an actual integer solution
  (with ¬(p | n ∧ p | e), from IsCoprime n e) makes some residue pair pass, so "excluded" can only mean
  "no local solution mod M" (correct direction). `realExcluded_eq_false` likewise for the real place.
  Coverage of all square classes is `descent_gen` (every (d₁,d₂) lies in classesOf S₁ × classesOf S₂
  given the coverage hypotheses), and `Cert.check` iterates over exactly that product. Not every bad
  prime needs a table (exclusion is only sufficient); the forge confirms the coverage hypothesis is
  load-bearing (S₂ without 7 passes `check` but cannot discharge `hd₂`).
- Independent Python re-implementation of `Cert.check` (Lean semantics, /tmp/rev_v2/pycheck.py):
  all four census certificates and E′ (generic form, table mod 64) leave no non-good class unexcluded;
  no good class is excludable; sq49/sq64 tables match. So the `decide +kernel` goals are true.

## (2) Rank bound
- Proof chain: Ψ∘μ lands in Tset good (mem_T / hmemT), Ψ injective (Lagrange at the 3 roots), so
  #im μ ≤ #good ≤ 4; Stoll `card_range_μ` (SelmerGroup.lean:967): #im μ = 2^finrank · #E[2];
  `card_ker_nsmul_two` + 3 roots gives #E[2] = 4; hence 2^r ≤ 1. Correct.
- E′ consistency: im μ = {(1,1),(3,−3),(2,−1),(6,3)} = images of O/(16,0), (0,0)/(12,0), (8,16), (24,48);
  E′(ℚ) ≅ Z/2×Z/4, #E/2E = 4. Matches goodClasses.

## (3) Census identification
- scripts/g1_full/global.py builds `w*f0*f8` as ellinit([0, m⁴+n⁴, 0, m⁴n⁴, 0]) (X = m²n²w) -- the SAME
  Weierstrass model as `curve (−m⁴) (−n⁴)`. PARI (chk_census.py): for (4,3),(7,1),(12,5),(24,7) the census
  curve w*f0*f8 has rank [0,0], torsion 16/8/16/16, and is a closing curve; for (24,7) it is the ONLY
  closing curve. (15,8): rank [1,1], correctly refused. None of the four pairs is in open_pairs_1000.txt.
- So the Lean theorems replace exactly the PARI ellrank step of the census for those pairs on that curve;
  the torsion enumeration / "no torsion w gives a solution" step is still PARI.

## (4) Escape hatches
- grep sorry/admit/native_decide/ofReduceBool/implemented_by/axiom/opaque/unsafe/skipKernelTC/macro:
  none in the 5 Lean modules or 2 forge files. Only `decide +kernel`, plain `decide`, norm_num.
- `#print axioms` lines present in each file; not re-run (see limits).

## (5) Stoll's library (/tmp/MichaelStollBayreuth_EllipticCurves at 2164fb5 = manifest pin)
- No sorry/admit/axiom/opaque/implemented_by/extern/native_decide in any .lean (only a comment mentions
  native_decide). Used API: fg_point_of_numberField, μ, card_range_μ, card_ker_nsmul_two, μX_of_eval_f_*,
  redHom, eq_zero_of_isOfFinAddOrder_of_red_eq_zero, Point.mapEquiv. Statements read: as used.

## FIX-FIRST
- None in the Lean code. Precondition before any doc says "kernel-checked": someone re-runs
  `lake build LeanRank.RankZeroPoints LeanRank.RankZeroCensus` with the WORKTREE lakefile and the
  `#print axioms` lines (the author built with a modified lakefile; reviewer could not build).

## NIT
- RankZeroCensus.lean:9-11: the identification w·f0·f8 ≅ Y² = X(X+m⁴)(X+n⁴) is only a docstring; a one-line
  `ring` lemma (y²=w f0 f8 → (m²n²w, m³n³y) on curve (−m⁴)(−n⁴)) would put it in Lean. Low risk:
  the census's own PARI model is literally this curve.
- finrank = 0 is also the value for non-finite modules; non-vacuity rests on Stoll's Mordell–Weil
  (sorry-free). A corollary `Finite (curve ..).Point` for the census curves (as done for E′, E_B) would make
  the statement self-evidently strong.
- Report v2.md lists (17,7) among produced certificates; it is not in RankZeroCensus.lean (stated later; fine).
- Remaining linter/deprecation warnings (dif_neg) per author.

## What a doc may claim
"In Lean (lean_rank/LeanRank/RankZero*.lean, Stoll's EllipticCurves 2164fb5, no sorry/native_decide; author
build log, standard axioms): rank E′(ℚ) = 0 and #E′(ℚ) = 8 for y² = x(x−12)(x−16), Euler's x⁴+14x²y²+y⁴=z²
and Lemma B, and rank 0 for Y² = X(X+m⁴)(X+n⁴) -- the model global.py uses for w·f0·f8 -- at (4,3), (7,1),
(12,5), (24,7); the torsion step of those census closures is not in Lean."
