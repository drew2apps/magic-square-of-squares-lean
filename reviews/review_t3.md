# review_t3 — hostile referee of `sess_t3` (TwoDescentMap / TwoDescentCertificates)

Written by an AI model (reviewer, did not write the code); not read by a person. Worktree
`/home/claude/conductor/t3`, head 684f663, base 8214be3. Scratch under `/tmp/rev_t3/`.

## Log

- Read `TwoDescentMap.lean` (955 lines) in full, the head/tail and list of `TwoDescentCertificates.lean`,
  `two_descent_certs.py`, `scripts/TwoDescentChecks.lean`, `UNIFORM_RANK.md` §1–§7, `sigma_exact.model`,
  the `gen*` definitions of `G1SigmaGeneric.lean`.
- Generator run (check mode, 3 s): matches the committed Lean file byte for byte; transcript matches
  `two_descent_certs.txt`; worktree left clean.
- Independent enumeration (one-liner, `m ≤ 1000`, `1 ≤ n < m`, gcd 1, `m² + n² = 2k²`): 127 ratios, equal
  as a list to `twiceSquareRatios`. The 254 theorem names, their `(m, ±n, k)` and the 127 pairs used in
  `theorem2_twiceSquare` agree with it.
- Independent cross-check (`/tmp/rev_t3/scr/xcheck.py`, my own transcription of `eP/aP/bP` and of `curve`):
  `(a₂, a₄, a₆)` of `Eplus m (±n)` = `sigma_exact.model(m, n, ±2)` on all 254 curves and on ~200 random
  coprime ratios (not twice-square); `model(7,1,−2) ≠ Eplus 7 1` (the check can fail); `R`, `Q`, `T` on the
  model on all 254; every bracket number `N` and root `s` of the 254 certificates recomputed (`H(m,n)`,
  `H(m,−n)`, product, `2mn`, `−2mn`, in `theorem2_of`'s order; negatives exactly where `N < 0`).
- Grep of the two library files and `scripts/TwoDescentChecks.lean`: no `sorry`, `admit`, `native_decide`,
  `axiom`, `implemented_by`, `extern`, `unsafe`, `opaque`, no `instance`/`attribute`/`macro`. One
  `decide +kernel` (`twiceSquareRatios_ok`: 127 Nat triples, gcd and products ≤ ~2·10⁶, kernel GMP).
  One `set_option maxRecDepth 4000`.
- Symbolic check (`/tmp/rev_t3/scr/sym.py`, sympy, `sigma_exact.model`'s own `pmul/peval/pder` re-run on
  symbols): `(a₂, a₄, a₆)` of `model(m, n, +2)` and `model(m, n, −2)` minus those of Lean's `curve eP aP bP` at
  `(m, +n)` and `(m, −n)` expand to `[0, 0, 0]` for general `m, n`. So "`E⁻` at `(m, n)` is `Eplus m (−n)`" holds
  identically (outside Lean); inside Lean it is pinned at `(7, 1)` by my `P5` below.
- Lean, Mac slice a1 (container sources sha256-identical to the worktree: TwoDescentMap a90dc126…,
  TwoDescentCertificates 29766c39…, G1SigmaGeneric 1a5aad40…, toolchain v4.33.1, Mathlib 0df444a):
  `lake build …TwoDescentCertificates` (traces match, replayed); fresh `lake env lean` of
  `TwoDescentMap.lean` (6 s, 0.24 GB) and of `TwoDescentCertificates.lean` (9 s, 0.50 GB): exit 0, no
  errors/warnings; `#print axioms` all `[propext, Classical.choice, Quot.sound]` except
  `control_bracket` `[propext, Quot.sound]`, `twiceSquareRatios_ok` `[propext]`. `scripts/TwoDescentChecks.lean`:
  exit 0.
- Forge run, `scripts/ReviewForge.lean` in the slice host dir (copy `/tmp/rev_t3/ReviewForge.lean`, sha256
  79c538ac…), one run, 6 s. Positive controls compiled: P1 verbatim copy of `t2_991_263_p`; P2 `rv_teeth`
  (with `Q` replaced by `R`, the `Theorem2At` predicate is *false* on `E⁺(7,1)`: `1·R + 1·R + 0·T = R + R`);
  P3 `2·T = O`; P4 `Theorem2At 7 1 5` is `rfl`-equal to the expected ∃∀ statement; P5
  `(Eplus 7 1).a₁ = a₃ = 0`, `a₂ = 300176`, `a₄ = 20572160000`, `a₆ = 2007040000000000` (the numbers of
  `sigma_exact.model(7,1,2)`). Forgeries, all refused: F1 last digit of `H·H⁻` 6→7 (l.22); F2 its root
  …823→…824 (l.31); F3 root of `H` 793256→793255 (l.38); F4 `2mn` 521266→521267 (l.50); F5 `H`, `H⁻` swapped
  (l.56, 57); F6 `k = 724` (l.64); F7 list check with `(23, 7, 16)` (`decide` proved it false); F8 bracket for
  `300² = 90000` (l.77). `AddSubgroup.index_range_nsmul_of_fg` and `WeierstrassCurve.IsCharNeTwoNF`: unknown
  constants in this Mathlib.

## Step 1 — statements

`curve e a b` = Weierstrass `a₁ = a₃ = 0, a₂ = a − e, a₄ = b − ae, a₆ = −eb` (no shift needed: `X − e` is used
directly), equation `y² = (X − e)(X² + aX + b)` (`curve_equation_iff`). `Eplus m n = curve (−2S³)
(16m²n²(m+n)²) (16m²n²(m+n)⁴S²)` over ℚ; = `model(m,n,2)`, and `Eplus m (−n)` = `model(m,n,−2)` (symbolic, above).
`Theorem2At m n k` (l.769): ∃ proofs that `R = (xR, yR)`, `Q = (0, yQ)`, `T = (e, 0)` are nonsingular on
`Eplus m n`, and ∀ `c₁ c₂ c₃ : ℤ`, ∀ `W : (Eplus m n).toAffine.Point` (Mathlib's group), `c₁R + c₂Q + c₃T = W + W
→ all cᵢ even`. That is exactly F₂-independence of the images in `E(ℚ)/2E(ℚ)` (stronger than the {0,1} form, and
covers differences, so the 8 combinations are distinct mod 2E). Not vacuous: P2 shows the same predicate is
false for a dependent triple; `2E` is `{W + W}` over all rational points, so membership is not trivially false
(P3, `R_add_R`, `control_7_1`). `theorem2_of` (l.787) hypotheses: `m²+n² = 2k²`, `m, n, m±n ≠ 0`, and
`¬IsSquare` of `H(m,n)`, `H(m,−n)`, `H·H⁻`, `2mn`, `−2mn` — all satisfiable (they are, 254 times).
`theorem2_twiceSquare`: ∀ `p` in the 127-list, `Theorem2At p.1 p.2.1 p.2.2 ∧ Theorem2At p.1 (−p.2.1) p.2.2`.
This is UNIFORM_RANK Theorem 2's displayed conclusion, minus "hence rank ≥ 2".

Q step re-derived: on `y² = (x − e)((x − r)² + s²)`, with `θ = r + is`, `x(2W) − θ = (u + iv)²`,
`u + iv = ((x−r)² + s² − 2is(x − e))/(2y)` (`gauss_coord` k1, k2 are exactly `4y²(x(2W) − r) = U² − V²`,
`4y²(−s) = 2UV`). At `x(Q) = 0`: `(u²+v²)² = r² + s² = b`, so `u² + v² = ±4mn(m+n)²S`, giving `u² = 2mn(m+n)⁴`
or `u² = −2mn(m+n)²(m−n)²`, i.e. `±2mn` a square (`gauss_family`). Equivalent to UNIFORM_RANK's "`i·mn` a
square in ℚ(i) ⇒ `mn = ±2r²`". Correct; `u, v ∈ ℚ` because `W ∈ E(ℚ)`.

## Step 2 — α

`alpha` (l.214): `O ↦ 1`, `(e, y) ↦ q(e)`, else `x − e`. `alpha_add` (q(e) ≠ 0): `α(P₁+P₂)α(P₁)α(P₂)` square for
all points (all cases of Mathlib's `add_some`/`add_of_Y_eq`); `alphaHom : E(K) →+ Additive (Kˣ ⧸ (powMonoidHom
2).range)`; `alphaHom_two_smul`, `alpha_double`: 2E ⊆ ker (only ⊆; enough). The independence proof uses
`alpha_add`/`alpha_double` directly, not `alphaHom`.

## Findings

FIX-FIRST: none in the Lean or the generator.

Integration (author already lists it): merging the two files without `declared_theorems.txt` 655 → 961
breaks CI `quick` (ci.yml:33 counts `^theorem` in `MagicSquareOfSquares/*.lean`; 50 + 256 = 306 new, checked).

NIT:
1. two_descent_certs.py:151–155 negative #5 never calls generator code (it re-tests `isqrt` arithmetic);
   #1, #2 (l.130–139) only hit the `m²+n² = 2k²` assert. No negative for `crosscheck` (e.g. assert
   `model(7,1,−2) ≠ Eplus 7 1`, which I checked holds). Lean F1–F8 above are the real negatives.
2. two_descent_certs.py:170–172: the cross-check compares `model` with a Python *transcription* of
   `eP/aP/bP`, not with Lean's terms. Lean's own link is `Eplus_cubic` (cubic = `genCub`); add a Lean `example`
   pinning coefficients (like my P5) or cite the symbolic check.
3. TwoDescentMap.lean:738, :784 and the proposed doc text say `E⁺(ℚ)[2] = {O, T}`; Lean proves only `⊆`
   (affine `P + P = O ⇒ P = T`). True and trivial the other way, but say "is contained in".
4. Completeness of the 127-list is not in Lean (stated honestly, Certificates header l.15–17).

## Step 6 — prior art (Stoll, EllipticCurves, clone `/tmp/MichaelStollBayreuth_EllipticCurves`, 2164fb5, v4.36.0)

Stoll has the full `x − T` map `μ : Multiplicative W.Point →* Units.modPow (K[X]/(f)) 2`
(WeakMordellWeil.lean:452 `μ₀`, :636 `μ`), `μ₀_two_nsmul` (2E ⊆ ker) and the stronger `ker_μ_eq` (:762,
ker = 2E), `Point.exists_eq_f_sub_sq_of_add_add_eq_zero` (line through three points), under
`[W.IsElliptic] [W.IsCharNeTwoNF] [DecidableEq K]`. Duplicated by TwoDescentMap.lean: `alpha`, `alpha_add`,
`alpha_double`, `alphaHom`, `alphaHom_two_smul`, `alphaHom_eq_zero_iff`, `SqClass` = μ projected to the factor
`K[X]/(X − e) ≅ K` of `K[X]/(f)`; `double_isSquare`, `gauss_coord`, `gauss_double` = `μ₀_two_nsmul` projected to
`K` and to `K[X]/((X−r)² + s²)` (≅ `K(i)`), written over `K`; `chord_*`, `tangent_coord` = the three-point
identity. Not duplicated: the family (`Eplus`, `R`, `Q`, `T`, `R_add_R`, `Eplus_isElliptic`,
`two_torsion_*`), `indep_of_alpha*`, `theorem2_of`, the certificates. (This file's map needs only `q(e) ≠ 0`,
any characteristic, no AdjoinRoot.)

Rank ≥ 2 by combination: yes, in a v4.36.0 project. `fg_point_of_numberField` (MordellWeil.lean:116; ℚ is a
`NumberField`; needs `IsElliptic` = `Eplus_isElliptic`) + Mathlib's `AddSubgroup.index_range_nsmul_of_fg`
(index of 2G = 2^rank·#G[2], as used in `card_range_μ`, SelmerGroup.lean:967; absent from Mathlib 0df444a,
confirmed) + `two_torsion_Eplus` (#E[2] = 2) + `Theorem2At` (8 distinct classes) give
8 ≤ 2^rank·2, so `2 ≤ Module.finrank ℤ (Eplus m (±n)).toAffine.Point` on the 254 curves. Needs either porting
TwoDescentMap + G1SigmaGeneric (elementary) to v4.36.0 or a bridge project on v4.36.0 requiring Stoll's repo;
back-porting Stoll to v4.33.1 is not realistic. `μ` itself is not needed for this step.

## What a doc may claim (one sentence)

"For each of the 127 coprime ratios `n < m ≤ 1000` with `m² + n² = 2k²` and each of the curves `Eplus m n`,
`Eplus m (−n)` (which are `sigma_exact.model(m, n, ±2)`, an identity checked by sympy, not in Lean), Lean
proves with standard axioms that every integer combination `c₁R + c₂Q + c₃T` lying in `2E(ℚ)` has all `cᵢ` even,
i.e. `R, Q, T` are independent in `E(ℚ)/2E(ℚ)`; that the list is complete, that `rank ≥ 2`, and Lemmas B and D
are not Lean statements."
