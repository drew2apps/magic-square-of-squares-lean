# Hostile review of slice v3 (branch r3_v3, f03f8af..c538a1c: LemmaB, Theorem2Family, Theorem1Family, LemmaDEquiv, Theorem1Order, scripts/FamilyForge.lean)

Reviewer: separate AI agent (Claude), 25 Nov 2026, 13:49–14:10 PT. Written by an AI model; not read by a person.
Worktree not edited. Scratch: `/tmp/rev_v3/` (here: `ReviewV3.lean`, `d2search.py`; in container rhreplay:
`/tmp/rev_v3/{ReviewV3.lean,forge.log,review.log}`).

## Verdict
**Confirmed; no mathematical or formal-statement error found.** Lemma B is proved unconditionally; Theorem 2 (R, Q, T
independent mod 2E, rank ≥ 2) at every coprime twice-square pair and Theorem 1 (rank ≥ 1, R and P = 2R of infinite
order) at every coprime pair are proved with the single Prop hypothesis `LemmaD`, which is exactly Lemma D. One
FIX-FIRST in the proposed doc text (a stale count), the rest NITs.

## 1. Statements in plain mathematics (from `#check`/`#print` in the compiled environment)
- `LemmaD` (Theorem2Family.lean:35): for all x, w ∈ ℚ, w² = 2(x²+1)(x⁴−2x³+2x²+2x+1) ⇒ x = ±1. The sextic is
  2(x²+1)·H(x,1) with H = genH (G1SigmaGeneric.lean:46), i.e. UNIFORM_RANK §4's D2. The leading coefficient 2 is not a
  square, so D2 has no rational point at infinity; x = ±1 forces w² = 16. Hence LemmaD ⇔ D2(ℚ) = {(±1, ±4)}: exactly
  Lemma D, not stronger or weaker. Not vacuous: all four points (±1, ±4) satisfy the premise and are admitted (A1
  below); "D2 has no affine point" and "x = 1 only" are both refuted in Lean (A2, `control_lemmaD_sharp`). Not
  trivially false: an exact search |p|, q ≤ 1500 finds only x = ±1, and the LMFDB entry 576.a.576.1 (J(ℚ) ≅ ℤ/10, four
  points) cited by UNIFORM_RANK makes it true. Lean cannot derive False from it at an obvious point (R9).
- `lemmaD_iff` (LemmaDEquiv.lean:171): LemmaD ⇔ [for all coprime m, n ∈ ℤ with m² + n² = 2k², m² ≠ n², H(m,n) is not a
  square in ℚ]. The m² ≠ n² clause is needed (H(1, ±1) = 4). Proof checked by hand (x = p/q reduced; p + q odd gives
  2·odd = square; p, q odd gives 16·S₁H₁ = z² with S₁, H₁ coprime and positive).
- `theorem2_family`/`rank_family` (Theorem2Family.lean:133, 149): LemmaD ⇒ for all m, n, k ∈ ℤ, coprime, 1 ≤ n < m,
  m² + n² = 2k²: Theorem2At on Eplus m n and on Eplus m (−n), and 2 ≤ Module.finrank ℤ of both point groups.
  `rank_family_all` (LemmaDEquiv.lean:183): the same for all signs and orders, with m² ≠ n². The quantification is
  genuinely universal (∀ m n k : ℤ), not over a list; instantiated at (7,1), (−7,1), (1,7), (239,1) (A3, A4).
- `rank_one_family` (Theorem1Family.lean:191): LemmaD ⇒ for coprime m, n with mn ≠ 0, m² ≠ n²: finrank ≥ 1 on
  Eplus m (±n). `R_infinite_order`/`P_infinite_order` (Theorem1Order.lean:37, 103): same hypotheses ⇒ N•R ≠ O and
  N•P ≠ O for all N ≥ 1 on Eplus m n (E⁻ by n ↦ −n, allowed since n ranges over ℤ).
- Unconditional: `T_not_double` (no W with 2W = T), `no_order_four` (4W = O ⇒ 2W = O) on Eplus m n for coprime m, n,
  mn ≠ 0, m² ≠ n² (Theorem2Family.lean:175, 189). `lemmaB` / `euler14` / `pocklington` as in the docstrings.
- Rank notion: `Module.finrank ℤ (Eplus ..).toAffine.Point`, finite generation from Stoll's `fg_point_of_numberField`
  through `fg_Eplus`, and Stoll's `AddSubgroup.index_range_nsmul_of_fg`: identical to Rank.lean.

## 2. Lemma B, the descent
Read line by line (LemmaB.lean:78–226). (a² − b², ab, c) primitive (a = 2a' even, b odd); Mathlib's
`coprime_classification'` gives a² − b² = M² − N², a'b = MN, c = M² + N²; four-number lemma (Nat) |a'| = rs, |b| = tu,
|M| = rt, |N| = su; then s²(4r² + u²) = t²(r² + u²), gcd(s,t) = gcd(r,u) = 1, so r² + u² = κs², 4r² + u² = κt²,
κ | 3·gcd(r²,u²) = 3, and κ = 3 forces 3 | r, u; so (u, 2r, t) is a primitive triple, u = P² − Q², r = PQ,
P⁴ − P²Q² + Q⁴ = s², s < c = (rt)² + (su)²: strong induction on c. Correct. euler14: x, y odd gives 16(a⁴ − a²b² + b⁴)
with a = (x+y)/2, b = (x−y)/2, a + b = x odd; opposite parity uses (x² − y², 4xy, z) and M⁴ − M²N² + N⁴ = (x² + y²)².
Correct. The excluded cases are exactly the trivial solutions: an exact search (coprime, x, y < 400) finds only
(0,1), (1,0), (1,1) for d = 14 and (0,1), (1,0), (1,1) for Pocklington, each excluded by mn ≠ 0, m² ≠ n² or
a + b odd. Use: `gen_H_mul` (H(m,n)H(m,−n) = m⁸ + 14m⁴n⁴ + n⁸) then `lemmaB_rat` in theorem2_family_of_H, T_not_double,
rank_one_of_lemmaD, rank_family_all. Correct.

## 3. Grep and axioms
No sorry/admit/native_decide/axiom/implemented_by/extern/unsafe/opaque/set_option/instance/notation/macro in the
six files. `#print axioms` (my file, compiled environment) on rank_family, rank_family_all, theorem2_family,
rank_one_family, R_infinite_order, P_infinite_order, T_not_double, no_order_four, lemmaD_iff, lemmaB, euler14,
pocklington: [propext, Classical.choice, Quot.sound] for each; LemmaD is a `def … : Prop`, used only as a hypothesis.

## 4. Compilation (container rhreplay, /work/mss/v3, Lean v4.36.0 via /work/vol/elan)
Sources in /work/mss/v3 are sha256-identical to the worktree (all 6 LeanRank files, FamilyForge, the 4 main-tree
modules). `lake build --no-build` on the five new modules: "All targets up-to-date", so the oleans are those of these
sources (nothing rebuilt; shared packages untouched). Peak RSS ~4.9 GB per `lake env lean` run.
- `scripts/FamilyForge.lean`: exit 1, exactly 8 errors at lines 20, 24, 29, 34, 42, 46, 64, 71 (examples 1–6, 9, 10);
  examples 7, 8, 11 accepted. Matches the report.
- My `ReviewV3.lean`: refused R1 (rank_family_all with IsCoprime deleted), R2 (rank ≥ 2 at non-coprime (14,2)),
  R3 (rank ≥ 2 at square-type (12,5), no k), R4 (rank ≥ 3 at (7,1)), R5 (rank ≥ 2 from Theorem 1 alone), R6 (Lemma B
  at (1,−1)), R7 (euler14 at x = y = 1), R8 (−1 = 1 from LemmaD at (−1,4)), R9 (False from LemmaD): 9 errors at lines
  43, 48, 54, 61, 67, 71, 75, 81, 84, one per block. Accepted A1 (the four points admitted), A2 (LemmaD is not
  "no point"), A3 ((7,1), (−7,1), (1,7)), A4 ((239,1)), A5 (T_not_double at (4,1), where H(4,1) = 13², no LemmaD).

## 5. Prior art
Mathlib as pinned (v4.36.0): only `not_fermat_42` (a⁴ + b⁴ = c², FLT/Four.lean); no x⁴ − y⁴ = z², no Pocklington
quartic, no d = 14 (grep; the only "Pocklington" is a TODO about the primality test, Tactic/Primality/Producers.lean:24).
Web (short): the Lean Zulip thread "Fermat's right triangle theorem" (Sep 2020, "Is there code for X?") talks about a
descent for x⁴ − y⁴ = z², not d = 14 or Pocklington; nothing else found. Not exhaustive; upstream Mathlib master and
the EllipticCurves repo past 2164fb5 were not checked.

## Findings
- **FIX-FIRST (doc text)** reports/v3.md:95: the proposed UNIFORM_RANK §7 paragraph says FamilyForge has "six
  refused, two controls accepted"; the committed script has 8 refused, 3 accepted (verified above). Fix before pasting.
- NIT FamilyForge.lean:31–34 (ex. 4): Theorem2At 14 2 10 is true (Eplus 14 2 is Eplus 7 1 scaled by λ = 2:
  X ↦ λ⁶X, y ↦ λ⁹y, and R, Q, T scale the same way); :66–71 (ex. 10): the hypothesis hR is uninhabited (at (1,−1)
  R = (0,0), the singular point of y² = (X+16)X²), so the claim is vacuously true; :44–46 (ex. 6) only shows that a
  lemma proving IsSquare does not prove its negation. Several "forged claims" therefore test that the hypotheses
  cannot be discharged, not that a false statement is refused. Doc wording: "inputs outside the hypotheses are
  refused", not "forged (false) claims refused". (R1/R2 above are of the same kind: true by scaling.)
- NIT Theorem1Order.lean:38 and Theorem2Family.lean:176: `R_infinite_order` takes the nonsingularity proof hR, and
  `T_not_double` takes hT, as hypotheses. They are provable at every admissible ratio (`R_add_R`; `nonsingular_T`
  with `qe_ne_zero`), so nothing is vacuous, but an ∃ form like `P_infinite_order` would make that visible.
- NIT LemmaB.lean:11–12, 228–230: `pocklington` is only the opposite-parity case (Pocklington's theorem also covers a,
  b both odd with a² ≠ b²); enough for euler14. The proposed doc text states it correctly; keep it that way.
- NIT the proposed header status line ("Theorem 2 and Theorem 1 are in Lean for every ratio") should say "for every
  coprime pair, on the models Eplus m (±n)": coprimality (the scaling reduction) and E⁻ = Eplus m (−n) (A2, sympy)
  are not in Lean. The §7 paragraph does say this.
- NIT deprecations `if_pos`/`if_neg` (Theorem1Family.lean:161,163; Theorem1Order.lean:59,62), cosmetic.
- Environment note: in rhreplay `sh -lc 'lake …'` resolves to the /root/.elan shim and starts downloading a v4.36.0
  toolchain; call /work/vol/elan/toolchains/leanprover--lean4---v4.36.0/bin/lake. My first attempt did this; I killed
  those 5 PIDs by number and removed the partial toolchain, the lock and the 609 MB download in /root/.elan
  (settings.toml and known-projects there were touched and left).

## One sentence a doc may claim
"Lean (standard axioms only) proves Lemma B for every coprime pair, and proves that Lemma D, which is stated as the
proposition that every rational point of D2 has x = ±1 and is shown equivalent to H(m, n) never being a square at a
coprime twice-square pair with m² ≠ n², implies Theorem 2 (R, Q, T independent in E(ℚ)/2E(ℚ), rank ≥ 2 on
Eplus m (±n)) at every coprime twice-square pair and Theorem 1 (rank ≥ 1, R and 2R of infinite order) at every coprime
pair with mn(m² − n²) ≠ 0; Lemma D itself and the identification of Eplus m (−n) with E⁻ are not in Lean."
