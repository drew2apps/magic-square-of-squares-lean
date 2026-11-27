# z7 — hostile referee of roadmap H, H7 (int7/experiments/roadmap-h/h7/)

Started: Fri Nov 27 13:19:31 UTC 2026
Written by an AI model; not read by a person.


## 1. Own complete 2-descent (written before reading rank_second.py or the bundle's scripts)

`/tmp/z7/my2desc.py`: Sel² = {d ∈ Q(S,2)² : res_v d ∈ W_v for v ∈ S ∪ {∞}}, W_v = image of the local Kummer map
P ↦ (x−e1, x−e2) computed by sampling points of E(Q_v) (exact rational x, f(x) a square in Q_v tested exactly) until
the span reaches the theoretical order 4 (odd p), 8 (p = 2), 2 (∞). Sanity: on y² = x³ − n²x for 31 values of n,
dim Sel² − 2 = R + s of PARI 2.17.4's ellrank in every case (e.g. n = 17: 4, PARI [0,0,2]; n = 205: 5, PARI [1,1,2]).

Curve, derived symbolically (`/tmp/z7/velu_check.py`, sympy, exact): from dm_g4k.g4k_curve's formula
(x = c3·s, y² = x³ + c2x² + c1c3x + c0c3², P(s) = g4(s)K(s)), e = −c3 is a root (residual 0); the shift X = x + c3 gives
y² = X(X² + aX + b) with a = −32m²n²S, b = 64m²n²S⁴; 16b = (32mnS²)², a square for every (m,n); Vélu's map
(x,y) ↦ (y²/x², y(b−x²)/x²) to Y² = X(X² − 2aX + a² − 4b) checked as a polynomial identity (residual 0), and
E′: Y² = X(X − 16Smn(m−n)²)(X + 16Smn(m+n)²). So E′ always has full rational 2-torsion (claim holds, for all m,n).

My dim Sel²(E′) (and the same on the model with e/4): (12201,7360) 3; (7735,4968) 3; (6881,6233) 7; (13348,3795) 5;
(6557,324) 5; controls (11165,7548) 4, (12204,9853) 6. All agree with REPORT_H7's table. For (12201,7360) and
(7735,4968) this alone gives rank ≤ 3 − 2 = 1 with no pairing. For (6881,6233) and (13348,3795) the bound 1 needs
the pairing (rank 4 and 2).

## 2. Reruns under PARI 2.17.4 (copy of int7 in /tmp/z7/repo, MSS_REPO/H1_DIR set there, python3 -I)
- h7_g4k_rank1.py: `# matches` (7.6 s); --negative: `# matches` (14.5 s).
- h7_e04_rank.py: `# matches` (65 s); --negative: `# matches` (1.3 s).
So the committed outputs (made with PARI 2.19.1) reproduce byte-for-byte with CI's PARI.

## 3. What was reviewed about the pairing (rank_second.py header, reviews/review_t2.md)
DFull.pairing / full_selmer2 were written for the curves E0L (y² = (x−e1)(x−e2)(x−e3), full 2-torsion) and reviewed
in review_t2.md: formula re-derived, independent re-implementation reproduced the 7×7 matrix of E0L(920,129)
(pairing rank 2) and the trivial matrix of E0L(72,65); dim Sel² cross-checked with ell2cover. So E′ (full 2-torsion,
e's integers) is the same setting. Differences: (i) bound_full (H1, bundle) checks only symmetry and alternation;
it does NOT do e0l_rank's two extra checks (point images in the kernel; same matrix with tangents at other conic
points, alt=True); (ii) the reviewed cases had pairing rank 0 and 2; (6881,6233) needs rank 4 (dim 7) and
rank0_ct's curve rank 4 (dim 6); never reviewed at that rank. (iii) roots_of returns sorted roots, so e1 is the least
root, not 0 as in e0l_rank; the formulas in DFull are written for general (e1,e2,e3), so this is not a problem by
reading, but it is a different input shape from what was reviewed.

## 4. Own Cassels–Tate pairing (claim 1)
`/tmp/z7/ct_mine.py`, written from my own derivation (below), with review_t2's *approach* (PARI p-adics to p^80,
pari.hilbert, conic points from qfsolve moved along 2 random lines so tangent points differ from rank_second's;
4 local x-values per place × all 8 sign choices; extra places 3, 5, 7, 11, 13 asserted to give +1), not its code
and not rank_second's code. Real place by mpmath at 120 digits (a cross-check, not a certificate).
Derivation: z_i ∘ (translation by T) = e₂(T,T_i)·z_i gives T1 = (u,−v,−w), T2 = (−u,v,−w); reflections σ_u, σ_w
compose to T2, so F_u − F_w ~ T2, F_v − F_w ~ T1; div(l_i/l_j) = 2(F_i − F_j); b = b2⊗T1 + b1⊗T2 for Kummer
coordinates (b1,b2) = (e₂(b,T1), e₂(b,T2)); so ⟨b, D_a⟩ = Π_v (b2, l_v/l_w)_v (b1, l_u/l_w)_v, the formula of
rank_second (line 837). Places outside {∞, 2, primes of a, b, Δ}: (b, g) is unramified on the smooth Z_p-model
(b a unit, l primitive so v_special(g) = 0, div g even), so the factor is 1.
Sanity: y² = x³ − n²x, n = 17, 73, 205, 34, 41, 113, 5: pairing ranks 2, 2, 2, 0, 0, 0, 0 = PARI's s ([0,0,2],
[0,0,2], [1,1,2], [2,2,0], [2,2,0], [0,2,0], [1,1,0]).

Results on E′ (same roots as h7_g4k_rank1.txt; my own Sel² basis):

| ratio | my dim Sel² | my pairing rank | bound | images of G (Vélu) and T in kernel | PARI 2.17.4 ellrank on E′ |
|---|---|---|---|---|---|
| (12201,7360) | 3 | 0 | 1 | yes | [1,1,0] |
| (7735,4968) | 3 | 0 | 1 | yes | [1,1,0] |
| (6881,6233) | 7 | 4 | 1 | yes | [1,1,4] |
| (13348,3795) | 5 | 2 | 1 | yes | [1,1,2] |
| (6557,324) | 5 | 0 | 3 | – | [1,3,0] |
| (11165,7548) | 4 | 0 | 2 | – | [2,2,0] |
| (12204,9853) | 6 | 2 | 2 | yes (PARI found 2 points on C this time) | [0,2,2] |
| (13899,1180) | 4 | 0 | 2 | yes | [2,2,0] |

Every number agrees with h7_g4k_rank1.txt, h7_g4k_rank1_negative.txt and h7_b13000_rank1.txt. Point of infinite
order: R.mul is the textbook chord-tangent law on y² = x³+ax²+bx (read), nG ≠ O for n ≤ 12 suffices by Mazur. 
**Claim 1 holds** (computed, exact; resting on quoted theory: Vélu — which I also checked as a polynomial identity —,
isogeny invariance of rank, Mazur, soundness of the 2-descent and the CT pairing formula).
Gap (minor): bound_full does not itself check the point images in the kernel or tangent-choice independence
(e0l_rank does); I did both checks independently (images: yes; my tangent points differ from rank_second's and
give the same ranks). Note: PARI's ellrank is randomised; the negative file's "(12204,9853): PARI finds 0
independent points" is run-state dependent (in my run PARI found 2 on C). Harmless for the control.

## 5. Claim 2 (F4's chain; condition exactly "rank 1")
g4k_exact.py header and output: "IF rank C(Q) = 1, every solution t is degenerate", C = dm_g4k.g4k_curve(m,n)[0]
(line 147), the same curve as h7_g4k_rank1.py's `shifted()` (I re-derived the shift: a = 3e + a2, b = 3e² + 2a2e + a4,
e = −c3 = −2S³, exact). review_x3.md: "F4 DM step: HOLDS ... Every step after 'rank C(Q) = 1' re-derived"; the only
open condition was rank 1 (plus quoted Silverman/CPS, tier E). Torsion ℤ/2 is proved inside the chain (not assumed).
Rerun under PARI 2.17.4 of g4k_exact.py on 12201,7360 and 13348,3795: identical to h7_g4k_exact.txt and to the first
8 lines of h7_b13000_g4k_exact.txt. F4's own g4k_exact.txt has the conclusion for (7735,4968), (6557,324), (6881,6233).
**Holds.** With claim 1, (12201,7360), (7735,4968), (6881,6233), (13348,3795): every solution degenerate, resting on
quoted theorems only (Vélu, isogeny invariance, Mazur, 2-descent/CT soundness, Silverman, CPS Thms 2 and 15).
Note: g4k_exact also prints "IF rank 1 ... degenerate" for (13899,1180), whose C has two independent points
(PARI [2,2], my bound 2): there the condition is false; the report correctly leaves it open.

## 6. Claim 3 (partner)
((6881+6233)/2, (6881−6233)/2) = (6557,324), gcd(13114,648) = 2 (both odd: halved partner). partner_forms /
partner_forms_half (G1Relations.lean, kernel-checked) give the five forms. I checked symbolically (sympy,
/tmp/z7/partner_cells.py) that all nine cells of (M/2,N/2) at t' = (t+1)/(t−1) are the nine cells of (m,n) at t,
permuted [0,2,1,3,7,8,6,4,5], times the one factor 1/(t−1)² (4/(t−1)² unhalved); so SOLUTION/degenerate/no are
preserved; t = 1 ↔ t' = ∞ (both degenerate: repeated cells). CENSUS9000.md already moves exclusions this way
("73 through their partner"; (7393,5783), (7409,2737) through partners). **Holds**, with the right strength as
long as it is stated as "(6557,324) is excluded through its partner (6881,6233)", same tier as (6881,6233), and not
as "rank 1 for (6557,324)'s curve" (its own E′ bound is 3: my Sel² 5, pairing 0; PARI [1,3,0]).

## 7. Claim 4 (E04 of tier B)
- Counts: dm_evidence.txt has 91 E04 lines (of 122: 91 E04, 25 D1, 6 E0L), 66 distinct run pairs (field 2), and
  66 distinct (run pair, G) — so 66 curves. Holds.
- Model: dm.py line 47 `EllipticCurve(QQ, [0, 2*S*S, 0, 4*S*S*m*m*n*n, 0])` with (m,n) the pair it was run on
  (field 2 of dm_evidence = "pair for which the curve was run"); h7_e04_rank.run_curve uses field 2 and the same
  a-invariants; roots 0, −2Sm², −2Sn² asserted. Holds.
- Generator: the line's G abscissa is checked to be on that curve (non-zero square), nG ≠ O for n ≤ 12. Forgery:
  the order-4 point x = −186116 of E04(23,7) (torsion ℤ/4×ℤ/2) is refused "the point has finite order"; the
  committed negative only tests x = 0, which is refused earlier (y = 0), not by the order test. Holds.
- My own descent + pairing (ct_mine.py) on all 66 E04 curves, on E04 itself: (dim, rank) = (3,0) ×46, (5,2) ×16,
  (7,4) ×2, (5,0) ×2; G's image and T's in the kernel on all 66. Per curve identical to h7_e04_rank.txt's first
  curve, except (320,289), where rank_second's local search gives up on E04 and mine gets (7,4), bound 1 on E04
  itself (the bundle gets 1 on the isogenous curve; consistent).
- Isogeny classes (PARI ellisomat as a finder; my descent on each full-2-torsion member): (320,289), (1015,961),
  (1015,713) each have 6 curves, of which exactly 2 have full rational 2-torsion. (1015,961): E04 (5,0) → 3, the
  2-isogenous one (7,4) → 1 (PARI [1,1,4]): holds. (1015,713): E04 (5,0) → 3, the other (7,2) → 3: neither gives
  1, so "(864,151) not covered" holds. (PARI gets 1 on a member with one 2-torsion point, [1,1,2]; out of reach here.)
- **Wording error:** "On the curve of (1015,713) all 6 curves give 3, or the pairing's local search gives up".
  The walk (least_bound) dedupes only by translation, not by scaling x → u²x, so its "6 curve(s)" are 2 curves
  each visited 3 times (I replayed the walk with the bundle's roots_of/isog/shift and compared minimal models:
  sequence [0,1,0,1,0,1]; same for (1015,961), (320,289) and (15,8)). So rank_second's pairing was never completed
  on the second curve of (1015,713) (I computed it: bound 3). Soundness unaffected (duplicates only waste slots).
- Forgery "rank ≥ 1 never gets bound 0": the bundle's least_bound on 24 random E04 curves (m ≤ 400) and 9
  congruent-number curves (n = 5,6,7 rank 1; 34,41,65,137,210 rank 2; 1254 rank 3): never below the exact lower
  bound from PARI's points (their images in Sel²), always equal to PARI's R; n = 5,6,7 give 1. Note the script does
  not assert bound ≥ 1 when G is known (a bound 0 would print "rank <= 0" and fall into "not covered", not abort);
  harmless here, worth an assert.
- Reruns match (section 2). **Claim 4 holds**, with the wording fix above. Status: the 90 ratios' upper rank bound
  1 becomes exact (rank_second, quoted theory); the closures still rest on PARI/Sage/eclib for G, heights and
  candidates, as the report says.

## 8. Claim 5 (numbers and status words)
Checked against the outputs (all hold): PARI ranks [1,1]/[2,2]/[2,2]/[2,4] and dm_g4k a ≤ 39, 79 candidates
(h7_g4k.jsonl); ĥ(G) ∈ [9.870024, 9.884853], saturation ≤ 199, a ≤ 18, 37 candidates, t = ±1 degenerate
(h7_g4k_exact.txt, rerun identical); F4 bounds 3, 2, 3, 5 (h7_g4k_rank_exact.txt; g4k_rank_exact.txt 2, 3, 5);
Sel²/pairing table (my own recomputation, section 4); (13348,3795) 5/2, (13899,1180) bound 2 and PARI [2,2];
E04 counts 91/66/65/90, 46/16/1 (section 7); 122 in tier B, 31 D1+E0L (25+6); G3 table 916 ratios, OPEN list
with (13348,3795), (13899,1180), route list with (13103,4271) (g3_table.py rerun: `# matches`); LEDGER row
"11000 < m ≤ 13000: 4 open" includes (12201,7360).
Wrong or loose:
1. REPORT Part 4 / bottom line: "On the curve of (1015,713) all 6 curves give 3, or the pairing's local search
   gives up" — 2 distinct curves (section 7).
2. HANDOFF "Ideas tried": "It turned three conditional rank-1 statements into exact ones" — the exact rank-1
   results are four curves; two of them were F4's conditional closures.
3. HANDOFF "Ideas tried": "the seven rank-1-route exclusions" — CENSUS9000.md calls the seven "conditional on a
   rank being 1"; "the rank-1 route" is the other condition set ((a), (b)), nine ratios. Mislabel.
4. LEDGER text and HANDOFF Label: "resting on Vélu, Mazur, Silverman and Cremona–Prickett–Siksek (quoted)" omits
   the one input the new step adds: the Cassels–Tate pairing and its Brauer–Manin formula (Cassels 1998;
   Skorobogatov 6.2), quoted, as implemented in rank_second.py (REPORT's own "What it rests on" lists it).
5. "no longer needs PARI, Sage or Magma": rank_second uses PARI's qfsolve (conic points) and factor (hints) as
   finders, each checked; say "no longer rests on".
6. "the numerical L-value behind that pair": since run 6 it is an exact interval enclosure conditional on PARI's
   ellap (LEDGER row "1000 < m ≤ 11000"); say "the L-value".
7. "finished at 12:23 UTC": times.txt gives restart 10:04 + 7829 s + 443 s ≈ 12:22, and block 14000 was started at
   12:15; immaterial.
Status words: "computed (exact)" is fair for the rank step (exact arithmetic, quoted theory). "unconditional" in the
bottom line ("there it is now unconditional") is too strong: it is conditional on the quoted statements (tier E)
and on the program; "no longer conditional on PARI's rank" is the right phrase. Tier: **D** (computed by a program of
the repository, not counted; the rule it implements — the CT formula — has a prose derivation in rank_second.py's
header, re-derived in review_t2 and here, but no tier-C written proof of the whole step). Not a census change.

## 9. Kernel = image of the Mordell–Weil group (extra consistency check, claim 1 and rank0_ct)
With my pairing, the kernel on Sel²(E′) equals the span of the images of T1, T2 and the known points (points found
by PARI on E′, checked exactly): (6881,6233) 3 = 3, (13348,3795) 3 = 3, (12201,7360) 3 = 3, +14 curve 2 = 2 (torsion
only). So each bound is attained by known points: no slack, and no rational point's image pairs non-trivially.
(Images of C's generator pushed through Vélu alone gave only the torsion span for two curves — φ(G) can be ≡ T mod
2E′(Q); that is consistent, not a failure.)

## 10. rank0_ct.py (conductor's request): rank 0 for curve +14 of (609,200)/(809,409)
- Reruns (copy in /tmp/z7/repo): `python3 -I rank0_ct.py` → `# matches rank0_ct.txt` (13 s); `--negative` →
  `# matches rank0_ct_negative.txt` (3 s).
- Shift and Vélu: x0 = 15821452025664321 is a root (asserted); A = 3x0 + a2, B = 3x0² + 2a2x0 + a4 (same formula I
  derived for g4k; recomputed independently with sympy's integer root: A = 47464356076992964, B =
  3408344719936636397319208503016324 = s², s integer); E′ roots 0, A ± 2s = [−69297704883520000, 0,
  164226417037505928], equal to the first full-2-torsion member listed by indep_rank.txt's class_full_2descent
  (dim Sel² 6 there too). Vélu: polynomial identity checked in section 1. Holds.
- Point counts: #E(F_p) = p + 1 + Σ Legendre(f(x)) for y² = f(x) is right; good primes = odd p ∤ B(A² − 4B),
  which is the set of odd bad primes of both models. A necessary check only (as the header says). Holds.
- Same curve for (809,409): magma_second_ranks.json has identical a-invariants for both pairs (that is how
  indep_rank.txt prints "same curve"); independently, second_123.jsonl's rank0_curve of (609,200) (twist +1) and of
  (809,409) (twist −1) are both isomorphic over Q to rank0_ct's curve (my exact c4/c6 check: u² = 1 and u² = 4);
  second_exact.py, which rebuilds E from m, n alone and asserts the record's rank-0 curve is E or E^(−1), rerun:
  `second_exact OK: 6 closures`, output identical to second_exact.txt. SECOND_DESCENT_CHAIN.md row (809,409):
  "E⁽⁻¹⁾ (the same curve as for (609,200))". Holds.
- Independent pairing (ct_mine.py, my tangent points, my local points, two seeds): E′ of +14: dim Sel² 6, pairing
  rank 4, bound 0 (both seeds); the other two full-2-torsion members of the class: (6,4) → 0 each; kernel = image of
  torsion (2 = 2). Control (861,620) E′: (4,2) → 0; (440,279) L-value curve y² = x(x − 719²)(x + 161²): (4,0) → 2.
  All equal to rank0_ct.txt. PARI 2.17.4: s = 4 on this curve (SECOND_DESCENT_CHAIN), consistent.
- Controls in --negative: x0+1, A+2, a non-square B ((223,167)), forged e3+2 refused; without pairing bound 4;
  congruent n = 5,6,7 → 1, n = 1,2,3 → 0. Fine. Missing (nit): bound() does not check that the torsion images lie in
  the pairing's kernel (I did: they span it).
- **Wrong sentence:** header "PARI: none in this script (rank_second factors with PARI as a finder and proves each
  prime)". With cypari2 blocked the script fails in rank_second.gen_conic_point (PARI's qfsolve for the conic
  points). Replace with: "PARI: only as a finder inside rank_second.py — qfsolve for the conic points and factor
  for factorisations — each result checked exactly; no rank, Selmer or pairing output of PARI is used."
- **Verdict: holds** (computed, exact; quoted: Vélu, isogeny invariance, the CT pairing and its Brauer–Manin
  formula). This makes SECOND_DESCENT_CHAIN.md step 7's "not by any 2-descent" true only of 2-descent without the
  pairing: suggest "no 2-descent alone gives rank 0 (2-Selmer bound 4); the complete 2-descent with the
  Cassels–Tate pairing on the 2-isogenous curve with full 2-torsion does (rank0_ct.py, tier D)". Tier D.

## 11. Further checks of "What failed" (scratch claims)
- D1: dm_evidence's 25 D1 lines use 16 curves; PARI ellisomat (finder): every class has 2 curves, none with full
  rational 2-torsion. σ: the 12 curves of sigma_evidence.txt: every class has 2 curves, none with full 2-torsion.
  Both "does not apply" statements hold. Not checked: the rank_indep bounds for D1 (388,287), (920,371), the
  "[1,3] on every curve" of the seven conditional exclusions, and (6557,324)'s other curves (all scratch/not run).
- dm_g4k.py on the six G3 ratios and g4k_rank_exact.py on (12201,7360), rerun under PARI 2.17.4: identical to
  h7_g4k.jsonl + h7_b13000_g4k.jsonl and h7_g4k_rank_exact.txt.
- rank_second.py in int7 is byte-identical to run 6 (last change d2bd5db, the review_t2 FIX-FIRST in DLam only;
  DFull/full_selmer2 as reviewed).

## Final
Bottom line: H7's rank step holds. An independent complete 2-descent and an independent Cassels–Tate pairing
(my code, my tangent and local points) reproduce every Sel² dimension and pairing rank: rank 1 exactly for the
g4·K curves of (12201,7360), (7735,4968), (6881,6233), (13348,3795); the kernel equals the image of the known
points. With F4's reviewed chain (condition exactly "rank C(Q) = 1"), all solutions of those four are degenerate;
(6557,324) follows through the partner map (all nine cells permute with one common factor). E04: 90/91 ratios,
65/66 curves, per-curve agreement; (864,151) not covered. rank0_ct.py holds (6, 4, so 0; same curve for (809,409)).
Verdicts: claim 1 holds; 2 holds; 3 holds; 4 holds with a wording error; 5 holds with gaps (wording, below);
rank0_ct holds with one wrong sentence. Tier D for all of it.

Replacement text:
- REPORT_H7 bottom line, "and there it is now unconditional." → "and there it no longer depends on PARI's rank
  (computed, exact; resting on the quoted statements of 'What it rests on')."
- REPORT_H7 Part 4, "On the curve of (1015,713) all 6 curves give 3, or the pairing's local search gives up." →
  "The isogeny class of the curve of (1015,713) has two curves with full rational 2-torsion (the walk visits each
  three times, as scaled copies): E04 gives 3, and the pairing's local search gives up on the other (a separate
  computation gives dim Sel² = 7 and a pairing of rank 2 there, so 3 again)."
- Proposed LEDGER text, last sentence → "Computed (exact; tier D), resting on Vélu's formula, isogeny invariance,
  Mazur, the Cassels–Tate pairing in its Brauer–Manin form (Cassels 1998; Skorobogatov 6.2) as computed by
  rank_second.py (reviews/review_t2.md), Silverman and Cremona–Prickett–Siksek Theorems 2 and 15 (quoted, tier E).
  Written by an AI model; not read by a person." Also in it: "three ratios" → "four ratios ((12201,7360),
  (7735,4968), (6881,6233), (13348,3795))"; "no longer needs PARI, Sage or Magma" → "no longer rests on PARI, Sage or
  Magma (PARI only finds conic points and factorisations, each checked)"; "the numerical L-value" → "the L-value".
- HANDOFF "Ideas tried": "It turned three conditional rank-1 statements into exact ones" → "It gave rank 1 exactly
  for four g4·K curves, two of them F4's conditional closures"; "on the seven rank-1-route exclusions" → "on the seven
  exclusions conditional on a rank being 1". HANDOFF H7 Label: add "the Cassels–Tate pairing as computed by
  rank_second.py" to the quoted list and "tier D". HANDOFF "Beyond H7" "needs no PARI, Sage or Magma" → "rests on no
  PARI, Sage or Magma output". H8: "numerical L-value" → "L-value (conditional on PARI's ellap)".
- rank0_ct.py header "PARI: none in this script (...)" → "PARI: only as a finder inside rank_second.py (qfsolve for
  the conic points, factor for factorisations), each result checked; no PARI rank, Selmer or pairing output is used."
Nits: bound_full/rank0_ct.bound do not check point/torsion images in the kernel or tangent independence (I did);
h7_e04_rank does not assert bound ≥ 1 when G is known; least_bound dedupes only by translation.
Could not check: the scratch claims listed in §11; anything needing Magma. Finished Fri Nov 27 13:50:18 UTC 2026.

## Packaging (conductor's request, done 14:36 UTC)
/home/claude/conductor/r7/referee_z7/ (to be committed as experiments/referee_z7/): my2desc.py, ct_mine.py,
velu_check.py, partner_cells.py, run_ct_g4k.py, ct_rank0.py, classes.py, run_ct_e04.py, forge_e04.py, each with its
.txt, z7_points.json (points found by PARI ellrank in --find modes, checked on use), README.md. Each output was
regenerated twice and is identical (deterministic); run_ct_g4k and ct_rank0 also checked from r7 with MSS_REPO/H7_DIR/
H1_DIR pointing at int7 (python3 -I -B; int7 untouched). Command: python3 -I X.py 2>&1 | grep -v '^PARI stack' | diff - X.txt.
