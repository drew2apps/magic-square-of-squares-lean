# Review z5 of run6_y5 (Fisher Theta_4 generalised) — hostile referee

Written by an AI model; not read by a person. Started 00:09 UTC 27 Nov 2026. Scratch: /tmp/z5/.

## 00:20 UTC - orientation (before reading reports/y5.md)
- Commits under review: 515c07b..493f0ba (5 commits; the doc changes in 419c91a..run6_y5 come from run-5 commits dabf7b4, acb95e1, not y5).
- All four computations run on the SAME orientation as (72,65): the lifted class is on the E' side (S_1 = S^(phi)(E/Q), classes d | a^2-4b), C3' is a torsor under E, C4 a torsor under E', g0 a pushout form for T', the paired y on the E side (S'_1, classes d | b), K2 classes b and x(G). For 920,129 and 760,39 this is rank_fisher's first-named side; for 631,71 it is the second side (--side=E'). The swapped-role branch of prepare() (side == 'E', working curve E') is exercised by NO output and NO CI line.
- Clean clone of run6_y5 in /tmp/z5/clone: rank_fisher_theta4_gen.py on 72,65 / 920,129 / 760,39 / 631,71 --side=E' (main and --negative), _order (both), _recheck reproduce byte for byte once the PARI stack notice (printed on STDOUT, not stderr) is filtered. Times: main 13-18 s each, negative 4-12 s, recheck 58 s, order < 3 s.
- _control (23 s) and _weil (50 s) also reproduce byte for byte. All 13 runs: exit 0, nothing on stderr.

## 04:20 UTC - resumed after the session restart (hard stop now 08:00 UTC)

### Item 2: independent local terms (referee script /tmp/z5/indep_local.py, outside the repo)
Uses the repo only to BUILD the object under test (QA, QB, f0, Xi, R3^0, N1^0). Independent: (e:g)-first sampler and a
(c:d)-first sampler written here; where both fail after 3 seeds, the repo's D_lam-guided (c:d) values are used ONLY as disk
centres, moved by a random multiple of p^j (j = 1..6), with (e:g), scale and W computed here (marked "seeded"); my own integer
Hensel certificate (v(F_i) >= 2 v(3x3 minor) + 1, class room +1 / +3 at 2); my own Hilbert-symbol formula, asserted equal to
PARI hilbert() on every term; my own place set S (own proved factorisations of b, a^2-4b, the pencil quartic's discriminant
AND its two end coefficients, primes of Xi and y, and a DIRECT degeneracy test: primes of the gcd of the 4x4 minors of
[QA, QB, Xi W^2 - f0, R3^0] and of the 3x3 minors of [QA, QB, Xi W^2 - f0] over the 15 quadratic monomials). As a sanity
check every prime 47..260 outside S was also evaluated: all terms +1.
| run | Xi (index) | lift | y-classes : product | -1 places for the open y |
|---|---|---|---|---|
| 760,39 | 219 (0) | repo k | 761: -1; b, x(G): +1 | 7, 761, 4903 (= repo) |
| 760,39 | 24758726262 (3, not used by repo) | k1*3, k3*5 | 761: -1; b, x(G): +1 | 7, 13, 4903 |
| 920,129 | -67071013 (0) | repo k | 156252226: -1; b, x(G): +1; y x(G) = 1758233776793 (= 1 mod 8): -1; y b x(G) = 253280314753626857 (= 1 mod 8): -1 | y: 2, 929, 6469 (= repo); y x(G): 3, 5, 929, 53299, 145585309 |
| 920,129 | -43154224789 (5, not used) | k1*3, k3*5 | 156252226: -1; b, x(G): +1; 1758233776793: -1 | y: 2, 43, 6469 |
| 631,71 E' | -15 (0) | repo k | 31836878641: -1; b: +1 | 2017, 8053, 1178063029213507 (= repo) |
| 631,71 E' | -23614549987199546592015 (3, not used) | k1 -> conj(k1) | 31836878641: -1; b: +1 | 2017, 15784273, 1178063029213507 |
| 72,65 (control) | 137 (0) | repo k | 2: -1; b, x(G): +1 | 2 (= repo) |
"Seeded" places: 760 at 761, 1949, 189661 (and 42227 for index 3); 920 at 159793; 631 at 2, 787, 1162771. Everywhere else the
points are the referee's own. Every term was constant over all certified points (8 to 368 per place).
Consequences: (i) every value -1 is reproduced with other local points, other level-4 lifts and other k_1, k_3; (ii) for 920,129
the -1 does NOT rest on the 2-adic symbol (two representatives = 1 mod 8 of the same class modulo the kernel give -1 with no
2-adic term); for 760,39 (y = 761) and 631,71 (y = 31836878641) y = 1 mod 8, so the 2-adic term is +1 identically; all y > 0,
so the real term is +1. (iii) My S contains the repo's S for every run; the extra primes (from the end coefficients of the
pencil quartic, which are NOT needed: the binary-form discriminant is the right invariant) all give +1.
(iv) The direct degeneracy test agrees with the repo's den(r) test: the primes where R3^0 mod p lies in the span of the
quadrics of C4 are exactly the repo's den(r)-degenerate primes plus primes already in B0 (920: 179, 53299, 6375647,
743810044013624383; 760: 4903, 11298841251693869044727293547; 631: 8053, 1178063029213507). So no degenerate prime is missed.
VERDICT item 2 (local terms and sums): HOLDS for all three curves and the control.

### Item 2 (cont.) and item 3: own exact checks of the inputs of the evenness argument (/tmp/z5/check_exact.py)
Repo used only to build l_1, l_3, q, c, QA, QB. Written here: (a) f = N(l_1), f_dag = N(l_3), and f f_dag - c q^2 written
exactly as sum m_i QA + n_i QB in degree 4 (PARI matinverseimage over Q, then re-expanded to 0); (b) tangency: the roots of
det(w A + B) in L, L' recomputed; the cone there has rank 3 and its restriction to {l_1 = 0} (resp. {l_3 = 0}) has rank 1
(a double line); (c) (I2) resultants of QA, QB on {alpha_1 = beta_1 = 0}, {alpha_3 = beta_3 = 0} and, over L, on
{l_1 = delta3 = 0} for two twists, all non-zero; (d) (I3') the ruling {l_1 = m = 0} lies on the cone and QA on it has non-zero
discriminant. Result: ALL PASS for 760,39, 920,129, 631,71 (E') and 72,65. (I1) itself is then automatic: N1 N1' - D1 delta1^2
= 4 f / N(k1) - 4 xi4 W^2 / N(k1) and R3 R3' - D3 W^2 delta3^2 = 4 (W^2 f_dag / N(k3) - tau3^2 q^2), both in I(C4) given (a);
the repo also checks them by linear algebra.
Item 3 (the written evenness argument). Re-derived line by line: it uses exactly (I1), (I2) (three resultants), tangency of
l_1 (div l_1 = 2 A_1, so div f = 2(A_1 + A_2), C4 -> C3' etale, ord_P W = m), q >= A_1 + A_2 (which follows from (a) and
tangency: div q = A_1 + A_2 + div(f_dag)/2), (I3') (m = 1), and b not a square (so E'(Q)[2] = {0, T'} and (I4) forces
[B] = T'). Nothing in it depends on L = Q(sqrt -2). Every curve-specific hypothesis is checked per curve and per twist by the
repo (check_g0 for (I1), (I2), (I3'); rank_fisher5.analyse for tangency) and independently by me as above; b is not a square
(sqfree b = 515996348257602069980161, 111314284150462675500481, 1574303332479969772801). (I4): my own certified points give
several classes of g0 at p = 2 on every curve (e.g. 760: (0,3),(0,5),(1,3),(1,5)). The argument carries over verbatim; no
missing check found. VERDICT item 3: HOLDS (a written argument, plus the mod-p Weil test; not a proof over Q by machine).

### Item 4: bad places
The good-reduction argument (p odd, p outside B0, p not dividing Xi) re-derived: C3~ smooth (pencil discriminant a unit; the
binary-form discriminant is the right invariant, the end coefficients are not needed); each rank-3 cone over L reduces to a
rank-3 cone, and the reduction of a tangent plane is tangent (the cone restricted to the plane is kappa m^2 with kappa a unit,
since a rank-3 cone contains no plane); so f0~ = unit x l_1' conj(l_1') has divisor 2(A_1~ + A_2~) and f0~ can NOT lie in
span(QA~, QB~) at such p. Hence the new 'f0 degenerates' check is a consistency check: its primes must lie in B0, and they do
(920: 2,3,23,43,113,929,159793,4883483573; 760: 2,3,5,19,103,761,1949,23761,42227,189661; 631: 2,3,5,7,13,449,787,1162771 -
all in B0; my own 3x3-minor test gives the same lists). N1^0, R3^0 not vanishing on C4~: the repo's den(r) argument is
confirmed by my direct 4x4-minor test (item 2 (iv)). primes(y), primes(Xi), all p < 47, B0 and the degenerate primes are in S.
B0's bad(E) part equals my own factoring of 2 b (a^2-4b) on every curve. Every member of S is a proved prime in my run
(pari isprime); the repo's primes_of() does not prove primality (BPSW factors), but pencil_primes does for the cofactor,
and my run covers the rest (NIT). Nothing is assumed about primes outside S beyond the written argument; 37-39 extra primes
in 47..260 give +1 in my runs. VERDICT item 4: HOLDS.

### Item 1: faithfulness and the side of (631,71)
- Code diff against rank_fisher_theta4.py: setup_gen = setup minus the assertion (D1, D3) = (-2, 2) (line for line);
  theta4_value = theta4_for plus the f0-degeneracy primes in S and generic labels; build_g0, check_g0 (I1)/(I2)/(I3'),
  normalise, degenerate_primes, local_terms, real_term, hilbert are imported unchanged and are generic in D1, D3.
  choose_Xi_lin replaces the enumeration by F_2 linear algebra; this only CHOOSES Xi: every Xi used is then shown ELS by
  local_terms (points at every p of S, asserted) and real_term (asserted), and the rest by the good-reduction argument. So
  a wrong choose_Xi_lin could only lose twists, never certify a wrong one. (It does lose some: 760,39 gives 8 twists while
  rank_fisher's bound for S'_1 there is 5-dimensional; a sampler limit, harmless.)
- Orientation. For all four runs, the lifted class y0 is in S_4 inside S_1 = S^(phi)(E/Q) (E' side), [C4] in
  S^(4)(E'/Q), g0 a pushout form for T' (generator of ker phihat), y in S'_4 inside S'_1 = S^(phihat)(E'/Q) (E side; the
  working curve is E, K2 classes x(T) = b and x(G)). This is literally the (72,65) set-up and Fisher's Sect. 3 formula
  Theta_4(y0, y) = <x_1, y> with x_1 any lift of y0 in S^(4)(E'), which needs y in S'_4 (Thm 3.1(i)). phi and phihat are
  NOT swapped anywhere. For (631,71), --side=E' is the SAME orientation as (72,65); what differs from rank_fisher's named
  side is only that on E' S_1 (dim 5) is bigger than S_2 = S_3 = S_4 (dim 3), so a given level-3 lift need not have a
  non-empty level-4 Selmer set; the code requires one (4 ELS twists found and checked), which is all that is needed.
- Rank bound for (631,71) on E': S_4 (E' side) = <2 known classes, y0 = -627214> (dim 3), S'_4 (E side) = <b, y> (dim 2;
  x(G) is a square: PARI core(x(G)) = 1, so the image of E(Q) is <b> only). The known classes lie in both kernels
  (Thm 2.1); Theta_4(y0, y) != 0 cuts y0 out of the left kernel and y out of the right one, so S_5 = known (dim 2),
  S'_5 = <b> (dim 1), and rank <= dim S_5 + dim S'_5 - dim E(Q)[phi] - dim E'(Q)[phihat] = 2 + 1 - 1 - 1 = 1. The printed
  "dim 2 and 1 ... 2 + 1 - 2 = 1" is right. For 920,129 and 760,39: (2-1) + (3-1) - 2 = 1, right. DIMS (2,3), (2,3), (3,2)
  agree with rank_fisher.txt lines "(5)". The open classes are those rank_fisher.txt names (y0 = -23, -799, -627214;
  y = 156252226, 761, 31836878641), and y's primes are bad primes of E (y unramified elsewhere), as S'_1 requires.
- Own PARI checks of the curves: G on E, ellorder(G) = 0 (infinite order; reductions mod 11, 17, 23 give the orders the
  _order script prints), torsion Z/2 on E and E' (so b and a^2-4b are not squares; a^2-4b = -1 mod squares on all three),
  ellrank(E) = [1, 3] on all three (the 2-descent bound 3 that Theta_4 lowers).
- Inherited, not re-audited (as in review_w4 G2): rank_indep's first descent, rank_second, rank_fisher's level-2/3 ELS
  lifts and dims (so y0 in S_4, y in S'_4), rank_fisher5's construction of C3'.
VERDICT item 1: HOLDS (for all three curves; tier C because Fisher is read through summaries).

### Item 5: negative tests and controls
- --negative (9 refusals per curve) runs through the GENERAL code (prepare/setup_gen/choose_Xi_lin and the shared check
  functions), not copies. Reasons checked in the outputs: conic point off its cone -> rank_fisher5 "P1 is not on the cone";
  q -> ideal test; sigma_1, tau_3 -> (I1); k_1 = 1 + sqrt D1 -> sigma_1 not rational; g0 c/e -> the term depends on the point
  (at a place where y is not a square, chosen so; good); precision 4 -> no usable 2-adic point. All for the intended reason.
  Weak ones: the "forged twist" refusal is "class not in the closure of the SAMPLED classes at 2", a statement about the
  sampler, not a certified non-existence (the certified pipeline is not run on the forged Xi); the "wrong y" test exercises
  rank_second.certified (inherited), not the new code. Missing (NIT): a mutation that drops a degenerate prime (e.g. 8053
  and 1178063029213507 on 631, which carry two of the three -1) and shows the product changes.
- Controls: _control re-runs the four BSD controls through setup_gen/choose_Xi_lin/theta4_value (monkeypatched into
  rank_fisher_theta4_control; the transcript has the new "f0 mod p" lines, so the general path is used). A, B give -1
  (non-zero predicted), C, D +1 (zero predicted): yes, two non-zero values are tested. Gaps: BSD-conditional (as before);
  all four, and the (72,65) positive control, have (D1, D3) = (-2, 2) and y1 = 1 mod 8, so NO control with a known answer
  exercises a general quadratic field or the 2-adic symbol. The new curves' own K2 checks (2 classes, 4-8 non-trivial
  local terms each, product +1, a theorem) and K3/K5 are the field-generic checks; my independent run adds other twists,
  other k_1/k_3 and (920,129) a representative with no 2-adic term.
VERDICT item 5: HOLDS WITH GAPS (no known-answer control in a general field; forged-twist test is about the sampler).

### Item 3 (cont.): Weil test at two other primes
The repo's test function at p = 30011 and 1000003 (seed 3, 15 pairs each) on 920, 760, 631: every ratio a square; the
odd-divisor control g0 c/e fails 4-8 of 15 each time. Code read: sections W = lin(c, d, e, g) are the right choice (for
sections not involving W, N(Z) = prod (-D3/D1)(delta3/delta1)^2 is a square by (I1) whatever the divisor, so they would
test nothing); chart g = 1 with shape-position check; norms by determinants of multiplication maps. Sound as a test.

### Item 4 (cont.): the real place
Own interval-arithmetic check (/tmp/z5/real_check.py, mpmath iv at 120 digits, (e:g) chosen first, exact formulas for rho,
mu, W; enclosures of QA, QB contain 0): a real point of C4 exists for the first three twists of each curve (920: Xi < 0,
760: Xi > 0, 631: Xi < 0). All paired y are positive, so the real term is +1 whatever the sign of g0.

### Item 6: the proposed CI lines, literally (reports/y5.md lines 141-151 plus the Weil line of line 185)
Run as one bash script with `set -eo pipefail` from the clean clone /tmp/z5/clone (run6_y5 = 493f0ba), sympy 1.14, PARI
2.17.4: exit 0, every diff empty, 201 s in all (the proposal says about 4.5 minutes). The PARI stack notice is on stdout,
so the `grep -v 'PARI stack'` in each line is needed (the order lines do not need it and do not have it). VERDICT: HOLDS.

### Item 7: wording
The script docstrings, outputs and the proposed LEDGER row stay at tier C (they name Fisher read through summaries,
written arguments, BSD-conditional controls, "not reviewed"). Points to fix:
- CHANGES draft: "An exact check was added that the (72,65) script lacked (the primes where the pushout form degenerates
  modulo p)" - it is f0 (the form of the double cover), not the pushout form g0, and by the good-reduction argument it can
  only fire at primes of B0; it is a consistency check, not a gap of the (72,65) script. Suggested: "A consistency check
  was added: the primes where f0 mod p lies in span(QA, QB) mod p (by the good-reduction argument they lie in B0, and do)."
- The docstring's (b) "for a reduction on the E side the roles of E and E' are exchanged" describes code that no output
  or CI line runs (no hints exist for it); say "unused, untested", or remove it. As written, its K2 list would hold only
  the class of T' (sqfree(a^2-4b)), incomplete for a curve whose E'-side image has dimension 2 (e.g. (631,71)).
- Doc text the proposal does not spell out: README line 119 also says G generates "on Magma alone for the curves of
  (760,39), (945,248) and (920,129), and on Magma and Fisher's Theta_4 for the curve of (72,65)"; and LEDGER line 53,
  third column, "(for the curves of (760,39), (945,248) and (920,129), so for four counted ratios with the partners, from
  Magma alone; ...)". Both must change with the count: (945,248) alone rests on Magma alone; (72,65), (760,39), (920,129)
  on Magma and Theta_4. The count itself is right: the six of README line 119 are (760,39), (799,721), (945,248),
  (920,129), (631,71), (351,280); minus five leaves (945,248).
- After this review the row's "Not reviewed" should become "one AI referee (reviews/review_y5.md): holds with gaps".
VERDICT item 7: HOLDS WITH GAPS (text fixes only).

## Verdict (FINAL, 04:30 UTC)
**HOLDS WITH GAPS, per curve: (920,129) E+ HOLDS WITH GAPS; (760,39) E+ HOLDS WITH GAPS; (631,71) D1 (E' side) HOLDS WITH
GAPS; (72,65) control reproduced; (945,248) correctly reported as not reached.** I found no error that changes any value.
Theta_4(y0, y) = -1 on all three new curves, reproduced by my own points, Hensel certificates, Hilbert symbols and place
set, for other twists, other k_1/k_3 and (920,129) a representative with no 2-adic term; so rank <= 1 + 2 - 2 = 1
(631: 2 + 1 - 2 = 1, arithmetic and sides right) and rank = 1 with G of infinite order, without Magma, at tier C.
Gaps (none an error): G1 Fisher read through summaries (Thm 2.1, Lemma 2.2, Sect. 3, Thm 3.1(i)); G2 evenness of div g0
and the good-reduction step are written arguments (re-derived, carry over verbatim; every curve-specific hypothesis
checked by the repo and by me; Weil test at four primes); G3 inherited rank_indep/rank_second/rank_fisher/rank_fisher5
not re-audited; G4 the known-answer controls are BSD-conditional and all in Q(sqrt -2), Q(sqrt 2) with y1 = 1 mod 8.
Items: 1 HOLDS; 2 HOLDS; 3 HOLDS (as a written argument); 4 HOLDS; 5 HOLDS WITH GAPS; 6 HOLDS; 7 HOLDS WITH GAPS.
FIX-FIRST: none. FIX: (1) the CHANGES sentence on the "exact check the (72,65) script lacked"; (2) mark/remove the unused
E-side branch; (3) README line 119 and LEDGER line 53 tails about Magma alone; (4) review status in the row.
NIT: (1) forged-twist negative test is about the sampler, run the certified pipeline on the forged Xi; (2) add a mutation
dropping a degenerate prime (on 631, 8053 and 1178063029213507 carry two of the three -1); (3) primes_of() does not prove
primality (every prime of S is proved prime in my run, so cosmetic); (4) a known-answer control in a general field.
Not checked: Fisher's printed text; the finders (inputs are re-checked exactly); the other-side computation for (631,71);
rank_fisher's level-3 lifts. Referee scripts: /tmp/z5/indep_local.py, check_exact.py, real_check.py, curves_pari.py,
weil_more.py, ci_literal.sh (outside the repo; no repository file edited).
