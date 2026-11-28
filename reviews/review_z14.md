# Review z14: Cassels-Tate pairing on the mixed pair curves (EXPLORATION_3 §1, §2, §4; EXPLORATION_4 §5)

Written by an AI model (reviewer z14), not read by a person. Scratch and scripts: experiments/referee_z14/ in

> Committed as written apart from this note. In this repository the reviewer's folder is `experiments/referee_z14/; the notes reviewed are `experiments/research-pairing/explore3/`, `explore4/` (section 5) and `explore5/` (section 1)`. Paths under `/home/claude/` are the reviewer's scratch folders.
worktree /home/claude/conductor/z14 (branch run11_z14). Findings are appended as they are made.

## Log
- Copied explore3's cassels.py, sel.py to experiments/referee_z14/e3_cassels.py, e3_sel.py; the only edits are the
  sys.path lines (/tmp/t10 removed, /tmp/mix replaced by the script's folder). Read before running: no file writes,
  no network; PARI qfsolve, factor, hilbert, p-adic sqrt; random local-point search.
- EXPLORATION_4_2026-11-27.md is NOT in /home/claude/untrusted/explore4/explore4/ (nor anywhere on this machine:
  `find / -name 'EXPLORATION_4*'` finds nothing). For its §5 I review only the scripts and outputs that are there
  (kappa.py, kappa.txt, tangent_norm.py, tangent_norm.txt) and the claims as quoted in my brief.
- Early structural observations (from the pilot files and t10_pairct.Pair on new ratios):
  * The four "mixed curves" E_0L, E_0U, E_8L, E_8U of a ratio are ONE curve: the pilot prints the same integral model
    e = (e1, e2, e3) for all four (e.g. (465,416): e = (-930968778390566402, -930968508713617920, -930968441908727040)
    on 0L, 0U, 8L, 8U), and Pair(m, n, j, k) on new ratios (1234,567), (7,4) gives e-differences (2R-, 2R+) with u = 1
    on all four. t10's NOTES.md already says C_0L = C_0U = C_8L = C_8U; the four test different classes (d_j d_k).
    So "on all four mixed curves" is one curve per ratio, and every per-curve count in §1 is 4 x a per-ratio count.
  * (m, n) and (m+n, m-n) (made primitive) give isomorphic mixed curves: (A', B', C') = 2(B, A, C), so R'+- = 4R+-,
    the e-differences scale by 4 = 2^2. The pilot contains such partners (e.g. (465,416)/(881,49), (1852,1155)/(3007,697),
    (3320,2303)/(5623,1017) all appear in task23's 60-ratio sample, with identical bases of {d}); so counts "per ratio"
    in §1 also double-count curves.
- Task 5 (claims; `claims_z14.py` -> `claims_z14.txt`): the pilot (parse.load on t10_15000_s100_[ab], t10_997_[ab]) has
  201 distinct ratios but only 160 distinct mixed curves (partners). Mixed kill lines 6,776 = 4 x 1,694 (one curve
  printed four times, each line a different tuple) -> 1,183 distinct (ratio, a, b), 1,009 distinct (curve up to
  isomorphism, a, b). E_08 1,280 lines -> 133 distinct (ratio, a, b); E_LU 874 -> 111. "2,154 lines" = 244 cases.
  First component of b: 2,760 / 488 / 3,528 lines = 400 / 79 / 530 distinct (curve, a, b). Every mixed kill line has
  a with first component 2, so no killed class of the pilot is of the form (1, d, d).
  task23 samples: seed 11 = 60 ratios = 57 curves; seed 12 = 201 ratios = 160 curves.
- Task 3 (symbolic.py -> symbolic.txt, sympy): Delta(E_08) = Delta(E_LU) = 256 (mn(m^2-n^2)(m^2+n^2))^4: six I4 at the
  zeros of mn(m^4 - n^4) incl. infinity (c4 != 0 there), sum 24, (deg a2, a4, a6) = (4, 8, 0): K3. Halves of all
  three 2-torsion points of E_08 lie over Q(i)(r) (x = +-AC, A^2 +- iAB, C^2 +- BC; y^2 checked squares), so E[4] is
  rational over Q(i)(r); j(E_08) = 256 W^3/(ABC)^4 has degree 24 in r: the r-line is X(4) (over Q(i)). The twist of
  E_08 by -1 is y^2 = x(x - B^2)(x - C^2) (E_LU up to 2^2). Delta(E_mix) = 2^14 (mn(m^2-n^2))^2 (R+ R-)^2: twelve I2 at
  r = 0, infinity, +-1 and the 8 zeros of W, smooth at r = +-i; K3. j(E_08) = 0 exactly at the zeros of W
  (W = A^4 - A^2C^2 + C^4 = C^4 Phi_12(A/C)); j(E_mix) = infinity there. The corner forms vanish on t = +-ir, +-i/r and
  the conjugates (checked). psi(iota t0) on E_mix has infinite order at 5 specializations (section_z14.py; PARI
  ellorder as a finder, 4P != O checked), so E_mix has a section of infinite order over Q(r).
  Two corrections: (i) "the eight zeros of W, the points with j = 0" must say j(E_08) = 0 (the modular j); on E_mix
  these are I2 fibres, j = infinity. (ii) "its automorphisms are the octahedral group": as a P^1 its automorphism group
  is PGL2; the octahedral group S4 = PSL2(Z/4) is the group permuting the six cusps (and only a dihedral group of
  order 8 of it is defined over Q: r -> +-r^{+-1}, +-(r+1)/(r-1)^{+-1}).
  Also: R+ and R- are irreducible over Q with root field Q(zeta_12), so every odd prime of W is = 1 mod 12, and
  W = R+^2 = 1 mod 8 (R+ - R- = 4mn(m^2-n^2)). Hence the "controls" S*J(-1), S*J(2) of task23.py (14,160/14,160,
  "J factor = -1 on 0") are identically true and test nothing.
- Scope extension received: EXPLORATION_5 §1 "Theorem B" (/home/claude/untrusted/explore5/explore5/). Read
  EXPLORATION_5 §1, redei_k.py, redei_k.txt, redei_dbg2.py (tangent_norm.py is byte-identical to explore4's).
  redei_k.py compares its right side with pipeline.cassels (explore4's copy of cassels.py, the same engine), not with
  rank_second; it needs /tmp/rn, /tmp/mix (not created here).
- EXPLORATION_5 §1, Theorem B, line by line (my own check; the algebra by hand, Step 1 also in symbolic.py):
  * Set-up and Cassels' formula: <a,b> = prod_v (L2 L3 (P_v), beta)_v for b = (1, beta, beta) is the tangent-line
    formula of cassels.py / rank_second.DFull (H_i omits u_i and goes with b_i; b1 = 1 drops out). Correct. It is the
    one external input (Cassels 1998); the per-place factors are constant in P_v (asserted at every place by
    rank_second).
  * Step 1: L_k(P) = -N(z_k - q_k)/2 is the tangent lemma with qT = T = 1 (cassels.py's primitive integral L_k is a
    fixed rational multiple: dropped by the product formula over Q). N(z_k) = c_k, so eps = z/q has norm 1;
    Hilbert 90 holds in the etale algebra K (x) Q_v in all three cases (split: eps = (t, 1/t) = conj(eta)/eta with
    eta = (1, t)); eps != 1 forces y != 0 and z - q invertible (split case: (t - 1, 1/t - 1) is a unit unless t = 1);
    N(eps - 1) = N(-2y sqrt d)/N(eta) = -4d y^2/N(eta). So L_k = 2 c_k d y^2/N(eta) = 2 c_k d N(eta) mod squares.
    Correct.
  * Step 2: (N eta, beta)_v = prod_{w|v} (eta, beta)_w (norm compatibility, beta in Q_v; true for split, inert,
    ramified v, v = 2, and v = infinity: d < 0 gives one complex place, both sides 1); bimultiplicativity; conjugation
    permutes the w | v and Hilbert symbols are Galois-equivariant, so prod (eta, pi-bar)_w = prod (eta-bar, pi)_w;
    eta eta-bar = eps eta^2. Correct (the chain is also just (N eta, N pi)_v = prod_w (N eta, pi)_w).
  * Step 3: prod_v (4 c2 c3 d^2, beta)_v = 1 (Hilbert reciprocity over Q, c2 c3 rational) and prod_w (q2 q3, pi)_w = 1
    (Hilbert reciprocity over K, q2 q3 and pi global). Correct. The right side is a finite product (at w outside the
    primes of 2 d beta z0 c2 c3 the factor is 1 at a point with unit u_i, and the per-place factor does not depend on
    the point because the Cassels factor does not).
  * Verdict on the derivation: correct, given Cassels' formula; I found no gap. It is a reformulation of Cassels'
    formula (no reciprocity law beyond Hilbert reciprocity over Q and K), valid on any curve with full 2-torsion.
  * "Reading as a Redei symbol" overstates. (gamma_w, pi)_w = (pi/w)^{w(gamma_w)} needs pi a unit at w, but
    pi = (x0 + y0 sqrt d)/z0 is NOT a unit above the primes l of z0 (valuations +-v_l(z0) at the two primes above l),
    and there (gamma unit) the factor is (N gamma / l)^{v_l(z0)} = (4W/l)^{v_l(z0)}; these cancel against the z0 in
    pi at the primes of W only through Jacobi reciprocity (W0 = 1 mod 4). Also the split of v_p(gamma) between the two
    primes above p | W depends on the local point; it is harmless only because beta is a square mod p (V in N). So the
    places not accounted for are those above 2, infinity, d, beta, z0 and the bad primes of AB, not only "2, d, beta,
    infinity". Theorem B is an intermediate identity; the Redei form prod_{p | W odd} ((x0 + y0 s)/p) is NOT derived
    until those factors are computed (the note's own "Not done" list says the factors above 2, d, beta, infinity are
    not written; the sentence "So the Redei form is now derived" contradicts it). [fix]
  * Coverage: in the pilot every mixed kill line has first component of a equal to 2 (6,776 of 6,776; 1,009 distinct
    cases), so Theorem B (both first components 1) covers none of the pilot's kills. The rational point psi(iota t0)
    has class (1, -1) at (e1, e2) on all four mixed models (7 ratios checked), so it cannot move first component 2
    to 1; the 2-torsion images have first components W, 2R-, 2R+ (mod squares).
- Stevenhagen, arXiv:1806.06250 (Math. Proc. Camb. Phil. Soc. 172 (2022) 627-654), checked on the paper (WebFetch of
  the arXiv PDF): Theorem 1.1 (= Thm 7.9): for a, b, c in Q*/Q*^2 with (1) (a,b)_p = (a,c)_p = (b,c)_p = 1 at every
  p <= infinity and (2) gcd(Delta(a), Delta(b), Delta(c)) = 1, the Redei symbol [a, b, c] of Def. 7.8 (Artin symbol
  of an ideal of norm |c| in a MINIMALLY RAMIFIED F/K, K = Q(sqrt ab), Defs 7.5-7.6) is well defined, trilinear and
  [a,b,c] = [b,a,c] = [a,c,b]. The local part at p | c is "essentially a Legendre symbol (beta/p)" in Q(sqrt a) at a
  degree-1 prime unramified in Q(sqrt a, sqrt beta); the paper does not state a product formula over p | c for an
  arbitrary solution, and nothing in it concerns Cassels-Tate pairings. Consequences for EXPLORATION_3 §1:
  * Condition (1) is symmetric: when x^2 - d y^2 = beta z^2 has no point, (d, beta)_v = -1 somewhere and NEITHER
    [d, beta, W] NOR [W, d, beta] is defined. "For 2,248 of the first 2,995 pairs ... only the second reading is
    defined" is wrong: on those pairs the pairing is not a Redei symbol in Stevenhagen's sense at all. [fix]
  * For the soluble pairs, conditions (1)-(2) hold here (d, beta in N(Q(sqrt W)) on every curve tested; Delta(W) = W0
    is odd, = 1 mod 4 and prime to d beta since gcd(AB, R+R-) = gcd(AB, C^2) = 1).
  * Well-definedness of the tested product (my argument, elementary, not checked by anyone else): for a primitive
    solution, gamma = x0 + y0 sqrt d has even valuation at every prime of K not above 2 d beta (an inert l | z0 would
    divide x0, y0; a split one divides only one of gamma, gamma-bar, to the order 2 v_l(z0)). So two primitive
    solutions differ by a rational t in <-1, 2, primes of d beta> modulo K*^2, and (t/W0) = 1 for all such t: (-1/W0) =
    (2/W0) = 1 (W0 = 1 mod 8, primes = 1 mod 12) and (q/W0) = (W0/q) = (W, d)_q = 1 for q | d (same for beta). Also
    (gamma/P)(gamma/P-bar) = (beta z0^2/p) = 1, so the sign of s does not matter, and the "zero" rule of task23.py
    (switch to x0 - y0 s when p | x0 + y0 s) is the value at the conjugate prime, the correct convention. Hence
    prod_{p | W odd} ((x0 + y0 s)/p) is a well-defined function of (d, beta); the same twist argument shows that the
    minimal-ramification normalisation of Stevenhagen's Def. 7.8 does not change it, so it equals [d, beta, W]. The
    "independence of the solution (14 solutions)" observed in explore3 is therefore forced, not evidence.
  * What remains for the identity <(1,d,d),(1,beta,beta)> = [d, beta, W] after EXPLORATION_4 §5's lemma: the lemma
    and its consequence <a,b> = prod_v (N(w2) N(w3), beta)_v are correct (constant -2qT: the two signs cancel in
    L2 L3, and prod_v (4 qT2 qT3, beta)_v = 1 over ALL places incl. 2 and infinity; cassels.py's finite place set
    contains 2, infinity and every prime of qT not dividing d, so the truncated product is also unchanged). But
    w2, w3 are LOCAL elements (values at P_v), so no product formula collapses the right side; one needs every local
    factor: at the places above 2, infinity, d, beta, AB and the primes of the conic points. EXPLORATION_5's Theorem B
    removes the conic points (correct), but the factors above 2, infinity, d, beta, z0 and AB are still not computed;
    reciprocity over Q (Jacobi, for the z0-terms) is needed in that last step. Not derived.
- Task 1 interim (mixtest.py, 127 of 134 curves; final numbers in mix_analyse.txt): the 134 test ratios (ratios.py,
  seed 14; m <= 15,000; 2-adic condition 4 | mn or m = +-n mod 8; none in the pilot and none a pilot partner; 80 nd
  opposite parity, 40 nd odd-odd, 8 square type, 6 twice-square type) are 134 distinct mixed curves. On every one the
  model y^2 = (x + 2C^2)(x + 2AB)(x - 2AB) is t10's Pair model for 0L, 0U, 8L, 8U up to translation (u = 1).
  (a) cassels.py = DFull (bilinear extension of DFull's basis matrix) on 15,441 of 15,441 pairs (a, b), a over the
  basis and random elements of Sel^2 (all of Sel^2 on 34 curves), b over the basis; no failure; cassels.py with a
  second seed: same value on 1,524 of 1,524; DFull direct on non-basis a = bilinear extension on 1,400 of 1,400.
  Note: cassels.py implements the same formula as DFull (prod (L_i, b_i) = prod (b1, L1/L3)(b2, L2/L3)), so this
  checks the implementation (conic points, local points, real place, place set), not the theory.
  (b) Redei form = DFull on 838 of 838 soluble pairs (169 with pairing -1), 6,704 primitive solutions, both signs of
  s; 388 solutions hit p | x0 + y0 s and the conjugate-prime rule gave the same value; never solution- or
  sign-dependent (forced, see above). 4,185 sampled pairs have no point on x^2 - d y^2 = beta z^2 (1,266 of them with
  pairing -1): there the Redei form says nothing.
  (c) V = {d : (1,d,d) in Sel^2} is inside N = {d : (d, W)_v = 1 for all v} on 127 of 127 curves, supported on -1 and
  the primes of 2AB, but V = N on 0 of 127 (dim V 1..8, dim N 3..10): "d a norm from Q(sqrt W)" is necessary, never
  sufficient; the rest are the local conditions at 2 and the primes of AB. The conic statement itself is right:
  qfsolve finds a point on R+ s^2 - R- D^2 = 2AB d w^2 exactly when d is in N (598 + 2,450 tests), and it is a theorem
  (2AB R+ = R+^2 - W is a norm).
- Task 1 final (mix_analyse.txt; 134 of 134 curves): (a) cassels.py = DFull on 16,401 of 16,401 (a, b); second seed
  1,608 of 1,608; DFull direct vs bilinear 1,480 of 1,480; no failure anywhere; DFull's basis matrix symmetric (402
  checks) with diagonal 1. (b) Redei form = DFull on 863 of 863 soluble pairs (823 with d != beta, by replaying the
  draws: redei_trivial.py; 169 with value -1), 6,904 primitive solutions, 417 with p | x0 + y0 s; DFull direct on 320
  Redei pairs = bilinear. 4,279 sampled pairs insoluble (1,266 with pairing -1). (c) V in N on 134 of 134, V = N on 0.
  No disagreement anywhere, so there is no smallest disagreement to report.
- EXPLORATION_4 §5 constants (kappa_z14.py): kappa_r = 2 and kappa* = sqf(mn(m^2 - n^2)) x (1 on 0L, 8U; 2 on 0U, 8L),
  e_r = e1, e* = e2, with sign, on 80 of 80 (ratio, curve) blocks of 20 non-pilot ratios (12 nd of both parities, 8
  square type), the constants read by t10's PairData from psi(iota t0), the 2-torsion and their sums (no ellrank). On
  the pilot (claims_z14.txt) the fit holds on all 6,776 tuple lines; the evidence unit is the (ratio, curve) block:
  784 blocks of 196 ratios (the kappa.txt counts 313/451/930 are lines, 4 x 62/59/75 blocks). Tested, not proved
  (a symbolic proof would compute (x - e*) o psi / D(t) over Q(r)(t)).
- Theorem B tested (theoremB_z14.py, vs DFull bilinear + DFull direct; my own code: local points by DFull's search,
  rational p-adic lifts, PARI nfhilbert at every prime of K over 2 d beta z0 N(pi) and the e_i - e_j, signs at the
  real places, each finite place at two local points): interim 375 of 375 tested pairs agree (42 with value -1,
  necessarily d != beta) on 119 curves of 30 ratios (E_mix, E_mix with the first component at e3, E_08, E_LU);
  1,538 drawn pairs skipped because beta is not a norm from Q(sqrt d) (80 per cent of the draws). The local factor
  never changed between two local points. Final counts in mix_analyse.txt.
- Claims editor (EXPLORATION_3 §1, §4; EXPLORATION_4 §5), lines vs distinct cases (claims_z14.txt; my reruns):
  * "804 curves of the pilot": 201 ratios x 4 printed models of one curve = 160 distinct curves.
  * "-1 on all 6,776 mixed kill lines": right as lines; = 1,183 distinct (ratio, a, b), 1,009 distinct (curve up to
    isomorphism, a, b). "1,694 per curve": the same curve each time (four different tuple classes per line set).
  * "1,280 + 874 lines of E_08, E_LU" (2,154): 133 + 111 = 244 distinct cases.
  * "2,995 + 2,784 pairs, no disagreement" (tangent form of the E' conic vs cassels.py): seed 11 = 60 ratios = 57
    curves, seed 12 = 201 ratios = 160 curves, and seed 11's 57 curves are among seed 12's 160; 370 + 529 of the pairs
    are diagonal (d, d) (the "diag <d,d>" lines of the logs). Compared with cassels.py, not with rank_second.
  * "14,160 of 14,160 and 16,638 of 16,638": 747 + 13,413 and 902 + 15,736 (Task-2 sample + "extra" pairs); 84 and 105
    of the Task-2 parts are trivial (d = 1, beta = 1 or d = beta); same 57 / 160 curves as above, so not two
    independent samples; compared with cassels.py only. "it did not depend on the sign of s or on the solution (14
    solutions per pair)": forced by the arithmetic (my argument above), not evidence. The S*J(-1), S*J(2) controls are
    identically 1.
  * "2,248 of the first 2,995 pairs": right (747 + 2,248 = 2,995); the inference drawn from it is wrong (no Redei
    symbol is defined there; see above).
  * "2,760 / 488 / 3,528": lines; distinct (curve, a, b) 400 / 79 / 530.
  * "dimension 2 to 8 on the 201 ratios": 160 curves; on my 134 non-pilot curves dim V runs from 1 to 8 (1 on 41).
  * "201 of 201" Sel^2 dimensions: 160 curves.
  * §4 "37 / 34 / 62" and "36 / 47 / 41 / 26": reproduced exactly by rank0.py with H = 3000 (e3_rank0_rerun.txt), all
    150 ratios distinct (53 odd-odd), 150 distinct curve pairs; recomputed with rank_second (below).
    "38 of 80" (square type): NOT reproduced: rank0.py with H = 3000 gives 40 (bound 0 on E_A 24, on the twist 26,
    parities 19/23/18/20), and H in {1000, 2000, 5000, 10000, 13000} do not reproduce rank0.txt's square block either;
    the H used is not recorded. [fix: record the command line]
  * EXPLORATION_4 §5 kappa.txt (313 / 451 / 930 per curve): lines, = 62 / 59 / 75 (ratio, curve) blocks per curve name.
- explore3's conic_pairing (tangent form of the conic of w^2 = d t^4 + 8C^2 t^2 + 16A^2B^2/d; function text copied
  verbatim; conic_test.py) = DFull on 558 of 558 pairs (d != beta; 73 with value -1) on the 134 test curves, soluble
  and insoluble x^2 - d y^2 = beta z^2 alike. [note] That quartic is the class-d covering of the 2-isogenous curve
  Y^2 = X(X^2 + 8C^2 X + 16A^2B^2) (a' = -2a = 8C^2, b' = a^2 - 4b = 16A^2B^2: symbolic.txt), not the 2-covering of
  (1, d, d) on E_mix itself; "their covering is ..." should say so (the agreement is then the compatibility of the
  pairing with the isogeny, here tested, not proved).
- Task 4 wording (EXPLORATION_3 §4): the bound b = dim Sel^2 - 2 - rank(pairing) is the rank bound from Sel^2 and the
  pairing. b = 0 proves rank 0 (no hypothesis on Sha). b and the 2^infinity-Selmer corank have the same parity
  unconditionally (Sha/div carries a non-degenerate alternating pairing), so b odd gives rank >= 1 if Sha[2^infinity]
  is finite: the (odd, odd) sentence is correctly conditional. "one of the two curves has even parity and should have
  rank 0 or a large Sha" is not right as stated: for even b >= 2 (finite Sha[2^infinity]) the alternatives are rank 0
  with Sha containing (Z/4)^2, or rank >= 2; "should have rank 0" is the minimalist heuristic, not a consequence, and
  "or a large Sha" is the same branch, not an alternative. [fix] Also rank0.py counts a failed pairing evaluation
  (cassels returns None) as -1 (G = 0 if v == 1 else 1), which would raise the pairing rank and lower the bound,
  i.e. could fake rank 0; no failure occurred here (identical to rank_second, below), but the code should refuse. [fix]
- Task 4 numbers (rank0_z14.py with full_selmer2 + DFull; mix_analyse.txt): on rank0.py's own 150 ratios the
  (dim Sel^2, bound) of both curves equal explore3's bound() on 150 of 150, so 37 / 34 / 62 and 36 / 47 / 41 / 26 are
  confirmed by the repository's descent (DFull found no local point for 4 of 10 entries of the twist of (2279,1):
  there the bound 1 is an upper bound from the other entries). Square type with H = 3000 (80 ratios): 80 of 80 equal;
  24 / 26 / 40, parities 19 / 23 / 18 / 20 (not the note's 38). Independent draw (seed 6, H = 3000, 150 ratios, 48
  odd-odd; one curve pair in common with seed 5): 150 of 150 equal (one DFull failure on E_A of (1777,1776));
  bound 0 on E_08 for 43, on E_LU for 47, on at least one for 77 (51 per cent); parities (even, even) 40,
  (even, odd) 37, (odd, even) 41, (odd, odd) 32 (21 per cent). So "about five in six" with an even bound and "a
  quarter" both-odd hold roughly on the second draw too (118/150 = 79 per cent; 32/150).

## Bottom line (z14)
On 134 mixed curves outside the pilot (and 380 + 148 further curves for §4 and Theorem B) nothing disagrees with the
repository's rank_second/DFull: cassels.py 16,401/16,401, Redei form 863/863 soluble pairs (823 with d != beta, 169
with value -1), conic_pairing 558/558, Theorem B 454/454, rank0 bounds 380/380 ratios. The algebra of EXPLORATION_4
§5 and of Theorem B is correct (given Cassels' formula). What is overstated: the Redei form is still tested, not
derived (local factors above 2, infinity, d, beta, z0, AB not computed; Theorem B covers none of the pilot's kills,
which all have first component 2); "only the second reading is defined" is wrong; solution/sign independence and the
J(-1), J(2) controls are forced, not evidence; the four mixed curves are one curve and the counts are lines (6,776
lines = 1,009 cases); §2 needs two wording fixes; §4's "rank 0 or a large Sha" is mis-stated and its square-type
numbers are not reproducible. Final summary numbers: experiments/referee_z14/mix_analyse.txt.
