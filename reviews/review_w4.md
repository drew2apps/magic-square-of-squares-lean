# Referee report on w4: Fisher Theta_4 certificate for E+ (72,65)

Branch run4_w4, worktree /home/claude/conductor/w4, commit 2215a85 (diff vs 0a52347).
Referee started 01:02 UTC. Findings written as they come; verdict at the top is provisional until marked FINAL.

## Verdict (FINAL, 01:33 UTC)

**HOLDS WITH GAPS.** As far as it can be checked here, the certificate is correct. Theta_4(y0, y1) != 0 for E+ of (72,65),
so dim S_5 = 1 and dim S'_5 = 2, and rank E+(Q) <= 1 + 2 - 2 = 1. With G of infinite order, rank = 1. No input comes from Magma.
I found no error that changes the value or how it should be read. The gaps (none of them is an error):

- G1. Fisher's paper was read only through WebFetch summaries, by the builder and by me (three targeted fetches of the PDF, plus
  review_v4's). The statements used are the standard ones and agree across the fetches: Thm 2.1(ii), the Sect. 3 pushout formula,
  and the definition of a pushout form (div f = 2b with b representing T, and no specific divisor required). Lemma 2.2 is not even
  needed: the rank formula gives rank = dim E'(Q)/phi E(Q) + dim E(Q)/phihat E'(Q) - 2 <= dim S_5 + dim S'_5 - 2 = 1.
- G2. The result rests on earlier exact, Magma-free scripts that were refereed as HOLDS in review_v4.md and were not re-audited
  here: rank_indep (dim S_1 <= 2, dim S'_1 <= 3), rank_second (Theta_2 = 0), rank_fisher (ELS level-3 lifts, hence S'_4 = S'_1 by
  alternation) and rank_fisher5 (C3' and its conic points).
- G3. The evenness of div g0 and the vanishing of the local terms outside S are written arguments. I re-derived both line by line;
  they are correct and complete. I also machine-checked evenness by Weil reciprocity over F_10007 and F_100003: 0 failures in 37 pairs
  of hyperplane sections, while the odd-divisor control g0*c/e fails 17 of 37. This is not a proof over Q, and none is in the repo.
- G4. The external controls depend on the BSD order of Sha (rank 0 itself is proved). No control has Theta_4 != 0 as a theorem.
  The y1 of all four controls is 1 mod 8, so they never exercise a non-trivial 2-adic symbol. The main result does not need one
  either (see below).
- G5. (137,7) is covered only through the repository's (72,65)~(137,7) partner relation, which I did not check.

**Answers to the brief.**
1. Theory. The code's set-up is Sect. 3 literally. C4 is a torsor under E' with [C4] in Sel^(4)(E'), and it lifts y0 in S_4. The
   pushout class is T' = ker(phihat), and y1 = 2 lies in S'_1 = S^(phihat)(E'/Q) = S'_4. "div g0 = 2B with [B] = T'" is SUFFICIENT:
   B is rational, E'(Q)[2] = {0, T'} because b is not a square, and (I4) gives [B] != 0. Any other representative or scaling changes
   the sum by sum_v (c, y)_v = 0. S is correct: every p < 47, B0 (the pencil), the primes of Xi, the primes where R3^0 degenerates,
   the primes of the y, and R and 2. Outside S the written argument holds for every odd p not in B0 and prime to Xi. The global sum
   is independent of the choices: of P_v (the local images are exact annihilators), of g0 (Hilbert reciprocity), and of the lift [C4]
   (Thm 3.1(i), which needs y1 in S'_4, and it is). The PLACE of the -1 is not invariant: Xi = 274 puts it at 2, 5, 13, and K6 at 3.
   Moreover x(G) = 2 * 97505521 pairs to 0 (a theorem), so Theta_4(y0, 2) = sum_v (g0(P_v), 97505521)_v. Since 97505521 = 1 mod 8,
   that sum has NO 2-adic term. I computed it independently: -1, from the terms at 353, 97505521 and 252169271. The result does not
   hinge on the 2-adic symbol.
2. Re-runs and attacks. All six outputs reproduce byte for byte. My own local machinery reproduces every term at all 23 places for
   Xi = 137, and also -1 for a fourth twist Xi = 13358256377 that the repo never used. That machinery is independent of the repo's
   point finder and of its hilbert(): a different parametrisation, my own integer Hensel test, and PARI hilbert() on integers. The
   repo's hilbert() agrees with PARI's hilbert() on 24000 random pairs. Low 2-adic precision is refused, never mis-classified.
   Perturbed 2-adic points change the class in 176 of 240 cases, and the Hensel test rejects all 240. Mutations of the conic point,
   k_1, sigma_1, tau_3 and q are refused, while valid alternatives (R3', N1', 3 k_1) give the same value, as theory predicts.
3. Controls. A, B (BSD Sha 16, non-zero) and C, D (BSD Sha 64, zero) are what they are claimed to be (PARI). They have the same
   2-torsion structure, the same side and the same cone fields. With my own local machinery they give -1, -1, +1, +1, as predicted.
   These four are the strong controls, though conditional on BSD. K2 is medium. K1, K3-K6 and the recheck are consistency checks only.
4. Magma. None. The point G is from PARI. Every finder is PARI and every finder output is re-checked exactly.

## Findings

FIX-FIRST: none.

FIX
1. rank_fisher_theta4_control.py:112-113 (printed at rank_fisher_theta4_control.txt:5, 40, 75, 111). The sentence
   "a lift with an empty one has <[C3], y1> != 0, Sect. 3" is wrong-sided and false here. The obstruction for C3' is <[C3'], y> with
   y in S_1 (the E' side). Because Theta_3 is alternating, every ELS level-3 lift of y0 has a NON-EMPTY level-4 Selmer set. Control
   D's "0" entries (lam = 2, xi = -1, -2) are misses of the real sampler (my re-run: classes at every finite place, no sampled point
   of C3'(R)). Replace with: "(the first lift for which the sampler found a level-4 twist; by Thm 2.1(i) every ELS level-3 lift of
   y0 has one; 0 means none was found by sampling)".
2. rank_fisher_theta4.py:824. Add y = 2 x(G) = 97505521 (mod squares) to ys and print its product: -1, with no 2-adic term. This makes
   it visible that the -1 does not rest on the single 2-adic symbol.
3. reports/w4.md:157-174, the proposed README and LEDGER text. Say that the external controls are BSD-conditional, name the inherited
   scripts (G2), and give the review status (one AI referee, this report) in place of "unreviewed".
4. rank_fisher_theta4.py:45-55. Add a machine check of evenness next to the written argument, for example the Weil-reciprocity test of
   /tmp/rw4/indep/weil.py at two primes with its odd-divisor negative test.
5. rank_fisher_theta4.py:399-436 (real_cds, real_cds_from_Dlam, real_signs). The real sampler only searches near D_lam's certified real
   point, so it misses C3'(R) for lifts on the other branch (control D). Sample the whole real locus; choosing (e:g) first, as my
   script does, finds real points at once.

NIT
1. rank_fisher_theta4.py:843 (and the last line of the .txt). "S_5 = <images of E(Q)>" should read "<image of E'(Q)> = <-1>", because
   S_5 lies in S_1 on the E' side.
2. rank_fisher_theta4.py:56. Add the rule at p = 2: (g0, 2)_2 = (-1)^((w^2-1)/8) with w the odd part. All four 2-adic classes found
   have w = 3 or 5 mod 8.
3. rank_fisher_theta4_recheck.py:100 uses T.hilbert, the main script's own symbol code. Use pari.hilbert or an independent formula.
   (The two agree on 24000 random pairs.)
4. rank_fisher_theta4.py:26 and :52. Say that B is Q-rational because it is half of a rational divisor. Say also that 4[C4] = 0
   follows from the chain of phi and phihat double covers (Kummer theory), which is all of Thm 7.3 that is used.

## Not checked
- The other-side computation, Theta_4 = <y0, eta_1> with a 4-covering of E lifting y1. Its cone conic lies over
  Q(sqrt 15513721419678722), and no solver for that is at hand.
- Fisher's printed text (only summaries were read). Fisher's Examples 10.2/10.3 as controls.
- The internals of rank_indep, rank_second, rank_fisher and rank_fisher5 (see review_v4.md).
- The (72,65) ~ (137,7) relation.
- K4 and K6 with my machinery. I re-checked only lam = -1 (Xi = 137 and 13358256377) and the four external controls.
- A proof over Q of evenness beyond the written argument plus the checks at two primes.

## Log

### 01:10 UTC - Item 1, the theory (Fisher read by WebFetch of the PDF, three targeted fetches; summaries, as for the build agent)

What the paper says (Sect. 2, 3, 7, as returned by the fetch):
- S_m (m = 2n) = image of S^(2^n)(E'/K) in S_1 = S^(phi)(E/K); S'_m the same with E, E' swapped, inside S'_1 = S^(phihat)(E'/K).
  Thm 2.1(ii): Theta_4 : S_4 x S'_4 -> F_2, left kernel S_5, right kernel S'_5. Lemma 2.2: rank <= dim S_{2n-1} + dim S'_{2n} - dim E(K)[phi] - dim E'(K)[phihat].
- Sect. 3: Theta_m(xi, eta) = <xi_1, eta> for a lift xi_1 in S^(2^n)(E'/K); independence of the lift is Thm 3.1(i) (needs eta in S'_m).
  <x, y> for x in S^(psi)(E'/K), y in S^(phihat)(E'/K): C = covering curve of x (a torsor under E'), T generator of E'(K)[phihat],
  b in Div^0(C) K-rational representing T, f with div f = 2b, <x, y> = sum over ALL places of (f(P_v), y)_v, P_v in C(K_v) off the
  zeros/poles of f; only finitely many terms non-zero.
- Sect. 7, Lemma 7.1: the tangent-plane product l_1 l_2 is "a" pushout form; the definition of pushout form is only
  div(f/x_1^2) = 2b with b representing T. No specific divisor is required.

Referee's derivation of what must hold (my own, checked against the code's notation at rank_fisher_theta4.py:12-26):
- The code's set-up matches Sect. 3 literally: C4 is a torsor under E' (4-covering of E'), x = [C4] in S^(4)(E'/Q) lifts y0 in S_4,
  the pushout class is T' (= Sect. 3's T, generator of ker phihat), and y = y1 = 2 in S'_1 = S^(phihat)(E'/Q). The twisting group of
  C4 -> C3' (S'_1, the 8-element coset 137*S'_1 in the output) is the same group as the pairing argument's group; I derived this
  independently (C4 -> C3' is compatible with phihat: E' -> E, kernel E'[phihat]), so the coset "137 * S'_1" is the expected shape.
- "div g0 = 2B with [B] = T'" IS SUFFICIENT. Reason: B is K-rational (half of a rational divisor), deg 0, 2[B] = 0, so [B] is in
  Pic^0(C4)(Q)[2] = E'(Q)[2] = {0, T'} (b is not a square, so X^2 - 2aX + a^2 - 4b is irreducible); [B] != 0 by (I4). Any other b' in the
  class of T' changes f by c*h^2, which changes each local term by (c, y)_v and the sum by sum_v (c, y)_v = 0 (Hilbert reciprocity).
  Fisher's Lemma 7.1 divisor is one construction, not a requirement.
- Independence of P_v: g0(P)/g0(Q) lies in Im(delta_phi,v) (the Kummer map of the isogeny dual to the one whose kernel is [B]); y1 in
  S'_1 lies in Im(delta_phihat,v); these are exact annihilators under the Hilbert symbol, so the term is constant on C4(Q_v). K1 tests this.
- Independence of the lift [C4] (the twists Xi): Thm 3.1(i), needs y1 in S'_4. y1 in S'_4: rank_fisher.txt:71-86: y1 has ELS level-3
  lifts (S'_3 = <b, x(G), y1>), and Theta'_3 is alternating with b, x(G) in its kernel, so Theta'_3 = 0 and S'_4 = S'_3. Correct argument
  (it rests on rank_fisher.py's level-3 ELS checks and rank_indep's first-descent upper bounds dim S_1 <= 2, dim S'_1 <= 3, not re-audited here).
- Rank bound: dim S_5 <= 2 - 1, dim S'_6 <= dim S'_5 <= 3 - 1, so rank <= 1 + 2 - 2 = 1. Correct.
- The "single -1 at p = 2" question: the LOCATION of the -1 is not invariant (Xi = 274 puts it at 2, 5, 13), only the product is. A
  stronger observation: since x(G) = 2 * 97505521 pairs to 0 (a theorem: images of E(Q) lie in the right kernel), Theta_4(y0, 2) =
  sum_v (g0(P_v), 97505521)_v, and 97505521 = 1 mod 8 is a 2-adic square and positive. So the answer does NOT depend on the 2-adic term at
  all: from the printed terms it is the product of the odd-prime terms at 353, 97505521, 252169271 (-1 each) = -1. The 2-adic term enters
  only through the K2 check for x(G) (whose 2-adic term equals the one for 2). I re-check this independently below.
- Evenness argument (docstring lines 45-55) re-derived line by line: correct, given (I1), (I2), (I3), (I3'). Note div q = A_1+...+A_4
  exactly (degree 8), which the W = 0 case uses.

### 01:17 UTC - Item 2, re-runs and an independent local re-computation

- All six outputs reproduce byte for byte (rank_fisher_theta4.py 41 s, --negative, _k4, _k6, _control, _recheck), apart from the
  PARI stack banner line on stderr. Negative tests 9/9 refused as recorded.
- Referee's own script /tmp/rw4/indep/ref_local.py (outside the repo): uses the repo only to BUILD C3', C4 (Xi = 137) and g0; finds
  local points by a different parametrisation ((e:g) chosen first, (c:d) from beta*PA - alpha*PB = 0; the repo chooses (c:d) first),
  validates each point by its own integer Hensel test (v(F) >= 2 v(J) + 1 over 3x3 minors), checks both (I1) identities at every point
  to the precision, and evaluates the symbol with PARI's hilbert(num*den, y, p) on exact integers (NOT the repo's hilbert()).
  Result at all 23 finite places of S, 4-40 certified points each, every term point-independent:
    y = 2:          -1 only at p = 2; product -1            (matches rank_fisher_theta4.txt)
    y = b-class:    -1 at 353, 1279, 336527, 97505521; +1  (matches)
    y = x(G)-class: -1 at 2, 353, 97505521, 252169271; +1  (matches)
    y = 97505521:   -1 at 353, 97505521, 252169271; product -1 (new: 97505521 = 2 * x(G) mod squares, a 2-adic square, so this
                    form of Theta_4 has NO 2-adic term; the -1 is carried by three odd primes)
  So the value -1 does not hinge on the 2-adic Hilbert symbol or the 2-adic precision. The repo's hilbert() also agrees with PARI's
  hilbert() on every term computed here (same terms).
- Mutations run by the referee through ref_local.py (same independent machinery):
    R3 -> R3' in g0 (g0' = 2 (delta3/N1)^2 / g0, another pushout for T'):  y = 2 product -1 (predicted: same value)
    N1 -> N1' in g0 (g0'' = -2 g0 * square, another pushout):             -1 (predicted: same)
    k_1 -> 3 k_1 (a different valid norm solution):                        -1 (predicted: same)
    a fourth twist Xi = 13358256377 (not used by the repo):                 -1, K2 classes +1
  (Xi = 274 with the referee's (e:g)-first sampler found no point at 710053 in 8 seeds; the repo finds them with its (c:d)-first
  sampler and its Hensel re-check certifies them. Sampler weakness, not a finding.)
- The repo's hilbert() (rank_fisher_theta4.py:302-313) agrees with PARI's hilbert() on 24000 random pairs at p = 2, 3, 5, 7, 11,
  13, 353, 97505521: 0 mismatches. The 2-adic formula (eps(u)eps(v) + a*om(v) + b*om(u)) is the standard one.
- Mutation coverage of the repo's own negative tests (re-run, all refused): conic point P1 / P3 moved off the cone, non-tangent line,
  f perturbed in step (8) (rank_fisher5 --negative, 9/9); q perturbed, sigma_1 doubled, tau_3 tripled, k_1 of wrong norm, forged
  twist, y = 3, non-pushout g0*c/e, 2-adic precision 4 (rank_fisher_theta4 --negative, 9/9). A mutated g0 is caught by (I1); a mutated
  conic point by rank_fisher5. Good coverage.
- Context check: PARI ellrank(E+) = [1, 3] (the 2-descent bound 3 that Theta_4 is meant to lower); b, a^2 - 4b non-squares, so
  E(Q)[2] = E'(Q)[2] = Z/2 and [B] in {0, T'}; bad primes of E = {2,3,5,7,13,97,137,79553041,97505521} = the list used to certify
  y1 = 2 in S'_1 (rank_fisher_theta4.py:813-818). 710053 is bad only for the model of C3' (non-minimal pencil), and is in S.
- The degenerate primes 353, 1279, 336527, 252169271 (R3^0 reduces into (QA, QB)) are primes of GOOD reduction of E (a, b, a^2-4b
  non-zero mod p) and recur for every lift and conic point tried (K4, K5, K6). They are in S and their terms are computed. Without
  them the answer would change (they carry the -1's of the 97505521-form above).

### 01:21 UTC - Item 1 (cont.), the bad set S and the written good-reduction argument

Re-derived (docstring lines 58-78): for odd p outside B0 and not dividing Xi, (i) C3~ smooth (pencil discriminant a unit; the
binary-quartic discriminant via poldisc of the formal degree-4 polynomial is the right invariant even if the leading coefficient
vanishes mod p); (ii) f0 = unit * l1' conj(l1') with l1' locally primitive at the primes of L over p (Gauss lemma), l1'~ tangent to
a reduced cone along a ruling, so div f0~ = 2(A1~ + A2~) with A1~ + A2~ - H~ ~ T~ != 0 and the cover is connected and etale;
(iii) N1^0 and R3^0 do not vanish identically on C4~ (the case split on v_p(r) is complete: >= 0 by the "f0 would be a constant times a
square" argument, < 0 reduces to Q0~ in span(QA~, QB~), and those primes are put in S); (iv) at most 8 + 16 + 8 = 32 bad points,
p + 1 - 2 sqrt p > 32 for p >= 47; points with W != 0 are smooth on the model (Jacobian rank 3), so they lift. So the term vanishes
outside S. The argument is correct and complete. ELS outside S follows from the same count.
Real place: all y used are positive, so the real term is +1 whatever the sign of g0; the interval-arithmetic sign is not needed.

### 01:22 UTC - Item 3, controls

Independent PARI (cypari2) check of the four external curves (script /tmp/rw4/indep/ctl.py):
| curve | a^2-4b | conductor | tors E, E' | ellrank E, E' | an. rank | L(E,1) | BSD Sha |
|---|---|---|---|---|---|---|---|
| A (128, 4097) | -4 | 262208 | Z/2, Z/2 | [0,2], [0,2] | 0 | 6.282898 | 16.000000000 |
| B (74, 17753) | -2^16 | 284048 | Z/2, Z/2 | [0,2], [0,2] | 0 | 9.711180 | 16.000000000 |
| C (1310, 429041) | -2^6 | 3432328 | Z/2, Z/2 | [0,2], [0,2] | 0 | 15.712169 | 64.000000000 |
| D (1706, 743993) | -2^16 | 11903888 | Z/2, Z/2 | [0,2], [0,2] | 0 | 13.711033 | 64.000000000 |
All claims about them hold. Same shape as E+ of (72,65): one rational 2-torsion point on each side, b and b' non-squares, cone
fields Q(sqrt -2), Q(sqrt 2), y0 = 2 on the E' side (S_1 = {+-1, +-2}) lifted to C3' -> C4 (a torsor under E'), y1 on the E side.
Same side of the isogeny and the same code path.
The predicted answers, re-derived: |Sha(E)[phi]| = |Sha(E')[phihat]| = 2, so |Sha(E)[2]| <= 4; BSD order 16 forces Sha(E)[2^oo] =
(Z/4)^2 and 64 forces (Z/8)^2 (an alternating non-degenerate pairing gives M x M). y0 in S_5 (image of Sel^(4phi)(E), via
[4]: E[4phi] -> E[phi]) iff iota(y0) is in 4 Sha(E). A, B: 4 Sha[2^oo] = 0 and iota(y0) != 0, so Theta_4(y0, .) != 0 on S'_4 = <b, y1>,
hence Theta_4(y0, y1) != 0. C, D: iota(y0) is in Sha[2] = 4 Sha[2^oo], so Theta_4(y0, y1) = 0. Correct.
Strength ranking:
- STRONG (an answer known for a reason independent of the code, both values occur): ext A, B (non-zero), C, D (zero). They are
  conditional only on the BSD order of Sha (rank 0 is a theorem). A code that returns a constant, or one that drops the degenerate
  primes (337, 239, 7/1193, 449/736823 there), fails at least one of them.
- MEDIUM (theorems, non-trivial local terms): K2 (b, x(G) pair to 0, with four -1 terms each). Through the identity
  Theta_4(y0,2) = Theta_4(y0, 2 x(G)) it ties the 2-adic term to the odd-prime terms, so it is not merely a "trivial pairing" check here.
- CONSISTENCY ONLY (same code path, different choices; would not catch a systematic error that is the same for all choices):
  K1, K3, K4, K5, K6 and the Hensel re-check (the re-check validates the PARI p-adic points independently, but uses the repo's
  hilbert(); the referee's run above replaces that too).
- None of the controls has Theta_4 != 0 predicted by a proved theorem: the non-zero predictions (A, B) are BSD-conditional. The
  agent's own search found no rank-1 control of the same shape. Fisher's Examples 10.2/10.3 were not used (other cone fields).

### 01:25 UTC - Item 2 (cont.), precision and perturbed local points (script /tmp/rw4/indep/perturb2.py)
- 2-adic precision sweep with the repo's own local_terms (Xi = 137): prec 5, 6, 8, 12 -> refused at p = 2 (no certifiable point);
  prec 20, 60, 400 -> the same classes {(0,3),(0,5),(1,3),(1,5)} and term -1. At p = 353 the term for 97505521 is -1 from prec 6 up.
  So low precision is refused, never mis-classified.
- 240 perturbed 2-adic points (one coordinate + 2^k, k = 1, 2, 3, 5): 176 change the class of g0, and 0 of 240 pass the integer
  Hensel test used by rank_fisher_theta4_recheck.py. The re-check would catch a bad local point. (The main script does not read stored
  points: it recomputes them with PARI p-adics each run, so "mutating a stored local point" is not possible there.)

### 01:24 UTC - Item 4, independence from Magma
- No script on the certified path calls Magma; grep finds "magma" only in prose. No JSON input mentions Magma.
- The known point G comes from sigma_evidence.txt:13 with "gen=PARI,_saturated". That line also records Magma's rank bound, but no
  code reads that field. G is used only as a point (lower bound, and the class x(G)).
- Finders are PARI: qfsolve for norm equations and conics, the Legendre descent for the cone-conic points, p-adic sqrt for local points.
  Every output is re-checked exactly (norms, membership in the ideal, tangency, Hensel in integers). PARI's matrank over Q is used for
  exact linear algebra (I1), and that is exact.
- Inherited inputs, all Magma-free and refereed in review_v4.md (HOLDS): rank_indep's first descent (dim S_1 <= 2, dim S'_1 <= 3),
  rank_second's Theta_2 = 0, rank_fisher's ELS level-3 lifts (hence S'_4 = S'_1 by alternation), rank_fisher5's C3' and conic points.
Verdict on item 4: independent of Magma. Not self-contained: it rests on those four earlier scripts.

### 01:25 UTC - Item 1/2, a machine check of the written evenness argument (referee's /tmp/rw4/indep/weil.py)
Weil reciprocity: if div g0 = 2B then for two hyperplane sections D1, D2 of C4 off the support, g0(D1)/g0(D2) = (l1/l2)(B)^2,
a square. Over F_10007 (a prime of good reduction outside S), g0(D) is computed as the norm of g0 in the degree-8 algebra
F_l[C4 n {l = 0}] (Groebner basis and the determinants of the multiplication maps; the chart c = 1, and sections that are not 8
affine points are skipped). Result: 12 random pairs, g0 gives 0 non-square ratios (evenness predicts 0); the odd-divisor control
g0*c/e gives 8 of 12 non-square ratios. This is independent of (I1)-(I3') and of the written argument. It is a check mod one prime,
not a proof over Q, but an odd multiplicity over Q that survives reduction would show up. A second prime is running.
  Second prime: F_100003, 25 pairs: g0 0 non-square ratios; g0*c/e 9 of 25. Evenness confirmed at two primes.

### 01:29 UTC - Item 3 (cont.), the Selmer-set sampler and a wrong-sided sentence in the control output
- Theory (referee): for every ELS level-3 lift C3' of y0, <[C3'], y0> = Theta_3(y0, y0) = 0 (Thm 2.1(i), alternating, any lift) and
  <[C3'], -1> = 0 (-1 is the image of T'), so the level-4 Selmer set S(C3', f) is NEVER empty (Fisher Sect. 3). The pairing group for
  C3' (a torsor under E) is S_1 (E' side), not S'_1.
- Control D lists lifts (lam, xi) = (2, -1) and (2, -2) with "0" level-4 twists. Referee's re-run (/tmp/rw4/indep/ctlD.py, 6 x 300
  samples per place): every finite place gives classes, but the real sampler finds NO point of C3'(R) ("real signs []"), so choose_Xi
  rejects every Xi. Since -1 is in S_1 = Im(delta_phi) everywhere, these twists are ELS; the "0" is a sampler miss at R, not an empty
  Selmer set. The same happens for K6's lam = -79553041 at p = 97 (the agent says so). Soundness is not affected: only points actually
  found and checked certify anything.
- rank_fisher_theta4_control.py:112-113 (and so control.txt:5, 40, ...) prints "(the lift whose level-4 Selmer set is non-empty; a lift
  with an empty one has <[C3], y1> != 0, Sect. 3)". This is wrong twice. It has the wrong side: the obstruction for C3' is <[C3'], y> with
  y in S_1 (E' side), not y1 (E side). And by alternation no ELS lift has an empty set. Text fix only.

### 01:33 UTC - small checks
- G (sigma_evidence.txt:13) is on E, of infinite order (PARI ellorder 0), and x(G) mod squares = 195011042. The class used is right.
- C4(R), Xi = 137: 20000 real points of C3' sampled by the referee's (e:g)-first parametrisation all have f0 > 0, so Xi f0 > 0 has
  real solutions (mpmath at 80 digits, not interval arithmetic; the repo's interval-arithmetic point is the certificate).
- Other-side cross-check (Theta_4 = <y0, eta_1> with a 4-covering of E lifting y1): NOT DONE. The E-side level-3 curve (lam = -7,
  xi = 35) has cone fields Q(sqrt 15513721419678722) = Q(sqrt(2*79553041*97505521)) and Q(sqrt 2). The repo's conic finder only
  handles Q(sqrt +-2), and the referee did not attempt a conic over the large real quadratic field.

### 01:36 UTC - external controls re-evaluated with the referee's machinery (/tmp/rw4/indep/ctl_indep.py)
The repo builds C3', C4, g0 for controls A and C. The referee's own (e:g)-first points, integer Hensel test and PARI hilbert() give:
A: product -1 for y1 = 17 and +1 for b (predicted -1, +1); C: +1 for y1 = 409 and +1 for b (predicted +1, +1). Points were found at
every place of S. All y are positive, so the real place is trivial. B and D are running.

### 01:33 UTC - controls B, D with the referee's machinery
B: y1 = 41 gives -1 and b gives +1; D: y1 = 449 gives +1 and b gives +1, as predicted. The referee's sampler found no 2-adic point on these two,
but every y there is 1 mod 8, so the 2-adic term is +1 whatever the point and the products are complete. All four external controls
reproduce independently. Note that all four y1 (17, 41, 409, 449) are 1 mod 8: no external control exercises a non-trivial 2-adic symbol.
