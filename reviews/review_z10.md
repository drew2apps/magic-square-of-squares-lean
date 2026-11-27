# z10 — hostile review of t9 (commit df1bc35, branch run7_t9): Θ₄ on E' of (945,248)

Written by an AI model; not read by a person.

Start: Fri Nov 27 16:04 UTC 2026. Hard stop 19:00 UTC. Scratch: /tmp/z10/. No repository file edited.

## 16:12 UTC - item 5 (regression), first pass
Clean clone of /home/claude/conductor/t9 at df1bc35 in /tmp/z10/clone; sympy 1.14.0, PARI 2.17.4 (cypari2).
Script /tmp/z10/ci/run_ci.sh: the 13 existing CI lines of .github/workflows (block "run 6: Theta_4 by the general
module") literally, plus 945,248 --side=E' main and --negative, plus rank_fisher_theta4_norm.py (default mode), each
diffed against the committed output (PARI stack notice filtered as CI does). Result (/tmp/z10/ci/ci.log): every
diff EMPTY, every exit 0. Times: 72,65 13/6 s; 920,129 17/6; 760,39 21/7; 631,71 E' 20/15; _control 24; _recheck 65;
_order 1; _weil 54; _order --negative 1; 945,248 E' 43 (main) / 44 (--negative); norm.py 1 s.
So the edits to rank_fisher_theta4_gen.py change no existing output (as far as these runs show).

## 16:22 UTC - item 1: own local terms (referee script /tmp/z10/work/indep.py, outside the repo)
The clone is used ONLY to build the objects under test: QA, QB, f0 of C3' (xi = 17794, via prepare()), and N1^0, R3^0
of g0 (build_g0 + normalise) for a GIVEN Xi and GIVEN k_1, k_3; check_g0 ((I1), (I2), (I3')) is re-run on each.
Written here: two samplers ((e:g)-first and (c:d)-first, random p-adic integers to ~10^-260, own p-adic square roots
by Tonelli-Shanks + Newton, own 2-adic root); an integer multivariable Hensel certificate (fix 2 coordinates, 3x3
minor of the Jacobian of (QA, QB, Xi W^2 - f0): v(F) > 2 v(det), point known mod p^(v(F)-v(det)), and the valuations
of R3^0, W, N1^0 at the point at least 1 (3 at p = 2) below that); own Hilbert-symbol formula, asserted equal to PARI
hilbert() on every term; own place set: bad(E) from own factorisation of 2 b (a^2-4b), pencil-quartic discriminant
(1086 digits) factored with every factor proved prime (isprime), f0-degenerate primes (gcd of 3x3 minors of
[QA, QB, f0]), a DIRECT R3^0-degeneracy test on C4 (gcd of the 4x4 minors of [QA, QB, Xi W^2 - f0, R3^0] over the
15 quadratic monomials), N1^0 and R3^0 primitive, primes of Xi and of the y's, every p < 47; plus every prime
47..300 outside S as a sanity check. Real place: a real point of C4 by mpmath interval arithmetic (QA, QB enclosures
contain 0, Xi f0 > 0 certified); all paired classes are > 0, so the real term is +1.
My place set equals the repo's S exactly (26 primes); the direct 4x4 test gives den(r)-primes plus B0 primes, as
the repo's den(r) test. 12-16 certified points at every place, terms constant on them.

| config | Xi | k_1, k_3 | y = 2: -1 at | product | b, x(G) |
|---|---|---|---|---|---|
| repo | 521375088721 = 349*1493911429 | repo's | 2, 29, 349, 32746027, 45237653, 1493911429, 100067394997 (= repo) | -1 | +1, +1 |
| ownk | 521375088721 | conj(k1)(3+2 sqrt D1)^2, conj(k3)(5 - sqrt D3)^2 | 2, 13, 29, 61, 349, 32746027, 45237653 | -1 | +1, +1 |
| xm | 2 (not found by the repo) | own (qfsolve with a random square factor, times a square) | 2, 3, 5, 61, 349, 32746027, 45237653 | -1 | +1, +1 |
| x2 | 1042750177442 = 2*Xi_repo | own | 3, 5, 29, 61, 349, 32746027, 45237653, 1493911429, 100067394997 | -1 | +1, +1 |
| x1 | 1 | own | 29, 61, 349, 32746027, 45237653, 1493911429, 100067394997 | -1 | +1, +1 |
| xb | 678244410280322944874881 (= b) | own | 29, 61, 349, 32746027, 45237653, 1493911429, 100067394997 | -1 | +1, +1 |
| xg | 2601752269922 (= x(G)) | own | 2, 3, 5, 13, 29, 61, 349, 32746027, 45237653 | -1 | +1, +1 |
Every one of these C4 is shown everywhere locally soluble by MY certified points at all 26 places of S and at R
(good reduction elsewhere: p odd, outside B0, p not dividing Xi). The value -1 is reproduced with 6 level-4 twists
(5 not used by the repo) and 7 choices of norm elements; with and without a -1 at p = 2 (y = 2 is not 1 mod 8, so
the 2-adic symbol is exercised). VERDICT item 1: HOLDS.

Finding (sampler, f0-degenerate primes). f0 mod p lies in span(QA, QB) at 2, 3, 5, 13, 31, 97, 157, 349, 977, 245411.
At 245411 a stratified search over all 2(p+1) residues (/tmp/z10/work/strat.py, 171 s) gives the second unit class
of f0 at 4 of 490828 sampled points, all at one residue (c:d) = 103578 mod p; uniform sampling never sees it. This is
why the repo's sampler finds only 2 of the 8 twists Xi in S'_1 = <2, 349*1493911429, 13*100067394997> (the twists
form the coset S'_1 itself, since Xi_repo is in S'_1): only those with Xi_repo's unit class at BOTH 97 and 245411
survive. My runs above use that residue as one seed at 245411 (marked in the script), and random points elsewhere.

FINDING (negative test, 945,248 --negative): the refusal "a forged twist Xi q outside the Selmer coset (no point of
1042750177442 W^2 = f0 over Q_97 (sampled classes))" is VACUOUS. (i) The repo's sampler T.f0_classes finds NO point
of C3'(Q_97) with seeds 1, 2, 3 (/tmp/z10/work/forged.py), so the same assertion also "refuses" the genuine Xi =
521375088721 (its class (0,1) is not in the empty closure). (ii) 2 is a square mod 97 (14^2 = 196 = 2 mod 97), so
Xi*2 and Xi have the same class at 97. (iii) Xi*2 = 1042750177442 is NOT forged: it is in S'_1 and C4 for it is
everywhere locally soluble by my certified points (row x2), with Theta_4 = -1. The forged-twist test (already noted
as about the sampler in review_y5 NIT 1) here tests nothing. Not an error in the result; the negative output's
claim "9 of 9 refused" should be read as 8 meaningful refusals. Fix (text only, or code): pick a forged Xi outside
S'_1 (e.g. Xi*3, 3 not in the local image at 2 or 3) and use the seed-retry sampler in forged().

## 16:25 UTC - item 3: the conic points (hints)
- rank_fisher5.analyse steps (5)-(8) (read): exact sympy arithmetic in Q(t), t^2 = dU (resp. dV): the singular
  member of the pencil at a root w0 in L, block structure and ranks; (6) P1 on the cone, not the vertex, the tangent
  plane l_1 = grad Q_w(P1), rank <= 1 on {l_1 = 0} (2x2 minors), f = l_1 conj(l_1) rational; (7) the same for P3;
  (8) q through the 8 points by an exact nullspace (3-dimensional, one vector outside span(QA, QB)) and c by an exact
  linear solve. No PARI call there; nothing GRH-dependent. gen.py imports no bnfinit/rnfisnorm/bnfisprincipal (grep);
  the GRH-dependent calls are all in rank_fisher_theta4_norm.py's finder (bnfinit flag 1, bnfisprincipal, rnfisnorm),
  whose output (two points over L, L') is re-checked as above. The default mode of norm.py prints the class groups
  "for information" only. So no GRH assumption enters the result. HOLDS.
- Own exact check (/tmp/z10/work/check_hints.py, 2 s): det(wA + B) splits over Q into two quadratics; the U point is
  on the cone Q_w for a root w in Q(sqrt 5820647), the V point on a cone over Q(sqrt -5820647); both cones have rank 3;
  the planes are tangent along a ruling (restricted form rank 1); N(l_1), N(l_3) are rational and EQUAL to the repo's
  f, f_dag (ratio 1); and an explicit certificate f f_dag - c q^2 = QA mA + QB mB (PARI matinverseimage, re-expanded
  exactly to 0 in sympy), c of 389/97 digits. (q through the 8 points: the repo's exact check, not redone.)
- `rank_fisher_theta4_norm.py --find` (clean clone, 24 s): "equal to the committed hint: True" for U (13.4 s) and V
  (9.4 s); nothing written without --write (git status clean). `--obstruction` (130 s) reproduces the committed
  rank_fisher_theta4_norm_obstruction.txt byte for byte. Default mode: reproduces rank_fisher_theta4_norm.txt.
VERDICT item 3: HOLDS.

## 16:28 UTC - item 2 (part): a second level-3 lift in the liftable coset (K4-type control, not in the repo)
In a scratch copy (/tmp/z10/k4, not the repository) t9's finder found conic points for the level-3 lift
xi'' = -17794 (= 17794 * (-1), in the coset 17794 S_2 that the obstruction table predicts liftable; U 15 s, V 24 s),
and rank_fisher_theta4_gen.py 945,248 --side=E' ran on it unchanged (the hint '|xi' = -17794): twists found
[9074, 298983650112551441426]; Theta_4(-5820647, 2) = 1/2 (product -1; -1 at 3, 5, 61, 32746027, 45237653,
1493911429, 100067394997), K2 +1, +1, K3 (-1 at 2, 3, 5, 32746027, 45237653), K5 -1. My own local terms on that C3'
(indep.py with Z10_REPO=/tmp/z10/k4/...): Xi = 9074 (repo k): -1 at the same 7 places, product -1; Xi = 4537
(= 9074/2, own k_1, k_3): -1 at 3, 32746027, 45237653, product -1; K2 classes +1; ELS shown at all 26 places and R.
So the value is the same through a second level-3 lift over the same D_lam (Fisher Thm 3.1(i) independence, as
observed for (72,65) by its K4). (Xi = 1 is not a twist of this C3': no solution of the norm equation.)

## 16:30 UTC - item 2: the level-3 lift xi = 17794 instead of the json's 1271
Claim restated: y0 = -5820647 in S_4 (E' side, S_1 = S^(phi)(E/Q), classes d | a^2-4b); the json's level-3 lift
C3(lam = -349, xi = 1271) has an empty level-4 Selmer set; C3(lam = -349, xi = 17794 = 1271 * 14) does not; any
everywhere locally soluble level-4 curve over any level-3 lift of y0 may be paired with y = 2 in S'_4.
- Fisher, Sect. 3 (ar5iv HTML of arXiv:1509.03234, read through a WebFetch summary, as before: gap G1): for m = 2n
  even, xi in S_m, eta in S'_m, "choose xi_1 in S^(p^n)(E'/K) lifting xi ... Theta_m(xi, eta) = <xi_1, eta>",
  independent of the choice of lift (Thm 3.1(i) with psi = p^(n-1) phi); Thm 2.1(ii): for m even a pairing
  S_m x S'_m with left kernel S_(m+1) and right kernel S'_(m+1); for m odd alternating pairings on S_m x S_m and on
  S'_m x S'_m. So for m = 4 ANY element of S^(4)(E'/Q) mapping to y0 may be used; it need not lie over the json's C3.
  The hypothesis on the paired class is y in S'_4 (needed for the independence), inherited from rank_fisher.py.
- Which lift is paired: [C4] with C4 = {QA = QB = 0, Xi W^2 = f0}, Xi = 521375088721, over C3'(xi = 17794), over
  D_(-349), over y0. Only the third equation of C3 changes (xi w^2 = l_X l_Y); lam, U, V, l_X, l_Y are the json's
  (read in prepare(): C = dict(C, xi=str(xi2))). 14 is in S_1 of this side (rank_indep's Sp = <-1, 14, 17, 97826>;
  the assertion in prepare() tests membership in RD['Sp'] of first_descent(a, b), which is the E'-side group:
  checked). C3'(17794) is everywhere locally soluble because C4 is (its local points map down).
- Is C4 checked locally soluble at every place of S and at R? Yes: local_terms asserts a usable point at every
  p in S (26 primes), real_term asserts a real point; my own certified points (item 1) at every place of S and R for
  this C4 and five other twists. Outside S: the good-reduction argument (p odd, outside B0, p not dividing Xi; S has
  every p < 47) as written in rank_fisher_theta4.py (inherited, gap G2); the f0-degenerate primes
  [2, 3, 5, 13, 31, 97, 157, 349, 977, 245411] are all in B0, as that argument predicts.
- The obstruction table (norm.py --obstruction, a sampler; reproduced byte for byte) is consistent with Thm 2.1(i):
  Theta_1 alternating of rank 2 on S_1/S_2 with Theta_1(14, 17) = -1 predicts exactly the four rows (I checked the
  products by hand: Ob_(xi s) = Ob_xi * Theta_1(s, .) for s = 14, 17, 238). The symbol Ob(eta) is computed from
  sampled classes only, but by local duality (eta in S_1 annihilates the local image of E(Q_p)) it is constant on
  the true coset, so the sampler limit at 97, 245411 does not affect the table. Nothing rests on it.
- Extra control (16:28 entry): a SECOND level-3 lift in the liftable coset, xi'' = -17794, gives Theta_4 = -1 by
  the repo pipeline and by my own local terms (two twists). Not in the repository.
VERDICT item 2: HOLDS (tier C: Fisher through summaries).

## 16:31 UTC - item 4: the dimensions (2 - 1) + (3 - 1) - 2
- DIMS['945,248'] = (2, 3) was already in gen.py before t9 (the diff does not touch it). Its source:
  rank_fisher.txt, 945,248 block: "(1) first descent: dim S'_1 <= 3 (E side, classes d | b), dim S_1 <= 4 (E' side);
  second descent: dim S'_2 <= 3, dim S_2 <= 2; known: 2, 1; open classes y1 = 2 (E side), y0 = -5820647 (E' side)"
  and "(5) dim S_2 = dim S_3 = dim S_4 = 2, dim S'_2 = dim S'_3 = dim S'_4 = 3".
- Bookkeeping for --side=E': lifted = sides["E'"].dd = -5820647 (S_4, E' side), y = sides['E'].dd = 2 (S'_4);
  working curve E (aw, bw = a, b), K2 classes sqfree(b) = 678244410280322944874881 = 13*349*1493911429*100067394997
  and sqfree(x(G)) = 2601752269922 = 2*13*100067394997; DIMS read as (dim S_4 (E' side), dim S'_4 (E side)) =
  (2, 3); printed kernels (1, 2). Right.
- My own PARI checks: G on E, ellorder(E, G) = 0; torsion Z/2 on E and E'; core(a^2 - 4b) = -1 (class of T'),
  core(x(phi(G))) = 1, so the image of E'(Q) in S_1 is <-1> (dim 1 = "known: 1"); the image of E(Q) in S'_1 is
  <b, x(G)> (dim 2); b, x(G), 2 independent mod squares; ellrank(E) = [1, 3]. y = 2 is in x(E(Q_p)) mod squares at
  every bad prime (own search; at 2 by x = 162, x(x^2+ax+b) = 2^12 u, u = 1 mod 8) and at R.
- The bound: S_4 = <-1, y0> (upper bound S_4 in S_2, dim <= 2), S'_4 = <b, x(G), y> (upper bound dim S'_1 <= 3);
  the images of the rational points lie in all S_m, S'_m and in the kernels; Theta_4(y0, y) != 0 gives S_5 = <-1>,
  S'_5 = <b, x(G)>. Lemma 2.2 (summary): the 2^n-descent bound is dim S_(2n-1) + dim S'_(2n) - dim E(Q)[phi]
  - dim E'(Q)[phihat]; for n = 3 this is dim S_5 + dim S'_6 - 2 <= 1 + 2 - 2 = 1 (S'_6 in S'_5; the other order
  likewise). The printed "Lemma 2.2 (n = 3) gives rank <= 1 + 2 - 2" is that bound with S'_6 replaced by S'_5 (an
  inequality in the right direction). y in S'_4 rests on rank_fisher.py (y in S'_3 by checked level-3 lifts; Theta'_3
  alternating on S'_3 with b, x(G) in its kernel, so y in S'_4) - inherited, not re-audited (gap G3).
VERDICT item 4: HOLDS.

## 16:32 UTC - item 5 (cont.): could the edits change existing outputs?
Read the diff: (a) load_hints returns (key, xi2); its only caller is prepare() (grep: no other script calls it);
xi2 is not None only for the key "945,248 E+|E'" + "|xi" (the only such entry). (b) The seed retry in choose_Xi_lin
runs only where seed 1 sampled no point; the old code returned [] there (no twist; main asserts, so no output was
ever written in that case), so every existing output is untouched; _control, _recheck, _weil call choose_Xi_lin and
their diffs are empty. (c) The print and the prepare() assertion are guarded by xi2. (d) DIMS unchanged.
Runs (16:12 entry): all 13 existing lines empty diffs; 945,248 main 43 s and --negative 44 s reproduce the committed
outputs; norm.py default 1 s, --find 24 s (both hints equal), --obstruction 130 s (byte for byte). HOLDS.
NIT: the gen docstring lines 15-16 ("for (631,71) by --side=E', for the other curves because that is rank_fisher.py's
side") should now read "for (631,71) and (945,248) by --side=E', for the other curves ...".

## 16:46 UTC - more checks
- Weil-reciprocity test (the repo's test function, run here on 945,248 E', Xi = 521375088721, seed 3, 15 pairs) at
  two primes t9 did not use: p = 30011: 14 of 14 ratios squares, odd-divisor control g0 c/e 6 of 14; p = 1000003:
  14 of 14, control 8 of 14 (/tmp/z10/work/weil945.out). A test, not a proof (gap G2).
- Numbers in t9's report checked against the committed outputs: Xi = 521375088721 with -1 at 2, 29, 349,
  32746027, 45237653, 1493911429, 100067394997; K3 Xi = 1300876134961 with -1 at 13, 61, 349, 32746027, 45237653;
  K2 +1, +1; K5 -1; hint digits U [54, 51, 96, 95, 50, 46], V [56, 52, 94, 91, 40, 35]; class groups [4, 2, 2] and
  [604, 2, 2] (h = 16, 2416); --negative 9 of 9; obstruction table. All match. (The scratch values Xi = -13 etc. in
  t9's log are from uncommitted points; not checkable, not used.)
- NITs in rank_fisher_theta4_norm.py's docstring: "Hints: U point 35-45 digits per coordinate (Z, g), V point 36-56;
  e, c about 90" describes the scratch set; the committed U hint has Z 54/51, g 50/46, e 96/95 digits. "(5-10 steps
  here)": --find logs 9 and 19 steps (U, seeds 1-2) and 12 (V). The claim "the json's C3 (xi = 1271) has an empty
  level-4 Selmer set" is printed as a fact by gen.py and the docstrings; it rests on the sampled obstruction of
  --obstruction (PARI p-adic points, not Hensel-certified; the local-duality argument makes the sampled symbol the
  true one once the points are genuine). Nothing rests on it; wording should say "by a sampled dual obstruction".
- Not covered by committed outputs or CI for (945,248): _recheck (Hensel re-derivation), _order (G of infinite
  order), _weil (fixed curve lists in those scripts). t9 ran them in scratch; I re-derived the terms independently
  (item 1), checked ellorder(E, G) = 0 with PARI, and ran the Weil test at two more primes. Gap (new, G6 below):
  the committed CI does not re-derive the 945,248 terms by Hensel; suggest adding 945,248@E' to the three scripts'
  lists in a later run (builder work; not done here).

## 16:48 UTC - item 6: status words and replacement text
t9's report and outputs stay at tier C ("tier C, as for the other Theta_4 curves; not released; NOT REVIEWED");
the gen output's last line ("rank = 1 ..., without Magma (written arguments as listed in the docstring)") is the
wording of the other three curves. t9's proposed LEDGER/README/CHANGES text is in its final message, which I do not
have; below is exact replacement text against int8 (= run 7 + run7_t9; LEDGER line 140 = t9's 134; README lines
100, 102, 119 as in t9). Use it in place of t9's proposal where they differ.

LEDGER row 53, third column. Replace
  "for the curves of (760,39) and (920,129) also from `Θ₄` since 27 November 2026 (after release 0.6.0, not
  released), `rank_fisher_theta4_gen.py`, tier C, reviewed with gaps, so of this row only (945,248) rests on Magma
  alone;"
by
  "for the curves of (760,39) and (920,129) also from `Θ₄` since 27 November 2026 (after release 0.6.0, not
  released), `rank_fisher_theta4_gen.py`, tier C, reviewed with gaps, and for the curve of (945,248) also from `Θ₄`
  since run 7 (27 November 2026, not released; on its `E′` side; tier C, reviewed with gaps,
  `reviews/review_z10.md`), so no curve of this row rests on Magma alone;"
and replace
  "(from PARI and Magma for four of them; from Magma alone for the curve of (945,248); for the curves of (72,65),
  (760,39) and (920,129) from Magma and from Fisher's `Θ₄`, tier C (the last two since 27 November 2026, not
  released))"
by
  "(from PARI and Magma for four of them; for the curves of (72,65), (760,39), (920,129) and (945,248) from Magma
  and from Fisher's `Θ₄`, tier C (the last three since 27 November 2026, not released))".

LEDGER row 65, first column. Replace "one counted ratio, (945,248), now rests on Magma alone for an upper rank bound)"
by "one counted ratio, (945,248), then rested on Magma alone for an upper rank bound, until run 7: next row)".

LEDGER, new row after row 65 (three cells):
| **`Θ₄ ≠ 0` on `E⁺` of (945,248): rank exactly 1 without Magma** (since 27 November 2026, run 7: not released; the curve of the counted closure (945,248), which keeps its status; with it no counted ratio rests on Magma alone for an upper rank bound) | `scripts/g1_full/rank_fisher_theta4_gen.py 945,248 --side=E'` (the general module, on the `E′` side as for (631,71); the level-3 lift ξ = 17794 = 1271·14 over the same `D_λ`, λ = −349, because the lift of `rank_fisher_curves.json`, ξ = 1271, has no level-4 lift by a sampled dual obstruction, `rank_fisher_theta4_norm.py --obstruction`; Fisher, Sect. 3: any level-4 lift may be paired); conic points over ℚ(√5820647) and ℚ(√−5820647) from `rank_fisher_theta4_norm.py` (a finder that uses PARI under GRH: class-group reduction, a Legendre-type descent, `rnfisnorm`, a quaternion reduction of the point; every point re-checked exactly, so nothing rests on GRH); `--negative`; CI job `descents` | `Θ₄(−5820647, 2)`: the product of the local symbols is −1; K2 (the classes of `b` and `x(G)` pair to 0), K3 (a second twist) and K5 (other norm elements) agree; the paired class 2 is not 1 mod 8, so the 2-adic symbol is exercised. Reviewed by a separate AI agent (`reviews/review_z10.md`): holds with gaps; recomputed with its own Hensel-certified local points, Hilbert symbols and place set for six level-4 twists (five not found by the repository's sampler) and seven choices of norm elements, and through a second level-3 lift (ξ = −17794): the same value. Gaps: those of the row above (Fisher's paper read through tool summaries, here also for the independence of the lift; the evenness of `div g0` and the good-reduction step are written arguments, tested by Weil reciprocity modulo four primes; the older scripts not re-audited, among them `rank_fisher.py`'s `y ∈ S′₄`; no known-answer control in another field); `_recheck`, `_order` and `_weil` do not include this curve; the "forged twist" of its negative test is a genuine twist, refused only because the sampler finds no point at 97, so that refusal tests nothing |

LEDGER line 140 (t9: 134). Replace the sentence beginning "The one curve whose rank 1 rests on Magma alone for a
counted ratio that no other method closes, `E⁺` of (945,248)" and ending "(the exact second descent of
`rank_indep.py` stops at 3)." by
  "No curve is left whose rank 1 rests on Magma alone for a counted ratio that no other method closes: the last,
  `E⁺` of (945,248), has `Θ₄` since run 7 (27 November 2026, not released; on its `E′` side, conic points over
  `ℚ(√±5820647)` from `rank_fisher_theta4_norm.py`; tier C). For these curves (those of (72,65), (760,39),
  (920,129), (945,248) and the `D1` curve of (631,71)) the bound 1 without Magma rests on Fisher's paper read
  through tool summaries and on written steps (the evenness of `div g0`, the good-reduction argument). The `D1`
  curve of (1389,59), which closes (724,665) in the 122, still rests on Magma alone, but (724,665) is closed again by
  the four-point argument."

README line 100. Replace "for the whole census, one counted ratio, (945,248), rests on Magma alone for an upper rank
bound (eight until 26 November 2026, six until 27 November 2026 and in release 0.6.0; the change is not released)"
by "for the whole census, no counted ratio rests on Magma alone for an upper rank bound since run 7 (eight until 26
November 2026, six in release 0.6.0, then one, (945,248); not released; each of these bounds now also comes from
Fisher's `Θ₄`, tier C)".

README line 102. Replace "so only (945,248) does;" by "so only (945,248) did, until run 7 (27 November 2026, not
released), when its curve got `Θ₄` as well (on the `E′` side; tier C), so none does;".

README line 119. Replace "so **one counted ratio, (945,248), rests on Magma alone** (tier C, with the gaps of
(72,65); reviewed by a separate AI agent: holds with gaps), and in the census" by "so one counted ratio, (945,248),
rested on Magma alone (tier C, with the gaps of (72,65); reviewed by a separate AI agent: holds with gaps); since
run 7 (27 November 2026, not released) `Θ₄ ≠ 0` on `E⁺` of (945,248) too (on its `E′` side, with conic points over
`ℚ(√±5820647)` from `rank_fisher_theta4_norm.py` and the level-3 lift ξ = 17794; tier C; reviewed by a separate AI
agent: holds with gaps), so **no counted ratio rests on Magma alone for an upper rank bound**, and in the census";
and replace "on Magma alone for the curve of (945,248), and on Magma and Fisher's `Θ₄` for the curves of (72,65),
(760,39) and (920,129) (the last two since 27 November 2026)." by "on Magma and Fisher's `Θ₄` for the curves of
(72,65), (760,39), (920,129) and (945,248) (the last three since 27 November 2026, not released)."

CHANGES (under the unreleased run that takes run7_t9):
    - *Fisher's `Θ₄` on `E⁺` of (945,248)* (`scripts/g1_full/rank_fisher_theta4_gen.py 945,248 --side=E'`; conic points over ℚ(√±5820647) from the new finder `rank_fisher_theta4_norm.py`: class-group reduction, a Legendre-type descent over the field, PARI's `rnfisnorm` under GRH and a quaternion reduction of the point; a finder only, every point re-checked exactly; tier C; reviewed by a separate AI agent, `reviews/review_z10.md`: holds with gaps, recomputed with the reviewer's own local points, Hilbert symbols and place set for six level-4 twists and seven choices of norm elements, and through a second level-3 lift): `Θ₄(−5820647, 2) ≠ 0`, so rank ≤ (2 − 1) + (3 − 1) − 2 = 1 without Magma; with it no counted ratio rests on Magma alone for an upper rank bound. The level-3 lift of `rank_fisher_curves.json` (ξ = 1271) has no level-4 lift (a sampled dual obstruction, `--obstruction`); its twist ξ = 17794 = 1271·14 over the same `D_λ` is used (Fisher, Sect. 3: any level-4 lift may be paired). `rank_fisher_theta4_gen.py` gains an optional ξ hint and a seed retry where the sampler finds no point; no existing output changes. No census count changes.

CI (job `descents`, in the Θ₄ block, about 2 minutes):
          python3 rank_fisher_theta4_gen.py 945,248 "--side=E'" 2>/dev/null | grep -v 'PARI stack' | diff - rank_fisher_theta4_gen_945_248_Eprime.txt
          python3 rank_fisher_theta4_gen.py 945,248 "--side=E'" --negative 2>/dev/null | grep -v 'PARI stack' | diff - rank_fisher_theta4_gen_945_248_Eprime_negative.txt
          python3 rank_fisher_theta4_norm.py 2>/dev/null | grep -v 'PARI stack' | diff - rank_fisher_theta4_norm.txt
          test "$(python3 rank_fisher_theta4_norm.py --find 2>/dev/null | grep -c 'equal to the committed hint: True')" = 2
(the first three ran here from a clean clone with empty diffs; the --find line prints timings, so it is grepped, not
diffed; it relies on PARI 2.17.4's bnfinit/rnfisnorm and a 60 s alarm, so on another PARI it may find other points
without any effect on the result). MANIFEST.sha256: gen.py and the hints file changed; five new files.
Docstring NIT for gen.py lines 15-16: "for (631,71) and (945,248) by --side=E'".

## 16:55 UTC - the json's lift xi = 1271 does not lift: own certified check
/tmp/z10/work/obstr_own.py (output obstr_own.out, 53 s): f0 of C3'(xi) from a U conic point (t9's finder for
xi = 1271, run here, nothing written; the committed hint for 17794), my own Hensel-certified points of C3'(Q_p)
(2x2 Jacobian minors) at all primes < 60 and B0 (26 places), own symbols (= PARI hilbert), the real signs of f0 on
C3'(R) by interval arithmetic:
  xi = 1271:  Ob(-1) = +1, Ob(14) = +1, Ob(17) = -1 (term -1 at 17), Ob(97826) = -1 (at 41), Ob(-5820647) = +1
  xi = 17794: Ob = +1 for all five eta.
Exactly t9's table. With Poitou-Tate (a global Xi with Xi in the local coset at every place exists iff the
product of symbols vanishes on S_1; written, standard), the json's C3 has an empty level-4 Selmer set, and this no
longer rests on uncertified points. (Remark: the repo's --obstruction prints "sign of f0 on C3'(R) -1" for
xi = 17794 while gen.py has f0 > 0 there; my run with the committed hint gives +1. The --obstruction run finds its
own U point after the xi = 1271 search (another PARI random state), so its f0 differs from gen's by a constant;
Ob is unchanged by a constant (product formula). Not an error; worth one sentence in the docstring.)
So the wording NIT of the 16:46 entry becomes: "(a dual obstruction from sampled local points, confirmed with
certified points by the review)".

## Final (16:58 UTC)
**Verdict: HOLDS WITH GAPS.** Theta_4(-5820647, 2) != 0 on E+ of (945,248), E' side, so rank <= (2-1)+(3-1)-2 = 1
without Magma, at tier C (not released). I found no error that changes the value.
Items: 1 HOLDS (own Hensel-certified points, own symbols = PARI, own place set = the repo's S; -1 for six level-4
twists, five not found by the repo, seven choices of norm elements, and through a second level-3 lift xi = -17794);
2 HOLDS (Fisher Sect. 3, via summary: any lift in S^(4)(E'/Q) may be paired, y in S'_4 needed and inherited; C4
locally soluble at all of S and R by the repo and by me; the json's xi = 1271 does not lift: confirmed with certified
points); 3 HOLDS (analyse is exact and PARI-free; own exact check of both points, tangency, f, f_dag, and an explicit
ideal certificate; --find reproduces the hints; nothing rests on GRH); 4 HOLDS (DIMS from rank_fisher.txt; side
bookkeeping right; own PARI checks of G, torsion, classes); 5 HOLDS (13 existing CI lines, 945,248 main/--negative,
norm.py default/--find/--obstruction: all reproduce from a clean clone, PARI 2.17.4); 6 HOLDS WITH GAPS (text).
Gaps: G1-G4 of review_y5 (Fisher through summaries, now also for the independence of the lift; written evenness and
good-reduction steps, Weil test at 10007, 100003, 30011, 1000003; inherited scripts incl. y in S'_4; no known-answer
control in another field). New: N1 the --negative "forged twist" for this curve is a genuine twist and is refused
only because the sampler finds no point at 97 (it refuses the true Xi too): 8 of 9 refusals mean something;
N2 _recheck/_order/_weil do not include this curve; N3 the sampler sees 2 of 8 twists (one unit class at 97,
245411); N4 docstring nits. FIX-FIRST: none. FIX: the doc text above (16:48 entry). Referee scripts: /tmp/z10/work/
(indep.py, strat.py, check_hints.py, obstr_own.py, forged.py, curve_checks.py, y2_at_2.py, weil945.py), /tmp/z10/ci/,
scratch copy /tmp/z10/k4 (xi = -17794). No repository file edited.
