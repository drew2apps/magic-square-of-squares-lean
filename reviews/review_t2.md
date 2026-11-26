# Review of sess_t2 (rank_second.py), hostile referee, AI-written (not read by a person)

## Log
- 19:25 UTC start. Read t2.md (discounted), rank_second.py in full.
- (c) Ran `python3 -I rank_second.py` in a copy (/tmp/rev_t2/run): 2m46s, exit 0, output byte-identical to
  rank_second.txt. `--negative`: exit 0, byte-identical to rank_second_negative.txt (7/7 refused). md5 of every
  file in the copy before/after: only `__pycache__/*.pyc` appear (Python bytecode, not script output); no
  .txt/.json written without --write. --write paths: rank_second.py:1025, :1030, :1045 all guarded by WRITE.
- (d) PARI in the code: factor (rank_indep._find_factor; result re-multiplied and every prime proved by trial
  division / Pocklington, rank_indep.factor_int), qfsolve for conic points (rank_second.py:157, :518, :713; every
  point checked on its conic exactly, :178, :531, :719). No ellrank output is used in the certified path. Numpy
  floats only place real test points (:387); signs are then exact.
- (c) Also reran with `sys.modules['cypari2'] = None` (PARI blocked): exit 0, output byte-identical. So every
  PARI result the run needs is in the hints files, and each one is checked on use.
- (b) dim Sel^2 = 7 for E0L(920,129), y^2 = x(x - 1095776342882)(x - 1883582727842): cross-checked with PARI 2.17.4
  `ell2cover` (7 ELS 2-covers, a separate algorithm) and `ellrank` = [1, 3, 2]. Same for E0L(72,65): 7 covers,
  [1, 5, 0]. The code review of full_selmer2 (:638-706) found no fault. The local image orders 4 / 8 / 2 at odd p / 2 / R
  are right. The images are spans of images of explicit local points (exact square tests). The Kummer map,
  including the T1 and T2 special cases (:649-652), is right. The reduction `red` is linear because the basis is in
  reduced echelon form. The kernel bookkeeping gives an independent basis.
- (a) Formula re-derived. Take the 2-covering x = e_i + a_i z_i^2. The sign changes of one coordinate are
  reflections. Since z_i maps to e_2(T_i, T) z_i under translation by T, the product (u,v,w) -> (u,-v,-w) is
  translation by T1, and (u,v,w) -> (-u,v,-w) is translation by T2. So F_v - F_w ~ T1 and F_u - F_w ~ T2. The
  tangent line l_v vanishes to order 2 at p_v on C_v, so div l_v = 2F_v on the degree-4 model and div(l_v/l_w) = 2(F_v - F_w).
  The Kummer coordinates (b1, b2) = (e(xi,T1), e(xi,T2)) give xi = b2 (x) T1 + b1 (x) T2. With the Brauer
  element (c, g), div g = 2 Delta, [Delta] = T, which represents c (x) T in H^1(Q, Pic^0 D) = H^1(Q, E)
  (Skorobogatov, Torsors and Rational Points, Thm 6.2.3; Cassels 1998, sec. 7, in the 2-isogeny case), this gives
  <xi, D_a> = prod_v (b2, l_v/l_w)(b1, l_u/l_w)(P_v), which is :837. Places outside S: if
  p does not divide 2 a1 a2 disc b, the conic reductions are smooth and the primitive tangent line reduces to the
  tangent line, so there is a Z_p-point off the three fibres with unit l's. The local factor there is 1, and the
  delta/det primes in S are harmless extras. Hilbert symbols (:89-106) match the textbook formulas (Serre, Cours
  d'arithmetique III.1.2). The real-sign routines _sign2 and _sign3 (:848-869) were checked case by case.
- (a) INDEPENDENT re-implementation (/tmp/rev_t2/indep/ctp.py). It uses PARI p-adics to precision p^80 and
  pari.hilbert. Each conic point is moved by 1 and by 3 random lines, so the tangent points differ from the
  author's. The local points are its own: 4 and 10 x-values per place, all 8 sign choices. The places include
  extra 3, 5, 7, 11 and 13. Result: the 7x7 matrix for E0L(920,129) is IDENTICAL to rank_second.txt:213-219 for
  two seeds. The matrix for E0L(72,65) is all +1 for both seeds. The evaluation was constant at every place
  tested. At 9 good primes (17 to 59), for 3 D_a x 7 b, every factor is +1 and constant. A class outside Sel^2
  ((3,1), (1,3), (7,1), (-1,1), (2,1)) gives a non-constant evaluation, at 2 or at R.
- (b) rank(M) = 2: rows 2, 3 and rows 6, 7 are equal, and the submatrix on rows (2, 6) x cols (6, 2) is
  [[-1, +1], [+1, -1]] in +-1 form, so it is invertible over F_2. So 7 - 2 - 2 = 3. The point images span 5
  (re-checked below).
- (b) Point images: my own Kummer map on the three G_T points of dm4_log.jsonl and on T1, T2 gives F_2-rank 5
  (/tmp/rev_t2/indep/pts.py). The points are on the curve, and PARI's height regulator is about 2970.5, nonzero.
  So 3 <= rank <= 3. HOLDS.
- (a) Mutation: swapping the roles of b1 and b2 in the formula makes the evaluation non-constant on all 49 entries.
  So the constancy check has teeth against a T1/T2 convention error.
- Claim (1), the isogeny side. Re-derived: the reflections (X,-Y) and (-X,Y) have product the translation by T, so
  F_Y - F_X ~ T and (d', l_Y/l_X) is the Brauer element of d' (x) T. The control (945,248) is consistent with the
  theory. The image of psi_* is the orthogonal complement of Sha(E)[phi], so the excluded d' must pair -1 with
  some iota(lambda); 6/6 do, and 2/2 kept classes are +1. The logic "bound 1 => <y0,y1> != 0" is sound.
  GAP: DLam.S (rank_second.py:259-263) has only 2, lam, d', delta and det. The paired lifts are not certified
  locally soluble at the other bad primes of D_lam (primes of b'N0, disc U, disc V and the pencil discriminant;
  rank_indep.py:504-507). rank_indep keeps "undetermined" lambdas (rank_indep.py:516, `if r is not False`). Without an
  ELS lift, x would not be in phi_*(Sel^2) and the Sel^2 bound would be < 3. Patch test
  (/tmp/rev_t2/indep/els_lifts.py): with those primes added to S, all 138 evaluated lifts get exact local points
  there, and the output is otherwise identical. With NPTS = 8 and another seed (els_more.py), every evaluated
  value is still +1. Some lifts are skipped only because the search budget runs out at large p.
- Doc text in t2.md: "Magma's bound needs descents by isogenies of degree 8 or more" is not justified by anything
  computed. It should say "needs information beyond Sel^2 and its Cassels-Tate pairing (e.g. Fisher's higher
  descents, Math. Comp. 86 (2017))".

## Verdict
FIX-FIRST
1. rank_second.py:259-263: add rank_indep's bad primes of D_lam to DLam.S, or assert RI.solvable(...) is True there.
   Until then the "computed exactly" sentence for claim (1) is not self-certified. (Reviewer's patch: passes, no
   change in results.)
NIT
2. t2.md proposed G1_FULL_SYSTEM text: drop "isogenies of degree 8 or more".
3. rank_second.py:535-616: dead code (lift_pairing, analyse_general, find_conic_point_any, e0l_curve,
   RATIONAL_LIFTS), never run by main.
4. rank_second.py:953-1010: negative f3 is refused trivially (point off the conic). No negative test targets the
   E0L matrix; add the b1<->b2 swap.
5. rank_second.py:48: "every lift with local points everywhere" should read "local points at every place of S",
   until fix 1 is in.
6. rank_second.py:773: `e1 * 0 + e3`, cosmetic.
HOLDS: (a) formula, Hilbert symbols, real signs, places outside S, symmetric/alternating, point images in kernel;
an independent re-implementation reproduces the 7x7 matrix. (b) dim Sel^2 = 7 (agrees with PARI ell2cover),
pairing rank 2, images span 5: rank E0L(920,129) = 3. (c) Both outputs reproduce byte-identically, also with PARI
blocked; nothing is written without --write. (d) PARI supplies only hints (factorisations proved prime, conic
points and curve points, each checked).
