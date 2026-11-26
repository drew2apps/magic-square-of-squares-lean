# Review of v4 (branch r3_v4, e0df800 vs f03f8af): rank_second.py fix + rank_fisher.py

Hostile referee, written by an AI model (Claude), 25 Nov 2026, 13:50-14:25 PT. I did not edit the worktree. Scratch
files are in /tmp/rev_v4/: the copies wt/ and wt2/, run_blocked.py, and in indep/ check_c3.py, check_d2.py,
pts_from_author.py, mutA.json and mutB.json.

## (1) FIX-FIRST of review_t2: the bad primes of D_lam. APPLIED.
- `rank_second.py:239-267`: `DLam(..., bad=...)` and `S(dp) = {2} | bad | primes(lam, d', delta, det)`. Both main
  call sites (:462, :503) and every `DLam` call in rank_fisher.py pass `bad=info['primes']`. That list is
  `rank_indep.second_descent`'s `ps` (:497-503): 2, b', N0, disc U, disc V and the pencil discriminant.
  `g_class_local` asserts NPTS exact local points at every p in S, and `g_sign_real` asserts a real point. The
  sites without `bad` are dead code (:551) and negative tests (:972-1005), so they are harmless.
- I checked that this bad set is complete. The pencil is a quartic of the form x y Dp(x, y). Its discriminant
  involves Res(x, Dp) = disc V, Res(y, Dp) = disc U and disc Dp. The primes of lam lie inside ps.
- Code diff against t2/: only those lines change. Output diff against t2's rank_second.txt: after deleting the
  `|S| = n` fields the diff is EMPTY. The same lifts are evaluated or skipped, the values are the same, and so are
  the places with factor -1.
- Independent data check (check_d2.py: my own square test and my own bad-prime computation). All 10 level-2 curves
  in rank_fisher_curves.json have an exact local point at every bad prime and a real point. My search missed 3
  small-disc cases (72,65 E' at 710053; 920,129 E' at 4883483573; D1 E at 2). There I verified the author's
  `local_points` output with my own test, and it passes. Every D used in the Theta_4 reductions was found by my
  own search.

## (2) Theory, from Fisher's paper
WebFetch could read the PDF arxiv.org/pdf/1509.03234v1 (the author could not). I read it through tool extraction,
not page images. These statements were quoted from the paper:
- Thm 2.1: for odd m, the pairings Theta_m, Theta'_m are "alternating" with kernels S_{m+1}, S'_{m+1}; for even m,
  S_m x S'_m has left and right kernels S_{m+1}, S'_{m+1}.
- Lemma 2.2: rank <= dim S_{2n-1} + dim S'_{2n} - dim E(K)[phi] - dim E'(K)[phihat].
- Sect. 4: Theta_m(a_1, b_1) = <a_m, b_1>, for any lift a_m.
- Thm 3.1: (i) <psi x, y> = <x, psihat y>; (ii) x is in im(psi) iff x is orthogonal to ker(psihat).

The dimension count re-derived. Theta_2 = 0 is computed, so S_3 = S_2 (dim 2) and S'_3 = S'_2 (dim 3); these are
exact because the y-lifts are ELS. An alternating form whose kernel contains a subspace of codimension <= 1 is
zero, so S_4 = S_3 and S'_4 = S'_3, and n = 2 gives 2 + 3 - 2 = 3. HOLDS. The author's own argument also holds:
Theta_3(a, a) = 2<a_3, a_3> = 0 by antisymmetry alone.

Theta_4 reduction re-derived. With S_1 = S_2 on the C3 side, a_3 is orthogonal to Sha(E)[phi] = <points, y0>,
because 2<a_3, a_3> = 0. By 3.1(ii), a_3 = phihat(a_4), with the Selmer-level correction coming from rational
points. Then <[C3], [D]>_E = <a_4, phi[D]>_E' = Theta_4(y0, y1). HOLDS on all five curves (and with the sides
printed).

## (3) Reproduction
- I ran both scripts, main and --negative, with `sys.modules['cypari2'] = None`. All four outputs are BYTE-IDENTICAL
  to the committed .txt files (negatives: rank_second 7/7, rank_fisher 8/8). A `--write` run in a scratch copy
  reproduced rank_fisher.txt, rank_fisher_curves.json and rank_fisher_hints.json byte for byte.
- Without --write, sha256 of every file in the copy was unchanged (I ran with -B).
- No PARI in the certified path. sympy is still used, in pre-existing rank_indep code: ground_roots, the resultant
  and the pencil expansion. I recomputed the pencil discriminants myself (see (1) and (4)).

## (4) Local points: all 178 rechecked
check_c3.py is my own code: PARI p-adics to O(p^150) and exact sign tests at R. It confirms three things on all 10
level-3 curves:
- p_X and p_Y are on W^2 = lam U and W^2 = lam V, and l_X, l_Y are exactly the tangent lines there.
- All 178 points give xi f(P) a nonzero square for some sign choice: 178/178.
- Each place set contains every prime it must: 2-13, lam, xi, disc U, disc V and the pencil discriminant. So the
  argument for p outside S3 (Hasse bound p+1-2sqrt p > 8 for p >= 17, an etale connected cover, Hensel) applies.
  The content primes of l_X and l_Y lie in {2, lam, disc U}, which are in S3.

The checker has teeth: with xi -> 3xi, 57 local points fail on all 10 curves; with l_X[2] + 1, tangency fails and
69 fail.

## Verdict
FIX-FIRST (text only; no number changes)
1. The equation of C3 and D is printed wrong. rank_fisher.py:40, 56, 413, 584 (and so rank_fisher.txt:42, 69, 88,
   111, 135) and rank_fisher_note.md:38, 60 say `lam X^2 = U, lam Y^2 = V`. The code, p_X and l_X all use
   X^2 = lam U, Y^2 = lam V, as in rank_fisher.py:34. Under the printed reading, p_X is not on the curve for every
   lam other than -1 (18 of the 20 stored curves), so the "explicit curves" for Theta_4 are mislabelled. Fix the
   label, and put the equations into rank_fisher_curves.json, which currently has none.
2. rank_fisher_note.md:97-105, rank_fisher.py:73-77 and the v4 report give the wrong reason why Theta_4 is not
   done. They quote "the limit is C_4, which fails hypothesis (19) ... degree-4 conic" as Fisher's words; that is
   not his text. Fisher Sect. 7 says that C_3 satisfies (19): one builds a phi-covering of C_3 (that is, C_4) and
   then computes Theta_4, solving conics over K and over quadratic extensions. Sect. 9's limit is the next step:
   C_4 fails (19). Theta_4 is Fisher's fifth descent, which works for any rational 2-torsion point, and it is
   presumably how Magma's TwoPowerIsogenyDescentRankBound reached 1. State it as "not implemented (it needs conics
   over a quadratic field)", not as blocked.

NIT
3. "Theta_4 (the 8-descent pairing)" (rank_fisher.py:66, 574, 712): Fisher says the sixth descent corresponds to
   8-descent, and Theta_4 is the fifth descent. Say "Theta_4 (fifth descent); Theta_4(y0, y1) != 0 gives the bound
   1". Also, "or deeper" is moot: with one rational 2-torsion point, Magma's method stops at S_5.
4. review_t2's NITs 3, 4 and 6 are still open: dead code at rank_second.py:535-616, no b1<->b2 negative test for
   E0L, and :777 `e1 * 0 + e3`.
5. The report gives 313 s for rank_fisher; my hinted runs took 55 s CPU. The 313 s presumably included finding
   the hints.

HOLDS
- The FIX-FIRST is applied, and every lift used is ELS at every bad prime.
- Theta_2 = 0; Theta_3 and Theta'_3 are alternating, hence 0; the 4-descent bound is exactly 3 on all five curves.
- The level-3 lifts are constructive, and the 178 points are rechecked.
- Theta_4 is reduced correctly to <[C3], [D]>.
- E0L(920,129) has rank exactly 3: unchanged output, and review_t2's check stands.
- Nothing is written without --write, and there is no PARI in the certified path.

Doc sentence:
> On E+ of (760,39), (945,248), (72,65), (920,129) and D1 of (631,71), exact arithmetic (rank_second.py,
> rank_fisher.py; AI-written, independently re-checked) shows that descent by 2-isogenies through the 4-descent
> gives exactly 3 (Theta_2 = 0 is computed; Theta_3 and Theta'_3 are alternating and so vanish). So Magma's bound 1
> needs Theta_4(y0, y1) != 0 (Fisher's fifth descent, arXiv:1509.03234 Sect. 7), which is not reproduced.
