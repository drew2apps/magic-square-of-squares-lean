# z11 — hostile referee of t10 (commit 714b45c, branch run9_t10, experiments/conductor-t10/)

Start: Fri Nov 27 19:22 UTC 2026. Hard stop 22:30 UTC. Written by an AI model; not read by a person.
Scratch: /tmp/z11/. I edit no repository file.

## ~19:30 UTC — read
- Read BRIEF_COMMON7.md, t10's NOTES.md, t10_summary.txt, reports/t10.md, t10_pairct.py (733 lines), z7's my2desc.py and
  ct_mine.py, LAST_PAIR.md §9–12, tuple_descent.magic. Forms of t10 (a0 = m + int, a8 = n + imt, aL, aU) checked by hand
  against LAST_PAIR §11 γ's: γ0 = v a0, γ8 = v conj(a8), γL = v aL, γU = v aU; so class(a_j a_k(t)) = d_j d_k at a
  solution (the factor v cancels in the product), and F_j(t) = N(a_j) is a square there. TD.magic's tuple order is
  (F0, F8, GL, GU) = (0, 8, L, U), as t10's IDX assumes.
- Base points checked by hand: F0 F_L(m/n) = 2m^2 * 2(m^2-n^2)^2/n^2; F0 F_U(-m/n), F8 F_L(n/m), F8 F_U(-n/m) likewise;
  F0 F8(0) = m^2 n^2; F_L F_U(0) = 4(m^2+n^2)^2. All non-zero squares for every m > n > 0: t0 is rational for every ratio.
- Plan: own script /tmp/z11/z11_pair.py: own forms, own Washington map (verified symbolically), root identification by
  a finite-field test independent of t10's resolvent formula, c(tau) by the divisor-evaluation form of the Kummer map
  (no constants read off points), Sel^2 and pairing by z7's my2desc.Curve / ct_mine.Cov (reused, read; primes from a
  PARI factorisation instead of sympy's), residual tuples by G1's exact_local (as G2/t10).

## ~19:38 UTC — independent recomputation, first results
Own script /tmp/z11/z11_pair.py (own forms, base points, Washington's map — the identity and its inverse verified
symbolically by /tmp/z11/wash_check.py, residue 0 — root identification by a finite-field test at p = 101, 109, 113:
chi(X - e_i) = chi(q_P(t) q_P(t0)) at every F_p point, for exactly one factorisation P per root; class of a tuple by
the divisor form of the Kummer map, delta_r = F_j(t)/F_j(t0), delta_* = +-a_j a_k(t)/a_j a_k(t0) mod Q(i)*^2, so
c(tau) = (F_j(t0), +-d_j d_k D(t0), product) with no constant read off points; Sel^2 and pairing by z7's
my2desc.Curve / ct_mine.Cov, reused after reading, primes by PARI factor (isprime checked); BOTH signs tested
separately; every PARI point found is checked against the divisor formula, in my Sel^2, and pairs +1 with all of Sel^2).
- (465,416) (in g2_997_all.txt): 6 tuples x 6 curves. My classes equal t10's printed classes on every line (t10's
  integral models differ from mine on E_0L..E_8U, the classes agree), the same killing b's ((1,65); (1,7), (2,195),
  (168240919681,1); (193,7)), the same verdicts: all 6 pairing-killed, LU kills exactly (1,2,1,1), (1,2,2,2). Both
  signs give the same verdict everywhere. The places of the local -1's differ from t10's (e.g. [3] vs [13]): the
  local decomposition depends on the local points chosen; only the product is intrinsic.
- (5623,1017): both tuples NOT pairing-killed on any of the six curves, both signs; E_08 dim Sel^2 = 2 and the class is
  (1,1) = delta(O) (sign +) or (1,-1) (sign -), both torsion classes. Matches t10.
- Root roles from my F_p test equal t10's resolvent roles (e_r = e1, e* = e3 on 08 and LU; e_r = e1, e* = e2 on the
  mixed curves). delta(psi(iota t0)) = (1, -1, -1) in (r, *, ') on the mixed curves; on 08 and LU psi(iota t0) is a
  2-torsion point with * and ' coordinates -1: the sign ambiguity is the class of a rational point, as t10 says.

### Counts (item 5), own parser /tmp/z11/recount.py
- 201 ratios (79 + 22 + 50 + 50); 1,766 tuples = 436 + 162 (= 598) + 492 + 676 (= 1,168); per-ratio counts equal
  G2's "N left" for all 201 (/tmp/z11/g2cmp.py); the seed-1 sample of 100 and the split 0..79 / 79..101 reproduced.
- pairing-killed 1,762, not 4; per curve 08 1280, 0L/0U/8L/8U 1694 each, LU 874: as stated. First killing primes as
  stated (17:420 ... 173:4).
- **"18,222 kill lines" is wrong wording.** There are 8,930 kill lines (one per killed (tuple, curve), = the sum of the
  per-curve counts); 18,222 is the number of local -1 terms on them (histogram of terms per line: 1 x 5182, 3 x 2886,
  5 x 826, 7 x 36; always odd, as it must be).
- "101 with m <= 997" (task wording) is wrong: the 101 are 101 of the 997 ratios of tuple_selmer_997.txt, with m from
  291 to 5000 (g2_997_all.txt); NOTES.md itself says "the 101 of g2_997_all.txt" (correct). The 836 have 5029 <= m <= 14977.

### The 4 exceptions are removed by a global argument after all (rank 0)
- On (5623,1017) E_08 and on (3320,2303) E_LU, dim Sel^2 = 2 = dim E[2] (t10's rank_second, my z7-based descent, and
  PARI's [0,0] agree), so **rank 0, proved by the complete 2-descent**. Torsion: 8 points (PARI's elltors as finder;
  gcd of #E(F_p), 3 <= p <= 400 good, = 8). The 8 points of C(Q) lie over t in {0, 1, -1, inf}
  (/tmp/z11/rank0.py), and at none of these are the four forms squares together (F_L, F_U fail at 0 and inf; F_0, F_8
  at +-1, as m^2 + n^2 is not a square). So **C_08(Q) (resp. C_LU(Q)) contains no image of a four-form point: these two
  ratios have no four-form point at all**, and the 4 tuples die by this rank-0 argument on the very curve where their
  class is a torsion class. NOTES.md's "no obstruction on E_08 can ever remove these tuples" and "E_08 (resp. E_LU)
  cannot" are misleading: E_08 removes them (and the whole ratio), not by the pairing but by finiteness.
- 5 of the 201 ratios have a pair curve with dim Sel^2 = 2 in t10's own lines: (3320,2303) LU, (5623,1017) 08,
  (9191,7688) LU, (9532,4445) LU, (14885,6076) LU; for all five C(Q) lies over {0, +-1, inf} and none is a four-form
  point (same script). These five ratios are excluded outright by a rank-0 pair curve (computed; tier D; m > 1000, so
  nothing for the census).

## ~19:44 UTC — more recomputation, controls, forging, the OOM output
- My code vs t10's lines (/tmp/z11/cmp.py: verdict per (tuple, curve), t10's printed class and killing b):
  (465,416) 36/36, (5623,1017) 12/12, (3320,2303) 12/12, (14885,6076) 48/48 agree; classes and killing b's identical
  wherever t10 prints them; both signs always give the same verdict. (14885,6076): all 8 tuples killed on
  08, 0L, 0U, 8L, 8U; E_LU there has dim Sel^2 = 2 (rank 0; see above). (3320,2303): not killed anywhere, LU class a
  torsion class (dim Sel^2 = 2). Outputs: /tmp/z11/run_*.txt.
- Degenerate controls, my code (pairing always computed, no "image" shortcut): (4,3), (12,5) (t10's) and my own
  (15,8) [17^2], (23,7) [2*17^2], (21,20) [29^2], (31,17) [2*25^2]: every known tuple survives every curve, both signs;
  (12,5)'s two extra tuples (good prime 29) killed on 0L, 0U, 8L, 8U, as in t10_negative.txt (24/24 lines agree).
  t10's committed script on my four new controls: known tuples survive (image; (31,17)'s second tuple "kernel" on
  the mixed curves) (/tmp/z11/t10_mycontrols.txt).
- Forging (try to kill a class of a rational point): (a) my code, every subset sum of PARI's points and the torsion on
  all six curves of (7,1), (4,3), (465,416) (e.g. all 32 classes of Sel^2 = image on E_0L of (7,1)), plus the same
  classes recomputed by the divisor formula on C (both signs): all in the kernel (/tmp/z11/forge_mine.txt). (b) t10's
  code with PARI's points withheld (PairData built on the 2-torsion and psi(iota t0) only), so the point classes and the
  known tuples go through t10's pairing: (7,1), (4,3), (465,416), (12,5), (627,364), (14885,6076): nothing killed
  (/tmp/z11/forge_t10.txt, forge_t10b.txt). The method never killed a class of a rational point in my tests.
- **The OOM output t10_997_a.txt**: I reran the committed script (no --write) on its LAST five ratios (--from 74 --to
  79: (4247,2193), (4333,508), (4348,3247), (4404,2197), (4424,2057)): every line of each block is byte-identical to
  t10_997_a.txt (/tmp/z11/t10rerun_997_74_79.txt). With t10's rerun of the first three, and my own code on (465,416)
  and (3320,2303) from the same file, I find the recovered output trustworthy.
(Section times corrected at 19:47 UTC from date -u; earlier headings were estimates.)

## ~20:05 UTC — wider recomputation; pairing matrices; middle of the OOM output
- My code (with z7's real-place search replaced by a draw in (e1, e2) and (e3, e3 + 3(e3 - e1)): z7's draw in
  [-10H, 10H] found no real point on a curve with clustered roots, (3163,3116)) on 38 of the 201 ratios: a seeded
  sample of 30 (random.Random(11), ratios with <= 40 tuples; /tmp/z11/batch.py, bat_*.txt) and (465,416), (5623,1017),
  (3320,2303), (14885,6076), (2999,1081), (3073,412), (3143,2567), (9751,5084) [216 tuples]. **2,472 (tuple, curve)
  verdicts (412 tuples), all equal to t10's; every class and killing b that t10 prints is equal to mine; both signs of
  every class give the same verdict.** No class is ever outside my Sel^2.
- (5623,1017), full pairing matrices with my code: E_0L dim Sel^2 6, matrix identically +1 (rank <= 4, as t10/PARI
  s = 0); E_LU dim 5, matrix rank 2 (rank <= 1, as PARI [1,1] s=2); E_08 dim 2 (rank 0). Root number of E_LU: -1.
- t10's committed script rerun on 997_a indices 35..40 ((2999,1081), (3007,697), (3073,412), (3136,1825),
  (3143,2567), (3163,3116)): byte-identical blocks. So 11 ratios of the recovered file (first 3 by t10, 6 middle and 5
  last by me) are reproduced by the committed script, and 10 of its 79 ratios by my own code.
- Distinct tests: the class on E_jk depends only on d_j d_k, so the 1,766 x 6 tuple-curve tests are 5,707 distinct
  (ratio, curve, d_j d_k), and the 8,930 kill lines are 4,976 distinct (ratio, curve, class) kills. Not an error, but
  "1,762 tuples" is not 1,762 independent events.
- "never above the largest good prime of G2's kill set": checked, the largest first killing prime EQUALS it for all
  201 (/tmp/z11/distinct.py).
- My H1 control (/tmp/z11/h1_mine.txt; base points m/n, -m/n, not t10's t = 0): (1109503,597671) — the two known
  tuples survive E_0L and E_0U, the four extra ones are killed on exactly the curve h1_pairs.txt names, both signs,
  same killing b (920312395136104482698617, 8302). Agrees with h1_pairs.txt.
- t10's own controls rerun with the committed script: `--h1control` and `--negative` print exactly
  t10_h1control.txt and t10_negative.txt (diff: only PARI's stack banner). Worktree left clean (git status empty).

## Item 1 — the logic, derived independently
Let C = C_jk: y^2 = F_j F_k (t), F = a conj(a), the four linear forms pairwise coprime (LAST_PAIR §9), P0 = (t0, y0)
rational with y0 != 0, psi: C -> E, P -> [P - P0] (Washington's map is such an isomorphism up to the automorphism -1,
which acts trivially on E[2]). E[2] = the three partitions of the roots, all Galois-stable: {F_j|F_k}, {a_j a_k|conj},
{a_j conj a_k|conj}; T_ab|cd = [W_a + W_b - D_inf]. The Kummer map is delta_r(Q) = f_r(D) for any rational divisor D in
the class of Q and any f_r with div f_r = 2 D_r, D_r rational in the class of T_r (independent of the choices by Weil
reciprocity; x - e_r is the choice on the Weierstrass model). With f_r = F_j(t), D = (P) - (P0):
delta_r(psi P) = F_j(t_P)/F_j(t0). Over K = Q(i), f = a_j a_k(t): delta_*(psi P) = a_j a_k(t_P)/a_j a_k(t0) mod K*^2, and
Q*/Q*^2 -> K*/K*^2 has kernel {+-1}: delta_* = +-D(t_P) D(t0) with D the rational class (Lemma S1). At a four-form
point with tuple tau: F_j(t_P) is a square and D(t_P) = d_j d_k (the factor v cancels; conj does not change a
rational class), so delta(psi P) = c(tau) = (F_j(t0), +-d_j d_k D(t0), product), for one of the signs. The signs differ
by (1, -1, -1); my code finds delta(psi(iota t0)) = (1, -1, -1) on the mixed curves and psi(iota t0) = T_r with
(*, ') = (-1, -1) on 08 and LU, so the difference is the class of a rational point, and by bilinearity both signs pair
alike (my code also tests both signs separately: same verdict on every line). delta(E(Q)) lies in the kernel of the
Cassels–Tate pairing on Sel^2; so <b, c(tau)> = -1 for some b in Sel^2 excludes tau. Base points: rational for every
ratio (checked by hand above). Roles e_r, e*: my finite-field test (independent of the resolvent formula) gives t10's
roles on every curve I ran. The four-form point is never P0, iota P0, a Weierstrass point or at t = inf for a
non-degenerate ratio (F_0(m/n) = 2m^2, F_L(0) = 2(m^2+n^2) are not squares; F_jF_k > 0). **The logic holds.**
What is not shown: NOTES' gloss that the obstruction lives on "the genus-one torsor a_j = d_j z_j^2, a_k = d_k z_k^2"
(that curve is never identified with the 2-covering of class c(tau); the computation is about the latter).

## Final (written ~20:07 UTC by date -u)
**Bottom line: t10's computation is right; its description of the 4 exceptions is wrong in effect, and two counts are
mislabelled. Verdict: holds with gaps.**

Verdicts.
1. Logic of the kill (c(tau), base points, signs, roles, kernel): **holds** (derived above).
2. 1,762 of 1,766 pairing-killed, per-curve counts, first killing primes, 199/201: **holds** (recount exact; 2,472
   (tuple, curve) verdicts on 38 ratios recomputed by my own code, all equal, classes and killing b's equal).
3. Controls: **hold** — t10's rerun identically; 4 new degenerate ratios (15,8), (23,7), (21,20), (31,17) survive in
   both codes; H1 reproduced by my code; rational-point classes never killed (my code: all subset sums of points on
   (7,1), (4,3), (465,416); t10's code with PARI's points withheld).
4. The 4 exceptions: not pairing-killed on any curve, torsion class on E_08 (resp. E_LU): **holds**; but "no
   obstruction on E_08 can ever remove these tuples" / "E_08 (resp. E_LU) cannot" **fails as a description**: E_08
   (resp. E_LU) has rank 0 (dim Sel^2 = 2), torsion 8, and its 8 points lie over t in {0, +-1, inf}, none a
   four-form point, so the two ratios have no four-form point at all. Same for (9191,7688), (9532,4445), (14885,6076)
   via E_LU.
5. OOM output t10_997_a.txt: **holds** (11 ratios byte-identical on rerun, 10 recomputed by my code).
6. Wording gaps: "18,222 kill lines" (8,930 kill lines carrying 18,222 local -1 terms); "every local -1 at a bad prime"
   is forced by the formula and depends on the local points (my choices put the -1's elsewhere), not evidence;
   "1,762" are 4,976 distinct (ratio, curve, class) kills; task wording "101 with m <= 997" is wrong (101 of the 997
   ratios of tuple_selmer_997.txt, 291 <= m <= 5000); "first ... mechanism" must not be used (brief rule 2).
   Tier D, "candidate", "computed (exact)" with its quoted inputs: appropriate.

Not checked: the earlier script version's agreement on 662 tuples (logs not committed); residual sets beyond reusing
G1's exact_local (as G2); the 163 ratios I did not recompute; E± curves (t10 did not try them either).

### Exact replacement text for NOTES.md
(a) Header line 3: "not reviewed." -> "reviewed by a separate AI agent (z11; `reports/z11.md`), not by a person."
(b) In **Method**, "(a Brauer–Manin obstruction on the genus-one torsor `a_j = d_j z_j², a_k = d_k z_k²`)" ->
"(the 2-covering of `E_jk` with class `c(τ)` is everywhere locally soluble and has no rational point; the pairing is
its Brauer–Manin obstruction from `Sel²`)".
(c) In the `Sel²` bullet, "the torsor `a_j = d_j z_j², a_k = d_k z_k²` is locally soluble" -> "the 2-covering of class
`c(τ)` is locally soluble".
(d) Replace the bullet "**Every −1 is at a bad prime of the ratio:** … are among those of the ratio." by:
"- **Where the −1's sit.** The 8,930 kill lines (one per killed tuple and pair curve, one `b` each; 4,976 distinct
  (ratio, curve, class)) carry 18,222 local terms −1, an odd number on each line, all at bad primes of the ratio, none
  at ℝ or at a prime only auxiliary for `rank_second` (the primes of `δ`, `det`). This is forced, not evidence: the
  bad primes of every `E_jk` (1,206 curves) are among those of the ratio, the local terms are 1 outside the bad
  primes of `E_jk`, 2 and ℝ, and which of those places carries a −1 depends on the local points and tangent lines
  chosen (review z11 found other places with other choices); only the product is intrinsic."
(e) Replace the sub-bullet "On `E_08`, `Sel²` is the image of the torsion (dim 2) and the class is that of a torsion
point, so no obstruction on `E_08` can ever remove these tuples." by:
"  - On `E_08`, `Sel²` is the image of the torsion (dim 2), so `E_08` has rank 0, and its torsion has order 8 (8 points;
    the gcd of `#E(F_p)`, `3 ≤ p ≤ 400` good, is 8: review z11). The 8 rational points of `C_08` are the obvious ones,
    over `t = 0, 1, −1, ∞`, where the four forms are not squares together. So the ratio has no four-form point: the
    tuples are removed by `E_08` after all, by rank 0, not by the pairing (their class is a torsion class, so no
    Selmer or Cassels–Tate condition on `E_08` sees them)."
and "The partner gives the same picture with `08 ↔ LU`." -> "The partner gives the same picture with `08 ↔ LU`
(`E_LU` of (3320,2303) has rank 0 likewise). Three other ratios of the 201, (9191,7688), (9532,4445), (14885,6076),
have `E_LU` of rank 0 and are excluded the same way (their tuples are also pairing-killed)."
(f) Replace the paragraph **What it means** by:
"**What it means.** For 199 of the 201 ratios, the Cassels–Tate pairing on `Sel²` of the genus-one quotients `C_jk`
removes every tuple that the good primes remove. For the partner pair (5623,1017) / (3320,2303) the pairing leaves
two tuples each, but there `E_08` (resp. `E_LU`) has rank 0 and its 8 points lie over `t ∈ {0, ±1, ∞}`, so these two
ratios have no four-form point at all (review z11). So on these 201 of the 937 ratios that need a good prime, every
tuple that a good prime removes is also removed by a computation on one pair curve over ℚ (the pairing, or rank 0).
That is a candidate global mechanism, not a proof of one: the pairing does not reproduce the good-prime kills of the
degenerate controls (627,364), (991,263), and the residual sets rest on `exact_local.py`. The "E±" curves of the
question were not tried."
(g) Answer line: after "and **4 are not**." add " (those 4 are removed by rank 0 of `E_08`, resp. `E_LU`; see below)".
Task/LEDGER wording: say "the 101 of the 997 ratios of `tuple_selmer_997.txt`" (not "m ≤ 997"); never "first".

Reports/outputs: /home/claude/conductor/reports/z11.md; scripts and outputs in /tmp/z11/ (z11_pair.py, wash_check.py,
recount.py, g2cmp.py, cmp.py, rank0.py, forge_t10.py, forge_mine.py, h1_mine.py, pmat.py, batch.py, run_*.txt,
bat_*.txt, t10rerun_*.txt).
(Side check stopped unfinished: a search for a point of infinite order on E_LU of (5623,1017) (root number -1, rank <= 1)
with ellrank effort 4 / ellheegner, to see whether the exceptional classes are point classes there too; not needed for
any verdict.)

## Packaging (conductor's request, 20:10–20:33 UTC): /home/claude/conductor/r7/referee_z11/ (to be experiments/referee_z11/)
Scripts print exactly their .txt (no timings, PARI messages kept off stdout; `--write` rewrites; stderr says
`# matches`/`# DIFFERS`); paths relative (../.., ../conductor-t10, ../referee_z7, or MSS_REPO/T10_DIR/Z7_DIR); checked
CI-style in place and in a simulated repository layout without environment variables. rank0.py no longer uses
elltors/ellap: torsion bound by exact #E(F_p) counts (gcd 8 over the good p <= 199) and the 8 explicit points of C.
Times (CPU): rank0 1 s, recount 0.1 s, wash_check 0.5 s, pmat 4 s, h1_mine 10 s, forge_t10 16 s, z11_compare 40 s,
z11_pair 42 s, z11_compare --set mid 66 s, forge_mine 70 s (all CI); z11_compare --set big 228 s, --set batch30 351 s,
oom_check 335 s (not CI). All outputs as in the review (batch30: 960/960, big 1296/1296, oom 11/11 identical).
