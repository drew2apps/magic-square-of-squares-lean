# Hostile review of run4_w6 (lball_exact53.py; commits 5bdd931, 50f335f)

Reviewer: AI referee, 01:39-02:16 UTC 26 Nov 2026 (about 15 CPU-minutes used). Repository not edited; reviewer scripts in /tmp/rev_w6/scr.

## Verdict: HOLDS WITH GAPS

The mathematics holds: for all 53 curves, L(E,1) > 0 is certified from N (Tate), the a_n (counting / order
certificate / Tate) and w = +1 (F != 0), with a rigorous enclosure; I found no error in the enclosure, the order
certificate, the a_n recursion or the bad-prime a_p, and every independent comparison I made agrees (N and w of all
53 against PARI; a_p at 3579 random primes against ellap; PARI's lfun(E,1) inside all enclosures checked; reruns
byte-identical). The gaps are (i) a scope problem: these 53 are tier-D second-descent closures, NOT census closures,
while the census's own L-value closures (G1_FULL_SYSTEM.md:205, e.g. (112,15), (127,97), (700,529)) are not
reproduced by global.py and are covered by neither lball.py nor this script; (ii) the negative test tests copies of
the guard and an unused function rather than the code that certifies.

### FIX-FIRST
1. Scope (outside the diff, but it is the premise of N1): G1_FULL_SYSTEM.md:205 counts 19 census closures that use
   L(E,1) != 0; global.py finds no rank-0 curve for exactly 18 of the 183 counted trivial-point closures (list
   below) and for (700,529) (G1_FULL_SYSTEM.md:211 says it reproduces every closing curve), and no script computes
   those L-values. Find them, certify them with
   lball_exact53's machinery, or label them. Do not describe the 53 as "counted" (LEDGER.md:54,56; SECOND_DESCENT.md:3,104).

### FIX
2. --negative (lball_exact53.py:491-497) re-implements the guard in `two()` instead of calling LX.evaluate (and
   omits the Fricke residual): route forged N / a_p through evaluate or certify.
3. 120 of the "147 refused" are verify_order (lball_exact53.py:226-239), which the certificate never calls; count
   them apart, and add mutations of ap_order's own steps (unique_multiple, Hasse radius) -- my mut.py shows
   --aptest catches those, so this is about the test's claim, not soundness.

### NIT
4. Docstring (:30-32) calls the guard "not part of the proof"; it is what refuses a wrong N (asserted at :350, :451). Say so.
5. run() gates on agreement with the PARI-derived lball.txt (:350): conservative, but "without PARI" acceptance still reads a PARI output.
6. Line 53 ((751,289)) rests on one Mac run with an earlier copy (certificate path identical; diff read); say so.
7. Build report's doc text: "141" -> 147; the "existing row for lball.py" in LEDGER does not exist.
8. count_ap fallback at p ~ 8.7e7 allocates a p-element list (~700 MB): a crash, not a wrong value.


## Findings (written as found)

### Static reading of lball_exact53.py and the parts of lball_exact.py it uses (01:40-01:50)

- Enclosure (lball_exact.py:430-544): series L(E,1) = sum a_n/n (q_y^n + w q_{1/y}^n); at y = 1, w = +1 this is
  2 sum a_n/n q^n (the classical Cremona formula). Checked by hand: pi by Machin with bracketed truncation
  (:434-451), sqrt N by isqrt with s <= sqrt(N) 2^P < s+1 used in the right directions (:469-473), exp(-t) by
  alternating partial sums stopped after a negative / positive term for lower / upper (:453-467; t < 1 asserted,
  so terms decrease), q^n by mul_lo / mul_hi, each term rounded outward by sign (:479-488), tail
  2 q^(M+1)/(1-q) for |a_n|/n <= 2 and 2 q^(M+1)((M+1) - M q)/(1-q)^2 for |a_n| <= 2n, both with q^(M+1) from
  above and 1-q from below (:489-494). |a_n| <= d(n) sqrt n <= 2n is right (d(n) <= 2 sqrt n). M is chosen by
  floats but the tail actually subtracted is the exact integer bound, so a float error can only widen, not
  falsify. No defect found.
- w = +1 from F = sum a_n e^{-2 pi n/sqrt N} excluding 0 (lball_exact53.py:313): valid given that N is the
  conductor and the a_n are those of E (Fricke at the fixed point i/sqrt N). If N were wrong, F != 0 would prove
  nothing; what catches a wrong N is only the guard (two cut-offs + Fricke residual at y = 6/5), which the
  docstring calls "a check, not part of the proof" but which run() does assert (lball_exact53.py:350). So in
  practice a wrong N or w is refused, not silently accepted (see the negative test below). NIT: the docstring
  should say that the guard is what protects against a wrong N (Tate's algorithm itself is trusted), rather than
  calling it inessential.
- Group law (ec_add :73-87) is Silverman III.2.3 verbatim; doubling and the -P test are right.
- Order certificate (bsgs_multiple :173-190, exact_order :192-199, unique_multiple :204-208, ap_order :215-224):
  the argument is sound. m with mP = O is re-verified by ec_mul (so BSGS collisions cannot produce a wrong m);
  the exact order is m stripped of primes l while (n/l)P = O, then each remaining l asserted to be needed;
  L = lcm divides the group exponent hence #E(F_p) (Lagrange); hasse(p) = [p+1-isqrt(4p), p+1+isqrt(4p)] is
  exactly the set of n with (n-p-1)^2 <= 4p; unique multiple of L in that set => #E(F_p). When E(F_p) has small
  exponent the lcm never exceeds the exponent, uniqueness fails, and after 12 points the prime is counted
  (count_ap). So small exponent can cost time, never correctness. factor_small needs m < 2^32; with M <= 8.7e7
  the largest m is < 1e8 (asserted). The only trusted pieces are ec_add/ec_mul and Python integers.
- Good-reduction model: at p | Delta(a) with good reduction, the p-minimal model from Tate is used, and
  model_mod asserts Delta(model) != 0 mod p at every prime where the certificate is used (:210-213), counting
  asserts the same (:247). Correct.
- a_n: an_from_ap (lball_exact.py:407-426) and an_array (:402-419 of lball_exact53) implement multiplicativity,
  a_{p^k} = a_p^k at bad p and the Hecke recursion at good p; the spf sieve in an_array (reversed primes, slice
  assignment) gives the least prime factor of every composite; checked by reading.

### Reruns (01:45-01:50, cloud, 2 cores shared with the build agent's own (41,1) rerun)
- `--only 0,1,2,6,37 --nopari` ((21,20), (15,8), (117,44), (340,189), (684,37) with N = 1.02e10): the five lines
  are byte-identical to lines 1, 2, 3, 7, 38 of lball_exact53.txt (41 s wall, 26 s CPU).
- `--negative`: byte-identical to lball_exact53_negative.txt (12 s). `--controls`: identical (0.6 s).

### The negative test: what it mutates
- Forged N (N/l, lN, 2N), forged w, forged a_p at a small good prime, at a good prime near sqrt(N)/4 and at the bad
  prime 2: each is refused because the two cut-offs disagree, the right reason. But the test re-implements the
  guard in a local `two()` (lball_exact53.py:491-497, copied from lball_exact.py:640-646) instead of calling
  `LX.evaluate`, the code the certificate runs, and it checks only the overlap, not the Fricke residual. A bug in
  `evaluate` (e.g. a guard that always passes) would not be seen by --negative. FIX: route the forged inputs
  through `LX.evaluate` (or `certify` with a forged LD / N) and assert `guard` is False.
- "w = -1" is not a pipeline input (w is derived from F); the real failure mode, a curve whose true w is -1, is
  covered by --controls (79a1, 5077a1: F contains 0, refused). Fine.
- Forged #E(F_p) + d: refused by `verify_order` (lball_exact53.py:226-239), a function that the certificate never
  calls (ap_order computes the order; nothing is "claimed"). These 120 of the 147 refusals therefore say nothing
  about the code that produced the a_p. FIX: either drop them from the count, or mutate ap_order's own pieces
  (Hasse radius isqrt(4p) -> isqrt(4p)+k or isqrt(p); unique_multiple returning the first multiple even when it
  is not unique; exact_order not stripping a prime) and require detection against counting. A mutation of
  unique_multiple would be the important one: it is the step on which uniqueness (hence soundness) rests.
- Broken group laws (two): these do go through ap_order (the global ec_add is swapped) and are detected only by
  comparison with counting (the self-check mP = O uses the same broken law, so it cannot see it). Fine as a test;
  it confirms that the soundness of the order certificate rests on ec_add being the group law, as stated.
- Not mutated: Tate's output itself (an f_p or a split/non-split flag inside local_data, rather than N afterwards);
  forging N afterwards covers the effect on the series, so this is a NIT.

### Independent cross-checks with PARI (reviewer's own scripts in /tmp/rev_w6/scr, not in the repo)
- indep.py on 9 curves (indices 0, 1, 2, 7, 25, 27, 37, 51, 52, i.e. the three smallest, (23,7), (399,40),
  (552,385), (684,37) with N = 1.02e10, (41,1) and (751,289)): N from local_data equals PARI's ellglobalred on the
  minimal model for all 9; PARI's ellrootno is +1 for all 9 (w_local +1 for 8, not implemented for (41,1)); Tate's
  a_p at every bad prime equals PARI's ellap; ap_order at 60 random primes in [1000, Mmax] per curve (for the two
  large curves up to 2.6e7 and 8.7e7) equals PARI's ellap: 540 primes, 0 mismatches, 0 unsettled.
- smallexp.py (adversarial for the order certificate): on (21,20), (23,7), (552,385) (torsion Z/2 x Z/2), every
  prime 1000 <= p <= 2e5 at which the exponent of E(F_p) (PARI ellgroup) has two or more multiples in the Hasse
  interval: 3, 1 and 1 such primes; ap_order returns None (falls back to counting) at all 5, never a value. This
  is the case where a weaker certificate would go wrong; this one does not.
- (41,1): the build's scratch a_n file (work/big51__an.bin, 25.6M terms, from which line 52 was merged) agrees with
  PARI ellak at 1500 random n and ellap at 1500 random primes; the build's first a_p run and its current rerun with
  the committed script agree on the 1.4M a_p computed so far by the rerun.
- PARI's lfun(E,1) (a different algorithm from lball.py's sum) lies inside the 12-digit enclosure for (21,20),
  (15,8), (23,7), (195,28) and (552,385), the smallest L-value (0.0266132796); the run on all 53 is below.
- The Mac copy that produced line 53 ((751,289); work/lball_exact53_macrun.py, md5 4667...) differs from the
  committed script only in docstrings, verify_order / negative, merge's summary, and the extra good-reduction
  assertion at p < 1000 (diff read): the certificate path is the same code. Line 53 rests on that single run (the
  PARI ellap check at all 5,042,213 primes is recorded in the jsonl). NIT: say so in the LEDGER row.

### Item 5: which closures these are, and whether they are counted (01:48-01:52) -- READ THIS
- The 53 are exactly lball.txt's 53, in the same order (asserted in merge, lball_exact53.py:359; I re-derived the
  selection: closed rows of second_123.jsonl, curve `by`, N <= 1e14), and each integral model has the same PARI
  minimal model as the `rank0_curve` it came from (53 of 53; /tmp/rev_w6/scr/isom.py). All 53 enclosures lie
  inside lball.txt's intervals (checked from the jsonl, not only "meet"); smallest lower end 0.026613279413
  at (552,385), as the build report says.
- BUT these 53 are second-descent closing curves, and the repository says that the second descent is tier D, not
  counted: LEDGER.md:56 ("Nothing in tiers C and D is counted in it"), LEDGER.md:54 (the second descent of
  SECOND_DESCENT.md is "tier D"), SECOND_DESCENT.md:3 and :104-105 ("Nothing here is counted in the census"),
  and SECOND_DESCENT.md:69 (81 of the 89 closures re-close ratios the census already closes by other methods). The
  same holds in int4 (LEDGER.md:58 there). So the premise of referee N1 / this brief ("the 53 counted L-ball
  closures") does not match the repository: an error in lball_exact53 would not move a census number.
- The census does have L-value closures of its own, and they are NOT these: G1_FULL_SYSTEM.md:205 says "19 of the
  closures use L(E,1) != 0 (Kolyvagin) where 2-descent gives only [0, 2]", e.g. (112,15)~(127,97) and (700,529),
  conductors below 1e14. Those three pairs are among the 425 congruence survivors (I re-ran the closed() test of
  gen_census1000.py: none of the congruences closes them) and are not in open_pairs_1000.txt, so the census counts
  them as closed by a rank-0 quotient. Yet global.py, which G1_FULL_SYSTEM.md:211 says "reproduces the closing
  curve for each closed pair", returns no closing curve for any of the three (I ran global.analyse: ranks [0, 2]
  for w*f0*f8, g4*K, f8*K, w*f0*K, w*g4*K of (127,97); no curve with upper bound 0). No script in the repository
  computes an L-value for the G1 curves (grep for lfun / ellL1 / ellanalyticrank / Kolyvagin in scripts/g1_full and
  archive: nothing). So the counted L-value closures appear to rest on a computation that is not in the repository
  and that neither lball.py nor lball_exact53.py touches. I checked 3 of the 19 only. This is outside the w6 diff,
  but it is exactly the trust the brief is about: FIX-FIRST for the conductor (find the 19, put their L(E,1) > 0
  through lball_exact53's certificate, or record them as resting on PARI's L-values).
- The build report's proposed doc text says "141 forged inputs refused" twice; after 50f335f it is 147 (NIT). Its
  proposed edit "in the existing row for lball.py / the census L-ball route" has no target: LEDGER has no tier-B row
  for lball.py; the only place is SECOND_DESCENT.md section 5 (NIT, but it propagates the "census" misnomer).

### Mutations of the order certificate's own pieces (reviewer's mut.py, 01:53)
Monkey-patched in a copy of the module namespace, run over every good prime 1000 <= p < 30000 of (21,20), (23,7),
(399,40) (the prime range of --aptest), compared with PARI ellap:
- unique_multiple returning the first multiple even when not unique: wrong a_p (or an assertion) at 472 of 9231
  primes -> --aptest would refuse it.
- Hasse radius isqrt(p) instead of isqrt(4p): 3658 of 9231 -> refused.
- Hasse radius isqrt(4p) + 3 (wider: still sound, only weaker): 0 of 9231 wrong, as it should be.
So the dangerous mutations are caught by the existing --aptest comparison (CI "descents"), even though --negative
does not contain them. That lowers the verify_order point above to FIX (the count of "147 refused" overstates what
--negative tests), not a soundness gap.
- (751,289): ap_order against PARI ellap at 400 further random primes in [5e7, 8.68e7], the top of the range
  used: 0 mismatches, 0 unsettled.

### The census's own L-value closures, sampled (01:55, /tmp/rev_w6/scr/g1sample.py)
Recomputed the 425 congruence survivors (gen_census1000.closed), 306 with a trivial point, minus the 123 listed
open = 183 counted closed by a rank-0 quotient. global.py (the census program for Result 3) on 30 random ones of the
183: no closing curve for 3, (575,48), (487,217), (697,696) (plus (112,15), (127,97), (700,529) checked before).
3 of 30 matches ~18 of 183, i.e. the "19 closures that use L(E,1) != 0" of G1_FULL_SYSTEM.md:205. These counted
closures are not reproduced by any program in the repository, contrary to G1_FULL_SYSTEM.md:211, and rest on
L-values (presumably PARI's) that neither lball.py nor lball_exact53.py computes. This is the real "counted
closures that trust PARI's L-values". (Not exhaustive: 6 pairs found; the full list needs global.py over the 183,
about 7 CPU-minutes.)

### PARI lfun on all 53 (01:54)
PARI's lfun(E,1) (30 digits requested) lies inside lball_exact53's 12-digit enclosure for all 53 curves, including
(41,1) (1.39295152987733, 75 s) and (751,289) (3.45121464818643, 208 s). (lfun uses PARI's own a_n and root number,
so this checks the arithmetic of the enclosure, not the PARI-freeness.)

## What I could not check
- I did not rerun lines 52-53 ((41,1), (751,289)) end to end (30 and 60 CPU-minutes of a_p); I checked their N, w,
  540 + 400 random a_p, the (41,1) a_n file against PARI at 3000 points, and their values against lfun.
- Tate's algorithm itself was read only where it sets a_p and f (multiplicative split test, additive a_p = 0);
  its correctness for these 53 rests on agreement with PARI (N for 53 of 53, a_p at every bad prime for 53 of 53),
  as in the earlier review of lball_exact.py.
- The theorems assumed (modularity with level = conductor, the Fricke eigenvalue relation w = -eps, Kolyvagin or
  Kato) are standard and stated correctly in lball_exact.py:34-56; not re-derived.
- Whether the closing curve `by` of each second_123 row is the right curve for the closure (upstream, the second
  descent) is outside this diff.

### The full list (02:05): global.py over all 183 counted rank-0 closures with a trivial point
No closing curve for exactly 18: (112,15), (127,97), (240,161), (260,69), (288,175), (329,191), (352,135),
(401,79), (463,113), (487,217), (527,336), (544,33), (575,48), (577,511), (623,527), (697,696), (863,191),
(943,576). That is exactly the "18 pairs with a trivial point" of G1_FULL_SYSTEM.md:205; with (700,529) (no
trivial point, checked separately) the 19. None of these 19 is among the 53 of lball_exact53 (the 53 are
second-descent curves of ratios that the rank-0 step left OPEN, e.g. (15,8), (21,20)). So: the counted census
closures that trust an L-value are these 19, their L-values are computed by no program in the repository, and
G1_FULL_SYSTEM.md:211 ("global.py reproduces the closing curve for each closed pair") is false for them.
- (02:15) global.py over all 85 counted rank-0 closures without a trivial point: no closing curve for exactly 1,
  (700,529). So the 19 of G1_FULL_SYSTEM.md:205 are confirmed one by one: 18 + (700,529).

## Reviewer artefacts (not in the repo)
/tmp/rev_w6/scr/{indep,smallexp,big51,lfun,isom,mut,g1sample,g1sample_nt}.py and their outputs in /tmp/rev_w6
(indep_all.txt, smallexp.txt, lfun_all.txt, g1all.txt, g1nt.txt, sample.txt, neg.txt, ctl.txt).
