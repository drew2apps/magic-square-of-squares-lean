# Review of run5_x2 (01e35dd..run5_x2) by review agent y2 (run 5)

Written by an AI model; not read by a person. Hostile-referee review; scratch under /tmp/y2/. No repository file edited.

Started 21:24 UTC. Commits under review: e8ffa87, 4f338a2, c75190c, 706b455, e0e4488.

## Runs from a clean clone of run5_x2 (e0e4488) at /tmp/y2/clean, 2 CPUs, timed
- `torsion_exact.py` 1.2 s, last line `self-test: all as expected`; `--negative` 0.7 s, `20 forged certificates, 20 refused`.
- `census_lvalue_exact.py --torsion | diff - census_lvalue_exact_torsion.txt`: identical, 2.5 s.
- `--torsion-negative`: identical to committed `_torsion_negative.txt`, 1.4 s (15 of 15 refused).
- `rank_fisher_theta4_order.py` 0.6 s, `--negative` 0.8 s: both identical to committed files.
- With `sys.modules['cypari2'] = None` (any cypari2 import raises), `--torsion` and `--torsion-negative` still run and
  give byte-identical output: the torsion step does not touch PARI.

## (a) torsion step -- findings so far (about 21:30 UTC)
Argument (b2) of torsion_exact.py re-read: correct. Minimal-order P in E(Q)[2^inf] outside T; R = 2P is in T; if R = 2S
with S in T (any S, not only 2-power ones -- the code's `dbl` uses all of T, fine) then P - S in E(Q)[2] in T; otherwise
R is "lonely" and R~ = 2P~ in 2E~(F_p) for every good odd p, contradicting the listed p. This does exclude a point of
order 8 over a listed order-4 point (and any longer chain, by minimality). P of order 2 outside T is excluded by the
exact E(Q)[2] count. Z/2xZ/4 vs Z/8 vs Z/2xZ/8: T is an explicit closed set, so the structure is read off T; my extra
forgeries (15a1 T = cyclic <P4>: refused, E(Q)[2] not in T; 15a1 T = E[2]: refused; 15a4 T = its Z/4: refused;
210e2 T = its Z/2xZ/4 and T = cyclic <P8>: refused) all behave correctly.
"Good" is decided by p odd and p not dividing Delta of the actual integral model (0,a2,0,a4,a6) (formula checked; my
own Delta = 16 disc(cubic) agrees on all 19). Primes dividing Delta of this model but good for the minimal model are
merely skipped (conservative). Reduction injective on torsion at good p >= 3: correct citation (E_1(Q_p) torsion-free
since 1 < p - 1). The model is integral (Curve casts int); torsion points need not be assumed integral: the search may
miss points but then certify() refuses, it cannot certify an incomplete T.
global.py by ast: the seven functions are compiled into a namespace with Fr, math, json, sys only; none references
`pari`; analyse's call to torsion_w resolves to the namespace entry torsion_step installs. Runs with cypari2 blocked are
byte-identical. Sound. torsion_ws reproduces global.torsion_w's (X -> w) lines; all 19 curves are the cubic case
(P = w*K and w*f0*f8 have degree 3), so the quartic branch is not exercised by the 19.

Independent re-derivation (/tmp/y2/indep/tors_indep.py, own model from (m,n), own Delta, own Legendre-symbol point
counts, p < 600): B agrees on all 19. Completeness by a DIFFERENT argument: w*f0*f8 is y^2 = x(x+m^4)(x+n^4) (full
2-torsion), so E(Q)[2^inf] is computed completely by the exact halving criterion (x0 - e_i all squares; Washington
Thm 8.14, quoted) -- gives Z/2xZ/4 (B = 8) or Z/2xZ/8 (B = 16); w*K is y^2 = x(x^2 + a x + S^4) with a^2 - 4S^4 =
-64 m^2 n^2 (m^2 - n^2)^2 < 0, one rational 2-torsion point, so E(Q)[2^inf] is cyclic, and E~(F_p)[2^inf] has
exponent 4 at some good p (17, 41, 97, ...), so it is Z/4 (the points x = S^2 exhibited). Odd part of B is 1 for all
19. My torsion x-multisets equal torsion_exact's on all 19, and PARI elltors (guide only) agrees on all 19.
Three curves in detail: (112,15) w*K Z/4 (exp 4 at p=17); (623,527) w*K Z/4 although B = 16 (exp 4 at p=41);
(575,48) w*f0*f8 Z/2xZ/8, complete by exact halving.

GAP (library robustness, not the 19 results): torsion_exact.bound() checks only "p odd and p does not divide Delta",
not primality, and accepts an empty list. Forged: 15a1 (true Z/2xZ/4) with primes = [] gives B = 0 and T = {O} is
CERTIFIED; with primes = [21, 33] (composites), B = 1 and T = {O} is CERTIFIED; prime powers (243, 289, 343) make
certify() raise ValueError, contradicting its docstring ("never raises on a forged input"). The census and
rank_fisher callers pass good_odd_primes() (a sieve), so no certificate of the 19 is affected. Fix: in bound(),
assert the list is non-empty and each p is prime (trial division), and add these three forgeries to --negative.
- Heavier lines (same clean clone): `kctp_isog_exact.py` 2.6 s and `--negative` 1.4 s, both identical to committed;
  `census_lvalue_exact.py --negative` 47 s, identical to `census_lvalue_exact_negative.txt`; `lball_exact.py --all`
  58 s, identical to `cat lball_exact.txt lball_exact_negative.txt`; `lball_exact53.py --negative` 10 s, identical to
  `lball_exact53_negative.txt`. All rc = 0.

## (b) lball negative tests through evaluate -- HOLDS
Old (01e35dd) vs new `_negative.txt`: after deleting the new ", residual [...]" field and the "(reason)" suffix the files
are line-for-line identical (lball_exact: 15 lines; lball_exact53: 42 lines); all 24 resp. 54 printed L(y=1), L(y=6/5)
enclosures are byte-identical. Reason: evaluate() uses the same M1 = terms_needed(N,1,True,tL) and
Mb = terms_needed(N,6/5,True,tL) as the old copy; only the array cap Mmax changed (now max_terms, which also covers the
F/residual term counts evaluate asserts). Every forged input (10 + 21 = 31) is refused by ev['guard'] = overlap and
R_has_0, the same expression main_curve()/certify() test; in every case BOTH sub-tests fire ("enclosures disjoint,
residual excludes 0"); true inputs pass with residual containing 0. Small caveat for the wording: forged N and w are
run with tails below 10^-3 (tL = NEG_TL), the certifying runs with TARGET; same function, different tail target.

## (c) G of infinite order on E+(72,65) -- HOLDS
Own code (/tmp/y2/indep/order_check.py): G on y^2 = x(x^2 + a x + b) with the printed a, b (exact); my Delta, my
Legendre counts give B = 2 over the 161 good odd p < 1000, and y(G) != 0, so (1) holds; my orders of G mod
11, 17, 19, ... are 10, 8, 6, 28, 6, 34, 22, 6, 52, 58, 6, 6 (identical to the script), so (2) holds; 2G integral,
3G has x = 2460449975530023264078722/2016607325329 (= .../1420073^2), so (3) holds (Nagell-Lutz on an integral model
y^2 = x^3 + a x^2 + b x, quoted). G = sigma_evidence's x = 168914884320 shifted by -e = 1665944009858: consistent.
PARI (guide): elltors 2, ellheight(G) = 3.107. The certificates are each sufficient; rank = 1 still rests on
rank <= 1 from rank_fisher_theta4 (not reviewed here).

## (a) continued
- My own classify (global.py's definitions re-typed) on my own torsion lists: for all 19 curves the only torsion w that
  are rational squares are 0 and 1 (plus w = inf from O), and the verdicts per t in {0, 1, -1} are exactly those
  printed in census_lvalue_exact_torsion.txt (e.g. (700,529): no/no/no; (112,15): no/degenerate/degenerate). No
  SOLUTION on any of the 19.
- The --torsion run does not recompute rank 0: it reads the committed census_lvalue_exact.txt lines ("certified
  positive" + "guard passes", or "certified by isogeny") and asserts one exists per curve. Fine as a cross-reference,
  but the doc should say "rank 0 as recorded in census_lvalue_exact.txt (its own CI line)". Those L-value lines carry a
  PARI cross-check column (lfun, N, w agree) produced with --pari: a check, not used. The pair->curve table CLOSURES
  came from the PARI-guided census_lvalue_identify.py (a search, not part of the proof). So "PARI nowhere" is right for
  the certificate, not for the provenance of the table; suggested wording below.

## (d) kctp_isog_exact.py -- HOLDS (experiment tier, as labelled)
Argument re-read: area L(C')/area L(C) = |u|^2 with omega' = u omega (Silverman III.1: x = u^2 x' + r gives
omega' = u omega, Delta' = u^-12 Delta); lattices of one curve over C are homothetic and Delta(L) = g2^3 - 27 g3^2 has
weight -12, so |u|^2 = |Delta(C)/Delta(C')|^(1/6) = (N Delta(C) / N Delta(C'))^(1/12) for K = Q(i) imaginary
quadratic, independent of any software convention for u; a normalised Velu isogeny gives L(Emin) in L(C) of index 2,
area(C) = area(Emin)/2. So area(C')/area(Emin) = rho/2 with rho the exact 12th root; rho = 2 for all three, ratio 1.
The rejected 1/(2 N(u)) = 1/4 would make B(E3) = B(E). The u-consistency test is only a check (the ratio never uses u).
Independent cross-check with PARI over bnfinit(t^2+1) (/tmp/y2/indep/kctp_pari2.py, a guide): Emin is minimal
(change [1,0,0,0]), prod c = 4194304 = 2^22, |tors| = 8; for the three Velu codomains: u = 1 - i, N(u) = 2,
N(Delta C)/N(Delta C') = 4096 = 2^12, rho = 2, area ratio 1; prod c(C') = 4096, 65536, 16777216; |tors| = 4, 4, 8;
B(C')/B(E) = 1/256, 1/16, 4 -- identical to Sage's. So B(E3)/B(E) = 2^24/2^22 * 8^2/8^2 * 1 = 4 is forced.
Rests on: Sage's exact global minimal models over Q(i), Tate's algorithm (Tamagawa numbers), torsion over K, and
is_normalized() (PARI agrees on every number); the standard facts (i)-(iii); for the Sha reading, Cassels' isogeny
invariance with rank E(K) = 0 and Sha finite (hypotheses). Reading: |Sha(E3)|/|Sha(E)| = B(E)/B(E3) = 1/4, i.e. the
corner curve E = L1L3M1M2 has the LARGER Sha (4 times that of E3), conditional on rank 0 and finite Sha.
Small wording point: the header says "(Checked: ... so the ratio is N(u).)" -- that is area(C')/area(C); the printed
area(C')/area(Emin) is N(u)/2. Clear enough in context.

## More attacks (about 21:35 UTC)
- Same library gap reaches rank_fisher_theta4_order.certificates(): with primes = [] (B = 0), certificate (1)
  "k G != O for every k <= B" HOLDS for the torsion point (0,0) (wrong); (2) and (3) still refuse. Real run unaffected
  (good_odd_primes(E, 1000)). Add primes = [] (and a composite list) to its --negative after fixing bound().
- TE.Curve casts coefficients with int(): Curve(Fraction(1,2), ...) silently becomes a different curve. Both callers
  guard it (census: assert tuple(a) == E.model(); rank_fisher: a.denominator == 1), so no effect; an assert on
  integrality in Curve would be cleaner.
- Determinism: --torsion, --torsion-negative, torsion_exact --negative, rank_fisher_theta4_order --negative are
  byte-identical under PYTHONHASHSEED 0, 3, 999, random.
- Base census_lvalue_exact.py had two identical copies of c46/isomorphic/two_isogenous/isogeny_chain (checked by AST);
  removing one changes nothing.

## Verdicts
(a) HOLDS WITH GAPS. The 19 torsion groups and the no-solution verdicts are right (re-derived by a different exact
    completeness argument; PARI agrees), the halving argument is correct, PARI is not touched (cypari2 blocked: same
    bytes). Gaps: (1) certify()/bound() accept an empty or composite prime list and then certify a wrong T (15a1:
    T = {O} certified with [] and with [21, 33]); unreachable from the shipped callers, but the negative suite should
    cover it and bound() should check primality and non-emptiness; (2) the --torsion run reads rank 0 from the
    committed census_lvalue_exact.txt rather than recomputing it (fine, but say so); (3) (b2) and the injectivity
    theorem are written/quoted, not machine-checked (x2 says this). Quartic branch of torsion_ws not exercised (all 19
    are cubic), harmless.
(b) HOLDS. Enclosures byte-identical; every forgery refused by evaluate's guard, both sub-tests firing.
(c) HOLDS. All three certificates re-checked with my own code; rank = 1 still needs rank <= 1 from rank_fisher_theta4.
(d) HOLDS (experiment tier). B(E3)/B(E) = 4 is forced by discriminant norms + normalised Velu; PARI over Q(i) gives the
    same minimal models' discriminant norms, Tamagawa products, torsion orders and the ratios 1/256, 1/16, 4. Reading:
    under rank 0 and finite Sha, |Sha(E)| = 4 |Sha(E3)| (E the larger), via Cassels.

## Proposed wording changes
- LEDGER row 51 / G1_FULL_SYSTEM line 205 / CHANGES: after "so those 19 closures use PARI nowhere" add "(rank 0 as
  certified in census_lvalue_exact.txt, with modularity and Kolyvagin; the injectivity of reduction on torsion is
  quoted and the 2-part completion step of torsion_exact.py is a written argument checked by its tests, not by a
  machine; the pair-to-curve table came from a PARI-guided search, which is not part of the proof)".
- torsion_exact.certify docstring: "never raises on a forged input" -> true only for prime lists of good odd primes;
  until bound() checks primality, say "the prime list must be primes (callers pass good_odd_primes)".
- x2's proposed LEDGER row 86 text (item 4, not in my brief, not checked by me beyond the ratio): its "cannot prove
  rank 0 on four of them ... only if n = 3" is conditional on rank E(K) = 0 and Sha finite and uses an n not defined
  in that sentence; keep "if rank 0 and Sha finite" and define n (Sha(E)[2^inf] = (Z/2^n)^2) in the row.
- `census_lvalue_exact.py --torsion --pari` (check only): 19/19 "the same torsion w: True", 11 s. Self-test controls:
  PARI conductors 15, 15, 15, 210, 14, 30, 90, 11, 26, 54, 66, 11 and torsion match the labels claimed.

## Final (y2)
All proposed CI lines reproduce byte-for-byte from a clean clone; timings above (all under 1 min; slowest are the
existing lball_exact.py --all 58 s and census_lvalue_exact.py --negative 47 s). Verdicts: (a) HOLDS WITH GAPS (library
accepts empty/composite prime lists -- unreachable from the shipped callers), (b) HOLDS, (c) HOLDS, (d) HOLDS
(experiment tier). Scratch code: /tmp/y2/indep/ (tors_indep.py, cmp.py, wcheck.py, order_check.py, kctp_pari2.py,
forge*.py). No repository file edited.
