# Review w7: the 19 census L-value closures (hostile referee)

## Verdict: HOLDS WITH GAPS

The mathematics holds. The 19 pairs are exactly the census's non-ellrank closures, the curves are correct, and the
isogeny step is valid. What needs fixing: one false sentence in the new doc text, a qualifier missing in four places
("given modularity"), and (700,529), which rests on a single run on another machine that this review did not reproduce.

### Ranked findings
1. FIX-FIRST, G1_FULL_SYSTEM.md:211: "(on each of their twelve curves `ellrank` gives an upper bound of at least 2)" is
   false. In census_lvalue_identify.json, 53 of the 228 curves have [1,1]. Replace with "(on none of their twelve curves,
   nor on any curve isogenous to one, does `ellrank` give an upper bound of 0)". The same false claim is in reports/w7.md
   ("Why global.py ...").
2. FIX, G1_FULL_SYSTEM.md:205 and :363, LEDGER.md:51, CHANGES.md:9: "L(E,1) > 0 certified without PARI" should read
   "... given modularity", because the series and the sign test w = +1 from F != 0 both use the functional equation of
   level N. Then "rank 0 by Kolyvagin (or Kato)". At :205, write "since isogenous curves have the same L-function and the
   same rank" (the L-function, not the rank, is what carries L(E,1) > 0 across). Say that (700,529) is one off-CI run on
   another machine with a patched copy of the code, not reproduced by a review. In LEDGER:51, give the review status.
3. FIX, G1_FULL_SYSTEM.md:211: "census_lvalue_identify.py reproduces those 19" overstates. The pairs are hard-coded from
   git history (identify.py:20-23); the script only finds their curves. Suggested text: "takes these 19 from the history
   of open_pairs_1000.txt (4ad9ae6, d185af6) and finds their curves; global.py re-run on the 268 closed census pairs leaves
   exactly these 19 open".
4. FIX, census_lvalue_big.py bigmerge: nothing structural catches a chunk that is duplicated or swapped. With K = 3,
   chunks 0 and 1 both hold 1842538 primes, so a copy passes the length assert and Hasse. In a simulation on (112,15),
   the guard refuses such a copy, and so would PARI's bigcheck. Fix: put K, J, the first and last prime, the count and a
   hash in each .bin.json and assert them; also assert sum(ap_stats) == len(ps).
5. NIT: census_lvalue_exact.py:64-106 is a byte-identical copy of :107-149. The --torsion code (torsion_step,
   torsion_all) is never called from __main__. In the negative test, "a_7 one more" changes a_7 at a bad prime
   (7 | N), so no good small prime is forged. merge() lets a later record for a pair overwrite an earlier one without
   warning. G1:417 --maxN 2e11 prints 13 series lines, 3 of which appear in the file as isogeny lines. CHANGES.md:9 has a
   stray space in "it does .".

### What the referee checked (all in /home/claude/rv_w7, python3 -I)
- (1) global.analyse re-run on all 268 closed census pairs: 249 closed by ellrank (165 + 84), and the 19 it leaves open
  are exactly OLD2 + NEW17. The git history (4ad9ae6 \ d185af6 = NEW17; OLD2 named in the doc at 4ad9ae6) and the tables
  agree.
- (2) Independent PARI check (ellfromeqn on y^2 = P(w) written out from m and n): N and w agree for (112,15), (127,97),
  (700,529), (943,576), and lfun(112,15) lies inside the enclosure. All 5 isogeny targets are in their source's
  ellisomat class.
- (3) The Velu formula and the c4/c6 isomorphism test with a square scale are correct and sufficient. The series run
  directly on the targets (487,217) and (527,336) gives enclosures identical to their sources'.
- (5) --pairs '112,15', --negative, --merge, and --pairs '260,69', '544,33' are identical to the committed files. The
  run with cypari2 blocked is also identical.
- (4) The (700,529) record is internally consistent (N, pi(Mmax) = 5527613 = sum of the stats, terms, the PARI lfun lies
  inside the enclosure). A re-run was started but stopped at about 45 minutes, past the hard stop, so it is NOT
  reproduced.

## Notes as written (17:44-18:08 UTC)

## Running notes (17:5x UTC)

(1) History. `git log -- open_pairs_1000.txt`: only 4ad9ae6 (174 open: 35+139) and d185af6 (157: 34+123). Set difference
4ad9ae6 \ d185af6 = exactly the 17 of NEW17 (16 trivial-point + (700,529)); nothing added. (112,15),(127,97) are named in
G1_FULL_SYSTEM.md at 4ad9ae6 line 72 as L-value closures; d185af6 line 72 says "19 ... 18 trivial + (700,529)" and the
commit body says "17 further closures by L-values". Tables: 167->183 (+16), 84->85 (+1): consistent with 2+17 = 19.
identify.txt: for all 19, "rank-0 by ellrank: none" on all 12 curves (incl. isogeny classes) -> global.py cannot close them.
Converse (no *other* pair closed outside ellrank): referee re-run of global.analyse on all 425-157 = 268 census pairs not in
the open list, in progress (rv_w7/census_all.py).

(3) Isogeny code (census_lvalue_exact.py:84-148): translation of a rational root to 0, Velu
y^2 = x(x^2 - 2a'x + a'^2 - 4b'), c4 = 16a^2-48b, c6 = -64a^3+288ab, test s = d6c4/(c6d4) a rational square with
d4 = s^2 c4, d6 = s^3 c6: all correct, and sufficient for a Q-isomorphism (c4c6 != 0 asserted). A chain of Velu 2-isogenies
from rational 2-torsion points is a Q-isogeny, so same L-function (V_l isomorphic as Galois modules) and same rank. The
target model is recomputed from (m,n) (merge asserts JSON model == rows[i][3]); the N-equality filter is harmless. HOLDS.
NIT: c46/isomorphic/two_isogenous/isogeny_chain are defined twice verbatim (lines ~84-117 and ~119-148).
NIT: torsion_step/torsion_all (--torsion) exist but are not wired into __main__ (dead code; output file never produced).

(5) Re-runs in /home/claude/rv_w7 (python3 -I): `--pairs '112,15'` -> line identical to census_lvalue_exact.txt:1 (stripped
of the PARI suffix), 15.7 s. `--negative` -> identical to census_lvalue_exact_negative.txt, 11/11 refused (78 s, 2 cores busy).
NIT: the forged "a_7 one more" is at a BAD prime (7 | N = 2^4*3*5*7*97*113*127; q5 picks 7 because 5 is bad, assuming 7 is
good). So no good small prime is forged; the test still refuses, but its label/intent is off (census_lvalue_exact.py negative(), q5).

**FIX-FIRST (false sentence).** G1_FULL_SYSTEM.md:211 "(on each of their twelve curves `ellrank` gives an upper bound of at
least 2)" is false, per the work's own census_lvalue_identify.json: 53 of the 228 curves have ellrank [1,1] (upper bound 1),
e.g. (112,15) w*f0*g4, w*f8*g4, w*f0*f8; (127,97) w*K. Distribution: (0,2) 59, (1,1) 53, (2,2) 47, (1,3) 43, (3,3) 10,
(0,4) 6, (1,5) 4, (4,4) 4, (2,4) 2. What is true: no curve has upper bound 0 (so global.py's `hi == 0` test never fires).
Replace by: "(on none of their twelve curves, nor on any curve isogenous to one, does `ellrank` give an upper bound of 0)".
The same false claim is in reports/w7.md "Why global.py ...": "all twelve curves have upper bound >= 2".

(2) Independent PARI (rv_w7/indep.py: ellfromeqn on y^2 = P(w) written from m, n as strings, not the repo's model):
(112,15) w*K N = 2338638960, w = +1, lfun 1.02590869756 (inside [1.025908697273, 1.025908697853]); (127,97) w*f0*f8 same N;
(700,529) w*f0*f8 N = 86826771232530, w = +1; (943,576) w*f0*f8 N = 497910195510, w = +1. All equal the claim.
The five isogeny targets are in the PARI isogeny class (ellisomat) of their sources, N equal: (352,135)->(487,217),
(575,48)->(623,527), (863,191)->(527,336), (260,69)->(329,191), (544,33)->(577,511). Partners check (m+n, m-n)/2^k. HOLDS.
`--merge census_lvalue_exact.jsonl` re-run: first 20 lines identical to census_lvalue_exact.txt, rest identical to _timing.txt.

Rizzo row (PARI-calibrated c6 sign, lball_exact.py:301) feeds only w_local; the certificate uses w from F != 0. No effect.

(1) converse, DONE (17:55): referee re-ran global.analyse (PARI 2.17.4, unchanged global.py) on all 268 census exceptions
(425 from gen_census1000.py's rule) not in open_pairs_1000.txt@d185af6. Result: 249 closed by ellrank rank 0 (165 with a
trivial point, 84 without = 183-18 and 85-1 of the tables), 19 left open, and those 19 are EXACTLY NEW17 + OLD2. No ALERT.
So the 19 are exactly the census's non-ellrank closures, and the table counts are consistent. HOLDS.
(m^2+n^2 square or twice a square: true for 18 of the 19; false only for (700,529), as the doc says.)

(4) (700,529) record (census_lvalue_big_700_529.jsonl) consistency, checked: N = 2*3*5*7*19*23*641*1201*1229 = 86826771232530
(and = my independent PARI N); pi(Mmax = 95696087) = 5527613 = pari.primes in the record = ap_stats sum
161+5527442+1+9 (161 = primes < 1000 minus the 7 bad ones below 1000; 9 bad primes); terms at y=1 55951440 = --list;
L_y2 inside L_y1; PARI lfun 4.362311899261 inside L_y1. So the record is internally consistent and the stats show three
complete, un-resumed chunks (a resumed chunk's .bin.json counts only the resumed part, so the sum would fall short).
Structural detection in bigmerge (census_lvalue_big.py:~70-80): a missing chunk -> open() fails; a truncated or
wrong-K chunk -> the length assert. A DUPLICATED or swapped chunk is NOT detected structurally: with K = 3 the lengths are
1842538, 1842538, 1842537, so chunk0.bin copied over chunk1.bin passes the length assert, and Hasse (x^2 <= 4p) passes too
(next prime is larger). It would be caught only by the analytic guard (very likely, cf. the a_p forgeries refused, but
not a proof) and by the PARI bigcheck. FIX: write (K, J, first prime, last prime, count, a hash) into each .bin.json and
assert them in bigmerge; assert sum(ap_stats) == len(ps).
Not re-runnable by the referee in the time: the record was made on another machine with a patched copy of lball_exact*.py
(build report: AST-md5 equal on the certificate path). A 2-chunk re-run was started here at 17:50 (rv_w7/b700_*.bin);
result below if it finishes.
(4, cont.) Duplicated-chunk simulation on the (112,15) curve (rv_w7/dupchunk.py, K = 3, chunk lengths 9551 = 9551): copying
chunk 0's a_p over chunk 1, or 1 over 2, is refused by the guard (enclosures disjoint, residual excludes 0; L_y1 becomes
2.9156 resp. -0.8702). So the guard does catch this failure mode in practice; the gap is only that bigmerge has no
structural check for it. The (700,529) re-run here was stopped at 17:55 (about 1 ms/prime on 2 cores: ~45 min, past the
hard stop); the (700,529) certificate is therefore NOT independently reproduced by this review.

(6) Doc text.
- G1_FULL_SYSTEM.md:211 — FIX-FIRST, see above (false "upper bound of at least 2"). Also "census_lvalue_identify.py
  reproduces those 19" overstates: the 19 are hard-coded from git history (OLD2 + NEW17, identify.py:20-23); the script
  finds their curves, it does not find the pairs. Suggested: "census_lvalue_identify.py takes these 19 pairs from the
  history of open_pairs_1000.txt (4ad9ae6, d185af6) and finds their closing curves; a re-run of global.py on the 268 closed
  census pairs leaves exactly these 19 open (review of 26 November 2026)."
- G1_FULL_SYSTEM.md:205 — FIX (wording, not substance): (a) "L(E,1) > 0 is certified without PARI" — the series and the
  root-number test themselves assume modularity (functional equation of level N, Carayol for N = conductor, Fricke
  eigenvalue); as written, modularity appears only in "Rank 0 then follows by modularity and Kolyvagin". Replace the last
  two sentences' frame by: "`L(E,1) > 0` is certified without PARI, given modularity (the series and the sign test use the
  functional equation) ... Rank 0 then follows by Kolyvagin (or Kato)." (b) "for the other 5 by an exact chain of two
  2-isogenies ..., since the rank is an isogeny invariant": what transfers L(E,1) > 0 is that isogenous curves have the
  same L-function; say "since isogenous curves have the same L-function and the same rank". (c) (700,529): add "one run,
  on another machine, not repeated in CI or by a review".
- G1_FULL_SYSTEM.md:363 (Limits) — FIX: add "given modularity" after "certified without PARI", and "(700,529): one run on
  another machine".
- G1_FULL_SYSTEM.md:417 (Reproduce) — NIT: `--maxN 2e11` sums 13 curves, 3 of which (487,217), (623,527), (527,336) are
  isogeny lines in census_lvalue_exact.txt, so its output is not comparable line-by-line with the file (the build report
  itself says to use --pairs). Use the CI's --pairs list, or say "13 lines; 3 of them appear in the file as isogeny lines".
- LEDGER.md:51 (Tier B, first row) — FIX: "rank 0 comes from L(E,1) != 0 (Kolyvagin), certified without PARI" -> "L(E,1) > 0
  certified without PARI given modularity; rank 0 by Kolyvagin"; add "(not reviewed before this review)" or this review's
  verdict, as the neighbouring tier rows do; add (700,529) single off-CI run.
- CHANGES.md:9 "Settled in this run" — mostly accurate (says torsion is PARI's, says not separately reviewed). FIX: same
  "given modularity" qualifier; NIT: stray space in "as that document says it does ." Same false ">= 2" claim is NOT in
  CHANGES (good). "with PARI's conductor, root number, a_n and L-value agreeing (a check)": true per the jsonl pari fields.
- OUTPUTS.md:166-167 — accurate (A for lines 1-11, B for 12-19; CI covers 1-11).

(3, cont.) Direct check of the isogeny step's conclusion: the referee ran the SERIES on two of the five isogeny targets:
`--pairs '487,217'` -> certified positive, L(E,1) in [3.341535197853, 3.341535198616] (byte-identical to (352,135));
`--pairs '527,336'` -> certified positive, [0.387584886368, 0.387584886673] (= (863,191)). Consistent with the chain.
The duplicated block census_lvalue_exact.py:64-106 and :107-149 is byte-identical (NIT: delete one).
(5, cont.) PARI-free claim tested: `--pairs '112,15;127,97'` run with `sys.modules['cypari2'] = None` (import blocked)
-> output identical to census_lvalue_exact.txt lines 1-2 (PARI suffix stripped). cypari2 is imported only lazily in
check functions (lball_exact53.get_pari, lball_exact.py:578ff), never on the certificate path.
NIT: merge() (census_lvalue_exact.py, `recs[i] = r`) silently lets a later record for the same pair overwrite an earlier one;
and the txt it writes re-renders stored booleans (`positive`, `guard`) — for lines 12-19 it is a record, not a re-check
(OUTPUTS.md already classes these B, correctly).

(5, cont.) `--pairs '260,69'` and `--pairs '544,33'` (no --pari) re-run here: identical to census_lvalue_exact.txt lines
12-13. These are the sources of the isogeny lines for (329,191) and (577,511), and are not in CI (OUTPUTS class B), so this
is their first reproduction. Not reproduced by this review: lines 3-11 (CI's descents job covers them) and (700,529).
