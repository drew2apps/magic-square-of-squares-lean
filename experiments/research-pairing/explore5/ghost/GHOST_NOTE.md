# Ghost tuples: prediction test, sampling beyond m = 15000, laws, and the centre form (2026-11-27)

> **Read `../NOTES.md` first.** This is the session's note, committed unchanged apart from this line; where they differ, `../NOTES.md` (after review) holds.


Written by an AI model (Claude). Not read by a person. "Tested" = experiment with the counts given; "proved" is used
nowhere below except where a derivation is written out. Everything is in /tmp/ghost; nothing was published.
About 45 minutes of work; several items are unfinished and listed at the end.

## 1. Prediction test and attempt to break "every survivor has the ghost shape"

Prediction (t1a.py, t1a.txt): the partner (24599,5207) of (14903,9696) has 8 residual tuples; 4 have the ghost shape;
exactly 2 survive all six pairings, (5207,24599,2,1) and (49198,10414,2,1), both of ghost shape (S2). Confirmed.

Sampling beyond m = 15000 (sampler.py; residual() of pipeline.py called directly; t1_count.txt):

| sample | ratios tried | with residual tuples | tuples | ghost-shape tuples | survivors | survivors not of ghost shape |
|---|---|---|---|---|---|---|
| biased: m+n, m-n built from primes = +-1 mod 8, half replaced by the partner; m up to about 150000 | at least 5,892 (*) | 269 | 1,096 | 622 | 6 (3 ratios) | 0 |
| random coprime, 15000 < m < 60000, with the mod-8 pattern of the stage-1 lists | 1,836 | 14 | 52 | 30 | 0 | 0 |
| m^2+n^2 with at most 2 odd-exponent odd primes, same range and mod-8 pattern | 2,244 | 38 | 218 | 56 | 0 | 0 |

(*) the biased run was stopped by hand; the number is the "tried" counter at the last ratio written. No error lines, no tuple outside Sel^2.

New ratios with surviving tuples (all of ghost shape):
- (90384,58537), m^2+n^2 = 5^2 * 89 * 109 * 137 * 349: (2,1,31847,148921), (2,1,297842,63694)   (S1)
- (126127,16889), m^2+n^2 = 2 * 5^2 * 173 * 313 * 5981: (2,4260317806,1,2), (2130158903,1,1,2)   (S2)
- (115513,18217), m^2+n^2 = 2 * 53 * 129010493: (36434,231026,2,1), (115513,18217,2,1)   (S2)

So the claim was not broken: 937 + 321 + 1 (the partner) = 1,259 ratios with residual tuples and 6,904 + 1,366 + 8 =
8,278 tuples now tested (samples not de-duplicated against each other), 10 + 6 + 2 = 18 survivors at 9 ratios, all of
ghost shape. Note m^2+n^2 can have five odd prime factors at a ghost
ratio, so "few prime factors" is not part of the description.

## 2. Laws

(a) Residual. Tested on the 9,736 stage-1 ratios (law_a.py): an S1-shape tuple (d0 d8 ~ 2, dL dU ~ m^2-n^2) is
residual at 333 ratios, an S2-shape tuple at 166, never both. At all 333 (resp. 166) every odd prime dividing m^2-n^2
(resp. mn) to an odd power is = +-1 mod 8: 333 of 333 and 166 of 166. Not sufficient: 1,554 (resp. 1,392) stage-1
ratios satisfy the condition and have no such residual tuple. Primes = 3 mod 8 do occur to an even power (94 and 46
ratios), e.g. 3^2 | m-n at (3320,2303). This condition was used to bias the sampler; hit rate about 4.6% (269 of at least 5,892) against 0.8% (14 of 1,836).
All 333 S1 ratios have m, n of opposite parity.

(b) Sel^2 membership on E_08. Tested on 3,000 random coprime ratios, m < 4000, C not a square (law_b.py):
- m, n of opposite parity: (1,2mn,2mn) is in Sel^2(E_08) if and only if no prime = 5 mod 8 divides m^2-n^2: 1,994 of 1,994.
- m, n both odd: never in Sel^2: 1,006 of 1,006 (use the partner ratio).
Outline of the reason at odd primes (not written out in full, so count it as tested): the class is that of the
4-torsion point times (C,1,C). At p | mn and p | C the class (C,1,C) lies in the local image automatically (C is a
square mod p | mn; at p | C the local image contains all classes (*,1,*)). At p | m^2-n^2 the reduction is
multiplicative at the node x = 0, split iff p = 1 mod 4; non-split: every unramified class is in the image; split:
(C,1,C) is in the image iff C = 2n^2 is a square mod p, i.e. p = 1 mod 8. The 2-adic part (parity) is only tested.
So (b) is strictly weaker than (a): (b) allows p = 3 mod 8, (a) does not (to an odd power).

(c) Unkilled on E_08: no law found (law_c.py, 333 S1 ratios, 250 killed on E_08, 83 not). Predicates tried, none fits:
"C is a 4th power mod every odd p | m^2-n^2" (unkilled 26 true / 57 false; killed 103 / 147); same for p | mn
(25/58; 62/188); "2mn is a 4th power mod every odd-exponent prime of C" (43/40; 136/114); "m^2-n^2 a square mod
those primes" (32/51; 151/99). "2mn is a square mod every odd-exponent prime of C" holds at all 333 (it is implied).
The number of odd-exponent primes of C does not separate either. On E_LU no S1-shape tuple is ever killed (0 of 333;
its class there is a torsion class, see 3c). Mixed curves: not attempted.

## 3. More quotients

(b) THE CENTRE FORM (main finding). The repository's residual test uses four forms only (tuple_descent.magic returns
F0, F8, GL, GU; exact_local and g2_killsets take those). "F4 is a square" is not imposed and is not implied locally.
I added gamma_4 = (1+i)(m+in)(u+iv), of norm F4 = 2(m^2+n^2)(u^2+v^2), as a fifth linear form with class d4 supported
on the primes of 2(m^2+n^2), and reran the same local code at the same bad primes (f4_form.py, f4_all.py).
Sanity: at the degenerate ratios (4,3), (12,5), (7,1), (17,7), (21,20), (41,1) exactly the two known solutions remain.
This is valid only if a solution must have F4 square (I read the task that way; F4 is the centre).

| population | ratios | keep a residual 5-tuple | 4-tuples | 4-tuples that extend | pairing survivors | survivors that extend |
|---|---|---|---|---|---|---|
| the 937 ratios, m <= 15000 | 937 | 211 | 6,904 | 1,370 | 10 (5 ratios) | 2 (1 ratio) |
| biased sample, 269 ratios | 269 | 46 | 1,096 | 118 | 6 (3 ratios) | 0 |

At (3320,2303), (5623,1017), (9113,3200), (12313,5913) and the three new ghost ratios, the centre form at the bad
primes alone leaves no 5-tuple at all: no good prime and no pairing needed. Only (14903,9696) and its partner
(24599,5207) keep ghost tuples after six pairings and the centre form: (1,2,10414,49198,d4) with d4 in {2, 2p} and
(1,2,24599,5207,d4) with d4 in {1, p}, p = 12644473; partner: (5207,24599,2,1,d4), d4 in {2, 2p}, and
(49198,10414,2,1,d4), d4 in {1, p}. There m^2+n^2 = 5^2 * p.
Not done: the pairing on the four new pair curves y^2 = F_j F4. They all have a rational point (t = +-1 for j = 0, 8;
t = 0, infinity for j = L, U), and by the same argument as for the six old curves the class of a 5-tuple should be
(1, kappa_j d_j d_4) with a constant kappa_j to be fitted (the function Tr(gamma_j gamma_4) + 2y is d_j d_4 times a
square). This is the obvious next test on (14903,9696).

(c) Ranks (ranks.py, ranks2.py, ranks2_clean.txt). At each of the six ghost ratios, on one of E_08, E_LU the ghost
class lies in the image of the rational torsion Z/4 x Z/2 (so nothing on that curve can kill it), and on the other it
is a non-torsion class in the kernel of the Cassels-Tate pairing:

| ratio | curve where the class is not torsion | dim Sel^2 | dim kernel | rank bound | ellrank | root number |
|---|---|---|---|---|---|---|
| (3320,2303) | E_08 | 5 | 3 | <= 1 | [1,1] | -1 |
| (5623,1017) | E_LU | 5 | 3 | <= 1 | [1,1] | -1 |
| (9113,3200) | E_08 | 4 | 4 | <= 2 | [0,2] | +1 |
| (12313,5913) | E_LU | 4 | 4 | <= 2 | [0,2] | +1 |
| (14903,9696) | E_08 | 4 | 4 | <= 2 | [0,2] | +1 |
| (24599,5207) | E_LU | 4 | 4 | <= 2 | [0,2] | +1 |

At (3320,2303)/(5623,1017): the kernel is torsion plus one class, and the ghost class is that class. If the rank is 1
(PARI's [1,1]; no point was found, so I take this as resting on parity, i.e. on finiteness of the 2-part of Sha), the
ghost class is the image of a point of infinite order and no argument on that curve alone can kill it. The point was
not found: no rational point of height <= 100000 on two quartic models of the ghost 2-covering. At the other four the
rank is 0 or 2 and undecided; same point search, nothing found. So: not decided unconditionally at any ratio.
PARI ellheegner / higher effort / 4-descent: not tried.

(a) E1, E2, E3 of the genus-3 curve: not done (no class map worked out).

(d) Remarks. At t = 1 the four norms are C, C, 4(m-n)^2, 4(m+n)^2 and at t = 0 they are m^2, n^2, 2C, 2C, so the
degenerate "solution" fails only because C (resp. 2C) is not a square; the ghost class is the class of that point,
and law (b) is the statement that C is a square at enough places. The centre form sees the same defect a second time
(F4 = 4C at t = 1), which is presumably why it removes most ghosts. Guess, untested: the ghosts that survive F4 are
those where C is a square times one prime and the local conditions cannot see that prime (d4 and d4 * p both pass).

## Not done / caveats
- 3(a); pairing on y^2 = F_j F4; law (c) on E_08 and on mixed curves; 4-descent; unconditional ranks.
- The biased sampler's exact number of ratios tried is not recorded (at least 5,892).
- The F4 test reuses the repository's local code with five forms; I did not run its selftest/oracle on five forms.
- The statement "F4 must be a square" is taken from the task text, not checked against the repository's documents.

## Files
t1a.py/.txt, sampler.py, samp_biased.jsonl, samp_random.jsonl, samp_fewp.jsonl, t1_count.py/.txt, law_a.py,
ghostdata.json, law_b.py/.txt, law_b.json, law_c.py/.txt, law_c_rows.json, f4_form.py/.txt, f4_all.py/.txt,
f4_all.jsonl, ranks.py/.txt, ranks2.py, ranks2.txt, ranks2_clean.txt, feat.py, time1.py.
