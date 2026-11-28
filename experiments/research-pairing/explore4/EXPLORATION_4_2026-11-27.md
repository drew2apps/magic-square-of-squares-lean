# Exploration 4: the pairing on all 937 ratios, the surviving shape, root numbers (2026-11-27)

> **Read `NOTES.md` first.** This is the session's note, committed unchanged apart from this line; where they differ, `NOTES.md` (after review) holds.


Tier D. Written by an AI model (Claude), not read by a person. Nothing here is proved unless a line says so.
No counts change. Nothing was published. Scripts and outputs: `explore4.tgz`.

## 1. The pairing run on every ratio the bad primes leave open (m <= 15000)

Population: the repo's lists of non-degenerate ratios that keep a class after the first stage:
997 with m <= 5000 and 8,739 with 5000 < m <= 15000 (9,736 in all).

| stage | m <= 5000 | 5000 < m <= 15000 | all |
|---|---|---|---|
| ratios | 997 | 8,739 | 9,736 |
| bad primes alone leave no tuple (repo lists) | 896 | 7,903 | 8,799 |
| bad primes leave tuples | 101 | 836 | 937 |
| ... and every tuple is killed by the pairing on a pair curve | 100 | 832 | 932 |
| ... a tuple survives all six pair curves | 1 | 4 | 5 |

`pipeline.py` takes the residual tuples from the repo's own code (`tuple_descent.magic`, `exact_local.allowed_at`,
`exact_local.join`) and applies the general Cassels pairing (`cassels.py`, `sel.py` of exploration 3) on the six
pair curves. The pilot had run 201 of these 937 ratios.

Checks:
- On the pilot's 201 ratios: 1,766 tuples, 10,596 (tuple, curve) decisions. Mine and the pilot's agree on all:
  8,930 killed by both, 1,666 killed by neither.
- All 6,904 residual tuples map into Sel^2 of all six curves (0 outside).
- Killed per curve: E_08 4,580; E_LU 3,122; each mixed curve 6,587. Killed on E_08 or E_LU alone: 5,864 of 6,904.

The five ratios with a surviving tuple (`surv.txt`):

| ratio | m^2+n^2 | surviving tuples | rank bound E_08 / E_LU (Sel^2 minus pairing; PARI agrees) | repo's kill |
|---|---|---|---|---|
| (3320,2303) | 2621 * 6229 | (2,1,1,635399), (2,1,1270798,2) | <=1 / 0 | good primes 17, 29 |
| (5623,1017) | 2 * 2621 * 6229 | (2,1270798,1,2), (635399,1,1,2) | 0 / <=1 | (partner) |
| (9113,3200) | 17 * 5487457 | (1,2,511,1759), (1,2,3518,1022) | <=2 / <=1 | good prime 37 |
| (12313,5913) | 2 * 17 * 5487457 | (1022,3518,2,1), (1759,511,2,1) | <=1 / <=2 | good prime 37 |
| (14903,9696) | 5^2 * 12644473 | (1,2,10414,49198), (1,2,24599,5207) | <=2 / 0 | good prime 37 |

The first two are the pilot's exceptions. The other three are new; the pilot had not sampled them.
(9113,3200) and (12313,5913) are partners under (m,n) -> (m+n, m-n). The partner of (14903,9696) is
(24599,5207), outside the range. The census already closes all five with good primes; this run only asks
whether the pairing can replace the good primes.

Note: at (9113,3200) and (12313,5913) neither E_08 nor E_LU has rank 0 (root numbers +1/-1, bounds 2 and 1).
So "a pair curve of rank 0" does not describe the exceptions. The description that fits all five is section 2.

## 2. The shape of a surviving tuple (answer to the question in the task, on this data)

For a tuple (d0,d8,dL,dU) put a = sqf(mn d0 d8) (its class on E_08 is (1,a,a)) and
b = sqf((m^2-n^2) dL dU) (its class on E_LU, after the 4-torsion shift, is (1,b,b)).

Tested (`analyse2.py`, 6,904 tuples at 937 ratios): **all 10 surviving tuples have
(a,b) = (sqf(2mn), 1) or (a,b) = (1, sqf(2(m^2-n^2))).** No tuple of any other shape survives.

These are exactly the E_08 and E_LU classes of the degenerate solutions (`ghost.py`, 13 degenerate ratios):
at a square-type ratio the solutions t = +1, -1 have (a,b) = (sqf(2mn), 1); at a twice-square-type ratio the
solutions t = 0, infinity have (a,b) = (1, sqf(2(m^2-n^2))). On E_08 the class (1,2mn,2mn) is the class of the
rational 4-torsion point times (C,1,C), C = m^2+n^2; it is a torsion class exactly when C is a square.

So a surviving tuple is a "ghost" of a degenerate solution: it has the E_08 and E_LU classes that a degenerate
solution would have, at a ratio where m^2+n^2 is not a square or twice a square. The (dL,dU) split is not
always the degenerate one (at (9113,3200) a factor 7 sits on the other side), so the mixed-curve classes differ.

Counts: 2,042 residual tuples at 499 ratios have the ghost shape; 2,032 are killed (1,556 on E_08 or E_LU,
476 only on mixed curves; see `analyse2.txt`), 10 survive.

Why this matters for (R): (R) says solutions exist only at degenerate ratios. The data say the only classes that
get past bad primes and all six pairings are the classes of the degenerate solutions themselves. A proof of
"a tuple with trivial pairing on all six curves has the ghost shape" would reduce stage 2 to one explicit class
per ratio. This is tested on 937 ratios, not proved, and I have no proof route yet.

## 3. The bad-prime laws are not linear (a negative result)

I had proposed to write the bad-prime laws as linear conditions. Tested (`analyse2.py`): the residual set of a
ratio is a coset of an F_2-subspace at 682 of 929 ratios (8 larger sets not tested) and is **not** a coset at 247.
Sizes 6, 10, 12, 14, 18 occur. Reason: the local image of the genus-5 curve in the tuple group is a set, not a
subgroup. Linearity holds only after mapping to a pair curve (Sel^2 is a group). So the full-tuple stage cannot be
a single matrix computation; the pair-curve stage can. Every residual set has even size.

Failure rate asked for by the panel, on these lists: after bad primes 101/997 = 10.1% and 836/8,739 = 9.6%;
after bad primes and the six pairings 1/997 = 0.10% and 4/8,739 = 0.046%.

## 4. Root numbers: no correlation (a negative result)

`rootno.py`, `rootno2.py`. The root numbers (w(E_08), w(E_LU)) take each of the four sign patterns about a
quarter of the time in every population tested: class-keeping ratios m <= 5000 (997), a sample of 2,500 with
5000 < m <= 15000, 2,500 random ratios, 2,500 random ratios with a 2-adic point, and the 937 ratios with residual
tuples (230 / 216 / 242 / 249). w(E_mix) is +1 at 493 and -1 at 444 of the 937 (1.6 standard deviations).
The local conditions for a solution do not select a parity. Parity gives no handle on (R).

## 5. Mixed curves: class constants, and the first step of a proof of the Redei form

Class constants (`kappa.py`, all 6,776 pilot lines, 1,694 per curve): the class of a tuple on a mixed curve is
(2, kappa d_j d_k) at (e1,e2), with kappa = sqf(mn(m^2-n^2)) for 0L and 8U and 2 sqf(mn(m^2-n^2)) for 0U and 8L.

Lemma (proved; one line of algebra, checked by `tangent_norm.py`). On the conic u1^2 - d u3^2 = c T^2 with a
point Q = (q1,q3,qT), the tangent form L(P) = q1 u1 - d q3 u3 - c qT T satisfies, for P on the conic,

    L(P) = - N( qT (u1 + sqrt(d) u3) - T (q1 + sqrt(d) q3) ) / (2 qT T),     N = norm from Q(sqrt d).

Consequence (given Cassels' theorem and the product formula; the constant -2 qT drops out): for classes
a = (1,d,d), b = (1,beta,beta) on E_mix,

    <a,b> = prod_v ( N(w_{2,v}) N(w_{3,v}), beta )_v ,

where w_{2,v}, w_{3,v} are the differences "local point minus global point" on the two norm conics
u1^2 - d u3^2 = 2 R+ and u1^2 - d u2^2 = 2 R-. The product (u1 + sqrt(d) u3)(u1 + sqrt(d) u2) has norm 4W,
W = R+ R-. This is the shape of the Redei symbol [d, W, beta]. The form tested in exploration 3
(the symbol [d, beta, W]) would then follow from Redei reciprocity (Stevenhagen, arXiv:1806.06250).
Not done: turning the product of Hilbert symbols of norms into the Artin symbol, the places 2 and infinity,
and the minimal-ramification condition. Classes with first component not 1: not started.

## 6. What I did not do

- The point search on the 4-coverings for (627,364).
- Sel^2 as a full matrix for classes with b_r not 1.
- The pairing on the three genus-3 quotient curves E1, E2, E3 at the five ghost ratios (the class map from tuples
  to those curves is not worked out). This is the natural next test: does any quotient kill a ghost tuple?
- A 4-descent on E_08 at the ghost ratios.

## Files in explore4.tgz

pipeline.py, analyse.py, analyse2.py, surv.py, ghost.py, kappa.py, rootno.py, rootno2.py, tangent_norm.py,
their .txt outputs, out_997.jsonl, out_15000_0.jsonl, out_15000_1.jsonl, survivors.json, and copies of cassels.py,
sel.py, and the three ratio lists used as input (g2_*.txt, copied from the repo's stage-2 lists).
pipeline.py needs the repo (MSS_REPO); analyse.py needs the pilot files for the comparison.
