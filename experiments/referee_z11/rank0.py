#!/usr/bin/env python3
"""z11: the five ratios of t10's 201 with a pair curve of rank 0 have no four-form point at all.
Written by an AI model (review agent z11); not read by a person.  Experiment, not counted; tier D.

For (5623,1017) E_08 and (3320,2303), (9191,7688), (9532,4445), (14885,6076) E_LU:
  1. dim Sel^2(E) = 2 = dim E(Q)[2] by the complete 2-descent (z7's my2desc.Curve on z11_pair's integral model; t10's
     rank_second gives the same, see recount.txt), so rank E(Q) = 0;
  2. #E(F_p), counted exactly by listing x mod p, for the good primes 3 <= p <= 199: their gcd bounds #E(Q)_tors
     (reduction is injective on torsion at good p >= 3), and it is 8;
  3. C(Q) contains 8 explicit points, over t = 0, 1, -1, infinity (G(t), resp. the leading coefficient, a non-zero
     square, so two points over each); since C(Q) and E(Q)
     are in bijection (Washington's map, wash_check.py), C(Q) is exactly these 8;
  4. at t = 0, 1, -1, infinity the four forms F_0, F_8, F_L, F_U are not all squares (the ratio is not degenerate).
So no rational t makes the four forms squares: these ratios, and the tuples t10 could not pairing-kill, are excluded
by rank 0 of one pair curve (computed; tier D; m > 1000, so nothing for the census).
PARI: factorisation of integers (primes checked) and of the cubic (roots checked exactly); no elltors, ellap or ellrank.
    python3 -I rank0.py      (prints rank0.txt)"""
import sys, os
from math import gcd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _z11common import Out
from z11_pair import PairCurve, Pairing, Fcoef, Fr, issq, pari
say = Out('rank0.txt')
say('rank0.py (review z11; written by an AI model): pair curves of dim Sel^2 = 2 among t10\'s 201 ratios')
def count(C, p):
    A2, A4, A6 = [c % p for c in C.cub]; n = 1
    sq = [0] * p
    for y in range(p): sq[y * y % p] += 1
    for x in range(p): n += sq[(x * x * x + A2 * x * x + A4 * x + A6) % p]
    return n
def G_at(m, n, jk, t):
    if t == 'inf':
        return Fcoef(m, n, jk[0])[2] * Fcoef(m, n, jk[1])[2]
    return sum(c * Fr(t) ** i for i, c in enumerate(Fcoef(m, n, jk[0]))) * sum(c * Fr(t) ** i for i, c in enumerate(Fcoef(m, n, jk[1])))
def F_at(m, n, j, t):
    if t == 'inf': return Fr(Fcoef(m, n, j)[2])
    return sum(c * Fr(t) ** i for i, c in enumerate(Fcoef(m, n, j)))
allok = True
for m, n, jk in [(5623, 1017, '08'), (3320, 2303, 'LU'), (9191, 7688, 'LU'), (9532, 4445, 'LU'), (14885, 6076, 'LU')]:
    C = PairCurve(m, n, jk); PR = Pairing(C.es)
    e1, e2, e3 = C.es; disc = ((e1 - e2) * (e1 - e3) * (e2 - e3)) ** 2
    g = 0; used = 0
    for p in pari.primes([3, 199]):
        p = int(p)
        if disc % p == 0: continue
        g = gcd(g, count(C, p)); used += 1
    pts = []
    for t in (0, 1, -1, 'inf'):
        v = G_at(m, n, jk, t); assert v != 0 and issq(v)
        pts += [(t, +1), (t, -1)]                      # y = +-sqrt(G(t)) (at infinity: +-sqrt of the leading coefficient)
    four = {str(t): [issq(F_at(m, n, j, t)) for j in '08LU'] for t in (0, 1, -1, 'inf')}
    none = not any(all(v) for v in four.values())
    ok = PR.dim == 2 and g == 8 and len(set(pts)) == 8 and none
    allok &= ok
    say('(%d,%d) E_%s: dim Sel^2 = %d, so rank 0; gcd of #E(F_p) over %d good primes 3 <= p <= 199: %d; C(Q) = the 8 points'
        ' over t = 0, 1, -1, inf (G(t) a non-zero square at each); four forms squares there (F_0, F_8, F_L, F_U): %s; a four-form'
        ' point: %s' % (m, n, jk, PR.dim, used, g, ' '.join('%s:%s' % (t, ''.join('y' if x else 'n' for x in v)) for t, v in four.items()),
                        'NONE' if none else 'YES'))
say('ALL FIVE RATIOS HAVE NO FOUR-FORM POINT (rank 0 of one pair curve)' if allok else 'A CHECK FAILED')
say.finish()
