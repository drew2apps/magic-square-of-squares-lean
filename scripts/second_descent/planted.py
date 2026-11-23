#!/usr/bin/env python3
"""Planted points: the equation step of the second descent on real, non-degenerate points.

For any coprime m > n and any rational t0 with F0(t0) F8(t0) a square (a rational point of the classical curve
Y^2 = (m^2 + n^2 t^2)(n^2 + m^2 t^2)), the number g = Q_b(-i t0), Q_b = (m - n s)(n + m s)(1 - s^2), has norm a square, so
g = d * square in K = Q(i) for a rational d (Hilbert 90).  Then (s0, Y), s0 = -i t0, is a K-point of Y^2 = d Q_b(s): the
situation of SECOND_DESCENT.md, sections 2 and 3, with a point that is known to exist and is not degenerate.  The first
descent is not used.  If E(Q) or E^(-1)(Q) has rank 0 for this d, the solver must return t0.

usage: planted.py M U [k]     ratios m <= M, t0 = u/v with u, v <= U; with k > 1, a negative control: d is replaced by k d,
                              the point is then not on the curve, and the solver should miss it whenever a summand has rank 0"""
import sys
from math import gcd, isqrt
from fractions import Fraction as Fr
from sage.all__sagemath_schemes import *
from second_descent import K, I, sq, try_curve

def rational_part(g):
    """d squarefree with g/d a square in K, for g in K* whose norm is a rational square."""
    r = QQ(g.norm()).sqrt(); h = g + r
    if h == 0: h = g - r
    q = g/h**2; assert q[1] == 0
    d = QQ(q[0]).squarefree_part(); d = abs(d); assert (g/d).is_square()
    return d

def main(M, U, wrong=1):
    pts = 0; norank0 = 0; tested = 0; found = 0; seen = {}
    for m in range(2, M + 1):
        for n in range(1, m):
            if gcd(m, n) != 1: continue
            Q = (m - n*sq)*(n + m*sq)*(1 - sq**2)
            for v in range(1, U + 1):
                for u in range(1, U + 1):
                    if gcd(u, v) != 1 or u == v: continue
                    N = (m*m*v*v + n*n*u*u)*(n*n*v*v + m*m*u*u)
                    if isqrt(N)**2 != N: continue
                    t0 = QQ(u)/v; pts += 1
                    d = (rational_part(K(Q(-I*t0)))*wrong).squarefree_part()
                    if (m, n, d) not in seen: seen[(m, n, d)] = try_curve(m, n, '1356', Q, QQ(m)/n, d_given=d)
                    o = seen[(m, n, d)]
                    if o['status'] == 'no summand of rank 0': norank0 += 1; continue
                    tested += 1
                    got = set(o['rational_t'])
                    ok = str(Fr(u, v)) in got and str(Fr(-u, v)) in got
                    found += ok
                    print('m:n = %d:%d  t0 = %s  d = %s  rank bounds %s  torsion over K %d  rational t returned: %s  %s' % (m, n, t0, d, o['rank_bounds'], o['torsion_K'], ' '.join(sorted(got)), 'found' if ok else 'MISSED'), flush=True)
    if wrong != 1: print('NEGATIVE CONTROL: d multiplied by %d, so the point is not on the curve' % wrong)
    print('points (m, n, t0) with F0 F8 a square, t0 not 0, 1: %d on %d curves;  no summand of rank 0: %d;  a summand of rank 0: %d, of these t0 and -t0 returned: %d, missed: %d' % (pts, len(seen), norank0, tested, found, tested - found))
    assert wrong != 1 or tested == found

if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 1)
