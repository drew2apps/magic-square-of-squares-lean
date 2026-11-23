#!/usr/bin/env python3
"""The Lemma of GAUSS_DESCENT.md and the local classes, tested on real points.

For a ratio m:n and each of the five forms separately, take every coprime (u, v) with |u|, |v| <= U at which that one form
satisfies its condition (N(l(u, v)) in c Q^2: a conic, so there are many such points).  Then l(u, v) = delta * beta^2 in Z[i].
(1) Lemma: the squarefree part delta must be, up to a square, in the set D_j that gauss_descent.py uses.  Tested with a
    square test in Z[i] that uses integer square roots only (no class labels, no factorisation): exactly one element d of D_j
    must make l(u, v) * d a square.
(2) Local classes: for that d, and every prime p of a list, the class of l(u, v) at the places above p must equal the class
    of d, wherever both are determined.  This tests the labels (at 1+i: valuation parity and the unit modulo pi^5) on numbers
    that are known to differ by a square.

usage: real_points.py U m,n [m,n ...]"""
import sys
from math import gcd, isqrt
import gauss_descent as gd

def is_gsquare(z):
    """z = (a, b) in Z[i], z != 0, is a square in Z[i] (equivalently in Q(i)): w = x + iy with x^2 = (r + a)/2, y^2 = (r - a)/2, r = |z|."""
    a, b = z; N = a*a + b*b; r = isqrt(N)
    if r*r != N or (r + a) % 2: return False
    x2, y2 = (r + a)//2, (r - a)//2; x, y = isqrt(x2), isqrt(y2)
    return x*x == x2 and y*y == y2 and abs(2*x*y) == abs(b)

def is_sq(n): return n >= 0 and isqrt(n)**2 == n

PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]

def run(m, n, U):
    S = m*m + n*n; consts = [1, 1, 2*S, 2, 2]; D = gd.dsets(m, n); F = gd.forms(m, n)
    locs = [gd.Local(p) for p in PRIMES]; cD = [[[L.cls(d) for d in Dj] for Dj in D] for L in locs]
    pts = [0]*5; bad_lemma = 0; cmp = 0; bad_cls = 0; cmp2 = 0
    for u in range(0, U + 1):
        for v in range(-U, U + 1):
            if gcd(u, v) != 1 or (u == 0 and v < 0): continue
            for j in range(5):
                z = gd.value(F[j], u, v)
                if z == (0, 0) or not is_sq(consts[j]*gd.gnorm(z)): continue
                pts[j] += 1
                hits = [i for i, d in enumerate(D[j]) if is_gsquare(gd.gmul(z, d))]
                if len(hits) != 1: bad_lemma += 1; print('LEMMA FAILS', (m, n), j, (u, v), hits); continue
                for L, cd in zip(locs, cD):
                    cz = L.cls(z); cdl = cd[j][hits[0]]
                    for a, b in zip(cz, cdl):
                        if a is None or b is None: continue
                        cmp += 1; cmp2 += L.p == 2
                        if a != b: bad_cls += 1; print('CLASS DIFFERS', (m, n), j, (u, v), L.p, cz, cdl)
    return pts, bad_lemma, cmp, cmp2, bad_cls

if __name__ == '__main__':
    U = int(sys.argv[1]); tot = [0]*5; bl = 0; c = 0; c2 = 0; bc = 0; k = 0
    for s in sys.argv[2:]:
        m, n = (int(x) for x in s.split(',')); pts, a, b, b2, d = run(m, n, U)
        tot = [x + y for x, y in zip(tot, pts)]; bl += a; c += b; c2 += b2; bc += d; k += 1
    print('%d ratios, |u|, |v| <= %d: points on the single conditions (F0, F8, G4, GL, GU): %s, in all %d' % (k, U, tot, sum(tot)))
    print('Lemma: squarefree part outside the set D_j, or not unique: %d' % bl)
    print('local classes compared at the primes up to 37: %d (at 1+i: %d); different: %d' % (c, c2, bc))
    assert bl == 0 and bc == 0
