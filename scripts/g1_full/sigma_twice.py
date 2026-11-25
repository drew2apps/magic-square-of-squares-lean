#!/usr/bin/env python3
"""Why the curves of F0 F8 GL never close a ratio with m^2 + n^2 twice a square (G1_FULL_SYSTEM.md, Result 4).
Exact rational arithmetic, no Sage, no PARI (functions of sigma_exact.py).

For every ratio, t = 1 gives F0 F8 GL = S^2 (2(m - n))^2, a square, so E+ has the point P with z = 2, and t = -1 gives
the point with z = -2 on E- (G1Sigma.lean: sigma_plus, sigma_minus_U at t = 1).  When 2S is a square, t = 0 gives the
point X = 0 on both curves: the constant term of the model is A^2 times the leading coefficient 2 m^2 n^2 S of q.
This script proves, for each of the 127 coprime ratios m <= 1000 with m^2 + n^2 twice a square and both curves, that
the two points are independent, so the rank is at least 2 and the argument of Demjanenko and Manin with one generator
cannot apply: the determinant hhat(P) hhat(Q) - <P, Q>^2 of the height pairing is enclosed in rational arithmetic
(each canonical height from 2^k times the point, with h - hhat <= B and hhat - h <= U as in sigma_exact.py) and its
lower end is positive.  For the ratios with m^2 + n^2 a square it checks that the point of t = +-1 has infinite order
(its canonical height has a positive lower end); that point is the degenerate solution itself."""
import os, sys
from fractions import Fraction as Fr
from math import gcd, isqrt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sigma_exact import model, Curve, hnaive, ln_bounds, lnplus_hi, LN2


def bounds(C):
    j = Fr(C.c4**3, C.D); jsize = max(abs(j.numerator), j.denominator)
    common = ln_bounds(abs(C.D))[1]/6 + lnplus_hi(abs(j))/6 + lnplus_hi(Fr(abs(C.b2), 12)) + (LN2[1] if C.b2 else 0)
    assert jsize >= 32
    return common + ln_bounds(jsize)[1]/12 + Fr(1922, 1000), common + Fr(214, 100)


def hhat(C, P, B, U, k=6):
    x = P[0]
    for _ in range(k): x = C.dup_x(x)
    h = hnaive(x)
    return (h[0] - B) / 4**k, (h[1] + U) / 4**k


def main():
    deg = [(m, n) for m in range(2, 1001) for n in range(1, m) if gcd(m, n) == 1 and
           ((isqrt(m*m + n*n)**2 == m*m + n*n) or ((m*m + n*n) % 2 == 0 and isqrt((m*m + n*n)//2)**2 == (m*m + n*n)//2))]
    tw = [p for p in deg if isqrt(p[0]**2 + p[1]**2)**2 != p[0]**2 + p[1]**2]
    least = None; inf_order = 0; sq_curves = 0
    for (m, n) in deg:
        for s in (2, -2):
            M = model(m, n, s); C = Curve(M['a2'], M['a4'], M['a6']); A, z0 = M['A'], M['z0']
            B, U = bounds(C)
            zp = 2 if s > 0 else -2
            X = Fr(A, zp - z0); y2 = X**3 + C.a2*X*X + C.a4*X + C.a6
            assert y2 > 0 and isqrt(y2.numerator)**2 == y2.numerator and isqrt(y2.denominator)**2 == y2.denominator
            P = (X, Fr(isqrt(y2.numerator), isqrt(y2.denominator))); assert C.on(P)
            hP = hhat(C, P, B, U)
            if (m, n) not in tw:
                sq_curves += 1; inf_order += hP[0] > 0
                continue
            r = isqrt(C.a6); assert r * r == C.a6
            Q = (Fr(0), Fr(r)); assert C.on(Q)
            hQ = hhat(C, Q, B, U); hS = hhat(C, C.add(P, Q), B, U)
            assert hP[0] > 0 and hQ[0] > 0
            pq = ((hS[0] - hP[1] - hQ[1]) / 2, (hS[1] - hP[0] - hQ[0]) / 2)        # <P, Q>
            pq2 = max(pq[0]**2, pq[1]**2)
            det_lo = hP[0] * hQ[0] - pq2
            assert det_lo > 0, (m, n, s, float(det_lo))
            least = det_lo if least is None else min(least, det_lo)
    print('twice a square: %d ratios, %d curves: the points of t = +-1 and t = 0 are independent on every one '
          '(least lower bound for the regulator %.3f), so every rank is at least 2' % (len(tw), 2*len(tw), float(least)))
    assert inf_order == sq_curves
    print('a square: %d curves: the point of t = +-1 has infinite order on %d' % (sq_curves, inf_order))
    print('sigma_twice OK')


if __name__ == '__main__':
    main()
