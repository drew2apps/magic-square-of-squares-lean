#!/usr/bin/env python3
"""A numerical test of the height relations behind dm4.py, over number fields (absolute canonical heights from Sage).

For a rational t0 the points Q1 = (t0, sqrt(F0 GL)), Q2 = (-t0, sqrt(F0 GU)), Q3 = (1/t0, sqrt(F8 GL)/t0^2) of C are
defined over a quadratic or biquadratic field.  The script prints
    N(a1) - H = hhat(2 Q1 - R0)/4 - h(t0),   hhat(Q1 + Q2 - R0) - 2 h(t0),   hhat(Q1 - Q2) - 2 h(t0),
and the same for Q3.  According to dm4.py all of these are bounded independently of t0.  Over these fields the two
points Q1 + Q2 - R0 and Q1 - Q2 are conjugate, so their heights are equal: the test shows that the differences
stay bounded as h(t0) grows, not the size of <a1, a2> for rational points.  It is a plausibility check, not a proof.
usage: dm4_validate.py m,n t0 [t0 ...]     (needs passagemath-schemes; a few seconds per value)"""
import sys, math
from sage.all__sagemath_schemes import *
from sage.rings.number_field.number_field import NumberField
m, n = (int(v) for v in sys.argv[1].split(','))
S = m*m + n*n
R = PolynomialRing(QQ, 'T'); T = R.gen()
def qpoly(x): return (n*n*x**2 + m*m) * (2*S*x**2 - 8*m*n*x + 2*S)
t1 = QQ(m)/n; q0 = QQ(2*m*(m*m - n*n))/n
qs = qpoly(t1 + T); e0, d, c, b, a = [qs[i] for i in range(5)]
a1 = d/q0; a2 = c - d**2/(4*q0**2); a3 = 2*q0*b; a4 = -4*q0**2*a; a6 = a2*a4
E0 = EllipticCurve(QQ, [a1, a2, a3, a4, a6])
def toE0(tt, y):
    tau = tt - t1
    X = (2*q0*(y + q0) + d*tau) / tau**2
    Y = (4*q0**2*(y + q0) + 2*q0*(d*tau + c*tau**2) - d**2*tau**2/(2*q0)) / tau**3
    return X, Y
def sqfree(x):
    x = QQ(x); v = x.numerator() * x.denominator(); return ZZ(v).squarefree_part()
for s in sys.argv[2:]:
    t0 = QQ(s)
    A = qpoly(t0); B = qpoly(-t0); C3 = (m*m*t0**2 + n*n) * (2*S*t0**2 - 8*m*n*t0 + 2*S)
    for name, tt2, rad, scale in (('pm', -t0, B, 1), ('inv', 1/t0, C3, 1/t0**2)):
        dA, dB = sqfree(A), sqfree(rad)
        x = polygen(QQ)
        if dA == 1 or dB == 1 or dA == dB: print(t0, name, 'degenerate field'); continue
        Krel = NumberField([x*x - dA, x*x - dB], names=('r1', 'r2')); K = Krel.absolute_field('z'); fr, to = K.structure(); r1, r2 = to(Krel.gens()[0]), to(Krel.gens()[1])
        yA = r1 * (A / dA).sqrt(); yB = r2 * (rad / dB).sqrt() * scale
        EK = E0.change_ring(K)
        P1 = EK(toE0(K(t0), yA)); P2 = EK(toE0(K(tt2), yB)); iP2 = EK(toE0(K(tt2), -yB)); R0 = EK([-a2, a1*a2 - a3])
        H = math.log(max(abs(t0.numerator()), abs(t0.denominator())))
        hN = float((2*P1 - R0).height()) / 4; hU = float((P1 - iP2).height()); hV = float((P1 - P2).height())
        assert P1 - iP2 == P1 + P2 - R0
        print('t0=%s %s: h(t0)=%.3f  N(a1)-H=%.3f  hhat(U)-2H=%.3f  hhat(V)-2H=%.3f  <a1,a2>=%.4f' % (t0, name, H, hN - H, hU - 2*H, hV - 2*H, (hU - hV)/4), flush=True)
