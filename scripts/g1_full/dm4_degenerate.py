#!/usr/bin/env python3
"""The relations of dm4.py on the only rational configurations available: the degenerate solutions.

When m^2 + n^2 is a square, t0 = 1 makes the five forms squares (a degenerate square), and C: y^2 = F0 GL has
rational points Q1 above t = 1 and Q2 above t = -1.  For these the script prints N(a1), N(a2) (to compare with
H = h(1) = 0 and the interval [-eL, eU] of dm4.py), the heights of U = Q1 + Q2 - R0 and V = Q1 - Q2 (to compare
with 2H), the product <a1, a2> = (hhat(U) - hhat(V))/4 (to compare with kappa), and checks the consequence
x(U) x(V) = (4 m n S)^2 of the identity (I+-) for one of the points of order 2.
usage: dm4_degenerate.py m,n [m,n ...]     (needs passagemath-schemes)"""
import sys
from sage.all__sagemath_schemes import *
def run(m, n):
    S = m*m + n*n; R = PolynomialRing(QQ, 'T'); T = R.gen()
    assert ZZ(S).is_square()
    def qpoly(x): return (n*n*x**2 + m*m) * (2*S*x**2 - 8*m*n*x + 2*S)
    t1 = QQ(m)/n; q0 = QQ(2*m*(m*m - n*n))/n
    qs = qpoly(t1 + T); e0, d, c, b, a = [qs[i] for i in range(5)]
    a1 = d/q0; a2 = c - d**2/(4*q0**2); a3 = 2*q0*b; a4 = -4*q0**2*a; a6 = a2*a4
    E0 = EllipticCurve(QQ, [a1, a2, a3, a4, a6]); R0 = E0([-a2, a1*a2 - a3])
    def toE0(tt, y):
        tau = tt - t1
        return E0([(2*q0*(y + q0) + d*tau) / tau**2, (4*q0**2*(y + q0) + 2*q0*(d*tau + c*tau**2) - d**2*tau**2/(2*q0)) / tau**3])
    h = lambda P: float(P.height())
    Q1 = toE0(QQ(1), qpoly(QQ(1)).sqrt()); Q2 = toE0(QQ(-1), qpoly(QQ(-1)).sqrt())
    U = Q1 + Q2 - R0; V = Q1 - Q2
    e = [P[0] for P in E0.torsion_points() if P.order() == 2]
    ok = any((U[0] - ee)*(V[0] - ee) == (4*m*n*S)**2 for ee in e); assert ok
    print('(%d,%d): N(a1) = %.3f, N(a2) = %.3f, hhat(U) = %.3f, hhat(V) = %.3f, <a1,a2> = %.3f, x(U) x(V) = (4mnS)^2: %s' % (m, n, h(2*Q1 - R0)/4, h(2*Q2 - R0)/4, h(U), h(V), (h(U) - h(V))/4, ok))
for s in sys.argv[1:]:
    m, n = (int(v) for v in s.split(',')); run(m, n)
