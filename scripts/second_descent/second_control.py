#!/usr/bin/env python3
"""Control for the equation step of the second descent: the polynomial.
For a random x1 in K and y1 = sqrt(f(x1)), the product (x(P+T) - g(x1)) (x(P-T) - g(x1)), computed with Sage's group law over
K(y1), must equal pol(x1) / (D(x1)^2 (x1 - x0)^4)   [T not of order 2].
The other control is in the output of second3.py: for every ratio with m^2 + n^2 a square that is closed, the equations return
t = 1 and t = -1, the degenerate solutions, which are not put in by hand (field degenerate_found_by_equations)."""
import sys, json
from math import gcd
from fractions import Fraction as Fr
from sage.all__sagemath_schemes import *

from second_descent import K, I, sq, weierstrass

def polynomial(m, n, trials=3):
    Q = (m - n*sq)*(n + m*sq)*(1 - sq**2); r = QQ(m)/n
    E, A = weierstrass(Q, r); EK = E.change_ring(K); R = PolynomialRing(K, 'x'); x = R.gen(); a2 = EK.a2(); f = x**3 + a2*x**2 + EK.a4()*x + EK.a6()
    n_checked = 0
    for T in EK.torsion_points():
        if T.is_zero() or T[1] == 0: continue
        x0, y0 = T[0], T[1]; Dn = 2*r*x + A; hD = (-A*x + (a2 + x + x0)*Dn)*(x - x0)**2
        pol = ((f + y0**2)*Dn - hD)**2 - 4*y0**2*f*Dn**2; ratios = []
        for k in range(trials):
            x1 = K(3 + 2*k) + I*(k + 1)/7
            Lf = K.extension(PolynomialRing(K, 'Y')([-f(x1), 0, 1]), 'y1'); y1 = Lf.gen()
            EL = EK.change_ring(Lf); P = EL([x1, y1]); TL = EL([x0, y0]); g = -A*x1/(2*r*x1 + A)
            lhs = ((P + TL)[0] - g)*((P - TL)[0] - g)
            ratios.append(lhs*Dn(x1)**2*(x1 - x0)**4/pol(x1))
        assert all(c == 1 for c in ratios), ratios
        n_checked += 1
    print('polynomial identity', (m, n), ': %d torsion points checked, ratio 1 each time' % n_checked)

if __name__ == '__main__':
    for p in [(4, 3), (91, 60), (15, 8), (572, 555)]: polynomial(*p)
