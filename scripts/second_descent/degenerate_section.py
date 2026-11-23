#!/usr/bin/env python3
"""The degenerate point has infinite order for every ratio (EXPERIMENT; LAST_PAIR.md, section 9, Theorem C).

Over Q(tau), tau = m/n (n = 1 after scaling, which changes d Q1 by a square): d = 2(m^2 - n^2), and
    E1: the Jacobian of Y^2 = d q+(s)(m - n s)(n + m s), in the model of weierstrass(), with coefficients in Q[tau].
At s = -i:  Q1(-i) = 2(m^2 - n^2)(m + i n)^2, so d Q1(-i) is a square and P = (s, Y) = (-i, 2(m^2 - n^2)(m + i n)) is a point
over Q(i)(tau); G = P + conj(P) has coordinates in Q(tau), computed here by the chord through P and conj(P) (slope -2(m^2 - n^2)).
A rational torsion point has order 1, ..., 10 or 12 (Mazur).  For a rational tau0 at which the curve is not singular,
specialisation is a homomorphism and a multiple of G specialises to O exactly when its x has a pole.  So if G(tau0) had
order k, one of these would vanish at tau0:   y(jG) for k = 2j;   x((j+1)G) - x(jG) for k = 2j + 1.
The program prints the rational roots of the discriminant and of the numerators of these functions for k = 2, ..., 10, 12.
If they are only 0, 1, -1 (m n (m^2 - n^2) = 0, no ratio), G has infinite order for every ratio, and E1 has rank at least 1.
usage: degenerate_section.py"""
from second_descent import *
from sage.rings.fraction_field import FractionField
Rt = PolynomialRing(QQ, 'tau'); tau = Rt.gen(); F = FractionField(Rt)
Rs = PolynomialRing(F, 's'); s = Rs.gen()
m, n = F(tau), F(1); D = m*m - n*n; S = m*m + n*n; d = 2*D
Q = d*(S*s**2 - 2*D*s + S)*(m - n*s)*(n + m*s); r = m/n
A, B, C, E0 = [Q.derivative(k)(r)/factorial(k) for k in range(1, 5)]
assert all(F(c).denominator() == 1 for c in (A, B, C, E0))
E = EllipticCurve(F, [0, B, 0, A*C, A*A*E0])
# P = (A w, A Y w^2) with w = 1/(-i - tau), Y = 2 D (tau + i):  x(P) = a + b i, y(P) = c + e i
a = -A*tau/(tau**2 + 1); b = A/(tau**2 + 1); c = 2*A*D*tau/(tau**2 + 1); e = -2*A*D/(tau**2 + 1)
assert c**2 - e**2 == a**3 - 3*a*b**2 + B*(a**2 - b**2) + A*C*a + A*A*E0 and 2*c*e == 3*a**2*b - b**3 + 2*B*a*b + A*C*b      # P is on E1
lam = e/b; xG = lam**2 - B - 2*a; yG = -(lam*(xG - a) + c); G = E(xG, yG)
print('x(G) =', F(xG).factor()); print('y(G) =', F(yG).factor())
print('discriminant of the model:', Rt(E.discriminant()).factor())
mult = {1: G}
for k in range(2, 7): mult[k] = mult[k - 1] + G
num = lambda f: Rt(F(f).numerator())
tests = {2: num(mult[1][1]), 3: num(mult[2][0] - mult[1][0]), 4: num(mult[2][1]), 5: num(mult[3][0] - mult[2][0]), 6: num(mult[3][1]),
         7: num(mult[4][0] - mult[3][0]), 8: num(mult[4][1]), 9: num(mult[5][0] - mult[4][0]), 10: num(mult[5][1]), 12: num(mult[6][1])}
bad = set(x for x, _ in Rt(E.discriminant()).roots(QQ))
print('rational roots of the discriminant:', sorted(bad))
for k, f in sorted(tests.items()):
    assert f != 0
    rts = sorted(x for x, _ in f.roots(QQ)); bad |= set(rts)
    print('order %d: a polynomial of degree %d in tau, rational roots %s' % (k, f.degree(), rts))
print('rational tau that are not excluded:', sorted(bad), '; G has infinite order for every other rational tau:', sorted(bad) == [-1, 0, 1])
