#!/usr/bin/env python3
"""The estimates behind the outline of a closure of the last pair, with every constant computed (EXPERIMENT; SECOND_DESCENT.md, 6).

Statements used, as read from arXiv:1602.04097 (Checcoli, Veneziano, Viada) by an automated reader, not checked by a person:
 (T)  Theorem 1.2, N = 2:  h^(P) <= C1 h2(C) deg C + C2(E) (deg C)^2 + C3(E) for rational P on a curve C of genus >= 2 in E^2, E(Q) of rank 1;
 (P)  Proposition 3.2:  -c2(E) <= h2(P) - h^(P) <= c1(E) for P in E (algebraic points: the general constants);
 (Z)  Zhang:  h2(C) <= 2 deg C mu2(C), mu2(C) = inf of the theta for which the points of C of height h2 <= theta are Zariski dense.
h^ is three times Silverman's canonical height, that is 1.5 times PARI's ellheight; h2 is the Weil height with the l2 norm at infinity.

The curve.  V: the (t, y1, y1', y2, y2') with y1^2 = d Q1(-i t), y1'^2 = d Q1(i t), y2^2 = d Q2(-i t), y2'^2 = d Q2(i t).  The product of the
four right sides is a constant times a square (checked below), so V has two components, each of degree 8 over the line of t;
no other product of them is a square, so each component is irreducible.  On a component, A = P1 + P1' and B = tau(R - R') are sums of
two maps of degree 4 to E1 (through the four quotients y1, y1', y2, y2', which belong to four different characters of the group of
the cover, so they are independent: the image C is not in a translate of a proper algebraic subgroup, and has genus >= 2).
deg(A) <= 2 (4 + 4) = 16, the same for B, and deg C <= 3 (16 + 16) = 96.

The height.  Above t = a root of unity, s = -i t is a root of unity, x(P1) = A0/(s - r) on the model of weierstrass(), and
X = alpha x + beta on the short model.  h(X) <= h(alpha A0) + h(r) + h(beta) + 2 log 2 =: HX;  hW(X:Y:1) <= 1.5 h(X) + 0.5 log(1+|A|+|B|);
h2 <= hW + 0.5 log 3;  h^(P1) <= h2 + c2gen;  h^(A) <= 2 h^(P1) + 2 h^(P1') (parallelogram law);  h2((A, B)) <= h^(A) + h^(B) + 2 c1gen.
usage: closure_estimate.py m,n d a_max"""
import sys
from math import log
from second_descent import *
import cypari2
pari = cypari2.Pari(); pari.allocatemem(2*10**9)
m, n = (int(v) for v in sys.argv[1].split(',')); d = int(sys.argv[2]); amax = int(sys.argv[3])
fam = family(m, n); (Q1, r1), (Q2, r2) = fam['+13'], fam['+24']
R = PolynomialRing(K, 't'); t = R.gen()
f = [(d*Q1.change_ring(K))(-I*t), (d*Q1.change_ring(K))(I*t), (d*Q2.change_ring(K))(-I*t), (d*Q2.change_ring(K))(I*t)]
def is_const_times_square(g): return all(e % 2 == 0 for _, e in g.factor())
import itertools
sq = [''.join(str(i + 1) for i in c) for k in (1, 2, 3, 4) for c in itertools.combinations(range(4), k) if is_const_times_square(prod(f[i] for i in c))]
print('products of f1..f4 that are a constant times a square in K[t]:', sq, '(expected: only 1234)')
hq = lambda q: log(max(abs(QQ(q).numerator()), abs(QQ(q).denominator())))
def one(Q, r):
    E, A0 = weierstrass(d*Q, r); Em = E.minimal_model(); A, B = -27*Em.c4(), -54*Em.c6(); Es = EllipticCurve([A, B]); iso = E.isomorphism_to(Es)
    u, rr = iso.u, iso.r                                    # x_short = (x - r)/u^2 for Sage's isomorphism (u, r, s, t)
    alpha, beta = 1/u**2, -rr/u**2
    HX = hq(alpha*A0) + hq(r) + hq(beta) + 2*log(2)
    j = Es.j_invariant(); D = -16*(4*A**3 + 27*B**2); hinf = max(log(abs(float(j))), 0.0)
    c1g = (log(abs(D)) + hinf)/4 + hq(j)/8 + (log(abs(A)) + log(abs(B)))/2 + 3.724
    c2g = (log(abs(D)) + hinf)/4 + (log(abs(A)) + log(abs(B)))/2 + 4.015
    hWW = log(max(1.0, float(abs(A))**0.5, float(abs(B))**(1/3)))
    c1q = min((log(abs(D)) + hinf)/4 + hq(j)/8 + log(abs(A) + abs(B) + 3)/2 + 2.919, 3*hWW + 4.709)
    c2q = min((log(abs(D)) + hinf)/4 + log(abs(A) + abs(B) + 3)/2 + 3.21, 3*hWW/2 + 2.427)
    HX, c1g, c2g, c1q, c2q = (float(v) for v in (HX, c1g, c2g, c1q, c2q))
    Kpt = float(1.5*HX + 0.5*log(float(1 + abs(A) + abs(B))) + 0.5*log(3) + c2g)
    return dict(E=E, Es=Es, iso=iso, A=A, B=B, HX=HX, c1g=c1g, c2g=c2g, c1q=c1q, c2q=c2q, Kpt=Kpt, A0=A0)
c = one(Q1, r1); c2_ = one(Q2, r2)
print('E1: HX = %.2f, c1 (general) = %.2f, c2 (general) = %.2f, c1 (over Q) = %.2f, c2 (over Q) = %.2f; bound for h^ of a point above a root of unity: %.1f' % (c['HX'], c['c1g'], c['c2g'], c['c1q'], c['c2q'], c['Kpt']))
print('E2: HX = %.2f, c1 (general) = %.2f, c2 (general) = %.2f; bound for h^ of a point above a root of unity: %.1f' % (c2_['HX'], c2_['c1g'], c2_['c2g'], c2_['Kpt']))
mu = 4*c['Kpt'] + 4*c2_['Kpt'] + 2*c['c1g']; deg = 96; h2C = 2*deg*mu
print('mu2(C) <= 4 * %.1f + 4 * %.1f + 2 * %.2f = %.1f;  deg C <= %d;  h2(C) <= 2 deg mu = %.4e' % (c['Kpt'], c2_['Kpt'], c['c1g'], mu, deg, h2C))
C1 = 72.251
for name, c1, c2 in (('constants over Q', c['c1q'], c['c2q']), ('general constants', c['c1g'], c['c2g'])):
    C2 = C1*(6.019 + 4*c1); C3 = 4*c2; bound = C1*h2C*deg + C2*deg*deg + C3
    print('(T) with the %s: C2 = %.1f, C3 = %.1f; every rational point of C has h^ <= %.4e' % (name, C2, C3, bound))
# real values: the point G of the degenerate solution (t = 1, a root of unity) against (P) and against the bound for mu
Em = c['E'].minimal_model(); em = pari.ellinit([pari(str(a)) for a in Em.a_invariants()]); g = pari.ellrank(em, 2)[3][0]
hG = float(pari.ellheight(em, g)); G = c['iso'](Em.isomorphism_to(c['E'])(Em(QQ(str(g[0])), QQ(str(g[1])))))
X, Y = QQ(G[0]), QQ(G[1]); den = lcm(X.denominator(), Y.denominator()); co = [abs(X*den), abs(Y*den), den]
h2G = log(float(sum(v*v for v in co))**0.5) - log(float(gcd(co)))
print('real point: G has h^ = 1.5 * %.4f = %.3f and h2 = %.3f; h2 - h^ = %.3f, to be between -c2 = %.2f and c1 = %.2f: %s; h^(G) <= 4 * %.1f: %s' % (
    hG, 1.5*hG, h2G, h2G - 1.5*hG, -c['c2q'], c['c1q'], -c['c2q'] <= h2G - 1.5*hG <= c['c1q'], c['Kpt'], 1.5*hG <= 4*c['Kpt']))
low = 1.5*hG*(amax + 1)**2
print('sieve: a non-degenerate solution has a or b above %d, so its point has h^ >= %.4e' % (amax, low))
print('ratio of the sieve bound to the largest bound of (T): %.3e' % (low/(C1*h2C*deg + C1*(6.019 + 4*c['c1g'])*deg*deg + 4*c['c2g'])))
