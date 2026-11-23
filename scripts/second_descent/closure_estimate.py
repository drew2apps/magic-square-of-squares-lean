#!/usr/bin/env python3
"""The estimates behind the outline of a closure of the last pair, with every constant computed (EXPERIMENT; SECOND_DESCENT.md, 6).

Statements used, as read from arXiv:1602.04097 (Checcoli, Veneziano, Viada) by an automated reader, not checked by a person:
 (T)  Theorem 1.2, N = 2:  h^(P) <= C1 h2(C) deg C + C2(E) (deg C)^2 + C3(E) for rational P on a curve C of genus >= 2 in E^2, E(Q) of rank 1;
 (P)  Proposition 3.2:  -c2(E) <= h2(P) - h^(P) <= c1(E) for P in E (algebraic points: the general constants);
 (Z)  Zhang:  h2(C) <= 2 deg C mu2(C), mu2(C) = inf of the theta for which the points of C of height h2 <= theta are Zariski dense.
h^ is three times Silverman's canonical height, that is 1.5 times PARI's ellheight; h2 is the Weil height with the l2 norm at infinity.

The curve.  V: the (t, y1, y1', y2, y2') with y1^2 = f1 = d Q1(-i t), y1'^2 = f2 = d Q1(i t), y2^2 = f3 = d Q2(-i t), y2'^2 = f4 = d Q2(i t).
The product f1 f2 f3 f4 is a constant times a square and no other product is (checked below, with the constant), so V has two
components W, each irreducible of degree 8 over the line of t, a cover with group (Z/2)^3.  The four f have 8 distinct roots in
all, each a root of exactly two of them (checked below), so the cover has 8 branch points with inertia of order 2 and
2g - 2 = 8 (-2) + 8 * 4: W has genus 9.  (An earlier version said 16 branch points and genus 25; that was wrong.)
On W, A = P1 + P1' and B = tau(R - R') are sums of two maps of degree 4 to E1, through the four quotients y1, y1', y2, y2'.
The pull-backs of the invariant differential are dt/y1, dt/y1', dt/y2, dt/y2' up to nonzero constants: four eigenvectors of the
group for four different characters, hence linearly independent, so no u A + v B with (u, v) != (0, 0) is constant, and the image
C is not in a translate of a proper algebraic subgroup of E1 x E1 (LAST_PAIR.md, Lemma 2).
deg(A) <= 2 (4 + 4) = 16 by the parallelogram law for degrees, the same for B, and deg C <= 3 (16 + 16) = 96.  The numbers are also
printed for deg C <= 192, twice that, to show how little depends on it.

Logarithms.  Every log of an integer that enters an upper bound is replaced by (number of bits) * log 2, which is larger, so the
printed bounds are upper bounds whatever the rounding of the few floating-point operations left (the margin is a factor 10^7).

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
def ulog(x):
    """An upper bound for log|x|, x a nonzero integer: bits * log 2."""
    return ZZ(abs(ZZ(x))).nbits()*0.6931471806
import itertools
sq = [''.join(str(i + 1) for i in c) for k in (1, 2, 3, 4) for c in itertools.combinations(range(4), k) if is_const_times_square(prod(f[i] for i in c))]
print('products of f1..f4 that are a constant times a square in K[t]:', sq, '(expected: only 1234)')
gall = prod(f); const = gall.factor().unit(); print('f1 f2 f3 f4 = c g^2 with c a square in K: %s' % K(const).is_square())
rad = prod(f).radical(); mult = sorted(sum(1 for fi in f if fi % g == 0) for g, _ in rad.factor() for _ in range(g.degree()))
print('distinct roots of f1 f2 f3 f4: %d; number of the f that each root divides: %s; genus of a component: 1 + (8 * (-2) + %d * 4)/2 = %d' % (rad.degree(), sorted(set(mult)), rad.degree(), 1 + (8*(-2) + rad.degree()*4)//2))
hq = lambda q: max(ulog(QQ(q).numerator()), ulog(QQ(q).denominator()))
def one(Q, r):
    E, A0 = weierstrass(d*Q, r); Em = E.minimal_model(); A, B = -27*Em.c4(), -54*Em.c6(); Es = EllipticCurve([A, B]); iso = E.isomorphism_to(Es)
    u, rr = iso.u, iso.r                                    # x_short = (x - r)/u^2 for Sage's isomorphism (u, r, s, t)
    alpha, beta = 1/u**2, -rr/u**2
    HX = hq(alpha*A0) + hq(r) + hq(beta) + 2*0.6931471806
    j = Es.j_invariant(); D = -16*(4*A**3 + 27*B**2); hinf = hq(j)                                    # hinf(j) <= hW(j)
    c1g = (ulog(D) + hinf)/4 + hq(j)/8 + (ulog(A) + ulog(B))/2 + 3.724
    c2g = (ulog(D) + hinf)/4 + (ulog(A) + ulog(B))/2 + 4.015
    c1q = (ulog(D) + hinf)/4 + hq(j)/8 + ulog(abs(A) + abs(B) + 3)/2 + 2.919                       # the first of the two expressions of Table 1; the minimum is not needed
    c2q = (ulog(D) + hinf)/4 + ulog(abs(A) + abs(B) + 3)/2 + 3.21
    HX, c1g, c2g, c1q, c2q = (float(v) for v in (HX, c1g, c2g, c1q, c2q))
    Kpt = float(1.5*HX + 0.5*ulog(1 + abs(A) + abs(B)) + 0.5*log(3) + c2g)
    return dict(E=E, Es=Es, iso=iso, A=A, B=B, HX=HX, c1g=c1g, c2g=c2g, c1q=c1q, c2q=c2q, Kpt=Kpt, A0=A0)
c = one(Q1, r1); c2_ = one(Q2, r2)
print('E1: HX = %.2f, c1 (general) = %.2f, c2 (general) = %.2f, c1 (over Q) = %.2f, c2 (over Q) = %.2f; bound for h^ of a point above a root of unity: %.1f' % (c['HX'], c['c1g'], c['c2g'], c['c1q'], c['c2q'], c['Kpt']))
print('E2: HX = %.2f, c1 (general) = %.2f, c2 (general) = %.2f; bound for h^ of a point above a root of unity: %.1f' % (c2_['HX'], c2_['c1g'], c2_['c2g'], c2_['Kpt']))
mu = 4*c['Kpt'] + 4*c2_['Kpt'] + c['c1g'] + c2_['c1g']
print('mu2(C) <= 4 * %.1f + 4 * %.1f + %.2f + %.2f = %.1f  (the points of C above the roots of unity)' % (c['Kpt'], c2_['Kpt'], c['c1g'], c2_['c1g'], mu))
C1 = 72.251
bounds = {}
for deg in (96, 192):
    h2C = 2*deg*mu
    for name, c1, c2 in (('constants over Q', c['c1q'], c['c2q']), ('general constants', c['c1g'], c['c2g'])):
        C2 = C1*(6.019 + 4*c1); C3 = 4*c2; bounds[(deg, name)] = C1*h2C*deg + C2*deg*deg + C3
        print('(T) with deg C <= %d, h2(C) <= %.4e and the %s: C2 = %.1f, C3 = %.1f; every rational point of C has h^ <= %.4e' % (deg, h2C, name, C2, C3, bounds[(deg, name)]))
# real values: the point G of the degenerate solution (t = 1, a root of unity) against (P) and against the bound for mu
conj = [a for a in K.automorphisms() if a(I) == -I][0]; E1 = c['E']; EK = E1.change_ring(K)
vd = d*Q1.change_ring(K)(-I); wd = 1/(-I - r1); Pd = EK(c['A0']*wd, c['A0']*vd.sqrt()*wd*wd); Gs = Pd + EK(conj(Pd[0]), conj(Pd[1])); G1 = E1(QQ(Gs[0]), QQ(Gs[1]))
Em = E1.minimal_model(); em = pari.ellinit([pari(str(a)) for a in Em.a_invariants()]); Gm = E1.isomorphism_to(Em)(G1)
hG = float(pari.ellheight(em, [pari(str(Gm[0])), pari(str(Gm[1]))])); G = c['iso'](G1)
X, Y = QQ(G[0]), QQ(G[1]); den = lcm(X.denominator(), Y.denominator()); co = [abs(X*den), abs(Y*den), den]
h2G = log(float(sum(v*v for v in co))**0.5) - log(float(gcd(co)))
print('real point: G has h^ = 1.5 * %.4f = %.3f and h2 = %.3f; h2 - h^ = %.3f, to be between -c2 = %.2f and c1 = %.2f: %s; h^(G) <= 4 * %.1f: %s' % (
    hG, 1.5*hG, h2G, h2G - 1.5*hG, -c['c2q'], c['c1q'], -c['c2q'] <= h2G - 1.5*hG <= c['c1q'], c['Kpt'], 1.5*hG <= 4*c['Kpt']))
low = 1.5*(hG - 0.001)*(amax + 1)**2
print('sieve: a non-degenerate solution has a or b above %d, so its point has h^ >= 1.5 * %.3f * %d^2 = %.4e' % (amax, hG - 0.001, amax + 1, low))
for deg in (96, 192): print('ratio of the sieve bound to the largest bound of (T) with deg C <= %d: %.3e' % (deg, low/max(v for k, v in bounds.items() if k[0] == deg)))
for deg in (96, 192): print('so every solution has a, b <= %d (deg C <= %d); the sieve covers 0 <= a, b <= %d: %s' % (int((max(v for k, v in bounds.items() if k[0] == deg)/(1.5*(hG - 0.001)))**0.5), deg, amax, int((max(v for k, v in bounds.items() if k[0] == deg)/(1.5*(hG - 0.001)))**0.5) <= amax))
low1 = 1.5*(hG - 0.001)*596310**2
if amax > 596309: print('the same with the first sieve (a, b <= 596309): %.3e and %.3e' % (low1/max(v for k, v in bounds.items() if k[0] == 96), low1/max(v for k, v in bounds.items() if k[0] == 192)))
