#!/usr/bin/env python3
"""The second descent's closures of the ratios open in the census, rechecked without Sage or PARI (sympy only).

    python3 second_exact.py          the five closures by a rank bound 0 from descent (second_123.jsonl)

For each closing curve (SECOND_DESCENT.md, sections 2 and 3) this script recomputes, from m, n alone:
 1. the quartic Q(s) (a product of the linear forms L1..L6 and q+ of section 2), the class d at the two degenerate
    solutions (checked: d*Q(s) is a square in Q(i) at both, and the record's d is this d), a rational root r of Q;
 2. the Weierstrass model E of Y^2 = d Q(s) by X = A/(s - r) (the map of quartic_root_weierstrass in G1Sigma.lean),
    and that the record's curve of rank 0 is isomorphic over Q to E or to its twist by -1, as the record says;
 3. the torsion of E(Q(i)): every point found by halving from the points of order 2 over Q(i), closed under
    addition, and the bound: the gcd of #E(F_p) over primes p = 1 (mod 4) of good reduction, 13 <= p < 400 (Q(i)
    embeds in Q_p, and torsion of order prime to p injects into E(F_p)); the two agree;
 4. the equations of section 3: with g(x) = -A x/(2 r x + A) (the image of s -> -s), x(P + T) = g(x(P)) when the
    twist by -1 has rank 0, x(P - T) = g(x(P)), that is x(-P + T), when E has rank 0, for every torsion point T;
    each is a polynomial equation in x(P) over Q(i) (squared to remove y), solved by factoring over Q(i);
 5. every root gives s0 = r + A/x and t0 = i s0; the rational t0 are tested in the five forms.
Inputs not recomputed here: that the rank is 0 (PARI's descent in the record; Magma, magma_second_ranks_raw.txt;
eclib for three), and the first descent (that a solution has, form by form, the class of a degenerate solution)."""
import os, sys, json
from fractions import Fraction as Fr
import sympy as sp

here = os.path.dirname(os.path.abspath(__file__))
I = sp.I
s, x = sp.symbols('s x')


def forms(m, n):
    S = m*m + n*n; D = m*m - n*n
    return {'1': m - n*s, '2': m + n*s, '3': n + m*s, '4': n - m*s, '5': 1 - s, '6': 1 + s, '+': S*s**2 - 2*D*s + S}


def quartic(m, n, name):
    F = forms(m, n)
    Q = sp.Integer(1)
    for ch in name: Q *= F[ch]
    return sp.expand(Q)


def is_gauss_square(z):
    """z in Q(i): is it a square in Q(i)?"""
    z = sp.nsimplify(sp.expand(z))
    if z == 0: return True
    a, b = sp.re(z), sp.im(z)
    N = a*a + b*b
    r = sp.sqrt(N)
    if not r.is_Rational: return False
    for c in ((a + r)/2, (a - r)/2):
        if c >= 0:
            u = sp.sqrt(c)
            if u.is_Rational:
                if u == 0: return (b == 0 and a <= 0 and sp.sqrt(-a).is_Rational)
                v = b/(2*u)
                if sp.expand((u + I*v)**2 - z) == 0: return True
    return False


def degenerate_s(m, n):
    S = m*m + n*n
    if sp.sqrt(S).is_Rational: return [-I*1, -I*(-1)]          # t = 1, -1 (s = -i t)
    return ['inf0', 'inf1']                                      # t = 0 and t = infinity


def squarefree(q):
    q = Fr(q); num = q.numerator * q.denominator
    sgn = -1 if num < 0 else 1; num = abs(num)
    f = sp.factorint(num)
    return sgn * int(sp.prod([p for p, e in f.items() if e % 2]))


def weierstrass(Qd, r):
    """Y^2 = Qd(s) with Qd(r) = 0: the model y^2 = X^3 + c2 X^2 + c1 c3 X + c1^2 c4, X = c1/(s - r)."""
    h = sp.symbols('h')
    P = sp.Poly(sp.expand(Qd.subs(s, r + h)), h)
    c = [P.coeff_monomial(h**k) for k in range(5)]
    assert c[0] == 0 and c[1] != 0
    return c[1], (c[2], c[1]*c[3], c[1]**2*c[4])


def c4c6(a2, a4, a6):
    b2, b4, b6 = 4*a2, 2*a4, 4*a6
    c4 = b2*b2 - 24*b4; c6 = -b2**3 + 36*b2*b4 - 216*b6
    return sp.Rational(c4), sp.Rational(c6)


def isomorphic(e1, e2):
    """Two models y^2 = x^3 + a x^2 + b x + c (general a-invariants for e2 allowed) isomorphic over Q?"""
    def inv(ai):
        a1, a2, a3, a4, a6 = [sp.Rational(v) for v in ai]
        b2 = a1*a1 + 4*a2; b4 = 2*a4 + a1*a3; b6 = a3*a3 + 4*a6
        return b2*b2 - 24*b4, -b2**3 + 36*b2*b4 - 216*b6
    c4, c6 = inv(e1); d4, d6 = inv(e2)
    if c4 == 0 or c6 == 0:
        return c4 == d4 == 0 or c6 == d6 == 0          # not needed here
    u2 = (d6/c6) / (d4/c4)
    return sp.sqrt(u2).is_Rational and d4 == u2**2*c4 and d6 == u2**3*c6


class EK:
    """y^2 = x^3 + a2 x^2 + a4 x + a6 over Q(i), points as (x, y) of sympy numbers, None = O."""
    def __init__(self, a2, a4, a6): self.a2, self.a4, self.a6 = a2, a4, a6
    def f(self, X): return X**3 + self.a2*X**2 + self.a4*X + self.a6
    def add(self, P, Q):
        if P is None: return Q
        if Q is None: return P
        if sp.expand(P[0] - Q[0]) == 0:
            if sp.expand(P[1] + Q[1]) == 0: return None
            l = (3*P[0]**2 + 2*self.a2*P[0] + self.a4) / (2*P[1])
        else:
            l = (Q[1] - P[1]) / (Q[0] - P[0])
        l = sp.nsimplify(sp.expand(sp.radsimp(l)))
        X = sp.expand(l*l - self.a2 - P[0] - Q[0]); Y = sp.expand(l*(P[0] - X) - P[1])
        return (sp.nsimplify(X), sp.nsimplify(Y))
    def neg(self, P): return None if P is None else (P[0], sp.expand(-P[1]))
    def roots_K(self, poly):
        """Roots in Q(i) of a polynomial in x with coefficients in Q(i)."""
        out = []
        for fac, _ in sp.factor_list(sp.expand(poly), x, extension=I)[1]:
            p = sp.Poly(fac, x)
            if p.degree() == 1:
                a, b = p.all_coeffs(); out.append(sp.nsimplify(sp.expand(sp.radsimp(-b/a))))
        return out
    def points_with_x(self, X):
        v = sp.expand(self.f(X))
        if v == 0: return [(X, sp.Integer(0))]
        if not is_gauss_square(v): return []
        y = sp.symbols('y'); ys = self.roots_y(v)
        return [(X, yy) for yy in ys]
    def roots_y(self, v):
        y = sp.symbols('y')
        out = []
        for fac, _ in sp.factor_list(y**2 - v, y, extension=I)[1]:
            p = sp.Poly(fac, y)
            if p.degree() == 1:
                a, b = p.all_coeffs(); out.append(sp.nsimplify(sp.expand(-b/a)))
        return out
    def halves(self, P):
        """All Q over Q(i) with 2Q = P."""
        b2, b4, b6 = 4*self.a2, 2*self.a4, 4*self.a6; b8 = 4*self.a2*self.a6 - self.a4**2
        poly = x**4 - b4*x**2 - 2*b6*x - b8 - (0 if P is None else P[0])*(4*x**3 + b2*x**2 + 2*b4*x + b6)
        if P is None:
            poly = 4*x**3 + b2*x**2 + 2*b4*x + b6          # 2Q = O: the points of order 2
        out = []
        for X in self.roots_K(poly):
            for Q in self.points_with_x(X):
                if self.equal(self.add(Q, Q), P): out.append(Q)
        return out
    def equal(self, P, Q):
        if P is None or Q is None: return P is None and Q is None
        return sp.expand(P[0] - Q[0]) == 0 and sp.expand(P[1] - Q[1]) == 0


def torsion(E):
    pts = [None]; frontier = [None]
    while frontier:
        new = []
        for P in frontier:
            for Q in E.halves(P):
                if not any(E.equal(Q, R) for R in pts): pts.append(Q); new.append(Q)
        frontier = new
    for P in list(pts):
        for Q in list(pts):
            R = E.add(P, Q); assert any(E.equal(R, T) for T in pts), 'not closed under addition'
    return pts


def count_mod(a2, a4, a6, p):
    a2, a4, a6 = int(a2 % p), int(a4 % p), int(a6 % p)
    def leg(v):
        v %= p
        return 0 if v == 0 else (1 if pow(v, (p - 1)//2, p) == 1 else -1)
    return p + 1 + sum(leg(X**3 + a2*X*X + a4*X + a6) for X in range(p))


def torsion_bound(a2, a4, a6, disc):
    from math import gcd
    g = 0
    den = sp.lcm([sp.Rational(c).q for c in (a2, a4, a6)])
    for p in sp.primerange(13, 400):
        if p % 4 != 1 or disc % p == 0 or den % p == 0: continue
        A2, A4, A6 = [int(sp.Rational(c).p * pow(int(sp.Rational(c).q), -1, p)) % p for c in (a2, a4, a6)]
        g = gcd(g, count_mod(A2, A4, A6, p))
    return g


def five_forms(m, n, t):
    t = Fr(t); S = Fr(m*m + n*n)
    F0 = m*m + n*n*t*t; F8 = n*n + m*m*t*t; F4 = S*(1 + t*t)/2; L = F4 - 2*m*n*t; U = F4 + 2*m*n*t
    def sq(q):
        from math import isqrt
        q = Fr(q); return q >= 0 and isqrt(q.numerator)**2 == q.numerator and isqrt(q.denominator)**2 == q.denominator
    if not all(sq(v) for v in (F0, F8, F4, L, U)): return 'no'
    cells = [(m*t + n)**2, (m*t - n)**2, (n*t + m)**2, (n*t - m)**2, F0, F8, F4, L, U]
    return 'SOLUTION' if len(set(cells)) == 9 and all(c != 0 for c in cells) else 'degenerate'


def check(m, n, rec):
    name = rec['curve']; Q = quartic(m, n, name)
    # 1. the class at the degenerate solutions
    d = int(rec['d'])
    for s0 in degenerate_s(m, n):
        if isinstance(s0, str):
            P = sp.Poly(Q, s); assert P.degree() == 4
            val = P.LC() if s0 == 'inf1' else Q.subs(s, 0)          # t = infinity: the leading coefficient; t = 0: Q(0)
            assert is_gauss_square(d * val), ('class at t = 0 or infinity', name, s0)
        else:
            assert is_gauss_square(d * Q.subs(s, s0)), ('class at a degenerate point', name, s0)
    roots = [r for r in sp.Poly(Q, s).ground_roots() if r.is_Rational]
    assert roots, 'no rational root'
    r = roots[0]
    A, (a2, a4, a6) = weierstrass(sp.expand(d*Q), r)
    # 2. the record's curve of rank 0 is E or its twist by -1
    tw = int(rec['rank0_twist'])
    mine = (0, a2, 0, a4, a6) if tw == 1 else (0, -a2, 0, a4, -a6)
    assert isomorphic(mine, [sp.Rational(v) for v in rec['rank0_curve']]), 'rank-0 curve is not E or E^(-1)'
    # 3. torsion over Q(i)
    E = EK(sp.Rational(a2), sp.Rational(a4), sp.Rational(a6))
    T = torsion(E)
    b2, b4, b6 = 4*a2, 2*a4, 4*a6; b8 = 4*a2*a6 - a4*a4
    disc = sp.Rational(-b2*b2*b8 - 8*b4**3 - 27*b6*b6 + 9*b2*b4*b6)
    bound = torsion_bound(a2, a4, a6, disc.p * disc.q)
    assert bound == len(T), ('torsion', len(T), bound)
    assert len(T) == int(rec['torsion_K']), ('torsion count differs from the record', len(T), rec['torsion_K'])
    # 4. the equations
    g = -A*x / (2*r*x + A)
    yy = sp.symbols('yy')
    cand = set()
    for Tp in T:
        if Tp is None:
            # x(+-P) = g(x): x (2 r x + A) + A x = 0
            for X in E.roots_K(sp.expand(x*(2*r*x + A) + A*x)): cand.add(X)
            continue
        Tq = Tp if tw == -1 else E.neg(Tp)     # sigma P = P + T  (twist of rank 0)  or  -P + T, x(-P + T) = x(P - T)
        xt, yt = Tq
        # x(P + Tq) = ((yy - yt)/(x - xt))^2 - a2 - x - xt, with yy^2 = f(x); set equal to g, clear denominators
        num = sp.expand(((yy - yt)**2 - (a2 + x + xt)*(x - xt)**2) * (2*r*x + A) + A*x*(x - xt)**2)
        num = sp.expand(num.subs(yy**2, E.f(x)))
        c1 = sp.expand(num.coeff(yy, 1)); c0 = sp.expand(num.coeff(yy, 0))
        poly = sp.expand(c0**2 - c1**2 * E.f(x))     # c0 + c1 yy = 0, squared
        if poly == 0: raise AssertionError('identity for a torsion point')
        for X in E.roots_K(poly): cand.add(X)
        cand.add(xt)
    cand.add(sp.Rational(-A, 2*r) if r != 0 else None)
    ts = set()
    for X in cand:
        if X is None or X == 0: continue
        s0 = sp.nsimplify(sp.expand(r + A/X)); t0 = sp.nsimplify(sp.expand(I*s0))
        if t0.is_Rational: ts.add(Fr(int(t0.p), int(t0.q)))
    ts |= {Fr(0)}
    verdict = {str(t): five_forms(m, n, t) for t in ts}
    return dict(curve=name, d=d, torsion=len(T), ncand=len(cand), t=verdict)


if __name__ == '__main__':
    OPEN = [(195, 28), (223, 167), (609, 200), (809, 409), (861, 620)]
    ok = 0
    for l in open(os.path.join(here, 'second_123.jsonl')):
        r = json.loads(l); p = tuple(r['pair'])
        if p not in OPEN: continue
        rec = [c for c in r['curves'] if c.get('status') == 'closed'][0]
        out = check(p[0], p[1], rec)
        bad = [t for t, v in out['t'].items() if v == 'SOLUTION']
        assert not bad
        ok += 1
        print('%d,%d curve %s: d = %d, torsion over Q(i) %d (equal to the bound from #E(F_p)), %d candidates; rational t: %s'
              % (p[0], p[1], out['curve'], out['d'], out['torsion'], out['ncand'],
                 ', '.join('%s %s' % (t, v) for t, v in sorted(out['t'].items()))), flush=True)
    # (440,279), closed with its partner (719,161) by an L-value (SECOND_DESCENT.md, section 6): the curve L1 L3 L5 L6,
    # its twist by -1 of rank 0 by L(E,1) != 0 (lproof_440_279.txt); the rest is rechecked as above.
    for l in open(os.path.join(here, 'second_123.jsonl')):
        r = json.loads(l)
        if tuple(r['pair']) != (440, 279): continue
        c = [c for c in r['curves'] if c['curve'] == '1356'][0]
        Q = quartic(440, 279, '1356'); rt = [v for v in sp.Poly(Q, s).ground_roots() if v.is_Rational][0]
        A, (a2, a4, a6) = weierstrass(sp.expand(int(c['d'])*Q), rt)
        # the curve whose L-value is summed in lproof_440_279.txt: y^2 = x (x - 719^2)(x + 161^2)
        lcurve = [0, 161**2 - 719**2, 0, -(719*161)**2, 0]
        assert isomorphic((0, -a2, 0, a4, -a6), lcurve), 'the L-value curve is not the twist by -1'
        rec = dict(curve='1356', d=c['d'], rank0_twist=-1, rank0_curve=lcurve, torsion_K=None)
        E = EK(sp.Rational(a2), sp.Rational(a4), sp.Rational(a6)); rec['torsion_K'] = len(torsion(E))
        out = check(440, 279, rec)
        assert not [t for t, v in out['t'].items() if v == 'SOLUTION']
        ok += 1
        print('440,279 curve 1356 (twist by -1 of rank 0 by an L-value): d = %d, torsion over Q(i) %d (equal to the bound), '
              '%d candidates; rational t: %s' % (out['d'], out['torsion'], out['ncand'],
              ', '.join('%s %s' % (t, v) for t, v in sorted(out['t'].items()))))
    print('second_exact OK: %d closures' % ok)
