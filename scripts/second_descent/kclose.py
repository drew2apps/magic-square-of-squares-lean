#!/usr/bin/env python3
"""Second descent on curves over K = Q(i) that are not defined over Q.  EXPERIMENT.

The ten linear factors over K of the five forms and their conjugates, in s = -i t:
    L1 = m - n s, L2 = m + n s, L3 = n + m s, L4 = n - m s, L5 = 1 - s, L6 = 1 + s            (rational)
    M1 = (m+in) s - (m-in)  [gL],  M2 = (m-in) s + (m+in)  [conj gL],  M3 = (m+in) s + (m-in)  [gU],  M4 = (m-in) s - (m+in)  [conj gU]
For a set J of four of them, Q_J = their product.  At a solution the class of Q_J(s0) in K*/K*^2 is the class at one of the two
degenerate solutions (first descent), say delta_1 or delta_2.  So the solution is a K-point of delta_k Y^2 = Q_J(s) for k = 1
or 2.  Both curves have all their points of order 2 over K.  If the 2-Selmer rank over K (kdescent.py) is 2 for both, both
have rank 0 over K, the point of the solution is a torsion point, and the finitely many torsion points are tested directly.

usage: kclose.py m,n [m,n ...]"""
import sys, json, itertools
from math import isqrt, gcd
from fractions import Fraction as Fr
from sage.all__sagemath_schemes import *
import kdescent as kd
from second_descent import classify_t

K = QuadraticField(-1, 'i'); I = K.gen()
RK = PolynomialRing(K, 's'); s = RK.gen()

def factors(m, n):
    p = m + I*n; pb = m - I*n
    return {'L1': m - n*s, 'L2': m + n*s, 'L3': n + m*s, 'L4': n - m*s, 'L5': 1 - s, 'L6': 1 + s,
            'M1': p*s - pb, 'M2': pb*s + p, 'M3': p*s + pb, 'M4': pb*s - p}

def gint(z): return (int(z[0]), int(z[1]))

def curve(Q, r):
    """delta-twisted quartic Q over K with root r: E: y^2 = x^3 + B x^2 + A C x + A^2 E0, s = r + A/x; returns E, A and integral roots e_i (scaled)."""
    T = RK.gen(); cub = RK(sum(Q.derivative(k)(r)/factorial(k)*T**(4 - k) for k in range(1, 5)))
    A, B, Cc, E0 = cub[3], cub[2], cub[1], cub[0]
    return EllipticCurve(K, [0, B, 0, A*Cc, A*A*E0]), A

def sqfree_rep(g):
    """A Gaussian integer in the class of g in K*/K*^2 (small): strip square factors."""
    den = lcm(g[0].denominator(), g[1].denominator()); z = g*den*den
    a, b = int(z[0]), int(z[1]); out = K(1)
    fa = kd.pari('factor(%d + %d*I)' % (a, b))
    unit = z
    for q, e in zip(fa[0], fa[1]):
        qq = K(int(q.real())) + I*int(q.imag()); e = int(e)
        unit = unit/qq**e
        if e % 2: out *= qq
    return out*unit                                    # unit of Z[i] times squarefree part

def run(m, n, verbose=True):
    S = m*m + n*n; Ssq = isqrt(S)**2 == S; assert Ssq or (S % 2 == 0 and isqrt(S//2)**2 == S//2)
    F = factors(m, n); names = sorted(F); res = dict(pair=[m, n], closed_by=[], tried=0, bound0=0)
    roots = {k: -F[k][0]/F[k][1] for k in names}
    for J in itertools.combinations(names, 4):
        Q = prod(F[k] for k in J)
        vals = [Q(-I), Q(I)] if Ssq else [Q(0), Q.leading_coefficient()]
        if any(v == 0 for v in vals): continue
        deltas = []
        for v in vals:
            d = sqfree_rep(v)
            if not any((d/x).is_square() for x in deltas): deltas.append(d)
        r = roots[J[3]]; ok = True; cand = set(); info = []
        for d in deltas:
            # delta Y^2 = Q  <=>  (delta Y)^2 = delta Q
            E, A = curve(d*Q, r)
            xs = [A/(roots[k] - r) for k in J[:3]]
            den = lcm([x[j].denominator() for x in xs for j in (0, 1)]); e = [gint(x*den*den) for x in xs]
            b = kd.selmer_bound(e)
            info.append(b)
            if b != 0: ok = False; break
            for Tp in E.torsion_points():
                if Tp.is_zero(): cand.add(r); continue
                if Tp[0] != 0: cand.add(r + A/Tp[0])
        res['tried'] += 1
        if not ok: continue
        res['bound0'] += 1
        ts = {}
        for sv in cand:
            t = I*sv
            if t[1] != 0: continue
            tq = Fr(int(QQ(t[0]).numerator()), int(QQ(t[0]).denominator())); ts[str(tq)] = classify_t(m, n, tq)
        rec = dict(J=''.join(J), classes=len(deltas), rational_t=dict(sorted(ts.items())), status='ALERT' if 'SOLUTION' in ts.values() else 'closed')
        res['closed_by'].append(rec)
        if verbose: print('   ', (m, n), rec, flush=True)
    res['status'] = 'closed' if res['closed_by'] else 'open'
    return res

if __name__ == '__main__':
    for a in sys.argv[1:]:
        m, n = (int(v) for v in a.split(','))
        r = run(m, n); print(json.dumps(r), flush=True)
