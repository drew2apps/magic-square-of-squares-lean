#!/usr/bin/env python3
"""Second descent for the ratios with a degenerate solution.  EXPERIMENT, not reviewed.  See SECOND_DESCENT.md.

Ratio m:n with m^2 + n^2 = k^2 or 2 k^2.  With s = -i t (t = u/v), D = m^2 - n^2, S = m^2 + n^2, the factors
    L1 = m - n s = a0/v     L2 = m + n s = conj(a0)/v     L3 = n + m s = conj(a8)/v     L4 = n - m s = a8/v
    L5 = 1 - s,  L6 = 1 + s      (L5 L6 = N(g4)/v^2)
    q+ = S s^2 - 2 D s + S = -gL conj(gU)/v^2          q- = S s^2 + 2 D s + S = -conj(gL) gU/v^2
are the values of the five forms of GAUSS_DESCENT.md and their conjugates.  A quartic Q is a product of four of the L's or of
one q and two L's.  The first descent leaves only the classes of the two degenerate solutions; if Q has the same class d in
K*/K*^2 (K = Q(i)) at both, and d is rational, every solution gives a K-point (s0, Y) of  C: Y^2 = d Q(s), s0 = -i t0.
C is a curve over Q with a rational point (a root r of Q); E its Jacobian, x = A/(s - r).  If E(Q) or E^(-1)(Q) has rank 0,
then  sigma(P) = -P + T  or  P + T  with T torsion in E(K), and sigma(s0) = -s0 gives  x(P + T) = -A x/(2 r x + A)  for some
torsion point T: a polynomial equation for x = x(P), solved over K for every T.  Every rational t found is tested in the
five-form system.

usage: second_descent.py [--all] m,n [m,n ...]      one JSON line per ratio (--all: do not stop at the first curve that closes)
       second_descent.py --summary FILE [FILE ...]   summary of JSON lines"""
import sys, json
from math import isqrt, gcd
from itertools import combinations
from fractions import Fraction as Fr

def issq(q):
    q = Fr(q)
    if q < 0: return False
    a, b = q.numerator, q.denominator
    return isqrt(a)**2 == a and isqrt(b)**2 == b

def classify_t(m, n, t):
    """The five forms at t (same test as scripts/g1_full/dm.py)."""
    s = Fr(m*m + n*n, 2)
    F0 = m*m + n*n*t*t; F8 = n*n + m*m*t*t; F4 = s*(1 + t*t); L = F4 - 2*m*n*t; U = F4 + 2*m*n*t
    if not all(issq(v) for v in (F0, F4, F8, L, U)): return 'no'
    cells = [(m*t + n)**2, (m*t - n)**2, (n*t + m)**2, (n*t - m)**2, F0, F8, F4, L, U]
    if len(set(cells)) == 9 and all(c != 0 for c in cells): return 'SOLUTION'
    return 'degenerate'

def partner(m, n):
    a, b = m + n, m - n
    return (a//2, b//2) if a % 2 == 0 else (a, b)

def summary(files):
    R = {}
    for f in files:
        for l in open(f):
            j = json.loads(l); R[tuple(j['pair'])] = j
    closed = sorted(p for p in R if R[p]['status'] == 'closed')
    full = sorted(p for p in R if p in closed or partner(*p) in closed)
    sq = [p for p in R if isqrt(p[0]**2 + p[1]**2)**2 == p[0]**2 + p[1]**2]
    print('%d ratios: %d with m^2+n^2 a square, %d with twice a square' % (len(R), len(sq), len(R) - len(sq)))
    print('statuses:', ', '.join('%s %d' % (k, sum(r['status'] == k for r in R.values())) for k in sorted(set(r['status'] for r in R.values()))))
    print('closed directly: %d (square type %d, other type %d)' % (len(closed), sum(p in closed for p in sq), len(closed) - sum(p in closed for p in sq)))
    by = {}; ts = {}; deg = 0; rn = {}; tw = {}
    for p in closed:
        c = R[p]['curves'][-1]; by[c['curve']] = by.get(c['curve'], 0) + 1; rn[c['rank0_rootno']] = rn.get(c['rank0_rootno'], 0) + 1
        tw[c['rank0_twist']] = tw.get(c['rank0_twist'], 0) + 1
        for t, v in c['rational_t'].items(): ts[t + ': ' + v] = ts.get(t + ': ' + v, 0) + 1
        deg += bool(c.get('degenerate_found_by_equations'))
    print('closing curve:', ', '.join('%s %d' % kv for kv in sorted(by.items())))
    print('summand of rank 0: E(Q) in %d, E^(-1)(Q) in %d; root number of that curve:' % (tw.get(1, 0), tw.get(-1, 0)), ', '.join('%+d in %d' % kv for kv in sorted(rn.items())))
    print('rational t returned by the equations, with the result of the five-form test:', ', '.join('%s (%d)' % kv for kv in sorted(ts.items())))
    print('square type: t = 1 and t = -1 both returned by the equations in %d of %d closed' % (deg, sum(p in closed for p in sq)))
    print('closed, directly or through the partner: %d of %d' % (len(full), len(R)))
    print('not closed (%d):' % (len(R) - len(full)), ' '.join('%d,%d' % p for p in sorted(R) if p not in full))
    print('of these, partner outside the list (m > 1000, not run):', ' '.join('%d,%d' % p for p in sorted(R) if p not in full and partner(*p) not in R))

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == '--summary':
    summary(sys.argv[2:]); sys.exit(0)

from sage.all__sagemath_schemes import *
from cypari2 import Pari
pari = Pari()
K = QuadraticField(-1, 'i'); I = K.gen()
Rq = PolynomialRing(QQ, 's'); sq = Rq.gen()

def bounds(E):
    """Upper bound for the rank of E(Q): the least of PARI's ellrank upper bounds over the curves isogenous to E."""
    hi = 10**9
    for F in pari.ellisomat(E.minimal_model().pari_curve(), 0, 1)[0]:
        hi = min(hi, int(pari.ellrank(pari.ellinit(F), 1)[1]))
    return hi

def family(m, n):
    """The quartics, one from each orbit of s -> -s, each with a rational root r.  Key: the factors, e.g. '1356', '+14'."""
    D = m*m - n*n; S = m*m + n*n
    L = {1: (m - n*sq, QQ(m)/n), 2: (m + n*sq, QQ(-m)/n), 3: (n + m*sq, QQ(-n)/m), 4: (n - m*sq, QQ(n)/m), 5: (1 - sq, QQ(1)), 6: (1 + sq, QQ(-1))}
    quad = {'+': S*sq**2 - 2*D*sq + S, '-': S*sq**2 + 2*D*sq + S}
    swap = {1: 2, 2: 1, 3: 4, 4: 3, 5: 6, 6: 5, '+': '-', '-': '+'}
    out = {}; seen = set()
    cands = [tuple(c) for c in combinations(range(1, 7), 4)] + [(sg,) + tuple(c) for sg in '+-' for c in combinations(range(1, 7), 2)]
    first = [(1, 3, 5, 6), ('+', 5, 6), ('+', 1, 4), ('+', 1, 3), ('+', 2, 4), ('+', 2, 3), (1, 4, 5, 6)]
    for key in first + [c for c in cands if c not in first]:
        norm = tuple(sorted(key, key=str)); img = tuple(sorted((swap[k] for k in key), key=str))
        if norm in seen or img in seen: continue
        seen.add(norm); Q = 1; r = None
        for k in key:
            if k in quad: Q = Q*quad[k]
            else:
                Q = Q*L[k][0]
                if r is None: r = L[k][1]
        out[''.join(str(k) for k in key)] = (Q, r)
    return out

def rational_class(g, m, n):
    """A squarefree positive integer d with g/d a square in K, or None.  (-1 is a square in K, so the sign does not matter.)"""
    primes = sorted(set(p for p, _ in ZZ(2*m*n*(m*m - n*n)*(m*m + n*n)).factor()))
    N = g.norm(); cands = [QQ(1)]
    for p in primes:
        if N.valuation(p) != 0 or p == 2: cands = cands + [c*p for c in cands]
    for d in cands:
        if (g/d).is_square(): return d
    return None

def weierstrass(Q, r):
    """Y^2 = Q(s), Q(r) = 0.  With s = r + 1/w: (Y w^2)^2 = A w^3 + B w^2 + C w + E; x = A w, so s = r + A/x."""
    T = Rq.gen(); cub = Rq(sum(Q.derivative(k)(r)/factorial(k)*T**(4 - k) for k in range(1, 5)))
    A, B, Cc, E0 = cub[3], cub[2], cub[1], cub[0]
    return EllipticCurve(QQ, [0, B, 0, A*Cc, A*A*E0]), A

def try_curve(m, n, name, Q, r, d_given=None):
    out = dict(curve=name)
    Ssq = isqrt(m*m + n*n)**2 == m*m + n*n
    # the two degenerate solutions: t = 1, -1 (s = -i, i) or t = 0, infinity (s = 0, infinity)
    cls = [rational_class(K(v), m, n) for v in ((Q(-I), Q(I)) if Ssq else (Q(0), Q.leading_coefficient()))]
    if d_given is None and (cls[0] is None or cls[0] != cls[1]): out['status'] = 'class not rational or not the same at the two degenerate solutions'; return out
    d = cls[0] if d_given is None else QQ(d_given); out['d'] = str(d)
    E, A = weierstrass(d*Q, r); Em, _ = weierstrass(-d*Q, r)
    rp, rm = bounds(E), bounds(Em); out['rank_bounds'] = [rp, rm]
    if rp != 0 and rm != 0: out['status'] = 'no summand of rank 0'; return out
    E0 = E if rp == 0 else Em; out['rank0_curve'] = [str(c) for c in E0.a_invariants()]; out['rank0_twist'] = 1 if rp == 0 else -1
    out['rank0_rootno'] = int(pari(E0).ellrootno()); out['rank0_torsion'] = [int(c) for c in E0.torsion_subgroup().invariants()]
    EK = E.change_ring(K); R = PolynomialRing(K, 'x'); x = R.gen(); a2 = EK.a2(); f = x**3 + a2*x**2 + EK.a4()*x + EK.a6()
    tors = EK.torsion_points(); out['torsion_K'] = len(tors) + 1 if not any(T.is_zero() for T in tors) else len(tors)
    cand = set([K(r), K(0), K(1), K(-1)])                # s-values of special points (P = O; not needed otherwise, harmless)
    def add_x(xv):
        if xv != 0: cand.add(r + A/xv)                   # x = 0 is s = infinity, t = infinity: never a non-degenerate solution
    for T in tors:
        if T.is_zero():
            add_x(K(-A)/r)                               # x = g(x): 2 r x^2 + 2 A x = 0
            continue
        x0, y0 = T[0], T[1]; add_x(x0)                   # P = T or -T
        Dn = 2*r*x + A                                   # s(P') = -s(P)  <=>  x(P') = g(x) = -A x/(2 r x + A)
        hD = (-A*x + (a2 + x + x0)*Dn)*(x - x0)**2       # x(P + T) = lambda^2 - a2 - x - x0, lambda = (y - y0)/(x - x0)
        pol = (f*Dn - hD) if y0 == 0 else (((f + y0**2)*Dn - hD)**2 - 4*y0**2*f*Dn**2)
        if pol == 0: out['status'] = 'equation identically satisfied'; return out
        for rt, _ in pol.roots(): add_x(rt)
    ts = {}
    for s in cand:
        t = I*s
        if t[1] != 0: continue
        tq = Fr(int(QQ(t[0]).numerator()), int(QQ(t[0]).denominator())); ts[str(tq)] = classify_t(m, n, tq)
    out['rational_t'] = dict(sorted(ts.items())); out['status'] = 'ALERT' if 'SOLUTION' in ts.values() else 'closed'
    if Ssq: out['degenerate_found_by_equations'] = ('1' in ts and '-1' in ts)
    return out

def close(m, n, stop=True):
    S = m*m + n*n; assert gcd(m, n) == 1 and (isqrt(S)**2 == S or (S % 2 == 0 and isqrt(S//2)**2 == S//2))
    res = dict(pair=[m, n], curves=[], status='open')
    for name, (Q, r) in family(m, n).items():
        o = try_curve(m, n, name, Q, r); res['curves'].append(o)
        if o['status'] in ('closed', 'ALERT'):
            res['status'] = o['status']; res.setdefault('by', name)
            if stop: break
    return res

if __name__ == '__main__':
    args = sys.argv[1:]; stop = True
    if args and args[0] == '--all': stop = False; args = args[1:]
    for sarg in args:
        m, n = (int(v) for v in sarg.split(','))
        print(json.dumps(close(m, n, stop), default=str), flush=True)
