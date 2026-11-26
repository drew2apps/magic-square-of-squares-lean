"""Port of scripts/second_descent/second_descent.py (and wide.py's first-descent check) from Sage to PARI.
EXPERIMENT, written by an AI model (Claude), not read by a person.  Same family of 24 quartics, same classes, same
rank bounds (PARI ellrank, effort 1, least over the isogeny class), torsion over K = Q(i) by PARI elltors over nfinit(y^2+1),
roots over K by nfroots.  usage: sd_pari.py m,n [m,n ...]   (one JSON line per ratio)"""
import sys, json, os
from math import isqrt, gcd
from itertools import combinations
from fractions import Fraction as Fr
REPO = '/tmp/claude-0/-home-user/9ec7a8d3-42dc-53c8-8518-1a663ac809f1/scratchpad/repo/scripts/'
sys.path.insert(0, REPO + 'gauss_descent')
import cypari2
pari = cypari2.Pari(); pari.allocatemem(10**9)
nf = pari('nfinit(y^2+1)')
I = pari('Mod(y, y^2+1)')
s = pari('s'); X = pari('x')

def issq(q):
    q = Fr(q)
    return q >= 0 and isqrt(q.numerator)**2 == q.numerator and isqrt(q.denominator)**2 == q.denominator
def classify_t(m, n, t):
    h = Fr(m*m + n*n, 2)
    F0 = m*m + n*n*t*t; F8 = n*n + m*m*t*t; F4 = h*(1 + t*t); L = F4 - 2*m*n*t; U = F4 + 2*m*n*t
    if not all(issq(v) for v in (F0, F4, F8, L, U)): return 'no'
    cells = [(m*t + n)**2, (m*t - n)**2, (n*t + m)**2, (n*t - m)**2, F0, F8, F4, L, U]
    return 'SOLUTION' if len(set(cells)) == 9 and all(c != 0 for c in cells) else 'degenerate'
def fr(z): z = pari(z); return Fr(int(pari.numerator(z)), int(pari.denominator(z)))
def ksq(g):
    """Is g (an element of K, as Mod(a + b y, y^2+1) or rational) a square in K?"""
    g = pari.lift(pari(g) * I**0) if True else g
    return len(pari.nfroots(nf, X**2 - g)) > 0
def knorm(g):
    g = pari.lift(pari(g) + 0*I)
    a = pari.polcoef(g, 0, 'y'); b = pari.polcoef(g, 1, 'y'); return fr(a*a + b*b)

def bounds(ainv):
    E = pari.ellinit(ainv); hi = 10**9
    for F in pari.ellisomat(pari.ellminimalmodel(E), 0, 1)[0]:
        hi = min(hi, int(pari.ellrank(pari.ellinit(F), 1)[1]))
    return hi

def family(m, n):
    D = m*m - n*n; S = m*m + n*n
    L = {1: (m - n*s, Fr(m, n)), 2: (m + n*s, Fr(-m, n)), 3: (n + m*s, Fr(-n, m)), 4: (n - m*s, Fr(n, m)), 5: (1 - s, Fr(1)), 6: (1 + s, Fr(-1))}
    quad = {'+': S*s**2 - 2*D*s + S, '-': S*s**2 + 2*D*s + S}
    swap = {1: 2, 2: 1, 3: 4, 4: 3, 5: 6, 6: 5, '+': '-', '-': '+'}
    out = {}; seen = set()
    cands = [tuple(c) for c in combinations(range(1, 7), 4)] + [(sg,) + tuple(c) for sg in '+-' for c in combinations(range(1, 7), 2)]
    first = [(1, 3, 5, 6), ('+', 5, 6), ('+', 1, 4), ('+', 1, 3), ('+', 2, 4), ('+', 2, 3), (1, 4, 5, 6)]
    for key in first + [c for c in cands if c not in first]:
        norm = tuple(sorted(key, key=str)); img = tuple(sorted((swap[k] for k in key), key=str))
        if norm in seen or img in seen: continue
        seen.add(norm); Q = pari(1); r = None
        for k in key:
            if k in quad: Q = Q*quad[k]
            else:
                Q = Q*L[k][0]
                if r is None: r = L[k][1]
        out[''.join(str(k) for k in key)] = (Q, r)
    return out

def rational_class(g, m, n):
    primes = sorted(set(int(p) for p in pari.factor(2*m*n*(m*m - n*n)*(m*m + n*n))[0]))
    N = knorm(g); cands = [1]
    for p in primes:
        if N.numerator % p == 0 or N.denominator % p == 0 or p == 2: cands = cands + [c*p for c in cands]
    for d in cands:
        if ksq(pari(g) / d): return d
    return None

def weierstrass(Q, r):
    """Y^2 = Q(s), Q(r) = 0; s = r + 1/w: cubic A w^3 + B w^2 + C w + E; model [0, B, 0, A C, A^2 E], s = r + A/x."""
    rr = pari(r.numerator) / r.denominator
    c = []; D = Q
    for k in range(1, 5):
        D = pari.deriv(D, 's'); c.append(fr(pari.subst(D, 's', rr)) / [1, 2, 6, 24][k-1])
    A, B, C, E0 = c[0], c[1], c[2], c[3]     # coefficient of w^(4-k) is Q^(k)(r)/k!
    ain = [Fr(0), B, Fr(0), A*C, A*A*E0]
    return ain, A

def try_curve(m, n, name, Q, r):
    out = dict(curve=name)
    Ssq = isqrt(m*m + n*n)**2 == m*m + n*n
    if Ssq: vals = (pari.subst(Q, 's', -I), pari.subst(Q, 's', I))
    else: vals = (pari.subst(Q, 's', 0), pari.pollead(Q))
    cls = [rational_class(v, m, n) for v in vals]
    if cls[0] is None or cls[0] != cls[1]: out['status'] = 'class not rational or not the same at the two degenerate solutions'; return out
    d = cls[0]; out['d'] = d
    ain, A = weierstrass(d*Q, r); ainm, _ = weierstrass(-d*Q, r)
    def pint(a): return [pari(x.numerator) / x.denominator for x in a]
    rp, rm = bounds(pint(ain)), bounds(pint(ainm)); out['rank_bounds'] = [rp, rm]
    if rp != 0 and rm != 0: out['status'] = 'no summand of rank 0'; return out
    out['rank0_twist'] = 1 if rp == 0 else -1
    EK = pari.ellinit(pint(ain), nf)
    T = pari.elltors(EK); tors = [pari('[0]')]
    for g, k in zip(list(T[2]), [int(v) for v in T[1]]):
        tors = [pari.elladd(EK, P, pari.ellmul(EK, g, j)) for P in tors for j in range(k)]
    out['torsion_K'] = len(tors)
    a2 = pint(ain)[1]; a4 = pint(ain)[3]; a6 = pint(ain)[4]
    f = X**3 + a2*X**2 + a4*X + a6
    Ap = pari(A.numerator) / A.denominator; rp_ = pari(r.numerator) / r.denominator
    cand = set(); identity = 0
    def add_x(xv):
        xv = pari.lift(pari(xv) + 0*I)
        if xv != 0: cand.add(str(pari.lift(rp_ + Ap / xv + 0*I)))
    for P in [rp_, 0, 1, -1]: cand.add(str(P))
    for Tp in tors:
        if len(Tp) == 1:
            add_x(-Ap / rp_); continue
        x0 = pari.lift(Tp[0] + 0*I); y0 = pari.lift(Tp[1] + 0*I); add_x(x0)
        Dn = 2*rp_*X + Ap
        hD = (-Ap*X + (a2 + X + x0)*Dn)*(X - x0)**2
        pol = (f*Dn - hD) if y0 == 0 else (((f + y0**2)*Dn - hD)**2 - 4*y0**2*f*Dn**2)
        pol = pari.lift(pol)
        if pol == 0: identity += 1; continue
        for rt in pari.nfroots(nf, pol): add_x(rt)
    if identity:
        out['identity_torsion'] = identity; out['status'] = 'equation identically satisfied'; return out
    ts = {}
    for sv in cand:
        t = pari.lift(pari('Mod(y, y^2+1)') * pari(sv))      # t = i s
        if pari.poldegree(t, 'y') >= 1 and pari.polcoef(t, 1, 'y') != 0: continue
        tq = fr(pari.polcoef(t, 0, 'y')); ts[str(tq)] = classify_t(m, n, tq)
    out['rational_t'] = dict(sorted(ts.items())); out['status'] = 'ALERT' if 'SOLUTION' in ts.values() else 'closed'
    return out

def close(m, n, stop=True):
    res = dict(pair=[m, n], curves=[], status='open')
    for name, (Q, r) in family(m, n).items():
        try: o = try_curve(m, n, name, Q, r)
        except Exception as ex: o = dict(curve=name, status='error ' + str(ex)[:80])
        res['curves'].append(o)
        if o['status'] in ('closed', 'ALERT'):
            res['status'] = o['status']; res.setdefault('by', name)
            if stop: break
    return res

def first(m, n):
    import scan
    c, last = scan.descent_fast(m, n, pmax=200)
    S = m*m + n*n
    sols = [(1, 1), (1, -1)] if isqrt(S)**2 == S else [(0, 1), (1, 0)]
    forms = [(I*n, m), (I*m, n), (1, I), (m + I*n, -n - I*m), (m + I*n, n + I*m)]
    vals = [[A*u + B*v for A, B in forms] for u, v in sols]
    differ = not all(ksq(a / b) for a, b in zip(*vals))
    return c, differ

if __name__ == '__main__':
    for a in sys.argv[1:]:
        m, n = map(int, a.split(','))
        print(json.dumps(close(m, n), default=str), flush=True)
