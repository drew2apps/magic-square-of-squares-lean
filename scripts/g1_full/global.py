#!/usr/bin/env python3
"""G1 family: pairs with a trivial rational point. Try to prove that every rational solution of the
five-form system is degenerate, using elliptic quotients of rank 0.

With w = t^2:  f0 = m^2 + n^2 w, f8 = n^2 + m^2 w, g4 = 2(m^2+n^2)(w+1) (= 4 F4),
K = (m^2+n^2)^2 (w+1)^2 - 16 m^2 n^2 w (= 4L * 4U / 4).  A solution t gives a rational point with
first coordinate w = t^2 on y^2 = P(w) for each product P below (a factor w is allowed: multiply y by t).
If such a curve is an elliptic curve of rank 0, w is the coordinate of a torsion point."""
import json, math, sys
from fractions import Fraction as Fr
import cypari2
pari = cypari2.Pari(); pari.allocatemem(2 * 10**9)

def pmul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): r[i + j] += x * y
    return r
def issq_int(x):
    if x < 0: return False
    r = math.isqrt(x); return r * r == x
def issq(q):
    q = Fr(q); return issq_int(q.numerator) and issq_int(q.denominator)
def sqrt_fr(q):
    return Fr(math.isqrt(q.numerator), math.isqrt(q.denominator))

def curves(m, n):
    W = [0, 1]; f0 = [m*m, n*n]; f8 = [n*n, m*m]; S = m*m + n*n; g4 = [2*S, 2*S]
    K = [S*S, 2*S*S - 16*m*m*n*n, S*S]
    P = {'w*f0*g4': pmul(W, pmul(f0, g4)), 'w*f8*g4': pmul(W, pmul(f8, g4)), 'w*f0*f8': pmul(W, pmul(f0, f8)),
         'f0*g4*f8': pmul(f0, pmul(g4, f8)), 'f0*K': pmul(f0, K), 'g4*K': pmul(g4, K), 'f8*K': pmul(f8, K),
         'w*K': pmul(W, K),
         'w*f0*K': pmul(W, pmul(f0, K)), 'w*g4*K': pmul(W, pmul(g4, K)), 'w*f8*K': pmul(W, pmul(f8, K)),
         'w*f0*g4*f8': pmul(W, pmul(f0, pmul(g4, f8)))}
    return P

def torsion_w(poly):
    """poly: coefficients c0..c3 (cubic in w) or 0,c1..c4 (quartic with root w = 0).
    Returns (rank_lower, rank_upper, list of w-values of all torsion points; None means w = infinity)."""
    if len(poly) == 5:                      # w = 1/u: v^2 = c1 u^3 + c2 u^2 + c3 u + c4
        assert poly[0] == 0
        c = [poly[4], poly[3], poly[2], poly[1]]; inv = True
    else:
        c = poly; inv = False
    c0, c1, c2, c3 = c
    E = pari.ellinit([0, c2, 0, c1 * c3, c0 * c3 * c3])
    rk = pari.ellrank(E)
    lo, hi = int(rk[0]), int(rk[1])
    if lo == 0 and hi > 0:                 # the rank is the same for isogenous curves; the descent bound is not
        for item in pari.ellisomat(pari.ellminimalmodel(E), 0, 1)[0]:
            hi = min(hi, int(pari.ellrank(pari.ellinit(item))[1]))
            if hi == 0: break
    T = pari.elltors(E)
    orders = [int(x) for x in T[1]]; gens = list(T[2])
    pts = [pari('[0]')]
    for g, d in zip(gens, orders):
        pts = [pari.elladd(E, P, pari.ellmul(E, g, k)) for P in pts for k in range(d)]
    assert len(pts) == int(T[0])
    ws = []
    for P in pts:
        if len(P) == 1:                      # point at infinity of the cubic
            ws.append(Fr(0) if inv else None)
        else:
            X = Fr(int(pari.numerator(P[0])), int(pari.denominator(P[0])))
            x = X / c3
            if inv: ws.append(None if x == 0 else 1 / x)
            else: ws.append(x)
    return lo, hi, ws, int(T[0])

def classify_t(m, n, t):
    """'no' if the five forms are not all squares at t; 'degenerate' if they are but the nine cells are
    not distinct nonzero; 'SOLUTION' otherwise."""
    s = Fr(m*m + n*n, 2)
    F0 = m*m + n*n*t*t; F8 = n*n + m*m*t*t; F4 = s * (1 + t*t); L = F4 - 2*m*n*t; U = F4 + 2*m*n*t
    if not all(issq(v) for v in (F0, F8, F4, L, U)): return 'no'
    cells = [(m*t + n)**2, (m*t - n)**2, (n*t + m)**2, (n*t - m)**2, F0, F8, F4, L, U]
    if len(set(cells)) == 9 and all(c != 0 for c in cells): return 'SOLUTION'
    return 'degenerate'

def analyse(m, n):
    out = []
    for name, poly in curves(m, n).items():
        lo, hi, ws, T = torsion_w(poly)
        verdicts = set()
        if hi == 0:
            for w in ws:
                if w is None: verdicts.add('degenerate'); continue      # t = infinity
                if w < 0 or not issq(w): continue
                t = sqrt_fr(w)
                for tt in {t, -t}: verdicts.add(classify_t(m, n, tt))
        out.append(dict(curve=name, rank=[lo, hi], torsion=T, closes=(hi == 0 and 'SOLUTION' not in verdicts),
                        alert=('SOLUTION' in verdicts)))
    return out

if __name__ == '__main__':
    pairs = json.load(open(sys.argv[1]))['open']
    res = {}; closed = 0
    for (m, n) in pairs:
        r = analyse(m, n); res[f'{m},{n}'] = r
        cl = [x['curve'] for x in r if x['closes']]
        closed += bool(cl)
        print(m, n, 'CLOSED by' if cl else 'open  ', cl if cl else [tuple(x['rank']) for x in r],
              'ALERT' if any(x['alert'] for x in r) else '', flush=True)
    print(f'closed {closed} of {len(pairs)}')
    json.dump(res, open(sys.argv[2], 'w'))
