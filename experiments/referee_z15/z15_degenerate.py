"""z15: the (a, b) classes of the degenerate solutions, derived independently.  For every coprime m > n >= 1, m <= 3000,
with m^2 + n^2 a square (solutions t = 1, -1) or twice a square (t = 0, infinity): evaluate the Gaussian forms of
t10/h1_pairs (a0 = m + i n t, a8 = n + i m t, aL = (m+n)(t-1) + i(m-n)(t+1), aU = (m-n)(t-1) + i(m+n)(t+1); at infinity
the leading coefficients), take the rational class of each product z = a_j a_k of square norm as sqf(2(Re z + |z|))
(z = d (p + iq)^2 gives 2(Re z + |z|) = 4 d p^2; the class is defined up to -1, a square in Q(i)), and compare
a = sqf(mn d0d8), b = sqf((m^2-n^2) dLdU) with the shapes S1 = (sqf(2mn), 1) and S2 = (1, sqf(2(m^2-n^2))), and the
products with those of tuple_descent.magic's known tuples.  Also the shapes of all residual tuples at t10's degenerate
controls.  Written by an AI model (reviewer z15), not read by a person."""
import os, sys, json
from math import gcd, isqrt
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
os.environ.setdefault('MSS_REPO', REPO)
sys.path.insert(0, os.path.join(REPO, 'experiments', 'roadmap-g', 'g1'))
import cypari2
P = cypari2.Pari()
def sqf(x): return int(P.core(abs(int(x))))
def gm(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def gclass(z):
    N = z[0] ** 2 + z[1] ** 2; r = isqrt(N); assert r * r == N, 'norm not a square'
    x = z[0] + r
    if x == 0: x = z[0] - r
    return sqf(2 * x)
def forms_at(m, n, t):
    if t == 'inf': return {'0': (0, n), '8': (0, m), 'L': (m + n, m - n), 'U': (m - n, m + n)}
    return {'0': (m, n * t), '8': (n, m * t), 'L': ((m + n) * (t - 1), (m - n) * (t + 1)), 'U': ((m - n) * (t - 1), (m + n) * (t + 1))}
def issq(x): r = isqrt(x); return r * r == x
print('z15_degenerate.py: written by an AI model (reviewer z15), not read by a person.')
cnt = {'S1': 0, 'S2': 0, 'other': 0}; agree = 0; ratios = 0; bad = []
import exact_local as X
TD = X.TD
for m in range(2, 3001):
    for n in range(1, m):
        if gcd(m, n) != 1: continue
        C = m * m + n * n; A = m * m - n * n
        ts = [1, -1] if issq(C) else ['inf', 0] if C % 2 == 0 and issq(C // 2) else None
        if ts is None: continue
        ratios += 1
        kn = TD.magic(m, n)[2]
        kprod = sorted((sqf(t[0] * t[1]), sqf(t[2] * t[3])) for t in kn)
        mine = []
        for t in ts:
            f = forms_at(m, n, t)
            p08 = gclass(gm(f['0'], f['8'])); pLU = gclass(gm(f['L'], f['U']))
            a = sqf(m * n * p08); b = sqf(A * pLU)
            s = 'S1' if (a, b) == (sqf(2 * m * n), 1) else 'S2' if (a, b) == (1, sqf(2 * A)) else 'other'
            cnt[s] += 1; mine.append((p08, pLU))
            if s != ('S1' if issq(C) else 'S2'): bad.append((m, n, t, s))
        agree += sorted(mine) == kprod
print('degenerate ratios m <= 3000:', ratios, '; solutions by shape:', cnt, '; wrong shape:', bad[:5])
print('products (d0d8, dLdU) of my evaluation equal those of tuple_descent.magic\'s known tuples at %d of %d ratios' % (agree, ratios))
# residual tuples at t10's degenerate controls
exec(open(os.path.join(HERE, 'z15_residual_all.py')).read().split('def read_list')[0].split('import cypari2')[1].split('TD.factorint = fact')[1].join(['', '']) if False else '')
sys.argv = [sys.argv[0]]
src = open(os.path.join(HERE, 'z15_residual_all.py')).read()
ns = {'__file__': os.path.join(HERE, 'z15_residual_all.py'), '__name__': 'z15res'}
exec(src.split('def read_list')[0], ns)
for m, n in [(4, 3), (7, 1), (12, 5), (17, 7), (627, 364), (991, 263)]:
    T, known = ns['residual'](m, n); A = m * m - n * n; out = []
    for t in T:
        a = sqf(m * n * t[0] * t[1]); b = sqf(A * t[2] * t[3])
        s = 'S1' if (a, b) == (sqf(2 * m * n), 1) else 'S2' if (a, b) == (1, sqf(2 * A)) else 'other'
        out.append('%s%s:%s' % (t, '*' if t in known else '', s))
    print('  (%d,%d) residual tuples (* = known): %s' % (m, n, ' '.join(out)))
