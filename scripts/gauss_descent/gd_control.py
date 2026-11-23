"""Positive control for gauss_descent.py: twisted systems that have a solution by construction.
For random coprime (m, n) and (u, v), let c_j be the squarefree part of the norm of the j-th form at (u, v).  The system
"N(form_j) is c_j times a square" has the solution (u, v), so the descent set computed with the general rule for the
sets D_j must contain the tuple of that solution: it must never come out empty, and the tuple must survive every prime."""
import sys, random, itertools
import numpy as np
from math import gcd
import gauss_descent as gd

def sqfree(n):
    f = gd.factor(n); r = 1
    for q, e in f.items():
        if e % 2: r *= q
    return r

def dset_general(c, res):
    """Squarefree delta, up to squares, with N(delta) = c times a square, supported above 2, the primes of c and of res."""
    fc = gd.factor(c); primes = sorted(set(q for q in list(fc) + list(gd.factor(res)) if q != 2))
    opts = []
    for q in primes:
        odd = fc.get(q, 0) % 2 == 1
        if q % 4 == 3:
            if odd: return []                       # the norm of a Gaussian integer has even valuation at an inert prime
            opts.append([(1, 0), (q, 0)] if res % q == 0 else [(1, 0)])
        else:
            g = gd.gprime(q)
            if odd: opts.append([g, (g[0], -g[1])])
            else: opts.append([(1, 0), (q, 0)] if res % q == 0 else [(1, 0)])
    b = fc.get(2, 0) % 2; out = []
    for a in (0, 1):
        for ch in itertools.product(*opts):
            x = (0, 1) if a else (1, 0)
            if b: x = gd.gmul(x, (1, 1))
            for y in ch: x = gd.gmul(x, y)
            out.append(x)
    return out

def resultants(m, n):
    return [2*m*n, 2*m*n, 2, 2*(m*m - n*n), 2*(m*m - n*n)]

def run(m, n, cs, primes, sol=None, kmax=5):
    D = [dset_general(c, r) for c, r in zip(cs, resultants(m, n))]
    if any(len(d) == 0 for d in D): return D, None, None
    T = np.ones(tuple(len(d) for d in D), dtype=bool); M = np.ones(T.shape, dtype=bool)
    for p in primes:
        k = kmax if p > 50 else max(kmax, 8 if p == 2 else 6)
        L, W = gd.local_set(m, n, p, k); cd = [[L.cls(x) for x in d] for d in D]
        allowed = np.zeros(T.shape, dtype=bool)
        for w in W:
            idx = [[i for i, c in enumerate(cd[j]) if gd.match(c, w[j])] for j in range(5)]
            if all(idx): allowed[np.ix_(*idx)] = True
        T &= allowed
        if sol:
            c = [L.cls(gd.value(f, *sol)) for f in gd.forms(m, n)]
            ind = np.zeros(T.shape, dtype=bool); idx = [[i for i, cc in enumerate(cd[j]) if cc == c[j]] for j in range(5)]
            if all(idx): ind[np.ix_(*idx)] = True
            M &= ind
    return D, T, M

if __name__ == '__main__':
    random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1); ntest = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    primes = [q for q in gd.PR if q <= 60]; done = 0; fails = 0; sizes = []
    # the general rule reproduces the sets of gauss_descent.py for the real system
    for (m, n) in ((820, 799), (889, 320), (72, 65), (15, 8), (137, 7)):
        S = m*m + n*n; A = [sorted(d) for d in gd.dsets(m, n)]; B = [sorted(dset_general(c, r)) for c, r in zip([1, 1, sqfree(2*S), 2, 2], resultants(m, n))]
        assert A == B, (m, n)
    print('the general rule gives the sets D_j of gauss_descent.py for 5 ratios')
    while done < ntest:
        m = random.randint(2, 60); n = random.randint(1, m - 1); u = random.randint(1, 40); v = random.randint(1, 40)
        if gcd(m, n) != 1 or gcd(u, v) != 1: continue
        vals = [gd.value(f, u, v) for f in gd.forms(m, n)]
        if any(x == (0, 0) for x in vals): continue
        cs = [sqfree(gd.gnorm(x)) for x in vals]
        bad = sorted(set(gd.factor(2*m*n*(m*m - n*n)*(m*m + n*n))) | set(q for c in cs for q in gd.factor(c)))
        if max(bad) > 3000 or any(len(gd.factor(c)) > 5 for c in cs): continue
        D, T, M = run(m, n, cs, sorted(set(primes + bad)), sol=(u, v))
        done += 1
        if T is None or not (T & M).any():
            fails += 1; print('FAIL', m, n, u, v, cs, None if T is None else (int(T.sum()), int(M.sum())))
        else: sizes.append((int(np.prod(T.shape)), int(T.sum())))
    print('%d twisted systems with a known solution: %d failures; nontrivial pi/pibar choices in %d of them; tuples left: min %d, max %d' % (done, fails, sum(1 for s in sizes if s[0] > 0), min(s[1] for s in sizes), max(s[1] for s in sizes)))
