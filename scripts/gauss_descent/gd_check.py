"""A second implementation of the local step of gauss_descent.py, written differently, as a cross-check.
For a prime p and a fixed precision k it runs through all primitive (u, v) modulo p^k (no recursion, no class labels)
and decides directly, for each form j and each delta in D_j, whether  lam * form_j(u, v) * delta  is a square at every
place above p, or cannot be decided at this precision (then it counts as possible).  A tuple is allowed at p if some
(u, v) and some rational unit lam allow all five of its entries.  The set of allowed tuples is a superset of the true one."""
import sys, numpy as np
from math import isqrt
import gauss_descent as gd

def is_square_local(x, p, k, r=None, mod=None):
    """True / False / None (undecided) for the Gaussian integer x, known modulo p^k, at all places above p."""
    a, b = x
    if p == 2:
        N = a*a + b*b
        if N == 0: return None
        e = (N & -N).bit_length() - 1
        if e + 5 > 2*k: return None
        if e % 2: return False
        for _ in range(e): a, b = (a + b)//2, (b - a)//2
        return (a % 8, b % 8) in gd.SQ8
    if p % 4 == 3:
        if a == 0 and b == 0: return None
        e = 0
        while a % p == 0 and b % p == 0:
            a //= p; b //= p; e += 1
            if e >= k: return None
        if e % 2: return False
        return pow((a*a + b*b) % p, (p - 1)//2, p) == 1
    res = True
    for rr in (r, mod - r):
        y = (a + b*rr) % mod
        if y == 0: return None
        e = 0
        while y % p == 0: y //= p; e += 1
        if e >= k: return None
        if e % 2 or pow(y % p, (p - 1)//2, p) != 1: res = False
    return res

def allowed_at(m, n, D, p, k):
    F = gd.forms(m, n); shape = tuple(len(d) for d in D); out = np.zeros(shape, dtype=bool)
    r = mod = None
    if p % 4 == 1: K = k + 6; mod = p**K; r = gd.sqrt_m1(p, K)
    lams = [1, 3, 5, 7] if p == 2 else ([1] if p % 4 == 3 else [1, next(g for g in range(2, p) if pow(g, (p - 1)//2, p) == p - 1)])
    pk = p**k
    pts = [(1, v) for v in range(pk)] + [(u, 1) for u in range(0, pk, p)]
    seen = set()
    for (u, v) in pts:
        vals = [gd.value(f, u, v) for f in F]
        for lam in lams:
            idx = []
            for j in range(5):
                x = (lam*vals[j][0], lam*vals[j][1])
                idx.append(tuple(i for i, d in enumerate(D[j]) if is_square_local(gd.gmul(x, d), p, k, r, mod) is not False))
            key = tuple(idx)
            if key in seen or not all(idx): continue
            seen.add(key); out[np.ix_(*[list(t) for t in idx])] = True
    return out

if __name__ == '__main__':
    prec = {2: 9, 3: 5, 5: 4, 7: 3}
    for s in sys.argv[1:]:
        m, n = (int(v) for v in s.split(','))
        D, T, hist, bad = gd.descent(m, n, 100)
        T2 = np.ones(T.shape, dtype=bool); T1 = np.ones(T.shape, dtype=bool); same = True
        for p, _ in hist:
            if p > 60: break
            T2 &= allowed_at(m, n, D, p, prec.get(p, 2))
        print('(%d,%d): first implementation %d tuples, second implementation (primes up to %d) %d tuples' % (m, n, int(T.sum()), min(60, hist[-1][0]), int(T2.sum())), flush=True)
