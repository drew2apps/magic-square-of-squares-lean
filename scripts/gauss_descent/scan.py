"""Every coprime pair m > n with M0 < m <= M1 that avoids the 484 residue classes (scripts/g1_full/modp_classes.txt) and
has m^2 + n^2 neither a square nor twice a square: the descent set of gauss_descent.py.  Prints the pairs whose set is
not empty, and a summary.  usage: scan.py M0 M1 [classes file [list file]]
(list file: one line per ratio that was run: m,n, the last prime used, the number of tuples left)"""
import sys, time
from math import gcd, isqrt
import numpy as np
import gauss_descent as gd

def small_factors(x, bound=20000):
    return set(p for p in gd.PR if p <= bound and x % p == 0)

def descent_fast(m, n, pmax=150, kmax=5):
    S = m*m + n*n; A = m**4 - 2*m**3*n + 2*m*m*n*n + 2*m*n**3 + n**4; B = m**4 + 2*m**3*n + 2*m*m*n*n - 2*m*n**3 + n**4
    bad = set(gd.factor(2*m*n)) | set(gd.factor(m*m - n*n)) | set(gd.factor(S)) | small_factors(A) | small_factors(B)
    D = gd.dsets(m, n); T = np.ones(tuple(len(d) for d in D), dtype=bool); last = 0
    for p in sorted(set([q for q in bad if q <= 20000] + [q for q in gd.PR if q <= pmax])):
        k = kmax if p > 50 else max(kmax, 8 if p == 2 else 6)
        L, W = gd.local_set(m, n, p, k); cd = [[L.cls(x) for x in d] for d in D]
        allowed = np.zeros(T.shape, dtype=bool)
        for w in W:
            idx = [[i for i, c in enumerate(cd[j]) if gd.match(c, w[j])] for j in range(5)]
            if all(idx): allowed[np.ix_(*idx)] = True
        T &= allowed; last = p
        if not T.any(): break
    return int(T.sum()), last

if __name__ == '__main__':
    M0, M1 = int(sys.argv[1]), int(sys.argv[2]); path = sys.argv[3] if len(sys.argv) > 3 else '../g1_full/modp_classes.txt'
    cl = {}
    for l in open(path):
        if l.startswith('#') or not l.strip(): continue
        a = l.split(); cl.setdefault(int(a[0]), set()).update(int(x) for x in a[1:])
    listf = open(sys.argv[4], 'w') if len(sys.argv) > 4 else None             # one line per ratio run: m,n and the last prime used
    t0 = time.time(); npass = ntriv = nempty = 0; left = []; hist = {}
    for m in range(M0 + 1, M1 + 1):
        for n in range(1, m):
            if gcd(m, n) != 1: continue
            if any(n % p and (m * pow(n, -1, p)) % p in rs for p, rs in cl.items()): continue
            npass += 1; S = m*m + n*n
            if isqrt(S)**2 == S or (S % 2 == 0 and isqrt(S//2)**2 == S//2): ntriv += 1; continue
            c, last = descent_fast(m, n)
            if c == 0: nempty += 1; hist[last] = hist.get(last, 0) + 1
            if listf: listf.write('%d,%d %d %d\n' % (m, n, last, c))
            else: left.append((m, n, c)); print('NOT EMPTY', m, n, c, flush=True)
    print('%d < m <= %d: %d pairs avoid the classes; %d with a trivial point (not run); %d without, of which %d have an empty descent set and %d do not; last prime used: %s  [%.0f s]'
          % (M0, M1, npass, ntriv, npass - ntriv, nempty, len(left), ' '.join('%d:%d' % kv for kv in sorted(hist.items())), time.time() - t0), flush=True)
