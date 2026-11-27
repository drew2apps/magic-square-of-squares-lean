"""Review z8 of roadmap H (written by an AI model; not read by a person; experiment, not counted).
Explicit 5-adic points at p = 5 (a prime of m^2 + n^2, not in the supports) for the known and extra tuples of the
true phantom ratios with 5 | m^2 + n^2: (u, v) mod 5^8 primitive and lam in {1, 2} with every value lam d_j l_j^+-(u, v)
of even valuation v < 7 and with a unit part that is a square mod 5; then every lift is a point.  Rests on
octic_cut.comps (iota mod 5^8).
usage: witness5.py [--write]   (witness5.txt; about 7 s)"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from z8common import emit, load_phantoms, TRUE_PHANTOMS
import octic_cut as O
p, K = 5, 8; Mo = p ** K; I = O.sqrtm1(p, K)
xs = np.arange(Mo, dtype=np.int64)
U = np.concatenate([xs, np.ones(Mo // p, dtype=np.int64)]); V = np.concatenate([np.ones(Mo, dtype=np.int64), p * np.arange(Mo // p, dtype=np.int64)])
def cls(vals):
    vals = vals % Mo; v = np.zeros(vals.shape, dtype=np.int64); w = vals.copy(); z = vals == 0
    for _ in range(K):
        d = (w % p == 0) & ~z; w = np.where(d, w // p, w); v = v + d
    ok = (~z) & (v % 2 == 0) & (v < K - 1); u = w % p
    return np.where(ok, np.where((u == 1) | (u == 4), 0, 1), 2)
R = load_phantoms(); nf = 0; nc = 0
out = ['z8: explicit 5-adic points (u, v mod 5^8, lam class) for the true phantom ratios with 5 | m^2 + n^2. Written by an AI model; not read by a person.']
for (m, n) in TRUE_PHANTOMS:
    if (m*m + n*n) % 5: continue
    assert m % 5 and n % 5 and (m*m - n*n) % 5
    known, extra = R[(m, n)]
    C = [cls((a % Mo) * U + (b % Mo) * V) for (_, _, a, b) in O.comps(m, n, I, Mo)]
    for kind, Ls in (('known', known), ('extra', extra)):
        for T in Ls:
            hit = None
            for lam in (0, 1):
                want = [(lam + O.chi(d, p)) % 2 for d in T]
                ok = np.ones(U.shape, dtype=bool)
                for k in range(8): ok &= C[k] == want[k // 2]
                idx = np.nonzero(ok)[0]
                if len(idx): hit = (int(U[idx[0]]), int(V[idx[0]]), lam); break
            if kind == 'extra': nc += 1; nf += hit is not None
            out.append('(%d,%d) %s %s: %s' % (m, n, kind, T, hit))
out.append('extra tuples: %d; with an explicit 5-adic point: %d' % (nc, nf))
emit(out, 'witness5.txt')
