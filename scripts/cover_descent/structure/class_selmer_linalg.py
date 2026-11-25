#!/usr/bin/env python3
"""The criterion of COMPLETE.md, section 7, as linear algebra over F_2 (written after a review of 0b34613).

class_selmer.py lists the classes d = 2 * prod(T) by going through every subset T of the primes = +-1 (mod 8) of
mn(m^2 - n^2).  The conditions are linear in the exponents x_q of these primes: a condition "d has the residue symbol s
modulo r" reads  sum_q x_q [ (q/r) = -1 ] = [ (2/r) s = -1 ]  over F_2.  So the set of classes is empty or a coset of
the kernel of a matrix A, and its size is 2^(number of primes - rank A).  This program computes, for every coprime ratio
with m <= H, whether the system is solvable and the dimension of the kernel, and checks the size against class_selmer.py
for m <= H0.  The conditions are those of class_selmer.py (the function `empty` and the list `cond` of `criterion`).

usage: class_selmer_linalg.py H [H0] [--list FILE] [--from A]     (default H0 = min(H, 1000); with --list, the ratios
       with m <= H that are not of the two degenerate types and keep a class are written to FILE, one "m,n dim" per
       line; with --from, only the ratios with A < m <= H are counted)
"""
import sys, os
from math import gcd, isqrt, prod
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import class_selmer as C

def system(m, n):
    """None if some prime has no local point; otherwise (primes q, list of (r, s) conditions)."""
    if C.empty(m, n): return None
    X = m*n*(m*m - n*n); NW = (m*m + n*n)*(m**8 + 14*m**4*n**4 + n**8)
    P8 = [q for q in C.fac(X) if q > 2 and q % 8 in (1, 7)]
    cond = [(r, 1) for r in C.fac(NW) if r > 2]
    if C.L.exact5(m, n): cond.append((5, -1))
    for p in C.TP:
        if (X*NW) % p:
            tau = (m*pow(n, -1, p)) % p
            if tau in C.TAB[p][1]: cond.append((p, 1))
            elif tau in C.TAB[p][2]: cond.append((p, -1))
    return P8, cond

def solve(P8, cond):
    """(solvable, dimension of the kernel), by Gaussian elimination over F_2 on bit masks."""
    k = len(P8); rows = []
    for r, s in cond:
        mask = sum(1 << i for i, q in enumerate(P8) if C.chi(q, r) == -1)
        rhs = 1 if C.chi(2, r)*s == -1 else 0
        rows.append(mask | (rhs << k))
    rank = 0; piv = []
    for col in range(k):
        sel = next((i for i in range(rank, len(rows)) if rows[i] >> col & 1), None)
        if sel is None: continue
        rows[rank], rows[sel] = rows[sel], rows[rank]
        for i in range(len(rows)):
            if i != rank and rows[i] >> col & 1: rows[i] ^= rows[rank]
        rank += 1
    solvable = all(row != (1 << k) for row in rows[rank:])
    return solvable, k - rank

if __name__ == '__main__':
    argv = sys.argv[1:]; out = None; A = 1
    if '--list' in argv: i = argv.index('--list'); out = open(argv[i + 1], 'w'); del argv[i:i + 2]
    if '--from' in argv: i = argv.index('--from'); A = int(argv[i + 1]); del argv[i:i + 2]
    H = int(argv[0]); H0 = int(argv[1]) if len(argv) > 1 else min(H, 1000)
    checked = mismatch = 0; dims = {}
    for m in range(max(2, A + 1), H + 1):
        for n in range(1, m):
            if gcd(m, n) != 1: continue
            sy = system(m, n)
            if sy is None: size = 0; key = None
            else:
                ok, dim = solve(*sy); size = 2**dim if ok else 0
                N = m*m + n*n; deg = isqrt(N)**2 == N or isqrt(2*N)**2 == 2*N
                key = ('degenerate' if deg else 'other', dim if ok else None)
                dims[key] = dims.get(key, 0) + 1
                if out and not deg and ok: out.write('%d,%d %d\n' % (m, n, dim))
            if m <= H0:
                checked += 1; mismatch += size != len(C.criterion(m, n)[0])
    print('coprime ratios with m <= %d: the size 2^dim (or 0) against the classes listed by class_selmer.py: %d ratios, %d different' % (H0, checked, mismatch))
    print('ratios with %sm <= %d and a local point everywhere, by type and dimension of the kernel (None = no solution):' % ('%d < ' % A if A > 1 else '', H))
    for kind in ('degenerate', 'other'):
        row = sorted((d for t, d in dims if t == kind), key=lambda d: -1 if d is None else d)
        print('  %s: %s' % (kind, ', '.join('%s: %d' % ('none' if d is None else 'dim %d' % d, dims[(kind, d)]) for d in row)))
