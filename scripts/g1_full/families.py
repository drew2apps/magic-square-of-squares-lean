#!/usr/bin/env python3
"""Sister families of G1. Two of the four pairs of opposite cells are squares identically
(roots ma+nb, na-mb and ma-nb, na+mb); the other five cells must be squares.
G1: base = the two edge pairs.   G2: base = the two corner pairs.   G3: one corner pair, one edge pair.
Forms are 4 x cell, as integer polynomials in (m, n, a, b), quadratic in each pair.

usage: families.py G1|G2|G3 PRIME_LIMIT N      (needs numpy)
Result for PRIME_LIMIT = 700, N = 200: G1 reproduces the 34 primes and the 60 survivors;
G2 and G3 have no obstructed class at any prime and no excluded pair."""
import math, sys, json
import numpy as np
def forms(fam, m2, mn, n2, a2, ab, b2, M):
    g4 = 2 * (m2 + n2) % M * ((a2 + b2) % M) % M
    d1 = (m2 - n2) % M * ((a2 - b2) % M) % M          # (m^2-n^2)(a^2-b^2) = 2 d1
    d2 = mn * ab % M                                    # mnab = d2 / 2
    if fam == 'G1': return [(4 * (m2 * b2 + n2 * a2)) % M, (4 * (n2 * b2 + m2 * a2)) % M, g4, (g4 - 8 * d2) % M, (g4 + 8 * d2) % M]
    if fam == 'G2': return [g4, (g4 + 4 * d1) % M, (g4 - 4 * d1) % M, (g4 + 16 * d2) % M, (g4 - 16 * d2) % M]
    if fam == 'G3': return [g4, (g4 + 16 * d2) % M, (g4 - 16 * d2) % M, (g4 + 2 * d1 + 24 * d2) % M, (g4 - 2 * d1 - 24 * d2) % M]
_c = {}
def setup(M, p):
    if (M, p) not in _c:
        sq = np.zeros(M, dtype=bool); x = np.arange(M, dtype=np.int64); sq[(x * x) % M] = True
        a = np.concatenate([np.ones(M, dtype=np.int64), np.arange(0, M, p, dtype=np.int64)])
        b = np.concatenate([np.arange(M, dtype=np.int64), np.ones(M // p, dtype=np.int64)])
        _c[(M, p)] = (sq, a * a % M, b * b % M, a * b % M)
    return _c[(M, p)]
def obstructed(fam, m, n, M, p):
    sq, a2, b2, ab = setup(M, p)
    ok = np.ones(len(a2), dtype=bool)
    for f in forms(fam, m * m % M, m * n % M, n * n % M, a2, ab, b2, M): ok &= sq[f]
    return not ok.any()
def primes(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]
if __name__ == '__main__':
    fam = sys.argv[1]; PL = int(sys.argv[2]); N = int(sys.argv[3])
    table = {}
    for p in primes(PL):
        if p == 2: continue
        cls = [r for r in range(p + 1) if obstructed(fam, r if r < p else 1, 1 if r < p else 0, p, p)]
        if cls: table[p] = cls
    print(fam, 'primes with an obstructed class (below', PL, '):', {p: len(c) for p, c in table.items()})
    def ratio(m, n, p): return p if n % p == 0 else m * pow(n, -1, p) % p
    rest = 0; closed = 0
    for m in range(2, N + 1):
        for n in range(1, m):
            if math.gcd(m, n) != 1: continue
            if any(ratio(m, n, p) in set(c) for p, c in table.items()) or any(obstructed(fam, m, n, 2 ** k, 2) for k in (5, 6, 7, 8)): closed += 1
            else: rest += 1
    print(' N =', N, 'excluded by a congruence:', closed, '| not excluded:', rest)
