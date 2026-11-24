#!/usr/bin/env python3
"""The structure of the class d of Theorem D (EXPERIMENT; LAST_PAIR.md, section 17), tested on every coprime ratio in a range.

    (E)  an odd prime dividing d is = +-1 mod 8                         (Theorem E: proved by hand in section 17)
    (2)  2 divides d exactly once                                        (observed)
    (S)  (d/p) = +1 for every odd prime p dividing m^2 + n^2 to an odd power   (observed)

(2) and (S) are tested on the sets that cover_descent.py allows at the prime (upper bounds for the true local sets, so
"no class violates it" is proved for each ratio tested).  (E) is tested on the exhaustion of third/third.py at the odd
primes q = +-3 mod 8, q <= 61, dividing mn(m^2 - n^2): cover_descent.py is weaker than the exhaustion when q^2 divides
one of m, n, m - n, m + n, and so can the exhaustion at its depth (those pairs are skipped and counted).  usage: class_structure.py M0 M1"""
import sys, os
from math import gcd
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..')); sys.path.insert(0, os.path.join(here, '..', 'third'))
import cover_descent as cd
import third as bf
from sympy import factorint, legendre_symbol
M0, M1 = int(sys.argv[1]), int(sys.argv[2])
ratios = cE = vE = c2 = v2 = cS = vS = skipped = 0
for m in range(M0, M1 + 1):
    for n in range(1, m):
        if gcd(m, n) != 1: continue
        ratios += 1; S = m*m + n*n; P = bf.cover_primes(m, n); cl = bf.cands(P)
        def allowed(p):
            L = cd.Local(m, n, p, 14 if p == 2 else (8 if p < 50 else 5)); al = L.allowed()
            return cl if al == 'all' else [d for d in cl if L.dclass(d) in al]
        for d in allowed(2): c2 += 1; v2 += (d % 4 != 2)
        for q in P:
            if q % 8 in (3, 5) and q <= 61:
                N = bf.depth(q, m, n)
                if any(x % (q*q) == 0 for x in (m, n, m - n, m + n)): skipped += 1; continue      # the depth is too small to see it
                keep, aux = bf.cover_upper(m, n, q, N, cl)
                for j in np.nonzero(keep)[0]: cE += 1; vE += (cl[j] % q == 0)
        for p, e in factorint(S).items():
            if p > 2 and e % 2 == 1 and p < 400:
                for d in allowed(p): cS += 1; vS += (legendre_symbol(d, p) != 1)
print('coprime ratios with %d <= m <= %d: %d' % (M0, M1, ratios))
print('(E) classes kept by the exhaustion at an odd prime q = +-3 mod 8 (q <= 61) of mn(m^2 - n^2): %d; of these divisible by q: %d; pairs (ratio, q) skipped because q^2 divides one of m, n, m - n, m + n: %d' % (cE, vE, skipped))
print('(2) classes allowed at 2: %d; of these not = 2 mod 4: %d' % (c2, v2))
print('(S) classes allowed at an odd prime p < 400 dividing m^2 + n^2 to an odd power: %d; of these with (d/p) = -1: %d' % (cS, vS))
