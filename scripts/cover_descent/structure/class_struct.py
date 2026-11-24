#!/usr/bin/env python3
"""What the theorems on the class (LAST_PAIR.md, section 17: E, F, G and the octic law O) decide by themselves.

For the pairs of G1.exceptions1000 (119 with m^2 + n^2 neither a square nor twice a square, 306 of those two types): the candidates are
d = 2 * (a product of primes = +-1 mod 8 dividing m n (m^2 - n^2))  (E, F); G asks that d be a square modulo every odd
prime of m^2 + n^2; O asks the same modulo every odd prime of m^8 + 14 m^4 n^4 + n^8, (primes of the octic below 10^5 only, and a
cofactor below 10^10; using fewer primes can only keep more candidates).  Python only; a few seconds.
"""
import sys, re, itertools, os
from math import isqrt
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..')
src = open(os.path.join(root, 'MagicSquareOfSquares', 'G1Census1000.lean')).read()
body = src[src.index('def exceptions1000'):]; body = body[:body.index(']')]
pairs = [tuple(int(x) for x in t) for t in re.findall(r'\((\d+), (\d+)\)', body)]
SMALL = [q for q in range(2, 100000) if all(q % r for r in range(2, isqrt(q) + 1))]
def triv(m, n):
    S = m*m + n*n; return isqrt(S)**2 == S or (S % 2 == 0 and isqrt(S//2)**2 == S//2)
def leg(a, p): return pow(a % p, (p - 1)//2, p)
def fac(x, full=True):
    out = []
    for q in SMALL:
        if q*q > x: break
        if x % q == 0:
            out.append(q)
            while x % q == 0: x //= q
    if x > 1 and (full or x < 10**10): out.append(x)
    return out
tot = nG = nGO = n2G = n2GO = 0; first = {}; deg = {}; named = {}
for (m, n) in pairs:
    tot += not triv(m, n)
    T = [q for q in fac(m*n*(m*m - n*n)) if q % 8 in (1, 7)]
    pS = [p for p in fac(m*m + n*n) if p != 2]
    pO = [p for p in fac(m**8 + 14*m**4*n**4 + n**8, full=False) if p != 2]
    cG = []; cGO = []
    for r in range(len(T) + 1):
        for S in itertools.combinations(T, r):
            d = 2
            for q in S: d *= q
            if all(leg(d, p) == 1 for p in pS):
                cG.append(d)
                if all(leg(d, p) == 1 for p in pO): cGO.append(d)
    if (m, n) in ((627, 364), (991, 263)): named[(m, n)] = cGO
    if triv(m, n): deg[len(cGO)] = deg.get(len(cGO), 0) + 1; continue
    nG += not cG; nGO += not cGO; n2G += 2 not in cG; n2GO += 2 not in cGO
    if 2 in cG and 2 not in cGO:
        p = min(p for p in pO if leg(2, p) != 1); first[p] = first.get(p, 0) + 1
print('pairs with m^2 + n^2 neither a square nor twice a square:', tot)
print('class 2 excluded by G: %d; by G or O: %d; smallest prime of the octic that excludes it (prime: pairs): %s' % (n2G, n2GO, sorted(first.items())))
print('no candidate passes E, F, G: %d; no candidate passes E, F, G, O: %d' % (nG, nGO))
print('the %d pairs of the two degenerate types (a solution exists, so a candidate must pass): number of candidates that pass E, F, G, O (number: pairs): %s' % (sum(deg.values()), sorted(deg.items())))
print('candidates that pass for (627,364): %s; for (991,263): %s' % (named[(627, 364)], named[(991, 263)]))
