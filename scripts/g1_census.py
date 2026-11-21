#!/usr/bin/env python3
"""Census of the G1 family for coprime 1 <= n < m <= N.

Local obstructions counted: (i) 2-adic (table modulo 32), (ii) a prime p = 5 mod 8 dividing
m^2 + n^2 to an odd power, (iii) (m:n) = +-2, +-1/2 modulo 7 or 23, or +-5, +-1/5 modulo 47.
With a second argument, pairs without one get a PARI rank bound for
E_B: Y^2 = X (X - (m^2-n^2)^2)(X - (m^2+n^2)^2)   (needs cypari2).

usage: g1_census.py N [rank]
"""
import math, sys, time
from collections import Counter
N = int(sys.argv[1]); DO_RANK = len(sys.argv) > 2
def table(M):
    sq = [False] * M
    for x in range(M): sq[x * x % M] = True
    prim = [(U, V) for U in range(M) for V in range(M) if math.gcd(math.gcd(U, V), M) == 1]
    return {(a, b): not any(sq[(a*a*V*V + b*b*U*U) % M] and sq[(2*(a*a+b*b)*(U*U+V*V)) % M]
                            and sq[(b*b*V*V + a*a*U*U) % M] for (U, V) in prim)
            for a in range(M) for b in range(M)}
T32 = table(32)
def factor(x):
    f = {}; d = 2
    while d * d <= x:
        while x % d == 0: f[d] = f.get(d, 0) + 1; x //= d
        d += 1
    if x > 1: f[x] = f.get(x, 0) + 1
    return f
def reasons(m, n):
    r = []
    if T32[(m % 32, n % 32)]: r.append('2-adic')
    if any(p % 8 == 5 and e % 2 == 1 for p, e in factor(m*m + n*n).items()): r.append('5 mod 8')
    for p, c in ((7, 4), (23, 4), (47, 25)):
        if ((m*m - c*n*n) * (c*m*m - n*n)) % p == 0 and not (m % p == 0 and n % p == 0):
            r.append(f'mod {p}')
    return r
pairs = [(m, n) for m in range(2, N + 1) for n in range(1, m) if math.gcd(m, n) == 1]
rows = [(m, n, reasons(m, n)) for (m, n) in pairs]
loc = [r for r in rows if r[2]]; rest = [r for r in rows if not r[2]]
print(f'N={N}: pairs {len(pairs)} | locally excluded {len(loc)} ({100*len(loc)/len(pairs):.1f}%)'
      f' | no obstruction of the three kinds {len(rest)}')
kinds = Counter(k for r in loc for k in set(r[2]))
print('  by kind (overlapping):', dict(sorted(kinds.items())))
if DO_RANK:
    import cypari2
    pari = cypari2.Pari(); pari.allocatemem(3 * 10**9)
    c = Counter(); t0 = time.time()
    for (m, n, _) in rest:
        a, b = (m*m - n*n) ** 2, (m*m + n*n) ** 2
        rk = pari.ellrank(pari.ellinit([0, -(a + b), 0, a * b, 0]))
        c[(int(rk[0]), int(rk[1]))] += 1
    print('  rank bounds (lower, upper):', dict(sorted(c.items())), f'[{time.time()-t0:.0f}s]')
    print(f'  rank 0 proved: {c[(0, 0)]} | remaining open: {len(rest) - c[(0, 0)]}')
