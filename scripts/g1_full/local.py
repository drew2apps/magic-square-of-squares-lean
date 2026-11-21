#!/usr/bin/env python3
"""G1 family, full five-form system: local obstructions for coprime 1 <= n < m <= N.
Forms (t = a/b): F0 = m^2 b^2 + n^2 a^2, F8 = n^2 b^2 + m^2 a^2, G4 = 2(m^2+n^2)(a^2+b^2),
GL = G4 - 8mnab, GU = G4 + 8mnab  (G4, GL, GU are 4 F4, 4 L, 4 U).
An obstruction modulo M = p^k: no (a:b) primitive at p with all five values squares modulo M."""
import math, os, sys, json
import numpy as np
from collections import Counter
N = int(sys.argv[1])
_cache = {}
def setup(M, p):
    if (M, p) in _cache: return _cache[(M, p)]
    sq = np.zeros(M, dtype=bool); x = np.arange(M, dtype=np.int64); sq[(x * x) % M] = True
    a = np.concatenate([np.ones(M, dtype=np.int64), np.arange(0, M, p, dtype=np.int64)])
    b = np.concatenate([np.arange(M, dtype=np.int64), np.ones(M // p, dtype=np.int64)])
    _cache[(M, p)] = (sq, a * a % M, b * b % M, a * b % M)
    return _cache[(M, p)]
def obstructed(m, n, M, p):
    sq, a2, b2, ab = setup(M, p)
    m2, n2, mn8 = m * m % M, n * n % M, 8 * m * n % M
    ok = sq[(m2 * b2 + n2 * a2) % M]
    ok &= sq[(n2 * b2 + m2 * a2) % M]
    g4 = 2 * (m2 + n2) % M * ((a2 + b2) % M) % M
    ok &= sq[g4]
    x = mn8 * ab % M
    ok &= sq[(g4 - x) % M]
    ok &= sq[(g4 + x) % M]
    return not ok.any()
def factor(x):
    f = {}; d = 2
    while d * d <= x:
        while x % d == 0: f[d] = f.get(d, 0) + 1; x //= d
        d += 1
    if x > 1: f[x] = f.get(x, 0) + 1
    return f
# mod-p table from classify.c (all primes < 9700; ratio r = m/n mod p, r = p means n = 0 mod p)
table = {}
for line in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'modp_classes.txt')):
    p, r = map(int, line.split()); table.setdefault(p, set()).add(r)
def ratio(m, n, p):
    return p if n % p == 0 else m * pow(n, -1, p) % p
def reasons(m, n):
    r = []
    for k in (5, 6, 7):
        if obstructed(m, n, 2 ** k, 2): r.append(f'2^{k}'); break
    for p in table:
        if not (m % p == 0 and n % p == 0) and ratio(m, n, p) in table[p]: r.append(f'p={p}')
    bad = set(factor(m)) | set(factor(n)) | set(factor(m * m - n * n)) | set(factor(m * m + n * n))
    for p in sorted(bad):
        if p == 2: continue
        for k in (2, 3):
            if p ** k <= 20000 and obstructed(m, n, p ** k, p): r.append(f'{p}^{k}'); break
    for p, e in factor(m * m + n * n).items():
        if p % 8 == 5 and e % 2 == 1 and not any(x.startswith(f'{p}^') for x in r): r.append(f'{p}:5mod8')
    return r
pairs = [(m, n) for m in range(2, N + 1) for n in range(1, m) if math.gcd(m, n) == 1]
rows = [(m, n, reasons(m, n)) for (m, n) in pairs]
loc = [x for x in rows if x[2]]; rest = [x for x in rows if not x[2]]
print(f'N={N}: pairs {len(pairs)} | locally obstructed (five forms) {len(loc)} | not obstructed {len(rest)}')
c = Counter()
for x in loc:
    for k in set('2-adic' if y.startswith('2^') else ('bad odd prime power' if '^' in y or '5mod8' in y else 'mod p') for y in x[2]): c[k] += 1
print('  by kind (overlapping):', dict(c))
pc = Counter(y for x in loc for y in x[2] if y.startswith('p='))
print('  mod-p primes used:', sorted(pc.items(), key=lambda z: int(z[0][2:])))
k2 = Counter(y for x in loc for y in x[2] if y.startswith('2^')); print('  2-adic first level:', dict(k2))
json.dump({'N': N, 'open': [(m, n) for (m, n, _) in rest], 'closed': {f'{m},{n}': r for (m, n, r) in loc}}, open(f'local5_{N}.json', 'w'))
if N <= 20: print('  not obstructed:', [(m, n) for (m, n, _) in rest])
