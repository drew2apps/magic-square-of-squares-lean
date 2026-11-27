"""Review z9 of roadmap H (run 7), written by an AI model; not read by a person.  Experiment, not counted.
Checks the degenerate counts of H2's exhaustive run (h2_ex_15000_*of3.txt, the 'degenerate=' field of each m) by a
separate route: the coprime m : n, 15000 < m <= 30000, 0 < n < m, with m^2 + n^2 a square (legs of a primitive
Pythagorean triple (a^2 - b^2, 2ab)) or twice a square ((m + n)/2, (m - n)/2 the legs of one), enumerated from the
parameters a > b, and each checked by isqrt.  Rests on: Euclid's parametrisation of primitive triples.
usage: python3 -I degcount.py [--write]   (degcount.txt)"""
import os, sys, collections
from math import gcd, isqrt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import *
O = Out('degcount.txt')
O('Review z9: degenerate ratios with 15000 < m <= 30000 from Pythagorean parameters, against H2. Written by an AI model; not read by a person.')
LO, HI = 15000, 30000
S = {}
for a in range(2, isqrt(2 * HI) + 3):                     # the hypotenuse a^2 + b^2 is at most sqrt(2) * 30000 < 42427
    for b in range(1, a):
        if (a - b) % 2 == 0 or gcd(a, b) != 1: continue
        x, y = a * a - b * b, 2 * a * b
        for (m, n), t in (((max(x, y), min(x, y)), 'square'), ((x + y, abs(x - y)), 'twice')):
            if LO < m <= HI and 0 < n < m and gcd(m, n) == 1: S[(m, n)] = t
for (m, n), t in S.items():
    N = m * m + n * n; assert (isqrt(N) ** 2 == N) if t == 'square' else (isqrt(N // 2) ** 2 * 2 == N)
h2 = collections.Counter()
for i in range(3):
    for l in open(os.path.join(H2, 'h2_ex_15000_%dof3.txt' % i)):
        if l.startswith('m='):
            w = dict(x.split('=') for x in l.split()[1:]); h2[(int(l.split()[0][2:]) - 1) // 1000 * 1000] += int(w['degenerate'])
mine = collections.Counter((m - 1) // 1000 * 1000 for m, n in S)
O('degenerate ratios: %d (square type %d, twice-square type %d); H2 total %d' % (len(S), sum(1 for t in S.values() if t == 'square'), sum(1 for t in S.values() if t == 'twice'), sum(h2.values())))
for b in sorted(mine): O('  %d-%d: here %d, H2 %d %s' % (b, b + 1000, mine[b], h2[b], 'equal' if mine[b] == h2[b] else 'DIFFERENT'))
O('per-m counts equal too: %s' % all(sum(1 for (m, n) in S if m == mm) == int(dict(x.split('=') for x in l.split()[1:])['degenerate'])
    for i in range(3) for l in open(os.path.join(H2, 'h2_ex_15000_%dof3.txt' % i)) if l.startswith('m=') for mm in [int(l.split()[0][2:])]))
O.finish()
