#!/usr/bin/env python3
"""Review x4 (run 6): the congruence lists (Result 1b of G1_FULL_SYSTEM.md: the 484 classes of
scripts/g1_full/modp_classes.txt and the 5 (mod 8) condition) regenerated for 1000 < m <= 11000 by separate code (numpy,
prime by prime) and compared with the six committed lists that census_table.py reads.  Prints only.  About 45 s.
Written by an AI model; not read by a person.   python3 congr.py"""
import os, sys, json, math, collections, numpy as np
R = os.environ.get('MSS_REPO', os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
A, B = 1000, 11000
cls = collections.defaultdict(set)
for line in open(R + '/scripts/g1_full/modp_classes.txt'):
    p, r = line.split()[:2]; cls[int(p)].add(int(r))
order = sorted(cls, key=lambda p: -len(cls[p]) / p)
tabs = {}
for p in order:
    bad = np.zeros(p, dtype=bool); bad[list(cls[p])] = True
    inv = np.zeros(p, dtype=np.int64); inv[1:] = [pow(i, -1, p) for i in range(1, p)]
    tabs[p] = (bad, inv)
def five_odd(N):           # some prime = 5 (mod 8) divides N to an odd power (trial division, my own)
    x = N
    for d in [2] + list(range(3, math.isqrt(N) + 2, 2)):
        if d * d > x: break
        if x % d == 0:
            k = 0
            while x % d == 0: x //= d; k += 1
            if d % 8 == 5 and k % 2: return True
    return x > 1 and x % 8 == 5
out = []
for m0 in range(A + 1, B + 1, 250):
    ms = []; ns = []
    for m in range(m0, min(m0 + 250, B + 1)):
        n = np.arange(1, m, dtype=np.int64); n = n[np.gcd(n, m) == 1]
        ms.append(np.full(len(n), m, dtype=np.int64)); ns.append(n)
    M = np.concatenate(ms); N = np.concatenate(ns)
    for p in order:
        bad, inv = tabs[p]
        np_ = N % p; t = ((M % p) * inv[np_]) % p
        keep = ~((np_ != 0) & bad[t]); M = M[keep]; N = N[keep]
    for m, n in zip(M.tolist(), N.tolist()):
        if not five_odd(m * m + n * n): out.append((m, n))
out.sort()
E = R + '/experiments'; L = set()
for f in [E + '/census2000/p2000.json', E + '/census-wide/table/pairs_2000_3000.json'] + [E + '/census-wide/pairs_%s.json' % r for r in ('3000_5000', '5000_7000', '7000_9000', '9000_11000')]:
    L |= set(tuple(x) for x in json.load(open(f))['open'])
L = {k for k in L if A < k[0] <= B}
print('%d < m <= %d: %d ratios pass the congruences (this code); the committed lists hold %d; the same set: %s' % (A, B, len(out), len(L), set(out) == L))
assert set(out) == L and len(out) == 20682
