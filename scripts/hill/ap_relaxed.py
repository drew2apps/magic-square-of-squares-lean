#!/usr/bin/env python3
"""A problem weaker than the magic square of squares (HILL_NOTE.md, section 5).

Three arithmetic progressions of three squares x^2 < y^2 < z^2, with differences D that need not be equal, such that the
three numbers w = y^2 / D are distinct and in arithmetic progression.  A progression is, up to a square factor, given by
coprime p > q >= 1 of different parity:  y = p^2 + q^2,  D = 4pq(p^2 - q^2),  so  w = (p^2 + q^2)^2 / (4pq(p^2 - q^2)).
A magic square of squares gives such a triple (its three progressions have the same D and their middle terms are in
arithmetic progression), so if there is none with p <= P, no magic square of squares has its three progressions of shape
p <= P, whatever their sizes.

The search: all values w with p <= P, sorted, as floats; for every pair (w1, w3) the mean is looked up among the neighbours of
its position, and a candidate is confirmed with exact fractions.  Self-test: three made-up fractions in arithmetic progression
are added to the list and must be found.
usage: ap_relaxed.py P [first row, last row as fractions of the list, for running in parts]"""
import math, sys, numpy as np
from fractions import Fraction as Fr
P = int(sys.argv[1]); lo, hi = (float(sys.argv[2]), float(sys.argv[3])) if len(sys.argv) > 3 else (0.0, 1.0)
vals = {}
for p in range(2, P + 1):
    for q in range(1, p):
        if math.gcd(p, q) == 1 and (p - q) % 2 == 1: vals.setdefault(Fr((p*p + q*q)**2, 4*p*q*(p*p - q*q)), (p, q))
planted = [Fr(1000003, 7919), Fr(1000003, 7919) + Fr(104729, 15485863), Fr(1000003, 7919) + 2*Fr(104729, 15485863)]
for w in planted: assert w not in vals; vals[w] = None
ws = sorted(vals); f = np.array([float(w) for w in ws]); n = len(f); out = set()
for i in range(int(lo*n), min(n, int(hi*n) + 1)):
    mid = (f[i] + f[i+1:])/2; j = np.searchsorted(f, mid)
    for d in (-2, -1, 0, 1):
        jj = np.clip(j + d, 0, n - 1)
        for k in np.nonzero(np.abs(f[jj] - mid) < 1e-11*mid)[0]:
            a, b, c = ws[i], ws[int(jj[k])], ws[i + 1 + int(k)]
            if a + c == 2*b and a != b: out.add((a, b, c))
real = sorted(t for t in out if not any(vals[x] is None for x in t))
print('shapes with q < p <= %d: %d values of w%s' % (P, n - 3, '' if (lo, hi) == (0.0, 1.0) else '; rows %g to %g of the list' % (lo, hi)))
if int(lo*n) <= ws.index(planted[0]) <= int(hi*n): print('self-test, the planted progression is found: %s' % (tuple(planted) in out)); assert tuple(planted) in out
else: print('self-test: the planted progression starts in another part of the list')
print('three distinct values in arithmetic progression: %d' % len(real))
for t in real[:20]: print('  shapes', [vals[x] for x in t], 'w =', ', '.join(str(x) for x in t))
