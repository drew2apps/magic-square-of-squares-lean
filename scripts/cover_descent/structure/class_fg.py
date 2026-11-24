#!/usr/bin/env python3
"""Brute-force tests of Theorems F and G of LAST_PAIR.md, section 17 (both are proved there and in Lean; this is a test).

For a modulus M and ALL residues (m, n, u, v) modulo M with the coprimality conditions and the four norms squares modulo M,
list the residues d such that d times each of the eight products P0 * P8^± * PL^± * PU^± is a Gaussian square modulo M.
  F: modulo 8 and 16, with (m, n) not both even and (u, v) not both even: every such d is even.
  G: modulo p = 5, 13, 17 and modulo 25, with p | m^2 + n^2 and p not dividing m: every such d that is not 0 modulo p
     is a square modulo p.  (Modulo 25 the pairs with 25 | m^2 + n^2 are included: the theorem has no condition on the
     power of p.)
Python only; under a minute.
"""
import sys, os, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from lean_cert import forms, Q, masks
FLAGS = list(itertools.product((0, 1), repeat=3))
def test(M, p, cond_mn, check):
    sq, g = masks(M); bad = tot = pts = 0; R = range(M)
    for m in R:
        for n in R:
            if not cond_mn(m, n): continue
            for u in R:
                for v in R:
                    if u % p == 0 and v % p == 0: continue
                    P = forms(m, n, u, v)
                    if not all((z[0]*z[0] + z[1]*z[1]) % M in sq for z in P): continue
                    pts += 1; Qs = [Q(P, f) for f in FLAGS]
                    for d in R:
                        if (p == 2 or d % p) and all(((d*X) % M, (d*Y) % M) in g for X, Y in Qs):
                            tot += 1; bad += not check(d)
    return pts, tot, bad
for M in (8, 16):
    print('F, modulo %d: points %d, classes that pass %d, odd among them %d' % ((M,) + test(M, 2, lambda m, n: (m % 2, n % 2) != (0, 0), lambda d: d % 2 == 0)))
for p, N in ((5, 1), (5, 2), (13, 1), (17, 1)):
    M = p**N; qr = set(x*x % p for x in range(1, p))
    r = test(M, p, lambda m, n: (m*m + n*n) % p == 0 and m % p != 0, lambda d: d % p in qr)
    print('G, modulo %d: points %d, classes (units at %d) that pass %d, non-residues among them %d' % (M, r[0], p, r[1], r[2]))
