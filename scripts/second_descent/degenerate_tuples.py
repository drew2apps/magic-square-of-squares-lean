#!/usr/bin/env python3
"""The two degenerate solutions of a ratio have different classes in K*/K*^2 for at least one of the five forms.
So when the first descent leaves two tuples and keeps the tuples of both degenerate solutions (gauss_descent_157.txt),
the two tuples left are exactly those.  usage: degenerate_tuples.py pairs123.txt"""
import sys
from math import isqrt
from sage.all__sagemath_schemes import *
K = QuadraticField(-1, 'i'); I = K.gen()
pairs = [tuple(int(v) for v in s.split(',')) for s in open(sys.argv[1]).read().split()]
def forms(m, n): return [(I*n, K(m)), (I*m, K(n)), (K(1), I), (m + I*n, -n - I*m), (m + I*n, n + I*m)]      # a0, a8, g4, gL, gU as A u + B v
same = []; pat = {}
for m, n in pairs:
    S = m*m + n*n
    sols = [(1, 1), (1, -1)] if isqrt(S)**2 == S else [(0, 1), (1, 0)]
    vals = [[A*u + B*v for A, B in forms(m, n)] for u, v in sols]
    eq = tuple(bool((a/b).is_square()) for a, b in zip(*vals))
    pat[eq] = pat.get(eq, 0) + 1
    if all(eq): same.append((m, n))
print('%d ratios; pattern of equal classes (a0, a8, g4, gL, gU) at the two degenerate solutions:' % len(pairs), pat)
print('ratios where the two degenerate solutions have the same tuple:', same)
