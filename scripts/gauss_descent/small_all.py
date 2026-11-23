import sys, numpy as np
from math import gcd, isqrt
import gauss_descent as gd
M = int(sys.argv[1]); emp = 0; non = []; tot = 0
for m in range(2, M + 1):
    for n in range(1, m):
        if gcd(m, n) != 1: continue
        tot += 1; S = m*m + n*n
        triv = isqrt(S)**2 == S or (S % 2 == 0 and isqrt(S//2)**2 == S//2)
        D, T, hist, bad = gd.descent(m, n, 100)
        c = int(T.sum())
        if c == 0: emp += 1
        else: non.append((m, n, c, 'trivial point' if triv else 'NO TRIVIAL POINT'))
        if c == 0 and triv: print('BUG: empty with trivial point', m, n)
print('coprime pairs m <= %d: %d; empty descent set: %d; non-empty: %d' % (M, tot, emp, len(non)))
for x in non: print(x)
