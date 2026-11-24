#!/usr/bin/env python3
"""Where the class d = 2 is removed (EXPERIMENT; LAST_PAIR.md, section 17).  For the 119 pairs of G1.exceptions1000 with
m^2 + n^2 neither a square nor twice a square (read from ../../../MagicSquareOfSquares/G1Census1000.lean): is d = 2
allowed by cover_descent.py at every prime of 2mn(m^2 - n^2) and at every prime dividing m^2 + n^2 to an odd power
(primes below 400), and which is the first prime of good reduction that rejects it?"""
import sys, os, re
from math import isqrt
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(here, '..'))
import cover_descent as cd
from sympy import factorint, primerange
from collections import Counter
src = open(os.path.join(here, '..', '..', '..', 'MagicSquareOfSquares', 'G1Census1000.lean'), encoding='utf-8').read()
allp = eval(re.search(r'def exceptions1000 : List \(ℕ × ℕ\) :=\s*(\[.*?\])', src, re.S).group(1))
sq = lambda x: isqrt(x)**2 == x
ratios = [(m, n) for m, n in allp if not sq(m*m + n*n) and not sq(2*(m*m + n*n))]
okbad = 0; killer = Counter(); nokill = []
for m, n in ratios:
    S = m*m + n*n
    P = sorted(set([2]) | set(q for x in (m, n, m - n, m + n) for q in factorint(abs(x))))
    bad = [p for p in sorted(set(P) | set(p for p, e in factorint(S).items() if e % 2)) if p < 400]
    def has2(p):
        L = cd.Local(m, n, p, 14 if p == 2 else (8 if p < 50 else 5)); al = L.allowed()
        return al == 'all' or L.dclass(2) in al
    if all(has2(p) for p in bad): okbad += 1
    else: print('d = 2 rejected at a bad prime or a prime of S:', m, n, [p for p in bad if not has2(p)])
    k = next((p for p in primerange(3, 1500) if p not in bad and S % p and not has2(p)), None)
    if k is None: nokill.append((m, n))
    else: killer[(k, k % 4)] += 1
print('ratios: %d; d = 2 is allowed at every prime of 2mn(m^2-n^2) and of S below 400: %d' % (len(ratios), okbad))
print('first good prime that rejects d = 2 (prime, prime mod 4): count', sorted(killer.items()))
print('not rejected by a good prime below 1500:', nokill)
