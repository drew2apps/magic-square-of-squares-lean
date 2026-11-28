"""Run 10 (conductor check of the research note of 27 Nov 2026, section 1): at each of the nine open ratios m:n,
the parameter t = +-(u/v)^(+-1) with u:v one of the nine open ratios. Count how many of the five forms
F0, F4, F8, L, U of G1Full.FullSquares are squares of rationals. Exact arithmetic."""
from fractions import Fraction as Fr
from math import isqrt
from collections import Counter
OPEN = [(195,28),(223,167),(440,279),(609,200),(627,364),(719,161),(809,409),(861,620),(991,263)]
def sq(q):
    if q < 0: return False
    a, b = q.numerator, q.denominator
    return isqrt(a)**2 == a and isqrt(b)**2 == b
def forms(m, n, t):
    F4 = Fr(m*m + n*n) * (t*t + 1) / 2
    return [m*m + n*n*t*t, F4, n*n + m*m*t*t, F4 - 2*m*n*t, F4 + 2*m*n*t]
hist = Counter(); best = []
for (m, n) in OPEN:
    for (u, v) in OPEN:
        for t in (Fr(u, v), Fr(-u, v), Fr(v, u), Fr(-v, u)):
            k = sum(sq(f) for f in forms(m, n, t))
            hist[k] += 1
            if k >= 1: best.append((m, n, t, [sq(f) for f in forms(m, n, t)]))
print('pairs tested:', sum(hist.values()))
print('number of square forms -> count:', dict(sorted(hist.items())))
for b in best: print('  ', b)
