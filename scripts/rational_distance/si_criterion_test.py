#!/usr/bin/env python3
"""A falsification test of a criterion from the literature (EXPERIMENT; RATIONAL_DISTANCE.md, section 4).

Y. Si, "Arithmetic constraints and covering curves for rational distances in a square", version 1.1.0 (14 November
2026), doi:10.5281/zenodo.24212887, Theorem "Independence criterion": on E_n: Y^2 = X^3 - n^2 X (n > 1 odd squarefree),
two rational points P, P' outside E_n[2] with
    delta(P') = delta(P)(n, 1, n),   u + u' != 0,   v_2(2(u - 1)(u + u')) even      (u = X/n)
are independent modulo torsion.  Test: on curves of rank exactly 1 (PARI ellrank, both bounds 1) every two points are
dependent, so the three hypotheses must never hold together.  Points tried: aM + T with M the point PARI returns,
a = 1..A, T in E_n[2].  No factoring: classes are compared by testing products for squares.  This does not prove the
criterion; a single pair with all three hypotheses would refute it.
usage: si_criterion_test.py [A]      (default 4; needs cypari2)"""
import cypari2, itertools, sys
from fractions import Fraction as F
from math import isqrt
p = cypari2.Pari(); p.allocatemem(2*10**9)
def fr(q): return F(int(p.numerator(q)), int(p.denominator(q)))
def issq(q):
    if q <= 0: return False
    a, b = isqrt(q.numerator), isqrt(q.denominator); return a*a == q.numerator and b*b == q.denominator
def v2(q):
    a, b, v = q.numerator, q.denominator, 0
    while a % 2 == 0: a //= 2; v += 1
    while b % 2 == 0: b //= 2; v -= 1
    return v
tested = viol = h1only = curves = 0
A = int(sys.argv[1]) if len(sys.argv) > 1 else 4
for n in range(3, 400, 2):
    if any(n % (q*q) == 0 for q in range(2, 20)): continue
    E = p.ellinit([0, 0, 0, -n*n, 0]); r = p.ellrank(E)
    if int(r[0]) != 1 or int(r[1]) != 1: continue
    if len(r[3]) < 1: continue
    M = r[3][0]
    if len(str(M[0])) > 60: continue
    curves += 1; pts = []
    for a in range(1, A + 1):
        aM = p.ellmul(E, M, a)
        for T in (None, p([0, 0]), p([n, 0]), p([-n, 0])):
            P = aM if T is None else p.elladd(E, aM, T)
            if len(P) == 2 and P[1] != 0: pts.append((a, fr(P[0])))
    for (a, X), (b, X2) in itertools.permutations(pts, 2):
        u, u2 = X/n, X2/n
        if u + u2 == 0: continue
        tested += 1
        h1 = issq(X*X2*n) and issq((X - n)*(X2 - n)) and issq((X + n)*(X2 + n)*n)
        h3 = v2(2*(u - 1)*(u + u2)) % 2 == 0
        h1only += h1
        if h1 and h3: viol += 1; print('VIOLATION n=%d a=%d b=%d' % (n, a, b))
print('curves E_n of rank 1 used (n odd squarefree below 400, x-coordinate of the point PARI returns at most 60 characters long): %d; ordered pairs of dependent points tested (multiples up to %d): %d; pairs with the Kummer hypothesis: %d; with all hypotheses (would contradict the criterion): %d' % (curves, A, tested, h1only, viol))
E = p.ellinit([0, 0, 0, -199*199, 0]); r = p.ellrank(E)
print('E_199 (the open line y = 26/225 of rd_scan_2_300.txt): PARI ellrank gives lower bound %d and upper bound %d for the rank' % (int(r[0]), int(r[1])))
