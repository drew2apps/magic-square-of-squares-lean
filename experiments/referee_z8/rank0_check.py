"""Review z8 of roadmap H (written by an AI model; not read by a person; experiment, not counted).
A second implementation for H1 part 6 (claim 6 of the review), with PARI as the tool: for the 7 ratios where H1 finds
y^2 = w f0 f8 of rank 0, the curve A: Y^2 = x(x + m^4)(x + n^4) (x = m^2 n^2 t^2, Y = m^2 n^2 t y0 y8), PARI's ellrank
bounds on every curve of its isogeny class (ellisomat), PARI's torsion, the gcd of #E(F_p) over primes 3 <= p < 400 of
good reduction (p not dividing 2mn(m^4 - n^4)), and, computed exactly from the torsion points, the t with t^2 = x/(m^2 n^2).
PARI's rank bounds are not certified here; this is a cross-check of h1_rank0_exact.py, not a proof.
usage: rank0_check.py [--write]   (rank0_check.txt; about 4.5 minutes)"""
import sys, os, itertools
from fractions import Fraction as Fr
from math import isqrt, gcd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from z8common import emit, pari
P = pari()
def issq(q): q = Fr(q); return q >= 0 and isqrt(q.numerator) ** 2 == q.numerator and isqrt(q.denominator) ** 2 == q.denominator
out = ['z8: PARI cross-check of H1\'s rank-0 quotients y^2 = w f0 f8. Written by an AI model; not read by a person.']
for m, n in [(561607, 264401), (1109503, 597671), (1992873, 999664), (175775580, 15972419), (345152671, 33026129), (304012873, 159269639), (349960367, 118678903)]:
    a, b = m ** 4, n ** 4
    E = P.ellinit([0, a + b, 0, a * b, 0]); tors = P.elltors(E)
    pts = []
    for ks in itertools.product(*[range(int(o)) for o in tors[1]]):
        Q = [0]
        for k, g in zip(ks, tors[2]): Q = P.elladd(E, Q, P.ellmul(E, g, k))
        pts.append(Q)
    ts = set()
    for Q in pts:
        if len(Q) == 1: ts.add('inf'); continue
        w = Fr(str(Q[0])) / (m * m * n * n)
        if issq(w): t = Fr(isqrt(w.numerator), isqrt(w.denominator)); ts |= {str(t), str(-t)}
    g = 0
    for p in range(3, 400):
        if all(p % q for q in range(2, isqrt(p) + 1)) and (m * n * (m ** 4 - n ** 4)) % p: g = gcd(g, int(P.ellcard(E, p)))
    rb = [(int(r[0]), int(r[1])) for r in (P.ellrank(P.ellinit([c[0][0], c[0][1]]), 2) for c in P.ellisomat(E)[0])]
    out.append('(%d,%d): torsion %s (%d points), gcd #E(F_p) %d; t over torsion: %s; ellrank over the isogeny class: %s' % (m, n, [int(x) for x in tors[1]], int(tors[0]), g, sorted(ts), rb))
emit(out, 'rank0_check.txt')
