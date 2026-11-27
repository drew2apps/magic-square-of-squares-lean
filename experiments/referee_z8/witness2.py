"""Review z8 of roadmap H (written by an AI model; not read by a person; experiment, not counted).
Explicit 2-adic points for the known and extra tuples of H1's 24 ratios (the repository's test at 2 is one-sided).
A tuple (d_j) is allowed at 2 if some coprime integers (u, v) and a 2-adic unit lam in {1, 3, 5, 7} make
lam d_j a_j(u, v) a square in Q_2(i) for all j (a_j: m v + i n u, n v + i m u, (m+n)(u-v) + i(m-n)(u+v),
(m-n)(u-v) + i(m+n)(u+v), H1's forms, same classes as tuple_descent's).  Square test: z = (1+i)^k w, k even, and w = x^2
mod 8 for some x in Z[i] (Hensel needs mod (1+i)^5).  Search |u|, v <= 64.  Negative control: three forged tuples of
(1109503,597671) get no point.  Rests on that test only.
usage: witness2.py [--write]   (witness2.txt; a few s)"""
import sys, os
from math import gcd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from z8common import emit, load_phantoms
def gm(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
SQ8 = set()
for x in range(8):
    for y in range(8):
        if (x + y) % 2 == 1: s = gm((x, y), (x, y)); SQ8.add((s[0] % 8, s[1] % 8))
def is_sq2(z):
    x, y = z
    if x == 0 and y == 0: return False
    k = 0
    while (x + y) % 2 == 0: x, y = (x + y) // 2, (y - x) // 2; k += 1      # z / (1 + i)
    return k % 2 == 0 and (x % 8, y % 8) in SQ8
def forms(m, n, u, v): return [(m*v, n*u), (n*v, m*u), ((m+n)*(u-v), (m-n)*(u+v)), ((m-n)*(u-v), (m+n)*(u+v))]
def witness(m, n, T, B=64):
    for u in range(-B, B + 1):
        for v in range(0, B + 1):
            if gcd(u, v) != 1 or (v == 0 and u != 1): continue
            F = forms(m, n, u, v)
            for lam in (1, 3, 5, 7):
                if all(is_sq2((lam * d * a, lam * d * b)) for d, (a, b) in zip(T, F)): return (u, v, lam)
    return None
out = ['z8: explicit 2-adic points (u, v, lam) for the tuples of H1\'s 24 ratios. Written by an AI model; not read by a person.']
nf = 0; nt = 0
for (m, n), (known, extra) in sorted(load_phantoms().items()):
    for kind, L in (('known', known), ('extra', extra)):
        for T in L:
            w = witness(m, n, T); nt += 1; nf += w is None
            out.append('(%d,%d) %s %s: %s' % (m, n, kind, T, w))
for T in [(1109503, 597671, 6, 1), (1109503, 597671, 2, 3), (3, 597671, 2, 1)]:
    out.append('negative control (1109503,597671) forged %s: %s' % (T, witness(1109503, 597671, T)))
out.append('tuples: %d; without a 2-adic point: %d' % (nt, nf))
emit(out, 'witness2.txt')
