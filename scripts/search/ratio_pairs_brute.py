"""Independent check of ratio_pairs.c for small H: every ordered pair of ratios is tested directly, with no grouping.
Usage: python3 ratio_pairs_brute.py H ../g1_full/modp_classes.txt   (prints the same line as ratio_pairs; H = 100 takes a few seconds)"""
import sys
from math import gcd, isqrt
H = int(sys.argv[1])
cls = [tuple(int(v) for v in l.split()[:2]) for l in open(sys.argv[2]) if l.strip() and not l.startswith('#')]
assert len(cls) == 484
def sq(x): return x >= 0 and isqrt(x) ** 2 == x
def oddpow_ok(S):
    x = S; p = 2
    while x % 2 == 0: x //= 2
    p = 3
    while p * p <= x:
        e = 0
        while x % p == 0: x //= p; e += 1
        if e % 2 and p % 8 != 1: return False
        p += 2
    return x == 1 or x % 8 == 1
def incls(m, n): return any(n % p and (m - r * n) % p == 0 for p, r in cls)
first = [(m, n) for m in range(1, H + 1) for n in range(1, m) if gcd(m, n) == 1]
second = [(m, n) for m in range(1, H + 1) for n in range(1, H + 1) if gcd(m, n) == 1]
F = [0] * 5; triv = 0
def trivial(S): return sq(S) or (S % 2 == 0 and sq(S // 2))
info2 = [(m, n, m * m + n * n, oddpow_ok(m * m + n * n), incls(m, n)) for m, n in second]
for m, n in first:
    S = m * m + n * n; o1 = oddpow_ok(S); c1 = incls(m, n)
    for m2, n2, S2, o2, c2 in info2:
        if not sq(2 * S * S2): continue
        F[0] += 1; triv += trivial(S) and trivial(S2)
        if not (o1 and o2): continue
        F[1] += 1
        if c1 or c2: continue
        F[2] += 1
        e2 = S * S2 // 2; k = 2 * m * n * m2 * n2
        ent = [n * n * m2 * m2 + m * m * n2 * n2, m * m * m2 * m2 + n * n * n2 * n2, e2 - k, e2 + k]
        if not all(sq(x) for x in ent): continue
        F[3] += 1
        allent = ent + [e2, (m * m2 + n * n2) ** 2, (m * m2 - n * n2) ** 2, (n * m2 + m * n2) ** 2, (n * m2 - m * n2) ** 2]
        if len(set(allent)) == 9 and 0 not in allent: F[4] += 1
print('H=%d ratios: first=%d second=%d pairs=%d | F1=%d (of these with both numbers a square or twice a square: %d) F2=%d F3=%d F4=%d F5=%d' % (H, len(first), len(second), len(first) * len(second), F[0], triv, F[1], F[2], F[3], F[4]))
