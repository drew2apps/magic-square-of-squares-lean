"""Reviewer z13: the job list of z13_indep.py (fixed seed; ratios outside the pilot).
Written by an AI model (reviewer z13), not read by a person.  Prints 'm n cv category' lines."""
import os, re, glob, random
from math import gcd, isqrt
HERE = os.path.dirname(os.path.abspath(__file__))
T10 = os.path.join(HERE, '..', 'conductor-t10')
pilot = set()
for f in glob.glob(os.path.join(T10, 't10_*.txt')):
    for l in open(f):
        mm = re.match(r'ratio \((\d+),(\d+)\)', l)
        if mm: pilot.add((int(mm[1]), int(mm[2])))
def v2(x):
    v = 0
    while x % 2 == 0: x //= 2; v += 1
    return v
def issq(x): r = isqrt(x); return r * r == x
def twoadic(m, n): return (m * n) % 4 == 0 or (m - n) % 8 == 0 or (m + n) % 8 == 0
def sqtype(m, n):
    S = m * m + n * n
    return issq(S) or (S % 2 == 0 and issq(S // 2))
rng = random.Random(20261127)
jobs = []; used = set(pilot)
def add(m, n, cat, curves=('08', 'LU')):
    used.add((m, n))
    for cv in curves: jobs.append((m, n, cv, cat))
def rand_ratio(H, pred):
    while True:
        m = rng.randint(2, H); n = rng.randint(1, m - 1)
        if gcd(m, n) != 1 or (m, n) in used: continue
        if pred(m, n): return m, n
# 1. random, m <= 20000, with a 2-adic point, not of square type (both parities come out of the draw)
for _ in range(40): add(*rand_ratio(20000, lambda m, n: twoadic(m, n) and not sqtype(m, n)), 'rand20k')
# 2. random, m <= 2000
for _ in range(30): add(*rand_ratio(2000, lambda m, n: twoadic(m, n) and not sqtype(m, n)), 'rand2k')
# 3. both odd (the triple divided by 2), random, m <= 20000
for _ in range(15): add(*rand_ratio(20000, lambda m, n: m % 2 == 1 and n % 2 == 1 and twoadic(m, n) and not sqtype(m, n)), 'bothodd')
# 4. A odd with v2(C -+ A) exactly 5: on E_08 (m, n of opposite parity, v2(mn) = 2); on E_LU (both odd, v2(m -+ n) = 3)
for _ in range(12): add(*rand_ratio(20000, lambda m, n: (m + n) % 2 == 1 and v2(m * n) == 2 and not sqtype(m, n)), 'v5_08')
for _ in range(12): add(*rand_ratio(20000, lambda m, n: m % 2 == 1 and n % 2 == 1 and (v2(m - n) == 3 or v2(m + n) == 3) and not sqtype(m, n)), 'v5_LU')
# 5. square type m^2 + n^2 = square (m, n the legs of a primitive triple) and twice-square type m^2 + n^2 = 2 square
sq = []; tw = []
for a in range(2, 120):
    for b in range(1, a):
        if gcd(a, b) != 1 or (a - b) % 2 == 0: continue
        p, q = a * a - b * b, 2 * a * b
        m, n = max(p, q), min(p, q)
        if m <= 20000 and (m, n) not in used: sq.append((m, n))
        m2, n2 = p + q, abs(p - q)
        if m2 <= 20000 and (m2, n2) not in used and gcd(m2, n2) == 1: tw.append((m2, n2))
for m, n in rng.sample(sorted(sq), 14): assert twoadic(m, n); add(m, n, 'square')
for m, n in rng.sample(sorted(tw), 14): assert twoadic(m, n); add(m, n, 'twicesquare')
# 6. all small ratios m <= 18 with a 2-adic point (not in the pilot)
for m in range(2, 19):
    for n in range(1, m):
        if gcd(m, n) == 1 and twoadic(m, n) and (m, n) not in used: add(m, n, 'small')
# 7. controls OUTSIDE the note's hypotheses (no 2-adic point; H2 fails on one curve)
for _ in range(8): add(*rand_ratio(20000, lambda m, n: (m + n) % 2 == 1 and v2(m * n) == 1), 'noH2_08', ('08',))
for _ in range(8): add(*rand_ratio(20000, lambda m, n: m % 2 == 1 and n % 2 == 1 and not twoadic(m, n)), 'noH2_LU', ('LU',))
for j in jobs: print(*j)
