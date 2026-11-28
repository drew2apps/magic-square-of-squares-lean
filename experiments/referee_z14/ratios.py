"""z14: the test ratios (fixed seed), none in the pilot of t10 and none whose mixed curve is a pilot curve
(partner (m+n, m-n), resp. ((m+n)/2, (m-n)/2)).  2-adic condition as rank0.py / class_selmer.py: 4 | mn or
m = +-n mod 8.  Kinds: nd (m^2+n^2 neither a square nor twice a square), sq (a square), tsq (twice a square).
Written by an AI model (reviewer z14); not read by a person.   usage: python3 -I ratios.py > ratios.txt"""
import os, re, glob, random
from math import gcd, isqrt
HERE = os.path.dirname(os.path.abspath(__file__))
pilot = set()
for f in glob.glob(os.path.join(HERE, '..', 'conductor-t10', 't10_*_[ab].txt')):
    for l in open(f):
        mm = re.match(r'ratio \((\d+),(\d+)\)', l)
        if mm: pilot.add((int(mm[1]), int(mm[2])))
def partner(m, n):
    a, b = m + n, m - n; g = gcd(a, b); return (a // g, b // g)
def curve_key(m, n):          # the opposite-parity member of {(m,n), partner}
    return (m, n) if (m - n) % 2 else partner(m, n)
pkeys = {curve_key(*r) for r in pilot}
def issq(x): r = isqrt(x); return r * r == x
def kind(m, n):
    S = m * m + n * n
    return 'sq' if issq(S) else ('tsq' if S % 2 == 0 and issq(S // 2) else 'nd')
def two_adic(m, n): return (m * n) % 4 == 0 or (m - n) % 8 == 0 or (m + n) % 8 == 0
rng = random.Random(14)
out = []; keys = set()
def take(m, n, why):
    k = curve_key(m, n)
    if (m, n) in pilot or k in pkeys or k in keys: return False
    keys.add(k); out.append((m, n, kind(m, n), 'odd-odd' if m % 2 and n % 2 else 'opposite', why)); return True
# 120 random non-degenerate ratios, m <= 15000: 80 of opposite parity, 40 odd-odd
c_op = c_oo = 0
while c_op < 80 or c_oo < 40:
    m = rng.randint(3, 15000); n = rng.randint(1, m - 1)
    if gcd(m, n) != 1 or not two_adic(m, n) or kind(m, n) != 'nd': continue
    oo = m % 2 == 1 and n % 2 == 1
    if oo and c_oo < 40 and take(m, n, 'random'): c_oo += 1
    elif not oo and c_op < 80 and take(m, n, 'random'): c_op += 1
# 8 square-type and 6 twice-square-type
sq = []; tsq = []
for a in range(2, 90):
    for b in range(1, a):
        if gcd(a, b) == 1 and (a - b) % 2:
            m, n = max(a * a - b * b, 2 * a * b), min(a * a - b * b, 2 * a * b)
            if m <= 15000 and two_adic(m, n): sq.append((m, n))
            p, q = partner(m, n)
            if p <= 15000 and two_adic(p, q): tsq.append((p, q))
rng.shuffle(sq); rng.shuffle(tsq)
i = 0
for r in sq:
    if i < 8 and take(*r, why='square type'): i += 1
i = 0
for r in tsq:
    if i < 6 and take(*r, why='twice-square type'): i += 1
print('# z14 test ratios (seed 14): m,n kind parity origin; %d ratios, %d distinct mixed curves; pilot ratios %d excluded with partners' % (len(out), len(keys), len(pilot)))
for m, n, k, par, why in out: print('%d,%d %s %s %s' % (m, n, k, par, why))
