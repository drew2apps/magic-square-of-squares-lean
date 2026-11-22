#!/usr/bin/env python3
"""Counts for Result 1b of G1_FULL_SYSTEM.md. Pure Python, no dependencies.

Height of a ratio m : n (coprime, 1 <= n < m) is m. Two conditions:
  A: m : n lies in one of the 484 obstructed classes of modp_classes.txt;
  B: some prime p = 5 (mod 8) divides m^2 + n^2 to an odd power.
Prints the exact counts for m <= 1000, the product of (1 - c_p/(p+1)), and, for a sample of
20000 pairs drawn uniformly from the coprime pairs 1 <= n < m <= H with the fixed seed below
(m and n uniform in [1, H], kept if n < m and coprime), the proportion that does not satisfy B."""
import collections, math, os, random
here = os.path.dirname(os.path.abspath(__file__))
cls = collections.defaultdict(set)
for line in open(os.path.join(here, 'modp_classes.txt')):
    p, r = line.split()[:2]; cls[int(p)].add(int(r))

def cond_b(N):
    x = N; d = 2
    while d * d <= x:
        if x % d == 0:
            k = 0
            while x % d == 0: x //= d; k += 1
            if d % 8 == 5 and k % 2: return True
        d += 1 if d == 2 else 2
    return x > 1 and x % 8 == 5

def cond_a(m, n):
    return any(n % p and (m * pow(n, -1, p)) % p in rs for p, rs in cls.items())

print('c_p:', ', '.join('%d:%d' % (p, len(cls[p])) for p in sorted(cls)), ' total', sum(len(v) for v in cls.values()))
prod = 1.0
for p in cls: prod *= 1 - len(cls[p]) / (p + 1)
print('product over the 34 primes of (1 - c_p/(p+1)) = %.7f' % prod)
tot = a = b = ab = 0
for m in range(2, 1001):
    for n in range(1, m):
        if math.gcd(m, n) != 1: continue
        tot += 1; x = cond_a(m, n); y = cond_b(m * m + n * n)
        a += x; b += y; ab += (x and y)
print('m <= 1000: pairs %d, A %d, B %d, both %d, A or B %d, neither %d' % (tot, a, b, ab, a + b - ab, tot - a - b + ab))
random.seed(1)
for H in (10**3, 10**4, 10**5, 10**6):
    k = s = 0
    while k < 20000:
        m = random.randint(1, H); n = random.randint(1, H)
        if n >= m or math.gcd(m, n) != 1: continue
        k += 1; s += not cond_b(m * m + n * n)
    print('H = %d: %.4f of the sample does not satisfy B' % (H, s / k))
