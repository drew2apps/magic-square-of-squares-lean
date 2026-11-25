#!/usr/bin/env python3
"""Checks written for the review of commit 0b34613 (COMPLETE.md, sections 5 and 7).  Python, and sympy for part 6; about ten seconds.

1. The first version of the Lean predicate `Passes` (ClassTable.lean; only the eight products) was weaker than "d is the class of a point of X_p": for p = 3 (mod 8)
   dividing m and d = p there is a point (u, v) (for p exactly dividing m, the point (0, 1)) at which d times each of the
   eight products is a square in Z_p[i], so that version held for every k, although the point is not in X_p and
   the class p is not allowed at such p (Theorem E).
2. Proposition A of section 7: if q = 5 or 13 divides m^2 + n^2 exactly once, there is no q-adic point.  A finite
   certificate: at every (u, v) modulo q^2, not both divisible by q, some form is a non-residue unit or has valuation 1.
3. Proposition B of section 7: at a prime q of the octic, the eight roots of the forms are four points, each a root of
   exactly two forms, and the number N of t in F_q at which the four forms are nonzero squares satisfies
   16 N >= 2q - 76 - 42 sqrt(q).  Checked for every (prime q < 3000, root of the octic); the bound is positive from
   q = 515, so for every prime q >= 521.
4. The rows E of class_good_table.txt (no point at a good prime) are all at primes = +-1 (mod 8), and at a good prime
   = +-3 (mod 8) one of the points t = 1, t = 0 is a point, for every ratio (p < 1200).
5. The rule at 2 used by class_selmer.py: a point modulo 64 exists exactly when 4 | mn or m = +-n (mod 8), for every
   (m, n) modulo 64 not both even and every (u, v) modulo 64 not both even.
6. The discriminant used for good reduction (COMPLETE.md, section 7): F0 F8 GL GU (t, 1) has degree 8 with leading
   coefficient 4 m^2 n^2 (m^2 + n^2)^2 and discriminant 2^48 (fW)^8, f = mn(m^4 - n^4), W the octic (sympy, symbolic).
   So at a prime not dividing 2fW the eight roots are distinct modulo p and none is at infinity.
"""
import os
from itertools import product
from math import isqrt
HERE = os.path.dirname(os.path.abspath(__file__))

def primes(N):
    s = bytearray([1]) * (N + 1); s[0] = s[1] = 0
    for i in range(2, isqrt(N) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(N + 1) if s[i]]
def pairs(m, n, u, v):
    return [(m * v, n * u), (n * v, m * u), ((m + n) * (u - v), (m - n) * (u + v)), ((m - n) * (u - v), (m + n) * (u + v))]
def forms(m, n, u, v): return [x * x + y * y for x, y in pairs(m, n, u, v)]
def gm(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def products(m, n, u, v):
    P0, P8, PL, PU = pairs(m, n, u, v); cj = lambda z, c: (z[0], -z[1]) if c else z
    return [gm(gm(gm(P0, cj(P8, a)), cj(PL, b)), cj(PU, c)) for a, b, c in product((0, 1), repeat=3)]
def vp(x, p):
    k = 0
    while x % p == 0: x //= p; k += 1
    return k
def leg(a, p): return pow(a % p, (p - 1) // 2, p) == 1

# 1. Passes is weaker than a local point
print('1. the eight products alone (the first version of Passes in ClassTable.lean) with d = p, for p = 3 (mod 8) dividing m, at a point (u, v):')
for p, m, n, u, v in [(3, 3, 1, 0, 1), (3, 3, 2, 0, 1), (3, 9, 1, 1, 3), (11, 11, 1, 0, 1), (19, 19, 2, 0, 1), (43, 43, 1, 0, 1)]:
    def gsq(z):                                   # z a square in Z_p[i], p = 3 (mod 4): even valuation, unit of square norm
        x, y = z; k = min(vp(x, p) if x else 99, vp(y, p) if y else 99)
        if k % 2: return False
        x //= p ** k; y //= p ** k
        return leg(x * x + y * y, p)
    ok = all(gsq((p * a, p * b)) for a, b in products(m, n, u, v))
    inX = all(N != 0 and vp(N, p) % 2 == 0 and leg(N // p ** vp(N, p), p) for N in forms(m, n, u, v))
    print('   p = %d, m = %d, n = %d, (u, v) = (%d, %d): p times each product a square in Z_p[i]: %s; the point in X_p: %s' % (p, m, n, u, v, ok, inX))

# 2. Proposition A
print('2. q = 5, 13 dividing m^2 + n^2 exactly once (n = 1, m modulo q^2), every (u, v) modulo q^2 not both divisible by q:')
for q in (5, 13):
    M = q * q; ts = [t for t in range(M) if (t * t + 1) % q == 0 and (t * t + 1) % M]
    bad = 0; pts = 0
    for t in ts:
        for u in range(M):
            for v in range(M):
                if u % q == 0 and v % q == 0: continue
                pts += 1
                fs = [f % M for f in forms(t, 1, u, v)]
                if not any((f % q and not leg(f, q)) or (f % q == 0 and f % M) for f in fs): bad += 1
    print('   q = %d: ratios %d, points %d; points where every form could be a q-adic square: %d' % (q, len(ts), pts, bad))
for q in (17, 29):
    M = q * q; ts = [t for t in range(M) if (t * t + 1) % q == 0 and (t * t + 1) % M]
    first = sum(any(all(f % q and leg(f, q) for f in forms(t, 1, u, 1)) for u in range(q)) for t in ts)
    print('   control, q = %d: ratios with a point of first order: %d of %d' % (q, first, len(ts)))

# 3. Proposition B
cases = viol = shape = 0; worst = None
for q in primes(3000)[1:]:
    for m in range(1, q):
        if (pow(m, 8, q) + 14 * pow(m, 4, q) + 1) % q or (2 * m * (pow(m, 4, q) - 1)) % q == 0: continue
        sq = [False] * q
        for a in range(1, q): sq[a * a % q] = True
        R = [set() for _ in range(4)]; N = 0; full_ok = True
        for t in range(q):
            fs = [f % q for f in forms(m, 1, t, 1)]
            for j in range(4):
                if fs[j] == 0: R[j].add(t)
            if all(fs):
                pr = fs[0] * fs[1] * fs[2] * fs[3] % q
                full_ok &= sq[pr]
                N += all(sq[f] for f in fs)
        pts = set().union(*R)
        shape += not (len(pts) == 4 and all(len(r) == 2 for r in R) and all(sum(t in r for r in R) == 2 for t in pts) and full_ok)
        cases += 1; viol += 16 * N < 2 * q - 76 - 42 * q ** .5
        if worst is None or N < worst[0]: worst = (N, q)
print('3. pairs (prime q < 3000, root of the octic): %d; not of the shape used (four points, two forms each, full product a square): %d;'
      ' with 16N below 2q - 76 - 42 sqrt(q): %d; fewest points: %d (at q = %d); the bound is positive from q = %d'
      % (cases, shape, viol, worst[0], worst[1], next(q for q in range(3, 10 ** 4) if 2 * q - 76 - 42 * q ** .5 > 0)))

# 4. rows E
E = []
for line in open(os.path.join(HERE, 'class_good_table.txt')):
    if line.startswith('#'): continue
    a = line.split(' ')
    if a[2]: E.append(int(a[0]))
miss = 0
for p in primes(1200):
    if p < 5 or p % 8 not in (3, 5): continue
    sq = [False] * p
    for a in range(1, p): sq[a * a % p] = True
    for m in range(2, p - 1):
        W = (pow(m, 8, p) + 14 * pow(m, 4, p) + 1) % p
        if (m * (m * m - 1) * (m * m + 1) * W) % p == 0: continue
        if not any(all(sq[f % p] for f in forms(m, 1, u, v)) for u, v in ((1, 1), (0, 1))): miss += 1
print('4. primes with a row E: %s; all = +-1 (mod 8): %s; good ratios at primes = +-3 (mod 8) below 1200 where neither t = 1 nor t = 0 is a point: %d'
      % (' '.join(map(str, E)), all(p % 8 in (1, 7) for p in E), miss))

# 5. the rule at 2
M = 64; sq = set(a * a % M for a in range(M)); bad = cnt = 0
for m, n in product(range(M), repeat=2):
    if m % 2 == 0 and n % 2 == 0: continue
    has = any(all(f % M in sq for f in forms(m, n, u, v)) for u, v in product(range(M), repeat=2) if u % 2 or v % 2)
    cnt += 1; bad += has != ((m * n) % 4 == 0 or (m - n) % 8 == 0 or (m + n) % 8 == 0)
print('5. (m, n) modulo 64 not both even: %d; a point modulo 64 differs from the rule "4 | mn or m = +-n (mod 8)": %d' % (cnt, bad))

# 6. the discriminant
from sympy import symbols, discriminant, factor, expand, Poly
ms, ns, t = symbols('m n t')
F = [(ms)**2 + (ns*t)**2, (ns)**2 + (ms*t)**2, ((ms + ns)*(t - 1))**2 + ((ms - ns)*(t + 1))**2, ((ms - ns)*(t - 1))**2 + ((ms + ns)*(t + 1))**2]
P = Poly(expand(F[0]*F[1]*F[2]*F[3]), t)
f = ms*ns*(ms**4 - ns**4); Wo = ms**8 + 14*ms**4*ns**4 + ns**8
print('6. F0 F8 GL GU (t, 1): degree %d; leading coefficient %s; discriminant / (2^48 (fW)^8) = %s'
      % (P.degree(), factor(P.LC()), factor(discriminant(P.as_expr(), t)/(2**48*(f*Wo)**8))))
