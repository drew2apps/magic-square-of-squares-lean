#!/usr/bin/env python3
"""The fibres of the G1 system with two colliding branch points.

Res(F0, GL) = Res(F0, GU) = Res(F8, GL) = Res(F8, GU) = 4 (m^8 + 14 m^4 n^4 + n^8) up to squares of m, n.
For a prime p dividing O(m, n) = m^8 + 14 m^4 n^4 + n^8 the curve over F_p is singular and the Weil bound
for a smooth curve of genus 49 does not apply. There F0 F8 GL GU = (2 S m n)^2 (square of a quartic), the
other 30 products of the five forms are not constant multiples of squares, and character sums give
    32 N >= 2 p - 100 - 270 sqrt(p),   N = number of t in F_p at which all five values are nonzero squares,
which is positive for p >= 18400. This script checks the remaining range directly:
  1. the identity F0 F8 GL GU = (2 S m n)^2 (t^4 + c t^2 + 1)^2 modulo O(m, n), symbolically (n = 1);
  2. for every prime 9700 < p < 20000 and every root r of r^8 + 14 r^4 + 1 mod p, a value of t exists at
     which the five forms of (m, n) = (r, 1) are nonzero squares mod p;
  3. none of the 484 obstructed classes in modp_classes.txt is a root of the octic (so the statement
     'all obstructed classes are at ratios of good reduction' stays true with the octic included).
No dependencies."""
import os, sys

def primes(a, b):
    s = bytearray([1]) * (b + 1); s[0] = s[1] = 0
    for i in range(2, int(b**.5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(a, b + 1) if s[i]]

def forms(m, n):
    S = m*m + n*n
    return [(n*n, 0, m*m), (m*m, 0, n*n), (2*S, 0, 2*S), (2*S, -8*m*n, 2*S), (2*S, 8*m*n, 2*S)]

def has_point(r, p):
    fs = forms(r, 1); e = (p - 1)//2
    for t in range(p):
        if all(pow((a*t*t + b*t + c) % p, e, p) == 1 for (a, b, c) in fs): return True
    return False

# 1. polynomial identity modulo the octic, with exact integer polynomial arithmetic in m (n = 1), coefficients of t
def pmul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): r[i + j] += x * y
    return r
def padd(a, b):
    n = max(len(a), len(b)); a = a + [0]*(n - len(a)); b = b + [0]*(n - len(b))
    return [x + y for x, y in zip(a, b)]
def pmod_octic(a):                      # reduce a polynomial in m modulo m^8 + 14 m^4 + 1
    a = a[:]
    for k in range(len(a) - 1, 7, -1):
        c = a[k]; a[k] = 0; a[k - 4] -= 14 * c; a[k - 8] -= c
    while len(a) > 1 and a[-1] == 0: a.pop()
    return a
M = [0, 1]; M2 = pmul(M, M); S = padd(M2, [1])
# forms as polynomials in t with coefficients polynomials in m: lists [c0, c1, c2]
F0 = [M2, [0], [1]]; F8 = [[1], [0], M2]
GL = [[2*x for x in S], [-8*x for x in M], [2*x for x in S]]; GU = [GL[0], [8*x for x in M], GL[2]]
def tmul(A, B):
    R = [[0] for _ in range(len(A) + len(B) - 1)]
    for i, x in enumerate(A):
        for j, y in enumerate(B): R[i + j] = padd(R[i + j], pmul(x, y))
    return R
P = tmul(tmul(F0, F8), tmul(GL, GU))
# candidate: (2 S m)^2 (t^4 + c t^2 + 1)^2 with c determined from the t^2 and t^6 coefficients: c = P2 / (2 (2Sm)^2)
lead = pmul(pmul([2*x for x in S], M), pmul([2*x for x in S], M))          # (2 S m)^2
assert pmod_octic(P[0]) == pmod_octic(lead) and pmod_octic(P[8]) == pmod_octic(lead)
assert all(pmod_octic(P[k]) == [0] for k in (1, 3, 5, 7))
# square test: P = lead (t^4 + c t^2 + 1)^2 iff lead P4 = lead^2 (c^2 + 2) and P2 = P6 = 2 lead c; eliminate c:
# 4 lead P4 = P2^2 + 8 lead^2  (mod the octic)
lhs = pmod_octic(pmul([4*x for x in lead], P[4])); rhs = pmod_octic(padd(pmul(P[2], P[2]), [8*x for x in pmul(lead, lead)]))
assert pmod_octic(P[2]) == pmod_octic(P[6]) and lhs == rhs, 'identity fails'
print('1. F0 F8 GL GU is (2Smn)^2 times the square of t^4 + c t^2 + 1 modulo the octic: OK')

# 2. direct check
count = 0; classes = 0
for p in primes(9700, 20000):
    roots = [r for r in range(1, p) if (pow(r, 8, p) + 14 * pow(r, 4, p) + 1) % p == 0]
    if roots: count += 1
    for r in roots:
        classes += 1
        assert has_point(r, p), (p, r)
print('2. primes in (9700, 20000) with roots of the octic: %d; classes checked: %d; every one has a point: OK' % (count, classes))

# 3. obstructed classes versus the octic
here = os.path.dirname(os.path.abspath(__file__))
hits = 0; total = 0
for line in open(os.path.join(here, 'modp_classes.txt')):
    w = line.split()
    if len(w) < 2 or not w[0].isdigit() or not w[1].isdigit(): continue
    p, r = int(w[0]), int(w[1]); total += 1
    if (pow(r, 8, p) + 14 * pow(r, 4, p) + 1) % p == 0: hits += 1
print('3. obstructed classes read: %d; on the octic: %d' % (total, hits))
assert hits == 0 and total == 484
# 4. the exact number N of t for every class of step 2, against the bound:
#    |32 N - 2 p| <= 100 + 270 sqrt(p). The lower side is the bound of G1_FULL_SYSTEM.md; the upper side follows from
#    the same character sums. This tests the main term p/16 (not p/32) and the size of the constants on real data.
if True:
    worst = 0.0; nclass = 0
    for p in primes(9700, 20000):
        roots = [r for r in range(1, p) if (pow(r, 8, p) + 14 * pow(r, 4, p) + 1) % p == 0]
        if not roots: continue
        qr = bytearray(p)
        for x in range(1, p): qr[x * x % p] = 1
        for r in roots:
            fs = forms(r, 1); N = 0
            for t in range(p):
                if all(qr[(a * t * t + b * t + c) % p] for (a, b, c) in fs): N += 1
            dev = abs(32 * N - 2 * p); bound = 100 + 270 * p ** 0.5; nclass += 1
            assert dev <= bound, (p, r, N)
            worst = max(worst, dev / bound)
    print('4. exact counts for %d classes: |32N - 2p| <= 100 + 270 sqrt(p) in every class; largest ratio %.3f' % (nclass, worst))
print('octic_check OK')
