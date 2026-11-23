#!/usr/bin/env python3
"""The local sets W_p of gauss_descent.py against actual integer points.

W_p is meant to contain the class tuple of every primitive point of Z_p^2.  Every coprime integer pair (u, v) is such a
point, and for an integer pair the five classes are computed exactly (no precision involved).  So for every coprime (u, v)
with all five values non-zero, some element of W_p must match its classes at every place.  A missing pattern is the unsafe
direction (a tuple could be removed wrongly), so this is the test that matters for the eliminations.

Points tried for each ratio and prime: all residues modulo p^e (p^(2e) about 4000) lifted by random higher digits, random
large pairs, and pairs at which one of the forms is divisible by a high power of a place above p (where the lifting and the
wildcards of the program are at work).

usage: local_points.py [seed] m,n [m,n ...]"""
import sys, random
from math import gcd
import gauss_descent as gd

def points(m, n, p, rng):
    e = 1
    while p**(2*(e + 1)) <= 4096: e += 1
    q = p**e; out = []
    for u in range(q):
        for v in range(q):
            if u % p == 0 and v % p == 0: continue
            out.append((u + q*rng.randint(0, 10**6), v + q*rng.randint(0, 10**6)))
    for _ in range(300): out.append((rng.randint(1, 10**9), rng.randint(1, 10**9)))
    # high valuation of a form: A u + B v = 0 modulo p^h in Z[i]/(p^h), solved by brute force over v = 1, u in a residue system mod p,
    # then Hensel by search (small p only), or directly at a split place
    F = gd.forms(m, n)
    if p % 4 == 1:
        L = gd.Local(p)
        for (A, B) in F:
            for r in (L.r, L.mod - L.r):
                a = (A[0] + A[1]*r) % L.mod; b = (B[0] + B[1]*r) % L.mod
                for h in (1, 2, 3, 4, 6):
                    ph = p**h
                    if a % p: out.append(((-b*pow(a, -1, ph)) % ph + ph*rng.randint(0, 50), 1 + ph*rng.randint(0, 50)))
                    if b % p: out.append((1 + ph*rng.randint(0, 50), (-a*pow(b, -1, ph)) % ph + ph*rng.randint(0, 50)))
    else:
        for (A, B) in F:                                 # inert or 2: points where the norm of the value is divisible by a high power of p
            best = []
            for u in range(1, 4*p*p + 60):
                for v in range(1, 60):
                    if gcd(u, v) != 1: continue
                    N = gd.gnorm(gd.value((A, B), u, v)); k = 0
                    while N and N % p == 0: N //= p; k += 1
                    if k >= 3: best.append((k, u, v))
            for k, u, v in sorted(best)[-40:]: out.append((u, v))
    return [(u, v) for u, v in out if gcd(u, v) == 1]

def run(m, n, rng, pmax=37):
    S = m*m + n*n; F = gd.forms(m, n); tested = 0; und = 0; missing = 0
    bad = sorted(set(gd.factor(2*m*n*(m*m - n*n)*S)))
    for p in sorted(set([q for q in bad if q <= 100] + [q for q in gd.PR if q <= pmax])):
        k = max(5, 8 if p == 2 else 6) if p <= 50 else 5                # the precision used by descent()
        L, W = gd.local_set(m, n, p, k)
        for (u, v) in points(m, n, p, rng):
            vals = [gd.value(f, u, v) for f in F]
            if any(z == (0, 0) for z in vals): continue
            c = [L.cls(z) for z in vals]
            if any(None in cj for cj in c): und += 1; continue
            tested += 1
            if not any(all(gd.match(c[j], w[j]) for j in range(5)) for w in W):
                missing += 1; print('MISSING', (m, n), p, (u, v), c, flush=True)
    return tested, und, missing

if __name__ == '__main__':
    args = sys.argv[1:]; seed = 1
    if args and ',' not in args[0]: seed = int(args[0]); args = args[1:]
    rng = random.Random(seed); T = U = M = 0
    for s in args:
        m, n = (int(x) for x in s.split(',')); t, u, mi = run(m, n, rng); T += t; U += u; M += mi
    print('%d ratios: %d integer points with exact classes tested against W_p (primes up to 37 and the bad primes up to 100); pattern missing from W_p: %d' % (len(args), T, M))
    assert M == 0
