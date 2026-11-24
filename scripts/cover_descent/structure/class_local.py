#!/usr/bin/env python3
"""The class at the bad primes: checks for the proof of Theorem C (COMPLETE.md).  Python only; about a minute.

For a prime p and integers u, v at which the four norms are nonzero squares in Z_p, the class of the point is
D = prod (w_j + x_j) in Q_p^* modulo squares and modulo -1, where the four pairs are (x_j, y_j) and w_j^2 = x_j^2 + y_j^2
(2(w + x)(x + iy) = (w + x + iy)^2).  `cls` computes it exactly for odd p (the norms are integers; a square root to
30 digits).  Labels: '1' unit and residue, 'e' unit and non-residue, 'p', 'ep'; for p = 3 mod 4 only '1' and 'p'.

1. Upper bounds: at random points, for every odd prime below 100 dividing mn(m^2-n^2)(m^2+n^2)(octic), the class is in
   the set of Theorem C.
2. Lower bounds: the points of the proof (first order; u = p u1; the points for p = 5 when 25 divides) have the
   classes that the proof says, for p^e exactly dividing m, e = 1..4, p up to 97 and three larger primes.
3. The conic: for primes p = 1 mod 4, both values of chi(2v) occur among the v with 2(1 + v^2) a nonzero square,
   except for p = 5 (primes below 5000; for p >= 17 this is the bound of Hasse for y^2 = v^3 + v).
4. p = 2: modulo 64, in normal form, the residues of d that pass at a point where the norms are squares.
5. The table of images against the test modulo p^K of lean_cert.py, 180 cases.
"""
import sys, os, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import lean_cert as L
from math import gcd
random.seed(20261124)
def chi(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1)//2, p) == 1 else -1)
def vp(x, p):
    k = 0
    while x % p == 0: x //= p; k += 1
    return k
def sqrt_unit(a, p, K):                                     # a unit and a residue: a root modulo p^K
    r = next(x for x in range(1, p) if (x*x - a) % p == 0); M = p
    for _ in range(K):
        M *= p; r = (r - (r*r - a)*pow(2*r, -1, M)) % M
    return r
K = 30
def cls(m, n, u, v, p):
    xs = [m*v, n*v, (m + n)*(u - v), (m - n)*(u - v)]; ys = [n*u, m*u, (m - n)*(u + v), (m + n)*(u + v)]
    val = 0; ch = 1
    for x, y in zip(xs, ys):
        N = x*x + y*y
        if N == 0: return None
        k = vp(N, p)
        if k % 2 or chi(N//p**k, p) != 1: return None
        if y == 0: s = 2*x
        else:
            M = p**(K + k); w = p**(k//2)*sqrt_unit(N//p**k % p**K, p, K) % M
            s = min(((w + x) % M, (x - w) % M), key=lambda z: vp(z, p) if z else 10**9)
            assert s and vp(s, p) <= vp(y, p)
        val += vp(s, p); ch *= chi(s//p**vp(s, p), p)
    if p % 4 == 3: return '1' if val % 2 == 0 else 'p'
    return {(0, 1): '1', (0, -1): 'e', (1, 1): 'p', (1, -1): 'ep'}[(val % 2, ch)]
def predicted(p, e, kind):                                  # kind: 'supp' (p divides mn(m^2-n^2), p^e exactly) or 'other'
    if kind == 'other': return {'1'}
    if p % 8 == 3: return {'1'}
    if p % 8 == 7: return {'1', 'p'}
    if p % 8 == 5: return {'e'} if (p == 5 and e == 1) else {'1', 'e'}
    return {'1', 'e', 'p', 'ep'}
PR = [q for q in range(3, 100, 2) if all(q % r for r in range(3, int(q**0.5) + 1, 2))]

# 1. upper bounds at random points
pts = bad = 0; seen = {}
for _ in range(4000):
    m = random.randint(2, 400); n = random.randint(1, m - 1)
    if gcd(m, n) != 1: continue
    S = m*n*(m - n)*(m + n); O = (m*m + n*n)*(m**8 + 14*m**4*n**4 + n**8)
    for p in PR:
        if S % p and O % p: continue
        kind = 'supp' if S % p == 0 else 'other'; e = vp(S, p) if kind == 'supp' else 0
        for _ in range(40):
            u = random.randint(-p**3, p**3); v = random.randint(-p**3, p**3)
            if u % p == 0 and v % p == 0: continue
            c = cls(m, n, u, v, p)
            if c is None: continue
            pts += 1; seen.setdefault((p % 8, kind, min(e, 2) if p == 5 else 0), set()).add(c)
            if c not in predicted(p, e, kind): bad += 1
print('1. random points with a class: %d; outside the set of Theorem C: %d' % (pts, bad))
print('   classes seen, by (p mod 8, kind, exponent if p = 5): %s' % sorted((k, sorted(s)) for k, s in seen.items()))

# 2. the points of the proof, p^e exactly dividing m
n1 = n2 = n3 = n5 = wrong = 0
for p in PR + [193, 1009, 2017]:
    for e in (1, 2, 3, 4):
        for _ in range(3):
            while True:
                m1 = random.randint(1, 50); n = random.randint(1, 60); m = p**e*m1
                if m1 % p and n % p and gcd(m, n) == 1 and m != n: break
            for v in random.sample(range(1, p), min(p - 1, 12)):          # first order: class chi(2uv), u = 1
                if chi(2*(1 + v*v), p) != 1: continue
                n1 += 1; want = '1' if p % 4 == 3 else ('1' if chi(2*v, p) == 1 else 'e')
                wrong += cls(m, n, 1, v, p) != want
            if p % 8 in (1, 7) and e >= 2:                                   # u = p u1, v = 1: class p u1
                for u1 in random.sample(range(1, p), min(p - 1, 8)):
                    n2 += 1; want = 'p' if (p % 8 == 7 or chi(u1, p) == 1) else 'ep'
                    wrong += cls(m, n, p*u1, 1, p) != want
            if p % 8 in (1, 7) and e == 1:                                   # u = p u1, n u1 = (t - m1^2/t)/2: class p n t
                for t in range(1, p):
                    if (t*t - m1*m1) % p == 0 or (t*t + m1*m1) % p == 0: continue
                    u1 = (t - m1*m1*pow(t, -1, p))*pow(2*n, -1, p) % p
                    n3 += 1; want = 'p' if (p % 8 == 7 or chi(n*t, p) == 1) else 'ep'
                    wrong += cls(m, n, p*u1, 1, p) != want
            if p == 5 and e >= 2:                                            # v^2 = -1 modulo 25: a point of class 1
                got = set(cls(m, n, 1, v, 5) for v in range(5**6) if (1 + v*v) % 25 == 0) - {None}
                n5 += 1; wrong += not ('1' in got and got <= {'1', 'e'})
print('2. points of the proof: first order %d, u = p u1 with e >= 2: %d, with e = 1: %d, p = 5 with 25 dividing: %d ratios; not as the proof says: %d' % (n1, n2, n3, n5, wrong))

# 3. the conic
one = []; cnt = 0
for p in range(5, 5000, 4):
    if any(p % q == 0 for q in range(3, int(p**0.5) + 1, 2)): continue
    cnt += 1; got = set(chi(2*v, p) for v in range(1, p) if chi(2*(1 + v*v), p) == 1)
    if got != {1, -1}: one.append(p)
print('3. primes = 1 mod 4 below 5000: %d; with one value of chi(2v) only: %s' % (cnt, one))

# 4. p = 2
M = 64; sq, g = L.masks(M); norm = [(1, t) for t in range(M)] + [(2*s, 1) for s in range(M//2)]; res = set(); np2 = 0
for m, n in norm:
    for u, v in norm:
        P = L.forms(m, n, u, v)
        if not all((z[0]*z[0] + z[1]*z[1]) % M in sq for z in P): continue
        np2 += 1; Qs = [L.Q(P, f) for f in L.FLAGS]
        res |= set(d % 16 for d in range(M) if all(((d*X) % M, (d*Y) % M) in g for X, Y in Qs))
print('4. modulo 64, pairs in normal form: %d points with the four norms squares; residues of d modulo 16 that pass somewhere: %s' % (np2, sorted(res)))

# 5. the table against the test modulo p^K
def nonres(p): return next(a for a in range(2, p) if chi(a, p) == -1)
def image(m, n, p, Kp):
    Mp = p**Kp; surv = L.survivors(m, n, Mp, p); eps = nonres(p)
    reps = {'1': 1, 'p': p} if p % 4 == 3 else {'1': 1, 'e': eps, 'p': p, 'ep': eps*p}
    return set(k for k, d in reps.items() if not L.killed(d, surv, Mp))
rows = mis = 0
for p, es, K0 in [(3, (1, 2, 3), 3), (5, (1, 2), 2), (7, (1, 2), 1), (11, (1,), 1), (13, (1,), 1), (17, (1,), 1), (23, (1,), 1), (29, (1,), 1), (31, (1,), 1), (37, (1,), 1), (41, (1,), 0)]:
    for e in es:
        Kp = e + K0
        if p**Kp > 1400: Kp = e + 1 if p**(e + 1) <= 1700 else e
        for which in ('m', 'n', 'm-n', 'm+n'):
            for _ in range(3):
                while True:
                    a = p**e*random.choice([x for x in range(1, 40) if x % p]); b = random.randint(1, 60)
                    m, n = {'m': (a, b), 'n': (b, a), 'm-n': (a + b, b), 'm+n': (a - b, b)}[which]
                    if m > 0 and n > 0 and m != n and gcd(m, n) == 1 and (m*m + n*n) % p: break
                rows += 1; mis += image(m, n, p, Kp) != predicted(p, e, 'supp')
print('5. cases %d; image modulo p^K different from the table of Theorem C: %d' % (rows, mis))
