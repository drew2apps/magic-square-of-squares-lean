#!/usr/bin/env python3
"""The local eliminations of LAST_PAIR.md, section 11, by exhaustion (EXPERIMENT).  Shares no code with cover_descent.py.

For a ratio m : n the candidates are the products d of the primes of 2 m n (m^2 - n^2).  For a prime p and M = p^N this
program goes through EVERY point (u : v) of the projective line over Z/M and keeps a candidate d when, at some point,
    F0, F8, GL, GU are squares modulo M   (they are in a list of all x^2 modulo M), and
    d A B and d A B' are squares modulo M in Z[i]   (they are in a list of all (x + i y)^2 modulo M).
If d is the class of a point over Q_p, take (u, v) primitive in Z_p: the four forms are squares of p-adic integers and
d A B, d A B' are squares of elements of Z_p[i] (they are integral), so all six are squares modulo M and d is kept.
So a candidate that is NOT kept is not the class of any point over Q_p: this is a proof by exhaustion, with no valuation,
no Legendre symbol and no Hensel lemma in it.  It is weaker than cover_descent.py at a fixed N (a value divisible by M
counts as a square) and needs no depth.
Output: for each prime, N, the number of points with the four forms square, and the candidates kept; then the
candidates kept at every prime.
usage: cover_bruteforce.py m,n [largest prime, default 61] [largest modulus M, default 5000]"""
import sys

def factor(x):
    out = set(); x = abs(x); q = 2
    while q*q <= x:
        while x % q == 0: out.add(q); x //= q
        q += 1
    if x > 1: out.add(x)
    return out
def primes_to(b): return [q for q in range(2, b + 1) if all(q % r for r in range(2, int(q**0.5) + 1))]

def kept(m, n, cands, p, M):
    S, D = m*m + n*n, m*m - n*n
    sqZ = bytearray(M)
    for x in range(M): sqZ[x*x % M] = 1
    sqG = bytearray(M*M)
    for x in range(M):
        xx = x*x
        for y in range(M): sqG[((xx - y*y) % M)*M + (2*x*y) % M] = 1
    out = set(); points = 0; left = set(cands)
    for u, v in [(t, 1) for t in range(M)] + [(1, t) for t in range(0, M, p)]:
        w = S*(u*u + v*v)
        if not (sqZ[(m*m*v*v + n*n*u*u) % M] and sqZ[(n*n*v*v + m*m*u*u) % M] and sqZ[(2*(w - 4*m*n*u*v)) % M] and sqZ[(2*(w + 4*m*n*u*v)) % M]): continue
        points += 1
        a0, a1 = S*(v*v - u*u), 2*D*u*v                     # A = -S u^2 + 2 D i u v + S v^2
        b0, b1 = m*n*(u*u + v*v), -D*u*v                    # B = (m v + i n u)(n v - i m u);  B' = its conjugate
        g1 = ((a0*b0 - a1*b1) % M, (a0*b1 + a1*b0) % M); g3 = ((a0*b0 + a1*b1) % M, (-a0*b1 + a1*b0) % M)
        for d in list(left):
            if sqG[(d*g1[0] % M)*M + d*g1[1] % M] and sqG[(d*g3[0] % M)*M + d*g3[1] % M]: out.add(d); left.discard(d)
        if not left: break
    return points, out

if __name__ == '__main__':
    m, n = (int(x) for x in sys.argv[1].split(',')); pmax = int(sys.argv[2]) if len(sys.argv) > 2 else 61
    Mmax = int(sys.argv[3]) if len(sys.argv) > 3 else 5000
    P = sorted(factor(2*m*n*(m*m - n*n))); cands = [1]
    for q in P: cands += [c*q for c in cands]
    print('(%d,%d): candidates: the %d products of %s' % (m, n, len(cands), P))
    alive = set(cands)
    for p in primes_to(pmax):
        N = 1
        while p**(N + 1) <= Mmax: N += 1
        points, out = kept(m, n, cands, p, p**N)
        before = len(alive); alive &= out
        print('prime %d, modulo %d^%d: points with the four forms square: %d; candidates kept: %d of %d; kept at every prime so far: %d%s'
              % (p, p, N, points, len(out), len(cands), len(alive), '' if len(alive) == before else ' ' + str(sorted(alive)) if len(alive) <= 8 else ''))
    print('kept at every prime: %s' % sorted(alive))
