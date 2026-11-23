#!/usr/bin/env python3
"""A 2-descent on one double cover, for one ratio (EXPERIMENT; LAST_PAIR.md, section 11).  Python and sympy only.

Y: the curve  F0, F8, GL, GU all squares  in t = u/v; it contains every magic square of squares with this ratio.
Over Z[i], with (u, v) coprime integers:
    A = -S u^2 + 2 D i u v + S v^2,   B = (m v + i n u)(n v - i m u),   B' = (m v - i n u)(n v + i m u),
    N(A) = GL GU / 4,   B B' = F0 F8      (S = m^2 + n^2, D = m^2 - n^2).
At a rational point of Y the number A B has square norm, so A B = d z^2 with d a squarefree positive integer (-1 is a
square in Z[i]); d is the class of the point, and the point lies on the double cover W_d of LAST_PAIR.md:
d A B = square, d A B' = square in Q(i).

Support.  Each of the four forms is the norm of a linear form over Z[i]:  F0 = N(m v + i n u),  F8 = N(n v - i m u),
GL = N((1 + i)(m u - n v + i (m v - n u))),  GU = N((1 + i)(m u + n v + i (m v + n u))),  and A B is, up to a unit and a
power of 1 + i, the product of two of these linear forms and the conjugates of the two others.  If F is a square and
gamma is its linear form, an odd prime p can divide the squarefree part of the class of gamma only if a prime above p
divides gamma and its conjugate, and then p | 2 m n (for F0, F8) or p | 2 (m^2 - n^2) (for GL, GU); these divisibilities
are in Lean (GaussDescent.lean).  So d is a product of primes dividing 2 m n (m^2 - n^2): a finite set P.

Local test.  For a prime p, d can be the class of a rational point only if there is a (u : v) in P^1(Q_p) with
    F0, F8, GL, GU squares in Q_p,   d A B and d A B' squares in Q_p(i) (in both factors, if p = 1 mod 4).
The program covers P^1(Z_p) by discs t0 + p^k Z_p in the two charts (t, 1) and (1, t), p | t.  On a disc a value g(t) is
congruent to g(t0) modulo p^k, so its square class is constant as soon as the valuation of g(t0) is small enough against k;
a disc where a constant class is not the right one is dropped, a disc where all six are constant and right is a point, and
any other disc is cut into p discs, down to a depth K.  At depth K what is not constant is counted as possible.  So the set
of classes allowed at p can only be too large, never too small: a d that is rejected is rejected by a proof.
Output: the primes of P, and the d that pass at every prime tested, with the classes of the known points marked
(t = 1, -1 if S is a square; t = 0, infinity if 2 S is a square).
usage: cover_descent.py m,n [largest good prime to test, default 150]
       cover_descent.py m,n selftest     random points against the allowed sets"""
import sys
from math import gcd
from sympy import factorint, primerange

def vp(x, p):
    if x == 0: return None
    v = 0
    while x % p == 0: x //= p; v += 1
    return v
def leg(a, p): return 0 if pow(a % p, (p - 1)//2, p) == 1 else 1
def gmul(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
SQ8 = sorted(set(((x*x - y*y) % 8, (2*x*y) % 8) for x in range(8) for y in range(8) if (x + y) % 2))
def key2(z):
    """The class of a nonzero Gaussian integer in Q_2(i)* modulo squares, and its valuation at 1 + i."""
    x, y = z; w = 0
    while (x + y) % 2 == 0: x, y = (x + y)//2, (y - x)//2; w += 1
    return w, (w % 2, min(((x*a - y*b) % 8, (x*b + y*a) % 8) for a, b in SQ8))

class Local:
    def __init__(self, m, n, p, K):
        self.m, self.n, self.p, self.K = m, n, p, K; self.S, self.D = m*m + n*n, m*m - n*n
        self.kind = 2 if p == 2 else (1 if p % 4 == 1 else 3)
        if self.kind == 1:
            r = next(x for x in range(2, p) if (x*x + 1) % p == 0); M = p**(K + 6)
            for _ in range(K + 8): r = (r - (r*r + 1)*pow(2*r, -1, M)) % M
            self.iota, self.M = r, M
    def vals(self, u, v):
        m, n, S, D = self.m, self.n, self.S, self.D
        A = (-S*u*u + S*v*v, 2*D*u*v); B = (m*n*(u*u + v*v), -D*u*v); Bc = (B[0], -B[1])
        w = S*(u*u + v*v)
        return (m*m*v*v + n*n*u*u, n*n*v*v + m*m*u*u, 2*(w - 4*m*n*u*v), 2*(w + 4*m*n*u*v)), gmul(A, B), gmul(A, Bc)
    def qclass(self, N, k):
        """(determined, class) of the integer N in Q_p* modulo squares, on a disc of radius p^k."""
        p = self.p; v = vp(N, p)
        if v is None: return False, None
        if p == 2: return (v + 3 <= k), (v % 2, (N >> v) % 8)
        return (v < k), (v % 2, leg(N//p**v, p))
    def kclass(self, z, k):
        """(determined, class) of the Gaussian integer z in (Q_p(i))* modulo squares, on a disc of radius p^k."""
        p = self.p
        if z == (0, 0): return False, None
        if self.kind == 2:
            w, key = key2(z); return (w + 6 <= 2*k), key
        if self.kind == 1:
            w = (z[0] + self.iota*z[1]) % self.M; v = vp(w, p)
            if v is None or v >= self.K + 3: return False, None
            return (v < k), (v % 2, leg(w//p**v, p))
        v = min(x for x in (vp(z[0], p), vp(z[1], p)) if x is not None) if z[0] and z[1] else vp(z[0] or z[1], p)
        a, b = (z[0]//p**v) % p, (z[1]//p**v) % p; e = (p*p - 1)//2; r = (1, 0); base = (a, b)   # (a + b i)^((p^2-1)/2) in F_p[i]
        while e:
            if e & 1: r = ((r[0]*base[0] - r[1]*base[1]) % p, (r[0]*base[1] + r[1]*base[0]) % p)
            base = ((base[0]*base[0] - base[1]*base[1]) % p, (2*base[0]*base[1]) % p); e >>= 1
        return (v < k), (v % 2, 0 if r == (1, 0) else 1)
    def dclass(self, d):
        """The class of the positive integer d in (Q_p(i))* modulo squares, with the same labels as kclass."""
        if self.kind == 2: return key2((d, 0))[1]
        v = vp(d, self.p)
        return (v % 2, leg(d//self.p**v, self.p) if self.kind == 1 else 0)
    def one(self, x):
        if self.p == 2: return x == (0, 1)
        return x == (0, 0)
    def allowed(self):
        """A set of classes (labels of dclass) that contains the class of every d possible at p, or 'all'."""
        p, K = self.p, self.K; out = set(); stack = [(c, t0, 1) for c in (0, 1) for t0 in range(p) if c == 0 or t0 == 0]
        while stack:
            c, t0, k = stack.pop(); u, v = (t0, 1) if c == 0 else (1, t0)
            Ns, g1, g3 = self.vals(u, v)
            qs = [self.qclass(N, k) for N in Ns]; e1, k1 = self.kclass(g1, k); e3, k3 = self.kclass(g3, k)
            if any(dd and not self.one(cc) for dd, cc in qs) or (e1 and e3 and k1 != k3): continue
            if all(dd for dd, cc in qs) and e1 and e3: out.add(k1); continue
            if k < K: stack += [(c, t0 + j*p**k, k + 1) for j in range(p)]; continue
            if e1: out.add(k1)                              # depth K: what is constant still binds, the rest counts as possible
            elif e3: out.add(k3)
            else: return 'all'
        return out

def run(m, n, pmax=150, K=None, verbose=True):
    S, D = m*m + n*n, m*m - n*n
    P = sorted(set([2]) | set(q for x in (m, n, m - n, m + n) for q in factorint(abs(x))))
    sqf = lambda x: int(__import__('math').prod(q for q, e in factorint(x).items() if e % 2))
    issq = lambda x: int(x**0.5 + 0.5)**2 == x
    known = {}
    if issq(S): known['t = 1, -1'] = sqf(2*abs(D))
    if issq(2*S): known['t = 0, infinity'] = sqf(S*m*n)
    tests = sorted(set(P) | set(primerange(3, pmax + 1)))
    cands = [1]
    for q in P: cands += [c*q for c in cands]
    info = []
    for p in tests:
        L = Local(m, n, p, K or (14 if p == 2 else (8 if p < 50 else 5))); al = L.allowed()
        if al == 'all': info.append((p, 'all')); continue
        before = len(cands); cands = [d for d in cands if L.dclass(d) in al]
        for name, d in known.items(): assert L.dclass(d) in al, ('a known point is rejected', p, name, d)
        if len(cands) < before: info.append((p, len(cands)))
    return P, known, sorted(cands), info

def selftest(m, n, primes, samples=20000, seed=1):
    """Random points against the allowed sets: for random coprime (u, v), whenever the four forms are squares in Q_p (decided
    from the exact integers), the class of A B, computed from the exact value, has to be in the set allowed at p, and equal to
    the class of A B'.  Returns (points with the four forms square, failures)."""
    import random
    rnd = random.Random(seed); hits = fails = 0
    for p in primes:
        L = Local(m, n, p, 14 if p == 2 else 8); al = L.allowed()
        for _ in range(samples):
            u, v = rnd.randrange(-10**6, 10**6), rnd.randrange(1, 10**6)
            if gcd(u, v) != 1 or u == 0: continue
            Ns, g1, g3 = L.vals(u, v)
            if any(N == 0 for N in Ns) or g1 == (0, 0) or g3 == (0, 0): continue
            if not all(L.one(L.qclass(N, 10**6)[1]) for N in Ns): continue
            hits += 1; c1 = L.kclass(g1, 10**6)[1]; c3 = L.kclass(g3, 10**6)[1]
            if c1 is None or c3 is None: continue
            if c1 != c3 or (al != 'all' and c1 not in al): fails += 1
    return hits, fails

if __name__ == '__main__':
    if len(sys.argv) > 2 and sys.argv[2] == 'selftest':
        m, n = (int(x) for x in sys.argv[1].split(',')); pr = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 53]
        hits, fails = selftest(m, n, pr)
        print('(%d,%d): random points at the primes %s: %d with the four forms square in Q_p; class outside the allowed set, or different for A B and A B\': %d' % (m, n, pr, hits, fails)); sys.exit(0)
    m, n = (int(x) for x in sys.argv[1].split(',')); pmax = int(sys.argv[2]) if len(sys.argv) > 2 else 150
    P, known, cands, info = run(m, n, pmax)
    print('(%d,%d): primes that can divide d: %s (%d classes)' % (m, n, P, 2**len(P)))
    print('classes of the known points: ' + ('; '.join('%s: d = %d' % kv for kv in known.items()) or 'none (m^2 + n^2 is neither a square nor twice a square)'))
    print('classes left after each prime that removes some: ' + ', '.join('%d: %s' % x for x in info))
    print('left: %s' % cands, ('(no class is possible: no rational point with the four forms square)' if not cands else '(exactly the classes of the known points)') if sorted(set(known.values())) == cands else '(more than the classes of the known points)')
