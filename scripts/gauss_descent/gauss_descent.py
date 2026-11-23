#!/usr/bin/env python3
"""Descent over the Gaussian integers for the five-form system of a ratio m : n.  EXPERIMENT.

A solution is a coprime integer pair (u, v), t0 = u/v, with
    N(a0) = square,  N(a8) = square,  2S N(g4) = square,  2 N(gL) = square,  2 N(gU) = square,     S = m^2 + n^2,
for the five Gaussian linear forms (N = norm from Z[i] to Z; these are v^2 F0, v^2 F8, v^2 G4 / (2S), v^2 GL / 2, v^2 GU / 2)
    a0 = m v + i n u,   a8 = n v + i m u,   g4 = u + i v,   gL = (m + i n) u - i (m - i n) v,   gU = (m + i n) u + i (m - i n) v.
Z[i] is a unique factorisation domain, so each value is  delta_j * (square)  with delta_j squarefree, and a prime that
divides a form to an odd power without dividing its norm condition divides the form and its conjugate, hence their
resultant.  So delta_j lies in a finite set D_j that depends on m, n only:
    delta_0, delta_8 :  i^a d,              d | odd part of mn, squarefree
    delta_4          :  i^a (1+i)^b prod over the primes p | S with odd exponent of one of the two primes above p,
                        b = 1 + v_2(S) mod 2
    delta_L, delta_U :  i^a (1+i) d,        d | odd part of m^2 - n^2, squarefree.
For every rational prime p a solution gives a point (u : v) of P^1(Q_p), and the classes of the five values in
K_v^* / squares at the places v | p of Q(i) must be those of (delta_0, ..., delta_U).  This script computes, for a list
of primes, a SUPERSET W_p of the class tuples that occur for (u, v) in Z_p^2 primitive, and the set T of tuples in
D_0 x ... x D_U that agree with some element of W_p for every p.  If T is empty, the ratio has no solution.
For a ratio with a degenerate solution the tuple of that solution must be in T (a test of the program).

usage: gauss_descent.py [--pmax=N] m,n [m,n ...]"""
import sys, itertools, time
from math import gcd, isqrt
import numpy as np

def small_primes(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b'\0\0'
    for i in range(2, isqrt(n) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]
PR = small_primes(2000000)

def factor(n):
    n = abs(n); f = {}
    for p in PR:
        if p * p > n: break
        while n % p == 0: f[p] = f.get(p, 0) + 1; n //= p
    if n > 1:
        assert n < PR[-1] ** 2, 'factor too large'
        f[n] = f.get(n, 0) + 1
    return f

def gmul(x, y): return (x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0])
def gnorm(x): return x[0]*x[0] + x[1]*x[1]

def gprime(p):
    """A Gaussian prime (a, b) with a^2 + b^2 = p, for p = 1 mod 4."""
    for a in range(1, isqrt(p) + 1):
        b2 = p - a*a; b = isqrt(b2)
        if b*b == b2: return (a, b)
    raise ValueError(p)

def sqrt_m1(p, K):
    """r with r^2 = -1 modulo p^K (p = 1 mod 4)."""
    g = 2
    while pow(g, (p - 1)//2, p) != p - 1: g += 1
    r = pow(g, (p - 1)//4, p); mod = p
    for _ in range(K):                                   # Newton
        mod = min(mod * mod, p ** K); r = (r - (r*r + 1) * pow(2*r, -1, mod)) % mod
    assert (r*r + 1) % p**K == 0
    return r

# ---- local classes.  A class is a small tuple; None means "not determined at this precision".
SQ8 = set()
for a in range(8):
    for b in range(8):
        if (a + b) % 2 == 1:                             # units of Z[i] at 1+i
            z = ((a*a - b*b) % 8, (2*a*b) % 8); SQ8.add(z); SQ8.add(((z[0] + 4) % 8, (z[1] + 4) % 8))
def unit_class2(x):
    """Class of a unit of Z_2[i] modulo squares: the smallest representative modulo 8 of x * (squares mod pi^5)."""
    return min(((x[0]*s[0] - x[1]*s[1]) % 8, (x[0]*s[1] + x[1]*s[0]) % 8) for s in SQ8)

class Local:
    def __init__(s, p, K=None):
        s.p = p
        if p == 2: s.kind = 'ram'
        elif p % 4 == 1:
            s.kind = 'split'; s.K = K or max(8, 200 // p.bit_length()); s.mod = p ** s.K; s.r = sqrt_m1(p, s.K)
        else: s.kind = 'inert'
    def lams(s):
        p = s.p
        if s.kind == 'ram': return [1, 3, 5, 7]
        if s.kind == 'inert': return [1]
        g = 2
        while pow(g, (p - 1)//2, p) != p - 1: g += 1
        return [1, g]
    def cls(s, x, k=None):
        """Class of the Gaussian integer x (x != 0) at the places above p; k: x is only known modulo p^k (None: exact).
        Returns a tuple with one entry per place; an entry is None if it is not determined."""
        p = s.p; a, b = x
        if s.kind == 'ram':
            N = a*a + b*b
            if N == 0: return (None,)
            e = (N & -N).bit_length() - 1                # v_pi(x) = v_2(N(x))
            if k is not None and e + 5 > 2*k: return (None,)
            for _ in range(e): a, b = (a + b)//2, (b - a)//2
            return ((e % 2, unit_class2((a % 8, b % 8))),)
        if s.kind == 'inert':
            if a == 0 and b == 0: return (None,)
            e = 0
            while a % p == 0 and b % p == 0:
                a //= p; b //= p; e += 1
                if k is not None and e >= k: return (None,)
            return ((e % 2, 1 if pow((a*a + b*b) % p, (p - 1)//2, p) == 1 else -1),)
        out = []
        for r in (s.r, s.mod - s.r):
            y = (a + b*r) % s.mod; e = 0
            if y == 0: out.append(None); continue
            while y % p == 0: y //= p; e += 1
            if (k is not None and e >= k) or e >= s.K - 2: out.append(None); continue
            out.append((e % 2, 1 if pow(y % p, (p - 1)//2, p) == 1 else -1))
        return tuple(out)

def forms(m, n):
    """The five forms as pairs (A, B) of Gaussian integers: value = A u + B v."""
    return [((0, n), (m, 0)), ((0, m), (n, 0)), ((1, 0), (0, 1)), ((m, n), (-n, -m)), ((m, n), (n, m))]
def value(f, u, v): return (f[0][0]*u + f[1][0]*v, f[0][1]*u + f[1][1]*v)

def local_set(m, n, p, kmax):
    """A superset of the class tuples of the five forms on primitive (u, v) in Z_p^2, as a set of 5-tuples of classes
    (entries may contain None = any class)."""
    L = Local(p); F = forms(m, n); out = set(); lams = L.lams()
    def emit(u, v, k):
        for lam in lams:
            out.add(tuple(L.cls((lam*x[0], lam*x[1]), k) for x in (value(f, u, v) for f in F)))
    def explore(u, v, k, chart):
        vals = [value(f, u, v) for f in F]
        cl = [L.cls(x, k) for x in vals]
        if k >= kmax or all(None not in c for c in cl): emit(u, v, k); return
        pk = p ** k
        for a in range(p):
            if chart == 0: explore(u, v + pk*a, k + 1, chart)
            else: explore(u + pk*a, v, k + 1, chart)
    for v in range(p): explore(1, v, 1, 0)
    for u in range(0, p, p): explore(u, 1, 1, 1)         # the point (0 : 1) modulo p and its lifts
    return L, out

def dsets(m, n):
    S = m*m + n*n
    def odd_primes(x): return sorted(q for q in factor(x) if q != 2)
    def rational(x):                                     # i^a d, d a squarefree product of the odd primes of x
        ps = odd_primes(x); res = []
        for a in (0, 1):
            for mask in range(1 << len(ps)):
                d = 1
                for i, q in enumerate(ps):
                    if mask >> i & 1: d *= q
                res.append((0, d) if a else (d, 0))
        return res
    D0 = rational(m*n); DLU = [gmul(x, (1, 1)) for x in rational(m*m - n*n)]
    fS = factor(S); ps = sorted(q for q, e in fS.items() if q != 2 and e % 2 == 1)
    b4 = (1 + fS.get(2, 0)) % 2; D4 = []
    for a in (0, 1):
        for ch in itertools.product((0, 1), repeat=len(ps)):
            x = (0, 1) if a else (1, 0)
            if b4: x = gmul(x, (1, 1))
            for q, c in zip(ps, ch):
                g = gprime(q); x = gmul(x, g if c == 0 else (g[0], -g[1]))
            D4.append(x)
    return [D0, list(D0), D4, DLU, list(DLU)]

def match(c, w):
    return all(wi is None or wi == ci for ci, wi in zip(c, w))

def descent(m, n, pmax=400, kmax=5, verbose=False, primes=None):
    assert gcd(m, n) == 1 and m > n >= 1
    S = m*m + n*n; A = m**4 - 2*m**3*n + 2*m*m*n*n + 2*m*n**3 + n**4; B = m**4 + 2*m**3*n + 2*m*m*n*n - 2*m*n**3 + n**4
    bad = sorted(set(factor(2*m*n*(m*m - n*n)*S)) | set(q for q in list(factor(A)) + list(factor(B))))
    D = dsets(m, n)
    T = np.ones(tuple(len(d) for d in D), dtype=bool)
    plist = primes or sorted(set([q for q in bad if q <= 20000] + [q for q in PR if q <= pmax]))
    hist = []
    for p in plist:
        k = kmax if p > 50 else max(kmax, 8 if p == 2 else 6)
        L, W = local_set(m, n, p, k)
        cd = [[L.cls(x) for x in d] for d in D]
        allowed = np.zeros(T.shape, dtype=bool)
        for w in W:
            idx = [[i for i, c in enumerate(cd[j]) if match(c, w[j])] for j in range(5)]
            if all(idx): allowed[np.ix_(*idx)] = True
        T &= allowed
        hist.append((p, int(T.sum())))
        if verbose: print('  p = %d: %d patterns, %d tuples left' % (p, len(W), int(T.sum())), flush=True)
        if not T.any(): break
    return D, T, hist, bad

def degenerate(m, n):
    """(u, v) of the degenerate solutions, if any."""
    S = m*m + n*n; out = []
    if isqrt(S)**2 == S: out += [(1, 1), (1, -1)]
    if S % 2 == 0 and isqrt(S//2)**2 == S//2: out += [(0, 1), (1, 0)]
    return out

if __name__ == '__main__':
    pmax = 400
    for s in sys.argv[1:]:
        if s.startswith('--pmax='): pmax = int(s[7:]); continue
        if s == '-v': continue
        m, n = (int(v) for v in s.split(',')); t0 = time.time()
        D, T, hist, bad = descent(m, n, pmax, verbose=('-v' in sys.argv))
        # test: the tuples of the degenerate solutions must survive
        ok = True
        for (u, v) in degenerate(m, n):
            M = np.ones(T.shape, dtype=bool)             # the tuples with the classes of this solution at every prime used
            for p, _ in hist:
                L = Local(p); c = [L.cls(value(f, u, v)) for f in forms(m, n)]
                cd = [[L.cls(x) for x in d] for d in D]
                ind = np.zeros(T.shape, dtype=bool); idx = [[i for i, cc in enumerate(cd[j]) if cc == c[j]] for j in range(5)]
                if all(idx): ind[np.ix_(*idx)] = True
                M &= ind
            if not (T & M).any(): ok = False
        print('(%d,%d): |D| = %s = %d, tuples left %d after primes up to %d (last change at p = %s)%s  [%.1f s]' % (m, n, 'x'.join(str(len(d)) for d in D), int(np.prod([len(d) for d in D])), int(T.sum()), hist[-1][0],
              max([p for i, (p, c) in enumerate(hist) if i == 0 or c != hist[i-1][1]] or [0]), '' if not degenerate(m, n) else ('; degenerate solutions present: ' + ('kept' if ok else 'LOST (BUG)')), time.time() - t0), flush=True)
