"""Review z8 of roadmap H (written by an AI model; not read by a person; experiment, not counted).
A third method for claim 1 of H3, with no Hensel lemma and no residue symbols beyond 'a unit is a square mod p^2':
for a ratio (m, n) and a prime p = 1 mod 4 that is good for it, go through every point (u, v) of P^1(Z/p^2) (as (x, 1)
and (1, p y)) and both unit classes lam; a label vector e is listed when some point makes lam eps^{e_j} l_j^+-(u, v) a
square mod p^2 (0 counts) for all eight components.  A vector not listed is not allowed at p (a necessary condition).
Integer forms of LAST_PAIR.md s.11 with the actual m, n.
usage: brute_p2.py [m n p] [--write]   (default: (5,1) and (17,1) at p = 1181 -> brute_p2.txt; about 10 s)"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from z8common import emit, ARGS, good, isprime
def run(m, n, p):
    assert isprime(p) and p % 4 == 1 and good(m, n, p)
    M = p * p
    g = 2
    while pow(g, (p - 1) // 2, p) != p - 1: g += 1
    io = pow(g, (p - 1) // 4, p)
    for _ in range(3): io = (io - (io * io + 1) * pow(2 * io, -1, M)) % M
    assert (io * io + 1) % M == 0
    def lin(j, I):
        if j == 0: return I * n, m
        if j == 1: return -I * m, n
        if j == 2: return (1 + I) * (m - I * n), (1 + I) * (-n + I * m)
        return (1 + I) * (m + I * n), (1 + I) * (n + I * m)
    L = [lin(j, s * io) for j in range(4) for s in (1, -1)]
    xs = np.arange(M, dtype=np.int64)
    U = np.concatenate([xs, np.ones(p, dtype=np.int64)]); V = np.concatenate([np.ones(M, dtype=np.int64), p * np.arange(p, dtype=np.int64)])
    def cls(vals):
        vals = vals % M
        c = np.full(vals.shape, 3, dtype=np.int64)          # 3: valuation exactly 1 (never lam * unit * square, never 0)
        c[vals == 0] = 2
        unit = vals % p != 0
        r = vals[unit] % p; e = np.ones_like(r); base = r.copy(); k = (p - 1) // 2
        while k:
            if k & 1: e = e * base % p
            base = base * base % p; k >>= 1
        c[unit] = np.where(e == 1, 0, 1)
        return c
    C = [cls((A % M) * U + (B % M) * V) for (A, B) in L]
    allowed = set()
    for lam in (0, 1):
        for e in range(16):
            ok = np.ones(U.shape, dtype=bool)
            for j in range(4):
                want = lam ^ ((e >> j) & 1)
                for s in range(2): ok &= (C[2 * j + s] == 2) | (C[2 * j + s] == want)
            if ok.any(): allowed.add(e)
    return 'ratio (%d,%d), p = %d: vectors not allowed even mod p^2: %s' % (m, n, p, [e for e in range(16) if e not in allowed])
cases = [tuple(int(a) for a in ARGS[:3])] if len(ARGS) >= 3 else [(5, 1, 1181), (17, 1, 1181)]
emit(['z8: brute force over P^1(Z/p^2) at a good prime above 1153. Written by an AI model; not read by a person.'] + [run(*c) for c in cases], 'brute_p2.txt')
