"""Review z8 of roadmap H (written by an AI model; not read by a person; experiment, not counted).
An independent local test of the four-class descent at a good prime p = 1 mod 4, written from LAST_PAIR.md s.11-12
(not from the bundle's code).  Forms (n = 1 after scaling, m = r mod p), as in tuple_descent.magic:
  g0 = m v + i n u, g8 = n v - i m u, gL = (1+i)((m u - n v) + i(m v - n u)), gU = (1+i)((m u + n v) + i(m v + n u)).
With iota^2 = -1 mod p each g_j splits into l_j^+ = g_j(iota), l_j^- = g_j(-iota).  A tuple of units (d_j) is allowed at
p iff some primitive (u, v) and unit lam make lam d_j l_j^+-(u, v) both squares; its label vector is (chi(d_j)), up to
adding (1,1,1,1).  Enumerate P^1(F_p): a non-root x is allowed iff chi(l_j^+) = chi(l_j^-) for every j, with label
(chi(l_j^+)); at the root x of exactly one l_{j0}^s (a good prime has 8 distinct roots, asserted) l_{j0}^s takes every
value of p Z_p on the disc, so component j0 is chi(l_{j0}^-s(x)).  Rests on: the definition; Legendre symbols; exact
integers (numpy int64, residues < p < 3.0e9).
Run without arguments: compare with H3's h3_weil.json on every good residue at 15 primes from 17 to 1153 (claim 3 of
the review).  Output my_local.txt."""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from z8common import emit, isprime, good, H3

def sqrtm1(p):
    for g in range(2, p):
        z = pow(g, (p - 1) // 4, p)
        if z * z % p == p - 1: return z

def forms(m, n, io, p):
    res = []
    for j in range(4):
        pair = []
        for s in (1, -1):
            I = (s * io) % p
            if j == 0: A, B = I * n, m
            elif j == 1: A, B = -I * m, n
            elif j == 2: A, B = (1 + I) * (m - I * n), (1 + I) * (-n + I * m)
            else: A, B = (1 + I) * (m + I * n), (1 + I) * (n + I * m)
            pair.append((A % p, B % p))
        res.append(pair)
    return res

def chitab(p):
    t = np.ones(p, dtype=np.int64); t[0] = 2
    t[(np.arange(1, p, dtype=np.int64) ** 2) % p] = 0
    return t

def allowed(m, n, p, chi):
    """16-bit mask of the allowed label vectors (bit sum_j e_j 2^j; j = 0, 8, L, U; 0 = square)"""
    io = sqrtm1(p); F = forms(m, n, io, p)
    x = np.arange(p, dtype=np.int64)
    V = [[(F[j][s][0] * x + F[j][s][1]) % p for s in range(2)] for j in range(4)]
    Cs = [[chi[V[j][s]] for s in range(2)] for j in range(4)]          # 0 square, 1 non-square, 2 zero
    nzero = sum((V[j][s] == 0).astype(int) for j in range(4) for s in range(2))
    assert nzero.max() <= 1, 'two roots agree: not a good prime'
    ok = np.ones(p, dtype=bool); lab = np.zeros(p, dtype=np.int64)
    for j in range(4):
        a, b = Cs[j]
        cj = np.where(a == 2, b, a)
        ok &= (a == 2) | (b == 2) | (a == b)
        lab |= cj << j
    vecs = set(int(c) for c in np.unique(lab[ok]))
    Ai = [[F[j][s][0] for s in range(2)] for j in range(4)]            # the point (1, 0)
    def c1(z): return 2 if z == 0 else (0 if pow(z, (p - 1) // 2, p) == 1 else 1)
    assert sum(1 for j in range(4) for s in range(2) if Ai[j][s] == 0) <= 1
    okk = True; l = 0
    for j in range(4):
        a, b = c1(Ai[j][0]), c1(Ai[j][1])
        if a == 2: cj = b
        elif b == 2: cj = a
        elif a == b: cj = a
        else: okk = False; break
        l |= cj << j
    if okk: vecs.add(l)
    vecs |= {15 ^ v for v in vecs}
    mask = 0
    for v in vecs: mask |= 1 << v
    return mask

if __name__ == '__main__':
    T = json.load(open(os.path.join(H3, 'h3_weil.json')))
    out = ['z8: my local test against H3 h3_weil.json (every good residue). Written by an AI model; not read by a person.']
    for p in (17, 29, 37, 41, 53, 61, 97, 101, 137, 181, 229, 409, 617, 1069, 1153):
        ch = chitab(p); d = T[str(p)]; bad = 0; k = 0; lost = 0
        for r in range(1, p):
            if not good(r, 1, p):
                assert str(r) not in d; continue
            a = allowed(r, 1, p, ch); k += 1; bad += a != d[str(r)]; lost += a != (1 << 16) - 1
        assert k == len(d)
        out.append('p = %d: %d good residues, %d lose a vector, %d disagree with h3_weil.json' % (p, k, lost, bad))
    emit(out, 'my_local.txt')
