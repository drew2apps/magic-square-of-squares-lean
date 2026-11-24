#!/usr/bin/env python3
"""A third implementation of the local test of LAST_PAIR.md, section 11 (Theorem D (c)) and section 12 (tuples).

Provenance: written on 23 November 2026 by a separate instance of Claude (an AI model made by Anthropic) acting as a
reviewer, from the prose of sections 11 and 14 only, before it had opened scripts/cover_descent; copied here unchanged
except for this paragraph.  It shares no code with cover_descent.py, tuple_descent.py or the two exhaustion programs.
It gives TWO bounds for the set of classes possible at a prime: an upper bound by exhaustion modulo p^N (as
cover_bruteforce.py does), and a certified lower bound: the classes for which an explicit integer point (u, v) makes
every condition hold exactly in Q_p and Q_p(i).  Where the two agree the local set is known exactly.  Exact integers.
No valuations / residue symbols / Hensel in the UPPER bound (pure table lookups modulo p^N).
The LOWER bound (certified local points) evaluates the forms at honest integer points and
decides their exact square class in Q_p or Q_p(i) with textbook facts (stated in the code).

Notation of the prose (S = m^2+n^2, D = m^2-n^2):
  g0 = m v + i n u          N = F0 = m^2 v^2 + n^2 u^2
  g8 = n v - i m u          N = F8 = n^2 v^2 + m^2 u^2
  gL = (1+i)((mu-nv)+i(mv-nu))   N = GL = 2((mu-nv)^2+(mv-nu)^2)
  gU = (1+i)((mu+nv)+i(mv+nu))   N = GU = 2((mu+nv)^2+(mv+nu)^2)
  A  = -S u^2 + 2 D i u v + S v^2,  B = g0*g8,  B' = conj(g0)*conj(g8)
"""
import sys, itertools, math
import numpy as np


def prime_factors(n):
    n = abs(n); out = []; q = 2
    while q * q <= n:
        if n % q == 0:
            out.append(q)
            while n % q == 0:
                n //= q
        q += 1
    if n > 1:
        out.append(n)
    return out


def cands(primes):
    """all squarefree positive integers supported on `primes`, index = bitmask"""
    out = []
    for mask in range(1 << len(primes)):
        d = 1
        for j, q in enumerate(primes):
            if mask >> j & 1:
                d *= q
        out.append(d)
    return out


_tab = {}


def tables(M):
    """(all squares of Z/M, all squares of Z[i]/M) by squaring every element."""
    if M in _tab:
        return _tab[M]
    x = np.arange(M, dtype=np.int64)
    sqZ = np.zeros(M, bool)
    sqZ[(x * x) % M] = True
    sqG = np.zeros((M, M), bool)
    x2 = (x * x) % M
    for a in range(M):
        re = (x2[a] - x2) % M
        im = (2 * a * x) % M
        sqG[re, im] = True
    _tab[M] = (sqZ, sqG)
    return _tab[M]


def gmul(a, b, M):
    return ((a[0] * b[0] - a[1] * b[1]) % M, (a[0] * b[1] + a[1] * b[0]) % M)


def proj_points(p, N):
    M = p ** N
    t = np.arange(M, dtype=np.int64)
    t2 = np.arange(0, M, p, dtype=np.int64)
    u = np.concatenate([t, np.ones(len(t2), np.int64)])
    v = np.concatenate([np.ones(M, np.int64), t2])
    return u, v


def forms(m, n, u, v, M):
    """the four Gaussian linear forms (with 1+i for L, U), A, and the 4 norms, mod M (arrays)"""
    m %= M; n %= M
    g0 = ((m * v) % M, (n * u) % M)
    g8 = ((n * v) % M, (-m * u) % M)
    lr, li = (m * u - n * v) % M, (m * v - n * u) % M
    gL = ((lr - li) % M, (lr + li) % M)
    ur, ui = (m * u + n * v) % M, (m * v + n * u) % M
    gU = ((ur - ui) % M, (ur + ui) % M)
    S = (m * m + n * n) % M; D = (m * m - n * n) % M
    A = ((S * ((v * v - u * u) % M)) % M, (2 * D * ((u * v) % M)) % M)
    F0 = (g0[0] * g0[0] + g0[1] * g0[1]) % M
    F8 = (g8[0] * g8[0] + g8[1] * g8[1]) % M
    GL = (gL[0] * gL[0] + gL[1] * gL[1]) % M
    GU = (gU[0] * gU[0] + gU[1] * gU[1]) % M
    return g0, g8, gL, gU, A, (F0, F8, GL, GU)


def cover_primes(m, n):
    return prime_factors(2 * m * n * (m * m - n * n))


def cover_upper(m, n, p, N, cl=None):
    """Upper bound for Loc_p: boolean array over candidates d (bitmask index over cover_primes).
    d kept iff at some point of P^1(Z/p^N): F0,F8,GL,GU squares mod M and d*A*B, d*A*B' squares in Z[i]/M."""
    M = p ** N
    sqZ, sqG = tables(M)
    if cl is None:
        cl = cands(cover_primes(m, n))
    u, v = proj_points(p, N)
    g0, g8, gL, gU, A, (F0, F8, GL, GU) = forms(m, n, u, v, M)
    ok = sqZ[F0] & sqZ[F8] & sqZ[GL] & sqZ[GU]
    u, v = u[ok], v[ok]
    g0 = (g0[0][ok], g0[1][ok]); g8 = (g8[0][ok], g8[1][ok]); A = (A[0][ok], A[1][ok])
    B = gmul(g0, g8, M)
    Bp = (B[0], (-B[1]) % M)
    AB = gmul(A, B, M); ABp = gmul(A, Bp, M)
    keep = np.zeros(len(cl), bool)
    wit = {}
    for j, d in enumerate(cl):
        dm = d % M
        hit = sqG[(dm * AB[0]) % M, (dm * AB[1]) % M] & sqG[(dm * ABp[0]) % M, (dm * ABp[1]) % M]
        if hit.any():
            keep[j] = True
            wit[j] = np.nonzero(hit)[0]
    return keep, (u, v, wit)


# ---------- exact local square classes (for certified LOWER bounds only) ----------
def vp(x, p):
    v = 0
    while x % p == 0:
        x //= p; v += 1
    return v


def is_sq_Qp(x, p):
    """x nonzero integer. Textbook: x = p^v e is a square in Q_p iff v even and e is a QR mod p
    (p odd) / e = 1 mod 8 (p = 2)."""
    assert x != 0
    v = vp(x, p)
    if v % 2:
        return False
    e = x // p ** v
    if p == 2:
        return e % 8 == 1
    return pow(e % p, (p - 1) // 2, p) == 1


_unit_sq8 = None
_fp2 = {}


def is_sq_Qpi(x, y, p, K=60):
    """x+iy nonzero Gaussian integer; is it a square in Q_p (x) Q(i)?  Returns True/False.
    inert p: Q_p(i) unramified field; p^v e square iff v even and e a square in F_p[i] (table of all squares).
    split p: both components x +- iota*y squares in Q_p (iota^2 = -1 computed mod p^K; K large so the
             components are known exactly enough; assert valuation < K-2).
    p = 2 : pi = 1+i; pi^w e square iff w even and e mod 8 is the square of a unit mod 8
            (1 + pi^5 Z_2[i] consists of squares; 8 = pi^6 up to a unit)."""
    assert (x, y) != (0, 0)
    if p == 2:
        global _unit_sq8
        if _unit_sq8 is None:
            _unit_sq8 = set()
            for a in range(8):
                for b in range(8):
                    if (a + b) % 2 == 1:
                        _unit_sq8.add(((a * a - b * b) % 8, (2 * a * b) % 8))
        w = 0
        while (x - y) % 2 == 0:
            x, y = (x + y) // 2, (y - x) // 2   # divide by 1+i
            w += 1
        if w % 2:
            return False
        return (x % 8, y % 8) in _unit_sq8
    if p % 4 == 3:
        v = min(vp(x, p) if x else 10 ** 9, vp(y, p) if y else 10 ** 9)
        if v % 2:
            return False
        if p not in _fp2:
            s = set()
            for a in range(p):
                for b in range(p):
                    if (a, b) != (0, 0):
                        s.add(((a * a - b * b) % p, (2 * a * b) % p))
            _fp2[p] = s
        return ((x // p ** v) % p, (y // p ** v) % p) in _fp2[p]
    # split
    P = p ** K
    io = next(a for a in range(2, p) if (a * a + 1) % p == 0)
    for _ in range(8):
        io = (io - (io * io + 1) * pow(2 * io, -1, P)) % P
    assert (io * io + 1) % P == 0
    for c in ((x + io * y) % P, (x - io * y) % P):
        assert c != 0
        v = vp(c, p)
        assert v < K - 2, "precision"
        if v % 2:
            return False
        if pow((c // p ** v) % p, (p - 1) // 2, p) != 1:
            return False
    return True


def exact_forms(m, n, u, v):
    g0 = (m * v, n * u)
    g8 = (n * v, -m * u)
    lr, li = m * u - n * v, m * v - n * u
    gL = (lr - li, lr + li)
    ur, ui = m * u + n * v, m * v + n * u
    gU = (ur - ui, ur + ui)
    S = m * m + n * n; D = m * m - n * n
    A = (S * (v * v - u * u), 2 * D * u * v)
    return g0, g8, gL, gU, A


def zmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def cover_point_ok(m, n, u, v, d, p):
    """exact: is integer point (u:v) a Q_p-point in the sense of Theorem D(c) for d?"""
    g0, g8, gL, gU, A = exact_forms(m, n, u, v)
    for g in (g0, g8, gL, gU):
        if not is_sq_Qp(g[0] ** 2 + g[1] ** 2, p):
            return False
    B = zmul(g0, g8); Bp = (B[0], -B[1])
    AB = zmul(A, B); ABp = zmul(A, Bp)
    return is_sq_Qpi(d * AB[0], d * AB[1], p) and is_sq_Qpi(d * ABp[0], d * ABp[1], p)


def cover_lower(m, n, p, N, keep, aux, cl):
    """certified subset of Loc_p: d with an honest integer point (lift of a point mod p^N)."""
    u, v, wit = aux
    low = np.zeros(len(cl), bool)
    for j in np.nonzero(keep)[0]:
        for idx in wit[j][:400]:
            if cover_point_ok(m, n, int(u[idx]), int(v[idx]), cl[j], p):
                low[j] = True
                break
    return low


# ---------------- tuples (section 12) ----------------
def tuple_primes(m, n):
    return prime_factors(2 * m * n), prime_factors(2 * (m * m - n * n))


def unit_reps(p):
    if p == 2:
        return [1, 3, 5, 7]
    nr = next(a for a in range(2, p) if pow(a, (p - 1) // 2, p) != 1)
    return [1, nr]


def tuple_upper(m, n, p, N):
    """boolean array [d0, d8, dL, dU]: tuple kept iff some primitive (u, v) mod p^N has
    d_j * gamma_j(u, v) a square of Z[i]/p^N for all four j.
    primitive pairs = projective representatives times units; a unit only matters modulo squares of
    (Z/p^N)^*, so representatives of (Z/M)^*/(Z/M)^*2 are enough (exact, no theory beyond that)."""
    M = p ** N
    sqZ, sqG = tables(M)
    P1, P2 = tuple_primes(m, n)
    c1, c2 = cands(P1), cands(P2)
    u0, v0 = proj_points(p, N)
    out = np.zeros((len(c1), len(c1), len(c2), len(c2)), bool)
    sigs = {}
    for lam in unit_reps(p):
        u, v = (lam * u0) % M, (lam * v0) % M
        g0, g8, gL, gU, A, _ = forms(m, n, u, v, M)
        masks = []
        for g, cc in ((g0, c1), (g8, c1), (gL, c2), (gU, c2)):
            mk = np.zeros((len(cc), len(u)), bool)
            for j, d in enumerate(cc):
                dm = d % M
                mk[j] = sqG[(dm * g[0]) % M, (dm * g[1]) % M]
            masks.append(mk)
        alive = masks[0].any(0) & masks[1].any(0) & masks[2].any(0) & masks[3].any(0)
        for idx in np.nonzero(alive)[0]:
            key = tuple(mk[:, idx].tobytes() for mk in masks)
            if key not in sigs:
                sigs[key] = (int(u[idx]), int(v[idx]))
    for key, pt in sigs.items():
        a, b, c, d = [np.frombuffer(k, bool) for k in key]
        out |= a[:, None, None, None] & b[None, :, None, None] & c[None, None, :, None] & d[None, None, None, :]
    return out, sigs, (c1, c2)


def tuple_point_ok(m, n, u, v, tup, p):
    g0, g8, gL, gU, A = exact_forms(m, n, u, v)
    for g, d in zip((g0, g8, gL, gU), tup):
        if not is_sq_Qpi(d * g[0], d * g[1], p):
            return False
    return True


def class_of(g):
    """Lemma S1: g Gaussian integer with square norm w^2 -> squarefree positive d with g = d z^2."""
    x, y = g
    w = math.isqrt(x * x + y * y)
    assert w * w == x * x + y * y
    r = 2 * (x + w) if (y != 0 or x + w != 0) else x
    if y == 0 and x + w == 0:
        r = x
    r = abs(r)
    d = 1
    for q in prime_factors(r):
        if vp(r, q) % 2:
            d *= q
    return d


def known_points(m, n):
    S = m * m + n * n
    pts = []
    if math.isqrt(S) ** 2 == S:
        pts += [(1, 1), (1, -1)]
    if S % 2 == 0 and math.isqrt(S // 2) ** 2 == S // 2:
        pts += [(1, 0), (0, 1)]
    return pts


def known_tuples(m, n):
    out = set()
    for (u, v) in known_points(m, n):
        g0, g8, gL, gU, A = exact_forms(m, n, u, v)
        out.add(tuple(class_of(g) for g in (g0, g8, gL, gU)))
    return out


def known_cover_classes(m, n):
    out = set()
    for (u, v) in known_points(m, n):
        g0, g8, gL, gU, A = exact_forms(m, n, u, v)
        AB = zmul(A, zmul(g0, g8))
        out.add(class_of(AB))
    return out


DEFAULT_N = {2: 8, 3: 6, 5: 4, 7: 4, 11: 3, 13: 3}
PRIMES61 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61]


def depth(p, m, n, big=False):
    if p in DEFAULT_N:
        N = DEFAULT_N[p]
    elif p <= 17:
        N = 3 if big else 2
    else:
        N = 2
    return N
