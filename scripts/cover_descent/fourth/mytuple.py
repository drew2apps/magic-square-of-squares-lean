#!/usr/bin/env python3
"""Independent four-class local test (panel 2, agent `tuple`).  Written from the prose only.

A system = four linear forms g_j = a_j*u + b_j*v over Z[i]; form j is ((ar,ai),(br,bi)).
H_p = Q_p^* / (Q_p^* cap squares of Q_p (x) Q(i)), with integer representatives REPS(p):
   inert p : [1, p]            split p : [1, n, p, p*n] (n = least non-residue)       p = 2 : [1, 5, 2, 10]
Local set S_p = { (c_1..c_4) in H_p^4 : exists primitive (u,v) in Z_p^2, c_j * g_j(u,v) a square, all j }.
  upper(...)  : exhaustion modulo p^N, pure table of all squares of Z[i]/p^N           (S_p subset U_p)
  lower(...)  : exact classes at explicit integer pairs (u,v)                           (L_p subset S_p)
"""
import itertools, sys
import numpy as np


def vp(x, p):
    v = 0
    while x % p == 0:
        x //= p; v += 1
    return v


def nonres(p):
    return next(a for a in range(2, p) if pow(a, (p - 1) // 2, p) == p - 1)


def reps(p):
    if p == 2:
        return [1, 5, 2, 10]
    if p % 4 == 3:
        return [1, p]
    n = nonres(p)
    return [1, n, p, p * n]


def unit_reps(p):
    """representatives of Z_p^* modulo squares of Z_p^* (so that lam = rep * mu^2); see proof in REPORT."""
    if p == 2:
        return [1, 3, 5, 7]
    return [1, nonres(p)]


def gpow(x, y, e, p):
    rx, ry = 1, 0
    while e:
        if e & 1:
            rx, ry = (rx * x - ry * y) % p, (rx * y + ry * x) % p
        x, y = (x * x - y * y) % p, (2 * x * y) % p
        e >>= 1
    return rx, ry


_US8 = {((a * a - b * b) % 8, (2 * a * b) % 8) for a in range(8) for b in range(8) if (a + b) % 2}
_iota = {}


def iota(p, K):
    key = (p, K)
    if key not in _iota:
        io = next(a for a in range(2, p) if (a * a + 1) % p == 0)
        P = p
        while P < p ** K:
            P = min(P * P, p ** K)
            io = (io - (io * io + 1) * pow(2 * io, -1, P)) % P
        assert (io * io + 1) % p ** K == 0
        _iota[key] = io
    return _iota[key]


def is_sq(x, y, p):
    """x+iy exact Gaussian integer with non-zero norm: is it a square in Q_p (x) Q(i)?"""
    assert x * x + y * y != 0
    if p == 2:
        w = 0
        while (x - y) % 2 == 0:
            x, y = (x + y) // 2, (y - x) // 2
            w += 1
        return w % 2 == 0 and (x % 8, y % 8) in _US8
    if p % 4 == 3:
        v = min(vp(x, p) if x else 10 ** 9, vp(y, p) if y else 10 ** 9)
        if v % 2:
            return False
        return gpow((x // p ** v) % p, (y // p ** v) % p, (p * p - 1) // 2, p) == (1, 0)
    K = vp(x * x + y * y, p) + 4
    io = iota(p, K)
    P = p ** K
    for c in ((x + io * y) % P, (x - io * y) % P):
        assert c != 0
        v = vp(c, p)
        assert v < K - 1
        if v % 2 or pow((c // p ** v) % p, (p - 1) // 2, p) != 1:
            return False
    return True


def exact_class(x, y, p):
    """index k with REPS[k]*(x+iy) a square in Q_p(x)Q(i), or None.  At most one k can work."""
    if x * x + y * y == 0:
        return None
    hits = [k for k, r in enumerate(reps(p)) if is_sq(r * x, r * y, p)]
    assert len(hits) <= 1, (x, y, p, hits)
    return hits[0] if hits else None


def hclass(d, p):
    k = exact_class(d, 0, p)
    assert k is not None
    return k


def evalform(f, u, v):
    (ar, ai), (br, bi) = f
    return ar * u + br * v, ai * u + bi * v


def exact_tuple(forms, u, v, p):
    out = []
    for f in forms:
        x, y = evalform(f, u, v)
        c = exact_class(x, y, p)
        if c is None:
            return None
        out.append(c)
    return tuple(out)


_tab = {}


def table(M):
    if M not in _tab:
        x = np.arange(M, dtype=np.int64)
        x2 = (x * x) % M
        sq = np.zeros((M, M), bool)
        for a in range(M):
            sq[(x2[a] - x2) % M, (2 * a * x) % M] = True
        _tab[M] = sq
    return _tab[M]


def choose_N(p, cap=5200):
    N = 1
    while p ** (N + 1) <= cap:
        N += 1
    if p == 2:
        N = min(N, 10)
    return N


def decode(code, nr):
    masks = [(code >> (4 * (3 - j))) & 15 for j in range(4)]
    sets = [[k for k in range(nr) if m >> k & 1] for m in masks]
    return list(itertools.product(*sets))


def points(p, N, full=False):
    """(u,v,chart,t,lam) arrays: charts (t,1), (1,t) p|t, times unit reps; or all primitive pairs if full."""
    M = p ** N
    if full:
        a = np.arange(M, dtype=np.int64)
        u, v = np.meshgrid(a, a, indexing='ij')
        u = u.ravel(); v = v.ravel()
        ok = (u % p != 0) | (v % p != 0)
        return u[ok], v[ok], None, None, None
    t = np.arange(M, dtype=np.int64)
    t2 = np.arange(0, M, p, dtype=np.int64)
    U = np.concatenate([t, np.ones(len(t2), np.int64)])
    V = np.concatenate([np.ones(M, np.int64), t2])
    T = np.concatenate([t, t2])
    C = np.concatenate([np.zeros(M, np.int64), np.ones(len(t2), np.int64)])
    us, vs, cs, ts, ls = [], [], [], [], []
    for lam in unit_reps(p):
        us.append((lam * U) % M); vs.append((lam * V) % M); cs.append(C); ts.append(T)
        ls.append(np.full(len(U), lam, np.int64))
    return tuple(np.concatenate(z) for z in (us, vs, cs, ts, ls))


def upper(forms, p, N=None, full=False):
    """returns (set of local tuples, witnesses dict code -> (chart,t,lam)), by table lookups only."""
    if N is None:
        N = choose_N(p)
    M = p ** N
    sq = table(M)
    R = reps(p)
    u, v, C, T, L = points(p, N, full)
    code = np.zeros(len(u), np.int64)
    for f in forms:
        (ar, ai), (br, bi) = f
        gr = ((ar % M) * u + (br % M) * v) % M
        gi = ((ai % M) * u + (bi % M) * v) % M
        mask = np.zeros(len(u), np.int64)
        for k, r in enumerate(R):
            mask |= sq[(r % M) * gr % M, (r % M) * gi % M].astype(np.int64) << k
        code = code * 16 + mask
    uq, idx = np.unique(code, return_index=True)
    S = set(); wit = {}
    for c, i in zip(uq.tolist(), idx.tolist()):
        d = decode(c, len(R))
        if d:
            S.update(d)
            if not full:
                wit[c] = (int(C[i]), int(T[i]), int(L[i]))
    return S, wit, N


def lower(forms, p, wit, N, depth=6, extra=()):
    """exact classes at integer pairs: every witness of the exhaustion, and for ambiguous codes a digit-by-digit
    walk towards the point where a form (or one component) vanishes.  Returns dict tuple -> (u,v)."""
    M = p ** N
    L = {}

    def pair(chart, t, lam):
        return (lam * t, lam) if chart == 0 else (lam, lam * t)

    def tryit(u, v):
        tp = exact_tuple(forms, u, v, p)
        if tp is not None and tp not in L:
            L[tp] = (u, v)
        return tp

    for (u, v) in extra:
        tryit(u, v)
    nr = len(reps(p))
    for c, (chart, t, lam) in wit.items():
        tp0 = tryit(*pair(chart, t, lam))
        masks = [(c >> (4 * (3 - j))) & 15 for j in range(4)]
        if all(bin(m).count('1') == 1 for m in masks) and tp0 is not None and tp0 in decode(c, nr):
            continue
        mod = M // p if N >= 2 else M
        frontier = [t % mod]   # start one level up: also moves sideways inside the residue class mod p^(N-1)
        for e in range(depth + 1):
            cand = []
            for t0 in frontier:
                for k in range(p if p > 2 else 4):
                    t1 = t0 + mod * k
                    u, v = pair(chart, t1, lam)
                    tryit(u, v)
                    sc = 0
                    for f in forms:
                        x, y = evalform(f, u, v)
                        nrm = x * x + y * y
                        sc = max(sc, vp(nrm, p) if nrm else 10 ** 6)
                    cand.append((sc, t1))
            best = max(s for s, _ in cand)
            frontier = [t1 for s, t1 in cand if s == best][:3] + [t1 for s, t1 in cand if s != best][:1]
            mod *= p
    return L


def local_sets(forms, p, N=None, depth=6, extra=()):
    U, wit, N = upper(forms, p, N)
    L = lower(forms, p, wit, N, depth, extra)
    return U, L, N


# ---------- systems ----------
def pf(n):
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


def magic(m, n):
    """g0 = m v + i n u ; g8 = n v - i m u ; gL=(1+i)((mu-nv)+i(mv-nu)) ; gU=(1+i)((mu+nv)+i(mv+nu))"""
    forms = [((0, n), (m, 0)), ((0, -m), (n, 0)),
             ((m + n, m - n), (-n - m, m - n)),      # (1+i)(m - i n) u + (1+i)(-n + i m) v
             ((m - n, m + n), (n - m, m + n))]       # (1+i)(m + i n) u + (1+i)( n + i m) v
    s1 = pf(2 * m * n); s2 = pf(2 * (m * m - n * n))
    return forms, [s1, s1, s2, s2]


def ratdist(r, s):
    forms = [((s, 0), (0, r)), ((s, 0), (-s, r)), ((s, 0), (0, r - s)), ((s, 0), (-s, r - s))]
    s1 = pf(2 * r * s); s2 = pf(2 * s * (r - s))
    return forms, [s1, s1, s2, s2]


def cands(primes):
    out = []
    for mask in range(1 << len(primes)):
        d = 1
        for j, q in enumerate(primes):
            if mask >> j & 1:
                d *= q
        out.append(d)
    return out


if __name__ == '__main__':
    # self-test of the 2-adic unit-square rule: unit squares mod 16, 32 are exactly the lifts of unit squares mod 8
    for M in (16, 32, 64):
        s = {((a * a - b * b) % M, (2 * a * b) % M) for a in range(M) for b in range(M) if (a + b) % 2}
        lifts = {(x, y) for x in range(M) for y in range(M) if (x + y) % 2 and (x % 8, y % 8) in _US8}
        assert s == lifts, M
    print("2-adic unit squares: mod 16,32,64 = lifts of mod 8: OK;  #unit squares mod 8 =", len(_US8))
    # unit reps vs all primitive pairs
    import random
    random.seed(1)
    bad = 0; tot = 0
    for trial in range(60):
        forms = [((random.randint(-6, 6), random.randint(-6, 6)), (random.randint(-6, 6), random.randint(-6, 6)))
                 for _ in range(4)]
        for p, N in ((2, 5), (3, 3), (5, 2), (7, 2), (13, 1), (11, 1)):
            A, _, _ = upper(forms, p, N)
            B, _, _ = upper(forms, p, N, full=True)
            tot += 1
            if A != B:
                bad += 1; print("UNIT-REP MISMATCH", forms, p, N, len(A), len(B))
    print("unit representatives vs all primitive pairs: %d systems x primes, %d mismatches" % (tot, bad))


# ---------- second upper bound: exhaustion modulo p^N without the M x M table (odd p), for larger moduli ----------
def _val(c, p, N):
    """valuation of residues c mod p^N (N for 0), and unit part mod p"""
    v = np.zeros(len(c), np.int64); cc = c.copy()
    for _ in range(N):
        z = (cc % p == 0)
        v += z; cc = np.where(z, cc // p, cc)
    return v, cc % p


def sqmod_fast(x, y, p, N):
    """is x+iy a square in Z[i]/p^N (odd p)?  inert: z=0, or min-valuation even and unit a square of F_p[i] (table p x p);
    split: each component is 0 or has even valuation and unit a QR (CRT Z[i]/p^N = Z/p^N x Z/p^N)."""
    M = p ** N
    if p % 4 == 3:
        vx, _ = _val(x, p, N); vy, _ = _val(y, p, N)
        v = np.minimum(vx, vy)
        pw = np.array([p ** k for k in range(N + 1)], dtype=np.int64)[v]
        ex = np.where(v < N, (x // pw) % p, 0); ey = np.where(v < N, (y // pw) % p, 0)
        key = ('f', p)
        if key not in _tab:
            a = np.arange(p, dtype=np.int64)
            A, B = np.meshgrid(a, a, indexing='ij')
            t = np.zeros((p, p), bool); t[(A * A - B * B) % p, (2 * A * B) % p] = True
            _tab[key] = t
        return (v >= N) | ((v % 2 == 0) & _tab[key][ex, ey])
    io = iota(p, N)
    qr = np.zeros(p, bool); qr[(np.arange(1, p, dtype=np.int64) ** 2) % p] = True
    out = np.ones(len(x), bool)
    for s in (1, -1):
        c = (x + s * ((io * y) % M)) % M
        v, e = _val(c, p, N)
        out &= (v >= N) | ((v % 2 == 0) & qr[e])
    return out


def upper_fast(forms, p, N):
    M = p ** N
    assert p > 2 and M * M < 2 ** 62
    R = reps(p)
    u, v, C, T, L = points(p, N)
    code = np.zeros(len(u), np.int64)
    for f in forms:
        (ar, ai), (br, bi) = f
        gr = ((ar % M) * u + (br % M) * v) % M
        gi = ((ai % M) * u + (bi % M) * v) % M
        mask = np.zeros(len(u), np.int64)
        for k, r in enumerate(R):
            mask |= sqmod_fast((r % M) * gr % M, (r % M) * gi % M, p, N).astype(np.int64) << k
        code = code * 16 + mask
    uq, idx = np.unique(code, return_index=True)
    S = set(); wit = {}
    for c, i in zip(uq.tolist(), idx.tolist()):
        d = decode(c, len(R))
        if d:
            S.update(d); wit[c] = (int(C[i]), int(T[i]), int(L[i]))
    return S, wit, N


def choose_Nfast(p, cap=1500000):
    N = 1
    while p ** (N + 1) <= cap:
        N += 1
    return N


def local_best(forms, p, depth=6, extra=()):
    """U from the table at its N; if L != U, refine U with upper_fast at a larger N (and take the union of the L's)."""
    U, wit, N = upper(forms, p, choose_N(p) if p <= 70 else 1)
    L = lower(forms, p, wit, N, depth, extra)
    how = 'table'
    if set(L) != U and p > 2:
        N2 = choose_Nfast(p)
        if N2 > N:
            U2, wit2, _ = upper_fast(forms, p, N2)
            L2 = lower(forms, p, wit2, N2, depth)
            assert U2 <= U, "deeper exhaustion not inside shallower one"
            L.update(L2); U = U2; N = N2; how = 'fast'
    assert set(L) <= U, ("LOWER NOT IN UPPER", p, set(L) - U)
    return U, L, N, how
