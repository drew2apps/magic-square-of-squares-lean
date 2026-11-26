"""Descent by a cyclic 4-isogeny over K = Q(i) (mu_4 in K, so it is Kummer theory in K*/K*^4).

psi: A -> B with ker psi = <Q>, Q in A(K) of order 4;  ker psihat = <T'>, T' in B(K) of order 4.
Since i in K, ker psi ~ Z/4 ~ mu_4 and ker psihat ~ mu_4, and the connecting maps are
    delta_psi:    B(K)/psi A(K)    -> K*/K*^4,  P -> g_{T'}(P)
    delta_psihat: A(K)/psihat B(K) -> K*/K*^4,  P -> g_Q(P)
with div g_T = 4(T) - 4(O):  g_T = l_T^2 / (x - x(2T)), l_T the tangent line at T (g_T o psi is a 4th power).
Local images: sampled from local points; sampling stops only at the exact product
    |im_v psi| * |im_v psihat| = 16 (v odd),  256 (v = 1+i)
(from 0 -> B_v[psihat]/psi(A_v[4]) -> B_v/psi A_v -> A_v/4A_v -> A_v/psihat B_v -> 0 and |A_v/4A_v| = |A_v[4]| |4|_v^-1),
and at odd v each is also checked against |im_v psi| = 4 c_v(B)/c_v(A).
Global: 4^r = |B/psiA| |A/psihat B| / 16, so  r <= log_4(|Sel^psi| |Sel^psihat| / 16).
Written by an AI model (Claude); not read by a person.  EXPERIMENT."""
import sys, random, itertools
from fractions import Fraction as F
import cypari2
pari = cypari2.Pari(); pari.allocatemem(4 * 10**9)

# ---------------- Gaussian rationals as (F, F)
def qa(x, y): return (x[0] + y[0], x[1] + y[1])
def qs(x, y): return (x[0] - y[0], x[1] - y[1])
def qm(x, y): return (x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0])
def qi(x):
    n = x[0]*x[0] + x[1]*x[1]; return (x[0]/n, -x[1]/n)
def qd(x, y): return qm(x, qi(y))
def qc(c): return (F(c), F(0))
def gq(z):
    z = pari(z)
    return (F(str(pari.polcoef(z, 0, 'y'))), F(str(pari.polcoef(z, 1, 'y'))))

def gprime(p):
    for a in range(1, 10**7):
        b2 = p - a*a
        if b2 < 0: break
        b = pari.sqrtint(b2)
        if int(b)**2 == b2: return (a, int(b))
    raise ValueError

UNK = 'UNK'

def hensel_r(p, a, b, N):
    M = p**N
    r = (-a * pow(b, -1, p)) % p
    mod = p
    while mod < M:
        mod = min(mod*mod, M); r = (r - (r*r + 1) * pow(2*r, -1, mod)) % mod
    return r

class Split:
    """K_v = Q_p, i -> r with a + b r = 0 (place of the Gaussian prime a + b i)."""
    def __init__(s, p, pi, bits=600):
        s.p = p; s.N = max(12, bits // p.bit_length()); s.M = p**s.N
        s.r = hensel_r(p, pi[0], pi[1], s.N); s.name = f"split {p} ({pi[0]}{pi[1]:+d}i)"; s.dims = [4, 4]
        s.rp = s.r % p; s.unif = (F(pi[0]), F(pi[1]))
        s.mu = {1: 0, s.rp: 1, p - 1: 2, (p - s.rp) % p: 3}
    def _vint(s, n, cap):
        e = 0
        while n % s.p == 0 and e < cap: n //= s.p; e += 1
        return e, n
    def from_q(s, x):
        if x[0] == 0 and x[1] == 0: return None
        A, D1 = x[0].numerator, x[0].denominator; B, D2 = x[1].numerator, x[1].denominator
        M2 = s.p**(2*s.N)
        n = (A*D2 + B*D1*s.r) % M2
        if n == 0: return UNK
        e1, u1 = s._vint(n, 2*s.N)
        e2, u2 = s._vint(D1*D2, 10**9)
        return (e1 - e2, (u1 * pow(u2, -1, s.M)) % s.M, min(s.N, 2*s.N - e1))
    def mul(s, x, y):
        if x is None or y is None: return None
        if x is UNK or y is UNK: return UNK
        return (x[0] + y[0], x[1]*y[1] % s.M, min(x[2], y[2]))
    def inv(s, x):
        if x is None or x is UNK: return UNK
        return (-x[0], pow(x[1], -1, s.M), x[2])
    def neg(s, x):
        if x is None or x is UNK: return x
        return (x[0], (-x[1]) % s.M, x[2])
    def add(s, x, y):
        if x is UNK or y is UNK: return UNK
        if x is None: return y
        if y is None: return x
        if x[0] > y[0]: x, y = y, x
        d = y[0] - x[0]
        pr = min(x[2], y[2] + d)
        if pr <= 0: return UNK
        m = s.p**pr
        t = (x[1] + y[1] * s.p**d) % m
        if t == 0: return UNK
        v, u = s._vint(t, pr)
        return (x[0] + v, u % s.M, pr - v)
    def sqrt(s, x):
        if x is None: return None
        e, u, pr = x
        if e % 2 or pow(u % s.p, (s.p - 1)//2, s.p) != 1: return False
        t = int(pari.lift(pari.sqrt(pari.Mod(u % s.p, s.p))))
        m = s.p
        while m < s.M:
            m = min(m*m, s.M); t = (t - (t*t - u) * pow(2*t, -1, m)) % m
        return (e // 2, t, pr)
    def cls(s, x):
        if x is None or x is UNK or x[2] < 1: return None
        return (x[0] % 4, s.mu[pow(x[1] % s.p, (s.p - 1)//4, s.p)])
    def rand(s, rng):
        j = rng.randint(-6, 8); t = rng.randrange(1, s.M)
        return s.mul((j, 1, s.N), s.from_q((F(t), F(0))))

def gmm(x, y, M): return ((x[0]*y[0] - x[1]*y[1]) % M, (x[0]*y[1] + x[1]*y[0]) % M)
def gpw(x, k, M):
    r = (1, 0)
    while k:
        if k & 1: r = gmm(r, x, M)
        x = gmm(x, x, M); k >>= 1
    return r
def ginv(x, M):
    n = (x[0]*x[0] + x[1]*x[1]) % M; ni = pow(n, -1, M)
    return (x[0]*ni % M, (-x[1])*ni % M)

class Inert:
    def __init__(s, p, bits=600):
        s.p = p; s.N = max(12, bits // p.bit_length()); s.M = p**s.N; s.name = f"inert {p}"; s.dims = [4, 4]
        s.unif = (F(p), F(0))
        s.mu = {(1, 0): 0, (0, 1): 1, (p - 1, 0): 2, (0, p - 1): 3}
    def _fsqrt(s, a, b):
        p = s.p
        def rt(z):
            r = pow(z, (p + 1)//4, p); return r if r*r % p == z % p else None
        if b == 0:
            r = rt(a)
            if r is not None: return (r, 0)
            r = rt(-a % p); return (0, r)          # (r i)^2 = -r^2 = a
        n = rt((a*a + b*b) % p)
        for nn in (n, (-n) % p):
            x2 = (a + nn) * pow(2, -1, p) % p
            x = rt(x2)
            if x is not None and x != 0:
                return (x, b * pow(2*x, -1, p) % p)
        raise ValueError("no sqrt")
    def _vg(s, a, b, cap):
        e = 0
        while a % s.p == 0 and b % s.p == 0 and e < cap: a //= s.p; b //= s.p; e += 1
        return e, a, b
    def from_q(s, x):
        if x[0] == 0 and x[1] == 0: return None
        A, D1 = x[0].numerator, x[0].denominator; B, D2 = x[1].numerator, x[1].denominator
        e1, a, b = s._vg(A*D2, B*D1, 10**9)
        e2, d, _ = s._vg(D1*D2, 0, 10**9)
        di = pow(d, -1, s.M)
        return (e1 - e2, (a*di % s.M, b*di % s.M), s.N)
    def mul(s, x, y):
        if x is None or y is None: return None
        if x is UNK or y is UNK: return UNK
        return (x[0] + y[0], gmm(x[1], y[1], s.M), min(x[2], y[2]))
    def inv(s, x):
        if x is None or x is UNK: return UNK
        return (-x[0], ginv(x[1], s.M), x[2])
    def neg(s, x):
        if x is None or x is UNK: return x
        return (x[0], ((-x[1][0]) % s.M, (-x[1][1]) % s.M), x[2])
    def add(s, x, y):
        if x is UNK or y is UNK: return UNK
        if x is None: return y
        if y is None: return x
        if x[0] > y[0]: x, y = y, x
        d = y[0] - x[0]; pr = min(x[2], y[2] + d)
        if pr <= 0: return UNK
        m = s.p**pr; pd = s.p**d
        a = (x[1][0] + y[1][0]*pd) % m; b = (x[1][1] + y[1][1]*pd) % m
        if a == 0 and b == 0: return UNK
        v, a, b = s._vg(a, b, pr)
        return (x[0] + v, (a % s.M, b % s.M), pr - v)
    def sqrt(s, x):
        if x is None: return None
        e, u, pr = x
        if e % 2: return False
        p = s.p
        if gpw((u[0] % p, u[1] % p), (p*p - 1)//2, p) != (1, 0): return False
        t = s._fsqrt(u[0] % p, u[1] % p)
        m = p
        while m < s.M:
            m = min(m*m, s.M)
            t2 = gmm(t, t, m); num = ((t2[0] - u[0]) % m, (t2[1] - u[1]) % m)
            corr = gmm(num, ginv(((2*t[0]) % m, (2*t[1]) % m), m), m)
            t = ((t[0] - corr[0]) % m, (t[1] - corr[1]) % m)
        return (e // 2, t, pr)
    def cls(s, x):
        if x is None or x is UNK or x[2] < 1: return None
        p = s.p
        return (x[0] % 4, s.mu[gpw((x[1][0] % p, x[1][1] % p), (p*p - 1)//4, p)])
    def rand(s, rng):
        j = rng.randint(-6, 8); t = (rng.randrange(s.M), rng.randrange(s.M))
        if t[0] % s.p == 0 and t[1] % s.p == 0: t = (t[0] + 1, t[1])
        return (j, t, s.N)

# ---------- the place above 2:  pi = 1+i, units mod 2^N (= pi^(2N) up to a unit); precision in pi-digits
def pidiv(a, b): return ((a + b) // 2, (b - a) // 2)          # (a+bi)/(1+i), exact when a+b even
def red9(a, b):
    t = (b - b % 16) // 16; return ((a - 16*t) % 32, b % 16)    # canonical rep mod (1+i)^9 = 16(1+i)

def _ram_table():
    G = sorted(set(red9(a, b) for a in range(32) for b in range(16) if (a + b) % 2))
    def m9(x, y): return red9(x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0])
    def p4(x): return m9(m9(x, x), m9(x, x))
    H = sorted(set(p4(g) for g in G))
    def coset(x): return min(m9(x, h) for h in H)
    cos = sorted(set(coset(g) for g in G)); assert len(cos) == 64, len(cos)
    def pw(x, k):
        r = (1, 0)
        for _ in range(k): r = m9(r, x)
        return r
    for g1, g2, g3 in itertools.combinations(cos, 3):
        tab = {}
        for a in range(4):
            for b in range(4):
                for c in range(4):
                    tab[coset(m9(m9(pw(g1, a), pw(g2, b)), pw(g3, c)))] = (a, b, c)
        if len(tab) == 64:
            return {g: tab[coset(g)] for g in G}
    raise RuntimeError
RAMTAB = None

class Ram:
    def __init__(s, bits=600):
        global RAMTAB
        if RAMTAB is None: RAMTAB = _ram_table()
        s.p = 2; s.N = bits; s.M = 2**bits; s.P = 2*bits; s.name = "ram 1+i"; s.dims = [4, 4, 4, 4]
        s.unif = (F(1), F(1))
    def _vg(s, a, b, cap):
        e = 0
        while (a + b) % 2 == 0 and e < cap:
            if a == 0 and b == 0: return cap, 0, 0
            a, b = pidiv(a, b); e += 1
        return e, a, b
    def from_q(s, x):
        if x[0] == 0 and x[1] == 0: return None
        A, D1 = x[0].numerator, x[0].denominator; B, D2 = x[1].numerator, x[1].denominator
        e1, a, b = s._vg(A*D2, B*D1, 10**9)
        D = D1*D2; k = 0
        while D % 2 == 0: D //= 2; k += 1
        u = gmm((a, b), (pow(D, -1, s.M), 0), s.M)
        u = gmm(u, gpw((0, 1), k, s.M), s.M)                     # 1/2^k = i^k pi^(-2k)
        return (e1 - 2*k, u, s.P)
    def mul(s, x, y):
        if x is None or y is None: return None
        if x is UNK or y is UNK: return UNK
        return (x[0] + y[0], gmm(x[1], y[1], s.M), min(x[2], y[2]))
    def inv(s, x):
        if x is None or x is UNK: return UNK
        return (-x[0], ginv(x[1], s.M), x[2])
    def neg(s, x):
        if x is None or x is UNK: return x
        return (x[0], ((-x[1][0]) % s.M, (-x[1][1]) % s.M), x[2])
    def add(s, x, y):
        if x is UNK or y is UNK: return UNK
        if x is None: return y
        if y is None: return x
        if x[0] > y[0]: x, y = y, x
        d = y[0] - x[0]; pr = min(x[2], y[2] + d)
        if pr <= 0: return UNK
        pd = gpw((1, 1), d, s.M)
        t = gmm(y[1], pd, s.M); a = (x[1][0] + t[0]) % s.M; b = (x[1][1] + t[1]) % s.M
        v, a, b = s._vg(a, b, pr)
        if v >= pr: return UNK
        return (x[0] + v, (a % s.M, b % s.M), pr - v)
    def sqrt(s, x):
        if x is None: return None
        e, u, pr = x
        if e % 2: return False
        M = s.M
        def vdiff(t):
            t2 = gmm(t, t, M); return s._vg((t2[0] - u[0]) % M, (t2[1] - u[1]) % M, pr)
        t = None
        for a in range(8):
            for b in range(8):
                if (a + b) % 2 and vdiff((a, b))[0] >= 5: t = (a, b); break
            if t: break
        if t is None: return False
        for _ in range(4000):
            k, da, db = vdiff(t)
            if k >= pr - 2: break
            # t <- t - (t^2-u)/(2t);  (t^2-u) = pi^k * (da + db i);  2t = -i pi^2 t
            cor = gmm(gmm((da, db), gpw((1, 1), k - 2, M), M), gmm((0, 1), ginv(t, M), M), M)
            t = ((t[0] - cor[0]) % M, (t[1] - cor[1]) % M)
        return (e // 2, t, pr - 2)
    def cls(s, x):
        if x is None or x is UNK or x[2] < 9: return None
        return (x[0] % 4,) + RAMTAB[red9(x[1][0], x[1][1])]
    def rand(s, rng):
        j = rng.randint(-8, 12); t = (rng.randrange(s.M), rng.randrange(s.M))
        if (t[0] + t[1]) % 2 == 0: t = (t[0] + 1, t[1])
        return (j, t, s.P)

# ---------------- the descent map on a curve Y^2 = X^3 + c2 X^2 + c4 X + c6 for a K-point T of order 4
class Kummer4:
    def __init__(s, c, T, tors):
        s.c = c; s.T = T
        xT, yT = T
        s.lam = qd(qa(qa(qm(qc(3), qm(xT, xT)), qm(qm(qc(2), c[0]), xT)), c[1]), qm(qc(2), yT))
        s.x2 = qs(qs(qm(s.lam, s.lam), c[0]), qm(qc(2), xT))
        s.tors = tors                          # global torsion points (affine), used as sampling centres
    def local(s, pl):
        L = {k: pl.from_q(v) for k, v in (('c2', s.c[0]), ('c4', s.c[1]), ('c6', s.c[2]), ('xT', s.T[0]), ('yT', s.T[1]),
                                        ('lam', s.lam), ('x2', s.x2))}
        L['ctr'] = [pl.from_q(P[0]) for P in s.tors]
        return L
    def value(s, pl, L, X):
        f = pl.add(pl.mul(pl.add(pl.mul(pl.add(X, L['c2']), X), L['c4']), X), L['c6'])
        if f is None or f is UNK: return None
        Y = pl.sqrt(f)
        if Y is False or Y is None or Y is UNK: return None
        cc = pl.add(L['yT'], pl.mul(L['lam'], pl.add(X, pl.neg(L['xT']))))
        d = pl.add(Y, pl.neg(cc))
        den = pl.add(X, pl.neg(L['x2']))
        if d is None or d is UNK or den is None or den is UNK: return None
        g = pl.mul(pl.mul(d, d), pl.inv(den))
        return pl.cls(g)

# ---------------- Z/4 linear algebra via PARI
def span_size(vecs, d):
    if not vecs: return 1
    M = pari.matconcat([pari.matrix(d, len(vecs), [vecs[j][i] for i in range(d) for j in range(len(vecs))]), 4 * pari.matid(d)])
    H = pari.mathnf(M)
    return 4**d // int(pari.matdet(H))
def reduce_basis(vecs, d):
    if not vecs: return []
    M = pari.matconcat([pari.matrix(d, len(vecs), [vecs[j][i] for i in range(d) for j in range(len(vecs))]), 4 * pari.matid(d)])
    H = pari.mathnf(M)
    return [[int(H[i, j]) % 4 for i in range(d)] for j in range(int(pari.matsize(H)[1]))]

def sample_images(pl, maps, target, rng, budget=60000, verbose=False):
    """maps: list of two Kummer4 (one per side). Returns list of two generator lists."""
    Ls = [m.local(pl) for m in maps]
    gens = [[], []]; sizes = [1, 1]; d = sum(1 for _ in pl.dims)
    for t in range(budget):
        if sizes[0] * sizes[1] == target: break
        side = t % 2; m = maps[side]; L = Ls[side]
        mode = (t // 2) % 3
        if mode == 0 or not L['ctr']:
            X = pl.rand(rng)
        else:
            ctr = L['ctr'][rng.randrange(len(L['ctr']))]
            if ctr is None or ctr is UNK: continue
            off = pl.rand(rng)
            if off is None: continue
            off = (off[0] + rng.randint(1, 12) + max(0, ctr[0]), off[1], off[2])
            X = pl.add(ctr, off)
            if X is UNK or X is None: continue
        c = m.value(pl, L, X)
        if c is None: continue
        cand = gens[side] + [list(c)]
        sz = span_size(cand, d)
        if sz > sizes[side]:
            gens[side] = reduce_basis(cand, d); sizes[side] = sz
            if verbose: print(f"     {pl.name} side {side}: {sz}", flush=True)
    return gens, sizes

def selmer(globs, places, gens, d_list):
    """globs: list of global basis elements (Gaussian rationals); returns |Sel| and its invariants."""
    g = len(globs)
    rows = []; blocks = []
    Mrows = []
    for pl, G, d in zip(places, gens, d_list):
        cl = [pl.cls(pl.from_q(z)) for z in globs]
        assert None not in cl, pl.name
        for i in range(d): Mrows.append([cl[k][i] for k in range(g)])
        blocks.append((d, G))
    R = len(Mrows)
    ncols = g + sum(len(G) for _, G in blocks) + R
    Mat = [[0]*ncols for _ in range(R)]
    for i in range(R):
        for k in range(g): Mat[i][k] = Mrows[i][k]
    r0 = 0; c0 = g
    for d, G in blocks:
        for j, v in enumerate(G):
            for i in range(d): Mat[r0 + i][c0 + j] = -v[i]
        r0 += d; c0 += len(G)
    for i in range(R): Mat[i][c0 + i] = -4
    PM = pari.matrix(R, ncols, [Mat[i][j] for i in range(R) for j in range(ncols)])
    K = pari.matkerint(PM)
    nk = int(pari.matsize(K)[1])
    lat = pari.matconcat([pari.matrix(g, nk, [K[i, j] for i in range(g) for j in range(nk)]), 4 * pari.matid(g)])
    H = pari.mathnf(lat)
    size = 4**g // int(pari.matdet(H))
    snf = [int(z) for z in pari.matsnf(H)]
    assert len(snf) == g, snf
    return size, sorted(4 // z for z in snf if z != 4)
