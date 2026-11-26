#!/usr/bin/env python3
"""The Cassels-Tate pairing on the 2-Selmer group over K = Q(i) of the four corner K-quartics of (627,364).  EXPERIMENT.

    python3 kctp.py              the four curves L1L3M1M2, L1L3M3M4, L2L4M1M2, L2L4M3M4 of (627,364) (labels of kclose.py),
                                 their Sel^2 over K, the pairing on it, the torsion over K, the rational t from torsion
    python3 kctp.py --controls   curves y^2 = x(x + a)(x + b) over Q, base-changed to K, against PARI (guide)
    python3 kctp.py --negative   forged inputs that must be refused
    (add --write to rewrite kctp.txt / kctp_controls.txt / kctp_negative.txt and kctp_hints.json; without it the
    script only prints.  --nopari: use only the hints file, never PARI.)

Written by an AI model (Claude); not read by a person.  Experiment, not counted.

THE CURVES.  With the ten linear factors of kclose.py in s = -i t and J one of the four sets above, Q_J is the product of
the four factors, r the root of the fourth, and E_J: y^2 = (x - e1)(x - e2)(x - e3) over K with e_k = den^2 A/(r_k - r),
the model of kclose.py (s = r + den^2 A / x).  For all four J the values of Q_J at both degenerate solutions s = -i, i
are squares in K (checked exactly), so the class of kclose.py is 1 and a solution of the ratio gives a point of E_J(K).
(kclose.py takes the representative -1 for two of them; -1 is a square in K, and its model differs by x -> -x.)

THE SELMER GROUP.  Sel^2(E/K) by the complete 2-descent over K, as in kdescent.py (whose local images and local class
labels are imported): pairs (b1, b2) of S-units modulo squares, the classes of (x - e1, x - e2), locally in the image of
E(K_v)/2E(K_v) at every place above 2 N((e1 - e2)(e1 - e3)(e2 - e3)).  The local images are found from local points until
they have their known order (4 at an odd place, 16 at 1 + i), so the group is exact.  The prime factorisations of the norms
are hints from PARI, checked here (product, and primality by a deterministic Miller-Rabin test valid below 3.3e24).

THE PAIRING.  For a = (a1, a2) in Sel^2, a3 = a1 a2 / g^2 (g = gcd(a1, a2)), the 2-covering
    D_a:  x = e1 + a1 u^2 = e2 + a2 v^2 = e3 + a3 w^2
and its three conics (forget one coordinate)
    C_u: a2 v^2 - a3 w^2 + (e2 - e3) z^2 = 0,  C_v: a1 u^2 - a3 w^2 + (e1 - e3) z^2 = 0,  C_w: a1 u^2 - a2 v^2 + (e1 - e2) z^2 = 0.
Each has a point over K (it is the image of D_a, which has points everywhere locally); a point is found by a lattice
search over Z[i] (LLL and enumeration in PARI as a search tool) and checked exactly.  l_u, l_v, l_w are the tangent lines
there, made primitive over Z[i].  On D_a, div(l_v / l_w) = 2 (F_v - F_w) with F_v - F_w in the class of T1, and
div(l_u / l_w) = 2 (F_u - F_w) with class T2 (the translations (u, v, w) -> (u, -v, -w), (-u, v, -w)).  For b = (b1, b2)
in Sel^2 the Cassels-Tate pairing is (Cassels 1998; the Brauer-Manin description; the same formula as the full-2-torsion
part of scripts/g1_full/rank_second.py over Q, which was checked there against PARI)
    <b, D_a> = sum over the places v of K of  [ (b2, l_v l_w (P_v))_v + (b1, l_u l_w (P_v))_v ]      (additively, in F_2)
with P_v any point of D_a(K_v) where the lines do not vanish.  The complex place contributes nothing.  The local term is
checked to be constant on D_a(K_v) at several points with all eight sign choices (a check, not a proof).
  The places that can contribute: those above 2 and the primes dividing a1 a2 a3 b1 b2 (e1 - e2)(e1 - e3)(e2 - e3), i.e. the
places of S.  At any other place v (odd, all of these are units): D_a is a smooth intersection of two quadrics over O_v
(the four singular members of its pencil are distinct modulo v), l_v and l_w are primitive, so neither vanishes
identically on the special fibre, and div(l_v / l_w) on this regular model is 2 times a Cartier divisor (no vertical
part).  Locally g = l_v / l_w = unit * phi^2, so v(g(P)) is even for every P in D_a(K_v), and (b, g(P))_v = 1 for a unit b.
This argument is written here, not checked by a program.
  Hilbert symbols.  Odd places: the tame formula; with the labels (valuation mod 2, unit not a square) it is
x0 y1 + x1 y0, because -1 is a square in every residue field of K.  The place 1 + i: the bilinear form on the 4-dimensional
K_v*/K_v*^2 is computed from the product formula (sum over all places is 0, the complex place is 0) on global representatives
of a basis, and checked against PARI's nfhilbert (guide) and the product formula on random pairs.
  Local points: x = e_k + a_k (T/D)^2 (T, D in Z[i], D a power of the uniformiser), so that x - e_k is a_k times a square;
the other two square roots by Hensel lifting in Z_p (split p), Z_p[i] (inert p) or Z_2[i], with the precision asserted.

WHAT IT GIVES.  rank E(K) <= dim Sel^2 - 2 - rank of the pairing (the kernel of the pairing contains the image of
E(K)/2E(K), whose torsion part has dimension 2).  The theory used (Cassels-Tate pairing, its Brauer-Manin description,
constancy of the local term for Selmer classes, the purity argument above) is quoted, not proved here.

TORSION.  E(K)[2^inf] by halving (P in 2E(K) iff x(P) - e_i are squares in K; the halves from x(P) + r1 r2 + r1 r3 + r2 r3,
each checked by doubling), the odd part from #E(F_p) at split primes of good reduction.  Every torsion point gives
s = r + den^2 A / x (s = r for O, s = infinity for x = 0), t = i s; a rational t is tested in the five-form system
(second_descent.classify_t).
"""
import os, sys, json, random, itertools
from fractions import Fraction as Fr
from math import gcd, isqrt

WRITE = '--write' in sys.argv
NOPARI = '--nopari' in sys.argv
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import kdescent as kd            # local images and class labels over K (existing, AI-written)
from second_descent import classify_t
gd = kd.gd

HINTS_FILE = os.path.join(here, 'kctp_hints.json')
HINTS = json.load(open(HINTS_FILE)) if os.path.exists(HINTS_FILE) else {'factor': {}, 'conic': {}}
OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(s)

_PARI = None
def PARI():
    global _PARI
    assert not NOPARI, 'PARI needed but --nopari given (hint missing)'
    if _PARI is None:
        import cypari2
        _PARI = cypari2.Pari(); _PARI.allocatemem(2 * 10 ** 9)
    return _PARI

# ------------------------------------------------------------------ integers
def is_prime(n):
    """Deterministic Miller-Rabin; valid for n < 3.3e24 (bases the primes up to 41)."""
    assert n < 3317044064679887385961981, 'too large for the deterministic test'
    if n < 2: return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
        if n % p == 0: return n == p
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

def factor_int(N):
    """Prime factors of |N| (a hint from PARI, checked: the product, and each factor prime)."""
    N = abs(N); assert N > 0
    key = str(N)
    if key in HINTS['factor']: ps = [int(p) for p in HINTS['factor'][key]]
    else:
        ps = sorted(set(int(p) for p in PARI().factor(N)[0])); HINTS['factor'][key] = [str(p) for p in ps]
    M = N
    for p in ps:
        assert is_prime(p), ('not a prime', p)
        while M % p == 0: M //= p
    assert M == 1, 'the factorisation hint is incomplete'
    return ps

def vp(x, p):
    assert x != 0
    k = 0
    while x % p == 0: x //= p; k += 1
    return k

def legendre(a, p):
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1

def sqrt_mod_p(a, p):
    """Tonelli-Shanks; a a non-zero square modulo the odd prime p."""
    a %= p; assert legendre(a, p) == 1
    if p % 4 == 3: return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while legendre(z, p) != -1: z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1: t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    assert r * r % p == a
    return r

def sqrt_fp2(a, b, p):
    """(x, y) with (x + i y)^2 = a + i b in F_p[i] = F_(p^2) (p = 3 mod 4), or None if a + i b is not a square."""
    a %= p; b %= p
    if (a, b) == (0, 0): return (0, 0)
    Nu = (a * a + b * b) % p
    if legendre(Nu, p) != 1: return None
    nn = sqrt_mod_p(Nu, p)
    for ns in (nn, (p - nn) % p):
        h = (a + ns) * pow(2, -1, p) % p
        if h == 0:
            if b == 0 and legendre(-a, p) == 1: return (0, sqrt_mod_p(-a, p))
            continue
        if legendre(h, p) == 1:
            x0 = sqrt_mod_p(h, p); y0 = b * pow(2 * x0, -1, p) % p
            assert ((x0 * x0 - y0 * y0 - a) % p, (2 * x0 * y0 - b) % p) == (0, 0)
            return (x0, y0)
    raise AssertionError('sqrt in F_p^2 failed')

# ------------------------------------------------------------------ Gaussian integers (pairs) and Gaussian rationals
def gm(x, y): return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])
def ga(x, y): return (x[0] + y[0], x[1] + y[1])
def gs(x, y): return (x[0] - y[0], x[1] - y[1])
def gc(x): return (x[0], -x[1])
def gn(x): return x[0] * x[0] + x[1] * x[1]
def gsc(c, x): return (c * x[0], c * x[1])
def gpow(x, k):
    r = (1, 0)
    for _ in range(k): r = gm(r, x)
    return r
def gprod(xs):
    r = (1, 0)
    for x in xs: r = gm(r, x)
    return r
def gdivmod(x, y):
    """Euclidean division in Z[i]: x = q y + r with N(r) <= N(y)/2."""
    N = gn(y); t = gm(x, gc(y))
    q = (((2 * t[0] + N) // (2 * N)), ((2 * t[1] + N) // (2 * N)))
    return q, gs(x, gm(q, y))
def gdiv(x, y):
    q, r = gdivmod(x, y); assert r == (0, 0), 'inexact Gaussian division'; return q
def gdivides(y, x): return gdivmod(x, y)[1] == (0, 0)
def ggcd(x, y):
    while y != (0, 0): x, y = y, gdivmod(x, y)[1]
    return x
def gxgcd(x, y):
    """(g, s, t) with s x + t y = g = gcd(x, y)."""
    r0, r1, s0, s1, t0, t1 = x, y, (1, 0), (0, 0), (0, 0), (1, 0)
    while r1 != (0, 0):
        q, r = gdivmod(r0, r1)
        r0, r1 = r1, r; s0, s1 = s1, gs(s0, gm(q, s1)); t0, t1 = t1, gs(t0, gm(q, t1))
    return r0, s0, t0
def gmod(x, y): return gdivmod(x, y)[1]
def ginv_mod(x, y):
    g, s, t = gxgcd(x, y); assert gn(g) == 1, 'not invertible'
    return gmod(gm(s, gc(g)), y)
def gsqrt_exact(z):
    """sqrt of a Gaussian integer in Z[i] if it is a square there, else None."""
    if z == (0, 0): return (0, 0)
    N = gn(z); n = isqrt(N)
    if n * n != N: return None
    a2, b2 = n + z[0], n - z[0]
    if a2 % 2 or b2 % 2: return None
    a2 //= 2; b2 //= 2; a, b = isqrt(a2), isqrt(b2)
    if a * a != a2 or b * b != b2: return None
    for w in ((a, b), (a, -b)):
        if gm(w, w) == z: return w
    return None

def q_(x): return (Fr(x[0]), Fr(x[1])) if isinstance(x, tuple) else (Fr(x), Fr(0))
def qm(x, y): return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])
def qa(x, y): return (x[0] + y[0], x[1] + y[1])
def qs(x, y): return (x[0] - y[0], x[1] - y[1])
def qinv(x):
    N = x[0] * x[0] + x[1] * x[1]; assert N != 0; return (x[0] / N, -x[1] / N)
def qdiv(x, y): return qm(x, qinv(y))
def q_is_square(x):
    """x a Gaussian rational: is it a square in K?  (exact)  Returns a square root or None."""
    if x == (0, 0): return (Fr(0), Fr(0))
    d = 1
    for c in x: d = d * c.denominator // gcd(d, c.denominator)
    z = (int(x[0] * d * d), int(x[1] * d * d)); w = gsqrt_exact(z)
    return None if w is None else (Fr(w[0], d), Fr(w[1], d))

# ------------------------------------------------------------------ the curves
def factors(m, n):
    """The ten linear factors of kclose.py as (c1, c0): c1 s + c0."""
    p, pb = (m, n), (m, -n)
    return {'L1': ((-n, 0), (m, 0)), 'L2': ((n, 0), (m, 0)), 'L3': ((m, 0), (n, 0)), 'L4': ((-m, 0), (n, 0)),
            'L5': ((-1, 0), (1, 0)), 'L6': ((1, 0), (1, 0)),
            'M1': (p, gsc(-1, pb)), 'M2': (pb, p), 'M3': (p, pb), 'M4': (pb, gsc(-1, p))}

def pmul(P, Q):
    R = [(Fr(0), Fr(0))] * (len(P) + len(Q) - 1)
    for i, a in enumerate(P):
        for j, b in enumerate(Q): R[i + j] = qa(R[i + j], qm(a, b))
    return R
def peval(P, x):
    r = (Fr(0), Fr(0))
    for c in reversed(P): r = qa(qm(r, x), c)
    return r
def ptaylor(P, r):
    """coefficients of P(r + T) in T."""
    out = []; Q = list(P)
    for k in range(len(P)):
        # k-th Taylor coefficient: synthetic division
        val = (Fr(0), Fr(0)); quo = []
        for c in reversed(Q): val = qa(qm(val, r), c); quo.append(val)
        out.append(quo[-1]); Q = list(reversed(quo[:-1]))
    return out

def build_curve(m, n, J):
    F = factors(m, n); names = [J[i:i + 2] for i in range(0, 8, 2)]
    Q = [(Fr(1), Fr(0))]
    for k in names: c1, c0 = F[k]; Q = pmul(Q, [q_(c0), q_(c1)])
    roots = {k: qdiv(q_(gsc(-1, F[k][1])), q_(F[k][0])) for k in names}
    for k in names: assert peval(Q, roots[k]) == (0, 0)
    S = m * m + n * n; assert isqrt(S) ** 2 == S, 'the square case only (degenerate solutions t = 1, -1, s = -i, i)'
    degen = {}
    for t in (1, -1):
        v = peval(Q, (Fr(0), Fr(-t)))
        assert v != (0, 0) and q_is_square(v) is not None, 'the value at a degenerate solution is not a square in K'
        degen[t] = v
    r = roots[names[3]]
    c = ptaylor(Q, r); assert c[0] == (0, 0)
    A, B, C, E0 = c[1], c[2], c[3], c[4]
    xs = [qdiv(A, qs(roots[k], r)) for k in names[:3]]
    den = 1
    for x in xs:
        for y in x: den = den * y.denominator // gcd(den, y.denominator)
    e = [(int(x[0] * den * den), int(x[1] * den * den)) for x in xs]
    for x, ee in zip(xs, e): assert (x[0] * den * den, x[1] * den * den) == q_(ee)
    # the cubic: x^3 + B x^2 + A C x + A^2 E0, scaled by x -> x / den^2, equals (x - e1)(x - e2)(x - e3)
    D2 = Fr(den * den)
    lhs = [qm(qm(A, A), E0), qm(A, C), B, (Fr(1), Fr(0))]
    lhs = [(c_[0] * D2 ** (3 - k), c_[1] * D2 ** (3 - k)) for k, c_ in enumerate(lhs)]
    rhs = [(Fr(1), Fr(0))]
    for ee in e: rhs = pmul(rhs, [q_(gsc(-1, ee)), (Fr(1), Fr(0))])
    assert lhs == rhs, 'the Weierstrass model is not (x - e1)(x - e2)(x - e3)'
    assert len(set(e)) == 3
    return dict(J=J, Q=Q, r=r, A=A, den=den, e=e, degen=degen)

# ------------------------------------------------------------------ S-units over Z[i]
class SUnits:
    """The S-units of Z[i] modulo squares for S = {1 + i} and the primes of a list: basis i, 1 + i, the Gaussian primes."""
    def __init__(s, primes):
        s.primes = sorted(set(primes) | {2})
        s.basis = [(0, 1), (1, 1)]; s.over = [None, 2]
        for p in s.primes:
            if p == 2: continue
            if p % 4 == 3: s.basis.append((p, 0)); s.over.append(p)
            else:
                g = gd.gprime(p); s.basis += [g, gc(g)]; s.over += [p, p]
        s.g = len(s.basis)
    def vec(s, z, allow_square_cofactor=False):
        """exponents mod 2 of the Gaussian integer z != 0 on the basis (bit k), and the check that z is an S-unit
        (or an S-unit times a square, if allowed)."""
        assert z != (0, 0)
        v = 0
        for k in range(1, s.g):
            pi = s.basis[k]; e = 0
            while gdivides(pi, z): z = gdiv(z, pi); e += 1
            if e % 2: v |= 1 << k
        if gn(z) != 1:
            assert allow_square_cofactor, 'not an S-unit'
            u = None
            for k_ in range(4):
                if gsqrt_exact(gm(z, gpow((0, 1), k_))) is not None: u = k_; break
            assert u is not None, 'not an S-unit times a square'
            z = gpow((0, -1), u)
        k_ = [(1, 0), (0, 1), (-1, 0), (0, -1)].index(z)
        if k_ % 2: v |= 1
        return v
    def elt(s, v):
        return gprod([s.basis[k] for k in range(s.g) if v >> k & 1])
    def qvec(s, x):
        """vec of a non-zero Gaussian rational (modulo squares), allowing square cofactors."""
        d = 1
        for c in x: d = d * c.denominator // gcd(d, c.denominator)
        return s.vec((int(x[0] * d * d), int(x[1] * d * d)), allow_square_cofactor=True)

def f2_span(vecs):
    B = []
    for v in vecs:
        for b in B: v = min(v, v ^ b)
        if v: B.append(v)
    return B
def f2_rank(vecs): return len(f2_span(vecs))

def selmer2(e, seed=1):
    """Sel^2(E/K) for y^2 = prod (x - e_i): list of pairs (v1, v2) of basis vectors, and the S-unit data."""
    rng = random.Random(seed)
    ps = set([2])
    for i in range(3):
        for j in range(i + 1, 3): ps |= set(factor_int(gn(gs(e[i], e[j]))))
    SU = SUnits(ps); g = SU.g
    cons = []; places = []
    for p in SU.primes:
        for pl in kd.places_above(p):
            img = kd.local_image(e, pl, rng)
            assert img is not None, ('local image not found', p)
            places.append(pl)
            d = pl.dim; M = [pl.vec(b) for b in SU.basis]; assert None not in M
            for phi in kd.nullspace(img, 2 * d):
                row = 0
                for k in range(g):
                    a = sum((phi >> t & 1) & M[k][t] for t in range(d)) & 1
                    c = sum((phi >> (d + t) & 1) & M[k][t] for t in range(d)) & 1
                    row |= a << k | c << (g + k)
                cons.append(row)
    sel = kd.nullspace(cons, 2 * g)
    assert len(sel) == 2 * g - len(kd.rank_basis(cons))
    return SU, [(v & ((1 << g) - 1), v >> g) for v in sel], cons

def in_span(v, B):
    for b in sorted(f2_span(B), reverse=True): v = min(v, v ^ b)
    return v == 0

# ------------------------------------------------------------------ local fields at the places of K
UNIT2 = kd.UNIT2
class Loc:
    """A finite place of K: 'split' (Z_p, i -> r, the embedding of kdescent's place idx), 'inert' (Z_p[i]) or 'ram' (Z_2[i]);
    elements modulo p^N (pairs for inert and ram)."""
    def __init__(s, p, idx, N=None):
        s.p, s.idx = p, idx
        s.kdp = kd.places_above(p)[idx]
        if p == 2: s.kind = 'ram'; s.N = N or 400
        elif p % 4 == 1:
            s.kind = 'split'; s.N = N or max(60, 1200 // p.bit_length())
            L = s.kdp.L; r0 = L.r if idx == 0 else L.mod - L.r
            mod = p ** s.N; r = r0 % p; mm = p
            while mm < mod:
                mm = min(mm * mm, mod); r = (r - (r * r + 1) * pow(2 * r, -1, mm)) % mm
            assert (r * r + 1) % mod == 0 and (r - r0) % p ** min(s.N, L.K) == 0
            s.r = r
        else: s.kind = 'inert'; s.N = N or max(60, 1200 // p.bit_length())
        s.mod = s.p ** s.N
        s.dim = 4 if p == 2 else 2
    def name(s): return '%d%s' % (s.p, '' if s.kind != 'split' else 'ab'[s.idx])
    def emb(s, z):
        if s.kind == 'split': return (z[0] + z[1] * s.r) % s.mod
        return (z[0] % s.mod, z[1] % s.mod)
    def mul(s, x, y):
        if s.kind == 'split': return x * y % s.mod
        return ((x[0] * y[0] - x[1] * y[1]) % s.mod, (x[0] * y[1] + x[1] * y[0]) % s.mod)
    def add(s, x, y):
        if s.kind == 'split': return (x + y) % s.mod
        return ((x[0] + y[0]) % s.mod, (x[1] + y[1]) % s.mod)
    def neg(s, x):
        if s.kind == 'split': return (-x) % s.mod
        return ((-x[0]) % s.mod, (-x[1]) % s.mod)
    def iszero(s, x): return x == 0 if s.kind == 'split' else x == (0, 0)
    def val(s, x):
        """valuation (normalised: v(uniformiser) = 1; uniformiser p, p, 1 + i); None if x = 0 modulo p^N."""
        if s.iszero(x): return None
        if s.kind == 'split': return vp(x, s.p)
        if s.kind == 'inert': return min(vp(x[0], s.p) if x[0] else 10 ** 9, vp(x[1], s.p) if x[1] else 10 ** 9)
        return vp(x[0] * x[0] + x[1] * x[1], 2)
    def prec(s):
        """the valuation of p^N."""
        return 2 * s.N if s.kind == 'ram' else s.N
    def vec(s, x):
        """class in K_v*/K_v*^2 with the labels of kdescent.Place.vec; asserts that the precision determines it."""
        v = s.val(x); assert v is not None, 'p-adic precision (zero)'
        if s.kind == 'ram':
            assert v + 6 <= s.prec(), 'p-adic precision'
            c = gd.Local(2).cls(x, s.N)[0]; assert c is not None
            return (c[0],) + UNIT2[c[1]]
        assert v + 2 <= s.prec(), 'p-adic precision'
        if s.kind == 'split':
            u = (x // s.p ** v) % s.p; return (v % 2, 0 if legendre(u, s.p) == 1 else 1)
        u = (x[0] // s.p ** v, x[1] // s.p ** v)
        return (v % 2, 0 if legendre(u[0] * u[0] + u[1] * u[1], s.p) == 1 else 1)
    def gvec(s, z):
        """class of a Gaussian integer (exact), the same labels; equals kdescent's."""
        w = s.kdp.vec(z); assert w is not None; return w
    def sqrt(s, x):
        """y with y^2 = x modulo p^(N - v(x)/2 - 3) (precision tracked by the caller via prec_loss), or None."""
        v = s.val(x); assert v is not None and v + 8 <= s.prec(), 'p-adic precision (sqrt)'
        if v % 2: return None
        p, N = s.p, s.N
        if s.kind == 'split':
            u = (x // p ** v) % s.mod
            if legendre(u, p) != 1: return None
            y = sqrt_mod_p(u, p); mm = p
            while mm < s.mod:
                mm = min(mm * mm, s.mod); y = (y - (y * y - u) * pow(2 * y, -1, mm)) % mm
            return p ** (v // 2) * y % s.mod
        if s.kind == 'inert':
            u = (x[0] // p ** v, x[1] // p ** v)
            y0 = sqrt_fp2(u[0], u[1], p)
            if y0 is None: return None
            y = y0; mm = p
            while mm < s.mod:
                mm = min(mm * mm, s.mod)
                t = ((y[0] * y[0] - y[1] * y[1] - u[0]) % mm, (2 * y[0] * y[1] - u[1]) % mm)
                d = (2 * y[0], 2 * y[1]); Nd = (d[0] * d[0] + d[1] * d[1]) % mm; iN = pow(Nd, -1, mm)
                q = (((t[0] * d[0] + t[1] * d[1]) * iN) % mm, ((t[1] * d[0] - t[0] * d[1]) * iN) % mm)   # t / d
                y = ((y[0] - q[0]) % mm, (y[1] - q[1]) % mm)
            sc = p ** (v // 2)
            return ((y[0] * sc) % s.mod, (y[1] * sc) % s.mod)
        # ram: x = (2i)^(v/2) u
        j = v // 2; a, b = x
        assert a % (2 ** j) == 0 and b % (2 ** j) == 0
        u = (a >> j, b >> j)
        for _ in range(j % 4): u = (u[1], -u[0])                      # times (-i)^j
        M = 2 ** N
        y = None
        for c0 in range(8):
            for c1 in range(8):
                if (c0 + c1) % 2 == 0: continue
                dd = (c0 * c0 - c1 * c1 - u[0], 2 * c0 * c1 - u[1])
                if vp(dd[0] ** 2 + dd[1] ** 2, 2) >= 5 if dd != (0, 0) else True: y = (c0, c1); break
            if y: break
        if y is None: return None
        k = 5
        while k < 2 * N + 4:
            t = (y[0] * y[0] - y[1] * y[1] - u[0], 2 * y[0] * y[1] - u[1])
            if t == (0, 0): break
            assert t[0] % 2 == 0 and t[1] % 2 == 0
            t = (t[0] // 2, t[1] // 2)
            Ny = y[0] * y[0] + y[1] * y[1]; iN = pow(Ny, -1, 4 * M)
            q = (((t[0] * y[0] + t[1] * y[1]) * iN) % (4 * M), ((t[1] * y[0] - t[0] * y[1]) * iN) % (4 * M))
            y = ((y[0] - q[0]) % (4 * M), (y[1] - q[1]) % (4 * M))
            k = 2 * k - 4
        sc = (1, 0)
        for _ in range(j): sc = gm(sc, (1, 1))
        y = gm(y, sc)
        return (y[0] % M, y[1] % M)

def check_sqrt(L, x, y, loss):
    """y^2 = x modulo uniformiser^(prec - loss)."""
    d = L.add(L.mul(y, y), L.neg(x))
    v = L.val(d)
    assert v is None or v >= L.prec() - loss, ('square root not precise', L.name(), v)

# ------------------------------------------------------------------ Hilbert symbols on the label vectors
def h_odd(x, y): return (x[0] * y[1] + x[1] * y[0]) % 2

def odd_places_of(z):
    out = []
    for p in factor_int(gn(z)):
        if p == 2: continue
        out += [Loc(p, i, N=20) for i in range(2 if p % 4 == 1 else 1)]
    return out

def hilbert_global_odd(a, b):
    """sum over the odd places of (a, b)_v for Gaussian integers a, b != 0 (additive)."""
    ps = set(factor_int(gn(a))) | set(factor_int(gn(b))); tot = 0
    for p in sorted(ps):
        if p == 2: continue
        for i in range(2 if p % 4 == 1 else 1):
            L = Loc(p, i, N=20); tot += h_odd(L.gvec(a), L.gvec(b))
    return tot % 2

def ram_basis():
    """Gaussian integers z_0..z_3 whose labels at 1 + i are the unit vectors."""
    P2 = kd.places_above(2)[0]; want = {(1, 0, 0, 0): (1, 1)}
    for a in range(-6, 7):
        for b in range(-6, 7):
            if (a + b) % 2 == 0: continue
            v = P2.vec((a, b))
            if v[0] == 0 and sum(v) == 1 and v not in want: want[v] = (a, b)
    assert len(want) == 4
    return [want[tuple(1 if j == k else 0 for j in range(4))] for k in range(4)]
RB = None; H2 = None
def setup_h2():
    global RB, H2
    RB = ram_basis()
    H2 = [[hilbert_global_odd(RB[j], RB[k]) for k in range(4)] for j in range(4)]   # product formula
    for j in range(4):
        for k in range(4): assert H2[j][k] == H2[k][j], 'the 2-adic Hilbert form is not symmetric'
    # non-degenerate over F_2
    rows = [sum(H2[j][k] << k for k in range(4)) for j in range(4)]
    assert f2_rank(rows) == 4, 'the 2-adic Hilbert form is degenerate'
def h_ram(x, y): return sum(x[j] * H2[j][k] * y[k] for j in range(4) for k in range(4)) % 2
def h_loc(L, x, y): return h_ram(x, y) if L.kind == 'ram' else h_odd(x, y)

def hilbert_selftest(nrand=300, seed=5):
    """The product formula on random pairs (the 2-adic form applied to their labels; a test of bilinearity and of the
    labels) and PARI's nfhilbert at every place (guide)."""
    rng = random.Random(seed); P2 = Loc(2, 0); nbad = 0; npari = 0
    nf = None if NOPARI else PARI().nfinit('y^2+1')
    for _ in range(nrand):
        a = (rng.randint(-300, 300), rng.randint(-300, 300)); b = (rng.randint(-300, 300), rng.randint(-300, 300))
        if a == (0, 0) or b == (0, 0): continue
        tot = (h_ram(P2.gvec(a), P2.gvec(b)) + hilbert_global_odd(a, b)) % 2
        assert tot == 0, ('product formula fails', a, b)
        if nf is not None and npari < 60:
            pr = PARI()
            A = pr('%d + %d*y' % a); Bb = pr('%d + %d*y' % b)
            for p in sorted(set(factor_int(gn(a))) | set(factor_int(gn(b))) | {2}):
                prs = pr.idealprimedec(nf, p)
                Ls = [Loc(p, i, N=30) for i in range(2 if p % 4 == 1 else 1)]
                for L in Ls:
                    mine = h_loc(L, L.gvec(a), L.gvec(b))
                    # match the PARI prime to the place: the one containing (a - r) if split, by valuation of pi
                    if L.kind == 'split':
                        g = gd.gprime(p)
                        cand = [P for P in prs if int(pr.idealval(nf, pr('%d + %d*y' % (g[0], g[1])), P)) > 0]
                        # place idx: kdescent's r; the prime above p containing y - r
                        P = [P for P in prs if int(pr.idealval(nf, pr('y - %d' % L.r), P)) > 0][0]
                    else: P = prs[0]
                    h = int(pr.nfhilbert(nf, A, Bb, P))
                    if (h == -1) != (mine == 1): nbad += 1
            npari += 1
    assert nbad == 0, ('disagreement with PARI nfhilbert', nbad)
    return nrand, npari

# ------------------------------------------------------------------ conics over K
def qconic(co, P): return ga(ga(gm(co[0], gm(P[0], P[0])), gm(co[1], gm(P[1], P[1]))), gm(co[2], gm(P[2], P[2])))

def primitive(P):
    """a Gaussian-rational vector scaled to a primitive vector of Z[i]^3."""
    d = 1
    for c in P:
        for y in c: d = d * y.denominator // gcd(d, y.denominator)
    Q = [(int(c[0] * d), int(c[1] * d)) for c in P]
    g = (0, 0)
    for c in Q: g = ggcd(g, c) if g != (0, 0) else c
    return tuple(gdiv(c, g) for c in Q)

def reduce_conic(co, SU):
    """(red, g): red square-free, pairwise coprime coefficients; a point P' of red gives the point (g_i P'_i) of co."""
    c = list(co); g = [(Fr(1), Fr(0))] * 3
    for k in range(1, SU.g):
        pi = SU.basis[k]; piq = q_(pi)
        while True:
            e = [0, 0, 0]
            for i in range(3):
                z = c[i]
                while gdivides(pi, z): z = gdiv(z, pi); e[i] += 1
            if max(e) >= 2:
                i = e.index(max(e)); c[i] = gdiv(c[i], gm(pi, pi)); g[i] = qdiv(g[i], piq); continue
            if sum(e) == 3:
                c = [gdiv(x, pi) for x in c]; continue
            if sum(e) == 2:
                i0, i1 = [i for i in range(3) if e[i]]; j = 3 - i0 - i1
                c[i0] = gdiv(c[i0], pi); c[i1] = gdiv(c[i1], pi); c[j] = gm(c[j], pi)
                g[i0] = qdiv(g[i0], piq); g[i1] = qdiv(g[i1], piq); continue
            break
    for i in range(3): SU.vec(c[i])                     # still S-units
    return c, g

def crt_sqrt(t, ks, SU, rng):
    """x with x^2 = t modulo the product M of the Gaussian primes basis[k], k in ks (random signs); (x, M)."""
    x = (0, 0); M = (1, 0)
    for k in ks:
        pi = SU.basis[k]; p = SU.over[k]
        if p == 2: r_ = (1, 0)
        elif p % 4 == 3:
            r_ = sqrt_fp2(t[0], t[1], p); assert r_ is not None, ('not a square modulo', p)
        else:
            ii = (-pi[0] * pow(pi[1], -1, p)) % p           # i modulo pi
            tv = (t[0] + t[1] * ii) % p; assert legendre(tv, p) == 1, ('not a square modulo', pi)
            r_ = (sqrt_mod_p(tv, p), 0)
        if rng.random() < 0.5: r_ = gsc(-1, r_)
        g_, s_, t_ = gxgcd(M, pi); assert gn(g_) == 1
        gi = gc(g_)
        x = gmod(ga(gm(gm(x, gm(t_, gi)), pi), gm(gm(r_, gm(s_, gi)), M)), gm(M, pi)); M = gm(M, pi)
        assert gmod(gs(gm(x, x), t), pi) == (0, 0)
    return x, M

def solve_conic(co, SU, rng):
    """A primitive point of co[0] X^2 + co[1] Y^2 + co[2] Z^2 = 0 over Z[i] (co S-units times squares); hints cached.
    Search only (PARI's qflllgram and qfminim on the lattice where the reduced form is divisible by abc); the point is
    checked exactly on the original conic."""
    key = str([list(c) for c in co])
    if key in HINTS['conic']:
        P = tuple(tuple(int(t) for t in c) for c in HINTS['conic'][key])
        assert P != ((0, 0),) * 3 and qconic(co, P) == (0, 0), 'a conic hint is not a point of the conic'
        return P
    (a, b, c), g = reduce_conic(co, SU)
    ks = [[k for k in range(1, SU.g) if gdivides(SU.basis[k], z)] for z in (a, b, c)]
    pr = PARI()
    def mulmat(w): return [[w[0], -w[1]], [w[1], w[0]]]
    for attempt in range(16):
        la, Ma = crt_sqrt(gsc(-1, gm(b, c)), ks[0], SU, rng)
        lb, Mb = crt_sqrt(gsc(-1, gm(a, c)), ks[1], SU, rng)
        lc, Mc = crt_sqrt(gsc(-1, gm(a, b)), ks[2], SU, rng)
        # b Y = la Z (a), c Z = lb X (b), a X = lc Y (c)  <=>  Y = mu1 Z (Ma), Z = mu2 X (Mb), X = mu3 Y (Mc)
        mu1 = gmod(gm(la, ginv_mod(b, Ma)), Ma) if ks[0] else (0, 0)
        mu2 = gmod(gm(lb, ginv_mod(c, Mb)), Mb) if ks[1] else (0, 0)
        mu3 = gmod(gm(lc, ginv_mod(a, Mc)), Mc) if ks[2] else (0, 0)
        F = [[0] * 12 for _ in range(6)]
        def put(r0, c0, mat, sg=1):
            for i_ in range(2):
                for j_ in range(2): F[r0 + i_][c0 + j_] += sg * mat[i_][j_]
        I2 = [[1, 0], [0, 1]]
        put(0, 2, I2); put(0, 4, mulmat(mu1), -1); put(2, 4, I2); put(2, 0, mulmat(mu2), -1)
        put(4, 0, I2); put(4, 2, mulmat(mu3), -1)
        for blk, Mm in enumerate((Ma, Mb, Mc)): put(2 * blk, 6 + 2 * blk, mulmat(Mm))
        ker = pr.matkerint(pr.matrix(6, 12, [x for row in F for x in row]))
        ncol = int(pr.matsize(ker)[1]); assert ncol == 6
        Bv = [[int(pr.component(pr.component(ker, j_ + 1), i_ + 1)) for i_ in range(6)] for j_ in range(6)]
        for v_ in Bv:
            X, Y, Z = (v_[0], v_[1]), (v_[2], v_[3]), (v_[4], v_[5])
            assert gmod(gs(Y, gm(mu1, Z)), Ma) == (0, 0) and gmod(gs(Z, gm(mu2, X)), Mb) == (0, 0) and gmod(gs(X, gm(mu3, Y)), Mc) == (0, 0)
            assert gmod(qconic((a, b, c), (X, Y, Z)), gm(gm(a, b), c)) == (0, 0), 'the form is not divisible by abc on the lattice'
        wts = [isqrt(gn(a)) + 1] * 2 + [isqrt(gn(b)) + 1] * 2 + [isqrt(gn(c)) + 1] * 2
        G = [[sum(wts[t] * Bv[i_][t] * Bv[j_][t] for t in range(6)) for j_ in range(6)] for i_ in range(6)]
        T = pr.qflllgram(pr.matrix(6, 6, [x for row in G for x in row]))
        Tm = [[int(pr.component(pr.component(T, j_ + 1), i_ + 1)) for i_ in range(6)] for j_ in range(6)]  # Tm[j] = column j
        R = [[sum(Tm[j_][k] * Bv[k][t] for k in range(6)) for t in range(6)] for j_ in range(6)]
        Gr = [[sum(wts[t] * R[i_][t] * R[j_][t] for t in range(6)) for j_ in range(6)] for i_ in range(6)]
        # enumeration of small combinations of the reduced basis (a search only)
        order = sorted(range(6), key=lambda j_: Gr[j_][j_])
        R = [R[j_] for j_ in order]
        cands = sorted(itertools.product(range(-2, 3), repeat=6), key=lambda cv: sum(abs(x) for x in cv))
        for cv in cands:
            if True:
                v_ = [sum(cv[k] * R[k][t] for k in range(6)) for t in range(6)]
                Pr = ((v_[0], v_[1]), (v_[2], v_[3]), (v_[4], v_[5]))
                if Pr != ((0, 0),) * 3 and qconic((a, b, c), Pr) == (0, 0):
                    P = primitive([qm(g[i], q_(Pr[i])) for i in range(3)])
                    assert qconic(co, P) == (0, 0)
                    HINTS['conic'][key] = [list(t) for t in P]
                    return P
    raise AssertionError('no point found on a conic (search failed)')

def tangent(co, P):
    """the tangent line at P of the diagonal conic, primitive over Z[i]: coefficients (co_i P_i) / gcd."""
    l = [gm(co[i], P[i]) for i in range(3)]
    g = (0, 0)
    for c in l: g = ggcd(g, c) if g != (0, 0) else c
    assert g != (0, 0), 'singular point'
    return tuple(gdiv(c, g) for c in l)

# ------------------------------------------------------------------ the 2-coverings D_a and the pairing
SIGNS = list(itertools.product((1, -1), repeat=3))
NPTS = 6

class Cover:
    def __init__(s, e, a1, a2, SU, rng, forced_points=None, forced_lines=None):
        g = ggcd(a1, a2); a3 = gdiv(gm(a1, a2), gm(g, g))
        s.e, s.a, s.SU = e, (a1, a2, a3), SU
        e1, e2, e3 = e
        s.co = {'u': (a2, gsc(-1, a3), gs(e2, e3)), 'v': (a1, gsc(-1, a3), gs(e1, e3)), 'w': (a1, gsc(-1, a2), gs(e1, e2))}
        s.P = {}; s.l = {}
        for k in 'uvw':
            P = forced_points[k] if forced_points and k in forced_points else solve_conic(s.co[k], SU, rng)
            assert qconic(s.co[k], P) == (0, 0) and P != ((0, 0),) * 3, 'not a point of the conic C_%s' % k
            s.P[k] = P
            s.l[k] = forced_lines[k] if forced_lines and k in forced_lines else tangent(s.co[k], P)
        s.cache = {}
    def local_points(s, L, rng, want=None):
        """points (u', v', w', z') of D_a(K_v) as elements of the local ring, with the precision they are known to."""
        if want is None: want = NPTS
        if L.name() in s.cache: return s.cache[L.name()]
        e = s.e; a = s.a; out = []
        unif = (1, 1) if L.kind == 'ram' else (L.p, 0)
        a123 = gprod(a)
        for it in range(20000):
            if len(out) >= want: break
            k = rng.randrange(3)
            if it % 2 == 0:              # x = e_k + a_k (T/D)^2
                D = gpow(unif, rng.randrange(0, 5)); T = gm(gpow(unif, rng.randrange(0, 7)), (rng.randint(-999, 999), rng.randint(-999, 999)))
                if T == (0, 0): continue
                X = ga(gm(e[k], gm(D, D)), gm(a[k], gm(T, T)))
            else:                        # x = X / D^2 as in kdescent.local_image: near a root, random, or of negative valuation
                mode = rng.randrange(3); D = (1, 0); T = (0, 0)
                if mode == 0: X = ga(e[k], gm(gpow(unif, rng.randrange(0, 12)), (rng.randint(-99, 99), rng.randint(-99, 99))))
                elif mode == 1: X = (rng.randint(-10 ** 30, 10 ** 30), rng.randint(-10 ** 30, 10 ** 30))
                else:
                    D = gpow(unif, rng.randrange(1, 5)); X = (rng.randint(-999, 999), rng.randint(-999, 999))
                k = -1
            Ns = [gm(gs(X, gm(e[i], gm(D, D))), a[i]) for i in range(3)]
            if any(N_ == (0, 0) for N_ in Ns): continue
            roots = []; vmax = 0; ok = True
            for i in range(3):
                if i == k:
                    roots.append(L.emb(gm(a[k], T)))
                    vk = L.val(L.emb(Ns[i]))
                    if vk is None or vk + 16 > L.prec(): ok = False; break
                    vmax = max(vmax, vk); continue
                x_ = L.emb(Ns[i]); v_ = L.val(x_)
                if v_ is None or v_ + 16 > L.prec(): ok = False; break
                y_ = L.sqrt(x_)
                if y_ is None: ok = False; break
                check_sqrt(L, x_, y_, v_ // 2 + 4)
                roots.append(y_); vmax = max(vmax, v_)
            if not ok: continue
            others = [gm(a[1], a[2]), gm(a[0], a[2]), gm(a[0], a[1])]
            pt = [L.mul(roots[i], L.emb(others[i])) for i in range(3)] + [L.emb(gm(D, a123))]
            out.append((pt, L.prec() - vmax // 2 - 8, (k, D, T)))
        assert len(out) == want, ('no local points found', L.name())
        s.cache[L.name()] = out
        return out
    def lval(s, L, k, pt):
        idx = {'u': (1, 2, 3), 'v': (0, 2, 3), 'w': (0, 1, 3)}[k]
        r = L.emb((0, 0))
        for c, j in zip(s.l[k], idx): r = L.add(r, L.mul(L.emb(c), pt[j]))
        return r
    def check_point(s, L, pt, pr):
        """the two quadrics of D_a vanish at pt modulo uniformiser^pr."""
        e1, e2, e3 = s.e; a1, a2, a3 = s.a
        u, v, w, z = pt
        q1 = L.add(L.add(L.mul(L.emb(a1), L.mul(u, u)), L.neg(L.mul(L.emb(a2), L.mul(v, v)))), L.mul(L.emb(gs(e1, e2)), L.mul(z, z)))
        q2 = L.add(L.add(L.mul(L.emb(a1), L.mul(u, u)), L.neg(L.mul(L.emb(a3), L.mul(w, w)))), L.mul(L.emb(gs(e1, e3)), L.mul(z, z)))
        for q in (q1, q2):
            vq = L.val(q); assert vq is None or vq >= pr, ('a local point is not on D_a', L.name())
    def local_term(s, L, b, rng):
        """the set of values of (b2, l_v l_w)_v + (b1, l_u l_w)_v over the sample points and all sign choices."""
        vb1, vb2 = L.gvec(b[0]), L.gvec(b[1]); vals = set()
        for pt0, pr, _ in s.local_points(L, rng):
            for sg in SIGNS:
                pt = [pt0[0] if sg[0] == 1 else L.neg(pt0[0]), pt0[1] if sg[1] == 1 else L.neg(pt0[1]),
                      pt0[2] if sg[2] == 1 else L.neg(pt0[2]), pt0[3]]
                s.check_point(L, pt, pr)
                ls = {}
                for k in 'uvw':
                    x_ = s.lval(L, k, pt); vx = L.val(x_)
                    if vx is None or vx + 8 > pr: ls = None; break
                    ls[k] = L.vec(x_)
                if ls is None: continue
                gv = tuple((p_ + q_) % 2 for p_, q_ in zip(ls['v'], ls['w'])); gu = tuple((p_ + q_) % 2 for p_, q_ in zip(ls['u'], ls['w']))
                vals.add((h_loc(L, vb2, gv) + h_loc(L, vb1, gu)) % 2)
        return vals
    def pairing(s, b, places, rng):
        tot = 0; per = {}
        for L in places:
            vals = s.local_term(L, b, rng)
            assert len(vals) == 1, ('the local term is not constant on D_a(K_v)', L.name(), vals)
            per[L.name()] = min(vals); tot += per[L.name()]
        return tot % 2, per

def other_point(co, P):
    """a second point of the diagonal conic: the other intersection with the line through P and a fixed R off the conic."""
    for R in ((1, 0), (2, 1), (3, -1)), ((1, 1), (0, 1), (2, 0)), ((3, 0), (1, 2), (1, -1)):
        QR = qconic(co, R)
        if QR == (0, 0): continue
        B = ga(ga(gm(co[0], gm(P[0], R[0])), gm(co[1], gm(P[1], R[1]))), gm(co[2], gm(P[2], R[2])))
        P2 = tuple(gs(gm(QR, P[i]), gm(gsc(2, B), R[i])) for i in range(3))
        P2 = primitive([q_(z) for z in P2])
        if qconic(co, P2) == (0, 0) and any(gm(P2[i], P[j]) != gm(P2[j], P[i]) for i in range(3) for j in range(3)):
            return P2
    raise AssertionError('no second point')

def lift4(D, places, rng, npts=12):
    """Is there lam = (lam1, lam2) in K(S,2)^2 with (l_v l_w (P_v), l_u l_w (P_v)) = lam modulo squares for some P_v in
    D_a(K_v), at every place of S?  (Then D_a lifts to an everywhere locally soluble 4-covering, i.e. a is in the image of
    Sel^4; this is equivalent to <b, D_a> = 0 for all b in Sel^2, computed here without Hilbert symbols.)  The local sets
    Lambda_v are cosets of the image of E(K_v)/2E(K_v) (order 4, or 16 at 1 + i); they are collected from local points until
    the differences of the values found span a space of that order.  Returns (solvable, dimension of the solution space)."""
    SU = D.SU; g = SU.g; rows = []
    for L in places:
        d = L.dim; target = 16 if L.kind == 'ram' else 4
        Lam = set()
        for rnd in range(400 if L.kind == 'ram' else 40):
            D.cache.pop(L.name(), None)
            for pt0, pr, _ in D.local_points(L, rng, npts):
                for sg in SIGNS:
                    pt = [pt0[0] if sg[0] == 1 else L.neg(pt0[0]), pt0[1] if sg[1] == 1 else L.neg(pt0[1]),
                          pt0[2] if sg[2] == 1 else L.neg(pt0[2]), pt0[3]]
                    ls = {}
                    for k in 'uvw':
                        x_ = D.lval(L, k, pt); vx = L.val(x_)
                        if vx is None or vx + 8 > pr: ls = None; break
                        ls[k] = L.vec(x_)
                    if ls is None: continue
                    g1 = tuple((p_ + q_) % 2 for p_, q_ in zip(ls['v'], ls['w'])); g2 = tuple((p_ + q_) % 2 for p_, q_ in zip(ls['u'], ls['w']))
                    Lam.add(g1 + g2)
            if Lam:
                lam0 = min(Lam); W = f2_span([to_int_vec(tuple((x ^ y) for x, y in zip(l, lam0))) for l in Lam])
                assert len(W) <= target.bit_length() - 1, ('Lambda_v is too large', L.name())
                if len(W) == target.bit_length() - 1: break
        # Lambda_v is a coset of the image of E(K_v)/2E(K_v) (order 4, or 16 at 1 + i): once the differences of the values
        # found span a space of that order, Lambda_v = lam0 + that space
        assert len(W) == target.bit_length() - 1, ('Lambda_v not found', L.name(), len(Lam))
        # conditions phi(vec_v(lam)) = phi(lam0) for phi in the annihilator of W (vectors of length 2d)
        Mb = [L.gvec(b) for b in SU.basis]
        for phi in kd.nullspace(f2_span(W), 2 * d):
            row = 0
            for k in range(g):
                a_ = sum((phi >> t & 1) & Mb[k][t] for t in range(d)) & 1
                c_ = sum((phi >> (d + t) & 1) & Mb[k][t] for t in range(d)) & 1
                row |= a_ << k | c_ << (g + k)
            rhs = sum((phi >> t & 1) & lam0[t] for t in range(2 * d)) & 1
            rows.append((row, rhs))
    # Gaussian elimination over F_2 with right-hand side
    piv = []
    for row, rhs in rows:
        for pr_, pv, prhs in piv:
            if row >> pv & 1: row ^= pr_; rhs ^= prhs
        if row == 0:
            if rhs: return False, None
            continue
        pv = row.bit_length() - 1
        piv = [(r_ ^ row, p_, h_ ^ rhs) if r_ >> pv & 1 else (r_, p_, h_) for r_, p_, h_ in piv]
        piv.append((row, pv, rhs))
    return True, 2 * g - len(piv)

def to_int_vec(v): return sum(b << i for i, b in enumerate(v))

def all_places(SU):
    out = []
    for p in SU.primes: out += [Loc(p, i) for i in range(2 if p % 4 == 1 else 1)]
    return out

def ctp_matrix(e, rng, say=None):
    """Sel^2 over K, the covers, the pairing matrix (rows: D_a, columns: b), torsion images; checks."""
    SU, sel, cons = selmer2(e)
    g = SU.g; places = all_places(SU)
    e1, e2, e3 = e
    tors = [(SU.vec(gm(gs(e1, e2), gs(e1, e3))), SU.vec(gs(e1, e2))), (SU.vec(gs(e2, e1)), SU.vec(gm(gs(e2, e1), gs(e2, e3))))]
    selv = [v1 | v2 << g for v1, v2 in sel]
    for t in tors: assert in_span(t[0] | t[1] << g, selv), 'a 2-torsion image outside Sel^2'
    els = [(SU.elt(v1), SU.elt(v2)) for v1, v2 in sel]
    covers = [Cover(e, a[0], a[1], SU, rng) for a in els]
    M = [[covers[i].pairing(els[j], places, rng)[0] for j in range(len(els))] for i in range(len(els))]
    for i in range(len(els)):
        assert M[i][i] == 0, 'the pairing is not alternating'
        for j in range(len(els)): assert M[i][j] == M[j][i], 'the pairing is not symmetric'
    # the torsion images pair trivially with every D_a
    for i in range(len(els)):
        for t in tors:
            assert covers[i].pairing((SU.elt(t[0]), SU.elt(t[1])), places, rng)[0] == 0, 'a torsion image pairs non-trivially'
    rho = f2_rank([sum(M[i][j] << j for j in range(len(els))) for i in range(len(els))])
    return dict(SU=SU, sel=sel, els=els, covers=covers, M=M, rho=rho, tors=tors, places=places, cons=cons)

def controls_K(seed=7, N=40, B=24):
    """Curves over K that are not base changes (j not in Q): y^2 = x (x - e2)(x - e3), e2, e3 small Gaussian integers.
    Lower bound: points with x = X / D^2 (|re X|, |im X| <= B, D in a small set) found by search and checked exactly; their
    Kummer images with those of E[2] span W, and rank >= dim W - 2.  Upper bound: dim Sel^2 - 2 - rank of the pairing.
    Checked on every curve: the images of the points found pair to 0 with every D_a, lower <= upper, and a basis class
    lifts to an everywhere locally soluble 4-covering (lift4, no Hilbert symbols) exactly when its row of the pairing is 0."""
    rng = random.Random(seed); st = dict(curves=0, sound=0, rho_pos=0, sharp=0, sharp_rho_pos=0)
    Ds = [(1, 0), (1, 1), (2, 0), (2, 1), (1, 2), (3, 0)]
    tried = 0
    while st['curves'] < N and tried < 10 * N:
        tried += 1
        e = [(0, 0), (rng.randint(-40, 40), rng.randint(-40, 40)), (rng.randint(-40, 40), rng.randint(-40, 40))]
        if len(set(e)) < 3 or jinv(e)[1] == 0: continue
        try:
            if len(selmer2(e)[1]) < 4: continue          # only curves where the pairing can matter (dim Sel^2 - 2 >= 2)
            R = ctp_matrix(e, random.Random(1))
        except AssertionError as ex:
            if 'local image not found' in str(ex) or 'no local points' in str(ex): continue
            raise
        SU = R['SU']; g = SU.g
        e1, e2, e3 = [q_(z) for z in e]
        imgs = [t[0] | t[1] << g for t in R['tors']]; pts = 0
        for D in Ds:
            D2 = gm(D, D)
            for X0 in range(-B, B + 1):
                for X1 in range(-B, B + 1):
                    x = qdiv(q_((X0, X1)), q_(D2))
                    f = qm(qm(qs(x, e1), qs(x, e2)), qs(x, e3))
                    if f == (0, 0) or q_is_square(f) is None: continue
                    v = SU.qvec(qs(x, e1)) | SU.qvec(qs(x, e2)) << g
                    if not in_span(v, imgs):
                        imgs.append(v); pts += 1
                        b = (SU.elt(v & ((1 << g) - 1)), SU.elt(v >> g))
                        for Dd in R['covers']:
                            assert Dd.pairing(b, R['places'], random.Random(2))[0] == 0, 'the image of a K-point pairs non-trivially'
        for i, Dd in enumerate(R['covers']):
            ok_, _ = lift4(Dd, R['places'], random.Random(9))
            assert ok_ == (not any(R['M'][i])), 'the lift to a 4-covering disagrees with the pairing'
        lo = f2_rank(imgs) - 2; up = len(R['sel']) - 2 - R['rho']
        st['curves'] += 1; st['sound'] += lo <= up; st['rho_pos'] += R['rho'] > 0
        st['sharp'] += lo == up; st['sharp_rho_pos'] += (lo == up and R['rho'] > 0)
        say('  y^2 = x(x - (%d%+di))(x - (%d%+di)): dim Sel^2 %d, pairing rank %d, upper %d; points found give lower %d%s'
            % (e[1][0], e[1][1], e[2][0], e[2][1], len(R['sel']), R['rho'], up, lo, '' if lo <= up else '   ABOVE THE BOUND'))
    return st

def pari_ranks_Q(a, b, tw):
    pr = PARI()
    E = pr.ellinit([0, tw * (a + b), 0, a * b, 0]); lo, hi = 0, 10 ** 9
    for F in pr.ellisomat(pr.ellminimalmodel(E), 0, 1)[0]:
        r = pr.ellrank(pr.ellinit(F), 1); lo = max(lo, int(r[0])); hi = min(hi, int(r[1]))
    return lo, hi

def controls(seed=3, N=200):
    """y^2 = x(x + a)(x + b) over Q, base-changed to K: rank E(K) = rank E(Q) + rank E^(-1)(Q) (PARI, a guide).  The bound
    from Sel^2 and the pairing must never be below that rank; where the 2-Selmer bound exceeds it, the pairing should
    (usually) close the gap."""
    rng = random.Random(seed); stats = dict(curves=0, sound=0, excess2=0, closed=0, zero_ok=0, zero_curves=0)
    seen = set()
    for _ in range(N):
        a = rng.randint(-80, 80); b = rng.randint(-80, 80)
        if a == 0 or b == 0 or a == b or (a, b) in seen: continue
        seen.add((a, b)); seen.add((b, a))
        e = [(0, 0), (-a, 0), (-b, 0)]
        kb = kd.selmer_bound(e)
        if kb is None: continue
        l1, h1 = pari_ranks_Q(a, b, 1); l2, h2 = pari_ranks_Q(a, b, -1)
        if not (l1 == h1 and l2 == h2): continue
        rK = l1 + l2
        if kb == rK and rK == 0: continue                  # nothing to test
        if kb - rK not in (0, 2): continue
        if kb == rK and stats['zero_curves'] >= 12: continue
        R = ctp_matrix(e, random.Random(1))
        up = len(R['sel']) - 2 - R['rho']
        stats['curves'] += 1; stats['sound'] += up >= rK
        if kb == rK: stats['zero_curves'] += 1; stats['zero_ok'] += R['rho'] == 0
        else: stats['excess2'] += 1; stats['closed'] += up == rK
        say('  y^2 = x(x %+d)(x %+d): ranks over Q of E, E^(-1) (PARI) %d, %d; dim Sel^2(E/K) %d; pairing rank %d; bound %d%s'
            % (a, b, l1, l2, len(R['sel']), R['rho'], up, '' if up >= rK else '   BELOW THE RANK'))
    return stats

# ------------------------------------------------------------------ route (b): a search for points over K
def split_primes(lo, hi, avoid):
    out = []
    for p in range(lo, hi):
        if p % 4 == 1 and is_prime(p) and p not in avoid: out.append(p)
    return out

def quartic_points(curve, H, nprimes=10, say_=None):
    """s = sigma/tau with sigma, tau in Z[i], |re|, |im| <= H (tau up to units), Q_J(sigma, tau) a square in Z[i]: a
    numpy sieve modulo split primes (both embeddings), then the exact test.  Returns the list of s (Gaussian rationals)."""
    import numpy as np
    F = factors(627, 364) if curve.get('mn') is None else factors(*curve['mn'])
    names = [curve['J'][i:i + 2] for i in range(0, 8, 2)]
    lins = [F[k] for k in names]
    ps = split_primes(1000, 4000, set())[:nprimes]
    rs = [sqrt_mod_p(p - 1, p) for p in ps]
    rng_ = np.arange(-H, H + 1, dtype=np.int64)
    S0, S1 = np.meshgrid(rng_, rng_, indexing='ij'); S0 = S0.ravel(); S1 = S1.ravel()
    # quadratic residue tables
    qr = {}
    for p in ps:
        t = np.zeros(p, dtype=bool); t[(np.arange(1, p, dtype=np.int64) ** 2) % p] = True; t[0] = True; qr[p] = t
    found = []
    taus = [(t0, t1) for t0 in range(0, H + 1) for t1 in range(0, H + 1) if (t0 > 0) and (t0, t1) != (0, 0)] + [(0, 0)]
    for tau in taus:
        if tau == (0, 0):
            sig_list = [(1, 0)]
            cand = np.array([True])
            S0c, S1c = np.array([1]), np.array([0])
        else:
            mask = np.ones(S0.shape, dtype=bool)
            for p, r0 in zip(ps, rs):
                for r in (r0, p - r0):
                    sv = (S0 + S1 * r) % p; tv = (tau[0] + tau[1] * r) % p
                    val = np.ones(S0.shape, dtype=np.int64)
                    for (c1, c0) in lins:
                        a1 = (c1[0] + c1[1] * r) % p; a0 = (c0[0] + c0[1] * r) % p
                        val = val * ((a1 * sv + a0 * tv) % p) % p
                    mask &= qr[p][val]
            S0c, S1c = S0[mask], S1[mask]
        for s0_, s1_ in zip(S0c.tolist(), S1c.tolist()):
            sig = (s0_, s1_)
            if sig == (0, 0) and tau == (0, 0): continue
            if tau != (0, 0) and gn(ggcd(sig, tau) if sig != (0, 0) else tau) != 1: continue
            val = (1, 0)
            for (c1, c0) in lins: val = gm(val, ga(gm(c1, sig), gm(c0, tau)))
            if val == (0, 0): continue
            if gsqrt_exact(val) is not None:
                found.append((sig, tau))
    return found

def kummer_s(curve, SU, s0):
    """the Kummer image (classes of x - e1, x - e2) of the point of E with s = s0 (a K-point of Y^2 = Q_J(s), s0 != r)."""
    den = curve['den']; A = curve['A']; r = curve['r']
    x = qdiv((A[0] * den * den, A[1] * den * den), qs(s0, r))
    e1, e2, e3 = [q_(z) for z in curve['e']]
    f = qm(qm(qs(x, e1), qs(x, e2)), qs(x, e3))
    assert q_is_square(f) is not None, 'not a point of E'
    return (SU.qvec(qs(x, e1)), SU.qvec(qs(x, e2))), x

def ec_add(e, P, Q):
    """the group law on y^2 = (x - e1)(x - e2)(x - e3) over K (Gaussian rationals; None is O)."""
    if P is None: return Q
    if Q is None: return P
    e1, e2, e3 = [q_(z) for z in e]
    a2 = (-(e1[0] + e2[0] + e3[0]), -(e1[1] + e2[1] + e3[1])); a4 = qa(qa(qm(e1, e2), qm(e1, e3)), qm(e2, e3))
    if P[0] == Q[0]:
        if qa(P[1], Q[1]) == (0, 0): return None
        lam = qdiv(qa(qa(qm(q_(3), qm(P[0], P[0])), qm(q_(2), qm(a2, P[0]))), a4), qm(q_(2), P[1]))
    else: lam = qdiv(qs(Q[1], P[1]), qs(Q[0], P[0]))
    x3 = qs(qs(qs(qm(lam, lam), a2), P[0]), Q[0])
    return (x3, qs(qm(lam, qs(P[0], x3)), P[1]))

def order4_check(curve):
    """the points over s = -i, i have order 4 and twice them is a point of order 2 (exact group law)."""
    e = curve['e']; out = []
    for s0 in ((Fr(0), Fr(-1)), (Fr(0), Fr(1))):
        den = curve['den']; A = curve['A']; r = curve['r']
        x = qdiv((A[0] * den * den, A[1] * den * den), qs(s0, r))
        f = qm(qm(qs(x, q_(e[0])), qs(x, q_(e[1]))), qs(x, q_(e[2]))); y = q_is_square(f); assert y is not None
        P = (x, y); P2 = ec_add(e, P, P)
        assert P2 is not None and P2[1] == (0, 0) and ec_add(e, P2, P2) is None
        out.append([k for k in range(3) if P2[0] == q_(e[k])][0] + 1)
    return out

def torsion_image(curve, SU):
    """images of T1, T2 and of the points over the degenerate solutions s = -i, i (order 4, checked by doubling below)."""
    e1, e2, e3 = curve['e']
    vs = [(SU.vec(gm(gs(e1, e2), gs(e1, e3))), SU.vec(gs(e1, e2))), (SU.vec(gs(e2, e1)), SU.vec(gm(gs(e2, e1), gs(e2, e3))))]
    for s0 in ((Fr(0), Fr(-1)), (Fr(0), Fr(1))): vs.append(kummer_s(curve, SU, s0)[0])
    return vs

def cover_model(curve, SU, a, delta=0):
    """the 2-covering of class a = (v1, v2) as  y3^2 = al3 y1^2 + ga3 y2^2,  y4^2 = al4 y1^2 + ga4 y2^2  over Z[i]
    (from L_j(sigma, tau) = beta_j y_j^2, beta_4 = 1, beta_k = a_k kappa_k), and the map back to s."""
    F = factors(*curve['mn']); names = [curve['J'][i:i + 2] for i in range(0, 8, 2)]
    lins = [(q_(F[k][0]), q_(F[k][1])) for k in names]
    A = curve['A']; r = curve['r']; c14 = lins[3][0]
    rk = [qdiv((-lins[k][1][0], -lins[k][1][1]), lins[k][0]) for k in range(3)]
    kap = [SU.qvec(qm(qm(qm((-A[0], -A[1]), c14), lins[k][0]), qs(rk[k], r))) for k in range(3)]
    av = [a[0], a[1], a[0] ^ a[1]]
    beta = [q_(SU.elt(av[k] ^ kap[k] ^ delta)) for k in range(3)] + [q_(SU.elt(delta))]
    assert (av[0] ^ kap[0]) ^ (av[1] ^ kap[1]) ^ (av[2] ^ kap[2]) == 0, 'the product of the betas is not a square'
    # (sigma, tau) = Minv (beta1 Y1, beta2 Y2)
    (p1, q1), (p2, q2) = lins[0], lins[1]
    det = qs(qm(p1, q2), qm(q1, p2))
    Minv = [[qdiv(q2, det), qdiv((-q1[0], -q1[1]), det)], [qdiv((-p2[0], -p2[1]), det), qdiv(p1, det)]]
    def coef(k):
        P = qa(qm(lins[k][0], Minv[0][0]), qm(lins[k][1], Minv[1][0])); Q = qa(qm(lins[k][0], Minv[0][1]), qm(lins[k][1], Minv[1][1]))
        return qdiv(qm(P, beta[0]), beta[k]), qdiv(qm(Q, beta[1]), beta[k])
    model = []
    for k in (2, 3):
        al, ga_ = coef(k); d = 1
        for z in (al, ga_):
            for y in z: d = d * y.denominator // gcd(d, y.denominator)
        model.append(((int(al[0] * d * d), int(al[1] * d * d)), (int(ga_[0] * d * d), int(ga_[1] * d * d))))
    def to_s(y1, y2):
        Y1 = qm(beta[0], q_(gm(y1, y1))); Y2 = qm(beta[1], q_(gm(y2, y2)))
        sig = qa(qm(Minv[0][0], Y1), qm(Minv[0][1], Y2)); tau = qa(qm(Minv[1][0], Y1), qm(Minv[1][1], Y2))
        return None if tau == (0, 0) else qdiv(sig, tau)
    return model, to_s, beta

def cover_search(model, H, nprimes=12):
    """(y1, y2) in Z[i]^2, |re|, |im| <= H, y2 up to units, with al y1^2 + ga y2^2 a square in Z[i] for both equations:
    numpy sieve modulo split primes, then the exact test."""
    import numpy as np
    ps = split_primes(1000, 6000, set())[:nprimes]; rs = [sqrt_mod_p(p - 1, p) for p in ps]
    g = np.arange(-H, H + 1, dtype=np.int64)
    A0, A1 = np.meshgrid(g, g, indexing='ij'); A0 = A0.ravel(); A1 = A1.ravel()
    qr = {}
    for p in ps:
        t = np.zeros(p, dtype=bool); t[(np.arange(0, p, dtype=np.int64) ** 2) % p] = True; qr[p] = t
    pre = []
    for p, r0 in zip(ps, rs):
        for r in (r0, p - r0):
            y1 = (A0 + A1 * r) % p; y1sq = y1 * y1 % p
            terms = []
            for (al, ga_) in model:
                terms.append((((al[0] + al[1] * r) % p) * y1sq % p, (ga_[0] + ga_[1] * r) % p))
            pre.append((p, r, terms))
    out = []
    y2s = [(b0, b1) for b0 in range(0, H + 1) for b1 in range(0, H + 1) if b0 > 0] + [(0, 0)]
    for y2 in y2s:
        idx = np.arange(A0.shape[0])
        for p, r, terms in pre:
            if idx.size == 0: break
            v2 = (y2[0] + y2[1] * r) % p; v2sq = v2 * v2 % p
            keep = np.ones(idx.size, dtype=bool)
            for t1, gco in terms: keep &= qr[p][(t1[idx] + gco * v2sq) % p]
            idx = idx[keep]
        for j in idx.tolist():
            y1 = (int(A0[j]), int(A1[j]))
            if y1 == (0, 0) and y2 == (0, 0): continue
            ok = True; ys = []
            for (al, ga_) in model:
                w = ga(gm(al, gm(y1, y1)), gm(ga_, gm(y2, y2)))
                sq = gsqrt_exact(w)
                if sq is None: ok = False; break
                ys.append(sq)
            if ok: out.append((y1, y2, ys))
    return out

def resultant_classes(curve, SU):
    """squarefree classes (S-unit vectors) that beta_4 = class of L_4(sigma, tau) can take for coprime sigma, tau: the units
    {1, i} times products of the Gaussian primes dividing the resultants of L_4 with the other three factors."""
    F = factors(*curve['mn']); names = [curve['J'][i:i + 2] for i in range(0, 8, 2)]
    c4 = F[names[3]]; ks = set()
    for k in names[:3]:
        c = F[k]; R = gs(gm(c4[0], c[1]), gm(c4[1], c[0]))
        v = SU.vec(R, allow_square_cofactor=True)
        for j in range(1, SU.g):
            if gdivides(SU.basis[j], R): ks.add(j)
    ks = sorted(ks); out = []
    for mask in range(1 << (len(ks) + 1)):
        v = mask & 1
        for t, j in enumerate(ks):
            if mask >> (t + 1) & 1: v |= 1 << j
        out.append(v)
    return out

def admissible_deltas(curve, SU, a):
    """the delta = beta_4 for which every beta_j (j = 1..4) divides the product of the resultants of L_j with the
    other three factors (necessary for a point with sigma, tau coprime in Z[i] to sit at integral y in that model)."""
    F = factors(*curve['mn']); names = [curve['J'][i:i + 2] for i in range(0, 8, 2)]
    allowed = []
    for j in range(4):
        cj = F[names[j]]; ks = 1                      # bit 0: the unit i is always allowed
        for k in range(4):
            if k == j: continue
            ck = F[names[k]]; R = gs(gm(cj[0], ck[1]), gm(cj[1], ck[0]))
            for t in range(1, SU.g):
                if gdivides(SU.basis[t], R): ks |= 1 << t
        allowed.append(ks)
    out = []
    A = curve['A']; r = curve['r']; lins = [(q_(F[k][0]), q_(F[k][1])) for k in names]
    rk = [qdiv((-lins[k][1][0], -lins[k][1][1]), lins[k][0]) for k in range(3)]
    kap = [SU.qvec(qm(qm(qm((-A[0], -A[1]), lins[3][0]), lins[k][0]), qs(rk[k], r))) for k in range(3)]
    av = [a[0], a[1], a[0] ^ a[1]]
    for d in resultant_classes(curve, SU):
        bs = [av[k] ^ kap[k] ^ d for k in range(3)] + [d]
        if all(bs[j] & ~allowed[j] == 0 for j in range(4)): out.append(d)
    return out

# ------------------------------------------------------------------ the runs
CURVES = ['L1L3M1M2', 'L1L3M3M4', 'L2L4M1M2', 'L2L4M3M4']

def isog_full(e, i):
    """E/<T_i> as y^2 = x (x - e2')(x - e3') if all its points of order 2 are over K (then (e_j - e_i)(e_k - e_i) is a
    square w^2 and e' = 0, -(d_j + d_k) + 2w, -(d_j + d_k) - 2w), else None."""
    j, k = [t for t in range(3) if t != i]
    dj, dk = gs(e[j], e[i]), gs(e[k], e[i])
    w = gsqrt_exact(gm(dj, dk))
    if w is None: return None
    s_ = ga(dj, dk)
    return [(0, 0), ga(gsc(-1, s_), gsc(2, w)), gs(gsc(-1, s_), gsc(2, w))]

def jinv(e):
    """j-invariant of y^2 = (x - e1)(x - e2)(x - e3): 256 (l^2 - l + 1)^3 / (l^2 (l - 1)^2), l = (e3 - e1)/(e2 - e1)."""
    l = qdiv(q_(gs(e[2], e[0])), q_(gs(e[1], e[0]))); one = (Fr(1), Fr(0))
    t = qa(qs(qm(l, l), l), one); lm1 = qs(l, one)
    return qdiv(qm(q_(256), qm(qm(t, t), t)), qm(qm(l, l), qm(lm1, lm1)))

def iso_over_K(e, f):
    """is y^2 = prod (x - f_i) isomorphic over K to y^2 = prod (x - e_i)?  f_s(i) = u^2 e_i + t for a permutation s and
    u^2 a square in K (exact)."""
    for s_ in itertools.permutations(range(3)):
        u2 = qdiv(q_(gs(f[s_[1]], f[s_[0]])), q_(gs(e[1], e[0])))
        if qm(u2, q_(gs(e[2], e[0]))) != q_(gs(f[s_[2]], f[s_[0]])): continue
        if q_is_square(u2) is not None: return True
    return False

def coset_reps(SU, sel, tv):
    g = SU.g; T = f2_span([a | b << g for a, b in tv]); reps = []
    for a, b in sel:
        v = a | b << g
        if not in_span(v, T + reps): reps.append(v)
    return T, reps

def search_points(c, SU, sel, tv, H):
    """route (b): K-points on the 2-coverings of the non-torsion cosets, in every admissible model, |y| <= H."""
    g = SU.g; T, reps = coset_reps(SU, sel, tv)
    assert len(reps) == 2
    cos = [reps[0], reps[1], reps[0] ^ reps[1]]; found = []
    for v in cos:
        a = (v & ((1 << g) - 1), v >> g); ds = admissible_deltas(c, SU, a)
        for d in ds:
            model, to_s, beta = cover_model(c, SU, a, d)
            for y1, y2, ys in cover_search(model, H):
                s0 = to_s(y1, y2)
                if s0 is None: continue
                kv, x = kummer_s(c, SU, s0)
                assert in_span((kv[0] | kv[1] << g) ^ v, T), 'a point found on D_a has another Kummer image'
                found.append((v, s0, kv))
    # torsion classes as a test of the models: all eight torsion points must be found
    tors_s = set()
    for a in [(0, 0)] + tv:
        for d in admissible_deltas(c, SU, a):
            model, to_s, beta = cover_model(c, SU, a, d)
            for y1, y2, ys in cover_search(model, 4):
                s0 = to_s(y1, y2)
                if s0 is not None: tors_s.add(s0)
    return found, tors_s, len(cos)

def main(H=24):
    say('kctp.py: the Cassels-Tate pairing on Sel^2 over K = Q(i) for the four corner K-quartics of (627,364).')
    say('Written by an AI model; not read by a person.  Experiment, not counted.')
    say()
    setup_h2()
    say('Hilbert symbols at 1 + i (product formula on the basis %s of K_v*/K_v*^2): matrix %s' % (RB, H2))
    nr, npari = hilbert_selftest()
    say('  check: product formula on %d random pairs (all places)%s: agree' % (nr, '' if NOPARI else ', PARI nfhilbert at every place for %d of them' % npari))
    say()
    m, n = 627, 364; res = {}
    F = factors(m, n); S_ = m * m + n * n; D_ = m * m - n * n
    qp = pmul([q_(F['M1'][1]), q_(F['M1'][0])], [q_(F['M4'][1]), q_(F['M4'][0])])
    qmn = pmul([q_(F['M2'][1]), q_(F['M2'][0])], [q_(F['M3'][1]), q_(F['M3'][0])])
    assert qp == [q_(S_), q_(-2 * D_), q_(S_)] and qmn == [q_(S_), q_(2 * D_), q_(S_)]
    say('Labels: q+ = S s^2 - 2 D s + S = M1 M4 and q- = M2 M3 (checked), so L1L3M1M2 is not q+ L1 L3 (= L1 L3 M1 M4).')
    say()
    for J in CURVES:
        c = build_curve(m, n, J); c['mn'] = (m, n)
        say('== %s of (%d,%d): y^2 = (x - e1)(x - e2)(x - e3) over K, s = r + den^2 A / x, den = %d' % (J, m, n, c['den']))
        say('   e = %s' % (c['e'],))
        say('   Q_J at s = -i, i (t = 1, -1): squares in K (checked), so the class of kclose.py is 1')
        jj = jinv(c['e']); assert jj[1] != 0
        say('   j(E) has imaginary part %s != 0: j is not rational (exact)' % (('%.6e' % float(jj[1])),))
        rng = random.Random(1)
        R = ctp_matrix(c['e'], rng)
        SU = R['SU']; g = SU.g
        tv = torsion_image(c, SU)
        T = f2_span([a | b << g for a, b in tv])
        o4 = order4_check(c)
        say('   the points over s = -i, i have order 4, twice them T%d, T%d (exact group law)' % tuple(o4))
        say('   S: places above %s; dim Sel^2(E/K) = %d; torsion image (T1, T2, the points over s = -i, i) of dimension %d' % (SU.primes, len(R['sel']), len(T)))
        for a, b in R['sel']:
            say('     (%s, %s)' % ('*'.join(str(SU.basis[k]) for k in range(g) if a >> k & 1) or '1', '*'.join(str(SU.basis[k]) for k in range(g) if b >> k & 1) or '1'))
        selv = [a | b << g for a, b in R['sel']]
        for t in T: assert in_span(t, selv)
        # all torsion images pair trivially (also the order-4 points)
        for i, D in enumerate(R['covers']):
            for a, b in tv: assert D.pairing((SU.elt(a), SU.elt(b)), R['places'], rng)[0] == 0
        # the same matrix with the tangent lines at other conic points (D_a rebuilt with a second rational point)
        say('   pairing matrix <b_j, D_(a_i)>: %s (symmetric, alternating, torsion images in the kernel: checked)' % R['M'])
        alt = [Cover(c['e'], a[0], a[1], SU, rng, forced_points={k: other_point(D.co[k], D.P[k]) for k in 'uvw'}) for a, D in zip(R['els'], R['covers'])]
        M2 = [[alt[i].pairing(R['els'][j], R['places'], rng)[0] for j in range(len(alt))] for i in range(len(alt))]
        assert M2 == R['M'], 'the pairing depends on the choice of the conic points'
        say('   the same matrix with the tangent lines at second rational points of all twelve conics')
        good = [p for p in (3, 7, 11, 13, 17, 19, 23, 41, 53, 61) if p not in SU.primes]
        gl = [Loc(p, i) for p in good for i in range(2 if p % 4 == 1 else 1)]
        for D in R['covers'] + alt:
            for b in R['els']:
                for L in gl:
                    vals = D.local_term(L, b, rng); assert vals == {0}, ('a place outside S contributes', L.name())
        say('   check of the purity argument: the local term is 0 at the %d places above %s (outside S), for all D_a and b' % (len(gl), good))
        lf = [lift4(D, R['places'], random.Random(9)) for D in R['covers']]
        for (s, _), row in zip(lf, R['M']): assert s == (not any(row)), 'the lift to a 4-covering disagrees with the pairing'
        assert all(s for s, _ in lf)
        say('   lifts to 4-coverings (no Hilbert symbols: lam in K(S,2)^2 with (l_v l_w, l_u l_w)(P_v) = lam locally, every v in S):'
            ' solvable for all %d basis classes (solution spaces of dimension %s): every class of Sel^2 lifts to an everywhere'
            ' locally soluble 4-covering, as the zero pairing requires' % (len(lf), [d for _, d in lf]))
        up = len(R['sel']) - len(T) - R['rho']
        say('   rank of the pairing %d  =>  rank E(K) <= %d' % (R['rho'], up))
        # the 2-isogenous curves with all 2-torsion over K
        for i in range(3):
            e2 = isog_full(c['e'], i)
            if e2 is None: say('   E/<T%d>: not all points of order 2 over K' % (i + 1)); continue
            R2 = ctp_matrix(e2, random.Random(1))
            say('   E/<T%d> = y^2 = x(x - %s)(x - %s): dim Sel^2 = %d, pairing rank %d  =>  rank <= %d' % (i + 1, e2[1], e2[2], len(R2['sel']), R2['rho'], len(R2['sel']) - 2 - R2['rho']))
            nxt = [isog_full(e2, k) is None for k in (1, 2)]
            assert all(nxt)
            say('     its 2-isogenous curves other than E (kernels its T2, T3; T1 = (0, 0) is the kernel of the dual): not all points of order 2 over K')
            up = min(up, len(R2['sel']) - 2 - R2['rho'])
        found, tors_s, ncos = search_points(c, SU, R['sel'], tv, H)
        say('   route (b): points on the 2-coverings of the %d non-torsion cosets, all admissible models, |y| <= %d: %d found' % (ncos, H, len(found)))
        say('     (test of the models, |y| <= 4: the torsion classes give s = %s, the s of all 8 torsion points: O at r, T1..T3 at the roots, the order-4 points at -i, i)' % sorted('%s%+si' % (t[0], t[1]) for t in tors_s))
        assert len(tors_s) == 6
        res[J] = dict(sel=len(R['sel']), rho=R['rho'], up=up, found=len(found))
        say()
    es = {J: build_curve(m, n, J)['e'] for J in CURVES}
    conj1 = [gc(z) for z in es[CURVES[0]]]
    iso = (iso_over_K(es[CURVES[0]], es[CURVES[3]]), iso_over_K(es[CURVES[1]], es[CURVES[2]]), iso_over_K(conj1, es[CURVES[1]]))
    say('Isomorphisms over K (exact): %s ~ %s: %s; %s ~ %s: %s; %s ~ complex conjugate of %s: %s' % (CURVES[3], CURVES[0], iso[0], CURVES[2], CURVES[1], iso[1], CURVES[1], CURVES[0], iso[2]))
    assert all(iso)
    say('So the four curves have one rank over K between them.')
    say()
    say('RESULT: on all four curves dim Sel^2 = 4, the torsion image has dimension 2, and the Cassels-Tate pairing on Sel^2 is')
    say('identically zero (also on the 2-isogenous curve with full 2-torsion).  So rank E(K) <= 2, and the rank is 0 or 2 if')
    say('Sha(E/K)[2^inf] is finite (PARI\'s root number +1, a guide, agrees): the level-2 pairing does not decide.  Every class')
    say('of Sel^2 lifts to an everywhere locally soluble 4-covering.  No point of infinite order was found.')
    say('(627,364) is NOT closed by this computation.')
    assert all(r['rho'] == 0 and r['up'] == 2 for r in res.values())
    return res

def negative():
    """forged inputs that must be refused."""
    global NPTS
    say('kctp.py --negative: forged inputs; each must be refused (10 local points per place here, 6 in the main run).')
    setup_h2(); n_ok = 0; tests = 0
    c = build_curve(627, 364, CURVES[0]); c['mn'] = (627, 364)
    rng = random.Random(1)
    R = ctp_matrix(c['e'], rng); SU = R['SU']; g = SU.g; places = R['places']
    selv = [a | b << g for a, b in R['sel']]
    def refused(name, f):
        nonlocal n_ok, tests
        tests += 1
        try: f()
        except AssertionError as ex:
            n_ok += 1; say('  refused: %s  (%s)' % (name, str(ex)[:110])); return
        say('  NOT REFUSED: %s' % name)
    D = R['covers'][1]; a1, a2, a3 = D.a
    NPTS = 10
    for Dd in R['covers']: Dd.cache = {}
    # 1. a point off the conic
    P = D.P['w']; Pbad = (ga(P[0], (1, 0)), P[1], P[2])
    refused('a conic point off the conic C_w', lambda: Cover(c['e'], a1, a2, SU, rng, forced_points={'w': Pbad}))
    # 2. a forged conic hint
    key = str([list(z) for z in D.co['w']]); old = HINTS['conic'].get(key)
    HINTS['conic'][key] = [list(Pbad[0]), list(Pbad[1]), list(Pbad[2])]
    refused('a forged conic hint', lambda: solve_conic(D.co['w'], SU, rng))
    if old is None: del HINTS['conic'][key]
    else: HINTS['conic'][key] = old
    # 3. a line that is not tangent in place of l_w: the local term must stop being constant
    lw = D.l['w']; lbad = (ga(lw[0], (1, 0)), lw[1], lw[2])
    def f3():
        Db = Cover(c['e'], a1, a2, SU, random.Random(2), forced_lines={'w': lbad})
        for b in R['els']:
            Db.pairing(b, places, random.Random(3))
    refused('a non-tangent line in place of l_w', f3)
    # 4. classes b outside Sel^2: the local term must not be constant (each class separately)
    for k in range(SU.g):
        if in_span(1 << k, selv): continue
        b = (SU.basis[k], (1, 0))
        def f4(b=b):
            for Dd in R['covers'][:2]: Dd.pairing(b, places, random.Random(4))
        refused('the class (%s, 1) outside Sel^2' % (SU.basis[k],), f4)
    # 5. classes a outside Sel^2 used as 2-coverings: no local points somewhere (each separately)
    for k in range(2, SU.g, 3):
        if in_span(1 << k | 1 << (k + g), selv): continue
        def f5(k=k):
            Db = Cover(c['e'], SU.basis[k], SU.basis[k], SU, random.Random(5))
            for L in places: Db.local_points(L, random.Random(6))
        refused('the class (%s, %s) outside Sel^2 used as a 2-covering' % (SU.basis[k], SU.basis[k]), f5)
    # 6. a factorisation hint with a prime missing
    N = gn(gs(c['e'][0], c['e'][1])); keyf = str(N); oldf = HINTS['factor'].get(keyf)
    HINTS['factor'][keyf] = [str(p) for p in factor_int(N)][:-1]
    refused('an incomplete factorisation hint', lambda: factor_int(N))
    HINTS['factor'][keyf] = oldf
    # 7. a vector outside the Selmer group claimed to be in it: some local condition fails
    def f7():
        v = (1 << 2) | (1 << (2 + g))
        assert in_span(v, selv) or all(bin(row & (v)).count('1') % 2 == 0 for row in R['cons']), 'a local condition fails'
    refused('a forged Selmer element', f7)
    # 8. the lift to 4-coverings with a non-tangent line in place of l_w: the local values are no longer a coset of the
    #    right order (or the system has no solution)
    def f8b():
        Db = Cover(c['e'], a1, a2, SU, random.Random(2), forced_lines={'w': lbad})
        ok_, _ = lift4(Db, places, random.Random(9))
        assert ok_, 'no lift with a non-tangent line'
    refused('lift4 with a non-tangent line in place of l_w', f8b)
    # 9. the lift to 4-coverings on a control curve whose pairing has rank 2: a class with a non-zero row must not lift
    def f9():
        e9 = [(0, 0), (14, -33), (32, -25)]
        R9 = ctp_matrix(e9, random.Random(1))
        i = [i for i in range(len(R9['covers'])) if any(R9['M'][i])][0]
        ok_, _ = lift4(R9['covers'][i], R9['places'], random.Random(9))
        assert ok_, 'no lift: the class with a non-zero row of the pairing does not lift (as it must not)'
    refused('a lift claimed for a class of y^2 = x(x - (14-33i))(x - (32-25i)) whose pairing row is non-zero', f9)
    # 10. a forged pairing matrix entry: the symmetry check (pairing computed twice in opposite order)
    def f8():
        M = [row[:] for row in R['M']]; M[0][1] ^= 1
        for i in range(len(M)):
            for j in range(len(M)): assert M[i][j] == M[j][i], 'the pairing is not symmetric'
    refused('a pairing matrix with one entry flipped', f8)
    say('%d of %d forged inputs refused' % (n_ok, tests))
    assert n_ok == tests
    return n_ok, tests

if __name__ == '__main__':
    if '--negative' in sys.argv:
        negative(); fn = 'kctp_negative.txt'
    elif '--controls' in sys.argv:
        say('kctp.py --controls: y^2 = x(x + a)(x + b) over Q base-changed to K; ranks over Q from PARI (a guide).')
        setup_h2(); st = controls(3, 120); say('summary: %s' % st)
        assert st['sound'] == st['curves'] and st['zero_ok'] == st['zero_curves']
        say()
        say('Curves over K with j not in Q (not base changes): lower bound from points found by search, upper from Sel^2 and the pairing.')
        st2 = controls_K(7, 16); say('summary: %s' % st2)
        assert st2['sound'] == st2['curves']
        fn = 'kctp_controls.txt'
    elif any(a.startswith('--search=') for a in sys.argv):
        H = int([a for a in sys.argv if a.startswith('--search=')][0][9:])
        say('kctp.py --search=%d: points on the 2-coverings of the non-torsion cosets of Sel^2 for L1L3M1M2 of (627,364)' % H)
        say('(the other three curves are isomorphic over K to it or to its complex conjugate: checked in kctp.txt)')
        c = build_curve(627, 364, CURVES[0]); c['mn'] = (627, 364)
        SU, sel, cons = selmer2(c['e']); tv = torsion_image(c, SU)
        found, tors_s, ncos = search_points(c, SU, sel, tv, H)
        say('%d cosets x all admissible models, |re y|, |im y| <= %d: %d points of infinite order found' % (ncos, H, len(found)))
        for v, s0, kv in found: say('  coset %s: s = %s' % (hex(v), s0))
        fn = 'kctp_search.txt'
    else:
        Hs = [a for a in sys.argv if a.startswith('--H=')]
        main(int(Hs[0][4:]) if Hs else 24); fn = 'kctp.txt'
    if WRITE:
        open(os.path.join(here, fn), 'w').write('\n'.join(OUT) + '\n')
        if os.path.exists(HINTS_FILE):                 # merge with hints written by the other modes
            old = json.load(open(HINTS_FILE))
            for k in old: old[k].update(HINTS.get(k, {}))
            for k in HINTS: old.setdefault(k, HINTS[k])
            merged = old
        else: merged = HINTS
        json.dump(merged, open(HINTS_FILE, 'w'), indent=0, sort_keys=True)
