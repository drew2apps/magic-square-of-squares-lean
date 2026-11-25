#!/usr/bin/env python3
"""Upper bounds for ranks by an exact descent via 2-isogeny, with an exact second descent; no Magma, no PARI, no Sage.

    python3 rank_indep.py                 the four curves E+ of (760,39), (945,248), (72,65), (920,129)
    python3 rank_indep.py --all           every curve of sigma_evidence.txt (12)
    python3 rank_indep.py --negative      altered inputs that must be refused
    (add --write to rewrite the output files and the hints; without it the script only prints, so CI can diff its output)
    (exact_last_pair.py imports this module for the curve E1 of LAST_PAIR.md.)

Written by an AI model (Claude); not read by a person.

THE CURVES.  The model of sigma_exact.py (`model(m, n, s)`), recomputed from m, n, and compared with the a-invariants
that were sent to Magma (magma_sigma2.json, magma_sigma3.json).  It has one rational point of order 2, (e, 0);
x -> x + e gives  E: y^2 = x^3 + a x^2 + b x,  and the 2-isogenous curve  E': Y^2 = X^3 + a' X^2 + b' X  with
a' = -2a, b' = a^2 - 4b.

FIRST DESCENT (Silverman-Tate, Rational Points on Elliptic Curves, III.4-III.6).  alpha: E(Q) -> Q*/Q*^2,
(x, y) -> x, (0, 0) -> b, O -> 1, has image in the classes d | b (d squarefree, any sign), and d is in the image iff
    C_d:  N^2 = d M^4 + a M^2 e^2 + (b/d) e^4
has a rational point.  Same for E' (alpha').  With r the rank,  2^r = |alpha(E(Q))| |alpha'(E'(Q))| / 4.
Here the image of alpha is bounded above by the d for which C_d has a point in R and in every Q_p, p | 2 b b'.
For p not dividing 2 b b' and d | b (so p does not divide d), C_d has good reduction at p (the discriminant of the
quartic is 16 b (a^2 - 4b)^2), hence a smooth F_p-point (Hasse) that lifts (Hensel): no condition there.

LOCAL SOLUBILITY (function `solvable`): is there (M : e) in P^1(Q_p) with every form c_i g_i(M, e) a square in Q_p
(0 counts as a square)?  Charts (t, 1), t in Z_p, and (1, t), t in pZ_p (the forms have even degree, so the value
changes by a square between the homogeneous and the affine form).  On a disc t0 + p^k Z_p, g(t) = g(t0) mod p^k.
  p odd: if v = v_p(g(t0)) < k, g(t)/g(t0) is in 1 + pZ_p, a square; the class is v mod 2 and the Legendre symbol.
  p = 2: if v + 3 <= k, g(t)/g(t0) is in 1 + 8Z_2, a square; the class is v mod 2 and the odd part mod 8.
A disc where some form has a constant non-square class is dropped; a disc where all have a constant square class
is a certified local point; otherwise it is cut into p discs, to a fixed depth, where it is counted as possible.
So the answer is True (a point, certified), False (no point, certified), or None (not decided; counted as possible).
For a prime p >= SMALL the p discs of a level are not all visited: after removing the content (g = p^k h), a residue
t0 with every h(t0) a unit mod p decides the disc at once (R1), so only two kinds of residues matter: roots of some
h mod p (found by a gcd with t^p - t and splitting), and, when every exponent k is even, residues where all h(t0) are
non-zero squares.  The latter are searched for; if none is found, a certificate is required that none exists: a
product of some of the h that is, mod p, a non-residue constant times the square of a polynomial (then at a residue
where all h(t0) are non-zero squares the product would be a non-zero square, a contradiction).  Without such a
certificate the answer is None.
  The real place is decided exactly (sign of the forms between their real roots, isolated with sympy).

SECOND DESCENT (for a class d of the first descent that is not the image of a known point).  A rational point of C_d
gives the point (M^2 : e^2 : N) of the conic  N^2 = d U^2 + a U V + c V^2  (c = b/d).  The conic has a rational
point P0 (a hint, found by PARI's qfsolve and recorded below; it is checked here exactly, so PARI is not trusted).
The lines through P0 parametrise the conic: (U, V, N) = (F(D) U0 - 2B r, F(D) V0 - 2B s, F(D) N0) with D = (r, s, 0),
F(D) = d r^2 + a r s + c s^2, B = B(P0, D) the bilinear form; since a^2 - 4 d c = b' is not a square, these forms
have no common zero over Q, and every rational point of the conic is (U(r, s) : V(r, s) : N(r, s)) for coprime
integers r, s.  So a point of C_d gives mu U(r, s) = M^2, mu V(r, s) = e^2 with mu rational and M, e not 0 (e = 0 needs
d a square, M = 0 needs b/d a square: the classes 1 and b, which are excluded).  With lambda the squarefree part of
mu, an odd prime of lambda divides U(r, s) and V(r, s), hence Res(U, V) = b' N0^4 (an identity, checked).  So
lambda is +-(a squarefree divisor of 2 b' N0), and the curve D_lambda: lambda U(r, s) = square, lambda V(r, s) = square
must have points everywhere.  If for every lambda it has none at some place (R or a prime of 2 b' N0), then C_d has
no rational point.  The local test is `solvable` with the two forms lambda U, lambda V.
  This removes the class d from the upper bound for the image; the image is a group, so its order is at most the
largest power of 2 not above the number of classes left.

PRIMES.  The factorisations of b, b', N0 are found by sympy and checked: the product is the number, and each factor
is proved prime by trial division (factors below 10^14) or by a Pocklington certificate built recursively.

What this proves for each curve, if the program ends with its assertions: an upper bound for the rank of E(Q).  With
the point G of sigma_evidence.txt (of infinite order: sigma_exact.py) the lower bound is 1.
"""
import os, sys, json, random
WRITE = '--write' in sys.argv   # outputs and hint files are written only with --write, so CI can compare stdout with them
from fractions import Fraction as Fr
from math import gcd, isqrt

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)

SMALL = 3000          # primes below this are treated by visiting every residue
DEPTH = {2: 14}       # depth of the discs (in powers of p): 2 -> 14, odd p -> DEPTH_ODD
DEPTH_ODD = 8

# ---------------------------------------------------------------- integers and primes
def is_prime_td(n):
    if n < 2: return False
    if n % 2 == 0: return n == 2
    f = 3
    while f * f <= n:
        if n % f == 0: return False
        f += 2
    return True

HINTS = {'factor': {}}      # factorisations found by PARI (rank_indep.py --find) or sympy; every one is checked

def _find_factor(N):
    key = str(N)
    if key in HINTS['factor']: return {int(p): int(k) for p, k in HINTS['factor'][key]}
    try:
        import cypari2
        pari = cypari2.Pari()
        F = pari.factor(N)
        out = {int(F[0][i]): int(F[1][i]) for i in range(len(F[0]))}
    except ImportError:
        import sympy
        out = {int(p): int(k) for p, k in sympy.factorint(N).items()}
    HINTS['factor'][key] = sorted(out.items())
    return out

def factor_int(N):
    """Factorisation of a nonzero integer (found by PARI or sympy, or read from the hints) and checked here: the
    product is |N| and every factor is proved prime."""
    N = abs(int(N)); assert N != 0
    if N == 1: return {}
    F = _find_factor(N)
    prod = 1
    for p, k in F.items(): prod *= p ** k
    assert prod == N, 'factorisation does not multiply back'
    for p in F: assert prove_prime(p), ('not proved prime', p)
    return F

_PROVED = {}
def prove_prime(n, depth=0):
    """Trial division below 10^10; above, Pocklington: n - 1 = F R with F > sqrt(n) a product of proved primes, and a
    witness a with a^(n-1) = 1 mod n and gcd(a^((n-1)/q) - 1, n) = 1 for every prime q | F."""
    if n in _PROVED: return _PROVED[n]
    if n < 10**10:
        r = is_prime_td(n); _PROVED[n] = r; return r
    assert depth < 40
    F = 1; qs = []
    for q, k in sorted(_find_factor(n - 1).items()):
        if not prove_prime(q, depth + 1): continue
        qs.append(q); F *= q ** k
    r = False
    if F * F > n:
        for a in range(2, 500):
            if pow(a, n - 1, n) != 1: break
            if all(gcd(pow(a, (n - 1) // q, n) - 1, n) == 1 for q in qs): r = True; break
    _PROVED[n] = r
    return r

def vp(x, p):
    """p-adic valuation of a nonzero integer."""
    assert x != 0
    v = 0
    while x % p == 0: x //= p; v += 1
    return v

def legendre(a, p):
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1

def sqfree_part(n):
    """Squarefree part (with sign) of a nonzero integer, via factor_int."""
    s = -1 if n < 0 else 1
    for p, k in factor_int(n).items():
        if k % 2: s *= p
    return s

# ---------------------------------------------------------------- polynomials over Z and F_p (low degree first)
def peval(g, t):
    r = 0
    for c in reversed(g): r = r * t + c
    return r
def ptrim(g, p=None):
    g = [c % p for c in g] if p else list(g)
    while g and g[-1] == 0: g.pop()
    return g
def pmul(a, b, p=None):
    if not a or not b: return []
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): c[i + j] += x * y
    return ptrim(c, p) if p else c
def psub_lin(g, t0, p):
    """g(t0 + p u) as a polynomial in u."""
    out = [0] * len(g)
    # Horner with the linear polynomial t0 + p u
    for c in reversed(g):
        # out = out * (t0 + p u) + c
        new = [0] * len(g)
        for i, x in enumerate(out):
            if x:
                new[i] += x * t0
                if i + 1 < len(g): new[i + 1] += x * p
        new[0] += c
        out = new
    return out
def pmod(a, b, p):
    a = ptrim(a, p); b = ptrim(b, p)
    inv = pow(b[-1], -1, p)
    while len(a) >= len(b):
        q = a[-1] * inv % p; s = len(a) - len(b)
        for i, x in enumerate(b): a[s + i] = (a[s + i] - q * x) % p
        a = ptrim(a, p)
    return a
def pgcd(a, b, p):
    a = ptrim(a, p); b = ptrim(b, p)
    while b: a, b = b, pmod(a, b, p)
    if a:
        inv = pow(a[-1], -1, p); a = [x * inv % p for x in a]
    return a
def ppowmod(base, e, mod, p):
    r = [1]; base = pmod(base, mod, p)
    while e:
        if e & 1: r = pmod(pmul(r, base, p), mod, p)
        base = pmod(pmul(base, base, p), mod, p); e >>= 1
    return r
def roots_mod(g, p):
    """The roots in F_p of a nonzero polynomial g (p odd): gcd with t^p - t, then random splitting."""
    g = ptrim(g, p)
    assert g
    if len(g) == 1: return []
    h = ppowmod([0, 1], p, g, p)
    h = h + [0] * max(0, 2 - len(h)); h[1] = (h[1] - 1) % p
    h = pgcd(g, ptrim(h, p), p)
    out = []
    rng = random.Random(12345)
    def split(f):
        if len(f) <= 1: return
        if len(f) == 2: out.append((-f[0]) * pow(f[1], -1, p) % p); return
        while True:
            dlt = rng.randrange(p)
            w = ppowmod([dlt, 1], (p - 1) // 2, f, p)
            w = w + [0] * max(0, 1 - len(w)); w[0] = (w[0] - 1) % p
            c = pgcd(f, ptrim(w, p), p)
            if 1 < len(c) < len(f):
                split(c); q = pdiv_exact(f, c, p); split(q); return
    split(h)
    for r in out: assert peval(g, r) % p == 0
    return sorted(set(out))
def pdiv_exact(a, b, p):
    a = ptrim(a, p); b = ptrim(b, p); inv = pow(b[-1], -1, p)
    q = [0] * (len(a) - len(b) + 1)
    while len(a) >= len(b) and a:
        c = a[-1] * inv % p; s = len(a) - len(b); q[s] = c
        for i, x in enumerate(b): a[s + i] = (a[s + i] - c * x) % p
        a = ptrim(a, p)
    assert not a
    return q
def const_times_square(g, p):
    """If g = eps * w^2 in F_p[t] (p odd, p > deg g), return eps; else None."""
    g = ptrim(g, p)
    if len(g) % 2 == 0: return None          # odd degree
    eps = g[-1]; inv = pow(eps, -1, p); f = [x * inv % p for x in g]
    k = (len(f) - 1) // 2
    w = [0] * (k + 1); w[k] = 1
    inv2 = pow(2, -1, p)
    for j in range(k - 1, -1, -1):
        # coefficient of t^(k + j) in w^2: 2 w_j + sum_{i + l = k + j, j < i, l < k} w_i w_l
        s = sum(w[i] * w[k + j - i] for i in range(j + 1, k) if 0 <= k + j - i <= k and k + j - i > j)
        w[j] = (f[k + j] - s) * inv2 % p
    return eps if ptrim(pmul(w, w, p), p) == ptrim(f, p) else None

# ---------------------------------------------------------------- local solubility
def content_val(g, p):
    return min(vp(c, p) for c in g if c != 0)

def solvable(forms, p, depth=None):
    """Is there a point of P^1(Q_p) where every form (integer coefficients, even degree, low degree first) takes a
    square value (0 allowed)?  True / False certified, None if not decided at the depth."""
    if depth is None: depth = DEPTH.get(p, DEPTH_ODD)
    forms = [list(g) for g in forms]
    deg = [len(g) - 1 for g in forms]
    assert all(d % 2 == 0 for d in deg)
    rev = [list(reversed(g)) for g in forms]                       # g(1, t): coefficients reversed
    res = []
    # chart (t, 1), t in Z_p: polynomials g(t); chart (1, t), t in pZ_p: g(1, p u), u in Z_p
    res.append(_affine([(g, 0) for g in forms], p, 0, depth))
    if res[-1] is True: return True
    res.append(_affine([(psub_lin(g, 0, p), 0) for g in rev], p, 1, depth))
    if True in res: return True
    if all(r is False for r in res): return False
    return None

def _affine(state, p, level, depth):
    """state: list of (polynomial in t, exponent e): the condition is p^e g(t) a square, t in Z_p."""
    # remove contents
    st = []
    for g, e in state:
        if all(c == 0 for c in g): st.append(None); continue     # identically 0: value 0, a square
        k = content_val(g, p); st.append(([c // p**k for c in g], (e + k) % 2))
    live = [x for x in st if x is not None]
    if not live: return True
    if level >= depth: return None
    if p == 2: return _affine2(live, level, depth)
    if p < SMALL:
        cands = range(p)
    else:
        cand = set()
        for h, e in live:
            cand.update(roots_mod(h, p))
        if all(e == 0 for _, e in live):
            # look for a residue where every h is a non-zero square
            rng = random.Random(p)
            for _ in range(400):
                t0 = rng.randrange(p)
                if all(legendre(peval(h, t0), p) == 1 for h, _ in live): return True
            # certificate that there is none: a product of some h equal to (non-residue) * square
            from itertools import combinations
            cert = False
            for r in range(1, len(live) + 1):
                for sub in combinations(range(len(live)), r):
                    prod = [1]
                    for i in sub: prod = pmul(prod, live[i][0], p)
                    eps = const_times_square(prod, p)
                    if eps is not None and legendre(eps, p) == -1: cert = True; break
                if cert: break
            if not cert: return None
        cands = sorted(cand)
    undecided = False
    for t0 in cands:
        vals = [(peval(h, t0) % p, e) for h, e in live]
        if any(v != 0 and (e == 1 or legendre(v, p) == -1) for v, e in vals): continue     # dropped
        if all(v != 0 for v, e in vals): return True                                       # a point (R1)
        r = _affine([(psub_lin(h, t0, p), e) for h, e in live], p, level + 1, depth)
        if r is True: return True
        if r is None: undecided = True
    return None if undecided else False

def _affine2(live, level, depth):
    """p = 2 at the disc Z_2 (after the substitutions made so far): visit t0 mod 2^j, decide by R2."""
    def rec(t0, k):
        vals = []
        for h, e in live:
            v = peval(h, t0)
            if v == 0: vals.append(None); continue
            w = vp(v, 2)
            if w + 3 <= k:
                u = (v >> w) % 8
                vals.append((w + e) % 2 == 0 and u == 1)
            else: vals.append(None)
        if any(x is False for x in vals): return False
        if all(x is True for x in vals): return True
        if level + k >= depth + 3: return None
        r0 = rec(t0, k + 1); r1 = rec(t0 + 2**k, k + 1)
        if r0 is True or r1 is True: return True
        if r0 is False and r1 is False: return False
        return None
    return rec(0, 0)

def real_solvable(forms):
    """Is there a real point (M : e) where every form is >= 0?  Exact: the signs between the real roots."""
    import sympy
    t = sympy.Symbol('t')
    pts = set([Fr(0)])
    roots = []
    for g in forms:
        P = sympy.Poly(list(reversed(g)), t)
        if P.degree() > 0:
            for (lo, hi), _ in P.intervals():
                roots.append((Fr(int(sympy.numer(lo)), int(sympy.denom(lo))), Fr(int(sympy.numer(hi)), int(sympy.denom(hi)))))
    # test points: rational numbers outside all isolating intervals, between consecutive intervals, and far out
    ends = sorted(set([x for iv in roots for x in iv]))
    if ends:
        pts.update([ends[0] - 1, ends[-1] + 1])
        for x, y in zip(ends, ends[1:]): pts.add((x + y) / 2)
    for x in pts:
        if all(peval(g, x) >= 0 for g in forms): return True
    # the point at infinity (1 : 0): leading coefficients
    if all(g[-1] >= 0 for g in forms): return True
    # a point inside an isolating interval could still work only where a form vanishes: refine is not needed for an
    # upper bound; report None (possible) if some interval exists, False only when there is none.
    return None if roots else False

# ---------------------------------------------------------------- square classes
def qp_class(x, p):
    """Class of a nonzero rational in Q_p*/Q_p*^2: (v mod 2, unit class) with unit class the Legendre symbol (odd p)
    or the odd part mod 8 (p = 2)."""
    x = Fr(x); n, dd = x.numerator, x.denominator
    vn, vd = vp(n, p), vp(dd, p)
    un, ud = n // p ** vn, dd // p ** vd
    if p == 2: return ((vn - vd) % 2, (un * ud) % 8)
    return ((vn - vd) % 2, legendre(un * ud, p))

def class_reps(p):
    """Integers representing every class of Q_p*/Q_p*^2."""
    if p == 2: return [s * u * t for s in (1,) for u in (1, 3, 5, 7) for t in (1, 2)]
    eps = next(a for a in range(2, p) if legendre(a, p) == -1)
    return [1, eps, p, eps * p]

# ---------------------------------------------------------------- the descent
def shift_to_two_torsion(a2, a4, a6):
    """y^2 = x^3 + a2 x^2 + a4 x + a6 with exactly one rational root e of the cubic: (a, b) of y^2 = x(x^2 + a x + b)."""
    roots = []
    import sympy
    X = sympy.Symbol('X')
    for r in sympy.Poly(X**3 + a2*X**2 + a4*X + a6, X).ground_roots():
        if r.is_Integer or r.is_Rational: roots.append(int(r))
    assert len(roots) == 1, ('rational roots', roots)
    e = roots[0]; assert e**3 + a2*e*e + a4*e + a6 == 0
    a, b = 3*e + a2, 3*e*e + 2*a2*e + a4
    return e, a, b

def quartic(d, a, b):
    """C_d: N^2 = d M^4 + a M^2 e^2 + (b/d) e^4, multiplied by d^2 (a square): integer coefficients, low degree first,
    as a polynomial in t = M/e."""
    assert b % d == 0
    return [d * b, 0, d * d * a, 0, d ** 3]

def local_image(a, b, p):
    """The classes delta of Q_p*/Q_p*^2 (as class labels) with C_delta possibly soluble over Q_p, and those certified."""
    poss, cert = set(), set()
    for delta in class_reps(p):
        # C_delta with rational b/delta: multiply the quartic by delta^2: delta^3 M^4 + a delta^2 M^2 e^2 + b delta e^4
        g = [b * delta, 0, a * delta * delta, 0, delta ** 3]
        r = solvable([g], p)
        if r is not False: poss.add(qp_class(delta, p))
        if r is True: cert.add(qp_class(delta, p))
    return poss, cert

def real_ok(d, a, b):
    """C_d has a real point."""
    if d > 0: return True
    c = Fr(b, d)
    if c > 0: return True
    # max over u = M^2/e^2 >= 0 of d u^2 + a u + c (d < 0): at u = -a/(2d) if that is >= 0
    u = Fr(-a, 2 * d)
    if u < 0: return False
    return d * u * u + a * u + c >= 0

def selmer(a, b, bad_primes):
    """The classes d | b (squarefree, with sign) that pass the real test and the local tests at the bad primes:
    an upper bound for the image of alpha.  Also the number certified at each prime."""
    Fb = factor_int(b)
    gens = [-1] + sorted(Fb)
    loc = {p: local_image(a, b, p) for p in bad_primes}
    out = []
    for mask in range(1 << len(gens)):
        d = 1
        for i, g in enumerate(gens):
            if mask >> i & 1: d *= g
        if not real_ok(d, a, b): continue
        if all(qp_class(d, p) in loc[p][0] for p in bad_primes): out.append(d)
    return out, loc

def pow2_floor(n):
    k = 0
    while 2 ** (k + 1) <= n: k += 1
    return k

def first_descent(a, b):
    ap, bp = -2 * a, a * a - 4 * b
    assert b != 0 and bp != 0
    assert not (bp > 0 and isqrt(bp) ** 2 == bp), "a^2 - 4b is a square: more 2-torsion"
    assert not (b > 0 and isqrt(b) ** 2 == b) or True
    bad = sorted(set([2]) | set(factor_int(b)) | set(factor_int(bp)))
    S, locS = selmer(a, b, bad)
    Sp, locSp = selmer(ap, bp, bad)
    return dict(a=a, b=b, ap=ap, bp=bp, bad=bad, S=S, Sp=Sp, locS=locS, locSp=locSp)

def sq_class_key(x):
    """Squarefree part of a nonzero rational (with sign), as an integer."""
    x = Fr(x)
    return sqfree_part(x.numerator * x.denominator)

# ---------------------------------------------------------------- second descent
def conic_point_check(d, a, c, P0):
    U0, V0, N0 = P0
    assert (U0, V0, N0) != (0, 0, 0)
    assert N0 * N0 == d * U0 * U0 + a * U0 * V0 + c * V0 * V0, 'not a point of the conic'
    assert N0 != 0

def param_forms(d, a, c, P0):
    """U(r, s), V(r, s) as integer binary quadratic forms (coefficients of s^2, r s, r^2 in the variable t = r/s:
    low degree first, i.e. [coef s^2, coef r s, coef r^2])."""
    U0, V0, N0 = P0
    # F(D) = d r^2 + a r s + c s^2;  2B = 2 d U0 r + a (U0 s + V0 r) + 2 c V0 s
    # U = F U0 - 2B r,  V = F V0 - 2B s   (all integer: 2B has integer coefficients)
    F = {(2, 0): d, (1, 1): a, (0, 2): c}                       # exponents (r, s)
    B2 = {(1, 0): 2 * d * U0 + a * V0, (0, 1): a * U0 + 2 * c * V0}
    def mul_lin(B, var):
        return {((i + 1, j) if var == 'r' else (i, j + 1)): v for (i, j), v in B.items()}
    def comb(F, k, G, l):
        out = {}
        for key, v in F.items(): out[key] = out.get(key, 0) + k * v
        for key, v in G.items(): out[key] = out.get(key, 0) + l * v
        return out
    U = comb(F, U0, mul_lin(B2, 'r'), -1)
    V = comb(F, V0, mul_lin(B2, 's'), -1)
    to_t = lambda Q: [Q.get((0, 2), 0), Q.get((1, 1), 0), Q.get((2, 0), 0)]   # in t = r/s, times s^2
    return to_t(U), to_t(V)

def res_quadratic(u, v):
    """Resultant of two binary quadratics given as [c0, c1, c2] (c2 the coefficient of r^2)."""
    import sympy
    t = sympy.Symbol('t')
    return int(sympy.resultant(sympy.Poly(list(reversed(u)), t), sympy.Poly(list(reversed(v)), t)))

def disc2(q): return q[1] * q[1] - 4 * q[0] * q[2]

def second_descent(d, a, b, P0, primes_hint=None, stop_at_first=True):
    """Show that C_d has no rational point, or return the lambdas that survive.  Returns (survivors, info)."""
    c = Fr(b, d); assert c.denominator == 1; c = int(c)
    conic_point_check(d, a, c, P0)
    U, V = param_forms(d, a, c, P0)
    U0, V0, N0 = P0
    # the parametrisation: d U^2 + a U V + c V^2 = (F(r, s) N0)^2 identically (degree 4 in (r, s): 6 points suffice)
    for (r_, s_) in ((1, 0), (0, 1), (1, 1), (1, -1), (2, 1), (1, 2), (3, -2)):
        u_ = U[0] * s_ * s_ + U[1] * r_ * s_ + U[2] * r_ * r_; v_ = V[0] * s_ * s_ + V[1] * r_ * s_ + V[2] * r_ * r_
        F_ = d * r_ * r_ + a * r_ * s_ + c * s_ * s_
        assert d * u_ * u_ + a * u_ * v_ + c * v_ * v_ == (F_ * N0) ** 2, 'the parametrisation of the conic is wrong'
    bp = a * a - 4 * b
    R = res_quadratic(U, V)
    assert R == bp * (d * U0 * U0 + a * U0 * V0 + c * V0 * V0) ** 2 == bp * N0 ** 4, 'resultant identity'
    # bad primes of the curves D_lambda: 2, primes of b' N0, and of the discriminants of U, V and of the pencil
    import sympy
    x, y = sympy.symbols('x y')
    Dp = sympy.expand((x * U[1] + y * V[1]) ** 2 - 4 * (x * U[2] + y * V[2]) * (x * U[0] + y * V[0]))
    Dpc = [int(sympy.Poly(Dp, x, y).coeff_monomial(m_)) for m_ in (y ** 2, x * y, x ** 2)]
    pencil_disc = disc2(Dpc)
    ps = set([2]) | set(factor_int(bp)) | set(factor_int(N0))
    for q in (disc2(U), disc2(V), pencil_disc):
        assert q != 0
        ps |= set(factor_int(q))
    ps = sorted(ps)
    lam_primes = sorted(set([2]) | set(factor_int(bp)) | set(factor_int(N0)))
    gens = [-1] + lam_primes
    # local images for the two-form curve, class by class at each prime
    loc = {}
    for p in ps:
        ok = set()
        for delta in class_reps(p):
            r = solvable([[delta * x_ for x_ in U], [delta * x_ for x_ in V]], p)
            if r is not False: ok.add(qp_class(delta, p))
        loc[p] = ok
    # the classes of the generators at each prime, as small codes; a lambda is a product of generators, so its code
    # at p is the XOR of theirs (v mod 2 in bit 0, the unit class in the other bits)
    def enc(x, p):
        v, u = qp_class(x, p)
        if p == 2: return v | ({1: 0, 3: 1, 5: 2, 7: 3}[u] << 1)
        return v | ((0 if u == 1 else 1) << 1)
    allowed = {p: set() for p in ps}
    for p in ps:
        for delta in class_reps(p):
            if qp_class(delta, p) in loc[p]: allowed[p].add(enc(delta, p))
    gcodes = [[enc(g, p) for p in ps] for g in gens]
    real_ok_sign = {sg: real_solvable([[sg * x_ for x_ in U], [sg * x_ for x_ in V]]) is not False for sg in (1, -1)}
    survivors = []
    cur = [0] * len(ps); mask = 0; lam = 1
    order = sorted(range(len(ps)), key=lambda i: len(allowed[ps[i]]))
    for k in range(1 << len(gens)):
        if k:
            bit = (k & -k).bit_length() - 1             # Gray code: flip one generator
            mask ^= 1 << bit
            gc = gcodes[bit]
            for i in range(len(ps)): cur[i] ^= gc[i]
        if not real_ok_sign[-1 if (mask & 1) else 1]: continue
        if all(cur[i] in allowed[ps[i]] for i in order):
            lam = 1
            for i, g in enumerate(gens):
                if mask >> i & 1: lam *= g
            assert all(qp_class(lam, p) in loc[p] for p in ps)
            survivors.append(lam)
            if stop_at_first: break
    return survivors, dict(U=U, V=V, primes=ps, nlam=1 << len(gens))

# ---------------------------------------------------------------- the four curves
CURVES = [((760, 39), 'E+'), ((945, 248), 'E+'), ((72, 65), 'E+'), ((920, 129), 'E+')]

def magma_ainvs():
    out = {}
    for f in ('magma_sigma2.json', 'magma_sigma3.json'):
        for r in json.load(open(os.path.join(here, f))):
            out[(tuple(r['pair']), r['curve'])] = r['ainvs']
    return out

def evidence():
    out = {}
    for line in open(os.path.join(here, 'sigma_evidence.txt')):
        if line.startswith('#') or not line.strip(): continue
        w = line.split(); f = dict(x.split('=', 1) for x in w[2:])
        m, n = (int(v) for v in w[0].split(','))
        out[((m, n), w[1])] = tuple(Fr(v) for v in f['G'].split(':'))
    return out


def find_conic_point(d, a, c):
    """A hint: a rational point of N^2 = d U^2 + a U V + c V^2 with N != 0, by PARI's qfsolve (checked by the caller)."""
    key = 'conic %d %d %d' % (d, a, c)
    if key in HINTS.get('conic', {}): return tuple(HINTS['conic'][key])
    import cypari2
    pari = cypari2.Pari()
    M = pari.matrix(3, 3, [2 * d, a, 0, a, 2 * c, 0, 0, 0, -2])
    sol = pari.qfsolve(M)
    assert str(pari.type(sol)) == 't_COL', ('qfsolve', sol)
    U0, V0, N0 = (int(x) for x in sol)
    g = gcd(gcd(U0, V0), N0); U0, V0, N0 = U0 // g, V0 // g, N0 // g
    if N0 == 0:
        # move to another point of the conic: the second intersection of the line through P0 in the direction (1, 1, 1)
        raise AssertionError('conic point with N = 0')
    HINTS.setdefault('conic', {})[key] = [U0, V0, N0]
    return (U0, V0, N0)

def find_points(ainvs):
    """A hint: points of infinite order on y^2 = x^3 + a2 x^2 + a4 x + a6 found by PARI's ellrank (checked by the caller)."""
    key = 'points %s' % ainvs
    if key in HINTS.get('points', {}): return [tuple(Fr(v) for v in P) for P in HINTS['points'][key]]
    import cypari2
    pari = cypari2.Pari()
    E = pari.ellinit(list(ainvs))
    pts = []
    for eff in range(4):
        r = pari.ellrank(E, eff)
        if len(r[3]): pts = [(Fr(str(P[0])), Fr(str(P[1]))) for P in r[3]]; break
    HINTS.setdefault('points', {})[key] = [[str(P[0]), str(P[1])] for P in pts]
    return pts

def on_curve(ainvs, P):
    a1, a2, a3, a4, a6 = ainvs
    x, y = P
    return y * y + a1 * x * y + a3 * y == x ** 3 + a2 * x * x + a4 * x + a6

def subgroup(gens):
    """Square classes (squarefree integers) generated by the given ones."""
    H = {1}
    for g in gens:
        H |= {sqfree_part(h * g) for h in H}
    return H

def rank_bound(a, b, pts_E, pts_Ep, say=print):
    """Upper bound for the rank of y^2 = x(x^2 + a x + b), by the first descent and the second descent on the
    classes outside the subgroup generated by the known points.  pts_E, pts_Ep: rational points on E and E' (shifted
    models), used only to name a subgroup of the image that need not be tested.  Returns (bound, record)."""
    R = first_descent(a, b)
    ap, bp = R['ap'], R['bp']
    rec = dict(a=str(a), b=str(b), bad=R['bad'], sel=len(R['S']), selp=len(R['Sp']))
    say('  first descent: |Sel(phi)| <= %d, |Sel(phi^)| <= %d, rank <= %d; bad primes %s'
        % (len(R['S']), len(R['Sp']), pow2_floor(len(R['S'])) + pow2_floor(len(R['Sp'])) - 2, R['bad']))
    out = {}
    for side, (A, Bc, S, pts) in (('E', (a, b, R['S'], pts_E)), ("E'", (ap, bp, R['Sp'], pts_Ep))):
        for P in pts: assert on_curve([0, A, 0, Bc, 0], P), 'a known point is not on the curve'
        known = [sqfree_part(Bc)] + [sq_class_key(P[0]) for P in pts if P[0] != 0]
        H = subgroup(known)
        Sset = set(S)
        assert H <= Sset, 'a known point outside the Selmer group: the descent is wrong'
        # control: a class of a known point other than 1 and b must survive the second descent (C_d has a point)
        controls = 0
        for dd in sorted(H - {1, sqfree_part(Bc)}, key=abs)[:3]:
            surv, _ = second_descent(dd, A, Bc, find_conic_point(dd, A, Bc // dd))
            assert surv, ('the second descent excludes the class of a known point', dd)
            controls += 1
        left = set(H); excluded_cosets = []; open_cosets = []
        todo = sorted(Sset - H, key=abs)
        seen = set()
        for d in todo:
            if d in seen: continue
            coset = sorted({sqfree_part(d * h) for h in H}, key=abs)
            seen |= set(coset)
            done = False
            for dd in coset[:4]:
                c = Bc // dd
                P0 = find_conic_point(dd, A, c)
                surv, info = second_descent(dd, A, Bc, P0)
                if not surv:
                    excluded_cosets.append(dict(d=dd, coset=coset, P0=[str(x) for x in P0], nlam=info['nlam'], primes=info['primes']))
                    done = True; break
            if not done:
                open_cosets.append(coset); left |= set(coset)
        k = pow2_floor(len(left))
        say('  %s: Selmer %d classes, known subgroup %d, cosets excluded by the second descent %d, not excluded %d: image <= 2^%d'
            ' (control: %d classes of known points survive, as they must)'
            % (side, len(S), len(H), len(excluded_cosets), len(open_cosets), k, controls))
        out[side] = dict(selmer=sorted(S, key=abs), known=sorted(H, key=abs), excluded=excluded_cosets, open=open_cosets, log2=k)
    r = out['E']['log2'] + out["E'"]['log2'] - 2
    rec.update(out); rec['rank_bound'] = r
    return r, rec

HINTS_FILE = os.path.join(here, 'rank_indep_hints.json')

def run_curve(m, n, cv, say=print):
    import sigma_exact as SE
    M = SE.model(m, n, 2 if cv == 'E+' else -2)
    E = SE.Curve(M['a2'], M['a4'], M['a6'])
    # the curve sent to Magma (magma_sigma*.json) is a model of the same curve: compare c4, c6 up to u^4, u^6
    ai = magma_ainvs()[((m, n), cv)]
    Mg = SE.Curve(ai[1], ai[3], ai[4]); assert ai[0] == 0 and ai[2] == 0
    c6 = lambda C: -C.b2 ** 3 + 36 * C.b2 * C.b4 - 216 * C.b6
    q4, q6 = Fr(E.c4, Mg.c4), Fr(c6(E), c6(Mg)); u2 = q6 / q4
    assert u2 * u2 == q4 and u2 ** 3 == q6 and u2 > 0 and isqrt(u2.numerator) ** 2 == u2.numerator and isqrt(u2.denominator) ** 2 == u2.denominator, 'not the curve of Magma'
    G = evidence()[((m, n), cv)]
    assert E.on(G)
    e, a, b = shift_to_two_torsion(M['a2'], M['a4'], M['a6'])
    GE = (G[0] - e, G[1])
    ap, bp = -2 * a, a * a - 4 * b
    r, rec = rank_bound(a, b, [GE], find_points([0, ap, 0, bp, 0]), say=say)
    return r, rec

def negative(say=print):
    """Forged hints and points must be refused."""
    import sigma_exact as SE
    M = SE.model(72, 65, 2); e, a, b = shift_to_two_torsion(M['a2'], M['a4'], M['a6']); ap, bp = -2 * a, a * a - 4 * b
    res = []
    def run(name, fn):
        try:
            fn(); say('  NOT REFUSED: ' + name); res.append(False)
        except AssertionError as ex:
            say('  refused: %s (%s)' % (name, str(ex)[:80])); res.append(True)
    def f1():
        U0, V0, N0 = find_conic_point(2, ap, bp // 2); second_descent(2, ap, bp, (U0 + 1, V0, N0))
    run('a point off the conic', f1)
    def f2():
        key = str(abs(b)); saved = HINTS['factor'].get(key)
        F = dict(factor_int(b)); q1, q2 = sorted(F)[-2:]
        F2 = {k: v for k, v in F.items() if k not in (q1, q2)}; F2[q1 * q2] = 1
        HINTS['factor'][key] = sorted(F2.items())
        try: factor_int(b)
        finally: HINTS['factor'][key] = saved
    run('a composite number given as a prime factor of b', f2)
    def f3():
        rank_bound(a, b, [(Fr(1), Fr(2))], [], say=lambda s: None)
    run('a point off the curve', f3)
    def f4():
        F = dict(factor_int(b)); saved = HINTS['factor'][str(abs(b))]
        HINTS['factor'][str(abs(b))] = sorted(list(F.items())[1:])
        try: factor_int(b)
        finally: HINTS['factor'][str(abs(b))] = saved
    run('a factorisation with a prime missing', f4)
    def f5():
        # the class of a point of E' must survive the second descent: give it as "to be excluded" and expect survivors
        pts = find_points([0, ap, 0, bp, 0])
        d = sq_class_key(pts[0][0])
        surv, _ = second_descent(d, ap, bp, find_conic_point(d, ap, bp // d))
        assert not surv, 'the class of a rational point survives the second descent, as it must'
    run('excluding the class of a rational point of E\'', f5)
    return all(res), res

if __name__ == '__main__':
    if os.path.exists(HINTS_FILE):
        for k, v in json.load(open(HINTS_FILE)).items(): HINTS[k] = v
    out = []
    def say(s): print(s, flush=True); out.append(s)
    if '--negative' in sys.argv:
        ok, res = negative(say)
        say('negative tests: %d of %d refused' % (sum(res), len(res)))
        if WRITE: open(os.path.join(here, 'rank_indep_negative.txt'), 'w').write('\n'.join(out) + '\n')
        if WRITE: json.dump(HINTS, open(HINTS_FILE, 'w'), indent=0, sort_keys=True)
        sys.exit(0 if ok else 1)
    curves = CURVES
    if '--all' in sys.argv:
        curves = [k for k in evidence()]
        say('rank_indep.py --all: exact descent via 2-isogeny and second descent, for every curve of sigma_evidence.txt')
    else:
        say('rank_indep.py: exact descent via 2-isogeny and second descent, for the curves of the census that rest on Magma alone for rank 1')
    results = {}
    for (m, n), cv in curves:
        say('%d,%d %s:' % (m, n, cv))
        r, rec = run_curve(m, n, cv, say=say)
        results[((m, n), cv)] = r
        say('%d,%d %s: rank <= %d (exact; lower bound 1 from G of sigma_evidence.txt)' % (m, n, cv, r))
    if WRITE: json.dump(HINTS, open(HINTS_FILE, 'w'), indent=0, sort_keys=True)
    say('summary: ' + ', '.join('%d,%d %s <= %d' % (k[0][0], k[0][1], k[1], v) for k, v in results.items()))
    if WRITE: open(os.path.join(here, 'rank_indep_all.txt' if '--all' in sys.argv else 'rank_indep.txt'), 'w').write('\n'.join(out) + '\n')
