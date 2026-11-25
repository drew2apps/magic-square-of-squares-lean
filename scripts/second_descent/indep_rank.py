#!/usr/bin/env python3
"""
indep_rank.py -- rank 0 without Magma for the curves of magma_second_ranks.json,
by an exact descent via 2-isogeny written here (Silverman-Tate, Rational Points
on Elliptic Curves, ch. III.6; Silverman, AEC X.4), and exact torsion over Q and
over Q(i).

Method.  E: y^2 = x^3 + a2 x^2 + a4 x + a6 (integral, a1 = a3 = 0) with a rational
root x0 of the cubic.  Translate x -> x + x0:  E: y^2 = x(x^2 + A x + B), and
E': y^2 = x(x^2 - 2A x + (A^2 - 4B)).  The maps alpha: E(Q) -> Q*/Q*^2,
(x, y) -> x, (0, 0) -> B, and alpha': E'(Q) -> Q*/Q*^2 satisfy

        2^r = #alpha(E(Q)) * #alpha'(E'(Q)) / 4          (r = rank E(Q) = rank E'(Q)).

alpha(E(Q)) lies in the Selmer set Sel = { d squarefree, d | B : the quartic
N^2 = d M^4 + A M^2 e^2 + (B/d) e^4 has a point over R and over every Q_p }, and
likewise for E'.  Only p | 2 B (A^2 - 4B) need testing (elsewhere the curve has
good reduction, and a genus-1 curve over F_p has a point, which lifts).  So

        r <= log2( #Sel(E) * #Sel(E') / 4 ).

The local test at p (function `qp_solvable`) is an exact search over P^1(Z_p),
classes x0 mod p^k: g(x0 + p^k t) = sum c_i t^i, and if v(c_0) < min_{i>=1} v(c_i)
(odd p) or v(c_0) + 3 <= min v(c_i) (p = 2), every value in the class has the
square class of c_0, which decides the class; otherwise refine.  Each exclusion of
a d is a certificate (d, p) that the program re-verifies; the search depth is
capped and an undecided class counts as "has a point" (the safe side).

Torsion: #E(Q)_tors divides #E(F_p) for every odd prime of good reduction
(reduction is injective on torsion there), and #E(K)_tors divides #E(F_p) for
split p = 1 mod 4 and #E(F_{q^2}) for inert q, q odd of good reduction.  The
point counts are computed here by Legendre sums (no PARI).  Lower bounds come from
explicit points (2-torsion from the roots, 4-torsion by halving).

PARI (cypari2) is used only to factor integers and polynomials and to find roots modulo primes, each result checked (each factor proved prime with
isprime, and the product checked), and, as a separate comparison, for its own
ellrank / elltors / ellglobalred.

Usage:  python3 indep_rank.py magma_second_ranks.json      (writes to stdout; indep_rank.txt)
        python3 indep_rank.py --ainvs a1,a2,a3,a4,a6 [m,n]  (one curve; indep_rank_440_279.txt)
        python3 indep_rank.py --selftest                    (indep_rank_selftest.txt)
"""
import sys, json
from fractions import Fraction
from math import isqrt, gcd
import cypari2

pari = cypari2.Pari()
pari.allocatemem(2 * 10 ** 9)

# ------------------------------------------------------------- integers
def factor_proved(n):
    """Prime factorisation of |n| > 0 with every prime proved by PARI's isprime."""
    n = abs(int(n))
    assert n > 0
    F = pari.factor(n)
    out = {}
    prod = 1
    for p, e in zip(F[0], F[1]):
        p, e = int(p), int(e)
        assert bool(pari.isprime(p)), p
        out[p] = e
        prod *= p ** e
    assert prod == n
    return out

def vp(n, p):
    if n == 0:
        return 10 ** 9
    k = 0
    while n % p == 0:
        n //= p
        k += 1
    return k

def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1

def is_square_Qp_unit_class(c, p):
    """c a nonzero integer: is c a square in Q_p?"""
    v = vp(c, p)
    if v % 2:
        return False
    u = c // p ** v
    if p == 2:
        return u % 8 == 1
    return legendre(u, p) == 1

def is_rat_square(q):
    q = Fraction(q)
    if q < 0:
        return False
    a, b = q.numerator, q.denominator
    return isqrt(a) ** 2 == a and isqrt(b) ** 2 == b

# ------------------------------------------------------------- local solvability
def poly_shift(coeffs, x0, h):
    """coeffs of g(x) = sum coeffs[i] x^i; return coefficients c_i of g(x0 + h t) in t."""
    n = len(coeffs) - 1
    # Taylor expansion via repeated synthetic division
    c = list(coeffs)
    out = []
    for i in range(n + 1):
        # evaluate and divide by (x - x0)
        acc = 0
        newc = [0] * (len(c) - 1) if len(c) > 1 else []
        for j in range(len(c) - 1, -1, -1):
            acc = acc * x0 + c[j]
            if j > 0:
                newc[j - 1] = acc
        out.append(acc * h ** i)
        c = newc
    return out

LIMIT = 3000   # below this, residues mod p are enumerated; above, the structure of the reduction is used

def _pmod(f, p):
    f = [x % p for x in f]
    while f and f[-1] == 0:
        f.pop()
    return f

def _pmul(f, g, p):
    if not f or not g:
        return []
    out = [0] * (len(f) + len(g) - 1)
    for i, a in enumerate(f):
        for j, b in enumerate(g):
            out[i + j] = (out[i + j] + a * b) % p
    return _pmod(out, p)

def _pdivmod(f, g, p):
    f = list(f)
    q = [0] * max(len(f) - len(g) + 1, 1)
    inv = pow(g[-1], -1, p)
    while len(f) >= len(g) and f:
        c = f[-1] * inv % p
        s = len(f) - len(g)
        q[s] = c
        for i, b in enumerate(g):
            f[s + i] = (f[s + i] - c * b) % p
        f = _pmod(f, p)
    return _pmod(q, p), f

def _pgcd(f, g, p):
    f, g = _pmod(f, p), _pmod(g, p)
    while g:
        f, g = g, _pdivmod(f, g, p)[1]
    if f:
        inv = pow(f[-1], -1, p)
        f = [x * inv % p for x in f]
    return f

def _ppowmod(base, e, mod, p):
    r = [1]
    b = _pdivmod(base, mod, p)[1]
    while e:
        if e & 1:
            r = _pdivmod(_pmul(r, b, p), mod, p)[1]
        b = _pdivmod(_pmul(b, b, p), mod, p)[1]
        e >>= 1
    return r

def roots_mod_p(f, p):
    """All roots in F_p of f (nonzero mod p, increasing-degree coefficients), exactly:
    enumeration for small p; for large p, the number of roots is deg gcd(f, t^p - t)
    (computed here) and the roots are PARI's polrootsmod, each checked by substitution."""
    f = _pmod(f, p)
    assert f
    if len(f) == 1:
        return []
    if p < LIMIT:
        return [t for t in range(p) if sum(c * pow(t, i, p) for i, c in enumerate(f)) % p == 0]
    tp = _ppowmod([0, 1], p, f, p)
    tp_minus_t = _pmod([(tp[i] if i < len(tp) else 0) - (1 if i == 1 else 0) for i in range(max(len(tp), 2))], p)
    g = _pgcd(f, tp_minus_t, p) if tp_minus_t else _pgcd(f, f, p)
    nroots = len(g) - 1
    if nroots == 0:
        return []
    P = pari.Pol(list(reversed(f)))
    rts = sorted(set(int(pari.lift(r)) for r in pari.polrootsmod(P, p)))
    for t in rts:
        assert sum(c * pow(t, i, p) for i, c in enumerate(f)) % p == 0
    assert len(rts) == nroots, (len(rts), nroots)
    return rts

def _poly_sqrt_monic(f, p):
    """If the monic f (even degree) is the square of a polynomial over F_p, return it, else None."""
    n = len(f) - 1
    if n % 2:
        return None
    k = n // 2
    h = [0] * (k + 1)
    h[k] = 1
    inv2 = pow(2, -1, p)
    # match coefficients from the top: f_{n-j} = sum_{a+b = 2k-j} h_a h_b
    for j in range(1, k + 1):
        s = sum(h[a] * h[2 * k - j - a] for a in range(k - j + 1, k) if 0 <= 2 * k - j - a <= k)
        h[k - j] = (f[n - j] - s) * inv2 % p
    return h if _pmul(h, h, p) == _pmod(f, p) else None

def class_search(coeffs, p, start_classes, maxdepth):
    """Is there x in the given union of classes of Z_p with g(x) a square in Q_p
    (g(x) = 0 counts)?  start_classes: list of (x0, k): x ranges over x0 + p^k Z_p.
    Returns False only when proved (every class decided: no square); True when a
    square value is found, or (large p only) when the reduction is not a constant
    times a square, where Weil's bound gives one (this answer is on the safe side
    and is not used as a proof); None when the depth cap is reached."""
    stack = list(start_classes)
    undecided = False
    while stack:
        x0, k = stack.pop()
        c = poly_shift(coeffs, x0, p ** k)
        if c[0] == 0:
            return True
        if p == 2:
            v0 = vp(c[0], p)
            w = min(vp(ci, p) for ci in c[1:]) if any(c[1:]) else 10 ** 9
            if v0 + 3 <= w:
                if is_square_Qp_unit_class(c[0], p):
                    return True
                continue
            if k >= maxdepth:
                undecided = True
                continue
            for d in range(2):
                stack.append((x0 + d * 2 ** k, k + 1))
            continue
        mu = min(vp(ci, p) for ci in c if ci)
        H = _pmod([ci // p ** mu for ci in c], p)
        # residues t with H(t) != 0 mod p: the class of g is p^mu * H(t)
        if mu % 2 == 0:
            if p < LIMIT:
                if any(legendre(sum(h * pow(t, i, p) for i, h in enumerate(H)), p) == 1 for t in range(p)):
                    return True
                rts = roots_mod_p(H, p)
            else:
                lead = H[-1]
                monic = [h * pow(lead, -1, p) % p for h in H]
                sq = _poly_sqrt_monic(monic, p)
                if sq is None or legendre(lead, p) == 1:
                    return True          # safe side: not a proof, and not needed as one
                rts = roots_mod_p(sq, p)  # H = lead * sq^2 with lead a non-residue
        else:
            rts = roots_mod_p(H, p)
        if k >= maxdepth:
            if rts:
                undecided = True
            continue
        for t in rts:
            stack.append((x0 + t * p ** k, k + 1))
    return None if undecided else False

def qp_solvable(q4, p, maxdepth=40):
    """q4 = (a, b, c, d, e): binary quartic a X^4 + b X^3 Z + c X^2 Z^2 + d X Z^3 + e Z^4.
    Does Y^2 = q4(X, Z) have a point over Q_p with (X:Z) in P^1(Q_p)?
    Charts: (x, 1), x in Z_p, and (1, z), z in p Z_p (every point of P^1(Q_p) is one of
    these up to a scalar lambda, and q4 scales by lambda^4, a square)."""
    a, b, c, d, e = q4
    g1 = [e, d, c, b, a]          # g(x) = q4(x, 1), coefficients in increasing degree
    g2 = [a, b, c, d, e]          # q4(1, z)
    r1 = class_search(g1, p, [(0, 0)], maxdepth)
    if r1 is True:
        return True
    r2 = class_search(g2, p, [(0, 1)], maxdepth)
    if r2 is True:
        return True
    if r1 is None or r2 is None:
        return None
    return False

def real_solvable(d, A, Bd):
    """Is d X^4 + A X^2 + Bd >= 0 for some real (X : 1) or at X = infinity?"""
    if d > 0 or Bd >= 0:
        return True
    # d < 0, Bd < 0: max over X^2 = s >= 0 of d s^2 + A s + Bd, at s = -A/(2d) if A > 0
    if A <= 0:
        return False
    return A * A - 4 * d * Bd >= 0

# ------------------------------------------------------------- 2-isogeny descent
def squarefree_divisors(primes):
    out = [1]
    for p in primes:
        out = out + [x * p for x in out]
    return out + [-x for x in out]

def selmer_phi(A, B, bad, maxdepth=40):
    """The Selmer set of the isogeny with kernel (0,0) on y^2 = x(x^2 + A x + B)."""
    primesB = sorted(factor_proved(B))
    sel = []
    certs = {}
    for d in squarefree_divisors(primesB):
        Bd = B // d
        assert Bd * d == B
        q4 = (d, 0, A, 0, Bd)
        if not real_solvable(d, A, Bd):
            certs[d] = "inf"
            continue
        excluded = None
        for p in bad:
            r = qp_solvable(q4, p, maxdepth)
            if r is False:
                excluded = p
                break
        if excluded is None:
            sel.append(d)
        else:
            certs[d] = excluded
    return sel, certs

def verify_cert(A, B, d, p, maxdepth=40):
    """Re-verify the claim 'N^2 = d M^4 + A M^2 e^2 + (B/d) e^4 has no point over Q_p'."""
    Bd = B // d
    assert Bd * d == B
    if p == "inf":
        return not real_solvable(d, A, Bd)
    return qp_solvable((d, 0, A, 0, Bd), p, maxdepth) is False

def isogeny_bound(A, B, maxdepth=40):
    A1, B1 = -2 * A, A * A - 4 * B
    bad = sorted(set(factor_proved(2 * B * B1)))
    s, cs = selmer_phi(A, B, bad, maxdepth)
    s1, cs1 = selmer_phi(A1, B1, bad, maxdepth)
    for d, p in cs.items():
        assert verify_cert(A, B, d, p, maxdepth)
    for d, p in cs1.items():
        assert verify_cert(A1, B1, d, p, maxdepth)
    prod = len(s) * len(s1)
    assert prod >= 4 and prod & (prod - 1) == 0, (len(s), len(s1))
    r = prod.bit_length() - 1 - 2
    return r, s, s1, cs, cs1, bad

# ------------------------------------------------------------- the curves
def rational_roots_cubic(a2, a4, a6):
    """Integer roots of x^3 + a2 x^2 + a4 x + a6 (monic integral: rational roots are integers)."""
    roots = []
    if a6 == 0:
        roots.append(0)
    P = pari.Pol([1, a2, a4, a6])
    for r in pari.polrootsmod(P, 10007) if False else []:
        pass
    # exact: factor over Q with PARI, then check each linear factor by substitution
    F = pari.factor(P)
    for f in F[0]:
        if int(pari.poldegree(f)) == 1:
            c1, c0 = int(pari.polcoef(f, 1)), int(pari.polcoef(f, 0))
            assert c1 == 1 or c1 == -1
            x0 = -c0 * c1
            assert x0 ** 3 + a2 * x0 ** 2 + a4 * x0 + a6 == 0
            if x0 not in roots:
                roots.append(x0)
    return sorted(roots)

def shifted(a2, a4, a6, x0):
    A = 3 * x0 + a2
    B = 3 * x0 * x0 + 2 * a2 * x0 + a4
    # check: (x + x0)^3 + a2 (x + x0)^2 + a4 (x + x0) + a6 == x^3 + A x^2 + B x identically
    for x in (1, 2, 3, -5):
        assert (x + x0) ** 3 + a2 * (x + x0) ** 2 + a4 * (x + x0) + a6 == x ** 3 + A * x * x + B * x
    return A, B

def discriminant(a2, a4, a6):
    b2, b4, b6 = 4 * a2, 2 * a4, 4 * a6
    b8 = b2 * a6 - a4 * a4  # a1 = a3 = 0: b8 = 4 a2 a6 - a4^2
    return -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6

def jinv(A, B):
    return Fraction((16 * (A * A - 3 * B)) ** 3, 16 * B * B * (A * A - 4 * B))

def other_roots(A, B):
    """Rational roots e != 0 of x^2 + A x + B (the other points of order 2)."""
    D = A * A - 4 * B
    if D < 0 or isqrt(D) ** 2 != D:
        return []
    s = isqrt(D)
    out = []
    for e2 in (-A + s, -A - s):
        assert e2 % 2 == 0
        out.append(e2 // 2)
    return sorted(set(out))

def shift0(A, B, e):
    """Move the root e of x(x^2 + A x + B) to 0: x -> x + e."""
    A2, B2 = 3 * e + A, 3 * e * e + 2 * A * e + B
    for x in (1, 2, -3):
        assert (x + e) * ((x + e) ** 2 + A * (x + e) + B) == x * (x * x + A2 * x + B2)
    return A2, B2

def isogeny_class_edges(A, B, maxcurves=16):
    """Walk the 2-isogeny class by the explicit formula E: y^2 = x(x^2 + A x + B) ->
    E/<(0,0)>: y^2 = x(x^2 - 2A x + (A^2 - 4B)) (Velu; Silverman-Tate III.4).  Every
    curve reached is 2-power isogenous to the first, so has the same rank.  Returns
    the list of models (A, B) with (0,0) the kernel point of the edge."""
    seen = {}
    todo = [(A, B)]
    models = []
    while todo and len(seen) < maxcurves:
        a, b = todo.pop()
        j = jinv(a, b)
        if j in seen:
            continue
        seen[j] = (a, b)
        pts = [(a, b)] + [shift0(a, b, e) for e in other_roots(a, b)]
        for (a1, b1) in pts:
            models.append((a1, b1))
            todo.append((-2 * a1, a1 * a1 - 4 * b1))
    return models, len(seen)

# ------------------------------------------------------------- full 2-descent (all 2-torsion rational)
def qp_class(q, p):
    """Class of the nonzero rational q in Q_p*/Q_p*^2 as a bit tuple (p = 'inf': sign)."""
    q = Fraction(q)
    assert q != 0
    if p == "inf":
        return (1 if q < 0 else 0,)
    a, b = q.numerator, q.denominator
    v = vp(abs(a), p) - vp(b, p)
    u = (abs(a) // p ** vp(abs(a), p)) * (b // p ** vp(b, p)) * (1 if a > 0 else -1)  # same class as unit part
    if p == 2:
        u8 = u % 8
        return (v % 2, ((u8 - 1) // 2) % 2, ((u8 * u8 - 1) // 8) % 2)
    return (v % 2, 0 if legendre(u, p) == 1 else 1)

def is_local_square(q, p):
    return all(x == 0 for x in qp_class(q, p))

def f2_span(vectors):
    """Row-reduced basis (as ints) of the F2-span of the given bit-ints."""
    basis = []
    for v in vectors:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
            basis.sort(reverse=True)
    return basis

def bits_to_int(bits):
    r = 0
    for b in bits:
        r = (r << 1) | b
    return r

def full_two_descent(e, rng_seed=5, max_tries=200000):
    """2-Selmer group of y^2 = (x - e1)(x - e2)(x - e3), e_i distinct integers (Silverman,
    AEC X.1.4): kappa(P) = (x - e1, x - e2) in (Q(S,2))^2.  Local images are generated
    by actual local points until they reach the known order #E(Q_v)/2E(Q_v): 4 at odd p
    (full 2-torsion), 8 at p = 2, 2 at infinity (three real roots).  Returns dim Sel^2."""
    import random
    rnd = random.Random(rng_seed)
    e1, e2, e3 = e
    S = sorted(set(factor_proved(2 * (e1 - e2) * (e1 - e3) * (e2 - e3))))
    gens = [-1] + S            # basis of Q(S,2)
    places = ["inf"] + S
    f = lambda x: (x - e1) * (x - e2) * (x - e3)
    def kappa(x):
        if x == e1:
            return (Fraction((e1 - e2) * (e1 - e3)), Fraction(e1 - e2))
        if x == e2:
            return (Fraction(e2 - e1), Fraction((e2 - e1) * (e2 - e3)))
        return (Fraction(x - e1), Fraction(x - e2))
    # the linear map G -> local classes at v, G = Q(S,2)^2 with basis (g, 1), (1, g)
    nG = 2 * len(gens)
    conditions = []   # list of (matrix rows: for each basis element of G its local image as int, W_v basis)
    for v in places:
        target = 2 if v == "inf" else (8 if v == 2 else 4)
        nbits = 1 if v == "inf" else (3 if v == 2 else 2)
        def loc(pair):
            return bits_to_int(qp_class(pair[0], v) + qp_class(pair[1], v))
        imgs = [loc(kappa(x)) for x in (e1, e2, e3)] + [0]
        # random local points
        span = f2_span(imgs)
        tries = 0
        es = sorted(e)
        while (1 << len(span)) < target:
            tries += 1
            assert tries < max_tries, ("local image not reached", v, len(span), target)
            if v == "inf":
                x = Fraction(rnd.randint(es[0], es[1])) if rnd.random() < 0.5 else Fraction(es[2] + rnd.randint(1, 10 ** 6))
            else:
                k = rnd.randint(0, 12)
                base = rnd.choice(es + [0])
                t = rnd.randint(-10 ** 6, 10 ** 6)
                x = Fraction(base) + Fraction(t) * Fraction(v) ** k if rnd.random() < 0.7 else Fraction(t, v ** (2 * rnd.randint(1, 4)))
            fx = f(x)
            if fx == 0 or x in (e1, e2):
                continue
            if v == "inf":
                if fx < 0:
                    continue
            elif not is_local_square(fx, v):
                continue
            span = f2_span(span + [loc(kappa(x))])
        assert (1 << len(span)) == target, (v, len(span), target)
        # images of the basis of G
        rows = []
        for g in gens:
            rows.append(loc((Fraction(g), Fraction(1))))
        for g in gens:
            rows.append(loc((Fraction(1), Fraction(g))))
        conditions.append((v, nbits, rows, span))
    # Sel = { c in F2^nG : sum c_i rows_i in span_v for all v }.  Reduce each row modulo
    # span_v (a row-echelon basis) and solve the linear system: reduced(sum) = 0.
    def reduce(x, span):
        for b in span:
            x = min(x, x ^ b)
        return x
    # build the matrix M (nG columns) whose kernel is Sel: stack reduced images
    cols = []
    for i in range(nG):
        col = []
        for (v, nbits, rows, span) in conditions:
            r = reduce(rows[i], span)
            col.append((r, 2 * nbits))
        cols.append(col)
    # convert each basis element's reduced image into one big int
    vecs = []
    for col in cols:
        big = 0
        for (r, w) in col:
            big = (big << w) | r
        vecs.append(big)
    # kernel dimension = nG - rank(vecs)   (reduction mod span is linear since span is row-reduced
    # and min(x, x ^ b) with sorted basis computes the canonical remainder)
    rank = len(f2_span(vecs))
    dimSel = nG - rank
    return dimSel, S

# ------------------------------------------------------------- full 2-descent, one rational 2-torsion point, L = Q(i)
def _gmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])

def _unit2_coords():
    """An explicit isomorphism (Z_2[i]^* / squares) -> F2^3, through units modulo 8 (a unit = 1 mod pi^5
    is a square; the ideal (pi^5) is {0, 4+4i} modulo 8)."""
    units = [(a, b) for a in range(8) for b in range(8) if (a + b) % 2 == 1]
    sq = set()
    for z in units:
        w = _gmul(z, z)
        w = (w[0] % 8, w[1] % 8)
        sq.add(w)
        sq.add(((w[0] + 4) % 8, (w[1] + 4) % 8))
    def coset(x):
        return frozenset(((_gmul(x, s)[0] % 8, _gmul(x, s)[1] % 8) for s in sq))
    cosets = {coset(x) for x in units}
    assert len(cosets) == 8
    # greedy basis
    basis = []
    span = {coset((1, 0)): ()}
    for g in [(0, 1), (1, 2), (3, 0), (5, 0), (1, 4), (3, 2), (1, 6), (7, 0)]:
        cg = coset(g)
        if cg in span:
            continue
        new = {}
        for c, word in span.items():
            rep = next(iter(c))
            new[coset(_gmul(rep, g))] = word + (len(basis),)
        basis.append(g)
        span.update(new)
        if len(basis) == 3:
            break
    assert len(basis) == 3 and len(span) == 8
    table = {}
    for c, word in span.items():
        bits = [0, 0, 0]
        for j in word:
            bits[j] ^= 1
        for x in c:
            table[x] = tuple(bits)
    # homomorphism check
    for x in units:
        for y in units:
            xy = _gmul(x, y)
            xy = (xy[0] % 8, xy[1] % 8)
            assert table[xy] == tuple(a ^ b for a, b in zip(table[x], table[y]))
    return table

UNIT2_BITS = _unit2_coords()

def gauss_local_bits(g, p, r=None, K=80):
    """Bits of the Gaussian integer g != 0 in (L tensor Q_p)^* / squares, L = Q(i), as a homomorphism."""
    a, b = g
    assert (a, b) != (0, 0)
    if p == 2:
        v = 0
        while (a - b) % 2 == 0:
            a, b = (a + b) // 2, (b - a) // 2
            v += 1
        return (v % 2,) + UNIT2_BITS[(a % 8, b % 8)]
    if p % 4 == 3:
        v = min(vp(a, p) if a else 10 ** 9, vp(b, p) if b else 10 ** 9)
        a1, b1 = a // p ** v, b // p ** v
        return (v % 2, 0 if legendre(a1 * a1 + b1 * b1, p) == 1 else 1)
    out = ()
    pK = p ** K
    for rr in (r, -r):
        X = (a + b * rr) % pK
        assert X != 0
        v = vp(X, p)
        assert v < K // 2
        U = (X // p ** v) % p
        out += (v % 2, 0 if legendre(U, p) == 1 else 1)
    return out

def sqrt_m1_mod(p, K):
    r = None
    for g in range(2, p):
        if legendre(g, p) == -1:
            r = pow(g, (p - 1) // 4, p)
            break
    mod = p
    for _ in range(1, K):
        mod *= p
        r = (r - (r * r + 1) * pow(2 * r, -1, mod)) % mod
    return r

def gaussian_prime_of(p):
    """a + b i with a^2 + b^2 = p (p = 1 mod 4), via a sqrt of -1 and a Gaussian gcd (Euclid)."""
    r = sqrt_m1_mod(p, 1)
    x, y = (p, 0), (r, 1)
    def gdivmod(x, y):
        n = y[0] * y[0] + y[1] * y[1]
        num = _gmul(x, (y[0], -y[1]))
        q = (round(Fraction(num[0], n)), round(Fraction(num[1], n)))
        qy = _gmul(q, y)
        return (x[0] - qy[0], x[1] - qy[1])
    while y != (0, 0):
        x, y = y, gdivmod(x, y)
    assert x[0] ** 2 + x[1] ** 2 == p
    return x

def full_two_descent_K(A, B, rng_seed=11, max_tries=400000):
    """2-Selmer group of y^2 = x (x^2 + A x + B) when A^2 - 4B = -c^2 (the quadratic factor splits
    over L = Q(i)): x^2 + A x + B = (x - rho)(x - rhobar), rho = (-A + c i)/2.  The map
    P -> (x, x - rho) in (Q x L)^* / squares has norm y^2, so it is determined by beta = x - rho in
    L(S,2) (the Q-part is N(beta)).  Basis of L(S,2): i, 1 + i, the Gaussian primes above the split
    p in S, the inert q in S (Z[i] is a PID, units +-1, +-i).  Local images generated by local points
    to the order #E(Q_p)/2E(Q_p): 4 at split p (E(Q_p)[2] has order 4), 2 at inert p, 4 at p = 2
    (E(Q_2)[2] has order 2, times 2); at infinity no condition (N(beta) > 0, E(R) connected).
    Returns dim Sel^2(E)."""
    import random
    rnd = random.Random(rng_seed)
    D1 = A * A - 4 * B
    c = isqrt(-D1)
    assert D1 < 0 and c * c == -D1
    S = sorted(set(factor_proved(2 * B * D1)))
    gens = [(0, 1), (1, 1)]
    for p in S:
        if p == 2:
            continue
        if p % 4 == 1:
            g = gaussian_prime_of(p)
            gens += [g, (g[0], -g[1])]
        else:
            gens.append((p, 0))
    f = lambda x: x * (x * x + A * x + B)
    def beta(x):   # a Gaussian integer in the class of x - rho
        q = Fraction(x)
        P, Q = q.numerator, q.denominator
        g = (2 * P + A * Q, -c * Q)      # (x - rho) * 2Q
        return (g[0] * 2 * Q, g[1] * 2 * Q)
    rows_all = []
    for p in S:
        r = sqrt_m1_mod(p, 80) if p % 4 == 1 else None
        target = 4 if (p == 2 or p % 4 == 1) else 2
        def loc(g):
            return bits_to_int(gauss_local_bits(g, p, r))
        span = []
        tries = 0
        bases = [0]
        if p % 4 == 1:
            # the p-adic roots rho, rhobar of x^2 + A x + B, as integers modulo p^80
            pK = p ** 80
            inv2 = pow(2, -1, pK)
            bases += [((-A + c * r) * inv2) % pK, ((-A - c * r) * inv2) % pK]
        while (1 << len(span)) < target:
            tries += 1
            assert tries < max_tries, ("local image not reached", p, len(span), target)
            k = rnd.randint(0, 12)
            t = rnd.randint(-10 ** 6, 10 ** 6)
            if rnd.random() < 0.7:
                x = Fraction(rnd.choice(bases)) + Fraction(t) * Fraction(p) ** k
            else:
                x = Fraction(t, p ** (2 * rnd.randint(1, 4)))
            if x == 0:
                continue
            fx = f(x)
            if fx == 0 or not is_local_square(fx, p):
                continue
            span = f2_span(span + [loc(beta(x))])
        assert (1 << len(span)) == target
        rows_all.append(([loc(g) for g in gens], span))
    def reduce(x, span):
        for b in span:
            x = min(x, x ^ b)
        return x
    vecs = []
    for i in range(len(gens)):
        big = 0
        for (rows, span) in rows_all:
            big = (big << 6) | reduce(rows[i], span)
        vecs.append(big)
    return len(gens) - len(f2_span(vecs)), S

# ------------------------------------------------------------- torsion
def count_Fp(a2, a4, a6, p):
    """#E(F_p) by a Legendre sum (p odd)."""
    s = 0
    for x in range(p):
        s += legendre(x * x * x + a2 * x * x + a4 * x + a6, p)
    return p + 1 + s

class Fp2:
    """F_{q^2} = F_q[i], q = 3 mod 4."""
    def __init__(self, q):
        self.q = q
    def mul(self, x, y):
        q = self.q
        return ((x[0] * y[0] - x[1] * y[1]) % q, (x[0] * y[1] + x[1] * y[0]) % q)
    def add(self, x, y):
        return ((x[0] + y[0]) % self.q, (x[1] + y[1]) % self.q)
    def is_square(self, x):
        # x in F_{q^2}^* is a square iff its norm is a square in F_q
        q = self.q
        nn = (x[0] * x[0] + x[1] * x[1]) % q
        return 0 if nn == 0 and x == (0, 0) else legendre(nn, q)

def count_Fq2(a2, a4, a6, q):
    F = Fp2(q)
    s = 0
    for a in range(q):
        for b in range(q):
            x = (a, b)
            x2 = F.mul(x, x)
            x3 = F.mul(x2, x)
            f = F.add(F.add(x3, F.mul((a2 % q, 0), x2)), F.add(F.mul((a4 % q, 0), x), (a6 % q, 0)))
            s += F.is_square(f)
    return q * q + 1 + s

def torsion_bounds(a2, a4, a6, Delta, nprimes=12):
    """gcd of #E(F_p) over good odd p (bound for E(Q)), and of the counts over the
    residue fields of K at good odd places (bound for E(K))."""
    gq, gk = 0, 0
    used_q, used_k = [], []
    p = 3
    while len(used_q) < nprimes or len(used_k) < nprimes:
        if all(p % r for r in range(2, isqrt(p) + 1)) and Delta % p != 0:
            Np = count_Fp(a2, a4, a6, p)
            if len(used_q) < nprimes:
                gq = gcd(gq, Np)
                used_q.append(p)
            if len(used_k) < nprimes:
                if p % 4 == 1:
                    gk = gcd(gk, Np)
                    used_k.append(p)
                elif p < 120:
                    gk = gcd(gk, count_Fq2(a2, a4, a6, p))
                    used_k.append(p)
        p += 2
    return gq, used_q, gk, used_k

def sqrt_rational(q):
    q = Fraction(q)
    if q < 0:
        return None
    a, b = q.numerator, q.denominator
    ra, rb = isqrt(a), isqrt(b)
    if ra * ra == a and rb * rb == b:
        return Fraction(ra, rb)
    return None

def is_square_in_K(q):
    """Is the rational q a square in Q(i)?  iff q or -q is a rational square."""
    return sqrt_rational(q) is not None or sqrt_rational(-q) is not None

def halvable_2torsion(roots_all, e, field):
    """Is the 2-torsion point (e, 0) divisible by 2 over the field (Q or K)?  For
    y^2 = (x - e1)(x - e2)(x - e3) with all e_i rational, (e, 0) = 2R with R rational
    over F iff e - e_j is a square in F for both other roots e_j (AEC X.1.4 / the
    2-descent map); over K, -1 is a square, so the test is up to sign."""
    test = sqrt_rational if field == "Q" else (lambda q: True if is_square_in_K(q) else None)
    others = [r for r in roots_all if r != e]
    return all(test(Fraction(e - r)) is not None for r in others)

def gauss_sqrt_frac(a, b):
    """Is a + b i (a, b Fractions) a square in Q(i)?"""
    a, b = Fraction(a), Fraction(b)
    den = a.denominator * b.denominator
    A, Bq = int(a * den * den), int(b * den * den)   # (a + bi) den^2, Gaussian integer
    nn = A * A + Bq * Bq
    s = isqrt(nn)
    if s * s != nn:
        return False
    if (s + A) % 2:
        return False
    c2, d2 = (s + A) // 2, (s - A) // 2
    c, d = isqrt(c2), isqrt(d2)
    if c * c != c2 or d * d != d2:
        return False
    return any(2 * sc * c * sd * d == Bq for sc in (1, -1) for sd in (1, -1))

def origin_halvable(A, B, field):
    """Is (0, 0) in 2E(F) on y^2 = x(x^2 + A x + B)?  The tangent at R = (x1, y1) passes
    through (0,0) iff y = l x meets the curve doubly at x1: x^2 + (A - l^2) x + B has a
    double root, i.e. B = b^2 and l^2 = A + 2b (x1 = b) or A - 2b (x1 = -b).  So:
    B a square in F with root b, and A + 2b or A - 2b a square in F (some sign of b)."""
    if field == "Q":
        b = sqrt_rational(B)
        if b is None:
            return False
        return any(sqrt_rational(A + 2 * sb) is not None for sb in (b, -b))
    # over K = Q(i): b = sqrt(B) or i*sqrt(-B)
    b = sqrt_rational(B)
    if b is not None:
        return any(gauss_sqrt_frac(A + 2 * sb, 0) for sb in (b, -b))
    c = sqrt_rational(-B)
    if c is None:
        return False
    return any(gauss_sqrt_frac(A, 2 * sc) for sc in (c, -c))

# ------------------------------------------------------------- main
def analyse(c, verbose=True):
    a1, a2, a3, a4, a6 = [int(x) for x in c["ainvs"]]
    assert a1 == 0 and a3 == 0
    Delta = discriminant(a2, a4, a6)
    assert Delta != 0
    roots = rational_roots_cubic(a2, a4, a6)
    res = {"pair": c["pair"], "curve": c["curve"], "ainvs": c["ainvs"],
           "rational_2_torsion_x": roots}
    best = None
    for x0 in roots:
        A, B = shifted(a2, a4, a6, x0)
        r, s, s1, cs, cs1, bad = isogeny_bound(A, B)
        res.setdefault("isogenies", []).append(
            {"x0": x0, "A": A, "B": B, "primes_tested": bad, "sel_phi": s, "sel_phihat": s1,
             "bound": r, "excluded_phi": {str(k): v for k, v in cs.items()},
             "excluded_phihat": {str(k): v for k, v in cs1.items()}})
        best = r if best is None else min(best, r)
    res["rank_bound_E_only"] = best
    # the whole 2-isogeny class, walked from E
    A0, B0 = shifted(a2, a4, a6, roots[0])
    models, ncurves = isogeny_class_edges(A0, B0)
    edges = []
    for (a, b) in models:
        r, s, s1, cs, cs1, bad = isogeny_bound(a, b)
        edges.append({"A": a, "B": b, "bound": r, "sel_phi": len(s), "sel_phihat": len(s1),
                      "n_excluded": len(cs) + len(cs1), "primes": bad})
        best = min(best, r)
    full = []
    for (a, b) in models:
        ro = other_roots(a, b)
        if len(ro) == 2:
            dimS, S = full_two_descent((0, ro[0], ro[1]))
            full.append({"A": a, "B": b, "roots": [0] + ro, "dim_Sel2": dimS, "bound": dimS - 2})
            best = min(best, dimS - 2)
    # E itself when its quadratic factor splits over Q(i)
    D1 = A0 * A0 - 4 * B0
    if D1 < 0 and isqrt(-D1) ** 2 == -D1:
        dimS, S = full_two_descent_K(A0, B0)
        full.append({"A": A0, "B": B0, "algebra": "Q x Q(i)", "dim_Sel2": dimS, "bound": dimS - 1})
        best = min(best, dimS - 1)
    res["class_full_2descent"] = full
    res["class_curves"] = ncurves
    res["class_edges"] = edges
    res["rank_bound"] = best
    # torsion
    gq, uq, gk, uk = torsion_bounds(a2, a4, a6, Delta)
    res["torsion_Q_divides"] = gq
    res["torsion_K_divides"] = gk
    # lower bounds: 2-torsion over Q, over K (roots of the cubic in K), 4-torsion by halving
    tQ = len(roots) + 1
    if len(roots) == 3:
        if any(halvable_2torsion(roots, e, "Q") for e in roots):
            tQ *= 2
    # over K: the quadratic cofactor's roots in K
    x0 = roots[0]
    A, B = shifted(a2, a4, a6, x0)
    disc = A * A - 4 * B
    if len(roots) == 3:
        rootsK = roots
    elif sqrt_rational(-disc) is not None:   # roots (-A +- i sqrt(-disc))/2 in K
        rootsK = None   # non-rational roots in K: not needed for the curves here
    else:
        rootsK = roots
    if rootsK is None:
        tK = 4      # the quadratic factor splits over K: full 2-torsion over K (halving over K not tested)
    else:
        tK = (4 if len(rootsK) == 3 else 2)
    if rootsK is not None and len(rootsK) == 3:
        n_half = sum(1 for e in rootsK if halvable_2torsion(rootsK, e, "K"))
        # each halvable 2-torsion point gives 4 points of order 4 over K
        if n_half:
            tK = 4 * (1 + n_half) if n_half < 3 else 16
    if len(roots) == 1 and origin_halvable(A, B, "Q"):
        tQ *= 2
    if len(roots) == 1 and rootsK is not None and len(rootsK) == 1 and origin_halvable(A, B, "K"):
        tK *= 2
    res["torsion_Q_lower"] = tQ
    res["torsion_K_lower"] = tK
    res["torsion_Q_exact"] = tQ if tQ == gq else None
    if len(roots) == 1 and not origin_halvable(A, B, "Q"):
        # E(Q)[2] = Z/2 and (0,0) is not in 2E(Q): the 2-primary part is Z/2;
        # the odd part divides the odd part of the gcd of the point counts
        g_odd = gq
        while g_odd % 2 == 0:
            g_odd //= 2
        if g_odd == 1:
            res["torsion_Q_exact"] = 2
    res["torsion_K_exact"] = tK if tK == gk else None
    # PARI's own answers, as a comparison only
    E = pari.ellinit([a1, a2, a3, a4, a6])
    res["pari_elltors"] = int(pari.elltors(E)[0])
    res["pari_conductor_digits"] = len(str(int(pari.ellglobalred(E)[0])))
    try:
        rk = pari.ellrank(E)
        res["pari_ellrank"] = [int(rk[0]), int(rk[1])]
    except Exception as ex:
        res["pari_ellrank"] = str(ex)[:80]
    return res

def selftest():
    print("self-test")
    # quartic local solvability against brute force over residue classes
    import random
    rnd = random.Random(3)
    # 1. curves with known rank (congruent number curves y^2 = x^3 - N^2 x: N = 5, 6, 7 rank 1; N = 1, 2, 3 rank 0)
    for N in [1, 2, 3, 5, 6, 7, 34, 41]:
        r, s, s1, cs, cs1, bad = isogeny_bound(0, -N * N)
        lo, hi = [int(x) for x in pari.ellrank(pari.ellinit([0, 0, 0, -N * N, 0]))[:2]]
        print("  y^2 = x^3 - %d^2 x: bound %d (PARI ellrank: %d..%d)" % (N, r, lo, hi))
        assert r >= lo
    # 2. random small curves with a 2-torsion point: bound >= PARI's lower bound, and of the same parity as Sel
    nt = 0
    nimg = 0
    while nt < 40:
        A = rnd.randint(-60, 60); B = rnd.randint(-60, 60)
        if B == 0 or A * A - 4 * B == 0:
            continue
        r, s, *_ = isogeny_bound(A, B)
        E = pari.ellinit([0, A, 0, B, 0])
        rk = pari.ellrank(E)
        lo, hi = int(rk[0]), int(rk[1])
        assert r >= lo, (A, B, r, lo)
        # the image alpha(P) of every point PARI found must be in the Selmer set
        for P in rk[3]:
            x = Fraction(int(pari.numerator(P[0])), int(pari.denominator(P[0])))
            if x == 0:
                continue
            num = x.numerator * x.denominator
            sf = 1 if num > 0 else -1
            for q, e in factor_proved(abs(num)).items():
                if e % 2:
                    sf *= q
            assert sf in s, (A, B, P, sf, s)
            nimg += 1
        nt += 1
    print("  40 random curves y^2 = x(x^2 + A x + B): bound >= PARI's lower bound for the rank;"
          " the images of %d points found by PARI are all in the Selmer sets" % nimg)
    # 3. negative test: a forged exclusion certificate must be refused
    A, B = 0, -36   # rank 1; d = 1 always has the point (1 : 0 : 1)
    for p in [2, 3, "inf"]:
        assert not verify_cert(A, B, 1, p), "forged certificate accepted"
    # the point (-3, 9) of y^2 = x^3 - 36 x... y^2 = 27*... use (x, y) = (-2, 8): 64 = -8 + 72
    assert 8 ** 2 == (-2) ** 3 - 36 * (-2)
    d = -2  # alpha((-2, 8)) = -2 mod squares, so d = -2 is in the image: it cannot be excluded anywhere
    for p in [2, 3, "inf"]:
        assert not verify_cert(A, B, d, p), "forged certificate for an element of the image accepted"
    print("  negative test: forged exclusions (d = 1, and d = -2 = alpha of the point (-2, 8)) are refused")
    # 4. full 2-descent against PARI on random curves with all 2-torsion rational
    nt = 0
    while nt < 25:
        e2 = rnd.randint(-80, 80); e3 = rnd.randint(-80, 80)
        if len({0, e2, e3}) < 3:
            continue
        dimS, S = full_two_descent((0, e2, e3))
        E = pari.ellinit([0, -(e2 + e3), 0, e2 * e3, 0])
        rk = pari.ellrank(E)
        lo, hi = int(rk[0]), int(rk[1])
        assert dimS - 2 >= lo, (e2, e3, dimS, lo)
        # PARI's ellrank returns [r, R, s, L] with R + s the bound from the 2-Selmer group alone
        assert dimS - 2 == hi + int(rk[2]), (e2, e3, dimS, hi, int(rk[2]))
        nt += 1
    print("  25 random curves y^2 = x(x - e2)(x - e3): dim Sel^2 - 2 equals PARI's R + s (its 2-Selmer bound) every time")
    for N, lo in [(5, 1), (34, 2), (41, 2)]:
        dimS, S = full_two_descent((0, N, -N))
        assert dimS - 2 >= lo
    print("  congruent number curves N = 5, 34, 41 (ranks 1, 2, 2): full 2-descent bounds 1, 2, 2 not below the rank")
    # 5. full 2-descent over Q x Q(i) against PARI on curves y^2 = x((x + a)^2 + b^2)
    nt = 0
    while nt < 25:
        a = rnd.randint(-40, 40); b = rnd.randint(1, 40)
        A, B = 2 * a, a * a + b * b
        dimS, S = full_two_descent_K(A, B)
        E = pari.ellinit([0, A, 0, B, 0])
        rk = pari.ellrank(E)
        lo, hi = int(rk[0]), int(rk[1])
        assert dimS - 1 >= lo, (a, b, dimS, lo)
        assert dimS - 1 == hi + int(rk[2]), (a, b, dimS, hi, int(rk[2]))
        nt += 1
    print("  25 random curves y^2 = x((x + a)^2 + b^2): dim Sel^2 - 1 over Q x Q(i) equals PARI's R + s every time")
    # 6. negative test: a local image that cannot be completed (too few local points) is refused,
    #    never used: an incomplete image would make the Selmer group too small
    for fn, arg in [(full_two_descent_K, (2, 5)), (full_two_descent, ((0, 5, -5),))]:
        try:
            fn(*arg, max_tries=1)
            raise SystemExit("negative test failed: an incomplete local image was accepted")
        except AssertionError:
            pass
    print("  negative test: with too few local points the full descents refuse to return a bound")
    print("self-test passed")

def main(argv):
    if argv and argv[0] == "--selftest":
        selftest()
        return
    if argv and argv[0] == "--ainvs":
        # one curve given by its a-invariants, e.g. the L-value curve of (440,279):
        #   python3 indep_rank.py --ainvs 0,-491040,0,-13400146081,0 440,279   (y^2 = x(x - 719^2)(x + 161^2))
        ainvs = [int(x) for x in argv[1].split(",")]
        pair = [int(x) for x in argv[2].split(",")] if len(argv) > 2 else []
        res = analyse({"pair": pair, "curve": "given", "ainvs": ainvs})
        print(json.dumps(res))
        E = pari.ellinit(ainvs)
        rk = pari.ellrank(E)
        print("SUMMARY pair %s: rank bound %s from every 2-descent here (2-isogeny class of %d curves); torsion over Q %s,"
              " over K %s; PARI ellrank [r, R, s] = %s" % (pair, res["rank_bound"], res["class_curves"],
              res["torsion_Q_exact"], res["torsion_K_exact"], [int(rk[0]), int(rk[1]), int(rk[2])]))
        return
    curves = json.load(open(argv[0]))
    seen = {}
    allzero = True
    for c in curves:
        key = tuple(int(x) for x in c["ainvs"])
        if key in seen:
            print("pair %s: same curve as pair %s" % (c["pair"], seen[key]["pair"]))
            continue
        res = analyse(c)
        seen[key] = res
        print(json.dumps(res), flush=True)
        print("SUMMARY pair %s curve %s: rank bound %s (least of: 2-isogeny descents on every edge of the 2-isogeny class of %d curves, full 2-descents on its members with full rational 2-torsion, and on E over Q x Q(i) when possible; 2-isogeny descent on E alone: %s); torsion over Q: lower %d, divides %d -> %s;"
              " over K: lower %d, divides %d -> %s; PARI elltors %d, ellrank %s, conductor %d digits"
              % (res["pair"], res["curve"], res["rank_bound"], res["class_curves"], res["rank_bound_E_only"], res["torsion_Q_lower"], res["torsion_Q_divides"],
                 res["torsion_Q_exact"], res["torsion_K_lower"], res["torsion_K_divides"], res["torsion_K_exact"],
                 res["pari_elltors"], res["pari_ellrank"], res["pari_conductor_digits"]), flush=True)
        allzero &= res["rank_bound"] == 0
    print("ALL rank bounds 0: %s" % allzero)

if __name__ == "__main__":
    main(sys.argv[1:])
