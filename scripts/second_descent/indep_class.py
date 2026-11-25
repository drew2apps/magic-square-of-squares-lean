#!/usr/bin/env python3
"""
indep_class.py -- an independent re-implementation of the class step of the
first descent over Q(i) (GAUSS_DESCENT.md sections 1-2, SECOND_DESCENT.md
section 1), written from the prose only, before reading gauss_descent.py or
degenerate_tuples.py.  Exact integer arithmetic throughout (no floating point).

For a ratio m:n (coprime, m > n > 0) a solution is a coprime pair (u, v) at
which the five conditions hold:

    N(alpha0) square, N(alpha8) square, 2S*N(gamma4) square,
    2*N(gammaL) square, 2*N(gammaU) square,

    alpha0 = m v + i n u,  alpha8 = n v + i m u,  gamma4 = u + i v,
    gammaL = (m+in)u - i(m-in)v = (mu - nv) + i(nu - mv),
    gammaU = (m+in)u + i(m-in)v = (mu + nv) + i(nu + mv).

Global step (the Lemma of GAUSS_DESCENT.md section 1, re-derived here in the
function `support_set`): the class delta_j of l_j(u,v) in K*/K*^2 lies in a
finite set D_j.  Local step: for every prime p, some primitive (u,v) in Z_p^2
must realise the five classes at every place above p.  For each p we build a
set W_p of local class patterns that contains the pattern of every primitive
pair (u,v) of Z_p^2, with a wildcard where the working precision does not
decide a class; a tuple survives p if it matches some element of W_p.

Differences from the description of the original program, on purpose:
  * primitive pairs are covered as lambda*(x,1), x in Z_p, and lambda*(1,p*y),
    y in Z_p, with lambda over representatives of Z_p^*/Z_p^*2 (the classes
    of l_j are not invariant under scaling (u,v) by a p-adic unit, so the
    scaling has to be enumerated or the set W_p can be too small);
  * a pattern is matched against the tuples as a product set (the compatible
    elements of each D_j), with numpy boolean arrays;
  * the labels at 1+i are cosets of the unit squares modulo pi^5, computed
    from scratch (and the script asserts that there are 8 of them).

Usage:
    python3 indep_class.py [--pmax=200] [--depth=D] m,n ...
    python3 indep_class.py --selftest       # controls and negative tests
Output: one block per ratio: tuples left after each prime, the tuples left,
the degenerate tuples (computed by an exact square test in Z[i]), and the
classes of the three quartics of the second descent at the tuples left.
"""
import sys
from math import gcd, isqrt
import numpy as np

# ---------------------------------------------------------------- Z[i]
def gmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])

def gconj(x):
    return (x[0], -x[1])

def gnorm(x):
    return x[0] * x[0] + x[1] * x[1]

def gpow(x, e):
    r = (1, 0)
    for _ in range(e):
        r = gmul(r, x)
    return r

I_UNIT = (0, 1)
PI2 = (1, 1)  # 1 + i

def is_square_zi(x):
    """Exact test: is x (nonzero, in Z[i]) a square in Z[i] (equivalently in K)."""
    a, b = x
    nn = a * a + b * b
    s = isqrt(nn)
    if s * s != nn:
        return False
    # x = (c + d i)^2: c^2 - d^2 = a, 2cd = b, c^2 + d^2 = s
    if (s + a) % 2:
        return False
    c2, d2 = (s + a) // 2, (s - a) // 2
    c, d = isqrt(c2), isqrt(d2)
    if c * c != c2 or d * d != d2:
        return False
    for sc in (1, -1):
        for sd in (1, -1):
            if 2 * (sc * c) * (sd * d) == b:
                return True
    return False

def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]

def factor(n):
    n = abs(n)
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f

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

def sqrt_minus1(p, K):
    """r with r^2 = -1 mod p^K (p = 1 mod 4), by search and Hensel lifting."""
    assert p % 4 == 1
    r = None
    for g in range(2, p):
        if legendre(g, p) == -1:
            r = pow(g, (p - 1) // 4, p)
            break
    assert (r * r + 1) % p == 0
    mod = p
    for _ in range(1, K):
        mod2 = mod * p
        # r_new = r - (r^2+1)/(2r) mod mod2
        r = (r - (r * r + 1) * pow(2 * r, -1, mod2)) % mod2
        mod = mod2
    assert (r * r + 1) % (p ** K) == 0
    return r

def gaussian_prime(p):
    """a + b i with a^2 + b^2 = p, a > 0 odd, b > 0 even (p = 1 mod 4)."""
    for a in range(1, isqrt(p) + 1):
        b2 = p - a * a
        b = isqrt(b2)
        if b * b == b2 and a % 2 == 1:
            return (a, b)
    raise ValueError(p)

# ---------------------------------------------------------------- the forms
FORMS = ["alpha0", "alpha8", "gamma4", "gammaL", "gammaU"]

def form_coeffs(m, n):
    """Each form as (re_u, re_v, im_u, im_v): value = (re_u u + re_v v) + i (im_u u + im_v v)."""
    return [
        (0, m, n, 0),        # alpha0 = m v + i n u
        (n * 0, n, m, 0),    # alpha8 = n v + i m u
        (1, 0, 0, 1),        # gamma4 = u + i v
        (m, -n, n, -m),      # gammaL = (mu - nv) + i(nu - mv)
        (m, n, n, m),        # gammaU = (mu + nv) + i(nu + mv)
    ]

def form_value(c, u, v):
    return (c[0] * u + c[1] * v, c[2] * u + c[3] * v)

def condition_constants(m, n):
    S = m * m + n * n
    return [1, 1, 2 * S, 2, 2]

def resultants(m, n):
    """Determinant of the coefficients of (u, v) in l and conj(l), up to a unit,
    computed from the coefficients (and compared with the values of the prose)."""
    out = []
    for c in form_coeffs(m, n):
        A = (c[0], c[2])  # coefficient of u
        B = (c[1], c[3])  # coefficient of v
        R = (gmul(A, gconj(B))[0] - gmul(gconj(A), B)[0],
             gmul(A, gconj(B))[1] - gmul(gconj(A), B)[1])
        out.append(R)
    # R is purely imaginary: R = A Bbar - conj(A Bbar)
    vals = []
    for R in out:
        assert R[0] == 0
        vals.append(abs(R[1]))
    S, D = m * m + n * n, m * m - n * n
    assert vals == [2 * m * n, 2 * m * n, 2, 2 * abs(D), 2 * abs(D)], vals
    return vals

# ---------------------------------------------------------------- global sets
def support_set(c, R):
    """The finite set of squarefree delta (up to squares of units) allowed by the
    Lemma for a linear form l with N(l) in c*Q^2 and resultant R (a rational integer):
      inert q:  q | delta only if q | R;
      split p:  v_p(c) odd  -> exactly one of pi, pibar;
                v_p(c) even -> none or both (= p), both only if p | R;
      1 + i:    divides delta iff v_2(c) odd;
      units:    i^a, a in {0, 1}.
    """
    fc = factor(c)
    fR = factor(R)
    choices = []  # list of lists of Gaussian integers, one factor chosen from each
    # 1 + i
    choices.append([PI2] if fc.get(2, 0) % 2 == 1 else [(1, 0)])
    primes = sorted(set(fc) | set(fR))
    for p in primes:
        if p == 2:
            continue
        if p % 4 == 3:
            assert fc.get(p, 0) % 2 == 0  # c is a norm up to squares
            choices.append([(1, 0), (p, 0)] if p in fR else [(1, 0)])
        else:
            if fc.get(p, 0) % 2 == 1:
                g = gaussian_prime(p)
                choices.append([g, gconj(g)])
            else:
                choices.append([(1, 0), (p, 0)] if p in fR else [(1, 0)])
    choices.append([(1, 0), I_UNIT])
    out = [(1, 0)]
    for ch in choices:
        out = [gmul(x, y) for x in out for y in ch]
    return out

def global_sets(m, n, consts=None):
    if consts is None:
        consts = condition_constants(m, n)
    R = resultants(m, n)
    return [support_set(consts[j], R[j]) for j in range(5)]

# ---------------------------------------------------------------- local classes
def unit_classes_mod_pi5():
    """Labels of the unit classes of Z_2[i] modulo squares.  A unit = 1 mod pi^5 is
    a square (Hensel: v(y^2 - x) > 2 v(2y) = 4 at y = 1), and conversely the
    square class of a unit depends only on it modulo pi^5.  Work modulo 8 = pi^6
    (unit times), where the ideal (pi^5) is {0, 4+4i}."""
    units = [(a, b) for a in range(8) for b in range(8) if (a + b) % 2 == 1]
    sq = set()
    for z in units:
        w = gmul(z, z)
        w = (w[0] % 8, w[1] % 8)
        sq.add(w)
        sq.add(((w[0] + 4) % 8, (w[1] + 4) % 8))
    label = {}
    for x in units:
        coset = sorted(((gmul(x, s)[0] % 8, gmul(x, s)[1] % 8) for s in sq))
        label[x] = tuple(coset)
    names = sorted(set(label.values()))
    assert len(units) == 32 and len(sq) == 4 and len(names) == 8, (len(sq), len(names))
    idx = {nm: k for k, nm in enumerate(names)}
    return {x: idx[label[x]] for x in units}

UNIT2 = unit_classes_mod_pi5()

class Places:
    """The places of K above a rational prime p, with a label function that is
    given a value known modulo p^k (as an integer pair) and returns, per place,
    (valuation mod 2, unit class) or None if the precision does not decide it."""

    def __init__(self, p, K=60):
        self.p = p
        if p == 2:
            self.kind = "ramified"
            self.n = 1
        elif p % 4 == 3:
            self.kind = "inert"
            self.n = 1
        else:
            self.kind = "split"
            self.n = 2
            self.K = K
            self.r = sqrt_minus1(p, K)
            self.pK = p ** K

    def labels(self, x, k):
        """x = (a, b) known modulo p^k (k = None: exact)."""
        p = self.p
        a, b = x
        if self.kind == "split":
            out = []
            for r in (self.r, -self.r):
                if k is None:
                    X = a + b * r
                    if X == 0:
                        out.append(None)
                        continue
                    # exact value: valuation in Q_p of a + b r, r a p-adic sqrt(-1) mod p^K
                    X %= self.pK
                    v = vp(X, p)
                    assert v < self.K // 2, "precision"
                else:
                    pk = p ** k
                    X = (a + b * r) % pk
                    if X == 0:
                        out.append(None)
                        continue
                    v = vp(X, p)
                U = (X // p ** v) % p
                out.append((v % 2, legendre(U, p)))
            return tuple(out)
        if self.kind == "inert":
            if k is not None:
                pk = p ** k
                a, b = a % pk, b % pk
                if a == 0 and b == 0:
                    return (None,)
            elif a == 0 and b == 0:
                return (None,)
            v = min(vp(a, p), vp(b, p))
            a1, b1 = a // p ** v, b // p ** v
            return ((v % 2, legendre(a1 * a1 + b1 * b1, p)),)
        # p = 2, place pi = 1 + i; precision in pi-units is 2k
        prec = None if k is None else 2 * k
        if k is not None:
            pk = 2 ** k
            a, b = a % pk, b % pk
            if a == 0 and b == 0:
                return (None,)
        elif a == 0 and b == 0:
            return (None,)
        v = 0
        while (a - b) % 2 == 0:
            # divide by 1 + i: (a + bi)(1 - i)/2
            a, b = (a + b) // 2, (b - a) // 2
            v += 1
            if prec is not None and v >= prec:
                return (None,)
            if a == 0 and b == 0:
                return (None,)
        if prec is not None and prec - v < 5:
            return (None,)
        return ((v % 2, UNIT2[(a % 8, b % 8)]),)

def scaling_reps(p):
    """Representatives of Z_p^* modulo squares (as integers)."""
    if p == 2:
        return [1, 3, 5, 7]
    g = next(g for g in range(2, p) if legendre(g, p) == -1)
    return [1, g]

def local_patterns(m, n, p, depth):
    """A set of patterns (one entry per form, per place: label or None) that
    contains the pattern of every primitive pair of Z_p^2."""
    pl = Places(p)
    co = form_coeffs(m, n)
    pats = set()
    for lam in scaling_reps(p):
        # chart A: (u, v) = lam*(x, 1), x mod p^k;  chart B: (u, v) = lam*(1, p y), y mod p^k
        for chart in ("A", "B"):
            stack = [(x0, 1) for x0 in range(p)]
            while stack:
                x0, k = stack.pop()
                if chart == "A":
                    u, v, prec = lam * x0, lam, k
                else:
                    u, v, prec = lam, lam * p * x0, k + 1
                pat = tuple(pl.labels(form_value(c, u, v), prec) for c in co)
                undecided = any(l is None for lj in pat for l in lj)
                if undecided and k < depth:
                    pk = p ** k
                    for d in range(p):
                        stack.append((x0 + d * pk, k + 1))
                else:
                    pats.add(pat)
    return pl, pats

# ---------------------------------------------------------------- the descent
def bad_primes(m, n):
    S, D = m * m + n * n, m * m - n * n
    return sorted(set(factor(2 * m * n * S * D)))

def descent(m, n, pmax=200, depth=None, consts=None, sets=None, verbose=False, plant=None):
    Ds = sets if sets is not None else global_sets(m, n, consts)
    sizes = [len(D) for D in Ds]
    grids = np.meshgrid(*[np.arange(s) for s in sizes], indexing="ij")
    idx = [g.ravel() for g in grids]
    alive = np.ones(idx[0].shape[0], dtype=bool)
    plist = sorted(set(primes_upto(pmax)) | set(bad_primes(m, n)))
    history = [("start", int(alive.sum()))]
    for p in plist:
        if depth is not None:
            dep = depth
        elif p == 2:
            dep = 12
        elif p < 30:
            dep = 6
        elif p < 200:
            dep = 4
        else:
            dep = 3
        pl, pats = local_patterns(m, n, p, dep)
        labs = [[pl.labels(x, None) for x in D] for D in Ds]
        ok = np.zeros_like(alive)
        cache = {}
        for pat in pats:
            masks = []
            for j in range(5):
                key = (j, pat[j])
                if key not in cache:
                    cache[key] = np.array(
                        [all(pw is None or pw == lw for pw, lw in zip(pat[j], labs[j][e]))
                         for e in range(sizes[j])], dtype=bool)
                masks.append(cache[key])
            sel = masks[0][idx[0]]
            for j in range(1, 5):
                sel &= masks[j][idx[j]]
            ok |= sel
        alive &= ok
        history.append((p, int(alive.sum())))
        if plant is not None:
            assert alive[plant], ("planted tuple removed at p =", p)
        if verbose:
            print("   p=%d patterns=%d left=%d" % (p, len(pats), int(alive.sum())), flush=True)
        if not alive.any():
            break
    left = [tuple(Ds[j][idx[j][t]] for j in range(5)) for t in np.nonzero(alive)[0]]
    return Ds, history, left, idx

# ---------------------------------------------------------------- degenerate tuples
def class_in(x, D):
    """The unique element d of D with x*d a square in K (x nonzero); None if none."""
    hits = [d for d in D if is_square_zi(gmul(x, d))]
    assert len(hits) <= 1, (x, hits)
    return hits[0] if hits else None

def degenerate_points(m, n):
    S = m * m + n * n
    k = isqrt(S)
    if k * k == S:
        return [(1, 1), (-1, 1)]
    k = isqrt(S // 2)
    if S % 2 == 0 and 2 * k * k == S:
        return [(0, 1), (1, 0)]
    return []

def tuple_at(m, n, u, v, Ds):
    co = form_coeffs(m, n)
    return tuple(class_in(form_value(co[j], u, v), Ds[j]) for j in range(5))

# ---------------------------------------------------------------- quartic classes
def rational_class(x):
    """If the class of x in K*/K*^2 is that of a rational number, return the
    squarefree integer d (sign and 2 included; d and -d are the same class) with
    the smallest |d| and d > 0; else None."""
    primes = sorted(set(factor(gnorm(x))) | {2})
    best = None
    for mask in range(1 << len(primes)):
        d = 1
        for t, q in enumerate(primes):
            if mask >> t & 1:
                d *= q
        if is_square_zi(gmul(x, (d, 0))):
            if best is None or d < best:
                best = d
    return best

def quartic_classes(t):
    d0, d8, d4, dL, dU = t
    Qb = gmul(gmul(d0, gconj(d8)), gmul(d4, gconj(d4)))          # L1 L3 L5 L6
    Qp14 = gmul(gmul(dL, gconj(dU)), gmul(d0, d8))                # q+ L1 L4 (-1 is a square)
    Qp23 = gmul(gmul(dL, gconj(dU)), gmul(gconj(d0), gconj(d8)))  # q+ L2 L3
    return {"L1L3L5L6": rational_class(Qb), "q+L1L4": rational_class(Qp14),
            "q+L2L3": rational_class(Qp23)}

def fmt(x):
    a, b = x
    if b == 0:
        return "%d" % a
    if a == 0:
        return "%di" % b
    return "%d%+di" % (a, b)

def run(m, n, pmax, depth, verbose=False):
    Ds, hist, left, _ = descent(m, n, pmax, depth, verbose=verbose)
    degs = [tuple_at(m, n, u, v, Ds) for (u, v) in degenerate_points(m, n)]
    return Ds, hist, left, degs

def report(m, n, pmax, depth):
    Ds, hist, left, degs = run(m, n, pmax, depth)
    print("ratio %d,%d  S=%d  |D_j|=%s" % (m, n, m * m + n * n, [len(D) for D in Ds]))
    print("  left after: " + " ".join("%s:%d" % h for h in hist))
    print("  tuples left: %d" % len(left))
    for t in left:
        print("    (" + ", ".join(fmt(x) for x in t) + ")  quartic classes " + str(quartic_classes(t)))
    print("  degenerate points %s, tuples:" % degenerate_points(m, n))
    for t in degs:
        assert None not in t, "a degenerate value is not in its D_j"
        print("    (" + ", ".join(fmt(x) for x in t) + ")")
    same = sorted(left) == sorted(degs) and len(set(degs)) == len(degs)
    print("  left == the degenerate tuples (and these are distinct): %s" % same)
    return same, left, degs

# ---------------------------------------------------------------- self-test
def selftest():
    print("self-test")
    # 1. square test
    for x in [(3, 4), (-7, 24), (0, 2), (5, 12)]:
        assert is_square_zi(gmul(x, x)) and not is_square_zi(gmul((0, 1), gmul(x, x)))
    assert is_square_zi((0, 2)) and is_square_zi((-1, 0)) and not is_square_zi((0, 1)) and not is_square_zi((2, 0))
    # 2. labels: x and x*z^2 have the same labels at every place; x and x*nonsquare differ somewhere
    import random
    rnd = random.Random(1)
    for p in [2, 3, 5, 7, 13]:
        pl = Places(p)
        for _ in range(300):
            x = (rnd.randint(-999, 999), rnd.randint(-999, 999))
            z = (rnd.randint(-99, 99), rnd.randint(-99, 99))
            if x == (0, 0) or z == (0, 0):
                continue
            assert pl.labels(x, None) == pl.labels(gmul(x, gmul(z, z)), None)
            # a known value modulo p^k that decides gives the exact label
            for k in (3, 6):
                lk = pl.labels(x, k)
                le = pl.labels(x, None)
                for a, b in zip(lk, le):
                    assert a is None or a == b
    # 3. controls (GAUSS_DESCENT.md section 3 numbers): no trivial point -> 0 tuples; trivial point -> 2
    for (m, n, expect) in [(4, 3, 2), (15, 8, 2), (820, 799, 0), (889, 320, 0)]:
        Ds, hist, left, degs = run(m, n, 40, None)
        print("  control %d,%d: left %d (expected %d)  history %s" % (m, n, len(left), expect, hist))
        assert len(left) == expect
        if expect:
            assert sorted(left) == sorted(degs)
    # 4. planted systems: random m, n, (u0, v0); constants = squarefree part of the norms
    #    at (u0, v0), so (u0, v0) is a solution; its tuple must survive every prime.
    rnd = random.Random(7)
    planted = 0
    while planted < 12:
        m = rnd.randint(2, 40); n = rnd.randint(1, m - 1)
        u0 = rnd.randint(-30, 30); v0 = rnd.randint(1, 30)
        if gcd(m, n) != 1 or gcd(u0, v0) != 1:
            continue
        co = form_coeffs(m, n)
        vals = [form_value(c, u0, v0) for c in co]
        if any(gnorm(x) == 0 for x in vals):
            continue
        consts = []
        for x in vals:
            N = gnorm(x)
            f = factor(N)
            c = 1
            for q, e in f.items():
                if e % 2:
                    c *= q
            consts.append(c)
        Ds = global_sets(m, n, consts)
        t = tuple_at(m, n, u0, v0, Ds)
        assert None not in t, ("planted value outside its support set", m, n, u0, v0)
        pos = 0
        sizes = [len(D) for D in Ds]
        for j in range(5):
            pos = pos * sizes[j] + Ds[j].index(t[j])
        Ds2, hist, left, _ = descent(m, n, 40, None, sets=Ds, plant=pos)
        assert t in left
        planted += 1
    print("  planted systems: 12 of 12 tuples of a real solution survive every prime up to 40 and the bad primes")
    # 5. negative tests: mutated tuples must be removed (each is a member of the D_j, so the descent really runs)
    Dx = global_sets(4, 3)
    sizes = [len(D) for D in Dx]
    Ds, hist, left, degs = run(4, 3, 40, None)
    exercised = 0
    for j, mult in [(0, (0, 1)), (1, (3, 0)), (3, (0, 1)), (4, (0, 1))]:
        t0 = list(degs[0])
        # the class of delta_j * mult, as the member of D_j that represents it
        mutated = class_in(gmul(t0[j], mult), Dx[j])
        assert mutated is not None and mutated != t0[j], ("mutation not in D_j or not a change", j, mult)
        t0[j] = mutated
        assert tuple(t0) not in degs
        pos = 0
        for i in range(5):
            pos = pos * sizes[i] + Dx[i].index(t0[i])
        removed = False
        try:
            descent(4, 3, 40, None, plant=pos)
        except AssertionError as e:
            # only the descent's own "planted tuple removed" assertion counts as a removal
            removed = bool(e.args) and isinstance(e.args[0], tuple) and e.args[0][:1] == ("planted tuple removed at p =",)
            if not removed: raise
        assert removed, ("negative test failed: the mutated tuple survived", j, mult)
        exercised += 1
    assert exercised == 4, "negative test not exercised"
    # a forged claim (a mutated tuple listed as left) is refused by the comparison used in report()
    forged = sorted(degs[:1] + [tuple(t0)])
    assert sorted(left) != forged
    print("  negative test: 4 mutated degenerate tuples of (4,3) (delta_0 * i, delta_8 * 3, delta_L * i, delta_U * i,"
          " each a member of D_j) are each removed by the descent; a forged list of tuples is refused")
    print("self-test passed")

def main(argv):
    pmax, depth = 200, None
    pairs = []
    for a in argv:
        if a == "--selftest":
            selftest()
            return
        if a.startswith("--pmax="):
            pmax = int(a.split("=")[1])
        elif a.startswith("--depth="):
            depth = int(a.split("=")[1])
        else:
            m, n = map(int, a.split(","))
            pairs.append((m, n))
    allok = True
    for (m, n) in pairs:
        ok, _, _ = report(m, n, pmax, depth)
        allok &= ok
    print("ALL: left == degenerate tuples for every ratio: %s" % allok)
    assert allok

if __name__ == "__main__":
    main(sys.argv[1:])
