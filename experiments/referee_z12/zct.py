"""z12's own evaluation of the conic formula on y^2 = x(x - A^2)(x - C^2), a = (1, dl, dl), with a free assignment of the
components of b = (b_r, b_m, b_*) to the forms l_r = u_m - u_*, l_m = u - C, l_* = u - A (perm), random local points,
all 8 sign choices, a constancy check at every place, and the real place included.  Hilbert symbols by hand."""
import random, cypari2
from fractions import Fraction as Fr
from sympy import factorint
pari = cypari2.Pari()
def vp(n, p):
    v = 0
    while n % p == 0: n //= p; v += 1
    return v
def cls(x, p):
    """(v mod 2, unit residue data) of a non-zero rational in Q_p*/Q_p*^2: odd p -> legendre bit; p = 2 -> unit mod 8"""
    x = Fr(x); n, d = x.numerator, x.denominator
    v = vp(abs(n), p) - vp(d, p); u = (n // p ** vp(abs(n), p)) * (d // p ** vp(d, p))
    return v, (u % 8 if p == 2 else u % p)
def issq(x, p):
    if x == 0: return False
    v, u = cls(x, p)
    if v % 2: return False
    return u == 1 if p == 2 else pow(u, (p - 1) // 2, p) == 1
def hil_vu(a, b, p):
    (al, u), (be, w) = a, b
    al %= 2; be %= 2
    if p == 2:
        eps = lambda t: ((t - 1) // 2) % 2; om = lambda t: ((t * t - 1) // 8) % 2
        e = eps(u) * eps(w) + al * om(w) + be * om(u)
        return -1 if e % 2 else 1
    leg = lambda t: 1 if pow(t, (p - 1) // 2, p) == 1 else -1
    s = (-1) ** (al * be * ((p - 1) // 2)) * (leg(u) ** be) * (leg(w) ** al)
    return s
def padic_cls(f, p):
    """class data of a PARI p-adic f (non-zero, enough precision)"""
    v = int(pari.valuation(f, p)); u = f / pari(p) ** v
    m = 8 if p == 2 else p
    r = int(pari.lift(pari.Mod(pari.truncate(u), m))) if True else None
    return v, r, int(pari.padicprec(f, p)) - v
def local_points(A, C, dl, p, rng, want=3):
    out = []
    bases = [0, A, -A, C, -C]
    for it in range(60000):
        if len(out) >= want: break
        base = rng.choice(bases); k = rng.randint(-8, 10)
        w = rng.randrange(1, p ** 4) * rng.choice((1, -1))
        u = Fr(base) + Fr(w) * Fr(p) ** k
        if u == 0 and base == 0: continue
        if issq(Fr(u * u - A * A) / dl, p) and issq(Fr(u * u - C * C) / dl, p) and u not in out: out.append(u)
    return out
def local_forms(A, C, dl, p, rng, N=90):
    """class data of (l_r, l_m, l_*) at several local points and all sign choices"""
    pts = local_points(A, C, dl, p, rng)
    res = []
    for u in pts:
        qm = Fr(u * u - A * A) / dl; qs = Fr(u * u - C * C) / dl
        um = pari('sqrt((%d)/(%d)+O(%d^%d))' % (qm.numerator, qm.denominator, p, N))
        us = pari('sqrt((%d)/(%d)+O(%d^%d))' % (qs.numerator, qs.denominator, p, N))
        for s1 in (1, -1):
            for s2 in (1, -1):
                for s3 in (1, -1):
                    U = s1 * u
                    if U == A or U == C: continue
                    lr = s2 * um - s3 * us
                    if lr == 0: continue
                    d = padic_cls(lr, p)
                    if d[2] < 6: continue        # not enough relative precision: skip this sign choice
                    res.append((d[:2], cls(U - C, p), cls(U - A, p)))
    return res, len(pts)
def real_forms(A, C, dl):
    """signs (as class data at inf) of (l_r, l_m, l_*) at real points and sign choices"""
    res = []
    for u in ([Fr(3 * max(abs(A), abs(C)) + 1)] if dl > 0 else [Fr(min(abs(A), abs(C)), 3)]):
        P = (u * u - A * A) / dl; Q = (u * u - C * C) / dl
        assert P > 0 and Q > 0
        for s1 in (1, -1):
            for s2 in (1, -1):
                for s3 in (1, -1):
                    U = s1 * u
                    sr = s2 if s2 != -s3 * 1 and False else None
                    # sign of s2 sqrt(P) - s3 sqrt(Q)
                    if s2 == -s3: sr = s2
                    else: sr = s2 * (1 if P > Q else -1)
                    res.append((sr, 1 if U - C > 0 else -1, 1 if U - A > 0 else -1))
    return res
def hil_inf(x, b): return -1 if (x < 0 and b < 0) else 1
PERMS = {'id': (0, 1, 2), 'swap m*': (0, 2, 1), 'swap rm': (1, 0, 2), 'swap r*': (2, 1, 0), 'cyc1': (1, 2, 0), 'cyc2': (2, 0, 1)}
def sqf(n):
    s = -1 if n < 0 else 1; r = 1
    for q, e in factorint(abs(n)).items():
        if e % 2: r *= q
    return s * r
def pairing_all(A, C, dl, bs, seed=1):
    """for each b in bs (b = (b_r, b_m)) and each perm: the global value, or 'nonconst@p'.  b_* = b_r b_m."""
    rng = random.Random(seed)
    B2 = C * C - A * A
    primes = set([2])
    for x in [A, C, B2, dl] + [y for b in bs for y in b]: primes |= set(factorint(abs(x)))
    out = {(b, pn): 1 for b in bs for pn in PERMS}
    loc = {}
    for p in sorted(primes):
        F, npts = local_forms(A, C, dl, p, rng)
        if not F:
            for k in out: out[k] = 'nopoint@%d' % p
            continue
        for b in bs:
            br, bm = b; bstar = sqf(br * bm); bb = (cls(br, p), cls(bm, p), cls(bstar, p))
            for pn, pm in PERMS.items():
                vals = set()
                for f in F:
                    h = 1
                    for i in range(3): h *= hil_vu(f[i], bb[pm[i]], p)
                    vals.add(h)
                if len(vals) != 1:
                    if not isinstance(out[(b, pn)], str): out[(b, pn)] = 'nonconst@%d' % p
                elif not isinstance(out[(b, pn)], str):
                    out[(b, pn)] *= vals.pop()
                if pn == 'id': loc[(b, p)] = vals.pop() if len(vals) == 1 else None
    R = real_forms(A, C, dl)
    for b in bs:
        br, bm = b; bstar = sqf(br * bm); bb = (br, bm, bstar)
        for pn, pm in PERMS.items():
            vals = set()
            for f in R:
                h = 1
                for i in range(3): h *= hil_inf(f[i], bb[pm[i]])
                vals.add(h)
            if isinstance(out[(b, pn)], str): continue
            if len(vals) != 1: out[(b, pn)] = 'nonconst@inf'
            else: out[(b, pn)] *= vals.pop()
    return out, loc
