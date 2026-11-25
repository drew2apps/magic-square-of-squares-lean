#!/usr/bin/env python3
"""
indep_second.py -- a second descent via 2-isogeny (EXPERIMENT), for the curve of (609,200)/(809,409), whose
2-Selmer bound is 4 (indep_rank.py): exact arithmetic, written here.

E: y^2 = x(x^2 + A x + B), alpha: E(Q) -> Q*/Q*^2, x -> x, and 2^r = #alpha(E(Q)) #alpha'(E'(Q)) / 4 (Silverman-Tate
III.6), E': y^2 = x(x^2 - 2A x + (A^2 - 4B)).  For d in the Selmer set of alpha, a rational point of
C_d: N^2 = d M^4 + A M^2 e^2 + (B/d) e^4 gives, with X = 2d M^2 + A e^2, Y = e^2, Z = 2N, a point of the conic
        X^2 - D1 Y^2 = d Z^2,        D1 = A^2 - 4B
(4d N^2 = (2dM^2 + Ae^2)^2 - D1 e^4).  From a rational point of the conic, every rational point is
mu (Qx, Qy, Qz)(p, q) with binary quadratic forms Qx, Qy, Qz and mu in Q* (lines through the point).  Then Y = e^2 and
X - A Y = 2d M^2 say: mu Qy(p, q) is a square and 2 d mu (Qx - A Qy)(p, q) is a square.  So d is in alpha(E(Q)) only
if for some squarefree mu the curve
        D_{d,mu}:  mu Qy(p, q) = rho^2,   2 d mu (Qx - A Qy)(p, q) = sigma^2
has a point over every Q_l.  If l divides mu to an odd power then l divides Qy(p,q) and (Qx - A Qy)(p,q)
(or 2d), hence their resultant (or the contents), for (p, q) primitive: so mu runs over the signed products of
the primes of 2 d Res(Qy, Qx - A Qy) cont(Qy) cont(Qx - A Qy).  Local test: residue classes of (p : q) in
P^1(Z_l), both square classes decided as in indep_rank.class_search; at large l only exclusions by odd
valuation (every residue that is not a root gives an odd valuation) are used, the rest counts as solvable.
An excluded d carries a certificate: for each mu, a prime with no local point (re-verified).

Then #alpha(E(Q)) <= the largest power of 2 not above the number of d that survive (alpha(E(Q)) is a group
inside that set).  The same for E'.  Soundness tests: the images of points found by PARI on random curves always
survive; the bound is never below PARI's lower bound.

Result (indep_second.txt): for the curve of (609,200)/(809,409) nothing is removed: all 4 and all 16 elements of
the two Selmer sets survive, and the bound stays 4.  This is what one should expect, and it is a consistency check
rather than a failure: the curves D_{d,mu} are 2-coverings of the phi-coverings C_d, and since phihat o phi = [2]
they are 2-coverings of E, so this "second descent via 2-isogeny" sees no more than the full 2-descent of
indep_rank.py (dim Sel^2 = 5, bound 4).  Going below 4 needs a 4-descent or the Cassels-Tate pairing, which is what
PARI's ellrank (s = 4) and Magma's TwoPowerIsogenyDescentRankBound do.  Nothing here proves a rank.

Usage: python3 indep_second.py              (the curve of magma_second_ranks.json for (609,200); indep_second.txt)
       python3 indep_second.py --selftest   (indep_second_selftest.txt)
"""
import sys, json, random
from fractions import Fraction
from math import isqrt
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import indep_rank as R

pari = R.pari
vp, legendre = R.vp, R.legendre

def sq_class_det(c, w, p):
    """c = value at the centre, w = min valuation of the other Taylor coefficients: is the square class of every
    value in the class determined?  Returns True (all squares), False (all non-squares), None (undetermined)."""
    if c == 0:
        return None
    v0 = vp(c, p)
    need = 3 if p == 2 else 1
    if v0 + need <= w:
        return R.is_square_Qp_unit_class(c, p)
    return None

def two_poly_search(g1, g2, p, start, maxdepth=30):
    """Is there x in the union of classes `start` with g1(x) and g2(x) both squares in Q_p?
    True (found), False (proved none), None (undecided).  Small p only (all residues)."""
    stack = list(start)
    undecided = False
    while stack:
        x0, k = stack.pop()
        c1 = R.poly_shift(g1, x0, p ** k)
        c2 = R.poly_shift(g2, x0, p ** k)
        w1 = min((vp(c, p) for c in c1[1:] if c), default=10 ** 9)
        w2 = min((vp(c, p) for c in c2[1:] if c), default=10 ** 9)
        s1 = sq_class_det(c1[0], w1, p)
        s2 = sq_class_det(c2[0], w2, p)
        if s1 is False or s2 is False:
            continue
        if s1 is True and s2 is True:
            return True
        if c1[0] == 0 and (c2[0] == 0 or R.is_square_Qp_unit_class(c2[0], p)):
            return True          # an exact point with rho = 0
        if c2[0] == 0 and R.is_square_Qp_unit_class(c1[0], p):
            return True
        if k >= maxdepth:
            undecided = True
            continue
        for dd in range(p):
            stack.append((x0 + dd * p ** k, k + 1))
    return None if undecided else False

def _const_times_square(H, p):
    """If H (mod p, nonzero) is c * h^2 with c a constant, return (c, h); else None."""
    lead = H[-1]
    monic = [x * pow(lead, -1, p) % p for x in H]
    if len(monic) == 1:
        return lead, [1]
    h = R._poly_sqrt_monic(monic, p)
    return (lead, h) if h is not None else None

def large_prime_search(g1, g2, p, start, maxdepth=12):
    """For a large odd p: is there x in the classes `start` with g1(x), g2(x) both squares in Q_p?
    False only when proved: at each class, the residues that are not roots of the reductions H1, H2 are all bad
    when a content valuation is odd, or H1 or H2 is a non-residue constant times a square, or H1 H2 is a
    non-residue constant times a square; then only the roots of H1 and H2 (found exactly) are refined.
    Otherwise True (by Weil's bound a good residue exists; this answer is on the safe side, not a proof)."""
    stack = list(start)
    undecided = False
    while stack:
        x0, k = stack.pop()
        c1 = R.poly_shift(g1, x0, p ** k)
        c2 = R.poly_shift(g2, x0, p ** k)
        if not any(c1) or not any(c2):
            return True
        m1 = min(vp(c, p) for c in c1 if c)
        m2 = min(vp(c, p) for c in c2 if c)
        H1 = R._pmod([c // p ** m1 for c in c1], p)
        H2 = R._pmod([c // p ** m2 for c in c2], p)
        none_good = False
        if m1 % 2 or m2 % 2:
            none_good = True
        else:
            for H in (H1, H2):
                cs = _const_times_square(H, p)
                if cs is not None and legendre(cs[0], p) == -1:
                    none_good = True
            if not none_good:
                cs = _const_times_square(R._pmul(H1, H2, p), p)
                if cs is not None and legendre(cs[0], p) == -1:
                    none_good = True
        if not none_good:
            return True
        rts = set(R.roots_mod_p(H1, p)) | set(R.roots_mod_p(H2, p))
        if rts and k >= maxdepth:
            undecided = True
            continue
        for t in rts:
            stack.append((x0 + t * p ** k, k + 1))
    return None if undecided else False

def large_prime_excludes(g1, g2, p):
    r1 = large_prime_search(g1, g2, p, [(0, 0)])
    if r1 is not False:
        return False
    h1, h2 = list(reversed(g1)), list(reversed(g2))
    r2 = large_prime_search(h1, h2, p, [(0, 1)])
    return r2 is False

def local_solvable(g1, g2, p):
    """Y1^2 = g1(x), Y2^2 = g2(x) with (x : 1) or (1 : z), z in p Z_p (g1, g2 of degree <= 2: binary quadratic forms,
    so the chart x -> 1/x multiplies them by x^2, a square).  False only if proved."""
    if p >= R.LIMIT:
        return not large_prime_excludes(g1, g2, p)
    r1 = two_poly_search(g1, g2, p, [(0, 0)])
    if r1 is True:
        return True
    h1, h2 = list(reversed(g1 + [0] * (3 - len(g1)))), list(reversed(g2 + [0] * (3 - len(g2))))
    r2 = two_poly_search(h1, h2, p, [(0, 1)])
    if r2 is True:
        return True
    if r1 is None or r2 is None:
        return True
    return False

# ------------------------------------------------------------- conics
def conic_point(D1, d):
    """A rational point (X, Y, Z) of X^2 - D1 Y^2 = d Z^2, for D1 = -c^2 (d a sum of two squares) or D1 = s^2."""
    s = isqrt(abs(D1))
    if D1 > 0 and s * s == D1:
        # (X - sY)(X + sY) = d Z^2: X - sY = d, X + sY = 1
        X = Fraction(d + 1, 2); Y = Fraction(1 - d, 2 * s)
        return (X, Y, Fraction(1))
    assert D1 < 0 and s * s == -D1
    assert d > 0
    # d = a^2 + b^2 from its prime factors (d squarefree, primes 2 or 1 mod 4)
    g = (1, 0)
    for q in R.factor_proved(d):
        assert q == 2 or q % 4 == 1, ("d is not a sum of two squares", d, q)
        gq = (1, 1) if q == 2 else R.gaussian_prime_of(q)
        g = R._gmul(g, gq)
    a, b = g
    assert a * a + b * b == d
    return (Fraction(a), Fraction(b, s), Fraction(1))

def param_conic(D1, d, P0):
    """Binary quadratic forms (Qx, Qy, Qz) with F(Qx, Qy, Qz) = 0 identically, F = X^2 - D1 Y^2 - d Z^2, covering
    every rational point up to a scalar: P = F(V) P0 - 2 B(P0, V) V with V = (p, q, 0) (Z = 0 misses P0 if z0 != 0)."""
    x0, y0, z0 = P0
    assert z0 != 0 and x0 * x0 - D1 * y0 * y0 - d * z0 * z0 == 0
    # V = (p, q, 0): F(V) = p^2 - D1 q^2; B(P0, V) = x0 p - D1 y0 q
    # coefficient lists [c_pp, c_pq, c_qq]
    FV = [Fraction(1), Fraction(0), Fraction(-D1)]
    Bv = [Fraction(x0), Fraction(-D1 * y0)]            # linear: x0 p - D1 y0 q
    def lin_times_lin(l1, l2):
        return [l1[0] * l2[0], l1[0] * l2[1] + l1[1] * l2[0], l1[1] * l2[1]]
    Qx = [FV[i] * x0 for i in range(3)]
    t = lin_times_lin(Bv, [Fraction(1), Fraction(0)])   # B * p
    Qx = [Qx[i] - 2 * t[i] for i in range(3)]
    Qy = [FV[i] * y0 for i in range(3)]
    t = lin_times_lin(Bv, [Fraction(0), Fraction(1)])   # B * q
    Qy = [Qy[i] - 2 * t[i] for i in range(3)]
    Qz = [FV[i] * z0 for i in range(3)]
    # check the identity at a few points
    for (p, q) in [(1, 0), (0, 1), (1, 1), (2, -3), (5, 7)]:
        ev = lambda Q: Q[0] * p * p + Q[1] * p * q + Q[2] * q * q
        assert ev(Qx) ** 2 - D1 * ev(Qy) ** 2 - d * ev(Qz) ** 2 == 0
    return Qx, Qy, Qz

def to_int_forms(*forms):
    """Multiply all forms by the same square L^2 so that they are integral."""
    from math import lcm
    L = 1
    for f in forms:
        for c in f:
            L = lcm(L, Fraction(c).denominator)
    return [[int(Fraction(c) * L * L) for c in f] for f in forms]

def resultant2(f, g):
    """Resultant of the binary quadratics f = [a0, a1, a2] (a0 p^2 + a1 pq + a2 q^2) and g."""
    M = pari.matrix(4, 4, [f[0], f[1], f[2], 0, 0, f[0], f[1], f[2], g[0], g[1], g[2], 0, 0, g[0], g[1], g[2]])
    return int(pari.matdet(M))

def second_descent_d(A, B, d, extra_primes=()):
    """Does d (in the Selmer set of alpha on y^2 = x(x^2 + A x + B)) survive the second descent?
    Returns (True, mu) or (False, {mu: prime})."""
    D1 = A * A - 4 * B
    sD = isqrt(abs(D1))
    if sD * sD != abs(D1):
        return True, None        # no conic point found by this program: d is kept (the safe side)
    P0 = conic_point(D1, d)
    Qx, Qy, Qz = param_conic(D1, d, P0)
    Qx, Qy, Qz = to_int_forms(Qx, Qy, Qz)
    G1 = Qy
    G2 = [Qx[i] - A * Qy[i] for i in range(3)]
    Res = resultant2(G1, G2)
    assert Res != 0
    cont = lambda f: abs(int(pari.content(pari.Vec(f))))
    support = sorted(set(R.factor_proved(2 * d * Res * cont(G1) * cont(G2))))
    tests = sorted(set(support) | set(extra_primes) | set(q for q in range(3, 60) if all(q % r for r in range(2, q))))
    mus = R.squarefree_divisors(support)
    certs = {}
    for mu in mus:
        # g1(x) = mu G1(x, 1), g2(x) = 2 d mu G2(x, 1), increasing degree in x = p/q
        g1 = [mu * G1[2], mu * G1[1], mu * G1[0]]
        g2 = [2 * d * mu * G2[2], 2 * d * mu * G2[1], 2 * d * mu * G2[0]]
        bad = None
        for l in tests:
            if not local_solvable(g1, g2, l):
                bad = l
                break
        if bad is None:
            return True, mu
        certs[mu] = bad
    # re-verify the certificates
    for mu, l in certs.items():
        g1 = [mu * G1[2], mu * G1[1], mu * G1[0]]
        g2 = [2 * d * mu * G2[2], 2 * d * mu * G2[1], 2 * d * mu * G2[0]]
        assert not local_solvable(g1, g2, l)
    return False, certs

def second_bound(A, B):
    """Rank bound after the second descent on both sides of the 2-isogeny with kernel (0, 0)."""
    A1, B1 = -2 * A, A * A - 4 * B
    bad = sorted(set(R.factor_proved(2 * B * B1)))
    out = {}
    sizes = []
    for (a, b, name) in ((A, B, "phi"), (A1, B1, "phihat")):
        sel, _ = R.selmer_phi(a, b, bad)
        surv = []
        excl = {}
        for d in sel:
            ok, info = second_descent_d(a, b, d)
            if ok:
                surv.append(d)
            else:
                excl[str(d)] = {str(k): v for k, v in info.items()}
        size = 1 << (len(surv).bit_length() - 1)
        sizes.append(size)
        out[name] = {"selmer": sel, "survive": surv, "excluded": excl, "group_bound": size}
    prod = sizes[0] * sizes[1]
    rbound = max(prod.bit_length() - 1 - 2, 0) if prod >= 4 else 0
    out["rank_bound"] = rbound if prod >= 4 else "inconsistent: fewer than 4"
    return out

def selftest(n=12):
    rnd = random.Random(21)
    done = 0
    nimg = 0
    while done < n:
        # A = 2a, B = a^2 + b^2 = k^2 (so that both conics have a point found here): a, b, k Pythagorean
        s_, t_ = rnd.randint(1, 6), rnd.randint(1, 6)
        if s_ == t_:
            continue
        a, b = (s_ * s_ - t_ * t_) * rnd.choice([1, -1]), 2 * s_ * t_
        if rnd.random() < 0.5:
            a, b = b * rnd.choice([1, -1]), abs(a)
        mult = rnd.randint(1, 3)
        a, b = a * mult, b * mult
        A, B = 2 * a, a * a + b * b          # D1 = -4 b^2, B a square
        E = pari.ellinit([0, A, 0, B, 0])
        rk = pari.ellrank(E)
        lo = int(rk[0])
        out = second_bound(A, B)
        assert isinstance(out["rank_bound"], int) and out["rank_bound"] >= lo, (A, B, out["rank_bound"], lo)
        for P in rk[3]:
            x = Fraction(int(pari.numerator(P[0])), int(pari.denominator(P[0])))
            if x == 0:
                continue
            num = x.numerator * x.denominator
            sf = 1 if num > 0 else -1
            for q, e in R.factor_proved(abs(num)).items():
                if e % 2:
                    sf *= q
            assert sf in out["phi"]["survive"], ("an image of a point was excluded", A, B, P, sf)
            nimg += 1
        done += 1
    print("self-test: %d random curves y^2 = x((x + a)^2 + b^2): bound >= PARI's lower bound; the images of %d points"
          " found by PARI all survive the second descent" % (n, nimg))

def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--selftest":
        selftest()
        return
    curves = json.load(open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "magma_second_ranks.json")))
    c = [c for c in curves if c["pair"] == [609, 200]][0]
    a1, a2, a3, a4, a6 = [int(x) for x in c["ainvs"]]
    roots = R.rational_roots_cubic(a2, a4, a6)
    A, B = R.shifted(a2, a4, a6, roots[0])
    out = second_bound(A, B)
    print(json.dumps({"pair": c["pair"], "ainvs": c["ainvs"], "A": A, "B": B, "second_descent": out}))
    print("SUMMARY (609,200)/(809,409) curve: Selmer sets %d and %d; after the second descent %d and %d survive;"
          " rank bound %s" % (len(out["phi"]["selmer"]), len(out["phihat"]["selmer"]), len(out["phi"]["survive"]),
                              len(out["phihat"]["survive"]), out["rank_bound"]))

if __name__ == "__main__":
    main()
