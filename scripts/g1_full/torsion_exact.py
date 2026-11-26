#!/usr/bin/env python3
"""E(Q)_tors exactly, without PARI, for y^2 = x^3 + a2 x^2 + a4 x + a6 with integers a2, a4, a6 (a1 = a3 = 0).
Written by an AI model; not read by a person.  A library (census_lvalue_exact.py --torsion, rank_fisher_theta4_order.py);
`python3 torsion_exact.py` runs its self-test on curves of known torsion, `--negative` forged certificates that must be
refused (no output files; CI greps the last line).

QUOTED THEOREM.  Let p be an odd prime not dividing Delta of the model.  Reduction modulo p is a homomorphism
E(Q) -> E~(F_p), and it is injective on E(Q)_tors (Silverman, AEC, Prop. VII.3.1(b) for torsion of order prime to p;
for p-power torsion, E_1(Q_p) has no torsion because v(p) = 1 < p - 1 when p >= 3, AEC Thm IV.6.4(b) and
Prop. VII.2.2).  So #E(Q)_tors divides B = gcd #E~(F_p) over any set of such primes.

THE CERTIFICATE for "E(Q)_tors = T" (T a list of exact points), checked by certify():
  (a) every point of T is O or lies on the curve (exact); T is closed under the group law (exact, rational
      arithmetic) and contains O, so T is a finite subgroup of E(Q) and T is in E(Q)_tors; #T divides B.
  (b) completeness, for every prime l dividing B:
      (b1) v_l(#T) = v_l(B): the l-part of B is attained, so T[l^inf] = E(Q)[l^inf];
      (b2) only for l = 2, when v_2(#T) < v_2(B): E(Q)[2] is contained in T (the rational roots of the cubic are
           counted exactly: the roots found, then the cofactor has no rational root -- a linear cofactor is a root, a
           quadratic one has non-square discriminant, a cubic one has no root modulo some good prime), and for every
           R in T[2^inf] that is not 2S for any S in T, a good odd prime p with R~ not in 2E~(F_p) (all of E~(F_p)
           listed).  Proof that then E(Q)[2^inf] = T[2^inf]: otherwise take P in E(Q)[2^inf] outside T of least order;
           R = 2P has smaller order, so R is in T.  If R = 2S with S in T, then P - S is in E(Q)[2], which is in T, so P
           is in T: no.  So R has no half in T, and R~ = 2 P~ in 2E~(F_p) for every good p: contradicts the listed p.
      (b3) an odd l with v_l(#T) < v_l(B): refused (no claim).
  Primes l not dividing B do not divide #E(Q)_tors.  So (a) and (b) give E(Q)_tors = T.
FINDING T (a search; its correctness is not assumed, since certify() refuses an incomplete T): E[2](Q) from the
rational roots of the cubic, then repeated halving (the rational roots of x(2P) = x(R), sympy's exact factorisation),
and for odd l | B the rational roots of the division polynomial psi_n, n the l-part of B.
"""
import math, sys
from fractions import Fraction as Fr

# ---------------------------------------------------------------- the curve and its group law (exact)
class Curve:
    def __init__(self, a2, a4, a6):
        self.a2, self.a4, self.a6 = int(a2), int(a4), int(a6)
        b2, b4, b6 = 4*self.a2, 2*self.a4, 4*self.a6; b8 = 4*self.a2*self.a6 - self.a4*self.a4
        self.D = -b2*b2*b8 - 8*b4**3 - 27*b6*b6 + 9*b2*b4*b6
        assert self.D != 0, 'singular'
    def f(self, x): return x**3 + self.a2*x*x + self.a4*x + self.a6
    def on(self, P): return P is None or Fr(P[1])**2 == self.f(Fr(P[0]))
    def neg(self, P): return None if P is None else (P[0], -P[1])
    def add(self, P, Q):
        if P is None: return Q
        if Q is None: return P
        x1, y1 = Fr(P[0]), Fr(P[1]); x2, y2 = Fr(Q[0]), Fr(Q[1])
        if x1 == x2:
            if y1 + y2 == 0: return None
            lam = (3*x1*x1 + 2*self.a2*x1 + self.a4) / (2*y1)
        else: lam = (y2 - y1) / (x2 - x1)
        x3 = lam*lam - self.a2 - x1 - x2
        return (x3, -(y1 + lam*(x3 - x1)))
    def mul(self, k, P):
        R = None; Q = P
        if k < 0: k, Q = -k, self.neg(P)
        while k:
            if k & 1: R = self.add(R, Q)
            Q = self.add(Q, Q); k >>= 1
        return R
    def order(self, P, bound):
        """least k <= bound with kP = O, else None"""
        Q = P
        for k in range(1, bound + 1):
            if Q is None: return k
            Q = self.add(Q, P)
        return None
    def model(self): return (0, self.a2, 0, self.a4, self.a6)

def key(P): return None if P is None else (Fr(P[0]), Fr(P[1]))
def show(P): return 'O' if P is None else '(%s, %s)' % (P[0], P[1])

# ---------------------------------------------------------------- reduction modulo p (exact, brute force)
def primes_upto(M):
    s = bytearray([1]) * (M + 1); s[0:2] = b'\x00\x00'
    for i in range(2, math.isqrt(M) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(M + 1) if s[i]]

def good_odd_primes(E, limit=1000):
    return [p for p in primes_upto(limit) if p > 2 and E.D % p]

def points_mod(E, p):
    """all affine points of E~(F_p) (p odd, good)"""
    roots = {}
    for y in range(p): roots.setdefault(y*y % p, []).append(y)
    return [(x, y) for x in range(p) for y in roots.get((x**3 + E.a2*x*x + E.a4*x + E.a6) % p, [])]

def count_mod(E, p): return len(points_mod(E, p)) + 1

def add_mod(E, P, Q, p):
    if P is None: return Q
    if Q is None: return P
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2:
        if (y1 + y2) % p == 0: return None
        lam = (3*x1*x1 + 2*E.a2*x1 + E.a4) * pow(2*y1, -1, p) % p
    else: lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (lam*lam - E.a2 - x1 - x2) % p
    return (x3, (-(y1 + lam*(x3 - x1))) % p)

def reduce_pt(P, p):
    """P~ in E~(F_p); None if P = O or P reduces to O (p divides a denominator)"""
    if P is None: return None
    x, y = Fr(P[0]), Fr(P[1])
    if x.denominator % p == 0 or y.denominator % p == 0: return None
    return (x.numerator * pow(x.denominator, -1, p) % p, y.numerator * pow(y.denominator, -1, p) % p)

def doubles_mod(E, p):
    pts = points_mod(E, p)
    return {add_mod(E, P, P, p) for P in pts} | {None}

def bound(E, primes):
    """B = gcd #E~(F_p) over the given primes, each asserted an odd prime (trial division) of good reduction; the list
    must not be empty (an empty list would give B = 0, which every #T divides)"""
    assert len(primes) > 0, ('no primes given',)
    g = 0
    for p in primes:
        assert p > 2 and all(p % q for q in range(2, math.isqrt(p) + 1)), ('not an odd prime', p)
        assert E.D % p != 0, ('not a good odd prime', p)
        g = math.gcd(g, count_mod(E, p))
    return g

# ---------------------------------------------------------------- finding torsion points (a search; certified later)
def rational_roots(coeffs):
    """rational roots of sum coeffs[i] x^i (exact coefficients), via sympy's factorisation over Q (a finder)"""
    import sympy as sp
    x = sp.Symbol('x')
    P = sp.Poly(sum(sp.Rational(Fr(c).numerator, Fr(c).denominator) * x**i for i, c in enumerate(coeffs)), x)
    if P.is_zero: return []
    out = set()
    for fac, _ in P.factor_list()[1]:
        if fac.degree() == 1:
            c1, c0 = fac.all_coeffs(); r = -sp.Rational(c0) / sp.Rational(c1); out.add(Fr(int(r.p), int(r.q)))
    return sorted(out)

def points_with_x(E, xs):
    out = []
    for x in xs:
        v = E.f(Fr(x))
        if v == 0: out.append((Fr(x), Fr(0))); continue
        if v > 0:
            a, b = math.isqrt(v.numerator), math.isqrt(v.denominator)
            if a*a == v.numerator and b*b == v.denominator: out += [(Fr(x), Fr(a, b)), (Fr(x), -Fr(a, b))]
    return out

def two_torsion(E):
    return [None] + points_with_x(E, rational_roots([E.a6, E.a4, E.a2, 1]))

def halves(E, R):
    """all rational P with 2P = R: x(2P) = (x^4 - b4 x^2 - 2 b6 x - b8) / (4 f(x)), so the roots of
    x^4 - b4 x^2 - 2 b6 x - b8 - 4 x_R f(x); each candidate is checked with the group law"""
    if R is None: return two_torsion(E)
    b4, b6, b8 = 2*E.a4, 4*E.a6, 4*E.a2*E.a6 - E.a4*E.a4; xr = Fr(R[0])
    co = [-b8 - 4*xr*E.a6, -2*b6 - 4*xr*E.a4, -b4 - 4*xr*E.a2, -4*xr, 1]
    return [P for P in points_with_x(E, rational_roots(co)) if key(E.add(P, P)) == key(R)]

def division_poly(E, n):
    """psi_n for odd n as a sympy polynomial in x (recursion with psi_2 = 2y, (2y)^2 = 4 f)"""
    import sympy as sp
    x = sp.Symbol('x'); f = sp.Poly(x**3 + E.a2*x**2 + E.a4*x + E.a6, x)
    b2, b4, b6, b8 = 4*E.a2, 2*E.a4, 4*E.a6, 4*E.a2*E.a6 - E.a4**2
    one = sp.Poly(1, x)
    psi = {0: (sp.Poly(0, x), 0), 1: (one, 0), 2: (one, 1),
           3: (sp.Poly(3*x**4 + b2*x**3 + 3*b4*x**2 + 3*b6*x + b8, x), 0),
           4: (sp.Poly(2*x**6 + b2*x**5 + 5*b4*x**4 + 10*b6*x**3 + 10*b8*x**2 + (b2*b8 - b4*b6)*x + b4*b8 - b6**2, x), 1)}
    def mul(A, B):
        e = A[1] + B[1]; P = A[0] * B[0]
        for _ in range(e // 2): P = P * 4 * f
        return (P, e % 2)
    def sub(A, B): assert A[1] == B[1] or A[0].is_zero or B[0].is_zero; return (A[0] - B[0], max(A[1], B[1]))
    def get(k):
        if k in psi: return psi[k]
        m = k // 2
        if k % 2:
            r = sub(mul(get(m + 2), mul(get(m), mul(get(m), get(m)))), mul(get(m - 1), mul(get(m + 1), mul(get(m + 1), get(m + 1)))))
        else:
            br = sub(mul(get(m + 2), mul(get(m - 1), get(m - 1))), mul(get(m - 2), mul(get(m + 1), get(m + 1))))
            r = mul(get(m), br)                      # psi_m * bracket, then divide by psi_2 = (2y)
            if r[1] == 1: r = (r[0], 0)
            else:
                q, rem = sp.div(r[0], 4 * f); assert rem.is_zero; r = (q, 1)
        psi[k] = r; return r
    P, e = get(n); assert e == 0
    return [Fr(int(sp.Rational(c).p), int(sp.Rational(c).q)) for c in reversed(P.all_coeffs())]

def find_torsion(E, B):
    """a list of torsion points (a search): the 2-part by repeated halving, the odd part from psi_n"""
    T2 = {key(P): P for P in two_torsion(E)}
    for _ in range(max(0, (B & -B).bit_length() - 1)):
        new = {}
        for R in list(T2.values()):
            for P in halves(E, R): new.setdefault(key(P), P)
        if all(k in T2 for k in new): break
        T2.update(new)
    odd = B
    while odd % 2 == 0: odd //= 2
    To = {None: None}
    if odd > 1:
        for P in points_with_x(E, rational_roots(division_poly(E, odd))): To[key(P)] = P
    T = {}
    for P in T2.values():
        for Q in To.values():
            S = E.add(P, Q); T[key(S)] = S
    return list(T.values())

# ---------------------------------------------------------------- the certificate
def vl(n, l):
    k = 0
    while n % l == 0: n //= l; k += 1
    return k

def factor_small(n):
    out = []; d = 2
    while d * d <= n:
        while n % d == 0: out.append(d); n //= d
        d += 1
    if n > 1: out.append(n)
    return sorted(set(out))

def structure(E, T):
    """the abelian group T as Z/n1 x Z/n2 (n1 | n2), from the element orders"""
    N = len(T); orders = [E.order(P, N) for P in T]; e = max(orders)
    return (N // e, e) if N // e > 1 else (e,)

def two_torsion_exact(E, T, primes):
    """E(Q)[2] is contained in T: the 2-torsion roots in T, then the cofactor of the cubic has no rational root"""
    rts = sorted({Fr(P[0]) for P in T if P is not None and P[1] == 0})
    co = [Fr(E.a6), Fr(E.a4), Fr(E.a2), Fr(1)]
    for r in rts:                                            # synthetic division by (x - r), remainder must be 0
        q = [Fr(0)] * (len(co) - 1); acc = Fr(0)
        for i in range(len(co) - 1, 0, -1): acc = acc * r + co[i]; q[i - 1] = acc
        assert acc * r + co[0] == 0, ('a listed 2-torsion point is not a root of the cubic', r)
        co = q
    deg = len(co) - 1
    if deg == 0: return True, 'all three roots of the cubic are listed'
    if deg == 1: return False, 'the cofactor is linear (a rational root not listed)'
    if deg == 2:
        c0, c1, c2 = co; d = c1*c1 - 4*c0*c2
        sq = d >= 0 and math.isqrt(d.numerator)**2 == d.numerator and math.isqrt(d.denominator)**2 == d.denominator
        return (not sq), ('the quadratic cofactor has discriminant %s, %s' % (d, 'a square (NOT EXCLUDED)' if sq else 'not a square'))
    for p in primes:                                         # deg 3: a prime with no root mod p
        if all((x**3 + E.a2*x*x + E.a4*x + E.a6) % p for x in range(p)):
            return True, 'the cubic has no root modulo %d' % p
    return False, 'no prime found where the cubic has no root'

def certify(E, T, primes, ndouble=40):
    """the certificate (a), (b) of the module docstring; returns dict(ok, B, structure, reasons); never raises on a
    forged input (it is refused with a reason)"""
    r = dict(ok=False, reasons=[])
    try: B = bound(E, primes)
    except AssertionError as ex: r['reasons'].append('bound: %s' % (ex.args[0],)); return r
    r['B'] = B
    for P in T:
        if not E.on(P): r['reasons'].append('%s is not on the curve' % show(P)); return r
    ks = {key(P) for P in T}
    if None not in ks: r['reasons'].append('O is not listed'); return r
    for P in T:
        for Q in T:
            if key(E.add(P, Q)) not in ks: r['reasons'].append('not closed: %s + %s is not listed' % (show(P), show(Q))); return r
    if B % len(ks):
        r['reasons'].append('#T = %d does not divide B = %d' % (len(ks), B)); return r
    T = [P for P in T if key(P) in ks]; T = list({key(P): P for P in T}.values())
    r['structure'] = structure(E, T); r['order'] = len(T)
    for l in factor_small(B):
        if vl(len(T), l) == vl(B, l): r['reasons'].append('%d-part of B attained' % l); continue
        if l != 2: r['reasons'].append('odd %d-part of B (%d) not attained by T (%d): NOT CERTIFIED' % (l, l**vl(B, l), l**vl(len(T), l))); return r
        ok2, why = two_torsion_exact(E, T, primes)
        if not ok2: r['reasons'].append('E(Q)[2] not shown to lie in T: %s' % why); return r
        T2 = [P for P in T if (lambda o: o & (o - 1) == 0)(E.order(P, len(T)))]      # T[2^inf]
        dbl = {key(E.add(S, S)) for S in T}
        lonely = [R for R in T2 if key(R) not in dbl]
        cache = {}; excl = []
        for R in lonely:
            hit = None
            for p in primes[:ndouble]:
                if p not in cache: cache[p] = doubles_mod(E, p)
                Rp = reduce_pt(R, p)
                if Rp is not None and Rp not in cache[p]: hit = p; break
            if hit is None: r['reasons'].append('%s has no half in T, and it is a double modulo each of the first %d good primes: NOT EXCLUDED' % (show(R), ndouble)); return r
            excl.append((R, hit))
        r['reasons'].append('2-part of B (%d) not attained by T (%d): E(Q)[2] in T (%s); the %d points of T[2^inf] with no half in T are not doubles modulo %s'
                            % (2**vl(B, 2), 2**vl(len(T), 2), why, len(excl), sorted({p for _, p in excl})))
        r['excluded'] = excl
    r['ok'] = True
    return r

def fmt_structure(s): return ' x '.join('Z/%d' % k for k in s) if s != (1,) else 'trivial'

def torsion(E, primes):
    """find, then certify; returns (T, certificate)"""
    B = bound(E, primes)
    T = find_torsion(E, B)
    return T, certify(E, T, primes)

# ---------------------------------------------------------------- self-test on curves of known torsion
CONTROLS = [  # (label, (a2, a4, a6) of y^2 = x^3 + b2 x^2 + 8 b4 x + 16 b6 for Cremona's curve, torsion from the tables)
    ('15a1', (5, -152, -624), (2, 4)), ('15a4', (5, 568, -1776), (8,)), ('15a8', (5, -1272, 15504), (4,)),
    ('210e2', (1, -17120, 499968), (2, 8)), ('14a1', (1, 72, -368), (6,)), ('30a2', (1, -296, 1680), (2, 6)),
    ('90c3', (-3, -1944, 110160), (12,)), ('11a1', (-4, -160, -1264), (5,)), ('26b1', (-3, -40, 208), (7,)),
    ('54b3', (-3, -216, 1872), (9,)), ('66c1', (1, -720, 5184), (10,))]

def selftest():
    bad = 0
    for lab, (a2, a4, a6), want in CONTROLS:
        E = Curve(a2, a4, a6); T, c = torsion(E, good_odd_primes(E, 400))
        got = c.get('structure')
        print('%-6s %-28s B = %-3s found %-14s certified %s  (tables: %s)' % (lab, E.model(), c.get('B'), fmt_structure(got) if got else '-',
              c['ok'], fmt_structure(want)))
        bad += not (c['ok'] and got == want)
    E = Curve(-4, -125120, -16869104); T, c = torsion(E, good_odd_primes(E, 400))   # 11a2: trivial torsion, B = 5
    print('11a2   %-28s B = %s, found %d point(s): %s (expected: refused, the odd part of B is not attained)' % (E.model(), c['B'], len(T), c['reasons'][-1]))
    bad += c['ok']
    assert bad == 0
    print('self-test: all as expected')

def negative():
    """forged certificates on curves of known torsion: each must be refused"""
    ntot = nref = 0
    def chk(lab, E, T, pr):
        nonlocal ntot, nref
        c = certify(E, T, pr); ntot += 1; nref += not c['ok']
        print('forged %-62s %s' % (lab + ':', 'NOT REFUSED' if c['ok'] else 'refused (%s)' % c['reasons'][-1]))
    for lab, (a2, a4, a6), want in CONTROLS:
        if lab not in ('15a8', '210e2', '90c3'): continue
        E = Curve(a2, a4, a6); pr = good_odd_primes(E, 400); T, c = torsion(E, pr); assert c['ok'] and c['structure'] == want
        o = {key(P): E.order(P, len(T)) for P in T}; P = max(T, key=lambda Q: o[key(Q)])
        bad = next(p for p in range(3, 400, 2) if E.D % p == 0 and all(p % q for q in range(3, p, 2)))
        chk('%s: a point (x, y + 1) off the curve' % lab, E, T + [(P[0], P[1] + 1)], pr)
        chk('%s: a point of largest order left out' % lab, E, [Q for Q in T if Q != P], pr)
        chk('%s: the subgroup of points of order dividing 2' % lab, E, [Q for Q in T if o[key(Q)] <= 2], pr)
        chk('%s: O left out' % lab, E, [Q for Q in T if Q is not None], pr)
        chk('%s: the bad prime %d among the primes' % (lab, bad), E, T, sorted(pr[:20] + [bad]))
        chk('%s: p = 2 among the primes' % lab, E, T, [2] + pr[:20])
    E = Curve(-4, -125120, -16869104)
    chk('11a2 (trivial torsion, B = 5): E(Q)_tors claimed trivial', E, [None], good_odd_primes(E, 400))
    E = Curve(5, -152, -624)                                  # 15a1, torsion Z/2 x Z/4 (a review found these two)
    chk('15a1: E(Q)_tors claimed trivial, with no primes', E, [None], [])
    chk('15a1: E(Q)_tors claimed trivial, with the composites 21, 33', E, [None], [21, 33])
    E = Curve(0, -2, 0)
    chk('y^2 = x^3 - 2x: (2, 2) (infinite order) listed with O, (0, 0)', E, [None, (0, 0), (2, 2), (2, -2)], good_odd_primes(E, 400))
    print('negative tests: %d forged certificates, %d refused' % (ntot, nref))
    assert ntot == nref

if __name__ == '__main__':
    if '--negative' in sys.argv: negative()
    else: selftest()
