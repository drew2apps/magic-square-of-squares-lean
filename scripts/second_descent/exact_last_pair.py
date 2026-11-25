#!/usr/bin/env python3
"""Theorem A of LAST_PAIR.md for (627,364), rechecked in exact arithmetic with no Sage, no PARI and no Magma.

    python3 exact_last_pair.py              the whole check (about ten seconds with the stored hints); output exact_last_pair.txt
    python3 exact_last_pair.py --negative   forged inputs that must be refused; output exact_last_pair_negative.txt
    (add --write to rewrite the output files and the hints; without it the script only prints, so CI can diff its output)

Written by an AI model (Claude); not read by a person.  In the style of scripts/g1_full/sigma_exact.py: Python
integers and fractions, Gaussian rationals as pairs of fractions, logarithms enclosed by the atanh series
(sigma_exact.ln_bounds), sympy for exact polynomial algebra over Q and Q(i) and to factor integers.  PARI is used only
as a finder of hints (conic points and points for the second descent, factorisations), stored in
exact_last_pair_hints.json; every hint is checked exactly before it is used, and the run with the stored hints does not
import PARI.

What is checked, in the order of LAST_PAIR.md:
 S1  the curves C1: Y^2 = d Q1(s), C2: Y^2 = d Q2(s); the maps psi_j: C_j -> E_j of weierstrass() (an identity of
     rational functions, checked at more points than its degree); the isomorphism of E1 with the stated minimal model;
     that E2 is isomorphic over Q to the twist of E1 by -1; the isomorphism theta: E2 -> E1 over K = Q(i) and
     sigma theta = -theta.  No complex multiplication: j(E1) is not an integer.
 S2  the point G: P + sigma P at t = 1 is the stated G (up to sign); the point B at t = 1 and the points at t = -1 are
     +-G + torsion.
 S3  torsion: E1(Q)_tors = {O, T}.
 S4  rank E1(Q) <= 1 by an exact descent via 2-isogeny and an exact second descent (scripts/g1_full/rank_indep.py).
 S5  the canonical height of G enclosed (Silverman's bounds, exact doubling); G has infinite order.
 S6  saturation: for every prime l <= 199, a functional E1(Q) -> Z/l kills l E1(Q) and the torsion and not G.
 S6b a lower bound for the heights of points of infinite order, without Sage: the model is minimal, 24 P lies in E0
     at every prime, and for such points h(x) - hhat <= 4.444 (Cremona-Prickett-Siksek, J. Number Theory 116 (2006):
     Theorem 2, no non-archimedean term for points in E0 at every prime, and Theorem 15, the upper side of Silverman's
     archimedean bound (Silverman 1990, p. 737) in their normalisation; the statement height_enclosure.py uses);
     a search of the x = a/e^2 with |a|, e^2 <= 110
     finds no point, so hhat(P) >= hhat(G)/210^2 for every P of infinite order, and with S6 G generates.  (Sage's
     bound 0.045162 of LAST_PAIR.md is printed for comparison and not used.)
 S7  Lemma 2 (the structure of the cover: the four f_i, their common roots, which products are squares) over Q(i).
 S8  Lemmas 4 and 5: every quantity of the height bound as an exact rational upper bound, from the constants of
     Checcoli-Veneziano-Viada (arXiv:1602.04097, Theorem 1.2 for N = 2, Proposition 3.2 and Table 1, Zhang's
     inequality as their Theorem 2.2), which are inputs (CVV below); the bound on a, b.
 S9  Lemma 6: the Mordell-Weil sieve as an exact computation of the residue sets modulo four primes, the lifting by
     the Chinese remainder theorem, and a test of the lifting with planted pairs.
 S10 Lemma 7: the fibre a = 1 as an exact polynomial equation over Q, its roots in Q(i), and the rational t;
     the same for b = 1 on E2 over Q(i) (not needed by the proof).
 S11 the final enumeration: the t of Theorem A are 1 and -1; t = 0 and t = infinity are not points.
The ratio (991,263) is the partner of (627,364) (checked here); that a solution for one gives one for the other is the
Lean theorem partner_forms and is not rechecked here.
"""
import os, sys, json, itertools
WRITE = '--write' in sys.argv   # outputs and hint files are written only with --write, so CI can compare stdout with them
from fractions import Fraction as Fr
from math import gcd, isqrt

here = os.path.dirname(os.path.abspath(__file__))
G1 = os.path.join(here, '..', 'g1_full')
sys.path.insert(0, G1)
import sigma_exact as SE          # ln_bounds, lnplus_hi, LN2 (rational enclosures of logarithms)
import rank_indep as RI           # the exact descent

HINTS_FILE = os.path.join(here, 'exact_last_pair_hints.json')

# ------------------------------------------------------------------ the data of LAST_PAIR.md
M_, N_ = 627, 364
D_CLASS = 521266
EM1 = (0, 1, 0, 15223208213202973031750280786182975, 372517508636065850025943055293341065999858640587423)
G_STATED = (Fr(83214074135015341), Fr(47069367709023802879687500))
LAM_SAGE = Fr(45162, 10**6)       # Sage's height_function().min for E1 (LAST_PAIR.md, Lemma 1); NOT recomputed here
SIEVE_PRIMES = (10009, 20089, 30013, 40009)
CLAIMED_AB = 14175                # LAST_PAIR.md, Lemma 5 (deg C <= 96)
XMAX = 50000                      # largest |a|, e^2 the height search of S6b may visit
CVV = dict(C1=Fr(72251, 1000), k6=Fr(6019, 1000), c1g=Fr(3724, 1000), c2g=Fr(4015, 1000), c1q=Fr(2919, 1000), c2q=Fr(321, 100))
DEG_C = 96

OUT = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.append(s)

# ------------------------------------------------------------------ Gaussian rationals
class GQ:
    __slots__ = ('re', 'im')
    def __init__(s, re, im=0): s.re = Fr(re); s.im = Fr(im)
    def __add__(s, o): o = gq(o); return GQ(s.re + o.re, s.im + o.im)
    __radd__ = __add__
    def __sub__(s, o): o = gq(o); return GQ(s.re - o.re, s.im - o.im)
    def __rsub__(s, o): return gq(o) - s
    def __neg__(s): return GQ(-s.re, -s.im)
    def __mul__(s, o): o = gq(o); return GQ(s.re * o.re - s.im * o.im, s.re * o.im + s.im * o.re)
    __rmul__ = __mul__
    def norm(s): return s.re * s.re + s.im * s.im
    def conj(s): return GQ(s.re, -s.im)
    def inv(s): nn = s.norm(); assert nn != 0; return GQ(s.re / nn, -s.im / nn)
    def __truediv__(s, o): return s * gq(o).inv()
    def __rtruediv__(s, o): return gq(o) * s.inv()
    def __pow__(s, k):
        r = GQ(1)
        for _ in range(k): r = r * s
        return r
    def __eq__(s, o): o = gq(o); return s.re == o.re and s.im == o.im
    def __hash__(s): return hash((s.re, s.im))
    def is_zero(s): return s.re == 0 and s.im == 0
    def __repr__(s): return '(%s + %s i)' % (s.re, s.im)
def gq(x): return x if isinstance(x, GQ) else GQ(x)
I = GQ(0, 1)

def issq_q(q):
    q = Fr(q)
    return q >= 0 and isqrt(q.numerator) ** 2 == q.numerator and isqrt(q.denominator) ** 2 == q.denominator
def sqrt_q(q): q = Fr(q); return Fr(isqrt(q.numerator), isqrt(q.denominator))
def gsqrt(z):
    """A square root of z in Q(i), or None.  (x + iy)^2 = a + bi:  x^2 = (a + w)/2, y^2 = (w - a)/2, w = |z|."""
    z = gq(z)
    if z.is_zero(): return GQ(0)
    nn = z.norm()
    if not issq_q(nn): return None
    w = sqrt_q(nn)
    for x2 in ((z.re + w) / 2, (z.re - w) / 2):
        if issq_q(x2):
            x = sqrt_q(x2)
            if x != 0:
                y = z.im / (2 * x)
                if GQ(x, y) * GQ(x, y) == z: return GQ(x, y)
            else:
                if issq_q(-z.re):
                    y = sqrt_q(-z.re)
                    if GQ(0, y) * GQ(0, y) == z: return GQ(0, y)
    return None

# ------------------------------------------------------------------ polynomials (low degree first) over Q or Q(i)
def pmul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): c[i + j] = c[i + j] + x * y
    return c
def peval(g, t):
    r = 0
    for c in reversed(g): r = r * t + c
    return r
def pder(g): return [i * c for i, c in enumerate(g)][1:]

def family(m, n):
    S, D = m * m + n * n, m * m - n * n
    qp = [S, -2 * D, S]
    Q1 = pmul(qp, pmul([m, -n], [n, m])); Q2 = pmul(qp, pmul([m, n], [n, -m]))
    return S, D, Q1, Fr(m, n), Q2, Fr(-m, n)

# ------------------------------------------------------------------ elliptic curves y^2 = x^3 + a2 x^2 + a4 x + a6
class Ell:
    def __init__(s, a2, a4, a6):
        s.a2, s.a4, s.a6 = a2, a4, a6
        s.b2, s.b4, s.b6 = 4 * a2, 2 * a4, 4 * a6
        s.b8 = 4 * a2 * a6 - a4 * a4
        s.D = -s.b2 ** 2 * s.b8 - 8 * s.b4 ** 3 - 27 * s.b6 ** 2 + 9 * s.b2 * s.b4 * s.b6
        s.c4 = s.b2 ** 2 - 24 * s.b4
        s.c6 = -s.b2 ** 3 + 36 * s.b2 * s.b4 - 216 * s.b6
        assert s.D != 0
    def f(s, x): return x ** 3 + s.a2 * x * x + s.a4 * x + s.a6
    def on(s, P): return P is None or P[1] * P[1] == s.f(P[0])
    def neg(s, P): return None if P is None else (P[0], -P[1])
    def add(s, P, Q):
        if P is None: return Q
        if Q is None: return P
        if P[0] == Q[0]:
            if (P[1] + Q[1]) == 0: return None
            l = (3 * P[0] * P[0] + 2 * s.a2 * P[0] + s.a4) / (2 * P[1])
        else:
            l = (Q[1] - P[1]) / (Q[0] - P[0])
        x = l * l - s.a2 - P[0] - Q[0]
        return (x, l * (P[0] - x) - P[1])
    def mul(s, k, P):
        R = None
        if k < 0: k, P = -k, s.neg(P)
        while k:
            if k & 1: R = s.add(R, P)
            P = s.add(P, P); k >>= 1
        return R
    def dup_x(s, x):
        num = x ** 4 - s.b4 * x * x - 2 * s.b6 * x - s.b8; den = 4 * x ** 3 + s.b2 * x * x + 2 * s.b4 * x + s.b6
        assert den != 0
        return num / den

def weierstrass(Q, r):
    """Y^2 = Q(s) with Q(r) = 0: s = r + 1/w, (Y w^2)^2 = A w^3 + B w^2 + C w + E; x = A w: y^2 = x^3 + B x^2 + AC x + A^2 E."""
    Qf = [Fr(c) for c in Q]
    d1 = pder(Qf); d2 = pder(d1); d3 = pder(d2); d4 = pder(d3)
    assert peval(Qf, r) == 0
    A, B, C, E = peval(d1, r), peval(d2, r) / 2, peval(d3, r) / 6, peval(d4, r) / 24
    assert A != 0
    return Ell(B, A * C, A * A * E), A

def psi(A, r, s, Y):
    """The point of the elliptic curve attached to (s, Y) on Y^2 = Q(s): x = A/(s - r), y = A Y/(s - r)^2."""
    w = 1 / (s - r)
    return (A * w, A * Y * w * w)

def check_psi(Q, r, E, A):
    """Q(s) (A/(s - r)^2)^2 = f(A/(s - r)) as rational functions: times (s - r)^6 both sides are polynomials of degree <= 6
    in s, so equality at 12 points proves it."""
    for k in range(1, 13):
        s = r + Fr(k, 7) + 1
        x = A / (s - r)
        assert peval(Q, s) * (A / (s - r) ** 2) ** 2 == E.f(x), 'psi is not a map to the curve'

def iso_params(Ef, Et):
    """(u, r) with x_f = u^2 x_t + r, y_f = u^3 y_t an isomorphism over Q from Et to Ef (both with a1 = a3 = 0), or None."""
    if Ef.c4 == 0 or Et.c4 == 0 or Ef.c6 == 0 or Et.c6 == 0: return None
    q4, q6 = Fr(Ef.c4, 1) / Et.c4, Fr(Ef.c6, 1) / Et.c6
    u2 = q6 / q4
    if u2 * u2 != q4 or u2 ** 3 != q6 or not issq_q(u2): return None
    u = sqrt_q(u2)
    r = (u2 * Et.a2 - Ef.a2) / 3
    if u2 ** 2 * Et.a4 != Ef.a4 + 2 * r * Ef.a2 + 3 * r * r: return None
    if u2 ** 3 * Et.a6 != Ef.a6 + r * Ef.a4 + r * r * Ef.a2 + r ** 3: return None
    return u, r

def to_model(P, u, r):
    """The image in Et of a point of Ef, for x_f = u^2 x_t + r, y_f = u^3 y_t."""
    return None if P is None else ((P[0] - r) / (u * u), P[1] / u ** 3)
def from_model(P, u, r):
    return None if P is None else (u * u * P[0] + r, u ** 3 * P[1])

# ------------------------------------------------------------------ setup (S1, S2)
def setup(m=M_, n=N_, d=D_CLASS, Gst=G_STATED):
    S, D, Q1, r1, Q2, r2 = family(m, n)
    assert isqrt(S) ** 2 == S
    dQ1 = [d * c for c in Q1]; dQ2 = [d * c for c in Q2]
    E1w, A1 = weierstrass(dQ1, r1); E2w, A2 = weierstrass(dQ2, r2)
    check_psi(dQ1, r1, E1w, A1); check_psi(dQ2, r2, E2w, A2)
    Em = Ell(*EM1[1::2][:1], *EM1[3:]) if False else Ell(EM1[1], EM1[3], EM1[4])
    ip = iso_params(E1w, Em); assert ip is not None, 'E1 is not isomorphic over Q to the stated minimal model'
    u1, rr1 = ip
    tw = Ell(-E1w.a2, E1w.a4, -E1w.a6)                    # the twist of E1 by -1:  -y^2 = f(x), x -> -x
    ip2 = iso_params(E2w, tw); assert ip2 is not None, 'E2 is not the twist of E1 by -1'
    u2, rr2 = ip2
    j = Fr(Em.c4 ** 3, Em.D); assert j.denominator != 1, 'j(E1) is an integer'
    # theta: E2w -> tw (over Q) -> E1w (over K): (x, y) -> (-x, i y)
    def theta(R):
        if R is None: return None
        X, Y = to_model(R, u2, rr2)
        return (gq(-X), I * Y)
    # degenerate points at t = 1 and t = -1 (s = -i t)
    def conjP(P): return None if P is None else (gq(P[0]).conj(), gq(P[1]).conj())
    E1K = E1w; E2K = E2w
    pts = {}
    for t in (1, -1):
        s = GQ(0, -t)
        Y1 = gsqrt(peval([gq(c) for c in dQ1], s)); Y2 = gsqrt(peval([gq(c) for c in dQ2], s))
        assert Y1 is not None and Y2 is not None, 'not a degenerate solution'
        P = psi(A1, gq(r1), s, Y1); R = psi(A2, gq(r2), s, Y2)
        assert P[1] * P[1] == E1K.f(P[0]) and R[1] * R[1] == E2K.f(R[0])
        Apt = E1K.add(P, conjP(P)); Bpt = theta(E2K.add(R, E2K.neg(conjP(R))))
        for Z in (Apt, Bpt):
            if Z is not None: assert gq(Z[0]).im == 0 and gq(Z[1]).im == 0, 'not rational'
        Ar = None if Apt is None else (gq(Apt[0]).re, gq(Apt[1]).re)
        Br = None if Bpt is None else (gq(Bpt[0]).re, gq(Bpt[1]).re)
        assert E1w.on(Ar) and E1w.on(Br)
        # sigma theta = - theta on a non-rational point: R itself
        th = theta(R); sth = (gq(th[0]).conj(), gq(th[1]).conj())
        tsr = theta(conjP(R))
        assert sth[0] == tsr[0] and sth[1] == -tsr[1], 'sigma theta != - theta'
        pts[t] = (to_model(Ar, u1, rr1), to_model(Br, u1, rr1))
    Gm = pts[1][0]
    assert Em.on(Gm) and Gm[0] == Gst[0] and abs(Gm[1]) == abs(Gst[1]), 'G is not P + conj(P) at t = 1'
    return dict(S=S, D=D, Q1=Q1, Q2=Q2, r1=r1, r2=r2, dQ1=dQ1, dQ2=dQ2, E1w=E1w, E2w=E2w, A1=A1, A2=A2, Em=Em,
                u1=u1, rr1=rr1, u2=u2, rr2=rr2, j=j, pts=pts, G=Gst, d=d, m=m, n=n, theta=theta)

# ------------------------------------------------------------------ S3 torsion and quadratic residues mod q
def qr_table(q):
    t = bytearray(q)
    for x in range(1, q): t[x * x % q] = 1
    return t
def count_mod(E, q, tab=None):
    """#E(F_q) for an integral model, q odd prime of good reduction."""
    tab = tab or qr_table(q)
    a2, a4, a6 = E.a2 % q, E.a4 % q, E.a6 % q
    N = q + 1
    for x in range(q):
        v = (x * x * x + a2 * x * x + a4 * x + a6) % q
        if v: N += 1 if tab[v] else -1
    return N

def torsion(Em):
    roots = [e for e in RI.shift_to_two_torsion(Em.a2, Em.a4, Em.a6)[:1]]
    e = roots[0]; T = (Fr(e), Fr(0)); assert Em.on(T)
    g = 0; used = []
    for q in range(3, 400):
        if RI.is_prime_td(q) and Em.D % q:
            g = gcd(g, count_mod(Em, q)); used.append(q)
            if g == 2 and len(used) >= 5: break
    assert g == 2, ('torsion', g)
    return T, used

# ------------------------------------------------------------------ S5 heights (sigma_exact's bounds, PARI normalisation)
def height_bounds(E):
    """B, U with h(x(P)) - hhat(P) <= B and hhat(P) - h(x(P)) <= U for every rational point (Silverman 1990, Thm 1,
    doubled for the normalisation of PARI), for an integral model; as in sigma_exact.py."""
    j = Fr(E.c4 ** 3, E.D); jsize = max(abs(j.numerator), j.denominator)
    common = SE.ln_bounds(abs(E.D))[1] / 6 + SE.lnplus_hi(abs(j)) / 6 + SE.lnplus_hi(Fr(abs(E.b2), 12)) + (SE.LN2[1] if E.b2 else 0)
    return common + SE.ln_bounds(jsize)[1] / 12 + Fr(1946, 1000), common + Fr(214, 100)

def enclose(E, P, k=7):
    B, U = height_bounds(E)
    x = P[0]
    for _ in range(k): x = E.dup_x(x)
    hx = SE.hnaive(x)
    return (hx[0] - B) / 4 ** k, (hx[1] + U) / 4 ** k, B, U

# ------------------------------------------------------------------ S6 saturation
def saturation(Em, G, T, lmax, qmax=60000):
    """For every prime l <= lmax a homomorphism E(Q) -> Z/l that kills l E(Q) and the torsion and not G.
    l odd: a good prime q with l || #E(F_q) and (N/l) G != O mod q; P -> (N/l) P mod q is such a map (it kills the
    torsion, of order 2).  l = 2: a good prime q at which the 2-part of E(F_q) is cyclic (the quadratic factor of the
    cubic irreducible mod q), and f(P) = [(N/2) P != O mod q], the map E(Q) -> E(F_q)/2E(F_q) = Z/2.  If f(T) = 0 and
    f(G) = 1, G is not in 2E(Q) + {O, T}; if f(T) = 1, two such q with different values on G are needed (then G and
    G + T are both outside 2 E(Q))."""
    den = G[0].denominator * G[1].denominator
    e = int(T[0])
    a, b = 3 * e + Em.a2, 3 * e * e + 2 * Em.a2 * e + Em.a4          # y^2 = (x - e)((x - e)^2 + a (x - e) + b)
    bp = a * a - 4 * b
    ls = [l for l in range(2, lmax + 1) if RI.is_prime_td(l)]
    need = {l: [] for l in ls}
    done = set()
    for q in range(5, qmax):
        if len(done) == len(ls): break
        if not RI.is_prime_td(q) or Em.D % q == 0 or den % q == 0: continue
        N = count_mod(Em, q)
        for l in ls:
            if l in done or N % l: continue
            if l > 2:
                if (N // l) % l == 0: continue
                assert (N // l) % 2 == 0, ('N/l odd', q, l)   # T mod q has order 2, so N is even: (N/l) T = O
                v = mulmod(Em, N // l, red(G, q), q)
                if v is not None: need[l].append(q); done.add(l)
            else:
                if q > 3000: continue
                # brute force over E(F_q): the set 2 E(F_q) of doubles; the reduction is a homomorphism, so if neither G
                # nor G + T reduces into it, neither is in 2 E(Q)
                tab = qr_table(q); sqr = {}
                for y in range(q): sqr.setdefault(y * y % q, []).append(y)
                doubles = set([None])
                for x in range(q):
                    v = (x ** 3 + Em.a2 * x * x + Em.a4 * x + Em.a6) % q
                    for y in sqr.get(v, []): doubles.add(addmod(Em, (x, y), (x, y), q))
                GT = Em.add(G, T)
                if red(G, q) not in doubles and red(GT, q) not in doubles:
                    need[l].append((q, 'G and G + T not in 2E(F_q)')); done.add(l)
    missing = [l for l in ls if l not in done]
    assert not missing, ('saturation not shown at', missing)
    return need

def red(P, q):
    if P is None: return None
    x, y = P
    if x.denominator % q == 0: return None
    return (x.numerator * pow(x.denominator, -1, q) % q, y.numerator * pow(y.denominator, -1, q) % q)
def addmod(E, P, Q, q):
    if P is None: return Q
    if Q is None: return P
    if P[0] == Q[0]:
        if (P[1] + Q[1]) % q == 0: return None
        l = (3 * P[0] * P[0] + 2 * E.a2 * P[0] + E.a4) * pow(2 * P[1], -1, q) % q
    else:
        l = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, q) % q
    x = (l * l - E.a2 - P[0] - Q[0]) % q
    return (x, (l * (P[0] - x) - P[1]) % q)
def mulmod(E, k, P, q):
    R = None
    while k:
        if k & 1: R = addmod(E, R, P, q)
        P = addmod(E, P, P, q); k >>= 1
    return R

# ------------------------------------------------------------------ S6b a lower bound for the heights, without Sage
def minimal_at(E, p):
    """The integral model E (a1 = a3 = 0) is minimal at p: either v(c4) < 4, v(c6) < 6 or v(D) < 12, or (written for
    p = 2 only) no change of coordinates with u = 2 and integers r, s, t, which may be taken modulo 2^6, gives integral
    coefficients.  If an integral model with u = 2 existed, r, s, t would be 2-integral: r is a root of
    3 r^4 + b2 r^3 + 3 b4 r^2 + 3 b6 r + b8 - u^8 b8' = 0, whose leading coefficient 3 is a unit at 2, and then
    2 s = u a1' and 2 t = u^3 a3' - r a1 (here a1 = a3 = 0) are integral, with s = a1', t = 4 a3'."""
    v = lambda x: RI.vp(x, p) if x else 99
    if v(E.c4) < 4 or v(E.c6) < 6 or v(E.D) < 12: return True
    assert p == 2, 'minimality test by search written for p = 2 only'
    a1, a2, a3, a4, a6 = 0, E.a2, 0, E.a4, E.a6
    u = 2
    for r in range(64):
        for s_ in range(64):
            if (a1 + 2 * s_) % u: continue
            if (a2 - s_ * a1 + 3 * r - s_ * s_) % u ** 2: continue
            for t in range(64):
                if (a3 + r * a1 + 2 * t) % u ** 3: continue
                if (a4 - s_ * a3 + 2 * r * a2 - (t + r * s_) * a1 + 3 * r * r - 2 * s_ * t) % u ** 4: continue
                if (a6 + r * a4 + r * r * a2 + r ** 3 - t * a3 - t * t - r * t * a1) % u ** 6: continue
                return False
    return True

def lam_exact(Em, need):
    """A proof that every point P of infinite order of E1(Q) has hhat(P) >= need (PARI normalisation), without Sage.
    (1) The model is minimal at every bad prime (minimal_at).  (2) At a bad p with v(c4) = 0 the reduction is
    multiplicative of type I_n, n = v(D), and E(Q_p)/E0(Q_p) is cyclic of order n (split) or of order <= 2 dividing n
    (non-split); with v(c4) > 0 it is additive and E(Q_p)/E0(Q_p) has order <= 4, so exponent 1, 2, 3 or 4.  So
    Q = m P with m = lcm(n_p, 12) lies in E0(Q_p) for every p.  (3) For such Q, on a minimal model, the local term at
    every prime vanishes and h(x(Q)) - hhat(Q) <= log+|j|/4 + log+|b2/12| + log 2* + 1.946 (Cremona-Prickett-Siksek,
    J. Number Theory 116 (2006), Theorem 2 for the vanishing of the non-archimedean terms on E0, and Theorem 15, the
    upper side of Silverman's archimedean bound, Silverman 1990 p. 737, in their normalisation, where hhat is PARI's and
    h(x) = log max(|a|, e^2); the same statement as in height_enclosure.py, with every alpha_p equal to 0).  (4) If
    hhat(P) < need, then h(x(Q)) < m^2 need + that constant =: hb, and x(Q) = a/e^2 with |a|, e^2 <= exp(hb): the search
    below finds no point there, so no such P exists.  The search is refused when exp(hb) > XMAX (a forged generator
    k G with k > 199 has a large height and would otherwise make the search astronomically long).  That G then
    generates needs rank 1 (S4 here; for other curves the caller must supply it)."""
    bad = sorted(RI.factor_int(Em.D))
    for p in bad: assert minimal_at(Em, p), ('not minimal at', p)
    m = 12
    for p in bad:
        if Em.c4 % p: m = m * RI.vp(Em.D, p) // gcd(m, RI.vp(Em.D, p))
    j = Fr(Em.c4 ** 3, Em.D)
    Binf = SE.lnplus_hi(abs(j)) / 4 + SE.lnplus_hi(Fr(abs(Em.b2), 12)) + (SE.LN2[1] if Em.b2 else 0) + Fr(1946, 1000)
    hb = m * m * need + Binf
    assert SE.ln_bounds(XMAX)[0] > hb, ('the height search would go beyond %d: refused (need too large)' % XMAX)
    X0 = 1
    while SE.ln_bounds(X0 + 1)[0] <= hb: X0 += 1          # every integer of absolute value <= exp(hb) is <= X0
    assert X0 < XMAX
    found = []
    e = 1
    while e * e <= X0:
        e2 = e * e; e4 = e2 * e2; e6 = e4 * e2
        for a in range(-X0, X0 + 1):
            # f(a/e^2) is a square iff e^6 f(a/e^2) = a^3 + a2 a^2 e^2 + a4 a e^4 + a6 e^6 is the square of an integer
            v = a * a * a + Em.a2 * a * a * e2 + Em.a4 * a * e4 + Em.a6 * e6
            if v >= 0 and isqrt(v) ** 2 == v and gcd(a, e) == 1: found.append(Fr(a, e2))
        e += 1
    assert not found, ('points of small height found', found)
    return dict(m=m, Binf=Binf, hb=hb, X0=X0)

def sym_gq(z):
    """An exact sympy Gaussian rational (no floats) as a GQ; refused if it is not one."""
    import sympy
    z = sympy.expand(sympy.cancel(z))
    re_, im_ = sympy.re(z), sympy.im(z)
    assert re_.is_Rational and im_.is_Rational, ('not an exact Gaussian rational', z)
    return GQ(Fr(int(re_.p), int(re_.q)), Fr(int(im_.p), int(im_.q)))

# ------------------------------------------------------------------ S7 Lemma 2 over Q(i), with sympy
def lemma2(St):
    import sympy
    t = sympy.Symbol('t'); Iy = sympy.I
    def fpoly(Q, sg):
        s = sg * (-Iy) * t
        return sympy.Poly(sympy.expand(sum(St['d'] * c * s ** k for k, c in enumerate(Q))), t, domain='QQ_I')
    f = [fpoly(St['Q1'], 1), fpoly(St['Q1'], -1), fpoly(St['Q2'], 1), fpoly(St['Q2'], -1)]
    for g in f:
        assert g.degree() == 4 and sympy.gcd(g, g.diff(t)).degree() == 0, 'f_i not squarefree'
    deg = {}
    for i, k in itertools.combinations(range(4), 2):
        deg[(i + 1, k + 1)] = sympy.gcd(f[i], f[k]).degree()
    # the cycle 1-3-2-4-1: each f_i shares two roots with each of its two neighbours and none with the opposite one
    assert deg == {(1, 2): 0, (1, 3): 2, (1, 4): 2, (2, 3): 2, (2, 4): 2, (3, 4): 0}, deg
    prod = f[0] * f[1] * f[2] * f[3]
    rad = sympy.quo(prod, sympy.gcd(prod, prod.diff(t)))
    assert rad.degree() == 8, 'distinct roots'
    # f1 f2 f3 f4 = c g^2 with c a square in K: the product is lc * (rad/lc(rad))^2 ... check exactly
    g = sympy.Poly(rad.monic().as_expr(), t, domain='QQ_I')
    c = sympy.cancel(prod.as_expr() / (g.as_expr() ** 2))
    cg = sym_gq(c)
    assert gsqrt(cg) is not None, 'the constant is not a square in K'
    # every root of the product is a root of exactly two f_i: the eight roots are those of the four gcds of degree 2
    for i in range(4):
        nb = [k for k in range(4) if k != i and deg[tuple(sorted((i + 1, k + 1)))] == 2]
        q = sympy.gcd(f[i], f[nb[0]]) * sympy.gcd(f[i], f[nb[1]])
        assert sympy.rem(f[i], q).is_zero and q.degree() == 4
    return deg

# ------------------------------------------------------------------ S8 Lemmas 4 and 5
def hq_hi(q):
    q = Fr(q)
    return SE.hnaive(q)[1]
def log_hi(n): return SE.ln_bounds(n)[1] if n > 1 else Fr(0)

def curve_consts(Ew, A0, r, short):
    """For the weierstrass model Ew with x = A0/(s - r), and its integral short model `short` = (A, B): the bound K for
    hhat (CVV normalisation) of a point above s a root of unity, and the constants c1, c2 of CVV Table 1."""
    A, B = short
    Es = Ell(0, A, B)
    ip = iso_params(Ew, Es); assert ip is not None, 'short model not isomorphic over Q'
    u, rr = ip
    alpha, beta = 1 / (u * u), -rr / (u * u)                       # x_short = alpha x + beta
    HX = hq_hi(alpha * A0) + hq_hi(r) + hq_hi(beta) + 2 * SE.LN2[1]
    Dl = -16 * (4 * A ** 3 + 27 * B * B); j = Fr(-1728 * (4 * A) ** 3, Dl)
    hD = log_hi(abs(Dl)); hinf = SE.lnplus_hi(abs(j)); hj = hq_hi(j)
    hA, hB = log_hi(abs(A)), log_hi(abs(B))
    c1g = (hD + hinf) / 4 + hj / 8 + (hA + hB) / 2 + CVV['c1g']
    c2g = (hD + hinf) / 4 + (hA + hB) / 2 + CVV['c2g']
    c1q = (hD + hinf) / 4 + hj / 8 + log_hi(abs(A) + abs(B) + 3) / 2 + CVV['c1q']
    c2q = (hD + hinf) / 4 + log_hi(abs(A) + abs(B) + 3) / 2 + CVV['c2q']
    L3 = SE.ln_bounds(3)[1]
    Kpt = Fr(3, 2) * HX + log_hi(1 + abs(A) + abs(B)) / 2 + L3 / 2 + c2g
    return dict(HX=HX, c1g=c1g, c2g=c2g, c1q=c1q, c2q=c2q, Kpt=Kpt, alpha=alpha, beta=beta)

def lemma45(St, hlo, degC=DEG_C):
    Em = St['Em']
    A, B = -27 * Em.c4, -54 * Em.c6
    k1 = curve_consts(St['E1w'], St['A1'], St['r1'], (A, B))
    k2 = curve_consts(St['E2w'], St['A2'], St['r2'], (A, -B))     # E2 is the twist by -1: short model (A, -B)
    mu = 4 * k1['Kpt'] + 4 * k2['Kpt'] + 2 * k1['c1g']             # both coordinates of (alpha(w), beta(w)) are on E1
    h2C = 2 * degC * mu
    C1 = CVV['C1']
    bounds = []
    for c1, c2 in ((k1['c1q'], k1['c2q']), (k1['c1g'], k1['c2g'])):
        C2 = C1 * (CVV['k6'] + 4 * c1); C3 = 4 * c2
        bounds.append(C1 * h2C * degC + C2 * degC * degC + C3)
    top = max(bounds)
    ab = isqrt(int(top / (Fr(3, 2) * hlo)))
    while (ab + 1) ** 2 * Fr(3, 2) * hlo <= top: ab += 1
    while ab * ab * Fr(3, 2) * hlo > top: ab -= 1
    return dict(k1=k1, k2=k2, mu=mu, h2C=h2C, top=top, ab=ab)

# ------------------------------------------------------------------ S9 the sieve, exact
class EllP:
    def __init__(s, E, p):
        s.p = p
        s.a2, s.a4, s.a6 = (Fr(c).numerator * pow(Fr(c).denominator, -1, p) % p for c in (E.a2, E.a4, E.a6))
    def add(s, P, Q):
        p = s.p
        if P is None: return Q
        if Q is None: return P
        if P[0] == Q[0]:
            if (P[1] + Q[1]) % p == 0: return None
            l = (3 * P[0] * P[0] + 2 * s.a2 * P[0] + s.a4) * pow(2 * P[1], -1, p) % p
        else:
            l = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, p) % p
        x = (l * l - s.a2 - P[0] - Q[0]) % p
        return (x, (l * (P[0] - x) - P[1]) % p)
    def neg(s, P): return None if P is None else (P[0], (-P[1]) % s.p)
    def on(s, P): return P is None or (P[1] * P[1] - (P[0] ** 3 + s.a2 * P[0] ** 2 + s.a4 * P[0] + s.a6)) % s.p == 0

def rmod(q, p):
    q = Fr(q); assert q.denominator % p, 'not p-integral'
    return q.numerator * pow(q.denominator, -1, p) % p
def redP(P, p):
    if P is None: return None
    if Fr(P[0]).denominator % p == 0: return None
    return (rmod(P[0], p), rmod(P[1], p))
def nd(q): q = Fr(q); return q.numerator * q.denominator

def admissible(St, p):
    """The conditions of Lemma 6 on the prime p, all exact."""
    if not RI.is_prime_td(p): return 'not prime'
    if p % 4 != 1: return 'not 1 mod 4'
    E1w, E2w = St['E1w'], St['E2w']
    for E in (E1w, E2w):
        for c in (E.a2, E.a4, E.a6):
            if Fr(c).denominator % p == 0: return 'model not p-integral'
        if Fr(E.D).numerator % p == 0: return 'bad reduction'
    for q in (St['A1'], St['A2'], St['r1'], St['r2'], St['d'], St['u2'], St['rr2']):
        if nd(q) % p == 0: return 'divides a constant of the maps'
    for Q in (St['dQ1'], St['dQ2']):
        if Q[-1] % p == 0: return 'leading coefficient'
        import sympy
        s = sympy.Symbol('s')
        if int(sympy.discriminant(sympy.Poly(list(reversed(Q)), s))) % p == 0: return 'quartic not separable mod p'
    j = St['j']
    if j.denominator % p == 0: return 'j not p-integral'
    jp = rmod(j, p)
    if jp in (0, 1728 % p): return 'j = 0 or 1728 mod p'
    return None

def sqrt_table(p):
    t = {}
    for x in range(p): t.setdefault(x * x % p, x)
    return t

def sieve_local(St, p, Gw, Tw):
    why = admissible(St, p)
    assert why is None, ('prime %d not admissible: %s' % (p, why))
    e1, e2 = EllP(St['E1w'], p), EllP(St['E2w'], p)
    sq = sqrt_table(p)
    iota = sq[p - 1]
    assert iota * iota % p == p - 1
    A1, A2, r1, r2 = (rmod(v, p) for v in (St['A1'], St['A2'], St['r1'], St['r2']))
    q1 = [c % p for c in St['dQ1']]; q2 = [c % p for c in St['dQ2']]
    u2, rr2 = rmod(St['u2'], p), rmod(St['rr2'], p)
    iu2 = pow(u2, -1, p)
    def tau(R):                     # E2 -> twist (over Q, reduced) -> E1: (x, y) -> (-X, iota Y)
        if R is None: return None
        X = (R[0] - rr2) * iu2 * iu2 % p; Y = R[1] * iu2 * iu2 * iu2 % p
        return ((-X) % p, iota * Y % p)
    Gp, Tp = redP(Gw, p), redP(Tw, p)
    assert Gp is not None and Tp is not None, ('point reduces to O', p)
    assert e1.on(Gp) and e1.on(Tp)
    dlog = {}; X = None; a = 0
    while True:
        dlog[X] = a; X = e1.add(X, Gp); a += 1
        if X is None: break
    o = a
    def pts(e, A, q, r, s):
        if s is None:
            v = q[4]
            return [(0, A * y % p) for y in {sq[v], (-sq[v]) % p}] if v in sq else []
        if s == r: return [None]
        v = sum(c * pow(s, k, p) for k, c in enumerate(q)) % p
        if v not in sq: return []
        w = pow((s - r) % p, -1, p); y = sq[v]
        return [(A * w % p, A * yy * w * w % p) for yy in {y, (-y) % p}]
    def res(Z):
        out = set()
        for eps in (0, 1):
            Y = Z if eps == 0 else e1.add(Z, e1.neg(Tp))
            if Y in dlog: out.add(dlog[Y]); out.add((-dlog[Y]) % o)
        return out
    Sp = set()
    for t in list(range(p)) + [None]:
        s1, s2 = (None, None) if t is None else ((-iota * t) % p, (iota * t) % p)
        P1, P2 = pts(e1, A1, q1, r1, s1), pts(e1, A1, q1, r1, s2)
        if not (P1 and P2): continue
        R1, R2 = pts(e2, A2, q2, r2, s1), pts(e2, A2, q2, r2, s2)
        if not (R1 and R2): continue
        for X in P1 + P2 + R1 + R2: assert (X is None) or (e1.on(X) if X in P1 + P2 else e2.on(X))
        ra = set(); rb = set()
        for X in P1:
            for Y in P2: ra |= res(e1.add(X, Y))
        for X in R1:
            for Y in R2:
                Z = tau(e2.add(X, e2.neg(Y))); assert e1.on(Z); rb |= res(Z)
        for x in ra:
            for y in rb: Sp.add((x, y))
    return o, Sp

def lift(loc, primes, N):
    (o1, S1), (o2, S2) = loc[primes[0]], loc[primes[1]]
    if gcd(o1, o2) != 1:
        # the orders are not coprime: lift from the first prime only, in a range of at most 30000
        N = min(N, 30000)
        cand = set((a, b) for (a0, b0) in S1 for a in range(a0, N + 1, o1) for b in range(b0, N + 1, o1))
        n12 = len(cand)
        for p in primes[1:]:
            o, Sp = loc[p]; cand = set(c for c in cand if (c[0] % o, c[1] % o) in Sp)
        return sorted(cand), N, n12
    Mo = o1 * o2; inv = pow(o1, -1, o2); N = min(N, Mo - 1)
    crt = lambda x1, x2: x1 + o1 * (((x2 - x1) * inv) % o2)
    A1s, B2 = {}, {}
    for (a1, b1) in S1: A1s.setdefault(a1, []).append(b1)
    for (a2, b2) in S2: B2.setdefault(a2, []).append(b2)
    cand = set()
    for a1, b1s in A1s.items():
        for a2, b2s in B2.items():
            a = crt(a1, a2)
            if a > N: continue
            for b1 in b1s:
                for b2 in b2s:
                    b = crt(b1, b2)
                    if b <= N: cand.add((a, b))
    n12 = len(cand)
    for p in primes[2:]:
        o, Sp = loc[p]; cand = set(c for c in cand if (c[0] % o, c[1] % o) in Sp)
    return sorted(cand), N, n12

# ------------------------------------------------------------------ S10 the fibre a = 1, exact
def fibre(St, consts):
    import sympy
    x = sympy.Symbol('x')
    E = St['E1w']; A0, r = St['A1'], St['r1']
    sy = lambda q: sympy.Rational(Fr(q).numerator, Fr(q).denominator)
    f = x ** 3 + sy(E.a2) * x ** 2 + sy(E.a4) * x + sy(E.a6)
    xs = set([Fr(-A0) / (2 * r), Fr(0)])
    for c in consts:
        x0, y0 = sy(c[0]), sy(c[1]); xs.add(Fr(c[0]))
        Dn = 2 * sy(r) * x + sy(A0); hD = (-sy(A0) * x + (sy(E.a2) + x + x0) * Dn) * (x - x0) ** 2
        pol = sympy.expand((f * Dn - hD) if y0 == 0 else (((f + y0 ** 2) * Dn - hD) ** 2 - 4 * y0 ** 2 * f * Dn ** 2))
        P = sympy.Poly(pol, x, domain='QQ'); assert not P.is_zero, 'identically satisfied'
        _, facs = P.factor_list()
        prodchk = sympy.Poly(1, x, domain='QQ')
        for g, k in facs: prodchk = prodchk * g ** k
        assert sympy.Poly(P.monic(), x) == sympy.Poly(prodchk.monic(), x), 'factorisation'
        for g, _ in facs:
            co = [Fr(str(cc)) for cc in reversed(g.all_coeffs())]
            if len(co) == 2: xs.add(GQ(-co[0] / co[1]))
            elif len(co) == 3:
                disc = co[1] ** 2 - 4 * co[2] * co[0]
                if disc < 0 and issq_q(-disc):
                    w = sqrt_q(-disc)
                    for sg in (1, -1): xs.add(GQ(-co[1] / (2 * co[2]), sg * w / (2 * co[2])))
    ts = set()
    for xv in xs:
        xv = gq(xv)
        if xv.is_zero(): continue                                       # s = infinity, t = infinity: treated in S11
        s = gq(r) + gq(A0) / xv
        t = I * s
        if t.im == 0: ts.add(t.re)
    return sorted(ts), len(xs)

def fibre_K(E, A0, r, consts):
    """The same equation x(P - c) = g(x(P)) for points c of E over K = Q(i), factored over Q(i) (sympy, domain QQ_I):
    the roots in K are those of the linear factors.  Returns the rational t = i s found and the number of x-values."""
    import sympy
    x = sympy.Symbol('x'); Iy = sympy.I
    sy = lambda q: sympy.Rational(Fr(q).numerator, Fr(q).denominator)
    syk = lambda z: sy(gq(z).re) + Iy * sy(gq(z).im)
    f = x ** 3 + sy(E.a2) * x ** 2 + sy(E.a4) * x + sy(E.a6)
    xs = set([gq(Fr(-A0) / (2 * r)), GQ(0)])
    for c in consts:
        x0, y0 = syk(c[0]), syk(c[1]); xs.add(gq(c[0]))
        Dn = 2 * sy(r) * x + sy(A0); hD = (-sy(A0) * x + (sy(E.a2) + x + x0) * Dn) * (x - x0) ** 2
        pol = sympy.expand((f * Dn - hD) if y0 == 0 else (((f + y0 ** 2) * Dn - hD) ** 2 - 4 * y0 ** 2 * f * Dn ** 2))
        P = sympy.Poly(pol, x, domain='QQ_I'); assert not P.is_zero, 'identically satisfied'
        _, facs = P.factor_list()
        prodchk = sympy.Poly(1, x, domain='QQ_I')
        for g, k in facs: prodchk = prodchk * g ** k
        assert P.monic() == prodchk.monic(), 'factorisation'
        for g, _ in facs:
            if g.degree() == 1:
                a1, a0 = g.all_coeffs()
                xs.add(sym_gq(-a0 / a1))
    ts = set()
    for xv in xs:
        if xv.is_zero(): continue
        sv = gq(r) + gq(A0) / xv; t = I * sv
        if t.im == 0: ts.add(t.re)
    return sorted(ts), len(xs)

def is_solution(St, t):
    """d Q1(-i t) and d Q2(-i t) are both squares in K."""
    s = GQ(0, -Fr(t))
    return all(gsqrt(peval([gq(c) for c in Q], s)) is not None for Q in (St['dQ1'], St['dQ2']))

# ------------------------------------------------------------------ the whole check
def load_hints():
    if os.path.exists(HINTS_FILE):
        H = json.load(open(HINTS_FILE))
        for k, v in H.items(): RI.HINTS[k] = v
def save_hints():
    if WRITE: json.dump(RI.HINTS, open(HINTS_FILE, 'w'), indent=0, sort_keys=True)

def rank_part(Em, G, say=say):
    e, a, b = RI.shift_to_two_torsion(Em.a2, Em.a4, Em.a6)
    GE = (G[0] - e, G[1])
    ap, bp = -2 * a, a * a - 4 * b
    r, rec = RI.rank_bound(a, b, [GE], RI.find_points([0, ap, 0, bp, 0]), say=say)
    return r, rec

def main(lmax=199, sieve_N=596309, primes=SIEVE_PRIMES, claimed_ab=CLAIMED_AB, lam=LAM_SAGE, hG_claim_lo=None):
    load_hints()
    say('exact_last_pair.py: Theorem A of LAST_PAIR.md for (627,364), exact arithmetic (no Sage, PARI or Magma)')
    St = setup()
    Em, G = St['Em'], St['G']
    say('S1  models: psi_1, psi_2 are maps C_j -> E_j (identity of rational functions); E1 ~ stated minimal model over Q '
        '(u = %s); E2 ~ twist of E1 by -1 over Q (u = %s); theta = (x, y) -> (-x, i y) after that, sigma theta = -theta; '
        'j(E1) not an integer (no CM)' % (St['u1'], St['u2']))
    (A1p, B1p), (Am1, Bm1) = St['pts'][1], St['pts'][-1]
    say('S2  t = 1: A = P + conj(P) = %s on the minimal model, the stated G up to sign: True' % ('G' if A1p[1] == G[1] else '-G'))
    T, tq = torsion(Em)
    E = Em
    def same_class(Z):            # Z = +-G + eps T ?
        for sg in (1, -1):
            for eps in (0, 1):
                W = E.mul(sg, G)
                if eps: W = E.add(W, T)
                if W == Z: return '%sG%s' % ('' if sg == 1 else '-', ' + T' if eps else '')
        return None
    labs = [same_class(Z) for Z in (A1p, B1p, Am1, Bm1)]
    assert all(labs), labs
    say('    the points (A, B) at t = 1 and t = -1 are %s, %s and %s, %s: (a, b) = (1, 1) at both' % tuple(labs))
    say('S3  torsion: T = %s; gcd of #E1(F_q) over the good odd primes %s is 2, so E1(Q)_tors = {O, T}' % ((int(T[0]), 0), tq))
    r, rec = rank_part(Em, G)
    assert r == 1, ('rank bound', r)
    say('S4  rank E1(Q) <= %d by the exact descent via 2-isogeny and the second descent (rank_indep.py); G has infinite order (S5): rank 1' % r)
    lo, hi, Bh, Uh = enclose(Em, G, 7)
    assert lo > 0
    if hG_claim_lo is not None: assert hG_claim_lo <= hi, 'claimed lower bound for hhat(G) above the proved upper bound'
    say('S5  hhat(G) (PARI normalisation) in [%.6f, %.6f], from 2^7 G and h - hhat <= %.3f, hhat - h <= %.3f (Silverman 1990, Thm 1)'
        % (float(lo) - 5e-7, float(hi) + 5e-7, float(Bh), float(Uh)))
    assert lam <= hi, 'a lower bound for the heights of points of infinite order above hhat(G) is impossible'
    nxt0 = lmax + 1
    while not RI.is_prime_td(nxt0): nxt0 += 1
    need = hi / (nxt0 - 1) ** 2
    LX = lam_exact(Em, need)
    msg6b = ('S6b lower bound for heights without Sage: model minimal at every bad prime; m = %d puts every point in E0 at every '
        'prime; a point P of infinite order with hhat(P) < hhat(G)/%d^2 would give Q = mP with h(x(Q)) < %d^2 * %.2e + %.3f = %.3f, '
        'x(Q) = a/e^2 with |a|, e^2 <= %d: none of these is a point.  So hhat(P) >= %.2e for every P of infinite order.  '
        'If G = k G0 + torsion with k > 1, then k has no prime factor <= %d (S6), so k >= %d and hhat(G0) <= hhat(G)/%d^2 < %.2e: '
        'impossible.  So G generates E1(Q) modulo torsion, with no input from Sage.'
        % (LX['m'], nxt0 - 1, LX['m'], float(need), float(LX['Binf']), float(LX['hb']), LX['X0'], float(need), lmax, nxt0, nxt0, float(need)))
    lam_used = lam / 2
    K = 1
    while (K + 1) ** 2 * lam_used <= hi: K += 1
    assert lmax >= K, ('saturation range below the index bound', lmax, K)
    sat = saturation(Em, G, T, lmax)
    lmin_ok = 2
    while True:
        nxt = lmin_ok + 1
        while not RI.is_prime_td(nxt): nxt += 1
        if nxt > lmax: break
        lmin_ok = nxt
    # G generates as soon as lam_true >= hi / (next prime above lmax)^2
    nxt = lmax + 1
    while not RI.is_prime_td(nxt): nxt += 1
    say('S6  saturation: for every prime l <= %d a functional E1(Q) -> Z/l separates G from l E1(Q) + torsion (primes q used: %s ...).'
        % (lmax, ', '.join('%d:%s' % (l, sat[l][0] if l > 2 else sat[l][:2]) for l in list(sat)[:6])))
    say(msg6b)
    tst = []
    for k in (2, 4, 8, 12, 24):
        Q = Em.mul(k, G); hx = SE.hnaive(Q[0])
        up_ = hx[1] - k * k * lo
        assert up_ <= LX['Binf'], ('a multiple of G violates the bound used in S6b', k)
        tst.append('%dG: %.3f' % (k, float(up_)))
    say('    test on real points: h(x(kG)) - hhat(kG) <= %s, all <= %.3f as S6b says for points in E0 everywhere '
        '(for G itself it is 18.96: G is not in E0 at every prime)' % (', '.join(tst), float(LX['Binf'])))
    say('    (comparison only, not used: with Sage\'s lam = %s, halved, the index would be at most %d)' % (lam, K))
    deg = lemma2(St)
    say('S7  Lemma 2: f1..f4 squarefree of degree 4; degrees of the gcds %s (the cycle 1-3-2-4-1); 8 distinct roots, each a root '
        'of exactly two f_i; f1 f2 f3 f4 = c g^2 with c a square in K.  So only the empty and the full product are a '
        'constant times a square, and a component W has genus 1 + (8(-2) + 8*4)/2 = 9' % deg)
    L = lemma45(St, lo)
    say('S8  Lemma 4: HX(E1) <= %.3f, HX(E2) <= %.3f; CVV constants c1 (general) <= %.3f, c2 (general) <= %.3f, c1 (Q) <= %.3f, '
        'c2 (Q) <= %.3f; hhat bound above a root of unity <= %.2f (E1), %.2f (E2); mu <= %.2f; h2(C) <= %.4e'
        % (float(L['k1']['HX']), float(L['k2']['HX']), float(L['k1']['c1g']), float(L['k1']['c2g']), float(L['k1']['c1q']),
           float(L['k1']['c2q']), float(L['k1']['Kpt']), float(L['k2']['Kpt']), float(L['mu']), float(L['h2C'])))
    say('    Lemma 5 (CVV Thm 1.2, N = 2, deg C <= %d): hhat(A, B) <= %.4e; with hhat(A, B) = 1.5 hhat(G)(a^2 + b^2) >= 1.5 * %.6f (a^2 + b^2): a, b <= %d'
        % (DEG_C, float(L['top']), float(lo), L['ab']))
    assert L['ab'] <= claimed_ab, ('the bound of LAST_PAIR.md is below the one computed here', claimed_ab, L['ab'])
    L192 = lemma45(St, lo, 192)
    say('    with deg C <= 192: a, b <= %d' % L192['ab'])
    Gw = from_model(G, St['u1'], St['rr1']); Tw = from_model(T, St['u1'], St['rr1'])
    assert St['E1w'].on(Gw) and St['E1w'].on(Tw)
    loc = {}
    for p in primes:
        loc[p] = sieve_local(St, p, Gw, Tw)
        say('S9  p = %d (admissible): order of G mod p %d, pairs (a mod o, b mod o) allowed %d of %d' % (p, loc[p][0], len(loc[p][1]), loc[p][0] ** 2))
    left, Ncov, n12 = lift(loc, list(primes), sieve_N)
    say('    lifting: pairs 0 <= a, b <= %d allowed by the first two primes %d; after all %d primes: %s' % (Ncov, n12, len(primes), left))
    assert left == [(1, 1)], left
    assert Ncov >= L192['ab'] >= L['ab'], ('the sieve does not cover the bound', Ncov, L['ab'])
    plant = sorted(set([(2, 1), (Ncov, 3), (Ncov // 2 + 7, Ncov - 5), (31337 % Ncov, 4242421 % Ncov), (loc[primes[0]][0] + 1, Ncov - 2)]))
    loc2 = {p: (o, Sp | set((a % o, b % o) for (a, b) in plant)) for p, (o, Sp) in loc.items()}
    left2, _, _ = lift(loc2, list(primes), sieve_N)
    assert all(x in left2 for x in plant) and set(left2) - set(plant) - set(left) == set(), 'lifting test'
    say('    test of the lifting: %d planted pairs all returned, nothing else' % len(plant))
    consts = [Gw, St['E1w'].add(Gw, Tw)]
    ts, nx = fibre(St, consts)
    sols = [t for t in ts if is_solution(St, t)]
    say('S10 fibre a = 1 (A = +-G or +-G + T): %d values of x(P) in Q(i) (roots and special points); rational t among them: %s; '
        'of these, solutions of Theorem A: %s' % (nx, [str(t) for t in ts], [str(t) for t in sols]))
    assert sorted(sols) == [Fr(-1), Fr(1)]
    # the same for B (not needed: Lemma 7 uses a = 1 only), on E2 over K with c' = theta^-1(c)
    def theta_inv(P):
        Xt, Yt = -gq(P[0]), GQ(0, -1) * gq(P[1])
        return (Xt * St['u2'] ** 2 + St['rr2'], Yt * St['u2'] ** 3)
    cps = [theta_inv(c) for c in consts]
    for cp in cps: assert cp[1] * cp[1] == St['E2w'].f(cp[0])
    ts2, nx2 = fibre_K(St['E2w'], St['A2'], St['r2'], cps)
    sols2 = [t for t in ts2 if is_solution(St, t)]
    say('    the same for b = 1 (B = +-G or +-G + T, on E2 over K, factored over Q(i)): rational t %s, solutions %s'
        % ([str(t) for t in ts2], [str(t) for t in sols2]))
    assert sorted(sols2) == [Fr(-1), Fr(1)]
    z0 = [gsqrt(GQ(Q[0])) for Q in (St['dQ1'], St['dQ2'])]
    zi = [gsqrt(GQ(Q[-1])) for Q in (St['dQ1'], St['dQ2'])]      # leading coefficient in t: d lc(Q) (-i)^4 = d lc(Q)
    assert z0[0] is None and z0[1] is None and zi[0] is None and zi[1] is None
    say('S11 t = 0 and t = infinity are not points (d Q_j(0) and d lc(Q_j) not squares in K).  Partner: (991,263) = (%d + %d, %d - %d).'
        % (M_, N_, M_, N_))
    assert (M_ + N_, M_ - N_) == (991, 263)
    say('THEOREM A (627,364): the rational t with d Q1(-it), d Q2(-it) both squares in Q(i) are exactly 1 and -1, given the '
        'statements from the literature listed in the header (CVV Theorem 1.2/4.2, Proposition 3.2, Table 1, Theorem 2.2; '
        'Silverman 1990 Theorem 1; Cremona-Prickett-Siksek 2006 Theorems 2 and 15; the Kodaira-Neron bounds '
        'for component groups) and the paper arguments of LAST_PAIR.md for Lemmas 2, 3 and 6.')
    save_hints()
    return dict(St=St, L=L, loc=loc, lo=lo, hi=hi)

def negative():
    load_hints()
    say('negative tests: each must be refused')
    St = setup(); Em, G = St['Em'], St['G']
    T, _ = torsion(Em)
    tests = []
    def run(name, fn):
        try:
            fn(); say('  NOT REFUSED: %s' % name); tests.append(False)
        except AssertionError as ex:
            say('  refused: %s (%s)' % (name, str(ex)[:90])); tests.append(True)
    G2 = Em.mul(2, G); G3 = Em.mul(3, G)
    run('forged generator 2G (saturation at 2)', lambda: saturation(Em, G2, T, 7))
    run('forged generator 3G (saturation at 3)', lambda: saturation(Em, G3, T, 7))
    run('forged generator 2G + T (saturation at 2)', lambda: saturation(Em, Em.add(G2, T), T, 7))
    run('forged generator 2G as the point of t = 1', lambda: setup(Gst=G2))
    run('a point off the curve as G', lambda: setup(Gst=(G[0] + 1, G[1])))
    run('another class d = 526', lambda: setup(d=526))
    Gw = from_model(G, St['u1'], St['rr1']); Tw = from_model(T, St['u1'], St['rr1'])
    run('sieve modulus 10007 (3 mod 4)', lambda: sieve_local(St, 10007, Gw, Tw))
    run('sieve modulus 10011 (composite)', lambda: sieve_local(St, 10011, Gw, Tw))
    bad = [p for p in (5, 13, 29, 37, 181) if admissible(St, p) is not None]
    run('sieve modulus %d (not admissible: %s)' % (bad[0], admissible(St, bad[0])), lambda: sieve_local(St, bad[0], Gw, Tw))
    lo, hi, _, _ = enclose(Em, G, 7)
    def lowered_bound():
        L = lemma45(St, lo)
        claimed = 13000
        assert L['ab'] <= claimed, ('computed bound %d above the claimed %d' % (L['ab'], claimed))
    run('height bound for a, b lowered to 13000 (below the computed one)', lowered_bound)
    def hhat_too_high():
        claim = Fr(2001, 100)
        assert claim <= hi, 'claimed lower bound 20.01 for hhat(G) above the proved upper bound'
    run('lower bound for hhat(G) raised to 20.01', hhat_too_high)
    def lam_too_high():
        lam = Fr(25)
        assert lam <= hi, 'lam above hhat(G)'
    run('lower bound lam for heights of points of infinite order set to 25', lam_too_high)
    def sieve_short():
        L = lemma45(St, lo)
        Ncov = 10000
        assert Ncov >= L['ab'], ('the sieve range %d does not cover the bound %d' % (Ncov, L['ab']))
    run('sieve range 10000 below the bound', sieve_short)
    def bad_sieve_set():
        # a forged local set without the residues of (1, 1): the lifting must then lose the degenerate pair, and the
        # check "only (1, 1) is left, and it is left" refuses
        loc = {}
        for p in SIEVE_PRIMES[:2]:
            o, Sp = sieve_local(St, p, Gw, Tw); loc[p] = (o, Sp - {(1, 1), (1, o - 1), (o - 1, 1), (o - 1, o - 1)})
        left, _, _ = lift(loc, list(SIEVE_PRIMES[:2]), 2000)
        assert (1, 1) in left, 'the degenerate pair (1, 1) was removed: a forged local set'
    run('a local set of the sieve with the pair of t = 1 removed', bad_sieve_set)
    def forged_211():
        # 211 G passes the saturation of S6 at every l <= 199 (211 is prime); its height is 211^2 hhat(G), so the step
        # S6b receives need = 211^2 hi / 210^2 and must refuse rather than search to exp(11000)
        lo_, hi_, _, _ = enclose(Em, G, 7)
        lam_exact(Em, 211 ** 2 * hi_ / 210 ** 2)
    run('forged generator 211 G at the height step S6b (search capped)', forged_211)
    run('the height search on a non-minimal model (E1 scaled by u = 2)', lambda: lam_exact(Ell(Em.a2 * 4, Em.a4 * 16, Em.a6 * 64), Fr(1, 10**4)))
    run('the height search on y^2 = x^3 - x + 1, which has points of small height', lambda: lam_exact(Ell(0, -1, 1), Fr(1, 100)))
    # the descent: forged hints
    e, a, b = RI.shift_to_two_torsion(Em.a2, Em.a4, Em.a6)
    ap, bp = -2 * a, a * a - 4 * b
    def forged_conic():
        dd = 33; c = bp // dd
        U0, V0, N0 = RI.find_conic_point(dd, ap, c)
        RI.second_descent(dd, ap, bp, (U0, V0, N0 + 1))
    run('a point off the conic as the hint of the second descent', forged_conic)
    def forged_factor():
        key = str(abs(bp))
        saved = RI.HINTS['factor'].get(key)
        F = dict(RI.factor_int(bp)); q1, q2 = [q for q in sorted(F) if F[q] == 2][:2]
        F2 = {k: v for k, v in F.items() if k not in (q1, q2)}; F2[q1 * q2] = 2      # same product, a composite factor
        RI.HINTS['factor'][key] = sorted(F2.items())
        try:
            RI.factor_int(bp)
        finally:
            if saved is not None: RI.HINTS['factor'][key] = saved
    run('a factorisation of b\' with a composite "prime"', forged_factor)
    def forged_point():
        RI.rank_bound(a, b, [(Fr(1), Fr(1))], [], say=lambda s: None)
    run('a point off the curve given to the descent', forged_point)
    ok = all(tests)
    say('negative tests: %d of %d refused' % (sum(tests), len(tests)))
    return ok

if __name__ == '__main__':
    if '--negative' in sys.argv:
        ok = negative()
        if WRITE: open(os.path.join(here, 'exact_last_pair_negative.txt'), 'w').write('\n'.join(OUT) + '\n')
        sys.exit(0 if ok else 1)
    main()
    if WRITE: open(os.path.join(here, 'exact_last_pair.txt'), 'w').write('\n'.join(OUT) + '\n')
