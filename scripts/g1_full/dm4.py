#!/usr/bin/env python3
"""Rank 2 and 3 on the four-fold quotient E0L: a Demjanenko-Manin argument with four points.

C: y^2 = q(t) = F0(t) GL(t) = (n^2 t^2 + m^2)(2S t^2 - 8mn t + 2S), S = m^2 + n^2, with the rational point
t1 = m/n.  E = C with origin (t1, y1); the involution y -> -y is Q -> iota(Q) = R0 - Q.

A solution t0 of the five-form system gives four rational points of C, with t-coordinates
    Q1: t0,   Q2: -t0,   Q3: 1/t0,   Q4: -1/t0
(F0 GL, F0 GU, F8 GL, F8 GU are squares at t0).  Put a_i = Q_i - R0/2 in E(Q) (x) Q, with the canonical
height as a positive definite quadratic form N and <,> its bilinear form.

E has full rational 2-torsion.  Let x be the coordinate of the model  E_T: Y^2 = x (x - al)(x - be)  in which
one particular 2-torsion point is at x = 0.  Identities in the function fields (checked below, exactly, for the
pair at hand; t, y1 = sqrt(F0 GL), y2 = sqrt(F0 GU), y3 = sqrt(F8 GL)):
    (I0)   x(2Q - R0)        = 4 Pi(t)^2 / (F0 GL),   Pi = 2mn (n^2 t^2 - m^2) + (m^4 - n^4) t
    (I+-)  x(Q1 + Q2 - R0)   = ((y1 - y2) / (2t))^2,          x(Q1 - Q2) = ((y1 + y2) / (2t))^2
    (Iinv) x(Q1 + Q3 - R0)   = ((y1 - y3) / (t^2 - 1))^2,     x(Q1 - Q3) = ((y1 + y3) / (t^2 - 1))^2
With t0 = p/q in lowest terms, M = max(p^2, q^2), H = h(t0) = (log M)/2, Psi0 = m^2 q^2 + n^2 p^2,
Psi = S (p^2 + q^2) - 4mnpq, G = 2S (p^2 + q^2), r^2 = G^2 - 64 m^2 n^2 p^2 q^2, G' = S (p^2 + q^2),
r'^2 = G'^2 - (m^2 - n^2)^2 (p^2 - q^2)^2 (r, r' integers, of either sign):
    x(2Q - R0) = 2 Pi(p,q)^2 / (Psi0 Psi),  x(Q1 +- Q2 ...) = 32 m^2 n^2 Psi0 / (G +- r),
    x(Q1 +- Q3 ...) = 2 (m^2 - n^2)^2 Psi / (G' +- r').
Elementary bounds for the heights of these fractions (gcd(p, q) = 1, pq != 0, p^2 != q^2):
    lo4  <= h(x(2Q - R0)) - 4H <= up4,   lo4 = log(n^2 (m-n)^2) - log(2 Res(Pi,Psi0)^2 Res(Pi,Psi)^2),
                                         up4 = log(2 (2 m^3 n + m^4 - n^4)^2)
    -log(2 m^2)       <= h(x(Q1 +- Q2 ...)) - 2H <= log(32 m^2 n^2 S)
    -log(8 (m+n)^2)   <= h(x(Q1 +- Q3 ...)) - 2H <= log(4 (m^2 - n^2)^2 (m+n)^2)
and for every rational point P of E_T:  -U <= h(x(P)) - hhat(P) <= BUP  (Cremona-Prickett-Siksek, Silverman; computed
in exact arithmetic as in height_certify.py).  Hence
    N(a_i) in [H - eL, H + eU],  eL = (BUP - lo4)/4, eU = (up4 + U)/4,
    |<a1, a2>|, |<a3, a4>| <= kappa_pm,   |<a1, a3>|, |<a2, a4>| <= kappa_inv,   kappa = (BUP + U + width)/4.
If the rank is at most 2, a1, a2, a4 are linearly dependent, so their Gram determinant vanishes. With x = <a1, a4>:
    N2 (N1 N4 - x^2) = N1 <a2,a4>^2 + N4 <a1,a2>^2 - 2 <a1,a2> <a2,a4> x <= kappa^2 (N1 + N4 + 2|x|).
The number N1 + N4 - 2|x| is the height of z = a1 - a4 or of z = a1 + a4, a point of E(Q) (Q1 - Q4 or Q1 + Q4 - R0), and
    (N1 + N4 - 2|x|)(N1 + N4 + 2|x|) = (N1 - N4)^2 + 4 (N1 N4 - x^2),  |N1 - N4| <= e := eL + eU,  N1 + N4 + 2|x| >= 2 (H - eL),
so  hhat(z) <= cstar(H) = (e^2/2 + 4 kappa^2) / (H - eL)   for H > eL                    (Lean: G1.fourPoint_rank2_bound).
If the rank is 3, the Gram matrix of a1, a4, a2, a3 is singular; in blocks [[A, B], [B^T, D]] (A for a1, a4; D for
a2, a3) every entry of B is one of the four bounded products.  A kernel vector shows that A - 2 kappa and D - 2 kappa
are not both positive definite (Lean: G1.fourPoint_blocks_not_both_definite), so x^2 >= (N1 - 2 kappa)(N4 - 2 kappa)
or the same for y = <a2, a3>; then |x| >= H - eL - 2 kappa, and as above
    hhat(z) <= cstar(H) = 4 kappa + (e^2 + 16 kappa e) / (4 (H - eL - kappa))   for H - eL > 2 kappa   (Lean: G1.fourPoint_rank3_bound)
for z = a1 -+ a4 or z = a2 -+ a3.  For z = a2 -+ a3 the equation below holds for Q2, whose t is -t0; both signs of
every candidate are tested.  Both bounds visibly decrease with H.  So for any Hstar: either
    (A) H < Hstar:  hhat(2 Q1 - R0) <= 4 (Hstar + eU), finitely many Q1, all listed; or
    (B) hhat(z) <= cstar(Hstar): finitely many z, and for each the equation t(Q - z) t(Q) = -1 on E,
        which has finitely many solutions; its rational solutions are found from a polynomial in t.
Every candidate t0 is tested exactly in the five-form system.

The argument for rank 2 uses only that three vectors are dependent, so it also covers rank 1 (used as a cross-check of
the pairs closed by dm3.py).

Exact arithmetic.  The constants BUP, U, eL, eU, kappa, cstar are fractions (logarithms are enclosed in rational
intervals).  The height pairing matrix G of the basis is enclosed entrywise: hhat(P) lies between
(h(x(2^k P)) - BUP)/4^k and (h(x(2^k P)) + U)/4^k, with 2^k P computed exactly.  With a rational centre G0 of the
enclosure, delta its half-width, mu a rational number such that G0 - mu is positive definite (checked by Sylvester's
criterion) and eps = r delta / mu:  v^T G v >= (1 - eps) v^T G0 v for all real v.  So {v^T G v <= b} is contained in
{v^T G0 v <= b/(1 - eps)}, which is listed exactly (box from the diagonal of the exact inverse of G0).  The two lists
are therefore supersets of the sets of the argument.  Floating point is used only to choose (the basis, k, mu, Hstar)
and for the diagnostic comparison of N(a) - H with [-eL, eU]; each choice is verified exactly where it matters.

Taken from other software without a check here: PARI's upper bound for the rank (ellrank), its saturation
(ellsaturation) and group law; Sage's lower bound lambda for the canonical height of non-torsion points
(height_function().min), whose floating-point value is used as an exact number; Kodaira symbols and Tamagawa numbers
(local_data) inside BUP.  The rank must equal the number r = 1, 2 or 3 of independent points found (PARI's lower bound
is not used: it can rest on parity).  The index of the points in E(Q)/torsion is at most sqrt(gamma^r Reg' / lambda^r),
gamma^r = 1, 4/3, 2 (Hermite), with Reg' bounded above from the enclosure; saturation is done up to that bound.
usage: dm4.py [-v] [--Hstar=H] [--upper=N] m,n [m,n ...]      (needs passagemath-schemes, passagemath-symbolics, cypari2; 15 to 60 seconds per pair)"""
import sys, json, math, time, os
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sage.all__sagemath_schemes import *
from dm import classify_t, tight_bounds, issq
from hbounds import height_bounds_exact, ln_bounds, LN2     # exact BUP, U for an integral model; rational enclosures of logarithms

def ln_int_bounds(v):
    """Rational lower and upper bounds for log(v), v a positive integer of any size."""
    v = int(v); assert v > 0
    k = v.bit_length() - 1
    if k <= 200: return ln_bounds(Fr(v))
    top = v >> (k - 120)                                # v / 2^(k-120) lies in [top, top + 1)
    a = ln_bounds(Fr(top)); b = ln_bounds(Fr(top + 1)); sh = k - 120
    return a[0] + sh * LN2[0], b[1] + sh * LN2[1]

def sqrt_hi(q):
    """A rational number >= sqrt(q), q a non-negative fraction."""
    q = Fr(q); assert q >= 0
    return Fr(math.isqrt((q.numerator << 160) // q.denominator) + 1, 1 << 80)

class Iv:
    """Closed intervals with rational ends: +, -, *."""
    def __init__(s, lo, hi=None): s.lo = Fr(lo); s.hi = Fr(lo if hi is None else hi); assert s.lo <= s.hi
    def __add__(s, o): return Iv(s.lo + o.lo, s.hi + o.hi)
    def __sub__(s, o): return Iv(s.lo - o.hi, s.hi - o.lo)
    def __mul__(s, o):
        ps = [s.lo*o.lo, s.lo*o.hi, s.hi*o.lo, s.hi*o.hi]; return Iv(min(ps), max(ps))

Rt = PolynomialRing(QQ, 't'); Kt = FractionField(Rt); t = Kt.gen()

class F4:
    """a + b*Y1 + c*Y2 + d*Y1*Y2 over Q(t), Y1^2 = r1, Y2^2 = r2"""
    r1 = None; r2 = None
    def __init__(s, a=0, b=0, c=0, d=0): s.v = (Kt(a), Kt(b), Kt(c), Kt(d))
    @staticmethod
    def co(o): return o if isinstance(o, F4) else F4(o)
    def __add__(s, o): o = F4.co(o); return F4(*[x + y for x, y in zip(s.v, o.v)])
    __radd__ = __add__
    def __neg__(s): return F4(*[-x for x in s.v])
    def __sub__(s, o): return s + (-F4.co(o))
    def __rsub__(s, o): return (-s) + o
    def __mul__(s, o):
        o = F4.co(o); a, b, c, d = s.v; e, f, g, h = o.v; r1, r2 = F4.r1, F4.r2
        return F4(a*e + b*f*r1 + c*g*r2 + d*h*r1*r2, a*f + b*e + (c*h + d*g)*r2, a*g + c*e + (b*h + d*f)*r1, a*h + d*e + b*g + c*f)
    __rmul__ = __mul__
    def conj1(s): a, b, c, d = s.v; return F4(a, -b, c, -d)
    def inv(s):
        p = s * s.conj1(); assert p.v[1] == 0 and p.v[3] == 0
        nn = p.v[0]**2 - p.v[2]**2 * F4.r2
        return s.conj1() * F4(p.v[0]/nn, 0, -p.v[2]/nn, 0)
    def __truediv__(s, o): return s * F4.co(o).inv()
    def __eq__(s, o): o = F4.co(o); return all(x == y for x, y in zip(s.v, o.v))

def issq_q(x):
    x = QQ(x); return x >= 0 and x.numerator().is_square() and x.denominator().is_square()

def hQ(x):
    x = QQ(x); return math.log(max(abs(int(x.numerator())), abs(int(x.denominator())), 1))

# Points on E_T found by Magma (public calculator, 2026-11-22): TwoDescent, FourDescent with class groups under GRH, PointsQI
# on the 4-coverings.  GRH is used only to find the points; that they lie on the curve is checked here.
extra_points = {
    (740, 539): [('420672258585150753650026931613802135666238724/227474913897287186171270349613249',
                  '-798372405517165038182538424755939465150787768190901855040094178108/3430838413082061446972975149716135176448129340193')],
}

EXT_UPPER = None                                    # an upper bound for the rank from another program (--upper=N)

def dm4(m, n, Hstar=None, verbose=False):
    T0 = time.time(); out = dict(pair=[m, n], curve='E0L', method='four points')
    S = m*m + n*n; T = Rt.gen()
    def qpoly(x): return (n*n*x**2 + m*m) * (2*S*x**2 - 8*m*n*x + 2*S)
    F0 = n*n*t**2 + m*m; F8 = m*m*t**2 + n*n; GL = 2*S*t**2 - 8*m*n*t + 2*S; GU = 2*S*t**2 + 8*m*n*t + 2*S
    t1 = QQ(m)/n; q0 = QQ(2*m*(m*m - n*n))/n
    qs = qpoly(t1 + T); e0, d, c, b, a = [qs[i] for i in range(5)]
    a1 = d/q0; a2 = c - d**2/(4*q0**2); a3 = 2*q0*b; a4 = -4*q0**2*a; a6 = a2*a4
    E0 = EllipticCurve(QQ, [a1, a2, a3, a4, a6])
    def toE0(tt, y):                                   # the point (tt, y) of C on E0 (Connell), in any field
        tau = tt - t1
        X = (2*q0*(y + q0) + d*tau) * (1/tau**2)
        Y = (4*q0**2*(y + q0) + 2*q0*(d*tau + c*tau**2) - d**2*tau**2/(2*q0)) * (1/tau**3)
        return X, Y
    def neg(P): return P[0], -P[1] - a1*P[0] - a3
    def addX(P, R):
        lam = (R[1] - P[1]) / (R[0] - P[0]); return lam*lam + a1*lam - a2 - P[0] - R[0]
    def add(P, R):
        lam = (R[1] - P[1]) / (R[0] - P[0]); nu = P[1] - lam*P[0]
        X3 = lam*lam + a1*lam - a2 - P[0] - R[0]; return X3, -(lam + a1)*X3 - nu - a3
    Pi = 2*m*n*(n*n*t**2 - m*m) + (m**4 - n**4)*t
    # ---- (I0): which 2-torsion point
    F4.r1 = qpoly(t); F4.r2 = Kt(1)
    y1 = F4(0, 1, 0, 0)
    X0 = addX(toE0(t, y1), neg(toE0(t, -y1)))
    assert X0.v[1] == 0 and X0.v[2] == 0 and X0.v[3] == 0
    e1 = None
    for P in E0.torsion_points():
        if P.order() == 2 and X0.v[0] - P[0] == 4*Pi**2/(F0*GL): e1 = P[0]
    if e1 is None: out['status'] = 'identity I0 fails'; return out
    es = sorted(P[0] for P in E0.torsion_points() if P.order() == 2 and P[0] != e1)
    al, be = es[0] - e1, es[1] - e1
    if not (al in ZZ and be in ZZ and e1 in ZZ): out['status'] = 'model not integral'; return out
    al, be = ZZ(al), ZZ(be)
    ET = EllipticCurve(QQ, [0, -(al + be), 0, al*be, 0])
    def toET(P): return ET([P[0] - e1, P[1] + (a1*P[0] + a3)/2])            # E0(Q) -> E_T(Q)
    def fromET(P): X = P[0] + e1; return E0([X, P[1] - (a1*X + a3)/2])
    # ---- (I+-), (Iinv)
    F4.r1 = qpoly(t); F4.r2 = qpoly(-t); y1 = F4(0, 1, 0, 0); y2 = F4(0, 0, 1, 0)
    Q1 = toE0(t, y1); Q2 = toE0(-t, y2); iQ2 = toE0(-t, -y2)
    ok = (addX(Q1, neg(iQ2)) - e1 == ((y1 - y2) / (2*t)) * ((y1 - y2) / (2*t))) and (addX(Q1, neg(Q2)) - e1 == ((y1 + y2) / (2*t)) * ((y1 + y2) / (2*t)))
    F4.r2 = F8*GL; y3 = F4(0, 0, 1, 0)                 # the point above 1/t is (1/t, y3/t^2)
    Q3 = toE0(1/t, y3 * (1/t**2)); iQ3 = toE0(1/t, -y3 * (1/t**2))
    ok = ok and (addX(Q1, neg(iQ3)) - e1 == ((y1 - y3) / (t*t - 1)) * ((y1 - y3) / (t*t - 1))) and (addX(Q1, neg(Q3)) - e1 == ((y1 + y3) / (t*t - 1)) * ((y1 + y3) / (t*t - 1)))
    if not ok: out['status'] = 'identity I+- or Iinv fails'; return out
    out['identities'] = 'I0, I+-, Iinv hold in the function fields'; out['model'] = [0, str(-(al + be)), 0, str(al*be), 0]
    # ---- rank and a basis
    Em = ET.minimal_model(); phi = ET.isomorphism_to(Em); psi = Em.isomorphism_to(ET)
    lo, hi, _ = tight_bounds(Em)
    if EXT_UPPER is not None and EXT_UPPER < hi: hi = EXT_UPPER; out['upper'] = 'external: --upper=%d (Magma, see magma_e0l.txt)' % EXT_UPPER
    out['rank_bounds'] = [lo, hi]
    r = hi
    if r not in (1, 2, 3): out['status'] = 'upper bound for the rank is not 1, 2 or 3'; return out
    from cypari2 import Pari
    pari = Pari(); Ep = Em.pari_curve()
    R0 = phi(toET(E0([-a2, a1*a2 - a3])))
    cand = [R0]
    for pp in range(-60, 61):                           # rational points of C above t of small height
        for qq in range(1, 61):
            if math.gcd(pp, qq) != 1: continue
            tt = QQ(pp)/qq; v = qpoly(tt)
            if tt != t1 and issq_q(v):
                P = toE0(tt, QQ(v).sqrt()); cand.append(phi(toET(E0([P[0], P[1]]))))
    for xs, ys in extra_points.get((m, n), []):      # points found by another program, on the model E_T; checked to lie on the curve
        cand.append(phi(ET([QQ(xs), QQ(ys)])))
    isom = pari.ellisomat(Ep, 0, 0)[0]
    def independent(pts):                               # greedy: a subset with nonzero regulator, smallest heights first
        ch = []
        for P in sorted(pts, key=lambda P: float(P.height())):
            if P.has_finite_order(): continue
            if float(Em.height_pairing_matrix(ch + [P]).det()) > 1e-6: ch.append(P)
            if len(ch) == r: break
        return ch
    basis = independent(cand)
    for eff in (0, 1, 2, 3):
        if len(basis) == r: break
        for tr in isom:
            F = pari.ellinit(tr[0])
            for P in pari.ellrank(F, eff)[3]:
                Pm = pari.ellisogenyapply(tr[2], P)
                if len(Pm) == 2: cand.append(Em([QQ(Pm[0]), QQ(Pm[1])]))
        basis = independent(cand)
    if len(basis) != r: out['status'] = '%d independent points not found (found %d)' % (r, len(basis)); return out
    out['rank'] = r                                     # r independent points and the upper bound r
    lam = float(Em.height_function().min(0.01, 20))
    if lam <= 0: out['status'] = 'no height lower bound'; return out
    lamq = Fr(lam)                                      # the exact value of that float
    # ---- exact arithmetic from here on: constants, enclosures of height matrices, the two lists
    BUP, U = height_bounds_exact([0, -(al + be), 0, al*be, 0], ET)          # fractions
    Aq = QQ(-(al + be)); Bq = QQ(al*be); b2 = 4*Aq; b4 = 2*Bq; b8 = -Bq*Bq
    def enclose(Pm):                                    # an interval with rational ends that contains hhat(Pm)
        x = QQ(psi(Pm)[0]); hf = float(Pm.height()); k = 0
        while 4**k * hf < 2e6:
            x = (x**4 - b4*x*x - b8) / (4*x**3 + b2*x*x + 2*b4*x); k += 1
        hx = ln_int_bounds(max(abs(int(x.numerator())), int(x.denominator())))
        return (hx[0] - BUP) / 4**k, (hx[1] + U) / 4**k
    def gram_enclosure(pts):                            # entrywise enclosure [lo, hi] of the height pairing matrix
        nn = len(pts); dd = [enclose(P) for P in pts]
        lo = [[None]*nn for _ in range(nn)]; hi = [[None]*nn for _ in range(nn)]
        for i in range(nn):
            lo[i][i], hi[i][i] = dd[i]
            for j in range(i):
                ee = enclose(pts[i] + pts[j])
                lo[i][j] = lo[j][i] = (ee[0] - dd[i][1] - dd[j][1]) / 2; hi[i][j] = hi[j][i] = (ee[1] - dd[i][0] - dd[j][0]) / 2
        return lo, hi
    def det_iv(lo, hi):                                 # interval of the determinant, by interval arithmetic
        M = [[Iv(lo[i][j], hi[i][j]) for j in range(len(lo))] for i in range(len(lo))]
        if len(M) == 1: return M[0][0]
        if len(M) == 2: return M[0][0]*M[1][1] - M[0][1]*M[1][0]
        return (M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1]) - M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0]) + M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0]))
    lo0, hi0 = gram_enclosure(basis); reg0 = det_iv(lo0, hi0)
    if not reg0.lo > 0: out['status'] = 'independence of the points not certified'; return out
    herm = {1: Fr(1), 2: Fr(4, 3), 3: Fr(2)}[r]         # gamma_r^r (Hermite): Reg >= lam^r / gamma_r^r
    kmax = int(math.floor(math.sqrt(float(reg0.hi * herm) / lam**r))) + 1
    while not kmax * kmax * lamq**r >= reg0.hi * herm: kmax += 1          # index^2 <= Reg' gamma^r / lam^r <= kmax^2, exactly
    sat = pari.ellsaturation(Ep, [[P[0], P[1]] for P in basis], kmax)
    gens = [Em([QQ(P[0]), QQ(P[1])]) for P in sat]
    def gram_of(g): return [[float(v) for v in row] for row in Em.height_pairing_matrix(g).rows()]
    for _ in range(200):                                # pairwise size reduction (floating point; only chooses the basis)
        gens.sort(key=lambda P: float(P.height())); Gm = gram_of(gens); ch = False
        for i in range(r):
            for j in range(i):
                k = int(round(Gm[i][j] / Gm[j][j]))
                if k: gens[i] = gens[i] - k*gens[j]; ch = True; break
            if ch: break
        if not ch: break
    Gm = gram_of(gens); reg = float(matrix(RDF, Gm).det())
    lo1, hi1 = gram_enclosure(gens)
    one60 = 1 << 60
    G0 = [[Fr(int(math.floor((lo1[i][j] + hi1[i][j]) / 2 * one60)), one60) for j in range(r)] for i in range(r)]      # a dyadic centre
    for i in range(r):
        for j in range(i): G0[i][j] = G0[j][i]
    delta = max(max(hi1[i][j] - G0[i][j], G0[i][j] - lo1[i][j]) for i in range(r) for j in range(r))
    assert all(lo1[i][j] <= Fr(Gm[i][j]) <= hi1[i][j] for i in range(r) for j in range(r))       # the floating-point matrix lies in the enclosure
    mu = Fr(min(float(v) for v in matrix(RDF, Gm).eigenvalues()) * 0.99)
    def det_fr(M):
        if len(M) == 1: return M[0][0]
        if len(M) == 2: return M[0][0]*M[1][1] - M[0][1]*M[1][0]
        return (M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1]) - M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0]) + M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0]))
    Gmu = [[G0[i][j] - (mu if i == j else 0) for j in range(r)] for i in range(r)]
    if not (mu > 0 and all(det_fr([row[:k] for row in Gmu[:k]]) > 0 for k in range(1, r + 1))):     # G0 - mu is positive definite (Sylvester)
        out['status'] = 'centre matrix not positive definite'; return out
    eps = r * delta / mu                                # v^T G v >= (1 - eps) v^T G0 v for the true matrix G and every real v
    if not eps <= Fr(1, 100): out['status'] = 'Gram matrix not enclosed tightly enough'; return out
    G0inv = matrix(QQ, [[QQ(int(v.numerator)) / QQ(int(v.denominator)) for v in row] for row in G0]).inverse()
    G0inv_d = [Fr(int(G0inv[i, i].numerator()), int(G0inv[i, i].denominator())) for i in range(r)]
    out['gens'] = 'points saturated at primes <= %d (height lower bound %.4f, regulator of the points before saturation at most %.3f)' % (kmax, lam, float(reg0.hi)); out['kmax'] = kmax; out['lam'] = lam; out['reg0'] = float(reg0.hi)
    out['min_model'] = [str(v) for v in Em.ainvs()]; out['G'] = [[str(P[0]), str(P[1])] for P in gens]; out['gram'] = Gm; out['reg'] = reg
    out['R0'] = [str(R0[0]), str(R0[1])]
    out['gram_check'] = dict(max_deviation=float(delta), mu=float(mu), relative=float(eps))
    tors = Em.torsion_points(); out['ntors'] = len(tors)
    import itertools
    G0i = [[int(G0[i][j] * one60) for j in range(r)] for i in range(r)]       # G0 is dyadic: these are integers
    def box(bound, centre2=None):
        """All integer vectors c with (c - centre2/2)^T G (c - centre2/2) <= bound for the true height matrix G
        (centre2 an integer vector, default 0): a superset, listed with the rational centre matrix G0 and the bound
        divided by 1 - eps.  Integer arithmetic: with u = 2c - centre2 the condition is u^T G0 u <= 4 bound/(1 - eps)."""
        centre2 = centre2 or [0]*r; bp = Fr(bound) / (1 - eps); rng = []
        thr = (4 * bp * one60).numerator // (4 * bp * one60).denominator      # floor: u^T G0i u is an integer
        for i in range(r):
            sh = sqrt_hi(bp * G0inv_d[i])               # |c_i - centre2_i/2|^2 <= bp (G0^-1)_ii
            rng.append(range(int(math.floor(Fr(centre2[i], 2) - sh)), int(math.ceil(Fr(centre2[i], 2) + sh)) + 1))
        for cvec in itertools.product(*rng):
            u = [2*cvec[i] - centre2[i] for i in range(r)]
            if sum(G0i[i][j]*u[i]*u[j] for i in range(r) for j in range(r)) <= thr: yield cvec
    def comb(cvec):
        P = Em(0)
        for cc, g in zip(cvec, gens): P = P + cc*g
        return P
    rho = None
    for cvec in box(enclose(R0)[1]):
        if (R0 - comb(cvec)).has_finite_order(): rho = cvec
    if rho is None: out['status'] = 'R0 not in the group generated'; return out
    out['rho'] = list(rho)
    # ---- the constants of the argument, as fractions (upper bounds for eL, eU, kappa)
    octic = m**8 + 14*m**4*n**4 + n**8
    res0 = m*m*n*n*octic; res1 = (m*m - n*n)**2*octic   # the resultants of Pi with Psi0 and Psi (G1.fourPoint_res_psi0_dvd, G1.fourPoint_res_psi_dvd)
    Rpq = PolynomialRing(ZZ, 'p'); p = Rpq.gen()
    Pi_p = 2*m*n*(n*n*p**2 - m*m) + (m**4 - n**4)*p; Psi0_p = n*n*p**2 + m*m; Psi_p = S*p**2 - 4*m*n*p + S
    assert res0 == abs(int(Pi_p.resultant(Psi0_p))) and res1 == abs(int(Pi_p.resultant(Psi_p))) and res0 > 0 and res1 > 0
    lnlo = lambda v: ln_int_bounds(v)[0]; lnhi = lambda v: ln_int_bounds(v)[1]
    lo4 = lnlo(n*n*(m - n)**2) - lnhi(2*res0**2*res1**2)
    up4 = lnhi(2*(2*m**3*n + m**4 - n**4)**2)
    wpm = lnhi(32*m*m*n*n*S) + lnhi(2*m*m)
    winv = lnhi(4*(m*m - n*n)**2*(m + n)**2) + lnhi(8*(m + n)**2)
    eL = (BUP - lo4) / 4; eU = (up4 + U) / 4; ee = eL + eU
    kap = (BUP + U + max(wpm, winv)) / 4
    if r <= 2:                                          # three vectors in a space of dimension at most 2
        def cstar(H): return (ee*ee/2 + 4*kap*kap) / (H - eL)
        def valid(H): return H > eL
        target = Fr(9, 10) * min(lo1[i][i] for i in range(r))       # below the smallest height of a basis vector
        guess = eL + (ee*ee/2 + 4*kap*kap) / target
    else:
        def cstar(H): return 4*kap + (ee*ee + 16*kap*ee) / (4*(H - eL - kap))
        def valid(H): return H - eL > 2*kap
        guess = eL + kap + (ee*ee + 16*kap*ee) / (Fr(12, 5)*kap)    # then cstar = 4.6 kappa
    if Hstar is None: Hstar = int(math.ceil(guess)) + 1
    Hstar = Fr(Hstar)
    if not valid(Hstar): out['status'] = 'Hstar too small'; return out
    cs = cstar(Hstar)                                   # cstar decreases with H (see the header), so it bounds hhat(z) for every H >= Hstar
    fl = lambda v: float(v)
    out.update(BUP=fl(BUP), U=fl(U), res=[str(res0), str(res1)], lo4=fl(lo4), up4=fl(up4), wpm=fl(wpm), winv=fl(winv), eL=fl(eL), eU=fl(eU), kappa=fl(kap), Hstar=fl(Hstar), cstar=fl(cs))
    out['arithmetic'] = 'exact (rational) for the constants, the enclosures of the height matrices and the two lists'
    # ---- the inverse map E -> t
    A1c = 4*q0**2; B1c = 2*d; A2c = d**2 - 4*q0**2*c; B2c = -4*q0**2*b; C2c = -4*q0**2*a; denc = A2c*B1c - A1c*B2c
    assert denc != 0
    R0E0 = E0([-a2, a1*a2 - a3])
    def t_of(Pm):                                       # Pm on Em; returns t, or None for t = infinity
        if Pm.is_zero(): return t1
        P = fromET(psi(Pm))
        if P == R0E0: return t1
        Ps = R0E0 - P; S1 = P[0] + Ps[0]; S2 = P[0]*Ps[0]
        sg = (A2c*S1 - A1c*S2 + A1c*C2c) / denc        # sigma = 1/(t - t1), from the fibre {P, R0 - P} as in dm3.py
        assert A1c*sg**2 + B1c*sg == S1 and A2c*sg**2 + B2c*sg + C2c == S2
        if sg == 0: return None
        tt = t1 + 1/sg
        assert issq_q(qpoly(tt))                        # the point lies above this t
        return tt
    verdicts = {}
    def test_t(tt):
        if tt is None: return
        for sgn in (1, -1):
            tf = Fr(sgn*int(QQ(tt).numerator()), int(QQ(tt).denominator()))
            if str(tf) not in verdicts: verdicts[str(tf)] = classify_t(m, n, tf)
    for tt in (0, 1, -1, t1, 1/t1): test_t(QQ(tt))
    # ---- (A) all Q with hhat(2Q - R0) <= 4 (Hstar + eU), that is, (c - rho/2)^T G (c - rho/2) <= Hstar + eU
    obs = [0.0, 0.0]; na = 0
    for cvec in box(Hstar + eU, list(rho)):
        Q = comb(cvec); u = [2*cvec[i] - rho[i] for i in range(r)]
        Na = sum(Gm[i][j]*u[i]*u[j] for i in range(r) for j in range(r)) / 4
        for Tt in tors:
            na += 1; tq = t_of(Q + Tt); test_t(tq)
            if tq is not None and tq != 0:              # numerical check of the bounds for N(a) on actual points
                dv = Na - hQ(tq); obs[0] = min(obs[0], dv); obs[1] = max(obs[1], dv)
    out['nA'] = na; out['observed N(a) - H'] = obs
    if not (-float(eL) <= obs[0] and obs[1] <= float(eU)): out['status'] = 'observed N(a) - H outside the bounds'; return out
    # ---- (B) all z with hhat(z) <= cstar, and the equation t(Q - z) t(Q) = -1
    zs = [comb(cvec) + Tt for cvec in box(cs) for Tt in tors]
    out['nz'] = len(zs)
    qt = qpoly(t)
    class F2:                                   # alpha + beta*y in Q(t)[y]/(y^2 - q(t))
        def __init__(s, al, be=0): s.al = Kt(al); s.be = Kt(be)
        def __add__(s, o): o = o if isinstance(o, F2) else F2(o); return F2(s.al + o.al, s.be + o.be)
        __radd__ = __add__
        def __neg__(s): return F2(-s.al, -s.be)
        def __sub__(s, o): return s + (-(o if isinstance(o, F2) else F2(o)))
        def __rsub__(s, o): return (-s) + o
        def __mul__(s, o):
            o = o if isinstance(o, F2) else F2(o)
            return F2(s.al*o.al + s.be*o.be*qt, s.al*o.be + s.be*o.al)
        __rmul__ = __mul__
        def inv(s):
            nn = s.al**2 - s.be**2*qt; return F2(s.al/nn, -s.be/nn)
        def __truediv__(s, o): o = o if isinstance(o, F2) else F2(o); return s * o.inv()
    yy = F2(0, 1); Qg = toE0(t, yy)
    two = [P for P in tors if P.order() == 2]
    def solve(z, const):                               # candidates t for t(Q - z) t(Q) = const, Q on E
        zz = fromET(psi(z)); zc = (QQ(zz[0]), QQ(zz[1]))
        D = add(Qg, neg(zc))                            # Q - z on E0
        f = t1 + (2*q0*(D[0] + c) - d**2/(2*q0)) / D[1]             # t(Q - z) = alpha + beta*y
        A = f.al*t - const; Bq = f.be*t
        Nn = A*A - Bq*Bq*qt                             # norm of t(Q - z) t(Q) - const
        if Nn == 0: return None
        return [rr for pol in (Rt(Nn.numerator()), Rt(Nn.denominator()), Rt(f.al.denominator()), Rt(f.be.denominator())) for rr, _ in pol.roots(QQ)]
    # self-test of the solver: for known points Q and z, t(Q) must be among the candidates for the constant t(Q - z) t(Q)
    gl = gens[-1]                                       # the last basis vector (the same as the first in rank 1)
    for Qs, zs_ in ((gens[0] + gl + (gens[0] if r == 1 else 0*gl), gens[0]), (2*gens[0] - gl + tors[1] + (3*gens[0] if r == 1 else 0*gl), gl + tors[2 % len(tors)]), (3*gl + tors[3 % len(tors)], tors[1])):
        ta, tb = t_of(Qs), t_of(Qs - zs_)
        assert ta is not None and tb is not None and ta in solve(zs_, ta*tb), 'solver self-test'
    out['solver_selftest'] = 'ok'
    nroots = 0
    for z in zs:
        for w in [z, -z, z + R0, z - R0, R0 - z, -z - R0] + [z + P for P in two] + [P for P in two]: test_t(t_of(w))     # the special positions
        if z.is_zero(): continue                        # t(Q)^2 = -1 has no rational solution
        cands = solve(z, -1)
        if cands is None: out['status'] = 'equation identically satisfied for a z'; return out
        for rr in cands: nroots += 1; test_t(rr)
    out['nroots'] = nroots
    out['nt'] = len(verdicts)
    out['t'] = {k: v for k, v in verdicts.items() if v != 'no' and len(k) < 40}
    out['status'] = 'ALERT' if 'SOLUTION' in verdicts.values() else 'closed'
    out['sec'] = '%.1f' % (time.time() - T0)
    return out

if __name__ == '__main__':
    Hs = None
    for s in sys.argv[1:]:
        if s.startswith('--Hstar='): Hs = float(s[8:]); continue
        if s == '-v': continue
        if s.startswith('--upper='): EXT_UPPER = int(s[8:]); continue
        m, n = (int(v) for v in s.split(','))
        try: r = dm4(m, n, Hs, verbose=('-v' in sys.argv))
        except Exception as ex:
            import traceback; traceback.print_exc(); r = dict(pair=[m, n], curve='E0L', status='error ' + str(ex)[:120])
        print(json.dumps(r, default=str), flush=True)
