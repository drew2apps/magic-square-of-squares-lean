#!/usr/bin/env python3
"""Fisher's fifth descent pairing Theta_4(y0, y1) on E+ of (72,65) (it governs the census ratios (72,65), (137,7)),
computed in exact arithmetic following T. Fisher, "Higher descents on an elliptic curve with a rational 2-torsion
point", arXiv:1509.03234 (Sect. 2: Thm 2.1, Lemma 2.2; Sect. 3: the pairing formula; Sect. 7: Lemma 7.1, Thms 7.2-7.4).

    python3 rank_fisher_theta4.py              the computation (about a minute)
    python3 rank_fisher_theta4.py --negative   altered inputs that must be refused
    (add --write to rewrite rank_fisher_theta4.txt / rank_fisher_theta4_negative.txt)

Written by an AI model (Claude); not read by a person.  Reviewed by a separate AI agent on 26 November 2026: holds
with gaps (LEDGER.md, tier C).  Not counted in the census.

NOTATION.  E: y^2 = x(x^2 + a x + b), T = (0,0), phi: E -> E' with kernel <T>; E': Y^2 = X(X^2 - 2aX + a^2 - 4b),
T' = (0,0) = generator of ker(phihat).  S_1 = S^(phi)(E/Q) (rank_indep's "E' side"), S'_1 = S^(phihat)(E'/Q)
(rank_indep's "E side", classes d | b).  rank_fisher.py: dim S_4 = 2 = <known, y0>, dim S'_4 = 3 = <known, y1>,
y0 = 2, y1 = 2, 4-descent bound exactly 3.  Theta_4 is a pairing S_4 x S'_4 with the images of the rational points
in both kernels (Thm 2.1), so it is determined by the single value Theta_4(y0, y1); if it is non-zero the left and
right kernels S_5, S'_5 have dimensions 1 and 2 and Lemma 2.2 (n = 3) gives rank <= 1 + 2 - 2 = 1.

THE CURVES.  C3' (rank_fisher5.py): {QA = QB = 0} in P^3(c:d:e:g), a level-3 lift of y0 (a torsor under E, a
phi-covering of the level-2 curve D_lam, lam = -1, xi = 7), checked there to be isomorphic to rank_fisher's C3.
Its pencil has the cone pairs over L = Q(sqrt -2) and L' = Q(sqrt 2); l_1 (over L) and l_3 (over L') are the
tangent planes at points of the cone conics, f = l_1 conj(l_1), f_dag = l_3 conj(l_3), q the quadric through
a_1 + ... + a_4 and c with f f_dag = c q^2 modulo I(C3') (rank_fisher5 steps (6)-(8), re-run here).
C4 = {QA = QB = 0, Xi W^2 = f0} in P^4(c:d:e:g:W), f0 = f / s0 primitive integral; equivalently xi4 W^2 = f with
xi4 = Xi s0.  It is the double cover of C3' attached to the pushout form f (Fisher Sect. 3 Selmer set; Thm 7.3:
a phihat-covering, so C4 is a torsor under E' and a level-4 lift of y0 when it is everywhere locally soluble).

THE PUSHOUT ON C4 (Fisher Thm 7.4: H - H_dag represents T'; the explicit function is derived here and its divisor
is CHECKED, not taken from the paper):
    k_1 in L with N(k_1) = xi4 mod squares, sigma_1 = sqrt(xi4 N(k_1)) in Q;
    k_3 in L' with N(k_3) = c xi4 mod squares, tau_3 = sqrt(c / (xi4 N(k_3))) in Q;
    l_1 / k_1 = (Tr1 + delta1 sqrt(-2)) / 2, l_3 / k_3 = (Tr3 + delta3 sqrt 2) / 2 (Tr, delta rational linear forms);
    N1 = Tr1 + 2 xi4 W / sigma_1,  N1' = Tr1 - 2 xi4 W / sigma_1        (linear forms on C4)
    R3 = W Tr3 + 2 tau_3 q,        R3' = W Tr3 - 2 tau_3 q              (quadratic forms on C4)
    g0 = R3 / (W N1)                                                    (a rational function on C4 over Q)
(N1 = (u_1 + conj u_1)^2 with u_1^2 = l_1/k_1; R3 / W = (u_3 + conj u_3)^2 with u_3^2 = l_3/k_3: v9's formula with
the denominator W cleared.)  EXACT CHECKS:
    (I1) N1 N1' = -2 delta1^2 and R3 R3' = 2 W^2 delta3^2 modulo I(C4) = (QA, QB, xi4 W^2 - f)
         (the second uses f f_dag - c q^2 in (QA, QB), re-checked by exact linear algebra);
    (I2) l_1, conj(l_1) have no common zero on C3' (resultant of QA, QB on the line {l_1 = conj l_1 = 0} != 0); the
         same for l_3; delta3 has no zero on A_1 (resultant on the line {l_1 = delta3 = 0} over L, != 0);
    (I3) div(l_1 | C3') = 2 A_1 (tangency, rank_fisher5 (6)), q >= A_1 + A_2 + A_3 + A_4 (rank_fisher5 (8));
    (I3') A_1 is two distinct points (the ruling {l_1 = 0} n cone meets C3' with non-zero discriminant), so q
         vanishes simply at each point of A_1 + A_2 and the multiplicity m below is 1.
WRITTEN ARGUMENT (evenness).  Let P in C4 (normalisation) over x0 in C3'; C4 -> C3' is etale.  From (I1):
ord N1 + ord N1' = 2 ord_x0 delta1 and N1 + N1' = 2 Tr1; if Tr1(x0) != 0 one of them is a unit, so ord N1 is 0 or
2 ord delta1; if Tr1(x0) = 0 then delta1(x0) != 0 by (I2) and both are units.  So ord_P N1 is even.  For R3: if
W(P) != 0 the same argument with R3 R3' = 2 W^2 delta3^2, R3 + R3' = 2 W Tr3 and (I2).  If W(P) = 0 then x0 is in
A_1 + A_2 with multiplicity m (ord_P W = m by (I3), since div f = 2(A_1 + A_2)); R3 - R3' = 4 tau_3 q has order
>= m (I3), R3 + R3' = 2 W Tr3 has order >= m, so both have order >= m, and ord R3 + ord R3' = 2m + 2 ord delta3(x0)
= 2m by (I2); so ord R3 = m and ord (R3 / W) = 0.  Hence div g0 = 2 B with B a Q-rational divisor of degree 0,
2B ~ 0, so [B] is in E'(Q)[2] = {0, T'} (b is not a square).  (I4) [B] != 0: two points P, P' of C4(Q_p) with
g0(P) g0(P') not a square (if B = div k then g0 = const k^2).  So g0 is a pushout form for T' on C4 (Fisher
Sect. 3) and, with y1 in S'_1 = S^(phihat)(E'/Q),
    Theta_4(y0, y1) = <[C4], y1> = sum_v (g0(P_v), y1)_v      (Fisher Sect. 3; P_v in C4(Q_v)).
For y1 = 2 the term at v is (g0(P_v), 2)_v: 0 at R and at p = +-1 mod 8; (-1)^(v_p) at p = +-3 mod 8.

BAD PLACES.  C3' is smooth modulo p for p outside B0 = {2, 3, 5, 7, 13, 97, 137, 710053, 79553041, 97505521}: the
discriminant of the pencil quartic det(x A + y B) is checked to be a product of these primes (the invariant I of the
quartic is lambda^4 c4(E), lambda = 2^17 3^2 5^4 7^3 13 97^6 137^3 710053^3, the model not being minimal).
WRITTEN ARGUMENT (good p).  Let p be odd, outside B0, p not dividing Xi.  The pencil, its cones, the rulings through
the reduced cone points and the tangent planes reduce (the pencil discriminant is a unit), so f0 mod p is, up to a
unit, the product of two tangent planes to the reduced cones and has divisor 2(A1~ + A2~) on the smooth genus-one
curve C3~ with A1~ + A2~ - H~ in the class of T~ != 0 (p odd, good reduction).  So the normalisation of C4 mod p is
a connected etale double cover of C3~: a smooth genus-one curve, with >= p + 1 - 2 sqrt p points over F_p.
  * ELS: at most 8 of them have W = 0; for p >= 17 one has W != 0, it is a smooth point of the model in P^4 (C3~
    smooth there and d/dW (Xi W^2 - f0) = 2 Xi W != 0) and lifts by Hensel.  So C4(Q_p) is non-empty.
  * The local term: let N1^0, R3^0 be the primitive integral multiples of N1, R3 (g0^0 = R3^0 / (W N1^0) differs
    from g0 by a rational constant, which does not change the sum: product formula).  N1^0 does not vanish on C4~
    (a + b W = 0 on C4~ would make f0~ = Xi (a/b)^2 a square, or a = 0 on C3~); R3^0 = W T3 + r Q (T3, Q primitive)
    does not vanish on C4~ when v_p(r) >= 0 (squaring W T3 = -r Q on C4~ would make f0~ a constant times a square,
    or T3 = 0, on C3~), and for v_p(r) < 0 it reduces to Q~, which vanishes on C4~ iff Q~ lies in span(QA~, QB~):
    the primes of den(r) where it does (all of them here: 2, 7, 353, 1279, 336527, 252169271) are put into S, and the
    terms there are computed.  Their zeros and those of W number at most
    8 + 16 + 8 = 32 points; for p >= 47 a point of C4~(F_p) avoids them and lifts, and there g0^0 is a p-adic unit.
So the sum runs over S = {R, 2, every p < 47} u B0 u primes(Xi) u {p | den(r): Q0 mod p in span(QA, QB)} (and the
primes of the control classes y); at each p of S points of C4(Q_p) are found and the term computed with PARI's
p-adic numbers (precision asserted); rank_fisher_theta4_recheck.py re-derives every term by Hensel's lemma in integer
arithmetic.  WHAT IS A WRITTEN ARGUMENT (not machine-checked): the evenness of div g0 from (I1)-(I3'), the
good-reduction arguments above, and the use of Fisher's Thm 2.1, Lemma 2.2 and Sect. 3 formula (read from the paper
through WebFetch summaries).

CONTROLS ON THE SAME CURVE (known answers).  (K1) The term at each place is computed at several local points and must
not depend on the point (Sect. 3).  (K2) For z in the span of the images of E(Q) (classes b and x(G)), the same sum
with z in place of y1 must be 0 (Thm 2.1: the images of the rational points lie in the kernels).  (K3) The value must
not depend on the level-4 lift: the same computation for a second twist Xi' = Xi z (z in S'_1) with new k_1, k_3.
(K4, rank_fisher_theta4_k4.py) the three other level-3 lifts xi' = 14, -7, -14 over the same D_lam, each with new
conic points: the same value.  EXTERNAL CONTROLS (rank_fisher_theta4_control.py): four rank-0 curves of the same
shape where the answer is predicted by the BSD order of Sha (two with Theta_4 != 0, two with Theta_4 = 0).
"""
import os, sys, json, random
from fractions import Fraction as Fr
from math import isqrt, gcd
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import sympy as sp
import rank_fisher5 as R5
import rank_second as R2
import rank_indep as RI
import cypari2
pari = cypari2.Pari()
pari.allocatemem(1 * 10 ** 9)
from mpmath import iv
iv.dps = 300

WRITE = '--write' in sys.argv
OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(s)

c, d, e, g = R5.c, R5.d, R5.e, R5.g
W = sp.Symbol('W'); t = sp.Symbol('t')
VS = (c, d, e, g)
B0 = [2, 3, 5, 7, 13, 97, 137, 710053, 79553041, 97505521]
def prec(p): return 400 if p < 10 ** 4 else 200

def sqfree(n):
    """squarefree part of a non-zero rational (PARI core of numerator * denominator), sign kept."""
    n = Fr(n); a, b = n.numerator * n.denominator, 1
    return int(pari.core(a))

def is_rat_square(x):
    x = Fr(x)
    return x >= 0 and isqrt(x.numerator) ** 2 == x.numerator and isqrt(x.denominator) ** 2 == x.denominator

def rat_sqrt(x):
    x = Fr(x); assert is_rat_square(x), 'not a rational square'
    return Fr(isqrt(x.numerator), isqrt(x.denominator))

def split_t(expr, D, m):
    """expr = A + B t with t = m sqrt(D): return (A, m B) so that expr = A + (m B) sqrt(D)."""
    P = sp.Poly(sp.expand(expr), t)
    assert P.degree() <= 1
    A = sp.expand(P.coeff_monomial(1)); Bt = sp.expand(P.coeff_monomial(t))
    return A, sp.expand(m * Bt)

def poly_int(expr, vs):
    """primitive integral multiple of a form: (content, primitive poly) with expr = content * prim."""
    P = sp.Poly(sp.expand(expr), *vs)
    den = sp.ilcm(*[sp.Rational(k).q for k in P.coeffs()])
    num = sp.igcd(*[int(sp.Rational(k) * den) for k in P.coeffs()])
    prim = sp.Poly(sp.expand(expr * den / num), *vs)
    assert all(sp.Rational(k).q == 1 for k in prim.coeffs())
    return sp.Rational(num, den), prim

def in_ideal_C3(expr, QA, QB):
    """expr (a form in c,d,e,g of degree k >= 2) lies in (QA, QB): exact linear algebra (PARI matsolve over Q)."""
    P = sp.Poly(sp.expand(expr), *VS)
    if P.is_zero: return True
    k = P.total_degree()
    mons = sorted(sp.itermonomials(VS, k - 2, k - 2), key=str) if k > 2 else [sp.Integer(1)]
    mons = [m_ for m_ in mons if sp.Poly(m_, *VS).total_degree() == k - 2]
    cols = [sp.Poly(sp.expand(QA * m_), *VS) for m_ in mons] + [sp.Poly(sp.expand(QB * m_), *VS) for m_ in mons]
    tm = sorted(set(sum([list(C.monoms()) for C in cols], []) + list(P.monoms())))
    M = [[str(sp.Rational(C.coeff_monomial(mo_ if False else sp.prod([v_ ** k_ for v_, k_ in zip(VS, mo_)]))))
          for C in cols] for mo_ in tm]
    rhs = [str(sp.Rational(P.coeff_monomial(sp.prod([v_ ** k_ for v_, k_ in zip(VS, mo_)])))) for mo_ in tm]
    Mp = pari.matrix(len(tm), len(cols), [x for row in M for x in row])
    vp_ = pari('[' + ','.join(rhs) + ']~')
    # solvable iff rank(M) == rank([M | rhs])
    return int(pari.matrank(Mp)) == int(pari.matrank(pari.matconcat([Mp, vp_])))

def reduce_W(expr, xi4, f):
    """reduce a polynomial in c,d,e,g,W modulo xi4 W^2 - f: returns (A, B) with expr = A + W B mod that."""
    P = sp.Poly(sp.expand(expr), W)
    A, B = sp.Integer(0), sp.Integer(0)
    for (k,), co in zip(P.monoms(), P.coeffs()):
        term = sp.expand(co * (f / xi4) ** (k // 2))
        if k % 2: B += term
        else: A += term
    return sp.expand(A), sp.expand(B)

def line_resultant(lin1, lin2, QA, QB, red=None):
    """the line {lin1 = lin2 = 0} (lin over Q or over Q(t) with reduction red): parametrise it by two independent
    kernel vectors and return the resultant of QA, QB restricted (a binary pair); non-zero iff no common zero."""
    rd = red if red else (lambda z: sp.expand(z))
    rows = [[rd(sp.Poly(l_, *VS).coeff_monomial(v_)) for v_ in VS] for l_ in (lin1, lin2)]
    M = sp.Matrix(rows)
    if red:
        # work over Q(t): kernel by hand via 2x2 minors (rows independent assumed, checked)
        # choose two columns with a non-zero minor
        piv = None
        for i in range(4):
            for j in range(i + 1, 4):
                mnr = rd(M[0, i] * M[1, j] - M[0, j] * M[1, i])
                if mnr != 0: piv = (i, j, mnr); break
            if piv: break
        assert piv, 'the two forms are dependent'
        i, j, mnr = piv
        free = [k for k in range(4) if k not in (i, j)]
        basis = []
        for fk in free:
            v = [0] * 4; v[fk] = mnr
            # solve M[:, (i,j)] (vi, vj)^T = -M[:, fk] mnr  by Cramer
            b0, b1 = -M[0, fk], -M[1, fk]
            v[i] = rd(b0 * M[1, j] - M[0, j] * b1); v[j] = rd(M[0, i] * b1 - b0 * M[1, i])
            basis.append(v)
    else:
        ns = M.nullspace(); assert len(ns) == 2, 'the two forms are dependent'
        basis = [list(v) for v in ns]
    s_, u_ = sp.symbols('s_ u_')
    pt = {v_: s_ * basis[0][k] + u_ * basis[1][k] for k, v_ in enumerate(VS)}
    qa = sp.Poly(rd(sp.expand(QA.subs(pt, simultaneous=True))), s_, u_)
    qb = sp.Poly(rd(sp.expand(QB.subs(pt, simultaneous=True))), s_, u_)
    A_ = [rd(qa.coeff_monomial(m_)) for m_ in (s_ ** 2, s_ * u_, u_ ** 2)]
    B_ = [rd(qb.coeff_monomial(m_)) for m_ in (s_ ** 2, s_ * u_, u_ ** 2)]
    # resultant of a0 s^2 + a1 s u + a2 u^2 and b0 s^2 + b1 s u + b2 u^2
    a0, a1, a2 = A_; b0, b1, b2 = B_
    res = rd((a0 * b2 - a2 * b0) ** 2 - (a0 * b1 - a1 * b0) * (a1 * b2 - a2 * b1))
    return res

# ---------------------------------------------------------------- local points (PARI p-adics) and the real place
class Model:
    """C3' = {PA(c,d) + QAeg(e,g) = 0, PB(c,d) + QBeg(e,g) = 0} (block diagonal pencil) and C4: Xi W^2 = f0."""
    def __init__(self, QA, QB):
        self.QA, self.QB = QA, QB
        PA = sp.Poly(QA, *VS); PB = sp.Poly(QB, *VS)
        for P_ in (PA, PB):
            for mo in P_.monoms():
                assert not ((mo[0] or mo[1]) and (mo[2] or mo[3])), 'the pencil is not block diagonal'
        cof = lambda P_, m_: int(P_.coeff_monomial(m_))
        self.A_cd = [cof(PA, c ** 2), cof(PA, c * d), cof(PA, d ** 2)]
        self.A_eg = [cof(PA, e ** 2), cof(PA, e * g), cof(PA, g ** 2)]
        self.B_cd = [cof(PB, c ** 2), cof(PB, c * d), cof(PB, d ** 2)]
        self.B_eg = [cof(PB, e ** 2), cof(PB, e * g), cof(PB, g ** 2)]

    @staticmethod
    def bq(F, x, y): return F[0] * x * x + F[1] * x * y + F[2] * y * y

    def points(self, p, cd_list):
        """points of C3'(Q_p) over the given integer pairs (c, d): (e : g) from the quadratic
        PB QAeg - PA QBeg = 0, then the scaling.  Returns PARI p-adic points [c, d, e, g]."""
        N = prec(p)
        out = []
        for (c0, d0) in cd_list:
            if isinstance(c0, int): c0 = pari('%d + O(%d^%d)' % (c0, p, N))
            if isinstance(d0, int): d0 = pari('%d + O(%d^%d)' % (d0, p, N))
            PA = self.bq(self.A_cd, c0, d0); PB = self.bq(self.B_cd, c0, d0)
            if PA == 0: continue
            qa = [PB * self.A_eg[k] - PA * self.B_eg[k] for k in range(3)]   # qa0 e^2 + qa1 e g + qa2 g^2 = 0
            if qa[0] == 0: continue
            Dp = qa[1] ** 2 - 4 * qa[0] * qa[2]
            if Dp == 0 or int(pari.padicprec(Dp, p)) - int(pari.valuation(Dp, p)) < 4: continue
            if not pari.issquare(Dp): continue
            sD = pari.sqrt(Dp)
            for sg in (1, -1):
                rho = (-qa[1] + sg * sD) / (2 * qa[0])
                den = self.bq(self.A_eg, rho, 1)
                if den == 0: continue
                G = -PA / den
                if G == 0 or int(pari.padicprec(G, p)) - int(pari.valuation(G, p)) < 4: continue
                if not pari.issquare(G): continue
                gg = pari.sqrt(G)
                out.append([c0, d0, rho * gg, gg])
        return out

    def real_points(self, cd_list):
        """real points as interval enclosures (mpmath.iv, outward rounding): [c, d, e, g] intervals."""
        out = []
        for (c0, d0) in cd_list:
            PA = self.bq(self.A_cd, c0, d0); PB = self.bq(self.B_cd, c0, d0)
            if PA == 0: continue
            qa = [PB * self.A_eg[k] - PA * self.B_eg[k] for k in range(3)]
            if qa[0] == 0: continue
            Dl = qa[1] ** 2 - 4 * qa[0] * qa[2]
            if Dl <= 0: continue
            sD = iv.sqrt(iv.mpf(Dl))
            for sg in (1, -1):
                rho = (-qa[1] + sg * sD) / (2 * qa[0])
                den = self.bq(self.A_eg, rho, 1)
                G = -PA / den
                if not (G.a > 0): continue
                gg = iv.sqrt(G)
                out.append([iv.mpf(c0), iv.mpf(d0), rho * gg, gg])
        return out

def ev_poly(P_, pt):
    """evaluate an integral sympy Poly in c,d,e,g (and W) at a point (PARI p-adics or iv intervals)."""
    tot = 0
    for mo, co in zip(P_.monoms(), P_.coeffs()):
        term = int(co)
        for x_, k_ in zip(pt, mo):
            if k_: term = term * x_ ** k_
        tot = tot + term
    return tot

def padic_class(x, p):
    """the class of a non-zero p-adic number in Q_p^*/Q_p^*2 as (v mod 2, unit class); the precision is asserted."""
    v = int(pari.valuation(x, p))
    rp = int(pari.padicprec(x, p)) - v      # relative precision
    u = x / pari('%d' % p) ** v
    need = 3 if p == 2 else 1
    assert rp >= need, 'p-adic precision lost (relative precision %d at p = %d)' % (rp, p)
    ui = int(pari.lift(u + pari('O(%d^%d)' % (p, need))))
    if p == 2: return (v % 2, ui % 8)
    return (v % 2, 0 if pow(ui, (p - 1) // 2, p) == 1 else 1)

def int_class(n, p):
    n = Fr(n)
    return padic_class(pari('%d/%d + O(%d^%d)' % (n.numerator, n.denominator, p, prec(p) + 10)), p)

def hilbert(x, y, p):
    """(x, y)_p for rationals / p-adic numbers given as classes (v, u)."""
    (a, ua), (b, ub) = x, y
    if p == 2:
        eps = lambda u: ((u - 1) // 2) % 2
        om = lambda u: ((u * u - 1) // 8) % 2
        s = (eps(ua) * eps(ub) + a * om(ub) + b * om(ua)) % 2
        return -1 if s else 1
    # odd p: (p^a u, p^b v)_p = (-1)^(a b eps(p)) (u|p)^b (v|p)^a
    epsp = ((p - 1) // 2) % 2
    s = (a * b * epsp + b * ua + a * ub) % 2
    return -1 if s else 1

# ---------------------------------------------------------------- the Sect. 7 data on C3' (rank_fisher5, re-run)
def level3_data():
    js = json.load(open(os.path.join(here, 'rank_fisher_curves.json')))
    cv = js[2]; side = "E'"
    assert cv['name'] == '72,65 E+' and cv['theta4_reduction'][0] == side
    S = cv['sides'][side]
    R = R5.analyse(cv['name'], side, S['level3_curve'], dd=int(S['dd']))
    assert R['pushout'] is not None and R['pushout3'] is not None and R['comp'] is not None
    return cv, R

def setup(R, mut=None):
    """rational data of Sect. 7: QA, QB integral; l_1 = al1 + B1 sqrt(D1), l_3 = al3 + B3 sqrt(D3); f, f_dag, q, c."""
    QA, QB = sp.expand(R['QA']), sp.expand(R['QB'])
    for Q_ in (QA, QB):
        assert all(sp.Rational(k).q == 1 for k in sp.Poly(Q_, *VS).coeffs())
    dU, dV = int(R['dU']), int(R['dV'])
    D1, D3 = sqfree(dU), sqfree(dV)
    m1, m3 = isqrt(dU // D1), isqrt(dV // D3)
    assert m1 * m1 * D1 == dU and m3 * m3 * D3 == dV and (D1, D3) == (-2, 2)
    al1, B1 = split_t(R['pushout']['l1'], D1, sp.Rational(1, m1))
    al3, B3 = split_t(R['pushout3']['l3'], D3, sp.Rational(1, m3))
    # t = m sqrt(D): coefficient of t times m is the coefficient of sqrt(D)... (split_t multiplies by its m argument)
    al1, B1 = al1, sp.expand(B1 * m1 * m1)    # B t = B m sqrt(D)
    al3, B3 = al3, sp.expand(B3 * m3 * m3)
    f = sp.expand(R['pushout']['f']); fd = sp.expand(R['pushout3']['f'])
    assert sp.expand(al1 ** 2 - D1 * B1 ** 2 - f) == 0, 'f != N(l_1)'
    assert sp.expand(al3 ** 2 - D3 * B3 ** 2 - fd) == 0, 'f_dag != N(l_3)'
    q = sp.expand(R['comp']['q']); cc = Fr(str(R['comp']['c']))
    if mut == 'q': q = sp.expand(q + e * g)
    assert in_ideal_C3(sp.expand(f * fd - sp.Rational(cc.numerator, cc.denominator) * q * q), QA, QB), \
        'f f_dag - c q^2 is not in (QA, QB)'
    s0, f0 = poly_int(f, VS)
    return dict(QA=QA, QB=QB, D1=D1, D3=D3, al1=al1, B1=B1, al3=al3, B3=B3, f=f, fd=fd, q=q, c=cc,
                s0=Fr(int(s0.p), int(s0.q)), f0=f0)

def check_smooth_primes(D, known=None):
    """the primes of the discriminant of det(x A + y B) (A, B: twice the symmetric matrices of QA, QB), factored
    completely (the known primes divided out first; the cofactor must be 1 or proven primes).  Returns (digits,
    primes)."""
    A2 = sp.hessian(D['QA'], VS); B2 = sp.hessian(D['QB'], VS)
    x_, y_ = sp.symbols('x_ y_')
    F = sp.Poly(sp.expand((x_ * A2 + y_ * B2).det()), x_, y_)
    co = [int(F.coeff_monomial(x_ ** (4 - i) * y_ ** i)) for i in range(5)]
    Fp = pari('Pol([' + ','.join(str(k) for k in co) + '])')
    disc = int(pari.poldisc(Fp))
    assert disc != 0 and co[0] != 0
    n = abs(disc); ps = []
    for p in (known if known is not None else B0):
        if n % p == 0: ps.append(p)
        while n % p == 0: n //= p
    if n > 1:
        Fa = pari.factor(n, 10 ** 7)
        for i in range(len(Fa[0])):
            q_ = int(Fa[0][i])
            assert pari.isprime(q_), 'incomplete factorisation of the pencil discriminant'
            ps.append(q_)
    return len(str(abs(disc))), sorted(set(ps))


S_ELS = [2, 3, 5, 7, 11, 13, 97, 137, 710053, 79553041, 97505521]

def cd_samples(p, n, seed):
    rng = random.Random(seed * 1000003 + p)
    out = [(1, 0), (0, 1)]
    for _ in range(n // 2):
        k = rng.randrange(0, 16)
        out.append((rng.randrange(-10 ** k, 10 ** k + 1), rng.randrange(1, 100)))
    out += [(rng.randrange(0, p ** 4) * rng.choice((1, p, p * p)), rng.randrange(1, p ** 4)) for _ in range(n // 2)]
    return [(x, y) for (x, y) in out if gcd(x, y) == 1]

def c3_points(M, C, p, n=100, seed=1):
    return M.points(p, cd_from_Dlam(C, p, n, seed) + cd_samples(p, n, seed))

def f0_classes(M, D, C, p, n=100, seed=1):
    """classes of f0 at sampled points of C3'(Q_p): {class: point}."""
    seen = {}
    for P in c3_points(M, C, p, n, seed):
        v = ev_poly(D['f0'], P)
        if v == 0: continue
        try: k = padic_class(v, p)
        except AssertionError: continue
        seen.setdefault(k, P)
    return seen

def real_cds(n=400, seed=7):
    rng = random.Random(seed)
    cds = []
    for _ in range(n):
        k = rng.randrange(0, 16)
        cds.append((rng.choice((1, -1)) * rng.randrange(1, 10 ** k + 1), rng.randrange(1, 100)))
    return cds + [(1, 0), (0, 1)]

def real_cds_from_Dlam(C, n=40):
    """integer (c, d) close to the projection of rank_fisher's real point of D_lam (and perturbations): a sampler
    only (Model.real_points recomputes the point and every sign is an interval enclosure)."""
    lp = C.get('local_points', {}).get('inf')
    if not lp: return []
    import mpmath
    mpmath.mp.dps = 60
    lam = int(C['lam']); U = [int(x) for x in C['U']]; pX = [int(x) for x in C['p_X']]
    r0, s0_ = Fr(lp[0]), Fr(lp[1]); sx = int(lp[2])
    out = []
    rng = random.Random(5)
    for i in range(n):
        r = mpmath.mpf(r0.numerator) / r0.denominator * (1 + (mpmath.mpf(rng.random()) - 0.5) * 10 ** -rng.randrange(1, 12)) if i else mpmath.mpf(r0.numerator) / r0.denominator
        s_ = mpmath.mpf(s0_.numerator) / s0_.denominator
        v = lam * (U[0] * s_ * s_ + U[1] * r * s_ + U[2] * r * r)
        if v <= 0: continue
        X = sx * mpmath.sqrt(v)
        cc = pX[2] * r - X * pX[0]; dd = pX[2] * s_ - X * pX[1]
        if dd == 0: continue
        k = 10 ** 40
        out.append((int(mpmath.nint(cc / dd * k)), k))
    return out

def real_signs(M, D, n=400, C=None):
    signs = {}
    for P in M.real_points(real_cds(n) + (real_cds_from_Dlam(C) if C else [])):
        v = ev_poly(D['f0'], P)
        if v.a > 0: signs.setdefault(1, P)
        elif v.b < 0: signs.setdefault(-1, P)
    return signs

def cls_vec(k, p):
    v, u = k
    if p == 2: return (v, ((u - 1) // 2) % 2, ((u * u - 1) // 8) % 2)
    return (v, u)

def vec_cls(x, p):
    if p == 2:
        u = next(u_ for u_ in (1, 3, 5, 7) if (((u_ - 1) // 2) % 2, ((u_ * u_ - 1) // 8) % 2) == (x[1], x[2]))
        return (x[0], u)
    return (x[0], x[1])

def coset_closure(Kp, p):
    """the coset k0 + span(k - k0) of F_2-vectors spanned by the sampled classes."""
    ks = [cls_vec(k, p) for k in Kp]
    if not ks: return dict(Kp)
    k0 = ks[0]
    span = {tuple(0 for _ in k0)}
    for k in ks:
        dlt = tuple((x + y) % 2 for x, y in zip(k, k0))
        span |= {tuple((x + y) % 2 for x, y in zip(s_, dlt)) for s_ in span}
    out = dict(Kp)
    for s_ in span:
        out.setdefault(vec_cls(tuple((x + y) % 2 for x, y in zip(s_, k0)), p), None)
    return out

def local_image_classes(R, p):
    """Im(delta_phihat, p) (E side) as classes (v, u) of this script, from rank_indep's certified local images."""
    poss, cert = R['locS'][p]
    assert poss == cert, 'the local image at %d is not decided' % p
    import rank_indep as RI_
    return [int_class(dl, p) for dl in RI_.class_reps(p) if RI_.qp_class(dl, p) in cert]

def choose_Xi(M, D, C, say=say, places=None, images=None):
    """Xi in <-1, places> with a point of C4: Xi W^2 = f0 at every place of places and at R (sampled points).
    images: optional {p: Im(delta_phihat, p) as classes}; then the sampled classes are closed up to the coset
    k0 Im(delta_phihat, p) (Fisher Sect. 3), only to choose Xi (the points are found and checked afterwards)."""
    S_ELS = places if places is not None else globals()['S_ELS']
    K = {p: f0_classes(M, D, C, p) for p in S_ELS}
    if images:
        for p in S_ELS:
            if p in images and K[p]:
                K[p] = coset_closure(dict(K, **{})[p] | {vec_cls(tuple((x + y) % 2 for x, y in zip(cls_vec(next(iter(K[p])), p), cls_vec(im, p))), p): None for im in images[p]}, p)
    Rs = real_signs(M, D, C=C)
    for p in S_ELS:
        say('  f0 classes found on C3\'(Q_%d): %s' % (p, sorted(K[p])))
    # the classes of f0 on C3'(Q_p) form a coset of Im(delta_phihat, p) (Fisher Sect. 3): close the sampled set up to
    # the coset it spans (only to choose Xi; a point of C4 at every place is then found and checked in local_terms)
    for p in S_ELS:
        K[p] = coset_closure(K[p], p)
    say('  signs of f0 on C3\'(R): %s' % sorted(Rs))
    gens = [-1] + S_ELS
    sols = []
    for mask in range(1 << len(gens)):
        Xi = 1
        for i, q_ in enumerate(gens):
            if mask >> i & 1: Xi *= q_
        if (1 if Xi > 0 else -1) not in Rs: continue
        if all(int_class(Xi, p) in K[p] for p in S_ELS):
            sols.append(Xi)
    return sols, K, Rs

def cd_from_Dlam(C, p, n, seed):
    """(c : d) values (lifted to integers modulo p^N) of points of D_lam(Q_p): X^2 = lam U(r, s), Y^2 = lam V(r, s),
    projected from p_X onto X = 0 (rank_fisher5 (1)-(2): (r : s : X) = (p1 : p2 : p3)(c, d)).  Only a sampler: the
    points of C3' are then computed from (c, d) by Model.points (and checked there)."""
    lam = int(C['lam']); U = [int(x) for x in C['U']]; pX = [int(x) for x in C['p_X']]
    N = prec(p); rng = random.Random(seed * 7919 + p)
    base = []
    lp = C['local_points'].get(str(p))
    if lp:
        r0, s0_ = Fr(lp[0]), Fr(lp[1]); D_ = r0.denominator * s0_.denominator
        base.append((int(r0 * D_), int(s0_ * D_)))
    out = []
    for i in range(n):
        if base and i % 2 == 0:
            k = rng.randrange(1, 12)
            r, s = base[0][0] + rng.randrange(-p ** 2, p ** 2) * p ** k, base[0][1] + rng.randrange(-p, p) * p ** k
        else:
            k = rng.randrange(0, 8)
            r, s = rng.randrange(-p ** k, p ** k + 1), rng.randrange(1, p ** 3)
        Xv = pari('%d + O(%d^%d)' % (lam * R5.qv(U, r, s), p, N))
        if Xv == 0 or not pari.issquare(Xv): continue
        X = pari.sqrt(Xv)
        for sg in (1, -1):
            cc = pX[2] * r - sg * X * pX[0]; dd = pX[2] * s - sg * X * pX[1]
            if cc == 0 or dd == 0: continue
            out.append((cc, dd))
    return out

# ---------------------------------------------------------------- the pushout g0 on C4 and its checks
def norm_solution(D_, n):
    """kappa = (X + Y sqrt(D_)) / Z with N(kappa) = X^2 - D_ Y^2 / Z^2 = n (n a squarefree integer): PARI qfsolve
    (a finder; the norm is re-checked exactly)."""
    sol = pari.qfsolve(pari.matdiagonal([1, -D_, -n]))
    assert str(pari.type(sol)) == 't_COL', 'no solution of X^2 - %d Y^2 = %d Z^2' % (D_, n)
    X, Y, Z = [int(v) for v in sol]
    assert Z != 0 and X * X - D_ * Y * Y == n * Z * Z
    return Fr(X, Z), Fr(Y, Z)

def build_g0(D, Xi, kap1=None, kap3=None, mut=None):
    """N1, N1', R3, R3', g0 = R3 / (W N1) on C4: Xi W^2 = f0 (xi4 = Xi s0, xi4 W^2 = f).  Exact rational data."""
    xi4 = Fr(Xi) * D['s0']
    D1, D3 = D['D1'], D['D3']
    if kap1 is None: kap1 = norm_solution(D1, sqfree(xi4))
    if kap3 is None: kap3 = norm_solution(D3, sqfree(D['c'] * xi4))
    ka, kb = kap1; la, lb = kap3
    Nk1 = ka * ka - D1 * kb * kb; Nk3 = la * la - D3 * lb * lb
    sig1 = rat_sqrt(xi4 * Nk1)                    # sigma_1^2 = xi4 N(k_1)
    tau3 = rat_sqrt(D['c'] / (xi4 * Nk3))         # tau_3^2 = c / (xi4 N(k_3))
    R_ = lambda x: sp.Rational(x.numerator, x.denominator)
    if mut == 'sigma': sig1 = 2 * sig1
    if mut == 'tau': tau3 = 3 * tau3
    # l / k = l conj(k) / N(k) = (Tr + delta sqrt(D)) / 2
    Tr1 = sp.expand(2 * (D['al1'] * R_(ka) - D1 * D['B1'] * R_(kb)) / R_(Nk1))
    de1 = sp.expand(2 * (D['B1'] * R_(ka) - D['al1'] * R_(kb)) / R_(Nk1))
    Tr3 = sp.expand(2 * (D['al3'] * R_(la) - D3 * D['B3'] * R_(lb)) / R_(Nk3))
    de3 = sp.expand(2 * (D['B3'] * R_(la) - D['al3'] * R_(lb)) / R_(Nk3))
    N1 = sp.expand(Tr1 + 2 * R_(xi4) * W / R_(sig1)); N1p = sp.expand(Tr1 - 2 * R_(xi4) * W / R_(sig1))
    R3 = sp.expand(W * Tr3 + 2 * R_(tau3) * D['q']); R3p = sp.expand(W * Tr3 - 2 * R_(tau3) * D['q'])
    return dict(Xi=Xi, xi4=xi4, kap1=kap1, kap3=kap3, Nk1=Nk1, Nk3=Nk3, sig1=sig1, tau3=tau3, Tr1=Tr1, de1=de1,
                Tr3=Tr3, de3=de3, N1=N1, N1p=N1p, R3=R3, R3p=R3p)

def check_g0(D, G, say=say):
    """(I1), (I2) of the docstring, exactly.  Returns the integers whose primes need care (none expected)."""
    xi4 = sp.Rational(G['xi4'].numerator, G['xi4'].denominator)
    f = D['f']
    # (I1) N1 N1' - D1 delta1^2 and R3 R3' - D3 W^2 delta3^2 vanish modulo (QA, QB, xi4 W^2 - f)
    for (lab, P_, Q_) in (('N1 N1\' = -2 delta1^2', sp.expand(G['N1'] * G['N1p']), sp.expand(D['D1'] * G['de1'] ** 2)),
                          ('R3 R3\' = 2 W^2 delta3^2', sp.expand(G['R3'] * G['R3p']),
                           sp.expand(D['D3'] * W ** 2 * G['de3'] ** 2))):
        A_, B_ = reduce_W(sp.expand(P_ - Q_), xi4, f)
        assert B_ == 0 or in_ideal_C3(B_, D['QA'], D['QB']), 'the W-part of %s is not in the ideal' % lab
        assert in_ideal_C3(A_, D['QA'], D['QB']), 'the identity %s fails modulo I(C4)' % lab
        say('  (I1) %s modulo I(C4) = (QA, QB, xi4 W^2 - f): checked (exact linear algebra)' % lab)
    # (I2) no common zeros
    r1 = line_resultant(D['al1'], D['B1'], D['QA'], D['QB'])
    assert r1 != 0, 'l_1 and conj(l_1) have a common zero on C3\''
    r3 = line_resultant(D['al3'], D['B3'], D['QA'], D['QB'])
    assert r3 != 0, 'l_3 and conj(l_3) have a common zero on C3\''
    # delta3 on A_1: the line {l_1 = 0, delta3 = 0} over L = Q(t), t^2 = D1
    red = lambda z: sp.expand(sp.Poly(sp.expand(z), t).rem(sp.Poly(t ** 2 - D['D1'], t)).as_expr())
    l1t = sp.expand(D['al1'] + D['B1'] * t)
    r13 = line_resultant(l1t, G['de3'], D['QA'], D['QB'], red=red)
    assert r13 != 0, 'delta3 vanishes at a point of A_1'
    say('  (I2) resultants on the lines {l_1 = conj l_1 = 0}, {l_3 = conj l_3 = 0} (over Q) and {l_1 = delta3 = 0} '
        '(over L) are non-zero (%d, %d, %d digits): A_1, A_2 disjoint, A_3, A_4 disjoint, delta3 != 0 on A_1, A_2'
        % (len(str(sp.Rational(r1).p)), len(str(sp.Rational(r3).p)), len(str(sp.Poly(r13, t).coeffs()[0]))))
    # (I3') A_1 consists of two distinct points: the cone Q_w (w0 a root of det(w A + B) in L) restricted to the plane
    # l_1 = 0 is kappa m^2 (rank 1: the ruling {l_1 = m = 0}); QA on that line is a binary quadratic of non-zero
    # discriminant.  (A_2 is the conjugate.)
    nd = ruling_discriminant(D, l1t, red)
    say('  (I3\') the ruling of the cone through P1 meets C3\' in two distinct points (discriminant of QA on it != 0, %d digits)' % nd)
    return True

def ruling_discriminant(D, l1t, red):
    A2 = sp.hessian(D['QA'], VS) / 2; B2 = sp.hessian(D['QB'], VS) / 2
    w_ = sp.Symbol('w_')
    F = sp.Poly(sp.expand((w_ * A2 + B2).det()), w_)
    facs = sp.factor_list(F.as_expr(), w_)[1]
    cand = []
    for fa, _ in facs:
        P_ = sp.Poly(fa, w_)
        if P_.degree() != 2: continue
        a2, a1, a0 = [sp.Rational(x) for x in P_.all_coeffs()]
        disc = a1 * a1 - 4 * a2 * a0
        q_ = disc / D['D1']
        if is_rat_square(Fr(int(q_.p), int(q_.q))):
            k = rat_sqrt(Fr(int(q_.p), int(q_.q)))
            for sg in (1, -1):
                cand.append((-a1 + sg * sp.Rational(k.numerator, k.denominator) * t) / (2 * a2))
    assert cand, 'no pair of roots over L'
    for w0 in cand:
        Qw = sp.expand(w0 * D['QA'] + D['QB'])
        co = {v_: red(sp.Poly(l1t, *VS).coeff_monomial(v_)) for v_ in VS}
        piv = next(v_ for v_ in VS if co[v_] != 0)
        others = [v_ for v_ in VS if v_ != piv]
        sub = sp.expand(-sum(co[v_] * v_ for v_ in others))
        Qr = red(sp.expand(Qw.subs({**{v_: co[piv] * v_ for v_ in others}, piv: sub}, simultaneous=True)))
        M3 = sp.Matrix(3, 3, lambda i, j: red(sp.diff(Qr, others[i], others[j]) / 2))
        if all(red(M3[i, k] * M3[j, l_] - M3[i, l_] * M3[j, k]) == 0 for i in range(3) for j in range(3)
               for k in range(3) for l_ in range(3)):
            # rank <= 1: Qr = (row . v)^2 / M3[i,i] for a row with M3[i,i] != 0; the ruling is {l_1 = 0, row . v = 0}
            i0 = next(i for i in range(3) if M3[i, i] != 0)
            mlin = sp.expand(sum(M3[i0, j] * others[j] for j in range(3)))      # in the coordinates 'others'
            # back to c,d,e,g: on the plane, v_piv = sub / co[piv] after the scaling v_other -> co[piv] v_other
            # parametrise the plane by 'others' (scaled) and the line by m = 0
            s_, u_ = sp.symbols('s_ u_')
            # solve mlin = 0 for one of the others
            jpiv = next(j for j in range(3) if red(M3[i0, j]) != 0)
            rest = [others[j] for j in range(3) if j != jpiv]
            val = sp.expand(-sum(M3[i0, j] * others[j] for j in range(3) if j != jpiv))
            pt_others = {rest[0]: s_ * M3[i0, jpiv], rest[1]: u_ * M3[i0, jpiv]}
            pt_others[others[jpiv]] = red(val.subs({rest[0]: s_, rest[1]: u_}, simultaneous=True))
            assert red(sp.expand(mlin.subs(pt_others, simultaneous=True))) == 0
            full = {v_: co[piv] * pt_others[v_] for v_ in others}
            full[piv] = red(sub.subs(pt_others, simultaneous=True))
            assert red(sp.expand(l1t.subs(full, simultaneous=True))) == 0, 'the line is not in the plane l_1 = 0'
            assert red(sp.expand(Qw.subs(full, simultaneous=True))) == 0, 'the line is not on the cone'
            qa = sp.Poly(red(sp.expand(D['QA'].subs(full, simultaneous=True))), s_, u_)
            a_ = [red(qa.coeff_monomial(m_)) for m_ in (s_ ** 2, s_ * u_, u_ ** 2)]
            dsc = red(a_[1] ** 2 - 4 * a_[0] * a_[2])
            assert any(x != 0 for x in a_), 'QA vanishes on the ruling'
            assert dsc != 0, 'the ruling is tangent to C3\': A_1 is a double point'
            return len(str(max(abs(sp.Rational(x).p) for x in sp.Poly(dsc, t).coeffs())))
    raise AssertionError('l_1 = 0 is not tangent to a cone of the U pair')

def c4_points(M, D, C, Xi, p, n=60, seed=1):
    """points (c, d, e, g, W) of C4(Q_p): Xi W^2 = f0(c, d, e, g) on C3'(Q_p) (PARI p-adics)."""
    out = []
    for P in c3_points(M, C, p, n, seed):
        v = ev_poly(D['f0'], P)
        if v == 0 or int(pari.padicprec(v, p)) - int(pari.valuation(v, p)) < 6: continue
        w2 = v / Xi
        if not pari.issquare(w2): continue
        out.append(P + [pari.sqrt(w2)])
    return out

def normalise(G):
    """primitive integral N1^0, R3^0 (polys in c, d, e, g, W) and the ratio r of R3 = t3 (W T3 + r Q0)."""
    vs5 = VS + (W,)
    n1, N10 = poly_int(G['N1'], vs5)
    r3, R30 = poly_int(G['R3'], vs5)
    t3, T30 = poly_int(G['Tr3'], VS)
    s_, Q0 = poly_int(sp.expand(G['R3'] - W * G['Tr3']), VS)
    r = Fr(int(sp.Rational(s_ / t3).p), int(sp.Rational(s_ / t3).q))
    return dict(N10=N10, R30=R30, T30=T30, Q0=Q0, r=r, n1=n1, r3=r3)

def g0_value(Nz, pt):
    """g0^0 = R3^0 / (W N1^0) at a point (c, d, e, g, W)."""
    num = ev_poly(Nz['R30'], pt); den = pt[4] * ev_poly(Nz['N10'], pt)
    assert num != 0 and den != 0, 'the point is a zero or pole of g0'
    return num / den

def degenerate_primes(D, Nz, say=say):
    """for p | den(r): R3^0 reduces to a unit times Q0; it must not lie in span(QA, QB) modulo p."""
    den = Nz['r'].denominator
    F = pari.factor(den)
    ps = [int(F[0][i]) for i in range(len(F[0]))]
    bad = []
    mons = [VS[i] * VS[j] for i in range(4) for j in range(i, 4)]
    rows = [[int(sp.Poly(Q_, *VS).coeff_monomial(m_)) for m_ in mons] for Q_ in (D['QA'], D['QB'], Nz['Q0'].as_expr())]
    for p in ps:
        rk = int(pari.matrank(pari.Mod(pari.matrix(3, 10, [x for r_ in rows for x in r_]), p)))
        if rk < 3: bad.append(p)
    say('  R3 = t3 (W T3 + r Q0) with r = %s: primes of den(r): %s; Q0 mod p in span(QA, QB) for: %s'
        % ('%d digits / %d digits' % (len(str(abs(Nz['r'].numerator))), len(str(Nz['r'].denominator))), ps, bad))
    return ps, bad

# ---------------------------------------------------------------- the first Selmer group of the E side
def in_Sprime1(y, a, b, p, tries=400):
    """y in Im(delta_phihat, p) = x(E(Q_p)) mod squares (E: y^2 = x(x^2 + a x + b)): a point x = y u^2 found
    (x (x^2 + a x + b) a non-zero square in Q_p), or x = 0 for the class of b (the image of T)."""
    if p == 'inf':
        for k in range(-60, 61):
            x = Fr(y) * Fr(2) ** k
            if x * (x * x + a * x + b) > 0: return True
        return False
    rng = random.Random(p)
    N = prec(p)
    for i in range(tries):
        u = Fr(rng.randrange(1, p ** 4), 1) * Fr(p) ** rng.randrange(-8, 9)
        x = Fr(y) * u * u
        v = x * (x * x + a * x + b)
        if v == 0: continue
        vv = pari('%d/%d + O(%d^%d)' % (v.numerator, v.denominator, p, N))
        if pari.issquare(vv): return True
    return False

def local_terms(M, D, C, Nz, Xi, places, ys, npts=40, seeds=(1, 2), say=say):
    """for each place p: points of C4(Q_p), the classes of g0^0 there and the terms (g0^0(P), y)_p for each y in ys,
    which must not depend on the point.  Returns {p: {y: term}}, {p: classes}."""
    terms, classes = {}, {}
    for p in places:
        cls = []
        sds = list(seeds) + [7, 8, 9, 10, 11, 12]
        k = 0
        while k < len(sds) and (k < len(seeds) or len(cls) < 4):
            for P in c4_points(M, D, C, Xi, p, npts if k < len(seeds) else 150, sds[k]):
                try: cls.append(padic_class(g0_value(Nz, P), p))
                except AssertionError: continue
            k += 1
        assert cls, 'no usable point of C4(Q_%d) found' % p
        T_ = {}
        for y in ys:
            vals = {hilbert(k, int_class(y, p), p) for k in cls}
            assert len(vals) == 1, 'the term (g0(P), %d)_%d depends on the point P: %s' % (y, p, vals)
            T_[y] = vals.pop()
        terms[p] = T_; classes[p] = sorted(set(cls))
    return terms, classes

def real_term(M, D, Nz, Xi, ys, C=None):
    """C4(R): a real point (interval enclosure) with Xi f0 > 0, and the sign of g0^0 there; (g0, y)_R = -1 iff
    both are negative."""
    found = False
    for P in M.real_points(real_cds(400, 11) + (real_cds_from_Dlam(C) if C else [])):
        v = ev_poly(D['f0'], P) / Xi
        if not v.a > 0: continue
        found = True
        Wv = iv.sqrt(v)
        pt = P + [Wv]
        num = ev_poly(Nz['R30'], pt); den = Wv * ev_poly(Nz['N10'], pt)
        if (num.a > 0 or num.b < 0) and (den.a > 0 or den.b < 0):
            sg = (1 if num.a > 0 else -1) * (1 if den.a > 0 else -1)
            return {y: (-1 if (sg < 0 and y < 0) else 1) for y in ys}, sg
    assert found, 'no real point of C4 found'
    assert all(y > 0 for y in ys), 'the sign of g0 at a real point is not certified'
    return {y: 1 for y in ys}, None

# ---------------------------------------------------------------- the computation
def primes_of(n):
    n = abs(int(n))
    if n <= 1: return []
    F = pari.factor(n)
    return [int(F[0][i]) for i in range(len(F[0]))]

def theta4_for(M, D, C, Xi, ys, say=say, label='', B0=B0, kmul=None):
    """the whole computation for one level-4 lift C4: Xi W^2 = f0.  Returns {y: product of the local terms}."""
    G = build_g0(D, Xi)
    if kmul:
        # (K5) other solutions of the two norm equations: k_1 (1 + sqrt D1)^2 and k_3 (1 + sqrt D3)^2 (norms times squares)
        mul = lambda k, Dq: (k[0] * (1 + Dq) + Dq * k[1] * 2, k[0] * 2 + k[1] * (1 + Dq))
        G = build_g0(D, Xi, kap1=mul(G['kap1'], D['D1']), kap3=mul(G['kap3'], D['D3']))
    say('%sC4: Xi W^2 = f0 with Xi = %d (xi4 = %s); k_1 = %s + %s sqrt(-2) (N = %s), k_3 = %s + %s sqrt(2) (N = %s); '
        'sigma_1 = %s, tau_3 = %s/%s (%d/%d digits)'
        % (label, Xi, G['xi4'], G['kap1'][0], G['kap1'][1], G['Nk1'], G['kap3'][0], G['kap3'][1], G['Nk3'], G['sig1'],
           '...' , '...', len(str(G['tau3'].numerator)), len(str(G['tau3'].denominator))))
    check_g0(D, G, say=say)
    Nz = normalise(G)
    ps, bad = degenerate_primes(D, Nz, say=say)
    small = [p for p in range(2, 47) if pari.isprime(p)]
    extra = sorted(set(q_ for y in ys for q_ in primes_of(y)))
    S = sorted(set(small + B0 + primes_of(Xi) + bad + extra))
    say('  S = {R} u %s  (primes < 47, B0, primes of Xi, the primes where R3^0 reduces into (QA, QB), primes of the y)' % S)
    terms, classes = local_terms(M, D, C, Nz, Xi, S, ys)
    rt, sg = real_term(M, D, Nz, Xi, ys, C=C)
    say('  R: a point of C4(R) (interval enclosure, Xi f0 > 0), g0^0 there %s' % ('> 0' if sg == 1 else ('< 0' if sg == -1 else 'of uncertified sign')))
    # (I4): two points of C4(Q_p) with g0(P) g0(P') not a square
    p4 = next((p for p in S if len(classes[p]) > 1), None)
    assert p4 is not None, '(I4) not shown: g0 might be a constant times a square'
    say('  (I4) at p = %d g0^0 takes the classes %s on C4(Q_p): g0 is not a constant times a square, so [B] = T\''
        % (p4, classes[p4]))
    out = {}
    for y in ys:
        nontriv = [p for p in S if terms[p][y] == -1] + (['R'] if rt[y] == -1 else [])
        prod = rt[y]
        for p in S: prod *= terms[p][y]
        out[y] = prod
        say('  sum_v (g0(P_v), %d)_v: the term is -1 exactly at %s; product = %+d' % (y, nontriv if nontriv else 'no place', prod))
    for p in S:
        say('    p = %d (%d mod 8): classes of g0^0 on C4(Q_p) %s; terms %s'
            % (p, p % 8, classes[p], ' '.join('%d:%+d' % (y, terms[p][y]) for y in ys)))
    return out, S

def main():
    say('rank_fisher_theta4.py: Fisher\'s Theta_4(y0, y1) on E+ of (72,65), exact (arXiv:1509.03234 Sects. 2, 3, 7)')
    a, b, Gpt = R2.sigma_curve(72, 65, 'E+')
    say('E: y^2 = x(x^2 + a x + b), a = %d, b = %d; known point G with x(G) = %s' % (a, b, Gpt[0]))
    cv, R = level3_data()
    C = cv['sides']["E'"]['level3_curve']
    say('C3\' (level-3 lift of y0 = 2, E\' side, lam = %s, xi = %s) re-built and re-checked by rank_fisher5.analyse '
        '(steps (1)-(8))' % (C['lam'], C['xi']))
    D = setup(R)
    say('  l_1 = al1 + B1 sqrt(-2), l_3 = al3 + B3 sqrt(2); f = N(l_1), f_dag = N(l_3) (identities), '
        'f f_dag - c q^2 in (QA, QB) (exact linear algebra); f = s0 f0, s0 = %s' % D['s0'])
    nd, ps_ = check_smooth_primes(D)
    assert set(ps_) <= set(B0), 'the pencil discriminant has a prime outside B0'
    say('  the discriminant of det(x A + y B) (%d digits) is a product of primes of B0 = %s: C3\' is smooth modulo '
        'every other prime' % (nd, B0))
    M = Model(D['QA'], D['QB'])
    yb, yG = sqfree(b), sqfree(Gpt[0])
    y1 = 2
    say('y1 = %d, and the images of E(Q) in S\'_1: b = %d (from T) and x(G) = %d mod squares' % (y1, yb, yG))
    # y1 in S'_1: a point of E(Q_v) with x = y1 u^2 at the places of S (and the unit classes elsewhere)
    Sy = ['inf'] + [p for p in range(2, 47) if pari.isprime(p)] + B0[1:]
    Sy = sorted(set(Sy), key=lambda v: -1 if v == 'inf' else v)
    RD = RI.first_descent(a, b)
    for p in RD['bad']:
        assert R2.certified(y1, a, b, p, RD['locS']), 'y1 is not certified in Im(delta_phihat) at %s' % p
    assert RI.real_ok(y1, a, b)
    say('  y1 = 2 is in Im(delta_phihat, v) at R and at the bad primes %s (rank_indep\'s certified local images: the '
        'quartic C_2 is soluble there); elsewhere 2 is a unit at a prime of good reduction: y1 in S\'_1 = S\'_4' % RD['bad'])
    say('The level-4 Selmer set: Xi with a point of Xi W^2 = f0 at R and at every p of %s:' % S_ELS)
    sols, K, Rs = choose_Xi(M, D, C)
    sols = sorted(sols, key=abs)
    say('  %d solutions: %s (a coset of <2, 79553041, 97505521> = S\'_1, as Fisher Sect. 3 predicts)' % (len(sols), sols))
    assert len(sols) == 8
    ys = [y1, yb, yG]
    res, S = theta4_for(M, D, C, sols[0], ys)
    assert res[yb] == 1 and res[yG] == 1, 'CONTROL (K2) FAILED: a class of E(Q) pairs non-trivially'
    say('  (K2) the classes of E(Q) pair to 0, as they must')
    th = res[y1]
    others = {}
    for Xi2 in sols[1:3]:
        r2, _ = theta4_for(M, D, C, Xi2, ys, label='(K3) ')
        assert r2[yb] == 1 and r2[yG] == 1, 'CONTROL (K2) FAILED for the twist %d' % Xi2
        assert r2[y1] == th, 'CONTROL (K3) FAILED: the value depends on the level-4 lift'
        others[Xi2] = r2[y1]
    say('  (K3) the same value for the level-4 lifts Xi = %s' % sorted(others))
    r5, _ = theta4_for(M, D, C, sols[0], ys, label='(K5) other k_1, k_3: ', kmul=True)
    assert r5[y1] == th and r5[yb] == 1 and r5[yG] == 1, 'CONTROL (K5) FAILED: the value depends on k_1, k_3'
    say('  (K5) the same value with k_1 (1 + sqrt -2)^2, k_3 (1 + sqrt 2)^2 in place of k_1, k_3 (another g0)')
    say('')
    say('RESULT: Theta_4(y0, y1) = sum_v (g0(P_v), 2)_v = %s  (%s)' % ('1/2 (non-zero)' if th == -1 else '0',
        'the product of the local Hilbert symbols is %+d' % th))
    if th == -1:
        say('  => S_5 = <images of E(Q)> (dim 1), S\'_5 = <images of E(Q)> (dim 2) and Lemma 2.2 (n = 3) gives '
            'rank E(Q) <= 1 + 2 - 2 = 1; with G of infinite order, rank E(Q) = 1 for E+ of (72,65), independently of Magma.')
    else:
        say('  => Theta_4 = 0: the bound stays 3 at this level.')
    return th

def negative():
    say('rank_fisher_theta4.py --negative: altered inputs that must be refused')
    cv, R = level3_data()
    C = cv['sides']["E'"]['level3_curve']
    a, b, Gpt = R2.sigma_curve(72, 65, 'E+')
    n = ok = 0
    def refuse(label, fn):
        nonlocal n, ok
        n += 1
        try: fn()
        except AssertionError as ex:
            ok += 1; say('  refused: %s (%s)' % (label, str(ex)[:100])); return
        say('  NOT REFUSED: %s' % label)
    D = setup(R)
    M = Model(D['QA'], D['QB'])
    quiet = lambda s: None
    refuse('q replaced by q + e g (f f_dag - c q^2 leaves the ideal)', lambda: setup(R, mut='q'))
    for m_, what in (('sigma', 'sigma_1 doubled in N1 (a mutated g0)'), ('tau', 'tau_3 tripled in R3 (a mutated g0)')):
        refuse(what, lambda m_=m_: check_g0(D, build_g0(D, 137, mut=m_), say=quiet))
    refuse('k_1 = 1 + sqrt(-2) (norm 3, not xi4 modulo squares): sigma_1 not rational',
           lambda: build_g0(D, 137, kap1=(Fr(1), Fr(1))))
    # a forged twist: Xi = 137 * 3 is not in the level-4 Selmer set (no point of C4 at some p)
    def forged_Xi():
        K = {p: f0_classes(M, D, C, p) for p in S_ELS}
        for p in S_ELS:
            assert int_class(137 * 3, p) in K[p], 'no point of 411 W^2 = f0 over Q_%d' % p
    refuse('a forged twist Xi = 137 * 3 (outside the Selmer coset)', forged_Xi)
    RD = RI.first_descent(a, b)
    def wrong_y():
        for p in RD['bad']:
            assert R2.certified(3, a, b, p, RD['locS']), '3 is not certified in Im(delta_phihat) at %s' % p
    refuse('a wrong y1 = 3 (not in S\'_1): the quartic C_3 is not soluble at some bad prime', wrong_y)
    # a wrong y1 = 3 fed to the pairing anyway: the local term must depend on the point somewhere
    G = build_g0(D, 137); Nz = normalise(G)
    refuse('the pairing evaluated with y = 3 (not Selmer): the local term depends on the point',
           lambda: local_terms(M, D, C, Nz, 137, [3, 5, 7, 13, 137], [3], say=quiet))
    # a function that is not a pushout: g0 * c / e (odd divisor): the terms must depend on the point somewhere
    Nbad = dict(Nz, R30=sp.Poly(sp.expand(Nz['R30'].as_expr() * c), *(VS + (W,))),
                N10=sp.Poly(sp.expand(Nz['N10'].as_expr() * e), *(VS + (W,))))
    refuse('g0 replaced by g0 c / e (not a pushout form): the local term depends on the point',
           lambda: local_terms(M, D, C, Nbad, 137, [2, 3, 5, 13], [2], say=quiet))
    # too little p-adic precision: every class must be refused (relative precision asserted), so no point is usable
    def low_prec():
        g_ = globals(); old = g_['prec']
        g_['prec'] = lambda p: 4
        try: local_terms(M, D, C, Nz, 137, [2], [2], say=quiet)
        finally: g_['prec'] = old
    refuse('the p-adic precision lowered to 4 digits at p = 2 (classes not certifiable)', low_prec)
    say('negative tests: %d of %d refused' % (ok, n))
    assert ok == n
    return ok, n

if __name__ == '__main__':
    if '--negative' in sys.argv:
        negative()
        if WRITE: open(os.path.join(here, 'rank_fisher_theta4_negative.txt'), 'w').write('\n'.join(OUT) + '\n')
    else:
        main()
        if WRITE: open(os.path.join(here, 'rank_fisher_theta4.txt'), 'w').write('\n'.join(OUT) + '\n')
