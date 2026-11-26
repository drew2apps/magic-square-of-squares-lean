#!/usr/bin/env python3
"""Fisher's fifth descent (Theta_4; T. Fisher, "Higher descents on an elliptic curve with a rational 2-torsion
point", arXiv:1509.03234, Section 7), FIRST STEPS ONLY, in exact arithmetic, on the level-3 curves C3 of
rank_fisher.py (rank_fisher_curves.json).

    python3 rank_fisher5.py              the ten level-3 curves (both sides of the five Magma-alone curves)
    python3 rank_fisher5.py --negative   altered inputs that must be refused
    (add --write to write rank_fisher5.txt / rank_fisher5_negative.txt; without it the script only prints)

Written by an AI model (Claude); not read by a person.  Experiment, not counted.  It does NOT compute Theta_4.

WHAT IS COMPUTED, EXACTLY (integers and rationals; polynomial identities expanded in full with sympy):
 (1) Both conics of D_lam (X^2 = lam U(r, s) and Y^2 = lam V(r, s), U = [u0, u1, u2] = u0 s^2 + u1 r s + u2 r^2 as in
     rank_second.qv) are parametrised from their stored rational points p_X, p_Y by binary quadratic forms:
     (r, s, X) = (p1, p2, p3)(c, d), (r, s, Y) = (q1, q2, q3)(e, g)  (Fisher (22)).  Checked: p3^2 = lam U(p1, p2)
     and q3^2 = lam V(q1, q2) identically, and l_X(p1, p2, p3) = kX L1(c, d)^2, l_Y(q1, q2, q3) = kY L2(e, g)^2
     identically (the stored l_X, l_Y are the tangent lines at p_X, p_Y; kX, kY non-zero integers).
 (2) The level-3 curve as a quadric intersection in P^3(c : d : e : g) (Fisher (23), twisted):
         C3':  p1(c, d) = Xi q1(e, g),  p2(c, d) = Xi q2(e, g),     Xi = xi kX kY,
     with the morphism C3' -> C3 = {X^2 = lam U, Y^2 = lam V, xi w^2 = l_X l_Y} (rank_fisher.py),
         (r : s : X : Y ; w) = (p1 : p2 : p3 : Xi q3 ; kX kY L1 L2),
     checked by the identities of (1) (on C3' one has (r, s) = (p1, p2) = Xi (q1, q2); then X^2 - lam U = 0,
     Y^2 - lam V = Xi^2 (q3^2 - lam V(q1, q2)) = 0, xi w^2 - l_X l_Y = xi (kX kY L1 L2)^2 - kX L1^2 Xi kY L2^2 = 0).
     It is of degree 1: C3' -> D_lam has the fibres {(c:d:e:g), (c:d:-e:-g)} and w changes sign between them.
 (3) The pencil: det(x A + y B), A, B the symmetric matrices of the two quadrics, is checked to equal
     const * U(x, y) V(x, y) (in the orientation found and printed), with const != 0 and disc(U V) != 0, so C3' is a
     non-singular quadric intersection (Fisher Sect. 7, the non-singularity criterion after Thm 7.2) and its four
     singular quadrics lie over the roots of U and of V.
 (4') Fisher Sect. 5 check: disc V is the class dd being lifted (modulo squares), on all ten curves.
 (4) Hypothesis (19) of Fisher Sect. 7, Im(rho) in <(12), (34)>: disc U, disc V and disc U disc V are checked to be
     non-squares (exact isqrt), so the Galois group acts on the four singular quadrics through the full group
     <(12), (34)> of order 4; the partition {roots of U | roots of V} is the only Galois-stable one, so it is the one
     attached to the rational 2-torsion point T of Pic^0(C3) (the three partitions <-> the three points of order 2;
     exactly one is rational).  (19) holds; both pairs are non-split: L = Q(sqrt disc U), L' = Q(sqrt disc V).
 (5) The conic Gamma_1 over L: at the root w0 of U(x, 1) (x/y = w0), the (c, d)-block of the singular quadric is
     checked to have rank 1, h0 (c + rho d)^2, and Gamma_1 = { h0 Z^2 - Xi B(e, g) = 0 } over L,
     B = (x q1 + y q2)-part at the root.  Printed: the field and the sizes of the coefficients.
     Gamma_1 has an L-point: C3 is everywhere locally soluble (rank_fisher.py), a Q_v-point of C3' off the vertex
     projects to a point of Gamma_1 over L (x) Q_v, so Gamma_1 is locally soluble at every place of L, and a conic over
     a number field satisfies the Hasse principle.  (Written argument, not a computation.)
 (6) Where rank_fisher5_hints.json holds a point P1 of Gamma_1 over L (found by rank_fisher5_find.py: PARI and a
     Legendre descent in Z[sqrt -2]; only (72,65), E' side, so far), the point is checked on the cone, l_1 = the
     tangent plane there (gradient), tangency along the ruling checked (the quadric restricted to l_1 = 0 has rank
     <= 1, all 2x2 minors exactly 0 in L), and f = l_1 conj(l_1) is checked to be rational: by Lemma 7.1 a pushout
     form on C3' for T (divisor 2(a_1 + a_2 - H)).  Written to rank_fisher5_pushout.json (--write).
 (7) The same for the other pair (a point P3 of Gamma_3 over L' = Q(sqrt disc V), hint '<name>|<side>|V';
     Legendre descent in Z[sqrt 2] for (72,65)): l_3, f_dag = l_3 conj(l_3), a second pushout form for T
     (divisor 2(a_3 + a_4 - H)); Fisher's companion C4^dagger (Thm 7.4) is built from it.
 (8) With both: the quadric q through the eight points a_1 + a_2 + a_3 + a_4 (each a_i = C3' on the ruling of
     the cone through P_i: the restriction of q to the ruling is proportional to that of QA; 8 rational linear
     conditions of rank 7, nullspace span(QA, QB, q), checked), and the constant c with f f_dag - c q^2 in I(C3')
     (an exact linear system in c and two quadrics; solvable, c != 0).  So f / f_dag = c h^2 with h = q / f_dag
     (Fisher Sect. 7 before Thm 7.4: the companion D^dagger is xi c w^2 = f_dag when D is xi w^2 = f).
NOT COMPUTED: steps (6)-(8) for the other nine curves; the level-4 Selmer set S(C3', f), local points on C4, the companion C4^dagger (Thm 7.4) and Theta_4.
Formulas used from Fisher's paper (read through WebFetch extractions of arxiv.org/pdf/1509.03234v1, not the PDF
itself): (19); Lemma 7.1; (20)-(23) and Thm 7.2; Thm 7.3 with (24); Thm 7.4 with (25).
"""
import os, sys, json
from fractions import Fraction as Fr
from math import isqrt
import sympy as sp

WRITE = '--write' in sys.argv
here = os.path.dirname(os.path.abspath(__file__))
HINTS_FILE = os.path.join(here, 'rank_fisher5_hints.json')
HINTS = json.load(open(HINTS_FILE)) if os.path.exists(HINTS_FILE) else {}
OUT = []
PUSH = {}
NEXTRA = 0
def say(s=''):
    print(s, flush=True); OUT.append(s)

c, d, e, g, x, y, Z = sp.symbols('c d e g x y Z')

def qv(F, r, s): return F[0] * s * s + F[1] * r * s + F[2] * r * r        # rank_second.qv

def issq(n):
    n = sp.Rational(n)
    if n < 0: return False
    p, q = int(n.p), int(n.q)
    return isqrt(p) ** 2 == p and isqrt(q) ** 2 == q

def parametrise(lam, W, P, l, a, b):
    """the conic T^2 = lam W(r, s) through P = (r0, s0, T0) (T0 != 0), lines through P meeting T = 0 at (a : b : 0):
    (r, s, T) = F(Q) P - 2 B(P, Q) Q, Q = (a, b, 0).  Returns (p1, p2, p3), k, Lin with l(p) = k Lin^2, all checked."""
    r0, s0, T0 = P
    assert T0 ** 2 - lam * qv(W, r0, s0) == 0, 'the stored point is not on the conic'
    assert T0 != 0
    F = lambda r, s, T: T * T - lam * qv(W, r, s)
    FQ = sp.expand(F(a, b, 0))
    BPQ = sp.expand((F(r0 + a, s0 + b, T0) - F(r0, s0, T0) - FQ) / 2)
    p = [sp.expand(FQ * r0 - 2 * BPQ * a), sp.expand(FQ * s0 - 2 * BPQ * b), sp.expand(FQ * T0)]
    assert sp.expand(p[2] ** 2 - lam * qv(W, p[0], p[1])) == 0, 'parametrisation identity fails'
    lp = sp.expand(l[0] * p[0] + l[1] * p[1] + l[2] * p[2])
    # l is the tangent line at P iff l(p) is a constant times BPQ^2
    k = sp.cancel(lp / sp.expand(BPQ ** 2))
    assert k.is_number and k != 0, 'l is not the tangent line at the stored point'
    assert sp.expand(lp - k * BPQ ** 2) == 0
    return p, sp.Rational(k), BPQ

def sym_matrix(Q, vs):
    return sp.Matrix(len(vs), len(vs), lambda i, j: sp.Rational(sp.diff(Q, vs[i], vs[j])) / 2)


def cone_V(A, B, uv, V, dV):
    """the singular quadric of the pencil at a root of uv(V)(x, 1): exactly as step (5), with the roles of the
    blocks swapped ((e, g)-block of rank 1).  Returns (red, t, Mw, h0, Bq) over L' = Q(t), t^2 = dV:
    Gamma_3 = { h0 Z^2 + Bq(c, d) = 0 }, Z = e + (Mw[2, 3] / Mw[2, 2]) g."""
    Vx = sp.Poly(sp.expand(uv(V).subs(y, 1)), x)
    assert Vx.degree() == 2
    a2, a1, a0 = [sp.Rational(t_) for t_ in Vx.all_coeffs()]
    t = sp.Symbol('t')
    assert sp.expand(a1 ** 2 - 4 * a2 * a0 - dV) == 0
    red = lambda E_: sp.expand(sp.Poly(sp.expand(E_), t).rem(sp.Poly(t ** 2 - dV, t)).as_expr())
    w0 = (-a1 + t) / (2 * a2)
    assert red(a2 * w0 ** 2 + a1 * w0 + a0) == 0
    Mw = (w0 * A + B).applyfunc(red)
    assert all(red(Mw[i, j]) == 0 for i in range(2) for j in range(2, 4)), 'the pencil is not block diagonal'
    assert red(Mw[2, 2] * Mw[3, 3] - Mw[2, 3] ** 2) == 0, 'the (e, g)-block is not of rank 1'
    assert red(Mw[2, 2]) != 0
    assert red(Mw[0, 0] * Mw[1, 1] - Mw[0, 1] ** 2) != 0, 'the (c, d)-block is singular'
    h0 = Mw[2, 2]
    Bq = sp.expand(Mw[0, 0] * c ** 2 + 2 * Mw[0, 1] * c * d + Mw[1, 1] * d ** 2)
    return red, t, Mw, h0, Bq

def tangent_norm_form(Mw, red, t, P, vs, what):
    """P a point (dict) of the cone sum Mw[i,j] v_i v_j = 0 over Q(t), not its vertex: the tangent plane l there,
    tangency along the ruling (rank <= 1 on l = 0, all 2x2 minors exactly 0), and N(l) = l conj(l) rational."""
    Qw = sp.expand(sum(Mw[i, j] * vs[i] * vs[j] for i in range(4) for j in range(4)))
    assert red(Qw.subs(P)) == 0, '%s is not on the cone' % what
    l = sp.expand(sum(red(sp.diff(Qw, v_).subs(P)) * v_ for v_ in vs))
    co = {v_: red(sp.Poly(l, *vs).coeff_monomial(v_)) for v_ in vs}
    assert any(co[v_] != 0 for v_ in vs), '%s is the vertex' % what
    piv = next(v_ for v_ in vs if co[v_] != 0)
    others = [v_ for v_ in vs if v_ != piv]
    sub = sp.expand(-sum(co[v_] * v_ for v_ in others))
    Qr = red(sp.expand(Qw.subs({**{v_: co[piv] * v_ for v_ in others}, piv: sub}, simultaneous=True)))
    M3 = sp.Matrix(3, 3, lambda i, j: sp.diff(Qr, others[i], others[j]) / 2)
    for i in range(3):
        for j in range(i + 1, 3):
            for k in range(3):
                for l_ in range(k + 1, 3):
                    assert red(M3[i, k] * M3[j, l_] - M3[i, l_] * M3[j, k]) == 0, \
                        'the plane at %s is not tangent along a ruling' % what
    fq = red(sp.expand(l * l.subs(t, -t)))
    assert not fq.has(t), 'N(l) is not rational'
    assert sp.expand(fq) != 0
    return l, sp.expand(fq)

def quadric_through_line(qcoef, monos, vs, V0, P0, QA, red, t):
    """the L-linear conditions (as rational rows) that the quadric sum qcoef_k monos_k restricted to the line
    {s V0 + u P0} is proportional to QA restricted to it (i.e. vanishes at the two points of C3' on the line)."""
    s_, u_ = sp.symbols('s_ u_')
    pt = {v_: s_ * V0[i] + u_ * P0[i] for i, v_ in enumerate(vs)}
    al = sp.Poly(red(sp.expand(QA.subs(pt, simultaneous=True))), s_, u_)
    a = [red(al.coeff_monomial(m_)) for m_ in (s_ ** 2, s_ * u_, u_ ** 2)]
    assert any(x_ != 0 for x_ in a), 'QA vanishes on the ruling'
    rows = []
    mq = [sp.Poly(red(sp.expand(mo.subs(pt, simultaneous=True))), s_, u_) for mo in monos]
    M = [[red(P_.coeff_monomial(m_)) for P_ in mq] for m_ in (s_ ** 2, s_ * u_, u_ ** 2)]
    for (i, j) in ((0, 1), (0, 2), (1, 2)):
        row = [red(a[j] * M[i][k] - a[i] * M[j][k]) for k in range(len(monos))]
        for part in (lambda z: sp.Poly(z, t).coeff_monomial(1), lambda z: sp.Poly(z, t).coeff_monomial(t)):
            rows.append([sp.Rational(part(z)) for z in row])
    return rows

def analyse(name, side, C, mutate=None, dd=None):
    lam, xi = int(C['lam']), int(C['xi'])
    U = [int(t) for t in C['U']]; V = [int(t) for t in C['V']]
    lX = [int(t) for t in C['l_X']]; lY = [int(t) for t in C['l_Y']]
    pX = [int(t) for t in C['p_X']]; pY = [int(t) for t in C['p_Y']]
    if mutate == 'U': U = [U[0] + 1] + U[1:]
    if mutate == 'pX': pX = [pX[0] + 1] + pX[1:]
    if mutate == 'lX': lX = [lX[0] + 1] + lX[1:]
    # (1)
    pp, kX, L1 = parametrise(lam, U, pX, lX, c, d)
    qq, kY, L2 = parametrise(lam, V, pY, lY, e, g)
    Xi = xi * kX * kY
    # (2)
    QA = sp.expand(pp[0] - Xi * qq[0]); QB = sp.expand(pp[1] - Xi * qq[1])
    if mutate == 'Xi': QB = sp.expand(pp[1] - 3 * Xi * qq[1])
    # the morphism C3' -> C3: on C3', (r, s) = (p1, p2) = Xi (q1, q2)
    r_, s_, X_, Y_, w_ = pp[0], pp[1], pp[2], Xi * qq[2], kX * kY * L1 * L2
    assert sp.expand(X_ ** 2 - lam * qv(U, r_, s_)) == 0
    assert sp.expand(Y_ ** 2 - lam * qv(V, Xi * qq[0], Xi * qq[1])) == 0
    lXv = sp.expand(lX[0] * r_ + lX[1] * s_ + lX[2] * X_)
    lYv = sp.expand(lY[0] * Xi * qq[0] + lY[1] * Xi * qq[1] + lY[2] * Y_)
    assert sp.expand(xi * w_ ** 2 - lXv * lYv) == 0, 'xi w^2 = l_X l_Y fails on the model'
    # (3) the pencil
    vs = (c, d, e, g)
    A, B = sym_matrix(QA, vs), sym_matrix(QB, vs)
    quart = sp.Poly(sp.expand((x * A + y * B).det()), x, y)
    found = None
    for (lab, uv) in (('U(x,y) V(x,y)', lambda F: qv(F, x, y)), ('U(y,x) V(y,x)', lambda F: qv(F, y, x)),
                      ('U(-y,x) V(-y,x)', lambda F: qv(F, -y, x)), ('U(y,-x) V(y,-x)', lambda F: qv(F, y, -x))):
        target = sp.Poly(sp.expand(uv(U) * uv(V)), x, y)
        ratio = sp.cancel(quart.as_expr() / target.as_expr())
        if ratio.is_number and ratio != 0:
            found = (lab, sp.Rational(ratio), uv); break
    assert found, 'det(x A + y B) is not a constant times U V'
    lab, const, uv = found
    UVd = sp.discriminant(sp.expand(uv(U) * uv(V)).subs(y, 1), x)
    assert UVd != 0 and sp.Poly(sp.expand(uv(U) * uv(V)), x, y).total_degree() == 4
    # (4) hypothesis (19)
    dU = U[1] ** 2 - 4 * U[0] * U[2]; dV = V[1] ** 2 - 4 * V[0] * V[2]
    if mutate == 'square': dU = 4
    assert not issq(dU) and not issq(dV) and not issq(dU * dV), 'a pair of singular quadrics splits: not the case treated'
    if dd is not None:
        # Fisher Sect. 5: the cone pairs of the level-3 curve lie over Q(sqrt xi_1) and Q(sqrt(b / xi_1)), xi_1 the
        # class being lifted; here the V-pair field is Q(sqrt dd) (exact)
        assert issq(dV * dd), 'the V-pair field is not Q(sqrt dd)'
    # (5) Gamma_1 at a root of uv(U)(x, 1)
    Ux = sp.Poly(sp.expand(uv(U).subs(y, 1)), x)
    assert Ux.degree() == 2
    a2, a1, a0 = [sp.Rational(t) for t in Ux.all_coeffs()]
    t = sp.Symbol('t')                    # t^2 = dU: elements of L as a + b t
    assert sp.expand(a1 ** 2 - 4 * a2 * a0 - dU) == 0      # the discriminant of Ux is dU
    red = lambda E_: sp.expand(sp.Poly(sp.expand(E_), t).rem(sp.Poly(t ** 2 - dU, t)).as_expr())
    w0 = (-a1 + t) / (2 * a2)
    assert red(a2 * w0 ** 2 + a1 * w0 + a0) == 0
    # the singular member x A + y B at (x : y) = (w0 : 1)
    Mw = (w0 * A + B).applyfunc(red)
    assert all(red(Mw[i, j]) == 0 for i in range(2) for j in range(2, 4)), 'the pencil is not block diagonal'
    assert red(Mw[0, 0] * Mw[1, 1] - Mw[0, 1] ** 2) == 0, 'the (c, d)-block is not of rank 1'
    assert red(Mw[0, 0]) != 0
    assert red(Mw[2, 2] * Mw[3, 3] - Mw[2, 3] ** 2) != 0, 'the (e, g)-block is singular'
    h0 = Mw[0, 0]
    # Gamma_1 (cone over it): h0 Z^2 + (e, g) block (e, g)^T = 0 with h0 Z = h0 c + Mw[0, 1] d
    Bq = sp.expand(Mw[2, 2] * e ** 2 + 2 * Mw[2, 3] * e * g + Mw[3, 3] * g ** 2)
    rho = None
    hgt = max(len(str(abs(t_))) for t_ in (int(sp.Rational(Xi).p), dU, dV))
    pushout = None
    H = HINTS.get('%s|%s' % (name, side))
    if H and mutate in (None, 'P1', 'c'):
        # (6) a point P1 = (Z, 0, e1, g1) of the cone over Gamma_1 over L (a hint, found by PARI; checked here), the
        # tangent plane l_1 = grad Q_w (P1) . (c, d, e, g), and f = l_1 * conj(l_1) (Fisher Lemma 7.1)
        D0, m = int(H['D0']), int(H['m'])
        assert m * m * D0 == dU
        el = lambda v: sp.Rational(v[0]) + sp.Rational(v[1]) * t / m      # a + b sqrt(D0), sqrt(D0) = t / m
        Z1, e1, g1 = el(H['Z']), el(H['e']), el(H['g'])
        if mutate == 'P1': Z1 = Z1 + 1
        Qw = sp.expand(sum(Mw[i, j] * vs[i] * vs[j] for i in range(4) for j in range(4)))
        P1 = {c: Z1, d: 0, e: e1, g: g1}
        assert red(Qw.subs(P1)) == 0, 'P1 is not on the cone over Gamma_1'
        assert red(e1) != 0 or red(g1) != 0, 'P1 is the vertex'
        l1 = sp.expand(sum(red(sp.diff(Qw, v_).subs(P1)) * v_ for v_ in vs))
        # tangency: on the plane l_1 = 0 the quadric Q_w has rank <= 1 (Q_w = 0 there is the ruling through P1, doubled)
        # eliminate a variable with a non-zero coefficient (exactly in L) and check the restricted 3x3 form has rank <= 1
        co = {v_: red(sp.Poly(l1, *vs).coeff_monomial(v_)) for v_ in vs}
        piv = next(v_ for v_ in (e, g, c, d) if co[v_] != 0)
        # v_piv = -(sum of the others)/co[piv]; work with N = co[piv] times it to stay polynomial in t
        others = [v_ for v_ in vs if v_ != piv]
        sub = sp.expand(-sum(co[v_] * v_ for v_ in others))
        # Q_w(others, sub / co) co^2 = Q_w(co * others, sub) (homogeneous of degree 2): polynomial in t
        Qr = red(sp.expand(Qw.subs({**{v_: co[piv] * v_ for v_ in others}, piv: sub}, simultaneous=True)))
        # Qr is a ternary form over L; rank <= 1 iff every 2x2 minor of its matrix vanishes
        M3 = sp.Matrix(3, 3, lambda i, j: sp.diff(Qr, others[i], others[j]) / 2)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l_ in range(3):
                        if i < j and k < l_:
                            assert red(M3[i, k] * M3[j, l_] - M3[i, l_] * M3[j, k]) == 0, \
                                'l_1 is not tangent to the cone along a ruling'
        conj = lambda E_: E_.subs(t, -t)
        fq = red(sp.expand(l1 * conj(l1)))
        assert sp.Poly(fq, t).degree() <= 0 if fq.has(t) else True
        assert not fq.has(t), 'f = N(l_1) is not rational'
        fq = sp.expand(fq)
        assert fq != 0
        # the primitive integral form with the same square class up to the printed rational scalar
        Pf = sp.Poly(fq, *vs)
        den = sp.ilcm(*[sp.Rational(c_).q for c_ in Pf.coeffs()])
        num = sp.igcd(*[int(sp.Rational(c_) * den) for c_ in Pf.coeffs()])
        fprim = sp.expand(fq * den / num)
        pushout = dict(l1=l1, f=fq, fprim=fprim, scalar=sp.Rational(num, den))
    coneV = cone_V(A, B, uv, V, dV)
    pushout3 = None
    H3 = HINTS.get('%s|%s|V' % (name, side))
    if H3 and mutate in (None, 'P3', 'c'):
        # (7) the same for the other pair: a point P3 = (c3, d3, Z, 0) of the cone over Gamma_3 over L' (hint),
        # l_3 its tangent plane, f_dag = l_3 conj(l_3): a second pushout form for T (Fisher Thm 7.4 uses it)
        redV, tV, MwV, h0V, BqV = coneV
        D0, m = int(H3['D0']), int(H3['m'])
        assert m * m * D0 == dV
        elV = lambda v: sp.Rational(v[0]) + sp.Rational(v[1]) * tV / m
        Z3 = elV(H3['Z']) + (1 if mutate == 'P3' else 0)
        P3 = {c: elV(H3['c']), d: elV(H3['d']), e: Z3, g: 0}
        l3, f3 = tangent_norm_form(MwV, redV, tV, P3, vs, 'P3')
        Pf = sp.Poly(f3, *vs)
        den = sp.ilcm(*[sp.Rational(c_).q for c_ in Pf.coeffs()])
        num = sp.igcd(*[int(sp.Rational(c_) * den) for c_ in Pf.coeffs()])
        pushout3 = dict(l3=l3, f=f3, fprim=sp.expand(f3 * den / num), scalar=sp.Rational(num, den))
    comp = None
    if pushout is not None and pushout3 is not None and mutate in (None, 'c'):
        # (8) the quadric q through a_1 + a_2 + a_3 + a_4 (unique modulo I(C3'), since a_1 + ... + a_4 ~ 2H) and the
        # constant c with f f_dag = c q^2 modulo I(C3') (so f / f_dag = c (q / f_dag)^2: Fisher Sect. 7 before Thm 7.4)
        monos = [vs[i] * vs[j] for i in range(4) for j in range(i, 4)]
        V1 = [-Mw[0, 1], Mw[0, 0], 0, 0]; P1v = [Z1, 0, e1, g1]
        redV, tV, MwV, h0V, BqV = coneV
        V3 = [0, 0, -MwV[2, 3], MwV[2, 2]]; P3v = [P3[c], P3[d], P3[e], 0]
        for (V0, Mx, rd) in ((V1, Mw, red), (V3, MwV, redV)):
            assert all(rd(sum(Mx[i, j] * V0[j] for j in range(4))) == 0 for i in range(4)), 'not the vertex'
        rows = quadric_through_line(None, monos, vs, V1, P1v, QA, red, t) + \
               quadric_through_line(None, monos, vs, V3, P3v, QA, redV, tV)
        Mrows = sp.Matrix(rows)
        ns = Mrows.nullspace()
        assert len(ns) == 3, ('the 8 points do not impose 7 conditions', len(ns))
        inI = lambda vec: sp.Matrix.hstack(sp.Matrix([sp.Poly(QA, *vs).coeff_monomial(mo) for mo in monos]),
                                           sp.Matrix([sp.Poly(QB, *vs).coeff_monomial(mo) for mo in monos]), vec).rank()
        qv_ = next(v_ for v_ in ns if inI(v_) == 3)
        q = sp.expand(sum(qv_[k] * monos[k] for k in range(10)))
        # c: f f_dag - c q^2 = QA mA + QB mB, a linear system in c and the 20 coefficients of mA, mB
        cc = sp.Symbol('cc'); mA = sp.symbols('mA0:10'); mB = sp.symbols('mB0:10')
        f_used = pushout['f'] + (e ** 2 if mutate == 'c' else 0)
        expr = sp.expand(f_used * pushout3['f'] - cc * q ** 2 - QA * sum(mA[k] * monos[k] for k in range(10))
                         - QB * sum(mB[k] * monos[k] for k in range(10)))
        eqs = sp.Poly(expr, *vs).coeffs()
        sol = sp.linsolve(eqs, [cc] + list(mA) + list(mB))
        assert sol != sp.EmptySet and len(sol) > 0, 'f f_dag is not a constant times a square modulo I(C3\')'
        cval = list(sol)[0][0]
        assert cval.is_number and cval != 0, ('c not determined', cval)
        comp = dict(q=q, c=sp.Rational(cval))
    return dict(dd=dd, comp=comp, coneV=coneV, A=A, B=B, pushout=pushout, pushout3=pushout3, name=name, side=side, lam=lam, xi=xi, kX=kX, kY=kY, Xi=Xi, lab=lab, const=const, dU=dU, dV=dV,
                h0=h0, Bq=Bq, QA=QA, QB=QB, hgt=hgt)

def core_small(n):
    """the squarefree part D0 of n (sympy factorint), checked: n = m^2 D0."""
    n = int(n); D0, m = (-1 if n < 0 else 1), 1
    for p, k in sp.factorint(abs(n)).items():
        D0 *= p ** (k % 2); m *= p ** (k // 2)
    assert m * m * D0 == n
    return D0

def main():
    say('rank_fisher5.py: first steps of Fisher\'s fifth descent (arXiv:1509.03234 Sect. 7) on the level-3 curves '
        'of rank_fisher.py, exact')
    js = json.load(open(os.path.join(here, 'rank_fisher_curves.json')))
    for cv in js:
        red = cv['theta4_reduction']
        for side, S in cv['sides'].items():
            C = S.get('level3_curve')
            if not C: continue
            R = analyse(cv['name'], side, C, dd=int(S['dd']))
            tag = ' (the C3 of the Theta_4 reduction)' if red and side == red[0] else ''
            say('%s, %s side%s: lam = %d, xi = %d' % (R['name'], side, tag, R['lam'], R['xi']))
            say('  (1) conics parametrised from p_X, p_Y; l_X, l_Y tangent: l_X(p) = kX L1^2, l_Y(q) = kY L2^2 with '
                'kX = %s, kY = %s' % (R['kX'], R['kY']))
            say('  (2) C3\': p1(c,d) = Xi q1(e,g), p2(c,d) = Xi q2(e,g), Xi = xi kX kY (%d digits); '
                'the morphism to C3 checked (identities)' % len(str(abs(int(sp.Rational(R['Xi']).p)))))
            say('  (3) det(x A + y B) = %s * %s (const %d digits): non-singular quadric intersection'
                % ('c' if True else '', R['lab'], len(str(abs(int(R['const'].p))))))
            say('  (4) (19) holds: disc U, disc V, disc U disc V non-squares; L = Q(sqrt %d), L\' = Q(sqrt %d) '
                '(dU = m^2 D0 checked); disc V = dd = %d modulo squares (Fisher Sect. 5: one pair over Q(sqrt xi_1))'
                % (core_small(R['dU']), core_small(R['dV']), R['dd']))
            if R['pushout'] is None:
                say('  (5) Gamma_1 over L: h0 Z^2 + B(e, g) = 0, (c, d)-block of rank 1 checked; locally soluble '
                    'everywhere because C3 is (written argument); an L-point is NOT computed here')
            else:
                say('  (5) Gamma_1 over L: h0 Z^2 + B(e, g) = 0, (c, d)-block of rank 1 checked; an L-point P1 (hint) '
                    'checked on the cone')
                say('  (6) Lemma 7.1: l_1 = tangent plane at P1 (tangency along the ruling checked); the pushout form '
                    'on C3\' is f = l_1 conj(l_1) = s0 * f0, f0 primitive in Z[c, d, e, g] (written to '
                    'rank_fisher5_pushout.json), s0 = %s' % R['pushout']['scalar'])
                PUSH['%s|%s' % (R['name'], side)] = dict(
                    C3prime=dict(QA=str(R['QA']), QB=str(R['QB']), vars='c d e g'), f0=str(R['pushout']['fprim']),
                    s0=str(R['pushout']['scalar']), L='Q(sqrt(%d))' % core_small(R['dU']))
            if R['pushout3'] is not None:
                say('  (7) the other pair: a point P3 of the cone over Gamma_3 over L\' (hint) checked, l_3 tangent along '
                    'the ruling (checked), f_dag = l_3 conj(l_3) = s3 * f3 rational (rank_fisher5_pushout.json), '
                    's3 = %s' % R['pushout3']['scalar'])
                PUSH.setdefault('%s|%s' % (R['name'], side), {}).update(
                    f3=str(R['pushout3']['fprim']), s3=str(R['pushout3']['scalar']),
                    Lp='Q(sqrt(%d))' % core_small(R['dV']))
            if R['comp'] is not None:
                cv_ = R['comp']['c']
                say('  (8) q through a_1 + a_2 + a_3 + a_4 (8 points, 7 conditions, checked) and f f_dag = c q^2 mod '
                    'I(C3\') with c = %s (exact linear algebra)' % cv_)
                PUSH['%s|%s' % (R['name'], side)].update(q=str(R['comp']['q']), c=str(cv_))
    say('Theta_4 is not computed by this script.')
    if WRITE:
        open(os.path.join(here, 'rank_fisher5.txt'), 'w').write('\n'.join(OUT) + '\n')
        json.dump(PUSH, open(os.path.join(here, 'rank_fisher5_pushout.json'), 'w'), indent=1, sort_keys=True)

def negative():
    js = json.load(open(os.path.join(here, 'rank_fisher_curves.json')))
    C = js[2]['sides']["E'"]['level3_curve']          # (72,65), E' side: the smallest
    try:
        analyse(js[2]['name'], "E'", C, dd=-int(js[2]['sides']["E'"]['dd']))
        say('NOT REFUSED: the class dd replaced by -dd (Sect. 5 field check)')
    except AssertionError as ex:
        say('refused: the class dd replaced by -dd (%s)' % str(ex)[:80])
        global NEXTRA
        NEXTRA = 1
    n = 0
    for m, what in (('U', 'a mutated U (the stored point leaves the conic)'), ('pX', 'a point off the conic'),
                    ('lX', 'a line that is not the tangent at p_X'), ('Xi', 'the second quadric with Xi -> 3 Xi'),
                    ('square', 'disc U forged to a square (a split pair)'),
                    ('P1', 'the hinted point of Gamma_1 moved off the cone (Z -> Z + 1)'),
                    ('P3', 'the hinted point of Gamma_3 moved off the cone (Z -> Z + 1)'),
                    ('c', 'f replaced by f + e^2 in step (8) (no longer c q^2 / f_dag modulo I(C3\'))')):
        try:
            analyse(js[2]['name'], "E'", C, mutate=m)
            say('NOT REFUSED: ' + what)
        except AssertionError as ex:
            n += 1; say('refused: %s (%s)' % (what, str(ex)[:80]))
    n += NEXTRA
    say('%d of 9 refused' % n)
    assert n == 9
    if WRITE: open(os.path.join(here, 'rank_fisher5_negative.txt'), 'w').write('\n'.join(OUT) + '\n')

if __name__ == '__main__':
    if '--negative' in sys.argv: negative()
    else: main()
