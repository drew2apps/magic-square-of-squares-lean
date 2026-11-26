#!/usr/bin/env python3
"""The Cassels-Tate pairing at the level of the 2-Selmer group, for the curves whose rank bound 1 rests on Magma alone.

    python3 rank_second.py              the five curves (four E+ of sigma_evidence.txt and D1 of (631,71)), a control,
                                        and the curves E0L of the four-point tier for (72,65) and (920,129)
    python3 rank_second.py --negative   altered inputs that must be refused
    (add --write to rewrite rank_second.txt / rank_second_negative.txt and the hints; without it the script only prints)

Written by an AI model (Claude); not read by a person.  Experiment, not counted.  It imports rank_indep.py (the exact
descent via 2-isogeny and the exact second descent) and adds one step to it.

THE QUESTION.  For each curve E: y^2 = x(x^2 + a x + b), E': Y^2 = X(X^2 + a' X + b'), the exact descent of
rank_indep.py leaves one coset of classes d | b (quartics C_d, torsors under E') and one coset of classes d' | b'
(quartics C'_d', torsors under E) that the second descent does not exclude, hence the bound 3; Magma's
TwoPowerIsogenyDescentRankBound gives 1, i.e. both cosets are not in the image of the Mordell-Weil group.  If rank = 1
these classes are elements x of Sha(E')[psi] and y0 of Sha(E)[phi], and Sha(E)[2] = <y0, y1> with phi(y1) = x.  The
one part of the Cassels-Tate pairing on Sel^2(E) that the second descents do not already compute is <y0, y1>.  If it
were 1/2, both cosets would be excluded and the bound would be 1 (this is the "pairing on the 2-Selmer group" of item
B2).  This script computes it exactly.

THE FORMULA.  A lift y1 of C_d to Sel^2(E) is a curve D_lam: X^2 = lam U(r, s), Y^2 = lam V(r, s) of the second
descent (rank_indep.second_descent), a torsor under E, with the covering D_lam -> C_d of deck transformation
(X, Y) -> (-X, -Y), the translation by T = (0, 0) of E.  The conics C_X: W^2 = lam U(r, s) and C_Y: W^2 = lam V(r, s)
have rational points p_X, p_Y (found by PARI's qfsolve as hints; checked here); with l_X, l_Y the tangent lines there,
div(l_X) = 2 F_X and div(l_Y) = 2 F_Y on D_lam, where F_X, F_Y are the fibres of D_lam -> C_X, C_Y, and
F_Y - F_X is the class of T (the reflections (X, Y) -> (X, -Y), (-X, Y) have fibre classes F_X, F_Y and their
product is the translation by T).  So g = l_Y / l_X has divisor 2 Delta with Delta in the class of T, and for
d' in Sel^phi(E) = H^1(Q, E[phi]) the quaternion algebra (d', g) is unramified on D_lam; its Brauer-Manin pairing
    <d', D_lam> = prod_v (d', g(P_v))_v,   P_v in D_lam(Q_v),
is the Cassels-Tate pairing of the class of C'_d' with the class of D_lam (Manin; see Skorobogatov, Torsors and
rational points, section 6.2, and Cassels, J. reine angew. Math. 494 (1998), for 2-isogenies).  It does not depend on
the P_v because d' is in the Selmer group; the script checks this at several points of each D_lam(Q_v) with all four
sign choices of (X, Y) (a check, not a proof).  Outside the places S = {inf, 2, the primes of lam, d', delta and
det of the two tangent decompositions k Q = l m + delta n^2, and (since v4) the bad primes of D_lam}, every local factor is 1 (for a primitive integral point
an odd valuation of l(P) forces p | l, m, n, hence p | r, s, X and p | Y; written out in local_factor_one below).
  The theory (Cassels-Tate pairing, its Brauer-Manin description, and that the second descent computes the pairing on
Sha(E')[psi]) is used, not proved here.  The arithmetic is exact: integers, p-adic square roots by Hensel lifting
modulo p^N with the precision asserted, Hilbert symbols by their formulas, signs at the real place by exact
comparison of squares.

TWO CHECKS OF THE FORMULA.  (1) The same number computed on the other side, <d, D'_mu>, with D'_mu the lifts of
C'_d' and d the class of C_d (the pairing is <x, c> on E' with psi(c) = y0, equal to <y0, y1> on E); a different
curve, different forms.  (2) A control where the answer is known from rank_indep.py: on a curve where the second
descent on the E' side excludes a class d' of Sel^phi(E), the pairing of d' with the lifts D_lam of the trivial class
C_1 (these are the images of Sel^phi(E) in Sel^2(E)) must be -1 for some lam; for a class that survives it must be +1
for all.

WHAT IT GIVES.  For all five curves <y0, y1> = 0 (+1 below), on every lift with local points everywhere, on both
sides.  So the bound 3 is also the bound from the 2-Selmer groups of E and E' with the Cassels-Tate pairing, as
PARI's ellrank finds; if the rank is 1 (Magma), Sha(E)[2] = (Z/2)^2 lies in 2 Sha(E), and Magma's bound needs a
descent deeper than this one (Fisher's higher descents by 2-power isogenies).  The ten ratios stay "Magma alone".

THE FOUR-POINT CURVES E0L (full 2-torsion; dm4_log.jsonl).  (72,65), (137,7) and (920,129) are also closed by the
four-point argument on E0L, with rank 3 from PARI's and Magma's bounds.  For y^2 = (x - e1)(x - e2)(x - e3): Sel^2
by the complete 2-descent (local images of the right order found from explicit local points); on each 2-covering
D_a: x = e1 + a1 u^2 = e2 + a2 v^2 = e3 + a1 a2 w^2 the three conics obtained by forgetting u, v or w, their
tangent lines l_u, l_v, l_w, and <b, D_a> = prod_v (b2, l_v/l_w)_v (b1, l_u/l_w)_v (the classes F_v - F_w,
F_u - F_w are T1, T2; b = b2 (x) T1 + b1 (x) T2 under the Weil pairing).  Checked: symmetric, alternating, the
images of the rational points in the kernel, curves with a rational point pair trivially, a class outside Sel^2
gives a non-constant evaluation (--negative).  E0L of (920,129): dim Sel^2 = 7, the pairing has rank 2, the images
of the three points and the 2-torsion span dimension 5: rank = 3 exactly, with no PARI and no Magma for the bound
(PARI's ellrank on this curve: [1, 3], s = 2, agrees).  E0L of (72,65) = (137,7) (one curve): the pairing on its
Sel^2 is trivial (PARI: [1, 5], s = 0, agrees); PARI's bound 3 comes from the 2-isogenous curve B via (0, 0)
([3, 3], s = 2), where the non-trivial part is the pairing between two lifts of the second descent, not computed
here; so 3 <= rank <= 5 here, and that four-point closure still rests on PARI and Magma for the bound.
"""
import os, sys, json, random
WRITE = '--write' in sys.argv
from fractions import Fraction as Fr
from math import gcd, isqrt

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import rank_indep as RI
import sigma_exact as SE

HINTS_FILE = os.path.join(here, 'rank_second_hints.json')
HINTS = {'conic': {}}
NPREC = 40            # p-adic precision of the local points
NPTS = 3              # local points per place (each with the four sign choices)

OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(s)

# ---------------------------------------------------------------- Hilbert symbols
def vp(x, p): return RI.vp(x, p)

def hilbert(a, b, p):
    """(a, b)_p for nonzero rationals a, b; p a prime or 'inf'."""
    a = Fr(a); b = Fr(b)
    a = a.numerator * a.denominator; b = b.numerator * b.denominator
    assert a != 0 and b != 0
    if p == 'inf': return -1 if (a < 0 and b < 0) else 1
    al, be = vp(a, p), vp(b, p)
    u, w = a // p ** al, b // p ** be
    if p == 2:
        eps = lambda x: ((x % 8 - 1) // 2) % 2
        om = lambda x: (((x % 8) ** 2 - 1) // 8) % 2
        e = eps(u) * eps(w) + al * om(w) + be * om(u)
        return -1 if e % 2 else 1
    s = 1
    if (al * be) % 2 and p % 4 == 3: s = -s
    if be % 2: s *= RI.legendre(u, p)
    if al % 2: s *= RI.legendre(w, p)
    return s

def hilbert_selftest():
    """The formulas against the definition at p = 2, 3, 5: (a, b)_p = 1 iff a x^2 + b y^2 is a non-zero square of Q_p
    for some x, y (searched over integers 0 <= x, y < p^3, resp. 2^5; the square test on an integer is exact)."""
    n = 0
    for p in (2, 3, 5):
        mod = 32 if p == 2 else p ** 3
        for a in (1, -1, 2, -2, 3, 5, 6, 10, 15, -3, -5, 7, 14):
            for b in (1, -1, 2, 3, 5, 6, -6, 7, 10, -10, 21):
                sol = any(_is_padic_square(a * x * x + b * y * y, p) for x in range(mod) for y in range(mod)
                          if a * x * x + b * y * y != 0)
                assert (hilbert(a, b, p) == 1) == sol, ('hilbert symbol formula', a, b, p)
                n += 1
    return n

# ---------------------------------------------------------------- p-adic square roots
def sqrt_mod_p(u, p):
    import sympy
    return int(sympy.sqrt_mod(u % p, p))

def padic_sqrt(A, p, N):
    """X in Z with X^2 = A in Z_p modulo p^(N + v_p(A)/2) (A a nonzero integer, a square in Q_p)."""
    k = vp(A, p); assert k % 2 == 0
    u = A // p ** k
    if p == 2:
        assert u % 8 == 1
        x = 1
        for j in range(3, N + 3):
            if (x * x - u) % 2 ** (j + 1): x += 2 ** (j - 1)
        assert (x * x - u) % 2 ** (N + 2) == 0
        x %= 2 ** N
    else:
        assert RI.legendre(u, p) == 1
        x = sqrt_mod_p(u, p); mod = p
        while mod < p ** N:
            mod = mod * mod
            x = (x - (x * x - u) * pow(2 * x, -1, mod)) % mod
        x %= p ** N
        assert (x * x - u) % p ** N == 0
    return p ** (k // 2) * x

# ---------------------------------------------------------------- conics and tangent decompositions
def conic_point(lam, F):
    """A rational point (r, s, W) of W^2 = lam F(r, s), F = [coef s^2, coef r s, coef r^2] (hint from PARI's qfsolve,
    checked)."""
    key = '%d %s' % (lam, F)
    if key in HINTS['conic'] and HINTS['conic'][key] != 'none':
        P = [int(x) for x in HINTS['conic'][key]]
    else:
        if key not in HINTS['conic']:
            import cypari2
            pari = cypari2.Pari()
            M = pari.matrix(3, 3, [-2 * lam * F[2], -lam * F[1], 0, -lam * F[1], -2 * lam * F[0], 0, 0, 0, 2])
            sol = pari.qfsolve(M)
            if str(pari.type(sol)) == 't_COL':
                P = [int(x) for x in sol]
                g = gcd(gcd(P[0], P[1]), P[2]); P = [x // g for x in P]
                HINTS['conic'][key] = P
            else: HINTS['conic'][key] = 'none'
        if HINTS['conic'][key] == 'none':
            # checked exactly: W^2 = lam F(r, s) is a F2 (r + ...)^2 + b s^2 - W^2 = 0 with a = lam F2,
            # b = lam (4 F0 F2 - F1^2) / (4 F2); it has no point over Q_p iff (a, b)_p = -1
            assert F[2] != 0
            a_, b_ = Fr(lam * F[2]), Fr(lam * (4 * F[0] * F[2] - F[1] ** 2), 4 * F[2])
            ps = ['inf', 2] + sorted(set(RI.factor_int(a_.numerator * a_.denominator)) | set(RI.factor_int(b_.numerator * b_.denominator)))
            bad = [p for p in ps if hilbert(a_, b_, p) == -1]
            assert bad, 'a conic said to have no rational point is soluble everywhere locally'
            raise AssertionError(('no conic point', 'none over Q_%s' % bad[0]))
        P = [int(x) for x in HINTS['conic'][key]]
    r, s, W = P
    assert (r, s, W) != (0, 0, 0)
    assert W * W == lam * (F[0] * s * s + F[1] * r * s + F[2] * r * r), 'not a point of the conic'
    return P

def qform(lam, F):
    """W^2 - lam F(r, s) as a symmetric matrix A (Q(z) = z^T A z) over Fractions, variables (r, s, W)."""
    return [[Fr(-lam * F[2]), Fr(-lam * F[1], 2), Fr(0)],
            [Fr(-lam * F[1], 2), Fr(-lam * F[0]), Fr(0)],
            [Fr(0), Fr(0), Fr(1)]]

def qeval(A, z): return sum(A[i][j] * z[i] * z[j] for i in range(3) for j in range(3))
def lin(c, z): return sum(x * y for x, y in zip(c, z))
def det3(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))

def tangent_decomposition(A, P):
    """Integer linear forms l (the tangent at P), m, n (n(P) = 0) and integers k > 0, delta with
    k Q = l m + delta n^2 identically (checked on the six monomials), and det(l, m, n) != 0."""
    assert qeval(A, P) == 0, 'the point is not on the conic'
    l = [2 * sum(A[i][j] * P[j] for j in range(3)) for i in range(3)]
    assert any(l), 'singular point'
    # n: the cross product of P with a unit vector, independent of l
    for w in ([1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [0, 1, 1], [1, 0, 1], [1, 1, 1]):
        n = [P[1] * w[2] - P[2] * w[1], P[2] * w[0] - P[0] * w[2], P[0] * w[1] - P[1] * w[0]]
        if any(n) and any(l[i] * n[j] - l[j] * n[i] for i in range(3) for j in range(3)): break
    # a point R != P on the line l = 0 with n(R) != 0: then Q - delta n^2 vanishes on l = 0 iff it vanishes at R
    # (Q and n^2 both vanish doubly at P on that line)
    cands = [[l[1], -l[0], 0], [l[2], 0, -l[0]], [0, l[2], -l[1]]]
    R = next(c for c in cands if any(c) and lin(n, c) != 0 and any(c[i] * P[j] - c[j] * P[i] for i in range(3) for j in range(3)))
    delta = qeval(A, R) / Fr(lin(n, R)) ** 2
    # m = (Q - delta n^2) / l: solve for the coefficients of m by comparing the six monomials
    B = [[A[i][j] - delta * n[i] * n[j] for j in range(3)] for i in range(3)]       # symmetric matrix of Q - delta n^2
    # B = (l m^T + m l^T)/2: pick i0 with l[i0] != 0; m_j = 2 B[i0][j]/l[i0] - (l_j m_i0)/l[i0], m_i0 = B[i0][i0]/l[i0]
    i0 = next(i for i in range(3) if l[i] != 0)
    mi0 = B[i0][i0] / l[i0]
    m = [(2 * B[i0][j] - l[j] * mi0) / l[i0] if j != i0 else mi0 for j in range(3)]
    for i in range(3):
        for j in range(3):
            assert B[i][j] == (l[i] * m[j] + m[i] * l[j]) / 2, 'tangent decomposition fails (not a tangent line?)'
    den = 1
    for x in m + [delta]: den = den * x.denominator // gcd(den, x.denominator)
    k = den
    mI = [int(x * k) for x in m]; dI = int(delta * k); lI = [int(x) for x in l]
    for i in range(3):
        for j in range(3):
            assert k * A[i][j] == Fr(lI[i] * mI[j] + mI[i] * lI[j], 2) + dI * n[i] * n[j]
    dt = det3([lI, mI, n])
    assert dt != 0 and dI != 0
    return dict(P=P, l=lI, m=mI, n=n, k=k, delta=dI, det=dt)

def qv(F, r, s): return F[0] * s * s + F[1] * r * s + F[2] * r * r

def check_conic_forms(d, a, c, U, V, N0):
    """U, V parametrise the conic N^2 = d U^2 + a U V + c V^2 (identity of quartic forms, checked at 7 points)."""
    for (r_, s_) in ((1, 0), (0, 1), (1, 1), (1, -1), (2, 1), (1, 2), (3, -2)):
        u_, v_ = qv(U, r_, s_), qv(V, r_, s_)
        F_ = d * r_ * r_ + a * r_ * s_ + c * s_ * s_
        assert d * u_ * u_ + a * u_ * v_ + c * v_ * v_ == (F_ * N0) ** 2, 'the forms U, V do not parametrise the conic'

# ---------------------------------------------------------------- the pairing on one D_lam
class DLam:
    def __init__(s, lam, U, V, bad=()):
        s.lam, s.U, s.V = lam, U, V
        # v4 fix (review_t2 FIX-FIRST): the bad primes of D_lam (rank_indep.second_descent: 2, the primes of b' N0,
        # of disc U, disc V and of the pencil discriminant) are added to S, so that a local point of D_lam is
        # found and checked at each of them (pairing() asserts it); before, S had only 2, lam, d', delta, det.
        s.bad = sorted(set(bad))
        s.special_rational = False
        s.cache = {}
        s.X = tangent_decomposition(qform(lam, U), conic_point(lam, U))
        s.Y = tangent_decomposition(qform(lam, V), conic_point(lam, V))
        # the fibres of D_lam over p_X, p_Y and over the second points of n = 0 must have no rational point
        for T, F, G in ((s.X, U, V), (s.Y, V, U)):
            A = qform(lam, F); P = T['P']
            # the second point p' of the line n = 0 on the conic: P + t R with n(R) = 0, t = -2 B(P, R) / Q(R)
            Rs = [[T['n'][1], -T['n'][0], 0], [T['n'][2], 0, -T['n'][0]], [0, T['n'][2], -T['n'][1]]]
            R = next(c for c in Rs if any(c) and any(c[i] * P[j] - c[j] * P[i] for i in range(3) for j in range(3)) and qeval(A, c) != 0)
            Bpr = sum(A[i][j] * P[i] * R[j] for i in range(3) for j in range(3))
            t = -2 * Bpr / qeval(A, R)
            P2 = [P[i] + t * R[i] for i in range(3)]
            assert qeval(A, P2) == 0 and lin(T['n'], P2) == 0
            for (r, s_, W) in (P, P2):
                val = Fr(lam) * qv(G, Fr(r), Fr(s_))
                if val >= 0 and isqrt(val.numerator) ** 2 == val.numerator and isqrt(val.denominator) ** 2 == val.denominator:
                    s.special_rational = True      # a rational point of D_lam over p_X, p_Y or a second point of n = 0
    def S(s, dp):
        ps = {2} | set(s.bad)
        for q in (s.lam, dp, s.X['delta'], s.X['det'], s.Y['delta'], s.Y['det']):
            ps |= set(RI.factor_int(q))
        return sorted(ps)
    def g_class_local(s, dp, p, rng):
        """the set of values (dp, g(P))_p over sample points P of D_lam(Q_p) (each with the four sign choices)."""
        vals = set()
        if p not in s.cache:
            pts = local_points(s.lam, s.U, s.V, p, rng, NPTS)
            assert len(pts) == NPTS, ('no local point found', p, s.lam)
            s.cache[p] = [(r, s_, padic_sqrt(s.lam * qv(s.U, r, s_), p, NPREC), padic_sqrt(s.lam * qv(s.V, r, s_), p, NPREC)) for (r, s_) in pts]
        used = 0
        for (r, s_, Xr, Yr) in s.cache[p]:
            try:
                vv = set()
                for sx in (1, -1):
                    for sy in (1, -1):
                        lY = lin(s.Y['l'], (r, s_, sy * Yr)); lX = lin(s.X['l'], (r, s_, sx * Xr))
                        vv.add(_hilb_mod(dp, lY, p) * _hilb_mod(dp, lX, p))
            except AssertionError as ex:
                if 'precision' not in str(ex): raise
                continue                      # a point too close to a zero of l_X or l_Y: not used
            vals |= vv; used += 1
        assert used >= 1, ('no local point found', p, s.lam)
        return vals
    def g_sign_real(s, dp):
        vals = set()
        cands = _real_samples(s.U, s.V)
        n = 0
        for r in cands:
            X2, Y2 = s.lam * qv(s.U, r, 1), s.lam * qv(s.V, r, 1)
            if X2 > 0 and Y2 > 0:
                n += 1
                for sx in (1, -1):
                    for sy in (1, -1):
                        sg = _sign_lin(s.Y['l'], r, 1, Y2, sy) * _sign_lin(s.X['l'], r, 1, X2, sx)
                        vals.add(hilbert(dp, sg, 'inf'))
        assert n > 0, 'no real point found'
        return vals
    def pairing(s, dp, seed=0):
        rng = random.Random(seed + abs(s.lam) % 1000)
        tot = 1; per = {}
        for p in s.S(dp):
            v = s.g_class_local(dp, p, rng)
            assert len(v) == 1, ('the local evaluation is not constant', p, v)
            per[p] = min(v); tot *= per[p]
        v = s.g_sign_real(dp)
        assert len(v) == 1, ('the local evaluation is not constant at the real place', v)
        per['inf'] = min(v); tot *= per['inf']
        return tot, per

def local_points(lam, U, V, p, rng, want):
    """Coprime integers (r, s) with lam U(r, s) and lam V(r, s) non-zero squares in Q_p (tested exactly on the
    integers), so that (r : s : sqrt : sqrt) is a point of D_lam(Q_p).  Random pairs first; then, in the charts
    t = r/s in Z_p and t = s/r in pZ_p, a search of the discs t0 + p^k Z_p: the content p^c of each form on the disc is
    removed, and a disc is cut further only at residues where a form without its content vanishes mod p (at p = 2 also
    at every residue, to a fixed depth, since the class needs the value mod 8).  Floating point is not used; the
    search only proposes candidates, and every candidate is tested exactly."""
    out = []
    def good(r, s_):
        if r == 0 and s_ == 0: return False
        A, B = lam * qv(U, r, s_), lam * qv(V, r, s_)
        return A != 0 and B != 0 and _is_padic_square(A, p) and _is_padic_square(B, p)
    def add(r, s_):
        if gcd(r, s_) == 1 and (r, s_) not in out and good(r, s_): out.append((r, s_))
        return len(out) >= want
    for _ in range(3000):
        if add(rng.randrange(-p ** 4 - 50, p ** 4 + 50), rng.randrange(1, p ** 3 + 50)): return out
    if p == 2:
        # every disc of radius 2^-17 of both charts contains one of these t (rank_indep.solvable certifies points on
        # discs of radius at least 2^-17 at p = 2)
        for t in range(2 ** 17):
            if add(t, 1) or add(1, 2 * t): return out
    for chart in (0, 1):
        polys = [[lam * F[0], lam * F[1], lam * F[2]] if chart == 0 else [lam * F[2], lam * F[1], lam * F[0]] for F in (U, V)]
        stack = [(0, 1 if chart == 0 else p, 0)]
        nodes = 0
        while stack and nodes < 20000:
            nodes += 1
            t0, sc, depth = stack.pop()
            hs = []
            for g in polys:
                h = RI.ptrim(RI.psub_lin(g, t0, sc))
                k = min(RI.vp(c_, p) for c_ in h if c_)
                hs.append([c_ // p ** k for c_ in h])
            if p == 2 or p < 3000: cands = range(p)
            else:
                cands = set()
                for h in hs:
                    hh = RI.ptrim(h, p)
                    if len(hh) > 1: cands.update(RI.roots_mod(hh, p))
                cands.update(rng.randrange(p) for _ in range(30))
                cands = sorted(cands)
            for c_ in cands:
                t = t0 + sc * c_
                if add(*((t, 1) if chart == 0 else (1, t))): return out
                if depth < (22 if p == 2 else 12) and ((p == 2 and depth < 4) or any(RI.peval(h, c_) % p == 0 for h in hs)):
                    stack.append((t, sc * p, depth + 1))
    return out

def _is_padic_square(x, p):
    k = vp(x, p)
    if k % 2: return False
    u = x // p ** k
    return u % 8 == 1 if p == 2 else RI.legendre(u, p) == 1

def _hilb_mod(dp, x, p):
    """(dp, x)_p for x in Z_p known modulo p^NPREC (as an integer)."""
    y = x % p ** NPREC
    assert y != 0, 'p-adic precision (zero modulo p^N)'
    v = vp(y, p)
    assert v <= NPREC - 6, 'p-adic precision'
    u = (y // p ** v) % (8 if p == 2 else p)
    return hilbert(dp, p ** v * u, p)

def _sign_lin(c, r, s_, W2, sw):
    """sign of c0 r + c1 s + c2 sw sqrt(W2), exactly (W2 > 0 rational)."""
    A = c[0] * r + c[1] * s_; B = c[2] * sw
    if B == 0: return (A > 0) - (A < 0)
    if A == 0: return 1 if B > 0 else -1
    if (A > 0) == (B > 0): return 1 if A > 0 else -1
    if A * A > B * B * W2: return 1 if A > 0 else -1
    if A * A < B * B * W2: return 1 if B > 0 else -1
    raise AssertionError('a real zero of a tangent line')

def _real_samples(U, V):
    """rational test points between and beyond the real roots of U and V (floating point only places them)."""
    import numpy as np
    roots = sorted(set(float(x.real) for F in (U, V) for x in np.roots([F[2], F[1], F[0]]) if abs(x.imag) < 1e-9 * (1 + abs(x))))
    ends = [roots[0] - 10 - abs(roots[0])] + roots + [roots[-1] + 10 + abs(roots[-1])] if roots else [-1.0, 1.0]
    out = []
    for x0, x1 in zip(ends, ends[1:]):
        for j in range(1, 8): out.append(Fr(x0 + (x1 - x0) * j / 8).limit_denominator(10 ** 15))
    return out

def local_factor_one():
    """The statement used for the places outside S (no computation): let P = (r, s, X, Y) be a primitive integral
    point of D_lam not over p_X or the second point of n_X = 0 on C_X.  Then l_X(P) m_X(P) = -delta_X n_X(P)^2 with
    n_X(P) != 0.  If p does not divide 2 lam delta_X det_X and v_p(l_X(P)) is odd, then v_p(m_X(P)) is odd, so p
    divides l_X, m_X, n_X at (r, s, X), so p | r, s, X (det_X != 0 mod p); then p^2 | lam U(r, s) = X^2 and
    p^2 | lam V(r, s) = Y^2 give p | Y, against primitivity.  The same for Y.  So v_p(g(P)) is even, and for p not
    dividing d' either, (d', g(P))_p = 1."""
    return True

# ---------------------------------------------------------------- the curves
def certified(dd, AA, BB, p, loc):
    """C_dd: N^2 = dd M^4 + AA M^2 e^2 + (BB/dd) e^4 has a certified point over Q_p: from the local image that
    rank_indep.first_descent computed, or else by rank_indep.solvable to the depth 40 (needed at p = 2 for D1 of
    (631,71))."""
    if RI.qp_class(dd, p) in loc[p][1]: return True
    return RI.solvable([[BB * dd, 0, AA * dd * dd, 0, dd ** 3]], p, depth=40) is True

def sigma_curve(m, n, cv):
    M = SE.model(m, n, 2 if cv == 'E+' else -2)
    e, a, b = RI.shift_to_two_torsion(M['a2'], M['a4'], M['a6'])
    G = RI.evidence()[((m, n), cv)]
    return a, b, (G[0] - e, G[1])

def d1_curve(m, n):
    """The curve D1 of dm2.py: Y^2 = x^3 + c1 x^2 + c0 c2 x + c0^2 c3, with the point x = 0."""
    S = m * m + n * n; D2 = (m * m - n * n) ** 2
    c3 = 4 * S * S * m * m * n * n; c2 = -16 * S * m ** 3 * n ** 3; c1 = 4 * S * S * D2; c0 = -16 * S * m * n * D2
    a2, a4, a6 = c1, c0 * c2, c0 * c0 * c3
    y0 = c0 * 2 * S * m * n; assert y0 * y0 == a6
    e, a, b = RI.shift_to_two_torsion(a2, a4, a6)
    return a, b, (Fr(-e), Fr(y0))

def analyse(name, a, b, G, cross=True):
    """rank_indep's descent, then the pairing <y0, y1> on both sides."""
    say('%s:  y^2 = x(x^2 + a x + b),  a = %d,  b = %d' % (name, a, b))
    ap, bp = -2 * a, a * a - 4 * b
    r, rec = RI.rank_bound(a, b, [G], RI.find_points([0, ap, 0, bp, 0]), say=lambda s: None)
    say('  rank_indep.rank_bound: exact bound %d; open cosets E side %s, E\' side %s'
        % (r, rec['E']['open'], rec["E'"]['open']))
    R = RI.first_descent(a, b)
    res = dict(name=name, bound=r)
    if not rec['E']['open'] or not rec["E'"]['open']:
        say('  a side has no open coset: the pairing step is not needed'); return res
    assert len(rec['E']['open']) == 1 and len(rec["E'"]['open']) == 1
    d0 = rec['E']['open'][0][0]; dp0 = rec["E'"]['open'][0][0]
    # the two classes must be certified Selmer elements (C_d0 and C'_dp0 soluble at every bad prime and at R)
    for (dd, AA, BB, loc) in ((d0, a, b, R['locS']), (dp0, ap, bp, R['locSp'])):
        for p in R['bad']:
            assert certified(dd, AA, BB, p, loc), ('class not certified locally soluble', dd, p)
        assert RI.real_ok(dd, AA, BB)
    say('  d = %d (E side) and d\' = %d (E\' side) are certified in the Selmer groups (local points at every bad prime)' % (d0, dp0))
    vals = []
    for side, (dd, AA, BB, other) in (('E', (d0, a, b, dp0)), ("E'", (dp0, ap, bp, d0))):
        if side == "E'" and not cross: continue
        c = BB // dd
        P0 = RI.find_conic_point(dd, AA, c)
        surv, info = RI.second_descent(dd, AA, BB, P0, stop_at_first=False)
        U, V = info['U'], info['V']
        check_conic_forms(dd, AA, c, U, V, P0[2])
        say('  %s side: class %d, the curves D_lam that the second descent of rank_indep.py does not exclude: lam in %s' % (side, dd, sorted(surv, key=abs)))
        done = 0
        for lam in sorted(surv, key=abs):
            try:
                D = DLam(lam, U, V, bad=info['primes'])
            except AssertionError as ex:
                if 'no conic point' not in str(ex): raise
                say('    lam = %d: the conic W^2 = lam U or W^2 = lam V has no point over %s (a Hilbert symbol, exact), so D_lam is not everywhere locally soluble; skipped' % (lam, ex.args[0][1].split()[-1]))
                continue
            try:
                tot, per = D.pairing(other)
            except AssertionError as ex:
                if 'no local point found' not in str(ex): raise
                say('    lam = %d: no local point found at %s (not shown to be everywhere locally soluble; skipped)' % (lam, ex.args[0][1]))
                continue
            done += 1
            ctl, _ = D.pairing(b if side == "E'" else bp)       # the class of the rational 2-torsion point: must be +1
            assert ctl == 1, 'the class of the 2-torsion point pairs non-trivially: the formula is wrong'
            assert not D.special_rational, 'a rational point of D_lam over p_X, p_Y: D_lam is not empty'
            minus = [str(k) for k, v in per.items() if v == -1]
            say('    lam = %d: <%d, D_lam> = %+d  (places with factor -1: %s; |S| = %d; control with the torsion class: +1)'
                % (lam, other, tot, ', '.join(minus) if minus else 'none', len(per) - 1))
            vals.append(tot)
        assert done >= 1, 'no lift with local points everywhere'
    res['pairing'] = sorted(set(vals))
    say('  => <y0, y1> = %s' % ('0 (every value +1): the pairing on Sel^2 does not lower the bound' if set(vals) == {1}
                               else '1/2 (a value -1): both cosets are excluded'))
    return res

def control(m, n, cv, say=say):
    """On a curve where the second descent on the E' side excludes classes: the pairing of each class of the E' side
    with the lifts of the trivial class C_1 must reproduce the exclusions."""
    a, b, G = sigma_curve(m, n, cv)
    ap, bp = -2 * a, a * a - 4 * b
    r, rec = RI.rank_bound(a, b, [G], RI.find_points([0, ap, 0, bp, 0]), say=lambda s: None)
    R = RI.first_descent(a, b)
    excl = [c['d'] for c in rec["E'"]['excluded']]
    keep = [h for h in rec["E'"]['known'] if h not in (1,)] + [c[0] for c in rec["E'"]['open']]
    say('control %d,%d %s: E\' side, classes excluded by the second descent (one per coset): %s; kept: %s' % (m, n, cv, excl, keep))
    P0 = RI.find_conic_point(1, a, b)          # a point of N^2 = U^2 + a U V + b V^2 (a hint, checked by second_descent)
    surv, info = RI.second_descent(1, a, b, P0, stop_at_first=False)
    U, V = info['U'], info['V']
    check_conic_forms(1, a, b, U, V, P0[2])
    lams = sorted(surv, key=abs)
    say('  lifts of C_1: lam in %s' % lams)
    Ds = [DLam(lam, U, V, bad=info['primes']) for lam in lams]
    ok = True
    for dp in excl + keep:
        if not all(certified(dp, ap, bp, p, R['locSp']) for p in R['bad']):
            say('  %d: not certified locally soluble at every bad prime; skipped' % dp); continue
        v = [D.pairing(dp)[0] for D in Ds]
        expect_minus = dp in excl
        good = (-1 in v) == expect_minus
        ok &= good
        say('  d\' = %d: values %s; %s (%s)' % (dp, v, 'excluded' if expect_minus else 'kept', 'agrees' if good else 'DISAGREES'))
    assert ok, 'the pairing does not reproduce the second descent'
    return ok

def find_conic_point_any(d, a, c):
    """rank_indep.find_conic_point, but when PARI's point has N = 0 (the conic N^2 = d U^2 + a U V + c V^2 is
    tangent to... this happens when a^2 - 4 d c is a square, on the side of a curve with full 2-torsion), the second
    point of the conic on a line through it is taken.  A hint; second_descent checks it."""
    key = 'conic %d %d %d' % (d, a, c)
    if key in RI.HINTS.get('conic', {}): return tuple(RI.HINTS['conic'][key])
    import cypari2
    pari = cypari2.Pari()
    sol = pari.qfsolve(pari.matrix(3, 3, [2 * d, a, 0, a, 2 * c, 0, 0, 0, -2]))
    assert str(pari.type(sol)) == 't_COL', ('qfsolve', sol)
    P = [int(x) for x in sol]
    Q = lambda z: d * z[0] * z[0] + a * z[0] * z[1] + c * z[1] * z[1] - z[2] * z[2]
    k = 1
    while P[2] == 0:
        w = [1, k, 1 + k]; k += 1
        Bpw = 2 * d * P[0] * w[0] + a * (P[0] * w[1] + P[1] * w[0]) + 2 * c * P[1] * w[1] - 2 * P[2] * w[2]
        if Q(w) == 0: continue
        P = [Q(w) * P[i] - Bpw * w[i] for i in range(3)]
    g = gcd(gcd(P[0], P[1]), P[2]); P = [x // g for x in P]
    assert Q(P) == 0 and P[2] != 0
    RI.HINTS.setdefault('conic', {})[key] = P
    return tuple(P)

RATIONAL_LIFTS = []
def lift_pairing(dd, AA, BB, other, max_lifts=2):
    """<other, D_lam> for the lifts D_lam of the quartic C_dd of y^2 = x(x^2 + AA x + BB) (second descent), for up to
    max_lifts lifts with local points everywhere; they must agree.  Returns (value, lams used)."""
    c = BB // dd
    P0 = find_conic_point_any(dd, AA, c)
    surv, info = RI.second_descent(dd, AA, BB, P0, stop_at_first=False)
    U, V = info['U'], info['V']
    check_conic_forms(dd, AA, c, U, V, P0[2])
    vals = {}; tried = []
    for lam in sorted(surv, key=abs):
        if len(vals) >= max_lifts: break
        try: D = DLam(lam, U, V)
        except AssertionError as ex:
            if 'no conic point' not in str(ex): raise
            tried.append(lam); continue
        try: tot, per = D.pairing(other)
        except AssertionError as ex:
            if 'no local point found' not in str(ex): raise
            tried.append(lam); continue
        if D.special_rational:
            # D_lam has a rational point (over p_X or p_Y): then the value must be +1 (Hilbert reciprocity)
            assert tot == 1, ('a curve with a rational point pairs to -1: the formula is wrong', dd, lam)
            RATIONAL_LIFTS.append((dd, lam))
        vals[lam] = tot
    assert vals, ('no lift with local points everywhere', dd)
    assert len(set(vals.values())) == 1, ('two lifts give different values', dd, vals)
    return list(vals.values())[0], sorted(vals)

def analyse_general(name, a, b, ptsE, ptsEp):
    """As analyse, for any number of open cosets: the matrix of <y, lift(x)> (y: one class of each open coset of the
    E' side, x: one of each open coset of the E side); its rank rho over F_2 lowers the bound by 2 rho."""
    say('%s:  y^2 = x(x^2 + a x + b),  a = %d,  b = %d' % (name, a, b))
    ap, bp = -2 * a, a * a - 4 * b
    r, rec = RI.rank_bound(a, b, ptsE, ptsEp, say=lambda s: None)
    X = [c[0] for c in rec['E']['open']]; Y = [c[0] for c in rec["E'"]['open']]
    say('  rank_indep.rank_bound: exact bound %d; E side: Selmer %d, known %d, open cosets %d; E\' side: Selmer %d, known %d, open cosets %d'
        % (r, len(rec['E']['selmer']), len(rec['E']['known']), len(X), len(rec["E'"]['selmer']), len(rec["E'"]['known']), len(Y)))
    R = RI.first_descent(a, b)
    Ycert = [y for y in Y if all(certified(y, ap, bp, p, R['locSp']) for p in R['bad']) and RI.real_ok(y, ap, bp)]
    say('  E\' side classes certified locally soluble at every bad prime: %s of %s' % (Ycert, Y))
    M = []
    for y in Ycert:
        row = []
        for x in X:
            v, lams = lift_pairing(x, a, b, y)
            row.append(0 if v == 1 else 1)
            say('    <%d, lift of C_%d> = %+d  (lifts used: lam = %s)' % (y, x, v, lams))
        M.append(row)
    # rank over F_2
    rows = [sum(bit << j for j, bit in enumerate(rw)) for rw in M]
    rho = len(RS_span(rows))
    newb = r - 2 * rho
    say('  => rank of the pairing matrix %d: rank <= %d' % (rho, newb))
    return dict(name=name, bound=r, rho=rho, new=newb)

def RS_span(vecs):
    basis = []
    for v in vecs:
        for bb in basis: v = min(v, v ^ bb)
        if v: basis.append(v)
    return basis

def e0l_curve(pair, e_index):
    """The curve E0L of dm4_log.jsonl (y^2 = x(x - e1)(x - e2), full 2-torsion), its points G_T; returned as the
    curve B = E0L / <(e, 0)> with one rational 2-torsion point (shifted), so that rank_indep's descent applies, with
    the images of the points on B and on B' = E0L (scaled by 4)."""
    rec = {tuple(json.loads(l)['pair']): json.loads(l) for l in open(os.path.join(here, 'dm4_log.jsonl'))}[pair]
    a2, a4 = int(rec['model'][1]), int(rec['model'][3]); assert int(rec['model'][4]) == 0
    roots = [0] + [int(x) for x in rec['roots']]
    for e in roots: assert e ** 3 + a2 * e * e + a4 * e == 0
    G = [(Fr(x), Fr(y)) for x, y in rec['G_T']]
    for P in G: assert P[1] ** 2 == P[0] ** 3 + a2 * P[0] ** 2 + a4 * P[0], 'a point of dm4_log.jsonl is not on the curve'
    e = roots[e_index]
    a, b = 3 * e + a2, 3 * e * e + 2 * a2 * e + a4          # E0L shifted: y^2 = x(x^2 + a x + b), a^2 - 4b a square
    Gs = [(P[0] - e, P[1]) for P in G]
    assert isqrt(a * a - 4 * b) ** 2 == a * a - 4 * b
    A, Bc = -2 * a, a * a - 4 * b                            # B: Y^2 = X(X^2 + A X + Bc), B' = E0L scaled by 4
    phi = lambda P: (P[1] ** 2 / P[0] ** 2, P[1] * (b - P[0] ** 2) / P[0] ** 2)
    ptsB = [phi(P) for P in Gs if P[0] != 0]
    ptsBp = [(4 * P[0], 8 * P[1]) for P in Gs]
    return A, Bc, ptsB, ptsBp, G, (a2, a4), roots

# ---------------------------------------------------------------- full 2-torsion: the Cassels-Tate pairing on Sel^2
def sq_cls(x, p):
    """class of a non-zero rational in Q_p*/Q_p*^2 (or R*/R*^2) as a bit vector (an integer)."""
    x = Fr(x)
    if p == 'inf': return 1 if x < 0 else 0
    n, d = x.numerator, x.denominator
    v = vp(n, p) - vp(d, p)
    u = (n // p ** vp(n, p)) * (d // p ** vp(d, p))
    if p == 2:
        u %= 8
        return (v % 2) | ((1 if u in (3, 7) else 0) << 1) | ((1 if u in (5, 7) else 0) << 2)
    return (v % 2) | ((0 if RI.legendre(u, p) == 1 else 1) << 1)

def f2_span(vecs):
    B = []
    for v in vecs:
        for b in B: v = min(v, v ^ b)
        if v: B.append(v)
    return B

def full_selmer2(es, rng):
    """Sel^2 of y^2 = (x - e1)(x - e2)(x - e3) (e_i distinct integers) by the complete 2-descent: classes (d1, d2)
    of (x - e1, x - e2) in the group generated by -1 and the primes of 2 (e1 - e2)(e1 - e3)(e2 - e3), locally in the
    image of E(Q_v) at every bad place.  The local image has order 4 (p odd), 8 (p = 2), 2 (R) [|E(K)/2E(K)| =
    |E(K)[2]| |2|_K^-1]; it is found from the 2-torsion and from explicit local points (x rational with
    (x - e1)(x - e2)(x - e3) a square in Q_v, tested exactly) until it has that order."""
    e1, e2, e3 = es
    disc = (e1 - e2) * (e1 - e3) * (e2 - e3)
    bad = sorted(set([2]) | set(RI.factor_int(disc)))
    places = ['inf'] + bad
    W = {p: (1 if p == 'inf' else (3 if p == 2 else 2)) for p in places}
    def kum(x):
        if x == e1: return ((e1 - e2) * (e1 - e3), e1 - e2)
        if x == e2: return (e2 - e1, (e2 - e1) * (e2 - e3))
        return (x - e1, x - e2)
    def enc(pair, p): return sq_cls(pair[0], p) | (sq_cls(pair[1], p) << W[p])
    loc = {}
    H = max(abs(e) for e in es) + 1
    for p in places:
        target = {'inf': 1}.get(p, 3 if p == 2 else 2)
        B = f2_span([enc(kum(x), p) for x in (e1, e2, e3)])
        tries = 0
        while len(B) < target and tries < 300000:
            tries += 1
            if p == 'inf': x = Fr(rng.randrange(-4 * H, 4 * H), rng.randrange(1, 60))
            else:
                e = es[rng.randrange(3)]
                x = e + Fr(rng.randrange(1, p ** 6)) * Fr(p) ** rng.randrange(-4, 12) * rng.choice((1, -1))
            y2 = (x - e1) * (x - e2) * (x - e3)
            if y2 == 0: continue
            if p == 'inf':
                if y2 < 0: continue
            elif sq_cls(y2, p) != 0: continue
            B = f2_span(B + [enc(kum(x), p)])
        assert len(B) == target, ('local image not found', p)
        loc[p] = B
    gens = [-1] + bad
    glob = [(g, 1) for g in gens] + [(1, g) for g in gens]
    def red(v, B):
        for b in sorted(B, reverse=True): v = min(v, v ^ b)
        return v
    # echelon bases so that red() is a canonical reduction
    for p in places:
        E = []
        for v in loc[p]:
            for b in E: v = min(v, v ^ b)
            if v: E = [min(b, b ^ v) for b in E] + [v]
        loc[p] = E
    sh = {}; tot = 0
    for p in places: sh[p] = tot; tot += 2 * W[p]
    rows = []
    for g in glob:
        r = 0
        for p in places: r |= red(enc(g, p), loc[p]) << sh[p]
        rows.append(r)
    piv = []; kern = []
    for i, r in enumerate(rows):
        c = 1 << i
        for pv, pc in piv:
            if r ^ pv < r: r ^= pv; c ^= pc
        if r == 0: kern.append(c)
        else: piv.append((r, c)); piv.sort(reverse=True)
    sel = []
    for c in f2_span(kern):
        d1 = d2 = 1
        for i in range(len(glob)):
            if c >> i & 1: d1 *= glob[i][0]; d2 *= glob[i][1]
        sel.append((RI.sqfree_part(d1), RI.sqfree_part(d2)))
    return bad, loc, sel, kum, enc, places

def gen_conic_point(A):
    """A rational point of z^T A z = 0 (A a symmetric 3x3 matrix of integers) from PARI's qfsolve (a hint; checked)."""
    key = 'gen %s' % [[int(x) for x in row] for row in A]
    if key in HINTS['conic']: P = [int(x) for x in HINTS['conic'][key]]
    else:
        import cypari2
        pari = cypari2.Pari()
        sol = pari.qfsolve(pari.matrix(3, 3, [int(2 * A[i][j]) for i in range(3) for j in range(3)]))
        assert str(pari.type(sol)) == 't_COL', ('no conic point', sol)
        P = [int(x) for x in sol]; g = gcd(gcd(P[0], P[1]), P[2]); P = [x // g for x in P]
        HINTS['conic'][key] = P
    assert qeval(A, P) == 0 and any(P), 'not a point of the conic'
    return P

def _other_point(C, P):
    """Another rational point of the conic z^T C z = 0: the second intersection with the line through P and
    (1, 2, 3) (or another fixed point off the conic), as a primitive integer vector."""
    for R in ([1, 2, 3], [3, 1, 2], [2, 3, 5], [1, -1, 4]):
        if qeval(C, R) == 0: continue
        t = -2 * sum(C[i][j] * P[i] * R[j] for i in range(3) for j in range(3)) / qeval(C, R)
        if t == 0: continue
        Q = [P[i] + t * R[i] for i in range(3)]
        den = 1
        for x in Q: den = den * Fr(x).denominator // gcd(den, Fr(x).denominator)
        Q = [int(x * den) for x in Q]; g = 0
        for x in Q: g = gcd(g, x)
        Q = [x // g for x in Q]
        assert qeval(C, Q) == 0
        if any(Q[i] * P[j] != Q[j] * P[i] for i in range(3) for j in range(3)): return Q
    raise AssertionError('no second point')

class DFull:
    """The 2-covering x = e1 + a1 u^2 = e2 + a2 v^2 = e3 + a1 a2 w^2 (z = 1) of y^2 = (x - e1)(x - e2)(x - e3),
    with the three conics C_u (v, w, z), C_v (u, w, z), C_w (u, v, z) obtained by forgetting one coordinate, their
    tangent lines l_u, l_v, l_w at rational points (div l = 2 F, F the fibre class), and the classes
    F_v - F_w ~ T1, F_u - F_w ~ T2 (translations (u,v,w) -> (u,-v,-w), (-u,v,-w) fix x - e1, resp. x - e2)."""
    def __init__(s, es, a, alt=False):
        s.es = es; e1, e2, e3 = es
        a1, a2 = a; a3 = a1 * a2
        s.a = (a1, a2, a3)
        F = Fr
        # C_u: a2 v^2 - a3 w^2 + (e2 - e3) z^2 = 0 ; C_v: a1 u^2 - a3 w^2 + (e1 - e3) z^2 ; C_w: a1 u^2 - a2 v^2 + (e1 - e2) z^2
        s.Cu = [[F(a2), F(0), F(0)], [F(0), F(-a3), F(0)], [F(0), F(0), F(e2 - e3)]]
        s.Cv = [[F(a1), F(0), F(0)], [F(0), F(-a3), F(0)], [F(0), F(0), F(e1 - e3)]]
        s.Cw = [[F(a1), F(0), F(0)], [F(0), F(-a2), F(0)], [F(0), F(0), F(e1 - e2)]]
        s.T = {}
        for nm, C in (('u', s.Cu), ('v', s.Cv), ('w', s.Cw)):
            P = gen_conic_point(C)
            if alt: P = _other_point(C, P)       # a second rational point of the conic: the tangent there instead
            s.T[nm] = tangent_decomposition(C, P)
        s.cache = {}
        # rational points of D_a over the tangent points and over the second points of n = 0 (where g is 0, infinity,
        # or the argument for the places outside S does not apply)
        s.special_rational = False
        idx = {'u': (1, 2), 'v': (0, 2), 'w': (0, 1)}
        for nm, C in (('u', s.Cu), ('v', s.Cv), ('w', s.Cw)):
            T = s.T[nm]; P = T['P']
            Rs = [[T['n'][1], -T['n'][0], 0], [T['n'][2], 0, -T['n'][0]], [0, T['n'][2], -T['n'][1]]]
            R = next(c for c in Rs if any(c) and any(c[i] * P[j] - c[j] * P[i] for i in range(3) for j in range(3)) and qeval(C, c) != 0)
            t = -2 * sum(C[i][j] * P[i] * R[j] for i in range(3) for j in range(3)) / qeval(C, R)
            P2 = [P[i] + t * R[i] for i in range(3)]
            assert qeval(C, P2) == 0
            for Q in (P, P2):
                if Q[2] == 0: continue
                # x from the two known coordinates (coordinate k of the conic is coordinate idx[nm][0] of (u, v, w))
                k = idx[nm][0]; e = (e1, e2, e1 * 0 + e3)[k]; ak = s.a[k]
                x = e + ak * Fr(Q[0]) ** 2 / Fr(Q[2]) ** 2
                m_ = {'u': 0, 'v': 1, 'w': 2}[nm]
                val = (x - s.es[m_]) / s.a[m_]
                if val == 0 or (val > 0 and isqrt(val.numerator) ** 2 == val.numerator and isqrt(val.denominator) ** 2 == val.denominator):
                    s.special_rational = True
    def coords(s, x, p):
        """(U, V, W, Z) integers (p-adic approximations) of a point above x, scaled by a1 a2 a3: U = sqrt((x - e1) a1) a2 a3."""
        a1, a2, a3 = s.a; e1, e2, e3 = s.es
        x = Fr(x)
        out = []
        den = x.denominator  # x = n/den with den a square times... use x - e = (n - e den)/den; multiply by den^2
        for e, ai, other in ((e1, a1, a2 * a3), (e2, a2, a1 * a3), (e3, a3, a1 * a2)):
            N = (x.numerator - e * den) * den * ai           # (x - e) ai den^2
            if p == 'inf': out.append((N, other))
            else: out.append(padic_sqrt(N, p, NPREC) * other)
        return out, a1 * a2 * a3 * den
    def ok_x(s, x, p):
        a1, a2, a3 = s.a
        vals = [(Fr(x) - e) * ai for e, ai in zip(s.es, (a1, a2, a3))]
        if any(v == 0 for v in vals): return False
        if p == 'inf': return all(v > 0 for v in vals)
        return all(sq_cls(v, p) == 0 for v in vals)
    def lvals(s, x, p, signs):
        (U, V, W), Z = s.coords(x, p)
        U, V, W = signs[0] * U, signs[1] * V, signs[2] * W
        lu = lin(s.T['u']['l'], (V, W, Z)); lv = lin(s.T['v']['l'], (U, W, Z)); lw = lin(s.T['w']['l'], (U, V, Z))
        return lu, lv, lw
    def local_xs(s, p, rng, want=NPTS):
        if p in s.cache: return s.cache[p]
        out = []; H = max(abs(e) for e in s.es) + 1
        for _ in range(200000):
            if len(out) >= want: break
            if p == 'inf':
                x = Fr(rng.randrange(-8 * H, 8 * H), rng.randrange(1, 60))
            else:
                e = s.es[rng.randrange(3)]
                x = e + Fr(rng.randrange(1, p ** 5)) * Fr(p) ** rng.randrange(-6, 14) * rng.choice((1, -1))
                x = Fr(x.numerator, x.denominator)
                if x.denominator % p == 0 and vp(x.denominator, p) % 2: continue
            if s.ok_x(x, p) and x not in out: out.append(x)
        assert len(out) == want, ('no local point found', p)
        s.cache[p] = out
        return out
    def S(s, bs):
        ps = {2}
        e1, e2, e3 = s.es
        for q in list(s.a) + [e1 - e2, e1 - e3, e2 - e3] + [s.T[k][kk] for k in 'uvw' for kk in ('delta', 'det')] + list(bs):
            ps |= set(RI.factor_int(q))
        return sorted(ps)
    def pairing(s, b, seed=0):
        """<b, D_a> = prod_v (b2, l_v/l_w (P_v))_v (b1, l_u/l_w (P_v))_v for b = (b1, b2) in Sel^2."""
        b1, b2 = b
        rng = random.Random(seed)
        tot = 1; per = {}
        for p in s.S(b) + ['inf']:
            vals = set()
            for x in s.local_xs(p, rng):
                for signs in ((1, 1, 1), (1, -1, 1), (-1, 1, 1), (1, 1, -1), (-1, -1, -1)):
                    if p == 'inf':
                        su, sv, sw = _dfull_real(s, x, signs)
                        h = hilbert(b2, sv * sw, 'inf') * hilbert(b1, su * sw, 'inf')
                    else:
                        lu, lv, lw = s.lvals(x, p, signs)
                        try: h = _hilb_mod(b2, lv, p) * _hilb_mod(b2, lw, p) * _hilb_mod(b1, lu, p) * _hilb_mod(b1, lw, p)
                        except AssertionError as ex:
                            if 'precision' in str(ex): continue
                            raise
                    vals.add(h)
            assert len(vals) == 1, ('the local evaluation is not constant', p, vals)
            per[p] = min(vals); tot *= per[p]
        return tot, per

def _sgn(x): return (x > 0) - (x < 0)

def _sign2(r, q, C):
    """sign of r + q sqrt(C), C >= 0 rational, exactly."""
    if q == 0 or C == 0: return _sgn(r)
    if r == 0: return _sgn(q)
    if _sgn(r) == _sgn(q): return _sgn(r)
    d = r * r - q * q * C
    return _sgn(r) if d > 0 else (_sgn(q) if d < 0 else 0)

def _sign3(c0, c1, A, c2, B):
    """sign of c0 + c1 sqrt(A) + c2 sqrt(B), A, B >= 0 rationals, exactly."""
    st = _sign2(0, c1, A) if c2 == 0 else (_sign2(0, c2, B) if c1 == 0 else None)
    if st is None:
        # t = c1 sqrt(A) + c2 sqrt(B): sign of t = sign of c1 sqrt(A) + c2 sqrt(B) = _sign2 on t / sqrt(B)... compare squares
        if _sgn(c1) == _sgn(c2): st = _sgn(c1)
        else:
            dd = c1 * c1 * A - c2 * c2 * B
            st = _sgn(c1) if dd > 0 else (_sgn(c2) if dd < 0 else 0)
    if c0 == 0: return st
    if st == 0 or st == _sgn(c0): return _sgn(c0) if st == 0 else st
    # opposite signs: compare t^2 = c1^2 A + c2^2 B + 2 c1 c2 sqrt(AB) with c0^2
    sd = _sign2(c1 * c1 * A + c2 * c2 * B - c0 * c0, 2 * c1 * c2, A * B)
    return st if sd > 0 else (_sgn(c0) if sd < 0 else 0)

def _dfull_real(D, x, signs):
    """signs of l_u, l_v, l_w at the real point above x with the given signs of (U, V, W)."""
    a1, a2, a3 = D.a; e1, e2, e3 = D.es
    x = Fr(x)
    Ns = [(x - e) * ai for e, ai in zip(D.es, (a1, a2, a3))]
    others = [a2 * a3, a1 * a3, a1 * a2]
    Z = a1 * a2 * a3
    # coordinate k = signs[k] * others[k] * sqrt(Ns[k]); Z rational
    def sgn_form(c, idx):
        # c: coefficients on (coord idx[0], coord idx[1], Z)
        c1 = c[0] * signs[idx[0]] * others[idx[0]]; c2 = c[1] * signs[idx[1]] * others[idx[1]]
        r = _sign3(c[2] * Z, c1, Ns[idx[0]], c2, Ns[idx[1]])
        assert r != 0, 'a real zero of a tangent line'
        return r
    return sgn_form(D.T['u']['l'], (1, 2)), sgn_form(D.T['v']['l'], (0, 2)), sgn_form(D.T['w']['l'], (0, 1))

def e0l_rank(pair, say=say):
    """The curve E0L of the four-point tier (dm4_log.jsonl; full 2-torsion): Sel^2 by the complete 2-descent, the
    Cassels-Tate pairing on it by the formula of DFull.pairing, and the images of the three points G_T and of the
    2-torsion.  rank <= dim Sel^2 - 2 - rank(pairing); rank >= the rank of the images of the points modulo torsion."""
    rec = {tuple(json.loads(l)['pair']): json.loads(l) for l in open(os.path.join(here, 'dm4_log.jsonl'))}[pair]
    a2, a4 = int(rec['model'][1]), int(rec['model'][3]); assert int(rec['model'][4]) == 0 and int(rec['model'][0]) == 0
    es = tuple([0] + [int(x) for x in rec['roots']])
    for e in es: assert e ** 3 + a2 * e * e + a4 * e == 0
    G = [(Fr(x), Fr(y)) for x, y in rec['G_T']]
    for P in G: assert P[1] ** 2 == P[0] ** 3 + a2 * P[0] ** 2 + a4 * P[0], 'a point of dm4_log.jsonl is not on the curve'
    say('E0L of %d,%d (dm4_log.jsonl): y^2 = x(x - %d)(x - %d)' % (pair[0], pair[1], es[1], es[2]))
    bad, loc, sel, kum, enc, places = full_selmer2(es, random.Random(1))
    say('  complete 2-descent: dim Sel^2 = %d (bad primes %s); basis %s' % (len(sel), bad, sel))
    gens = [-1] + bad
    def gvec(pr):
        v = 0
        for j, d in enumerate(pr):
            d = RI.sqfree_part(Fr(d).numerator * Fr(d).denominator)
            for i, g in enumerate(gens):
                if (g == -1 and d < 0) or (g != -1 and d % g == 0): v |= 1 << (i + j * len(gens))
            rest = abs(d)
            for g in bad:
                if rest % g == 0: rest //= g
            assert rest == 1, 'a class outside the group of the bad primes'
        return v
    selv = [gvec(a) for a in sel]
    assert len(f2_span(selv)) == len(sel)
    imgs = [kum(P[0]) for P in G] + [kum(es[0]), kum(es[1])]
    imgv = [gvec(c) for c in imgs]
    assert len(f2_span(selv + imgv)) == len(sel), 'a point image outside Sel^2: the descent is wrong'
    ri = len(f2_span(imgv))
    say('  the images of the three points and of two points of order 2 span a subgroup of dimension %d of Sel^2' % ri)
    Ds = [DFull(es, a) for a in sel]
    for i, D in enumerate(Ds):
        if D.special_rational: say('    D_%s has a rational point over a tangent point (so it must pair trivially)' % (sel[i],))
    M = []
    for i, a in enumerate(sel):
        row = [Ds[i].pairing(b)[0] for b in sel]
        M.append(row)
        say('    <., D_%s>: %s' % (a, ' '.join('%+d' % v for v in row)))
    for i in range(len(sel)):
        if Ds[i].special_rational: assert all(v == 1 for v in M[i]), 'a curve with a rational point pairs non-trivially'
        assert M[i][i] == 1, 'the pairing is not alternating'
        for j in range(len(sel)): assert M[i][j] == M[j][i], 'the pairing is not symmetric'
    # the images of points pair trivially with every D_a
    pv = [Ds[i].pairing((RI.sqfree_part(Fr(c[0]).numerator * Fr(c[0]).denominator), RI.sqfree_part(Fr(c[1]).numerator * Fr(c[1]).denominator)))[0]
          for i in range(len(sel)) for c in imgs]
    assert all(v == 1 for v in pv), 'the image of a rational point pairs non-trivially'
    # the same matrix with the tangent lines at other rational points of the nine conics (g changes by a square
    # times a constant that is a norm locally; the pairing must not change)
    M2 = [[DFull(es, a, alt=True).pairing(b)[0] for b in sel] for a in sel]
    assert M2 == M, 'the pairing depends on the choice of the tangent points'
    say('  the same matrix with the tangent lines taken at other rational points of the conics')
    rho = len(f2_span([sum((0 if v == 1 else 1) << j for j, v in enumerate(rw)) for rw in M]))
    up = len(sel) - 2 - rho; lo = ri - 2
    say('  pairing: symmetric, alternating, the point images in its kernel; rank of the pairing %d' % rho)
    say('  => %d <= rank E0L(Q) <= %d%s' % (lo, up, ' (exact)' if lo == up else ''))
    return dict(pair=pair, sel=len(sel), rho=rho, lo=lo, up=up)

CURVES = [('760,39 E+', ('s', 760, 39, 'E+')), ('945,248 E+', ('s', 945, 248, 'E+')), ('72,65 E+', ('s', 72, 65, 'E+')),
          ('920,129 E+', ('s', 920, 129, 'E+')), ('631,71 D1', ('d', 631, 71))]

def curve_of(spec):
    if spec[0] == 's': return sigma_curve(*spec[1:])
    return d1_curve(*spec[1:])

def negative():
    res = []
    def run(nm, fn):
        try:
            fn(); say('  NOT REFUSED: ' + nm); res.append(False)
        except AssertionError as ex:
            say('  refused: %s (%s)' % (nm, str(ex)[:90])); res.append(True)
    a, b, G = sigma_curve(72, 65, 'E+')
    ap, bp = -2 * a, a * a - 4 * b
    d0 = 2; P0 = RI.find_conic_point(d0, a, b // d0)
    surv, info = RI.second_descent(d0, a, b, P0, stop_at_first=False)
    U, V = info['U'], info['V']; lam = sorted(surv, key=abs)[0]
    def f1():
        key = '%d %s' % (lam, U); saved = HINTS['conic'].get(key)
        P = list(conic_point(lam, U)); HINTS['conic'][key] = [P[0] + 1, P[1], P[2]]
        try: DLam(lam, U, V)
        finally:
            if saved is None: HINTS['conic'].pop(key, None)
            else: HINTS['conic'][key] = saved
    run('a forged point of the conic C_X', f1)
    def f2():
        U2 = list(U); U2[1] += 2
        check_conic_forms(d0, a, b // d0, U2, V, P0[2])
    run('a mutated coefficient of U (a forged homogeneous space D_lam)', f2)
    def f3():
        A = qform(lam, U); P = conic_point(lam, U)
        T = tangent_decomposition(A, P)
        # a line through P that is not the tangent: the decomposition must fail
        A2 = [row[:] for row in A]; A2[2][2] += 1
        tangent_decomposition(A2, P)
    run('a tangent decomposition at a point that is not on the (mutated) conic', f3)
    def f4():
        # a class d' that is not in Sel^phi(E): the local evaluation is not constant somewhere, or a local test fails
        D = DLam(lam, U, V)
        R = RI.first_descent(a, b)
        bad = None
        for q in [3, 5, 7, 13, -3, 6, 10]:
            if any(RI.qp_class(q, p) not in R['locSp'][p][0] for p in R['bad']): bad = q; break
        assert bad is not None
        D.pairing(bad)
    run('pairing with a class outside the Selmer group (evaluation must not be constant)', f4)
    def f5():
        # a lam that is not everywhere locally soluble: no local point at some place
        lam2 = next(l for l in (3, 5, 6, 10, 15, 21, 30) if l not in surv)
        DLam(lam2, U, V).pairing(2)
    run('a curve D_lam with lam outside the survivors (no local point somewhere)', f5)
    def f6():
        # a line through p_X that is not the tangent: g has a divisor that is not even, the evaluation must vary
        D = DLam(lam, U, V)
        D.X = dict(D.X); D.X['l'] = list(D.X['n'])
        D.pairing(2)
    run('a non-tangent line in place of l_X (the evaluation must not be constant)', f6)
    def f7():
        # full 2-torsion: a class b outside Sel^2 must give a non-constant evaluation somewhere
        es = (0, 1095776342882, 1883582727842)
        D = DFull(es, (1, -1))
        D.pairing((1, 3))
    run('full 2-torsion: pairing with a class outside Sel^2 (evaluation must not be constant)', f7)
    return all(res), res

if __name__ == '__main__':
    if os.path.exists(RI.HINTS_FILE):
        for k, v in json.load(open(RI.HINTS_FILE)).items(): RI.HINTS[k] = v
    if os.path.exists(HINTS_FILE):
        H = json.load(open(HINTS_FILE))
        HINTS['conic'].update(H.get('conic', {}))
        for k, v in H.get('rank_indep', {}).items():
            for kk, vv in v.items(): RI.HINTS.setdefault(k, {}).setdefault(kk, vv)
    def save_hints():
        if WRITE:
            base = json.load(open(RI.HINTS_FILE)) if os.path.exists(RI.HINTS_FILE) else {}
            extra = {k: {kk: vv for kk, vv in RI.HINTS.get(k, {}).items() if kk not in base.get(k, {})} for k in ('factor', 'conic', 'points')}
            json.dump({'conic': HINTS['conic'], 'rank_indep': extra}, open(HINTS_FILE, 'w'), indent=0, sort_keys=True)
    if '--negative' in sys.argv:
        ok, res = negative()
        say('negative tests: %d of %d refused' % (sum(res), len(res)))
        save_hints()
        if WRITE: open(os.path.join(here, 'rank_second_negative.txt'), 'w').write('\n'.join(OUT) + '\n')
        sys.exit(0 if ok else 1)
    say('rank_second.py: the Cassels-Tate pairing <y0, y1> on Sel^2 for the curves whose rank bound 1 rests on Magma alone')
    nh = hilbert_selftest(); say('Hilbert symbol formulas checked against a search at p = 2, 3, 5 for %d pairs (a, b)' % nh)
    control(945, 248, 'E+')
    results = []
    for name, spec in CURVES:
        a, b, G = curve_of(spec)
        results.append(analyse(name, a, b, G))
    say('summary: ' + '; '.join('%s: exact bound %d, <y0, y1> values %s' % (r['name'], r['bound'], r.get('pairing')) for r in results))
    say('')
    say('The curves E0L of the four-point tier for (72,65) = (137,7) (one curve) and (920,129): full 2-torsion')
    e0 = [e0l_rank(pr) for pr in ((72, 65), (920, 129))]
    say('summary E0L: ' + '; '.join('%d,%d: dim Sel^2 %d, pairing rank %d, %d <= rank <= %d' % (r['pair'][0], r['pair'][1], r['sel'], r['rho'], r['lo'], r['up']) for r in e0))
    save_hints()
    if WRITE: open(os.path.join(here, 'rank_second.txt'), 'w').write('\n'.join(OUT) + '\n')
