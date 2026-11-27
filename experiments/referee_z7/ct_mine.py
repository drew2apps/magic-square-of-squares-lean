#!/usr/bin/env python3
"""Review z7 of roadmap H (H7) and of scripts/second_descent/rank0_ct.py: an independent Cassels-Tate pairing on Sel^2.
Written by an AI model (review agent z7, 27 November 2026); not read by a person.  Experiment, not counted.

Written from the reviewer's own derivation, with the APPROACH of reviews/review_t2.md (PARI p-adics, pari.hilbert,
conic points moved along random lines), not its code and not scripts/g1_full/rank_second.py's code.  For
E: y^2 = (x - e1)(x - e2)(x - e3) and a in Sel^2 (my2desc.py), the 2-covering
    D_a : x - e1 = a1 u^2, x - e2 = a2 v^2, x - e3 = a3 w^2   (a3 = a1 a2),
its conics C_u (v,w,z): a2 v^2 - a3 w^2 - (e3-e2) z^2 = 0, C_v (u,w,z): a1 u^2 - a3 w^2 - (e3-e1) z^2 = 0,
C_w (u,v,z): a1 u^2 - a2 v^2 - (e2-e1) z^2 = 0, and tangent lines l_u, l_v, l_w at rational points.
Derivation: z_i o (translation by T) = e2(T, T_i) z_i, so T1 = (u,-v,-w), T2 = (-u,v,-w); the reflections give
F_u - F_w ~ T2, F_v - F_w ~ T1, div(l_i / l_j) = 2 (F_i - F_j); b = b2 (x) T1 + b1 (x) T2 for Kummer coordinates
(b1, b2) = (e2(b,T1), e2(b,T2)); so
    <b, D_a> = prod_v (b2, l_v/l_w (P_v))_v (b1, l_u/l_w (P_v))_v,   P_v in D_a(Q_v).
Outside {inf, 2, primes of a, b, disc} the factor is 1 ((b, g) is unramified on the smooth Z_p-model); the small
primes 3, 5, 7, 11, 13 are evaluated anyway and asserted to give +1.  Each local value is computed at 4 local points
with all 8 sign choices and asserted constant.  rank E(Q) <= dim Sel^2 - 2 - rank of the pairing matrix.
Rests on (quoted): the Cassels-Tate pairing and its Brauer-Manin description (Cassels, J. reine angew. Math. 494
(1998); Skorobogatov, Torsors and rational points, 6.2).
PARI (cypari2): qfsolve finds the conic points (each checked exactly), p-adic square roots of exact rationals to
precision p^80 and hilbert() on exact rationals built from them.  The real place is evaluated with mpmath at 120
digits with a guard of 10^-80 relative to the terms: a cross-check, not a certificate.

    python3 ct_mine.py        self-test on y^2 = x^3 - n^2 x (prints ct_mine.txt)
"""
import sys, os, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)

from fractions import Fraction as Fr
from math import gcd
import cypari2, mpmath
pari = cypari2.Pari()
from my2desc import Curve, loc_class, is_sq_local, span
mpmath.mp.dps = 120
PREC = 80

def qe(A, z): return sum(A[i][j] * z[i] * z[j] for i in range(3) for j in range(3))

def conic_point(A, rng, moves=2):
    M = pari.matrix(3, 3, [int(A[i][j]) for i in range(3) for j in range(3)])
    sol = pari.qfsolve(M)
    assert str(pari.type(sol)) == 't_COL', 'no conic point'
    P = [Fr(int(x)) for x in sol]
    assert qe(A, P) == 0
    for _ in range(moves):          # move the point along random lines: our tangent points differ from anyone else's
        while True:
            R = [Fr(rng.randint(-50, 50)) for _ in range(3)]
            if qe(A, R) == 0: continue
            t = -2 * sum(A[i][j] * P[i] * R[j] for i in range(3) for j in range(3)) / qe(A, R)
            Q = [P[i] + t * R[i] for i in range(3)]
            if t != 0 and any(Q[i] * P[j] != Q[j] * P[i] for i in range(3) for j in range(3)): break
        P = Q
    den = 1
    for x in P: den = den * x.denominator // gcd(den, x.denominator)
    P = [int(x * den) for x in P]; g = 0
    for x in P: g = gcd(g, x)
    P = [x // g for x in P]
    assert qe(A, P) == 0
    return P

def tangent(A, P):
    return [sum(P[i] * A[i][j] for i in range(3)) for j in range(3)]

def padic_to_rat(t, p):
    """a rational with the same class in Q_p*/Q_p*^2 as the p-adic t (needs enough precision)"""
    if t == 0: raise ValueError('zero')
    k = int(pari.valuation(t, p))
    u = t / pari(p) ** k
    pr = int(pari.padicprec(u, p))
    if pr < (4 if p == 2 else 2): raise ValueError('precision')
    return Fr(p) ** k * int(pari.lift(u))

def hilb(a, b, p):
    """Hilbert symbol of two non-zero rationals at p (p = 0: real)"""
    a = Fr(a); b = Fr(b)
    A = a.numerator * a.denominator; B = b.numerator * b.denominator
    if p == 0: return -1 if (A < 0 and B < 0) else 1
    return int(pari.hilbert(A, B, p))

class Cov:
    def __init__(s, es, a, rng):
        s.es = es; e1, e2, e3 = es; a1, a2 = int(a[0]), int(a[1]); a3 = a1 * a2
        s.a = (a1, a2, a3)
        s.Cu = [[a2, 0, 0], [0, -a3, 0], [0, 0, -(e3 - e2)]]
        s.Cv = [[a1, 0, 0], [0, -a3, 0], [0, 0, -(e3 - e1)]]
        s.Cw = [[a1, 0, 0], [0, -a2, 0], [0, 0, -(e2 - e1)]]
        s.l = {}
        for nm, C in (('u', s.Cu), ('v', s.Cv), ('w', s.Cw)):
            P = conic_point(C, rng); s.l[nm] = tangent(C, P)
        s.pts = {}
    def ok(s, x, p):
        for e, ai in zip(s.es, s.a):
            q = (x - e) / ai
            if q == 0: return False
            if p == 0:
                if q < 0: return False
            elif not is_sq_local(q, p): return False
        return True
    def local_points(s, p, rng, want):
        out = []
        H = max(abs(e) for e in s.es)
        for _ in range(400000):
            if len(out) >= want: break
            if p == 0:
                x = Fr(rng.randint(-10 * H, 10 * H), rng.randint(1, 97))
            else:
                mode = rng.randint(0, 2)
                if mode == 0:
                    e = rng.choice(s.es); x = e + Fr(rng.randint(1, p ** 6) * rng.choice((1, -1))) * Fr(p) ** rng.randint(-8, 16)
                elif mode == 1:
                    x = Fr(rng.randint(-p ** 10, p ** 10), p ** (2 * rng.randint(0, 10)))
                else:
                    x = Fr(rng.randint(-10 * H, 10 * H))
            if s.ok(x, p) and x not in out: out.append(x)
        assert len(out) == want, ('no local points', p)
        return out
    def lvals_padic(s, x, p, sg):
        a1, a2, a3 = s.a; e1, e2, e3 = s.es
        def sq(q):
            q = Fr(q)
            return pari.sqrt(pari('(%d)/(%d) + O(%d^%d)' % (q.numerator, q.denominator, p, PREC)))
        U = sg[0] * sq((x - e1) / a1); V = sg[1] * sq((x - e2) / a2); W = sg[2] * sq((x - e3) / a3)
        lu = sum(c * z for c, z in zip(s.l['u'], (V, W, 1)))
        lv = sum(c * z for c, z in zip(s.l['v'], (U, W, 1)))
        lw = sum(c * z for c, z in zip(s.l['w'], (U, V, 1)))
        return [padic_to_rat(t, p) for t in (lu, lv, lw)]
    def lvals_real(s, x, sg):
        a1, a2, a3 = s.a; e1, e2, e3 = s.es
        def sq(q): return mpmath.sqrt(mpmath.mpf(q.numerator) / q.denominator)
        U = sg[0] * sq((x - e1) / a1); V = sg[1] * sq((x - e2) / a2); W = sg[2] * sq((x - e3) / a3)
        res = []
        for nm, zz in (('u', (V, W, 1)), ('v', (U, W, 1)), ('w', (U, V, 1))):
            c = s.l[nm]; t = sum(mpmath.mpf(ci) * zi for ci, zi in zip(c, zz))
            scale = sum(abs(mpmath.mpf(ci) * zi) for ci, zi in zip(c, zz))
            if abs(t) < scale * mpmath.mpf(10) ** -80: raise ValueError('near zero')
            res.append(Fr(1 if t > 0 else -1))
        return res
    def local_value(s, b, p, rng, want=4):
        b1, b2 = b; vals = set()
        signs = [(i, j, k) for i in (1, -1) for j in (1, -1) for k in (1, -1)]
        for x in s.local_points(p, rng, want):
            for sg in signs:
                try:
                    lu, lv, lw = s.lvals_real(x, sg) if p == 0 else s.lvals_padic(x, p, sg)
                except ValueError: continue
                vals.add(hilb(b2, lv / lw, p) * hilb(b1, lu / lw, p))
        assert len(vals) == 1, ('not constant', p, vals)
        return vals.pop()

def places_of(E, extra=()):
    return [0] + sorted(set(E.primes) | set(extra))

def pairing_matrix(es, sel, seed=1, extra=(3, 5, 7, 11, 13), want=4):
    rng = random.Random(seed)
    E = Curve(*es)
    covs = [Cov(es, a, rng) for a in sel]
    pl = places_of(E, extra)
    M = []
    for i, D in enumerate(covs):
        row = []
        for b in sel:
            tot = 1; per = {}
            for p in pl:
                v = D.local_value(b, p, rng, want)
                per[p] = v; tot *= v
                if p not in E.primes and p != 0: assert v == 1, ('a good prime contributes', p)
            row.append(tot)
        M.append(row)
    return M

def f2rank(M):
    rows = [sum((0 if v == 1 else 1) << j for j, v in enumerate(r)) for r in M]
    B = []
    for v in rows:
        for b in B: v = min(v, v ^ b)
        if v: B.append(v)
    return len(B)

def sqf(q):
    q = Fr(q); n = q.numerator * q.denominator
    s = -1 if n < 0 else 1; n = abs(n)
    from sympy import factorint
    r = 1
    for p, k in factorint(n).items():
        if k % 2: r *= p
    return s * r

def kernel_classes(E, M, sel):
    """the kernel of the pairing matrix, as bit vectors of classes in Q(S,2)^2"""
    n = len(sel)
    rows = [sum((0 if v == 1 else 1) << j for j, v in enumerate(r)) for r in M]
    ker = [c for c in range(1, 1 << n) if all(bin(r & c).count('1') % 2 == 0 for r in rows)]
    selv = [E.classvec(a) for a in sel]
    out = []
    for _, c in span(ker):
        v = 0
        for i in range(n):
            if c >> i & 1: v ^= selv[i]
        out.append(v)
    return out

def analyse(es, seed=1, xs=None, label='', show_matrix=True):
    """dim Sel^2, the pairing matrix and its rank, the bound; with xs (x-coordinates of rational points of E, which the
    caller has checked), also: the images of T1, T2 and the points lie in Sel^2 and in the kernel of the pairing,
    and whether they span the whole kernel (then the bound is attained by the known points)."""
    E = Curve(*es)
    d = E.selmer(seed=seed + 100)
    sel = [(sqf(a), sqf(b)) for a, b in E.sel]
    M = pairing_matrix(es, sel, seed=seed)
    for i in range(len(sel)):
        assert M[i][i] == 1, 'not alternating'
        for j in range(len(sel)): assert M[i][j] == M[j][i], 'not symmetric'
    rho = f2rank(M)
    msg = '%s dim Sel^2 = %d, pairing rank %d, rank <= %d' % (label, d, rho, d - 2 - rho)
    res = dict(dim=d, rho=rho, bound=d - 2 - rho)
    if xs is not None:
        imgs = [E.kummer(Fr(es[0])), E.kummer(Fr(es[1]))] + [E.kummer(Fr(x)) for x in xs]
        for im in imgs: assert E.sel_contains(im), 'a point image outside Sel^2: the descent is wrong'
        iv = [E.classvec(im) for im in imgs]
        kv = kernel_classes(E, M, sel)
        si, sk, both = len(span(iv)), len(span(kv)), len(span(iv + kv))
        assert both == sk, 'the image of a rational point is not in the kernel of the pairing'
        res.update(img=si, ker=sk)
        msg += '; images of T1, T2 and %d point(s) span %d = %s the kernel (dim %d)' % (len(xs), si, 'all of' if si == sk else 'part of', sk)
        if si == sk: msg += ', so rank = %d exactly' % (si - 2)
    print(msg, flush=True)
    if show_matrix:
        for r in M: print('   ', ' '.join('%+d' % v for v in r))
    return res

if __name__ == '__main__':
    print('# ct_mine.py self-test (review z7; written by an AI model): Sel^2 and Cassels-Tate pairing on y^2 = x^3 - n^2 x')
    print('# reference: PARI 2.17.4 ellrank [r,R,s] recorded at review time (not computed here); s must equal the pairing rank, R the bound')
    REF = {17: (0, 0, 2), 73: (0, 0, 2), 205: (1, 1, 2), 34: (2, 2, 0), 41: (2, 2, 0), 113: (0, 2, 0), 5: (1, 1, 0)}
    for n, (r, R, s) in REF.items():
        res = analyse((-n, 0, n), label='n = %d (PARI [%d, %d, %d]):' % (n, r, R, s))
        assert res['rho'] == s and res['bound'] == R, n
    print('# all agree')
