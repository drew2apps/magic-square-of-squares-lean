#!/usr/bin/env python3
"""Step 9 of Fisher's fifth descent (preparation only): local points of the level-3 curve C3' (rank_fisher5.py) at
the finite places of S3, from rank_fisher's certified local points of C3, and the class of the pushout form
f0 (rank_fisher5_pushout.json) there.  These classes are the data f(P_v) of the level-4 Selmer set
S(C3', f) = { xi4 : xi4 f(P_v) in Im(delta_v) for all v } (Fisher Sect. 3); the local images are NOT computed here.

    python3 rank_fisher5_local.py [--write]     (72,65), E' side only

Exact p-adic arithmetic: integers modulo p^N with the precision asserted.  The real place is not treated.
Written by an AI model (Claude); not read by a person.  Experiment, not counted.
"""
import os, sys, json
from fractions import Fraction as Fr
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import sympy as sp
import rank_second as R2
import rank_indep as RI
import rank_fisher5 as F5

WRITE = '--write' in sys.argv
N = 240
OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(s)

def vp(x, p): return RI.vp(x, p) if x != 0 else 10 ** 9

def main(mutate=None):
    js = json.load(open(os.path.join(here, 'rank_fisher_curves.json')))
    cv = js[2]; side = "E'"
    C = cv['sides'][side]['level3_curve']
    R = F5.analyse(cv['name'], side, C, dd=int(cv['sides'][side]['dd']))
    push = json.load(open(os.path.join(here, 'rank_fisher5_pushout.json')))['%s|%s' % (cv['name'], side)]
    c, d, e, g = F5.c, F5.d, F5.e, F5.g
    f0 = sp.Poly(sp.sympify(push['f0']), c, d, e, g)
    QA = sp.Poly(R['QA'], c, d, e, g); QB = sp.Poly(R['QB'], c, d, e, g)
    assert all(sp.Rational(k).q == 1 for k in QA.coeffs() + QB.coeffs() + f0.coeffs())
    ev = lambda P_, pt: sum(int(k) * pt[0] ** m[0] * pt[1] ** m[1] * pt[2] ** m[2] * pt[3] ** m[3]
                            for m, k in zip(P_.monoms(), P_.coeffs()))
    lam = int(C['lam']); U = [int(t) for t in C['U']]; V = [int(t) for t in C['V']]
    pX = [int(t) for t in C['p_X']]; pY = [int(t) for t in C['p_Y']]
    pp, kX, L1 = F5.parametrise(lam, U, pX, [int(t) for t in C['l_X']], c, d)
    P1 = [sp.Poly(x_, c, d) for x_ in pp]
    say('rank_fisher5_local.py: local points of C3\' for %s, %s side, and the classes of f0 there' % (cv['name'], side))
    rows = []
    for v, pt in sorted(C['local_points'].items(), key=lambda kv: (kv[0] != 'inf', int(kv[0]) if kv[0] != 'inf' else 0)):
        if v == 'inf': continue
        p = int(v)
        r, s_, sx, sy = Fr(pt[0]), Fr(pt[1]), int(pt[2]), int(pt[3])
        D = r.denominator * s_.denominator // __import__('math').gcd(r.denominator, s_.denominator)
        ri, si = int(r * D), int(s_ * D)
        if mutate == 'point': ri += p ** 3
        X = sx * R2.padic_sqrt(lam * R2.qv(U, ri, si), p, N)
        Y = sy * R2.padic_sqrt(lam * R2.qv(V, ri, si), p, N)
        mod = p ** N
        assert (X * X - lam * R2.qv(U, ri, si)) % p ** (N - 2) == 0 and (Y * Y - lam * R2.qv(V, ri, si)) % p ** (N - 2) == 0
        # (c : d) and (e : g): projections from p_X and p_Y onto X = 0 and Y = 0
        cc = pX[2] * ri - X * pX[0]; dd_ = pX[2] * si - X * pX[1]
        ee = pY[2] * ri - Y * pY[0]; gg = pY[2] * si - Y * pY[1]
        Qs = QA if vp(ri, p) <= vp(si, p) else QB        # use the coordinate r or s of smaller valuation
        A = ev(Qs, [cc, dd_, 0, 0])                     # = p1(c, d) or p2(c, d) (QA = p1(c,d) - Xi q1(e,g), ...)
        B = -ev(Qs, [0, 0, ee, gg])                     # = Xi q1(e, g) or Xi q2(e, g)
        assert A % mod != 0 and B % mod != 0
        AB = A * B
        k = vp(AB % p ** (N - 4), p)
        assert k < N - 10, 'precision'
        u = (AB // p ** k) if k > 0 else AB
        ok = (k % 2 == 0) and (u % 8 == 1 if p == 2 else RI.legendre(u % p, p) == 1)
        assert ok, 'p1(c,d) / (Xi q1(e,g)) is not a square: the point does not lift to C3\''
        al = R2.padic_sqrt(AB % p ** (N + k), p, N) if True else None
        Pt = [B * cc, B * dd_, al * ee, al * gg] if mutate != 'swap' else [B * cc, B * dd_, al * gg, al * ee]
        # the point is on C3' modulo p^M, M = N minus the valuations used
        res = [ev(QA, Pt), ev(QB, Pt)]
        vmin = min(vp(x_, p) for x_ in Pt)
        Ms = [min(N - 4, vp(x_ % p ** (N - 4), p)) for x_ in res]
        assert all(M_ >= 2 * vmin + 12 for M_ in Ms), ('not on C3\' to the required precision', v, Ms, vmin)
        cls = []
        for sg in (1, -1):
            P2 = [Pt[0], Pt[1], sg * Pt[2], sg * Pt[3]]
            y = ev(f0, P2)
            vy = vp(y % p ** (N - 4), p)
            assert vy <= N - 12, 'precision of f0(P)'
            uy = (y // p ** vy) % (8 if p == 2 else p)
            # normalise by the scale p^(2 vmin) (P is a representative with min valuation vmin)
            vy2 = vy - 2 * vmin
            sq = (uy % 8 == 1) if p == 2 else (RI.legendre(uy, p) == 1)
            cls.append((vy2 % 2, (uy % 8) if p == 2 else ('QR' if sq else 'NQR')))
        rows.append((p, cls))
        say('  p = %d: a point of C3\' modulo p^%d (QA, QB checked; coordinates of minimal valuation %d); '
            'f0 at it and at (c:d:-e:-g): (v mod 2, unit class) = %s, %s'
            % (p, min(Ms), vmin, cls[0], cls[1]))
    say('(the classes give one sample f(P_v) per place; the level-4 Selmer set also needs Im(delta_v): not computed)')
    return rows

if __name__ == '__main__':
    if '--negative' in sys.argv:
        res = ['rank_fisher5_local.py --negative']
        for m_, what in (('point', 'a moved local point (r + p^3: not on D_lam, no p-adic square root)'),
                         ('swap', 'e and g swapped in the lifted point (not on C3\')')):
            try:
                main(mutate=m_); res.append('NOT REFUSED: ' + what)
            except AssertionError as ex:
                res.append('refused: %s %s' % (what, str(ex)[:60]))
        res.append('%d of 2 refused' % sum(r_.startswith('refused') for r_ in res))
        print('\n'.join(res))
        assert res[-1].startswith('2 of 2')
        if WRITE: open(os.path.join(here, 'rank_fisher5_local_negative.txt'), 'w').write('\n'.join(res) + '\n')
    else:
        main()
        if WRITE: open(os.path.join(here, 'rank_fisher5_local.txt'), 'w').write('\n'.join(OUT) + '\n')
