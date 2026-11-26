#!/usr/bin/env python3
"""Independent re-check of the local terms of rank_fisher_theta4.py (E+ of (72,65), Xi = 137), without PARI's p-adic
arithmetic: every local point is re-validated by a multivariate Hensel lemma in pure integer arithmetic, and the class
of g0^0 = R3^0 / (W N1^0) is recomputed from integers.

    python3 rank_fisher_theta4_recheck.py [--write]     (writes rank_fisher_theta4_recheck.txt)

For each place p of S and each sampled point, the PARI p-adic point is only a STARTING VALUE: its coordinates are
lifted to integers x (scaled to be primitive), one coordinate x_j that is a p-adic unit is fixed, and for the three
equations F = (QA, QB, Xi W^2 - f0) and three of the remaining coordinates the Jacobian determinant J is computed.
Hensel's lemma (Newton's method in several variables): if v(F_i(x)) >= 2 v(J) + 1 for every i, there is a point
x* of C4(Z_p) with x* = x modulo p^(k - v(J)), k = min v(F_i(x)).  Then for an integral form G, G(x*) = G(x) modulo
p^(k - v(J)), so the class of G(x*) in Q_p^*/Q_p^*2 is that of G(x) when v(G(x)) + 3 (p = 2) or + 1 (odd p) is at
most k - v(J): asserted for G = R3^0 and G = W N1^0.  The term (g0(x*), 2)_p and (g0(x*), y)_p for the control
classes are recomputed and compared with rank_fisher_theta4.py's.  Written by an AI model (Claude); not read by a
person.  Experiment, not counted.
"""
import os, sys, itertools
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import sympy as sp
import rank_fisher_theta4 as T
import rank_second as R2
from rank_fisher_theta4 import pari, c, d, e, g, W

WRITE = '--write' in sys.argv
OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(s)

def vp(n, p):
    n = int(n)
    if n == 0: return 10 ** 9
    k = 0
    while n % p == 0: n //= p; k += 1
    return k

def ev(P_, x):
    tot = 0
    for mo, co in zip(P_.monoms(), P_.coeffs()):
        term = int(co)
        for xi, k in zip(x, mo):
            if k: term *= xi ** k
        tot += term
    return tot

def int_cls(n, p):
    """class of a non-zero integer: (v mod 2, unit class) as in rank_fisher_theta4.padic_class."""
    v = vp(n, p); u = n // p ** v
    if p == 2: return (v % 2, u % 8)
    return (v % 2, 0 if pow(u % p, (p - 1) // 2, p) == 1 else 1)

def main():
    say('rank_fisher_theta4_recheck.py: the local terms of rank_fisher_theta4.py re-checked by Hensel\'s lemma in '
        'integer arithmetic (E+ of (72,65), Xi = 137)')
    a, b, Gpt = R2.sigma_curve(72, 65, 'E+')
    ys = [2, T.sqfree(b), T.sqfree(Gpt[0])]
    cv, R = T.level3_data()
    C = cv['sides']["E'"]['level3_curve']
    D = T.setup(R)
    M = T.Model(D['QA'], D['QB'])
    Xi = 137
    G = T.build_g0(D, Xi); Nz = T.normalise(G)
    V5 = (c, d, e, g, W)
    Fs = [sp.Poly(D['QA'], *V5), sp.Poly(D['QB'], *V5), sp.Poly(sp.expand(Xi * W ** 2 - D['f0'].as_expr()), *V5)]
    J = [[sp.Poly(sp.diff(F.as_expr(), v_), *V5) for v_ in V5] for F in Fs]
    R30 = sp.Poly(Nz['R30'].as_expr(), *V5); WN1 = sp.Poly(sp.expand(W * Nz['N10'].as_expr()), *V5)
    ps, bad = T.degenerate_primes(D, Nz, say=lambda s: None)
    small = [p for p in range(2, 47) if pari.isprime(p)]
    S = sorted(set(small + T.B0 + [137] + bad))
    prod = {y: 1 for y in ys}
    total_pts = 0
    for p in S:
        pts = []
        for sd in (1, 2): pts += T.c4_points(M, D, C, Xi, p, 40, sd)
        terms = {y: set() for y in ys}
        n_ok = 0
        for P in pts:
            vmin = min(int(pari.valuation(z, p)) for z in P)
            N = min(int(pari.padicprec(z, p)) for z in P) - vmin
            x = [int(pari.lift(z / pari('%d' % p) ** vmin + pari('O(%d^%d)' % (p, N)))) for z in P]
            Fv = [ev(F, x) for F in Fs]
            k = min(vp(f_, p) for f_ in Fv)
            j0 = next(j for j in range(5) if x[j] % p)
            Jv = [[ev(J[i][j], x) for j in range(5)] for i in range(3)]
            best = None
            for cols in itertools.combinations([j for j in range(5) if j != j0], 3):
                det = sp.Matrix(3, 3, lambda i, jj: Jv[i][cols[jj]]).det()
                if det != 0 and (best is None or vp(det, p) < best): best = vp(det, p)
            if best is None or k < 2 * best + 1: continue          # not a certified approximation: skip it
            prec_ok = k - best
            num, den = ev(R30, x), ev(WN1, x)
            if num == 0 or den == 0: continue
            need = 3 if p == 2 else 1
            if vp(num, p) + need > prec_ok or vp(den, p) + need > prec_ok: continue
            cn, cd_ = int_cls(num, p), int_cls(den, p)
            # class of num/den = class of num * den (squares)
            if p == 2: cls = ((cn[0] + cd_[0]) % 2, (cn[1] * cd_[1]) % 8)
            else: cls = ((cn[0] + cd_[0]) % 2, (cn[1] + cd_[1]) % 2)
            for y in ys: terms[y].add(T.hilbert(cls, T.int_class(y, p), p))
            n_ok += 1
        assert n_ok > 0, 'no Hensel-certified point at p = %d' % p
        tp = {}
        for y in ys:
            assert len(terms[y]) == 1, 'the term at %d depends on the point' % p
            tp[y] = terms[y].pop(); prod[y] *= tp[y]
        total_pts += n_ok
        say('  p = %d: %d points certified by Hensel (of %d sampled); terms %s'
            % (p, n_ok, len(pts), ' '.join('%d:%+d' % (y, tp[y]) for y in ys)))
    say('  products over S (the real place contributes +1: all three classes are positive): %s'
        % ', '.join('%d: %+d' % (y, prod[y]) for y in ys))
    assert prod[2] == -1 and prod[ys[1]] == 1 and prod[ys[2]] == 1
    say('RE-CHECKED: Theta_4(y0, y1) product -1 (non-zero), the classes of E(Q) +1; %d Hensel-certified points' % total_pts)
    if WRITE: open(os.path.join(here, 'rank_fisher_theta4_recheck.txt'), 'w').write('\n'.join(OUT) + '\n')

if __name__ == '__main__':
    main()
