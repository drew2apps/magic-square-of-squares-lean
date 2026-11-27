#!/usr/bin/env python3
"""Independent re-check of the local terms of rank_fisher_theta4_gen.py (the general Theta_4 module), as
rank_fisher_theta4_recheck.py does for (72,65), with the Hilbert symbol from PARI's hilbert() on integers in place of
rank_fisher_theta4's own formula (a point raised by reviews/review_w4.md, NIT 3).

    python3 rank_fisher_theta4_gen_recheck.py [names] [--write]     (default 920,129 760,39 631,71@E'; writes
                                                                      rank_fisher_theta4_gen_recheck.txt)

For the first level-4 twist Xi of each curve and every place p of S (the same S as rank_fisher_theta4_gen.py), the
PARI p-adic points of C4 are only STARTING VALUES: the coordinates are lifted to integers x (made primitive), a unit
coordinate is fixed, and for F = (QA, QB, Xi W^2 - f0) and three other coordinates the Jacobian determinant J is
computed; if v(F_i(x)) >= 2 v(J) + 1 for all i, Hensel's lemma gives x* in C4(Z_p) with x* = x mod p^(k - v(J)),
k = min v(F_i(x)), and the classes of R3^0(x*) and W N1^0(x*) are those of R3^0(x), W N1^0(x) when their valuations
leave room (+3 at p = 2, +1 at odd p): asserted.  The term (g0^0(x*), y)_p = hilbert(R3^0(x) W N1^0(x), y, p)
(PARI, integers) must not depend on the point, and the product over S (the real place gives +1: every y used is
positive) must be rank_fisher_theta4_gen.py's value.  Written by an AI model (Claude); not read by a person.
Experiment, not counted.
"""
import os, sys, itertools
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import sympy as sp
import rank_fisher_theta4 as T
import rank_fisher_theta4_gen as TG
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

def det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))

def recheck(name, expected):
    TG.SIDE = None
    if '@' in name: name, TG.SIDE = name.split('@')
    P = TG.prepare(name, quiet=True)
    D, C, y, imgs, RD = P['D'], P['C'], P['y'], P['imgs'], P['RD']
    ys = [y] + [z for z in imgs if z != 1]
    assert all(z > 0 for z in ys), 'a negative class: the real term would need the sign of g0'
    nd, ps = TG.pencil_primes(D, known=sorted(set([2] + RD['bad'])))
    B0 = sorted(set([2] + RD['bad'] + ps))
    M = T.Model(D['QA'], D['QB'])
    places = sorted(set([2, 3, 5, 7, 11, 13] + B0))
    sols, K, Rs = TG.choose_Xi_lin(M, D, C, places, say=lambda s: None)
    Xi = sols[0]
    G = T.build_g0(D, Xi); Nz = T.normalise(G)
    V5 = (c, d, e, g, W)
    Fs = [sp.Poly(D['QA'], *V5), sp.Poly(D['QB'], *V5), sp.Poly(sp.expand(Xi * W ** 2 - D['f0'].as_expr()), *V5)]
    J = [[sp.Poly(sp.diff(F.as_expr(), v_), *V5) for v_ in V5] for F in Fs]
    R30 = sp.Poly(Nz['R30'].as_expr(), *V5); WN1 = sp.Poly(sp.expand(W * Nz['N10'].as_expr()), *V5)
    _, bad = T.degenerate_primes(D, Nz, say=lambda s: None)
    small = [p for p in range(2, 47) if pari.isprime(p)]
    extra = sorted(set(q_ for z in ys for q_ in T.primes_of(z)))
    S = sorted(set(small + B0 + T.primes_of(Xi) + bad + extra + TG.f0_degenerate_primes(D)))
    say('%s: Xi = %d, S = {R} u %s, classes %s' % (P['cv']['name'], Xi, S, ys))
    prod = {z: 1 for z in ys}
    total = 0
    for p in S:
        terms = {z: set() for z in ys}
        n_ok = 0; npts = 0
        for rnd, sds in enumerate(((1, 2, 7, 8), tuple(range(9, 21)))):
          if rnd and n_ok: break
          pts = []
          for sd in sds: pts += T.c4_points(M, D, C, Xi, p, 40 if not rnd else 150, sd)
          npts += len(pts)
          for Pt in pts:
              vmin = min(int(pari.valuation(z_, p)) for z_ in Pt)
              N = min(int(pari.padicprec(z_, p)) for z_ in Pt) - vmin
              x = [int(pari.lift(z_ / pari('%d' % p) ** vmin + pari('O(%d^%d)' % (p, N)))) for z_ in Pt]
              Fv = [ev(F, x) for F in Fs]
              k = min(vp(f_, p) for f_ in Fv)
              j0 = next(j for j in range(5) if x[j] % p)
              Jv = [[ev(J[i][j], x) for j in range(5)] for i in range(3)]
              best = None
              for cols in itertools.combinations([j for j in range(5) if j != j0], 3):
                  dt = det3([[Jv[i][cols[jj]] for jj in range(3)] for i in range(3)])
                  if dt != 0 and (best is None or vp(dt, p) < best): best = vp(dt, p)
              if best is None or k < 2 * best + 1: continue
              prec_ok = k - best
              num, den = ev(R30, x), ev(WN1, x)
              if num == 0 or den == 0: continue
              need = 3 if p == 2 else 1
              if vp(num, p) + need > prec_ok or vp(den, p) + need > prec_ok: continue
              for z in ys: terms[z].add(int(pari.hilbert(num * den, z, p)))
              n_ok += 1
        assert n_ok > 0, 'no Hensel-certified point at p = %d' % p
        tp = {}
        for z in ys:
            assert len(terms[z]) == 1, 'the term at %d depends on the point' % p
            tp[z] = terms[z].pop(); prod[z] *= tp[z]
        total += n_ok
        say('  p = %d: %d points certified by Hensel (of %d sampled); terms %s'
            % (p, n_ok, npts, ' '.join('%d:%+d' % (z, tp[z]) for z in ys)))
    say('  products over S (real place +1): %s' % ', '.join('%d: %+d' % (z, prod[z]) for z in ys))
    assert prod[y] == expected and all(prod[z] == 1 for z in ys[1:]), 'THE RE-CHECK DISAGREES'
    say('  RE-CHECKED: Theta_4 product %+d (as rank_fisher_theta4_gen.py), the classes of E(Q) +1; %d Hensel-certified '
        'points; symbols by PARI hilbert()' % (prod[y], total))
    return prod[y]

def main(names):
    say('rank_fisher_theta4_gen_recheck.py: the local terms of rank_fisher_theta4_gen.py re-checked by Hensel\'s lemma '
        'in integer arithmetic, symbols by PARI hilbert()')
    for nm in names:
        recheck(nm, -1)
    if WRITE: open(os.path.join(here, 'rank_fisher_theta4_gen_recheck.txt'), 'w').write('\n'.join(OUT) + '\n')

if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    main(args if args else ['920,129', '760,39', "631,71@E'"])
