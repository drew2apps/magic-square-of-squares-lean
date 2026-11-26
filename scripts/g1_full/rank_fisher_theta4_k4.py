#!/usr/bin/env python3
"""Control (K4) for rank_fisher_theta4.py: Theta_4(y0, y1) on E+ of (72,65) recomputed from DIFFERENT level-3 lifts
of y0 = 2: the same level-2 curve D_lam (lam = -1) with the other twists xi' of its level-3 Selmer coset
7 <-1, 2> = {7, -7, 14, -14} (rank_fisher.py).  Each C3_xi' is a new quadric intersection C3' (new conics Gamma_1,
Gamma_3, new points on them over Q(sqrt -2), Q(sqrt 2), new l_1, l_3, f, f_dag, q, c), so the whole Sect. 7
construction is redone from scratch; Fisher's Thm 2.1 says Theta_4(y0, y1) does not depend on the lift, so every
C3_xi' that lifts to level 4 must give the same value as rank_fisher_theta4.py (-1 as a product of Hilbert symbols).

    python3 rank_fisher_theta4_k4.py [--write]     (writes rank_fisher_theta4_k4.txt, rank_fisher_theta4_k4_hints.json)

The conic points are found by the method of rank_fisher5_find.py (PARI reduction modulo squares, then the Legendre
descent of rank_fisher5_legendre.py in Z[sqrt -2], Z[sqrt 2]): FINDERS, not certifying; rank_fisher5.analyse re-checks
each point exactly (on the cone, tangency along the ruling, f rational, q and c by exact linear algebra), and
rank_fisher_theta4 does the rest exactly.  A stored hints file is used when present (re-checked the same way).
Written by an AI model (Claude); not read by a person; not reviewed.  Experiment, not counted.
"""
import os, sys, json, time, subprocess
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import sympy as sp
import rank_fisher5 as R5
import rank_fisher_theta4 as T
import rank_second as R2
from rank_fisher_theta4 import pari

WRITE = '--write' in sys.argv
HFILE = os.path.join(here, 'rank_fisher_theta4_k4_hints.json')
OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(s)

def find_hint(R, pair):
    """rank_fisher5_find.py's finder as a function: a point of Gamma_1 over L (pair 'U') or of Gamma_3 over L'
    (pair 'V'), as a hint for rank_fisher5.analyse."""
    t = sp.Symbol('t'); e, g = R5.e, R5.g
    dU = R['dU']; h0, Bq = R['h0'], R['Bq']
    if pair == 'V':
        red_, t, Mw_, h0, Bq = R['coneV']
        dU = R['dV']; e, g = R5.c, R5.d
    fa = pari.factor(abs(dU))
    D0 = -1 if dU < 0 else 1; m = 1
    for p_, ex in zip(fa[0], fa[1]):
        D0 *= int(p_) ** (int(ex) % 2); m *= int(p_) ** (int(ex) // 2)
    assert m * m * D0 == dU and D0 in (-2, 2)
    def topari(expr):
        P = sp.Poly(sp.expand(expr), t)
        a = P.coeff_monomial(1); b = P.coeff_monomial(t)
        return pari('Mod((%s) + (%s)*%d*y, y^2 - (%d))' % (a, b, m, D0))
    P = sp.Poly(Bq, e, g)
    b20, b11, b02 = [topari(P.coeff_monomial(mm)) for mm in (e ** 2, e * g, g ** 2)]
    h0p = topari(h0)
    A1 = -4 * b20 * h0p; B1 = b11 ** 2 - 4 * b20 * b02          # E^2 - A1 Z^2 = B1 g^2, E = 2 b20 e + b11 g
    bnf = pari.bnfinit(pari('y^2 - (%d)' % D0), 1)
    def sqred(z):
        z = z.lift(); s = pari('Mod(1, y^2 - (%d))' % D0)
        dd = pari.denominator(pari.content(z)); zz = pari('Mod(%s, y^2 - (%d))' % (z * dd * dd, D0)); s = s / dd
        fz = pari.idealfactor(bnf, zz.lift())
        for i in range(int(pari.matsize(fz)[0])):
            Pr, ex = fz[i, 0], int(fz[i, 1])
            if ex >= 2:
                ip = pari.bnfisprincipal(bnf, Pr)
                pi = pari.nfbasistoalg(bnf, ip[1])
                zz = zz / pi ** (2 * (ex // 2)); s = s * pi ** (ex // 2)
        return zz, s
    A2, sa = sqred(A1); B2, sb = sqred(B1)
    assert A2 == A1 / sa ** 2 and B2 == B1 / sb ** 2
    def ab(z):
        z = z.lift(); return int(pari.polcoef(z, 0, 'y')), int(pari.polcoef(z, 1, 'y'))
    out_ = subprocess.run([sys.executable, '-B', os.path.join(here, 'rank_fisher5_legendre.py')] +
                          [str(v) for v in ab(A2) + ab(B2)], capture_output=True, text=True, timeout=600,
                          env=dict(os.environ, LEG_D=str(-D0)))
    assert out_.returncode == 0, ('Legendre descent failed', out_.stderr[-300:])
    sol = json.loads(out_.stdout)
    el = lambda v: pari('Mod(%d + %d*y, y^2 - (%d))' % (v[0], v[1], D0))
    xx, yy, zz = el(sol['x']), el(sol['y']), el(sol['z'])
    assert xx ** 2 == A2 * yy ** 2 + B2 * zz ** 2
    E_ = xx; Z = yy / sa; gg = zz / sb
    assert E_ ** 2 - A1 * Z ** 2 == B1 * gg ** 2
    ee = (E_ - b11 * gg) / (2 * b20)
    assert h0p * Z ** 2 + b20 * ee ** 2 + b11 * ee * gg + b02 * gg ** 2 == 0
    def out(z):
        z = z.lift(); return [str(pari.polcoef(z, 0, 'y')), str(pari.polcoef(z, 1, 'y'))]
    if pair == 'U': return dict(D0=D0, m=m, Z=out(Z), e=out(ee), g=out(gg))
    return dict(D0=D0, m=m, Z=out(Z), c=out(ee), d=out(gg))

def main():
    say('rank_fisher_theta4_k4.py: control (K4), Theta_4(y0, y1) on E+ of (72,65) from other level-3 lifts of y0')
    a, b, Gpt = R2.sigma_curve(72, 65, 'E+')
    yb, yG = T.sqfree(b), T.sqfree(Gpt[0])
    ys = [2, yb, yG]
    cv, R0 = T.level3_data()
    name, side = cv['name'], "E'"
    C = cv['sides'][side]['level3_curve']
    saved = json.load(open(HFILE)) if os.path.exists(HFILE) else {}
    hints_out = {}
    values = {}
    for xi2 in (14, -7, -14):
        t0 = time.time()
        C2 = dict(C, xi=str(xi2))
        key = '%s|%s|xi=%d' % (name, side, xi2)
        R5.HINTS.clear()
        Rb = R5.analyse(name, side, C2, dd=int(cv['sides'][side]['dd']))
        hU = saved.get(key + '|U') or find_hint(Rb, 'U')
        hV = saved.get(key + '|V') or find_hint(Rb, 'V')
        hints_out[key + '|U'] = hU; hints_out[key + '|V'] = hV
        R5.HINTS['%s|%s' % (name, side)] = hU; R5.HINTS['%s|%s|V' % (name, side)] = hV
        R = R5.analyse(name, side, C2, dd=int(cv['sides'][side]['dd']))
        assert R['pushout'] is not None and R['pushout3'] is not None and R['comp'] is not None
        say('C3 with lam = %s, xi = %d: rank_fisher5 steps (1)-(8) re-checked with new conic points over Q(sqrt -2), '
            'Q(sqrt 2) (%d, %d digits)' % (C['lam'], xi2, max(len(x) for x in hU['Z']), max(len(x) for x in hV['Z'])))
        D = T.setup(R)
        nd, ps = T.check_smooth_primes(D)
        B0 = sorted(set(T.B0) | set(ps))
        S_ELS = sorted(set([2, 3, 5, 7, 11, 13] + B0))
        say('  the pencil discriminant (%d digits) has the primes %s: B0 here = %s' % (nd, ps, B0))
        M = T.Model(D['QA'], D['QB'])
        sols, K, Rs = T.choose_Xi(M, D, C2, say=lambda s: None, places=S_ELS)
        sols = sorted(sols, key=abs)
        say('  level-4 Selmer set over the places %s: %d twists %s' % (S_ELS, len(sols), sols[:8]))
        if not sols:
            say('  this level-3 lift does not lift to level 4 (Sect. 3: then <[C3], eta> != 0 for some eta in S\'_1); '
                'not usable for Theta_4')
            continue
        assert len(sols) == 8, 'the Selmer set is not a coset of S\'_1'
        res, S = T.theta4_for(M, D, C2, sols[0], ys, say=say, label='  ', B0=B0)
        assert res[yb] == 1 and res[yG] == 1, 'CONTROL (K2) FAILED'
        values[xi2] = res[2]
        say('  => from this lift: product of the local symbols %+d' % res[2])
    say('')
    say('(K4) values from the lifts xi\' = %s: %s; rank_fisher_theta4.py (xi = 7): -1' % (sorted(values), values))
    assert values and all(v == -1 for v in values.values()), 'CONTROL (K4) FAILED: the value depends on the lift'
    say('(K4) passed: every level-3 lift that lifts to level 4 gives Theta_4(y0, y1) != 0')
    if WRITE:
        open(os.path.join(here, 'rank_fisher_theta4_k4.txt'), 'w').write('\n'.join(OUT) + '\n')
        json.dump(hints_out, open(HFILE, 'w'), indent=1, sort_keys=True)

if __name__ == '__main__':
    main()
