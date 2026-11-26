#!/usr/bin/env python3
"""Control (K6) for rank_fisher_theta4.py: Theta_4(y0, y1) on E+ of (72,65) from level-3 lifts over OTHER level-2
curve D_lam (lam = -2 instead of rank_fisher_theta4's lam = -1; lam = -79553041 was also tried: its level-4 twist
needs a point of C4 over Q_97 in a class of odd valuation, which the samplers did not find, so it is not used).  The level-2 curve, the level-3
twist xi (rank_fisher.selmer_set_level3, certified by rank_fisher.certify_C3), the quadric intersection C3', the conic
points (finder of rank_fisher_theta4_k4.py, re-checked by rank_fisher5.analyse), the level-4 Selmer set, the pushout g0
and the local sum are all recomputed.  Fisher's Thm 2.1: the value cannot depend on the lift; it must be -1.

    python3 rank_fisher_theta4_k6.py [--write]     (rank_fisher_theta4_k6.txt, rank_fisher_theta4_k6_hints.json)

Written by an AI model (Claude); not read by a person; not reviewed.  Experiment, not counted.
"""
import os, sys, json
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import rank_indep as RI
import rank_second as R2
import rank_fisher as RF
import rank_fisher5 as R5
import rank_fisher_theta4 as T
import rank_fisher_theta4_k4 as K4

WRITE = '--write' in sys.argv
HFILE = os.path.join(here, 'rank_fisher_theta4_k6_hints.json')
OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(s)

def main():
    say('rank_fisher_theta4_k6.py: control (K6), Theta_4(y0, y1) on E+ of (72,65) over other level-2 lifts D_lam of y0')
    a, b, Gpt = R2.sigma_curve(72, 65, 'E+')
    ap, bp = -2 * a, a * a - 4 * b
    y0, y1 = 2, 2
    yb, yG = T.sqfree(b), T.sqfree(Gpt[0])
    ys = [y1, yb, yG]
    R = RI.first_descent(a, b)
    P0 = RI.find_conic_point(y0, ap, bp // y0)
    surv, info = RI.second_descent(y0, ap, bp, P0, stop_at_first=False)
    R2.check_conic_forms(y0, ap, bp // y0, info['U'], info['V'], P0[2])
    saved = json.load(open(HFILE)) if os.path.exists(HFILE) else {}
    hints_out, values = {}, {}
    name, side = '72,65 E+', "E'"
    for lam in (-2,):
        assert lam in surv
        D = R2.DLam(lam, info['U'], info['V'], bad=info['primes'])
        th2, _ = D.pairing(y1)
        assert th2 == 1
        res = RF.selmer_set_level3(D, R, side)
        assert res['nonempty'] is True
        RF.certify_C3(res)
        C3 = dict(lam=str(D.lam), U=[str(x) for x in D.U], V=[str(x) for x in D.V], l_X=[str(x) for x in D.X['l']],
                  l_Y=[str(x) for x in D.Y['l']], p_X=[str(x) for x in D.X['P']], p_Y=[str(x) for x in D.Y['P']],
                  xi=str(res['xi']), local_points={str(v): [str(x) for x in P] for v, P in res['cert'].items()})
        key = '%s|%s|lam=%d|xi=%d' % (name, side, lam, res['xi'])
        R5.HINTS.clear()
        Rb = R5.analyse(name, side, C3, dd=y0)
        hU = saved.get(key + '|U') or K4.find_hint(Rb, 'U')
        hV = saved.get(key + '|V') or K4.find_hint(Rb, 'V')
        hints_out[key + '|U'] = hU; hints_out[key + '|V'] = hV
        R5.HINTS['%s|%s' % (name, side)] = hU; R5.HINTS['%s|%s|V' % (name, side)] = hV
        R3 = R5.analyse(name, side, C3, dd=y0)
        assert R3['pushout'] is not None and R3['pushout3'] is not None and R3['comp'] is not None
        say('D_lam with lam = %d (Theta_2 = +1), level-3 twist xi = %d (certified); C3\' re-built, conic points over '
            'Q(sqrt %d), Q(sqrt %d) found and checked' % (lam, res['xi'], T.sqfree(R3['dU']), T.sqfree(R3['dV'])))
        Dd = T.setup(R3)
        nd, ps = T.check_smooth_primes(Dd, known=sorted(set(T.B0) | set(R['bad'])))
        B0 = sorted(set(T.B0) | set(ps))
        S_ELS = sorted(set([2, 3, 5, 7, 11, 13] + B0))
        say('  pencil discriminant (%d digits): primes %s' % (nd, ps))
        M = T.Model(Dd['QA'], Dd['QB'])
        sols, K, Rs = T.choose_Xi(M, Dd, C3, say=lambda s: None, places=S_ELS)
        if not sols:
            # the sampled classes of f0 may miss part of a coset: close them with the certified local images
            # Im(delta_phihat, p) of rank_indep (bad primes of E) and the unit classes (other odd p); points are
            # then found and checked for the chosen Xi at every place (local_terms)
            ims = {p: (T.local_image_classes(R, p) if p in R['bad'] else [(0, 0), (0, 1)]) for p in S_ELS if p != 2 or p in R['bad']}
            sols, K, Rs = T.choose_Xi(M, Dd, C3, say=lambda s: None, places=S_ELS, images=ims)
            say('  (sampled classes closed up to cosets of the certified local images)')
        sols = sorted(sols, key=abs)
        say('  level-4 Selmer set over %s: %d twists %s' % (S_ELS, len(sols), sols[:8]))
        assert len(sols) == 8, 'the level-4 Selmer set is not a coset of S\'_1'
        out, S = T.theta4_for(M, Dd, C3, sols[0], ys, say=say, label='  ', B0=B0)
        assert out[yb] == 1 and out[yG] == 1, 'CONTROL (K2) FAILED'
        values[lam] = out[y1]
        say('  => over lam = %d: product of the local symbols %+d' % (lam, out[y1]))
    say('')
    say('(K6) values over lam = %s: %s; rank_fisher_theta4.py (lam = -1): -1' % (sorted(values), values))
    assert all(v == -1 for v in values.values()), 'CONTROL (K6) FAILED'
    say('(K6) passed')
    if WRITE:
        open(os.path.join(here, 'rank_fisher_theta4_k6.txt'), 'w').write('\n'.join(OUT) + '\n')
        json.dump(hints_out, open(HFILE, 'w'), indent=1, sort_keys=True)

if __name__ == '__main__':
    main()
