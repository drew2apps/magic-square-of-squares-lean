#!/usr/bin/env python3
"""External controls for rank_fisher_theta4.py: the same Theta_4 code on curves whose answer is known beforehand.

    python3 rank_fisher_theta4_control.py [--write]     (rank_fisher_theta4_control.txt, ..._control_hints.json)

THE CONTROL CURVES (found by a search over y^2 = x(x^2 + a x + b) with a^2 - 4b = -(power of 2)^2, PARI ellrank
upper bound - lower bound >= 2, L(E, 1) != 0):
    A: y^2 = x(x^2 + 128 x + 4097)    (conductor 262208)     B: y^2 = x(x^2 + 74 x + 17753)    (conductor 284048)
    C: y^2 = x(x^2 + 1310 x + 429041) (conductor 3432328)    D: y^2 = x(x^2 + 1706 x + 743993) (conductor 11903888)
For both: L(E, 1) != 0 (6.28..., 9.71...; PARI lfun), so rank E(Q) = 0 (Kolyvagin); E(Q)_tors = E'(Q)_tors = Z/2;
the first descent gives S_1 = {+-1, +-2} (E' side, known <-1>, y0 = 2) and S'_1 of dimension 2 (E side, known <b>,
y1 = 17 resp. 41); Theta_2(y0, y1) = +1 (rank_second's pairing, recomputed here), so S_2 = S_1, S'_2 = S'_1, and as
in rank_fisher.py Theta_3 = 0 (alternating, kernel of codimension <= 1): S_4 = S_1, S'_4 = S'_1, the 4-descent
bound is 2 + 2 - 2 = 2 > 0 = rank.
THE KNOWN ANSWER (conditional on the Birch--Swinnerton-Dyer formula for the order of Sha, the only input not proved):
y0 is in S_5 (the image of Sel^(4 phi)(E) in S_1) iff iota(y0) lies in 4 Sha(E) (a in Sha(E) with 4a = iota(y0) lies
in Sha(E)[4 phi] and so comes from Sel^(4 phi)(E); its image is y0 up to a class of E'(Q)).  |Sha(E)[phi]| =
|S_1 / known| = 2 and |Sha(E')[phihat]| = 2 give |Sha(E)[2]| <= 4.  PARI's BSD quotient L(E,1) / ellbsd(E) is
  A, B: 16.000...  so Sha(E)[2^oo] = (Z/4)^2, 4 Sha(E)[2^oo] = 0, iota(y0) != 0 (y0 is not a class of E'(Q)):
        y0 not in S_5, and as S'_4 = <known, y1>, Theta_4(y0, y1) != 0: the code must return -1;
  C, D: 64.000...  so Sha(E)[2^oo] = (Z/8)^2, Sha(E)[2] = 4 Sha(E)[2^oo] contains iota(y0): y0 in S_5, the left
        kernel of Theta_4, so Theta_4(y0, y1) = 0: the code must return +1.
For all four the class b of T (the image of E(Q)) must give +1.
The level-2 and level-3 curves are built with rank_indep / rank_second / rank_fisher (their certified routines), the
cone-conic points with the finder of rank_fisher_theta4_k4.py (Legendre descent in Z[sqrt -2], Z[sqrt 2]), every
step re-checked exactly as for (72,65).
Written by an AI model (Claude); not read by a person; not reviewed.  Experiment, not counted.
"""
import os, sys, json, time
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import rank_indep as RI
import rank_second as R2
import rank_fisher as RF
import rank_fisher5 as R5
import rank_fisher_theta4 as T
import rank_fisher_theta4_k4 as K4
from rank_fisher_theta4 import pari

WRITE = '--write' in sys.argv
HFILE = os.path.join(here, 'rank_fisher_theta4_control_hints.json')
OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(s)

CURVES = [('control A', 128, 4097, 2, 17, 16), ('control B', 74, 17753, 2, 41, 16),
          ('control C', 1310, 429041, 2, 409, 64), ('control D', 1706, 743993, 2, 449, 64)]

def control(name, a, b, y0, y1, sha_expected, saved, hints_out):
    t0 = time.time()
    ap, bp = -2 * a, a * a - 4 * b
    E = pari.ellinit([0, a, 0, b, 0]); Ep = pari.ellinit([0, ap, 0, bp, 0])
    L1 = pari.lfun(E, 1); sha = L1 / pari.ellbsd(E)
    say('%s: E: y^2 = x(x^2 + %d x + %d), conductor %s; L(E,1) = %s (non-zero: rank 0); BSD quotient for Sha(E) '
        '= %s; torsion E %s, E\' %s' % (name, a, b, pari.ellglobalred(E)[0], pari.round(L1 * 10 ** 6) / 10 ** 6,
                                      pari.round(sha * 10 ** 6) / 10 ** 6, pari.elltors(E)[1], pari.elltors(Ep)[1]))
    assert abs(float(L1)) > 1 and abs(float(sha) - sha_expected) < 1e-6
    R = RI.first_descent(a, b)
    assert sorted(R['Sp']) == [-2, -1, 1, 2] and len(R['S']) == 4
    for dd, AA, BB, loc in ((y1, a, b, R['locS']), (y0, ap, bp, R['locSp'])):
        assert all(R2.certified(dd, AA, BB, p, loc) for p in R['bad']) and RI.real_ok(dd, AA, BB)
    yb = RI.sqfree_part(b)
    say('  first descent: S_1 = %s (E\' side; known <-1>, y0 = %d), S\'_1 = %s (E side; known <%d>, y1 = %d); '
        'y0, y1 certified (rank_indep local images)' % (sorted(R['Sp']), y0, sorted(R['S']), yb, y1))
    # level 2 (E' side) and Theta_2
    P0 = RI.find_conic_point(y0, ap, bp // y0)
    surv, info = RI.second_descent(y0, ap, bp, P0, stop_at_first=False)
    R2.check_conic_forms(y0, ap, bp // y0, info['U'], info['V'], P0[2])
    tried = []
    chosen = None
    for lam in sorted(surv, key=abs):
        try:
            D = R2.DLam(lam, info['U'], info['V'], bad=info['primes'])
            th2, _ = D.pairing(y1)
        except AssertionError:
            continue
        assert th2 == 1, 'Theta_2(y0, y1) != 0: not a Theta_4 control'
        res = RF.selmer_set_level3(D, R, "E'")
        if res['nonempty'] is not True: continue
        RF.certify_C3(res)
        # the level-3 lifts over this D_lam: xi' in xi <-1, 2> (rank_fisher: a coset of S_1); each C3_xi' is tried
        for z in (1, -1, 2, -2):
            xi2 = res['xi'] * z
            C3 = dict(lam=str(D.lam), U=[str(x) for x in D.U], V=[str(x) for x in D.V], l_X=[str(x) for x in D.X['l']],
                      l_Y=[str(x) for x in D.Y['l']], p_X=[str(x) for x in D.X['P']], p_Y=[str(x) for x in D.Y['P']],
                      xi=str(xi2), local_points={str(v): [str(x) for x in P] for v, P in res['cert'].items()})
            key = '%s|E\'|lam=%d|xi=%d' % (name, D.lam, xi2)
            R5.HINTS.clear()
            try:
                Rb = R5.analyse(name, "E'", C3, dd=y0)
                hU = saved.get(key + '|U') or K4.find_hint(Rb, 'U')
                hV = saved.get(key + '|V') or K4.find_hint(Rb, 'V')
            except AssertionError as ex:
                msg = ex.args[0] if ex.args else type(ex).__name__
                msg = msg[0] if isinstance(msg, tuple) else msg  # PARI's error context (its padding varies) is dropped
                tried.append((D.lam, xi2, 'no conic point: %s' % str(msg)[:40])); continue
            R5.HINTS['%s|%s' % (name, "E'")] = hU; R5.HINTS['%s|%s|V' % (name, "E'")] = hV
            R3 = R5.analyse(name, "E'", C3, dd=y0)
            assert R3['pushout'] is not None and R3['pushout3'] is not None and R3['comp'] is not None
            Dd = T.setup(R3)
            nd, ps = T.check_smooth_primes(Dd, known=R['bad'])
            B0 = sorted(set([2] + R['bad'] + ps))
            S_ELS = sorted(set([2, 3, 5, 7, 11, 13] + B0))
            M = T.Model(Dd['QA'], Dd['QB'])
            sols, K, Rs = T.choose_Xi(M, Dd, C3, say=lambda s: None, places=S_ELS)
            tried.append((D.lam, xi2, len(sols)))
            if sols:
                hints_out[key + '|U'] = hU; hints_out[key + '|V'] = hV
                chosen = (D, xi2, C3, R3, Dd, nd, ps, B0, M, sorted(sols, key=abs)); break
        if chosen: break
    say('  level-3 lifts tried (lam, xi, number of level-4 twists found): %s' % tried)
    assert chosen, 'no level-3 lift with a level-4 lift found'
    D, xi2, C3, R3, Dd, nd, ps, B0, M, sols = chosen
    say('  level 2: D_lam with lam = %d checked, Theta_2(y0, y1) = +1; level 3: C3 with xi = %d (the first lift with level-4 '
        'twists found by the real sampler; a count of 0 for another lift is a miss of the sampler, not an empty '
        'level-4 Selmer set)' % (D.lam, xi2))
    say('  rank_fisher5 steps (1)-(8) on C3\': cones over Q(sqrt %d), Q(sqrt %d); conic points found and checked'
        % (T.sqfree(R3['dU']), T.sqfree(R3['dV'])))
    say('  pencil discriminant (%d digits): primes %s; B0 = %s' % (nd, ps, B0))
    say('  level-4 Selmer set: %d twists %s (a coset of S\'_1, dimension 2)' % (len(sols), sols))
    assert len(sols) == 4
    out, S = T.theta4_for(M, Dd, C3, sols[0], [y1, yb], say=say, label='  ', B0=B0)
    exp = -1 if sha_expected == 16 else 1
    say('  => %s: Theta_4(y0, y1) product %+d (expected %+d), class b of T: %+d (expected +1)'
        % (name, out[y1], exp, out[yb]))
    return out[y1], out[yb], exp

def main():
    say('rank_fisher_theta4_control.py: rank_fisher_theta4\'s Theta_4 code on two curves with a known answer')
    saved = json.load(open(HFILE)) if os.path.exists(HFILE) else {}
    hints_out = {}
    results = {}
    for (name, a, b, y0, y1, sh) in CURVES:
        results[name] = control(name, a, b, y0, y1, sh, saved, hints_out)
    say('')
    for name, (v1, vb, ex) in results.items():
        say('%s: Theta_4(y0, y1) %s (predicted %s), image of E(Q) %s'
            % (name, 'non-zero' if v1 == -1 else 'zero', 'non-zero' if ex == -1 else 'zero',
               '0 (as required)' if vb == 1 else 'NON-ZERO (WRONG)'))
    assert all(v1 == ex and vb == 1 for (v1, vb, ex) in results.values()), 'A CONTROL FAILED'
    say('controls passed')
    if WRITE:
        open(os.path.join(here, 'rank_fisher_theta4_control.txt'), 'w').write('\n'.join(OUT) + '\n')
        json.dump(hints_out, open(HFILE, 'w'), indent=1, sort_keys=True)

if __name__ == '__main__':
    main()
