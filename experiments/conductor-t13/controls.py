#!/usr/bin/env python3
"""controls.py (conductor-t13): the whole pipeline of covers.py + search.py + certify.py on control curves over
K = Q(i) with all 2-torsion over K and two independent points known by other means, to check that the 2-covering
quartics, their minimisation and reduction, the sieve and the map back find points of infinite order and certify
rank >= 2.  EXPERIMENT; written by an AI model (build agent t13), not read by a person.

    python3 controls.py [--write]   -> controls.txt

Controls: (a) y^2 = x(x + a)(x + b) over Q, base changed to K, with rank E(Q) = rank E^-1(Q) = 1 (the controls of
kctp_points.py; rank E(K) = 2); (b) curves y^2 = x(x - a)(x - b), a, b in Z[i] small with j not in Q, where a naive
search over small x in Z[i] finds two points certified independent (kctp_points.certify), so rank >= 2.
For each control: Sel^2 (kctp_recheck.Selmer), the cosets modulo the torsion image, one quartic per coset, a search
with H = 2000, the hits mapped back and certified.  Pass: every coset with dim Sel^2 = rank + 2 (Sha[2] = 0) gives a
point of infinite order, and the points found certify rank >= 2.
"""
from sage.all__sagemath_schemes import *
import sys, os, json, random, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', 'scripts', 'second_descent'))
import kctp_recheck as R
import kctp_points as KP
import covers as CV
import search as SE
import certify as CE

K = R.K; I = R.I
WRITE = '--write' in sys.argv
OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(str(s))

def run_curve(es, rng, name, H=2000):
    sel, E, comp, tbasis, rt = CV.curve_data(es, rng, name)
    vecs = [comp[0], comp[1], [(x + y) % 2 for x, y in zip(comp[0], comp[1])]] if len(comp) >= 2 else [comp[0]]
    pts = []; percoset = []
    for v in vecs:
        dd1, dd2, dd3 = sel.pair(v)
        c = CV.build(es, dd1, dd2, dd3, rng, 'coset')
        g = [(int(QQ(K(x)[0])), int(QQ(K(x)[1]))) for x in CV.coeffs(c['G'])]
        zc, tc = CV.covariant(c['G'])
        r = SE.search_one(g, zc, tc, H)
        cl = {'M': [[CV.kstr(x) for x in row] for row in c['M']], 'd': [CV.kstr(x) for x in c['d']],
              'X': [[CV.kstr(Xk.monomial_coefficient(m)) for m in (CV.t1 ** 2, CV.t1 * CV.t2, CV.t2 ** 2)] for Xk in c['X']]}
        ninf = 0; best = None
        for t in sorted(set(SE.primitive(a, b) for a, b in r['hits'])):
            pt = CE.point_of_hit(cl, es, t)
            assert pt is not None
            if pt == 'O': continue
            P = E(list(pt))
            if P.order() == +Infinity:
                ninf += 1
                hp = (abs(complex(*t[0]) - zc * complex(*t[1])) ** 2 + tc ** 2 * abs(complex(*t[1])) ** 2) / tc
                if best is None or hp < best[0]: best = (hp, P)
        percoset.append(ninf)
        if best: pts.append(best[1])
        say('    coset: T = %.3g, %d pairs, %d hits of infinite order%s' % (tc, r['npairs'], ninf,
            '' if not best else '; smallest h_p(t) = %.1f, canonical height (Sage) %.3f' % (best[0], best[1].height())))
    k = 0
    if len(pts) >= 2:
        try:
            k, nh = KP.certify(es, [tuple(P.xy()) for P in pts[:2]], sel)
        except AssertionError:
            k = 0
    say('  %s: dim Sel^2 = %d; certified rank >= %d from the points found by the pipeline' % (name, sel.dim, k))
    return sel.dim, k

def main():
    rng = random.Random(20261127)
    say('controls.py (conductor-t13).  Written by an AI model (build agent t13); not read by a person.')
    say('(a) base changes y^2 = x(x + a)(x + b), rank 2 over K (kctp_points controls):')
    ok = 0; tot = 0
    for a, b in [(1, 7), (1, 11), (1, 22), (6, 7)]:
        es = [K(0), K(-a), K(-b)]
        d, k = run_curve(es, rng, 'a=%d b=%d' % (a, b))
        tot += 1; ok += (k >= 2)
    say('(b) y^2 = x(x - a)(x - b) over K with j not in Q and two independent points found by a naive search:')
    nb = 0
    cands = [(K(aa) + K(ab) * I, K(ba) + K(bb) * I) for aa, ab, ba, bb in itertools.product(range(-6, 7), repeat=4)]
    rng2 = random.Random(5)
    rng2.shuffle(cands)
    for A, B in cands:
        if nb >= 3: break
        if A == 0 or B == 0 or A == B: continue
        es = [K(0), A, B]
        E = R.curve_from_e(es)
        if E.j_invariant() in QQ: continue
        found = []
        for x0, x1 in itertools.product(range(-8, 9), repeat=2):
            x = K(x0) + K(x1) * I
            fx = x * (x - A) * (x - B)
            if fx != 0 and fx.is_square():
                P = E([x, fx.sqrt()])
                if P.order() == +Infinity: found.append(P)
        if len(found) < 2: continue
        sel = R.Selmer(es, rng)
        try:
            k0, _ = KP.certify(es, [tuple(P.xy()) for P in found[:6]], sel)
        except AssertionError:
            continue
        if sel.dim != 4: continue
        best = []
        for P in found:
            try:
                KP.certify(es, [tuple(Q.xy()) for Q in best + [P]], sel); best.append(P)
            except AssertionError:
                pass
        if len(best) != 2: continue
        d, k = run_curve(es, rng, 'a=%s b=%s' % (A, B))
        tot += 1; ok += (k >= 2); nb += 1
    say('controls passed (rank >= 2 certified from the pipeline\'s own points): %d of %d' % (ok, tot))
    if WRITE:
        open(os.path.join(HERE, 'controls.txt'), 'w').write('\n'.join(OUT) + '\n')

if __name__ == '__main__':
    main()
