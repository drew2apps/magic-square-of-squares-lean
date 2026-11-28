#!/usr/bin/env python3
"""calib_torsion.py (conductor-t13): h_p of the preimages of the torsion points on the torsion-class coverings of
covers.json (the non-uniformity of the naive height). EXPERIMENT; written by an AI model (build agent t13), not read
by a person.  python3 calib_torsion.py > calib_torsion.txt  (seconds; floating point only in h_p, a guide)."""
from sage.all__sagemath_schemes import *
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'scripts', 'second_descent'))
import kctp_recheck as R
K = R.K; I = R.I
def kq(p): return K(QQ(p[0])) + K(QQ(p[1])) * I
cov = json.load(open(os.path.join(HERE, 'covers.json')))['curves']
for cv in ('E', 'E3'):
    es = [kq(x) for x in cov[cv]['e']]; E = R.curve_from_e(es)
    tors = [E(0)] + [E([e,0]) for e in es]
    for T in list(tors[1:]):
        tors += list(T.division_points(2))
    tors = list(set(tors))
    for cl in cov[cv]['classes']:
        if not cl['label'].startswith('torsion'): continue
        d1, d2, d3 = [kq(x) for x in cl['d']]
        (a, b), (c, d) = [[kq(x) for x in row] for row in cl['M']]
        Minv = matrix(K, [[a, b], [c, d]]).inverse()
        Xc = [[kq(x) for x in co] for co in cl['X']]
        zc = complex(cl['covariant'][0], cl['covariant'][1]); Tc = cl['covariant'][2]
        hps = []
        for P in tors:
            if P == 0: continue
            x = P.xy()[0]
            if x in es: continue
            r1, r2 = (x - es[0]) / d1, (x - es[1]) / d2
            if not (r1.is_square() and r2.is_square()): continue
            Rp = [r1.sqrt(), r2.sqrt(), K(1)]
            for s1 in (1, -1):
                for s2 in (1, -1):
                    Rv = [s1 * Rp[0], s2 * Rp[1], K(1)]
                    # find (t1 : t2) with X(t) proportional to Rv: brute force solve via resultants
                    S = PolynomialRing(K, 'w'); w = S.gen()
                    Xw = [co[0] * w ** 2 + co[1] * w + co[2] for co in Xc]
                    g = gcd(Xw[0] * Rv[2] - Xw[2] * Rv[0], Xw[1] * Rv[2] - Xw[2] * Rv[1])
                    sols = [(rt, K(1)) for rt, m in g.roots()] if g.degree() > 0 else []
                    for t in sols:
                        tf = Minv * vector(K, t)
                        den = lcm([q.denominator() for q in tf]); tf = tf * den
                        gg = K.ideal(list(tf)).gens_reduced()[0]; tf = tf / gg
                        t1c = complex(CC(tf[0])); t2c = complex(CC(tf[1]))
                        hps.append((abs(t1c - zc * t2c) ** 2 + Tc ** 2 * abs(t2c) ** 2) / Tc)
        print(cv, cl['label'], 'T =', Tc, 'h_p of torsion preimages:', sorted(set(round(h, 3) for h in hps)))
