#!/usr/bin/env python3
"""certify.py (conductor-t13): map the exact hits of search.py back to E(K) (or E3(K) and then to E(K) by the dual
isogeny), re-check every point exactly, and certify independence with kctp_points.certify.  EXPERIMENT; written by an
AI model (build agent t13), not read by a person.

    python3 certify.py search_E_H50000.json [search_E3_H50000.json ...] [--write]   -> certify.txt, hints_t13.json

For a hit t = (t1, t2) of the class with data (d1, d2, d3, X, M, s) (covers.json): (u : v : z) = X(M t),
x = e1 + d1 (u/z)^2; the point is accepted only if f(x) = (x - e1)(x - e2)(x - e3) is a square in K (exact), and then
(x, sqrt f(x)) is checked on the curve and its Kummer image (x - e1, x - e2) is checked to be the class (d1, d2) of the
covering modulo squares (exact: S-unit vector).  z = 0 gives the points above O.  Torsion is decided exactly
(Sage's order, and the certifier of kctp_points, which needs Kummer images independent of the torsion image).
"""
from sage.all__sagemath_schemes import *
import sys, os, json, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'scripts', 'second_descent'))
import kctp_recheck as R
import kctp_points as KP
import kctp_isog as KI

K = R.K; I = R.I
WRITE = '--write' in sys.argv
OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(str(s))

def kq(p): return K(QQ(p[0])) + K(QQ(p[1])) * I

def point_of_hit(cl, es, t):
    (a, b), (c, d) = [[kq(x) for x in row] for row in cl['M']]
    T1 = K(t[0][0]) + K(t[0][1]) * I; T2 = K(t[1][0]) + K(t[1][1]) * I
    o1, o2 = a * T1 + b * T2, c * T1 + d * T2
    X = [kq(co[0]) * o1 ** 2 + kq(co[1]) * o1 * o2 + kq(co[2]) * o2 ** 2 for co in cl['X']]
    u, v, z = X
    d1 = kq(cl['d'][0])
    if z == 0:
        return 'O'
    x = es[0] + d1 * (u / z) ** 2
    fx = (x - es[0]) * (x - es[1]) * (x - es[2])
    if fx == 0: return (x, K(0))
    if not fx.is_square(): return None
    return (x, fx.sqrt())

def to_E_map():
    """the curves E, E3 (models of covers.py) and the dual isogeny E3 -> E, exact (Sage Velu + isomorphisms)."""
    Q, fs = R.quartic(627, 364, 'L1L3M1M2'); es = R.full2_model(Q, fs)
    E = R.curve_from_e(es)
    phi = E.isogeny(E([es[2], 0])); E3s = phi.codomain()
    es3 = KI.short_full2(E3s); E3 = R.curve_from_e(es3)
    iso = E3.isomorphism_to(E3s)
    psi = phi.dual()
    back = psi.codomain().isomorphism_to(E)
    def m(P):
        Q_ = back(psi(iso(P))); assert Q_ in E; return Q_
    # check on the points of order 4: psi o phi = [2] (up to the sign of the isomorphisms)
    hs = E([es[2], 0]).division_points(2); assert hs
    for H in hs:
        assert m((~iso)(phi(H))) in (2 * H, -2 * H)
    return E, es, E3, es3, m

def main():
    files = [a for a in sys.argv[1:] if not a.startswith('--')]
    rng = random.Random(20261127)
    cov = json.load(open(os.path.join(HERE, 'covers.json')))['curves']
    if os.path.exists(os.path.join(HERE, 'altmodels.json')):
        alt = json.load(open(os.path.join(HERE, 'altmodels.json')))['curves']
        for cv in cov: cov[cv]['classes'] = cov[cv]['classes'] + alt[cv]['classes']
    E, es, E3, es3, m = to_E_map()
    say('certify.py (conductor-t13).  Written by an AI model (build agent t13); not read by a person.')
    sel = R.Selmer(es, rng)
    found = []
    for fn in files:
        res = json.load(open(os.path.join(HERE, fn)))
        cv = 'E3' if '_E3_' in fn else 'E'
        ee, C = (es3, E3) if cv == 'E3' else (es, E)
        for cl in cov[cv]['classes']:
            if cl['label'] not in res: continue
            hits = res[cl['label']]['hits']
            nO = nt = ninf = 0
            for t in hits:
                pt = point_of_hit(cl, ee, t)
                assert pt is not None, 'a hit does not give a point of the curve: %s' % (t,)
                if pt == 'O':
                    nO += 1; continue
                x, y = pt
                assert y ** 2 == (x - ee[0]) * (x - ee[1]) * (x - ee[2])
                P = C([x, y])
                d1, d2 = R.kummer(x, ee)
                if cv == 'E':
                    v = sel.vec_of(d1 / kq(cl['d'][0])); w = sel.vec_of(d2 / kq(cl['d'][1]))
                    assert v is not None and w is not None and not any(v) and not any(w), 'Kummer image is not the class'
                PE = m(P) if cv == 'E3' else P
                if P.order() != +Infinity:
                    nt += 1
                else:
                    ninf += 1; found.append((cv, cl['label'], t, PE))
            say('  %s %s (%s): %d hits: %d above O, %d torsion, %d of infinite order' % (cv, cl['label'], fn, len(hits), nO, nt, ninf))
    say('points of infinite order found: %d' % len(found))
    if found:
        pts = [tuple(P.xy()) for _, _, _, P in found]
        try:
            k, nh = KP.certify(es, pts, sel)
            say('certified: rank E(K) >= %d (kctp_points.certify; %d halvings)' % (k, nh))
        except AssertionError as ex:
            # certify a maximal independent subset
            best = []
            for p in pts:
                try:
                    KP.certify(es, best + [p], sel); best.append(p)
                except AssertionError:
                    pass
            say('certified: rank E(K) >= %d (a maximal independent subset of the points found)' % len(best))
    if WRITE:
        json.dump({'note': 'Written by an AI model (build agent t13), not read by a person. Points of E(K) of infinite order '
                           'found by search.py (exact coordinates as strings in K = Q(i)); empty if none was found.',
                   'files': files, 'points_on_E': [[str(c) for c in tuple(P.xy())] for _, _, _, P in found]},
                  open(os.path.join(HERE, 'hints_t13.json'), 'w'), indent=1, sort_keys=True)
        open(os.path.join(HERE, 'certify.txt'), 'w').write('\n'.join(OUT) + '\n')

if __name__ == '__main__':
    main()
