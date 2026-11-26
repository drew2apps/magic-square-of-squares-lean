#!/usr/bin/env python3
"""kctp_points.py: the rank-2 route for the corner curve L1L3M1M2 of (627,364) over K = Q(i): an exact certifier for
points of E(K) (on the curve, of infinite order, independent), with controls of known rank 2 and forged inputs.
EXPERIMENT; written by an AI model (Claude), not read by a person.

    python3 kctp_points.py [--write]       -> kctp_points.txt

THE CERTIFIER.  For E: y^2 = (x - e1)(x - e2)(x - e3) over K (all 2-torsion over K) and points P_1, ..., P_k given by
exact coordinates in K: each is checked on the curve exactly; the 2-descent (Kummer) map
P -> (x - e1, x - e2) mod squares (the usual formulas at T1, T2) is linear from E(K)/2E(K) to K(S,2)^2.  If the images
of T1, T2, the points of order 4 (generating the image of E(K)_tors) and P_1..P_k span a space of dimension
dim(torsion image) + k, the P_j are independent in E(K)/E(K)_tors: a relation sum n_j P_j in E(K)_tors with some n_j
odd (divide a relation by 2 as long as all n_j are even: 2 Q torsion implies Q torsion) would give a relation among the
images.  So rank E(K) >= k (exact).  A point not on the curve, a torsion point, or a dependent set is refused.

CONTROLS.  E_ab: y^2 = x(x + a)(x + b) over Q with rank E_ab(Q) = rank E_ab^-1(Q) = 1 (PARI ellrank, a guide for the
ranks over Q; the points it returns are checked exactly here).  A point (X, Y) of E^-1: Y^2 = X(X - a)(X - b) gives the
point (-X, i Y) of E(K).  The certifier must give rank E(K) >= 2 on each control.

THE TARGET.  No point of infinite order is known on E(K).  kctp.py --search=60 (run 3) searched the 2-coverings with
|re y|, |im y| <= 60 and found none; this script adds the certifier only.  Why a larger search was not run: kctp_isog.py
(guide) puts the BSD scale at |Sha| ~ 4.6e5 L(E/K,1) (rank 0) or Reg |Sha| ~ 4.6e5 L''(E/K,1)/2 (rank 2), i.e. generators
of canonical height of order 10^3 if the rank is 2 and Sha is small -- beyond 2- and 4-covering searches.
"""
from sage.all__sagemath_schemes import *
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kctp_recheck as R

K = R.K; I = R.I
WRITE = '--write' in sys.argv
OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(str(s))

def certify(es, pts, sel, max_halvings=12):
    """exact: points on y^2 = prod (x - e_k); returns the number k of points certified independent modulo torsion
    (rank E(K) >= k).  If the images in E(K)/2E(K) are dependent, the dependent combination Q lies in 2E(K) (the Kummer
    map is injective on E(K)/2E(K)); Q is halved exactly (2R = Q checked) and R replaces a point of the combination
    (2-saturation).  Raises AssertionError on a point off the curve, or when no independent set is reached within
    max_halvings halvings (a dependent or torsion input is never certified: certification needs independent images)."""
    E = R.curve_from_e(es)
    for x, y in pts:
        assert y ** 2 == (x - es[0]) * (x - es[1]) * (x - es[2]), 'point not on the curve'
    P = [E([x, y]) for x, y in pts]
    tgen = [E([es[0], 0]), E([es[1], 0])]
    for t in [E([ee, 0]) for ee in es]:
        try: hs = t.division_points(2)
        except Exception: hs = []
        for H in hs:
            assert 2 * H == t; tgen.append(H)
    def img(Q):
        if Q == 0: return [0] * (2 * sel.n)
        d1, d2 = R.kummer(Q.xy()[0], es); v1, v2 = sel.vec_of(d1), sel.vec_of(d2)
        assert v1 is not None and v2 is not None and sel.contains(d1, d2), 'Kummer image not in Sel^2'
        return v1 + v2
    tor = [img(t) for t in tgen]; rt = R.f2_rank(tor)
    for step in range(max_halvings + 1):
        rows = [img(Q) for Q in P]
        if R.f2_rank(tor + rows) == rt + len(P):
            return len(P), step
        if step == max_halvings: break
        # a dependency: sum eps_j img(P_j) + sum del_t img(T_t) = 0 with some eps_j = 1
        Mx = matrix(GF(2), rows + tor)
        ker = Mx.left_kernel().basis()
        vec = [v for v in ker if any(v[j] for j in range(len(P)))][0]
        Q = sum([P[j] for j in range(len(P)) if vec[j]], E(0)) + sum([tgen[t] for t in range(len(tgen)) if vec[len(P) + t]], E(0))
        hs = Q.division_points(2)
        assert hs, 'a combination with trivial Kummer image is not in 2E(K): impossible, the Kummer map is injective'
        Rh = hs[0]; assert 2 * Rh == Q
        j = [j for j in range(len(P)) if vec[j]][0]
        P[j] = Rh
    raise AssertionError('the points are not certified independent modulo torsion (after %d halvings)' % max_halvings)

def controls(rng):
    import cypari2
    pari = cypari2.Pari()
    say('Controls: y^2 = x(x + a)(x + b), rank 1 over Q for E and for E^-1 (PARI ellrank, guide), so rank E(K) = 2.')
    found = 0
    for a, b in [(1, 7), (2, 7), (1, 11), (3, 7), (2, 13), (5, 7), (1, 14), (3, 10), (2, 9), (4, 13), (1, 22), (6, 7)]:
        rk = []; gens = []
        for sg in (1, -1):
            Eq = pari.ellinit([0, sg * (a + b), 0, a * b, 0])
            r = pari.ellrank(Eq)
            rk.append((int(r[0]), int(r[1]))); gens.append([(QQ(P[0]), QQ(P[1])) for P in r[3]])
        if rk != [(1, 1), (1, 1)] or not all(gens): continue
        es = [K(0), K(-a), K(-b)]
        sel = R.Selmer(es, rng)
        (X1, Y1), (X2, Y2) = gens[0][0], gens[1][0]
        P1 = (K(X1), K(Y1)); P2 = (-K(X2), I * K(Y2))
        k, nh = certify(es, [P1, P2], sel)
        say('  a=%d b=%d: P1 = %s (from E(Q)), P2 = %s (from E^-1(Q)); dim Sel^2_K = %d; certified rank E(K) >= %d (%d halvings)' % (
            a, b, P1, P2, sel.dim, k, nh))
        assert k == 2 and sel.dim >= 4
        found += 1
        # negative tests on this control
        neg = 0
        for name, pts in [('P1 with y + 1 (off the curve)', [P1, (P1[0], P1[1] + 1)]),
                          ('P1, P1 + T1 (dependent mod torsion)', [P1, (lambda Q: (Q[0], Q[1]))(tuple((R.curve_from_e(es)(list(P1)) + R.curve_from_e(es)([es[0], 0])).xy()))]),
                          ('P1, 3 P1 (dependent)', [P1, tuple((3 * R.curve_from_e(es)(list(P1))).xy())]),
                          ('P1, P2, P1 + P2 (dependent)', [P1, P2, tuple((R.curve_from_e(es)(list(P1)) + R.curve_from_e(es)(list(P2))).xy())]),
                          ('a torsion point (e1, 0)', [(es[0], K(0))])]:
            try:
                certify(es, pts, sel); say('    NOT REFUSED: %s' % name)
            except AssertionError as ex:
                neg += 1
        say('    forged inputs refused: %d of 5' % neg)
        assert neg == 5
        if found >= 4: break
    assert found >= 3
    return found

def main():
    rng = random.Random(20261127)
    say('kctp_points.py: exact certifier for points of E(K), K = Q(i); controls; the target.  Written by an AI model;')
    say('not read by a person.  Experiment, not counted.')
    n = controls(rng)
    say('%d controls of rank 2 over K certified (rank >= 2, exact); every forged input refused.' % n)
    say()
    Q, fs = R.quartic(627, 364, 'L1L3M1M2'); es = R.full2_model(Q, fs)
    sel = R.Selmer(es, rng)
    E = R.curve_from_e(es)
    # the target: the torsion points must be refused as points of infinite order; a forged point too
    neg = 0
    for name, pts in [('(e1, 0)', [(es[0], K(0))]),
                      ('a point of order 4', [tuple(E([es[2], 0]).division_points(2)[0].xy())]),
                      ('(e1 + 1, 1) (off the curve)', [(es[0] + 1, K(1))])]:
        try:
            certify(es, pts, sel); say('  NOT REFUSED on the target: %s' % name)
        except AssertionError:
            neg += 1
    say('Target L1L3M1M2 of (627,364): no point of infinite order known; the certifier refuses %d of 3 forged/torsion' % neg)
    say('inputs.  The search of run 3 (kctp.py --search=60: 2-coverings, |re y|, |im y| <= 60) found none.')
    assert neg == 3

if __name__ == '__main__':
    main()
    if WRITE:
        open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'kctp_points.txt'), 'w').write('\n'.join(OUT) + '\n')
