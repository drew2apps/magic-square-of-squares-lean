#!/usr/bin/env python3
"""z11: try to make MY method kill the class of a rational point (it must never).  Written by an AI model (review agent
z11); not read by a person.  For every pair curve of (7,1), (4,3) (degenerate, rank 3 on the mixed curves) and
(465,416): every subset sum of PARI's points (ellrank and elltors, FINDERS; every point checked on the curve), its
Kummer class, and the same class recomputed by the divisor formula on C (both signs): all must lie in my Sel^2 and pair
to +1 with all of Sel^2.  PARI's random state is fixed (setrand(1)).
    python3 -I forge_mine.py      (prints forge_mine.txt)"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _z11common import Out
from z11_pair import PairCurve, Pairing, PAIRS, Fr, sqf, rclass, pari
say = Out('forge_mine.txt')
say('forge_mine.py (review z11; written by an AI model): classes of rational points through my pairing')
def add(C, P, Q):
    if P is None: return Q
    if Q is None: return P
    A2, A4, A6 = C.cub
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2:
        if y1 + y2 == 0: return None
        l = (3 * x1 * x1 + 2 * A2 * x1 + A4) / (2 * y1)
    else: l = (y2 - y1) / (x2 - x1)
    x3 = l * l - A2 - x1 - x2; return (x3, -(y1 + l * (x3 - x1)))
allok = True
for m, n in [(7, 1), (4, 3), (465, 416)]:
    for jk in PAIRS:
        C = PairCurve(m, n, jk); role = C.ffield_identify([101, 109, 113]); PR = Pairing(C.es)
        E = pari.ellinit([0, C.cub[0], 0, C.cub[1], C.cub[2]])
        pari.setrand(1); r = pari.ellrank(E); tor = pari.elltors(E)
        G = [(Fr(str(P[0])), Fr(str(P[1]))) for P in list(r[3]) + list(tor[2])]
        for X, Y in G: assert Y * Y == (X - C.es[0]) * (X - C.es[1]) * (X - C.es[2])
        combos = [None]
        for g in G: combos = combos + [add(C, P, g) for P in combos]
        cls = set(); viaC = 0
        for P in combos:
            if P is None or P[1] == 0: continue
            X, Y = P
            cls.add((sqf(X - C.es[0]), sqf(X - C.es[1])))
            tv = C.from_X(X, Y)
            if tv is not None:
                t, v = tv; parts = C.parts(t); P0 = C.parts(C.t0)
                comp = [None] * 3
                comp[role['r']] = sqf(parts['r'] * P0['r'])
                for sg in (1, -1):
                    comp[role['*']] = sg * sqf(rclass(parts['*']) * rclass(P0['*']))
                    comp[role["'"]] = sqf(comp[role['r']] * comp[role['*']])
                    cls.add((comp[0], comp[1])); viaC += 1
        res = {c: PR.test(c)[0] for c in sorted(cls)}
        bad = sorted(c for c, v in res.items() if v != 'kernel')
        allok &= not bad
        say('(%d,%d) E_%s: dim Sel^2 %d; %d points found, %d subset sums; %d distinct classes tested (%d from the divisor formula): %s'
            % (m, n, jk, PR.dim, len(G), len(combos), len(res), viaC, 'all in the kernel' if not bad else 'NOT IN THE KERNEL: %s' % bad))
say('NO CLASS OF A RATIONAL POINT IS KILLED' if allok else 'A CLASS OF A RATIONAL POINT IS KILLED: THE METHOD IS WRONG')
say.finish()
