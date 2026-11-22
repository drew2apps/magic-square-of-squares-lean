#!/usr/bin/env python3
"""Second software for the rank bound of the E04 lines of dm_evidence.txt.

The bound 'rank <= 1' in dm_evidence.txt comes from PARI's ellrank (2-descent with the Cassels pairing),
applied to the curve and to the curves isogenous to it. This script asks eclib (mwrank, an independent
2-descent without the Cassels pairing) for its bound on the same curves and prints the lines for which
eclib alone does not reach 1. For those lines the bound rests on PARI's Cassels-pairing step.
The D1 and E0L curves are not treated: eclib stops on them with an NTL error.
Needs passagemath (sage.all__sagemath_schemes, passagemath-eclib) and cypari2."""
import os
from sage.all__sagemath_schemes import *
import cypari2
pari = cypari2.Pari()

def eclib_bounds(E):
    C = E.mwrank_curve(verbose=False); C.two_descent(verbose=False, second_limit=8)
    return int(C.rank()), int(C.rank_bound())

here = os.path.dirname(os.path.abspath(__file__))
lines = [l.split() for l in open(os.path.join(here, 'dm_evidence.txt')) if l.strip() and not l.startswith('#')]
direct = iso = 0; left = []; seen = {}
for r in lines:
    if r[2] != 'E04': continue
    if r[1] not in seen:
        m, n = (int(x) for x in r[1].split(',')); S = m*m + n*n
        E = EllipticCurve(QQ, [0, 2*S*S, 0, 4*S*S*m*m*n*n, 0]).minimal_model()
        lo, hi = eclib_bounds(E); how = 'direct'
        if hi > 1:
            for c in pari.ellisomat(E.pari_curve(), 0, 1)[0]:
                F = EllipticCurve(QQ, [QQ(str(x)) for x in list(c)]).minimal_model()
                lo2, hi2 = eclib_bounds(F); lo = max(lo, lo2); hi = min(hi, hi2)
            how = 'isogenous'
        seen[r[1]] = (lo, hi, how)
    lo, hi, how = seen[r[1]]
    if hi == 1 and how == 'direct': direct += 1
    elif hi == 1: iso += 1
    else: left.append('%s (run %s): eclib [%d, %d]' % (r[0], r[1], lo, hi))
print('E04 lines: %d; eclib gives rank <= 1 on the curve: %d; on an isogenous curve: %d; not reached: %d' % (direct + iso + len(left), direct, iso, len(left)))
for l in left: print('  ' + l)
