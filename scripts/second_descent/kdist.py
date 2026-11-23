#!/usr/bin/env python3
import sys, itertools, json
"""The 2-Selmer bounds over K = Q(i) of all 210 quartics of kclose.py for one ratio, not only whether they are 0.  usage: kdist.py m,n"""
from kclose import *
def dist(m, n):
    S = m*m + n*n; Ssq = isqrt(S)**2 == S
    F = factors(m, n); names = sorted(F); roots = {k: -F[k][0]/F[k][1] for k in names}; out = []
    for J in itertools.combinations(names, 4):
        Q = prod(F[k] for k in J)
        vals = [Q(-I), Q(I)] if Ssq else [Q(0), Q.leading_coefficient()]
        if any(v == 0 for v in vals): continue
        deltas = []
        for v in vals:
            d = sqfree_rep(v)
            if not any((d/x).is_square() for x in deltas): deltas.append(d)
        r = roots[J[3]]; bs = []
        for d in deltas:
            E, A = curve(d*Q, r); xs = [A/(roots[k] - r) for k in J[:3]]
            den = lcm([x[j].denominator() for x in xs for j in (0, 1)]); e = [gint(x*den*den) for x in xs]
            b = kd.selmer_bound(e); bs.append(b if isinstance(b, int) or b is None else b[0])
        out.append((''.join(J), bs)); print(''.join(J), bs, flush=True)
    return out
m, n = (int(v) for v in sys.argv[1].split(','))
o = dist(m, n)
from collections import Counter
c = Counter(max((b if b is not None else 99) for b in bs) for J, bs in o)
print('(%d,%d): quartics by the larger of the Selmer bounds over K of their curves:' % (m, n), dict(sorted(c.items())))
