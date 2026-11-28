"""Reviewer z13: recount of the numbers of PAIRING_FORMULA_2.md that come from the pilot (sections 4, 5, 6), per line
and per distinct case.  Written by an AI model (reviewer z13), not read by a person."""
import os, sys, json, collections
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'referee_z12'))
import zparse
R = zparse.load()
print('ratios', len(R), 'residual tuples', sum(len(v['tuples']) for v in R.values()))
killed = [(r, t) for r, v in R.items() for t in v['tuples'] if any(x == 'killed' for x in t['ver'].values())]
mixed_only = [(r, t) for r, t in killed if t['ver'].get('08') != 'killed' and t['ver'].get('LU') != 'killed']
print('tuples killed on some pair curve', len(killed), '; killed on mixed curves only', len(mixed_only))
print('tuples killed on 08 or LU', len(killed) - len(mixed_only))
# section 4 (two.py, two2.py): the kill lines on 08 / LU, per line and per distinct (ratio, curve, class, b)
from math import gcd
for cv in ('08', 'LU'):
    lines = [(r, k['cls'], k['b']) for r, v in R.items() for k in v['curves'].get(cv, {}).get('kills', [])]
    def Aodd(r):
        m, n = r
        A = m * m - n * n if cv == '08' else 2 * m * n
        C = m * m + n * n
        while A % 2 == 0 and C % 2 == 0: A //= 2; C //= 2
        return A % 2 == 1
    ev = [x for x in lines if not Aodd(x[0])]; od = [x for x in lines if Aodd(x[0])]
    print(cv, 'kill lines', len(lines), 'A even: lines', len(ev), 'distinct', len(set(ev)), '| A odd: lines', len(od), 'distinct', len(set(od)))
# section 5: the exceptional tuples
for r in ((5623, 1017), (3320, 2303)):
    m, n = r
    for t in R[r]['tuples']:
        if not any(x == 'killed' for x in t['ver'].values()):
            d0, d8, dL, dU = t['tup']
            print(r, t['tup'], 'd0*d8 =', d0 * d8, 'dL*dU =', dL * dU, t['ver'])
