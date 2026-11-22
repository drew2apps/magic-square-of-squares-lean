#!/usr/bin/env python3
"""Result 4: the canonical heights of the recorded generators, recomputed at higher precision.

For every line of dm_evidence.txt: build the model, lift the recorded abscissa of G to a point, and compute
its canonical height with Sage at 53 and at 200 bits. Reports the largest difference between the two and
between the 200-bit value and the recorded hG (four decimals). height_sides.py shows that the listed bound
on |a| leaves a margin of at least 9% in every line, so an error of this size in hG is harmless.
This is a stability check, not interval arithmetic: Sage does not return an error bound for the height.
Needs passagemath (sage.all__sagemath_schemes); a few seconds."""
import os, sys
from sage.all__sagemath_schemes import EllipticCurve, QQ

def model(curve, m, n, f):
    S = m*m + n*n
    if curve == 'E04': return [0, 2*S*S, 0, 4*S*S*m*m*n*n, 0]
    if curve == 'D1':
        c3 = 4*S*S*m*m*n*n; c2 = -16*S*m**3*n**3; c1 = 4*S*S*(m*m - n*n)**2; c0 = -16*S*m*n*(m*m - n*n)**2
        return [0, c1, 0, c0*c2, c0*c0*c3]
    return [int(x) for x in f['model'].split(',')]

here = os.path.dirname(os.path.abspath(__file__))
d1 = d2 = 0.0; count = 0
for line in open(os.path.join(here, 'dm_evidence.txt')):
    if line.startswith('#') or not line.strip(): continue
    w = line.split(); curve = w[2]; m, n = (int(x) for x in w[1].split(','))
    f = dict(x.split('=', 1) for x in w[3:] if '=' in x)
    E = EllipticCurve(QQ, model(curve, m, n, f)); G = E.lift_x(QQ(f['G']))
    h53 = G.height(precision=53); h200 = G.height(precision=200)
    d1 = max(d1, abs(float(h53 - h200)) / float(h200)); d2 = max(d2, abs(float(h200) - float(f['hG'])))
    count += 1
assert count == 122 and d1 < 1e-9 and d2 < 1e-4
print('lines: %d; largest relative difference between the 53-bit and the 200-bit height: %.1e;' % (count, d1))
print('largest difference between the 200-bit height and the recorded hG: %.1e' % d2)
print('height_precision OK')
