#!/usr/bin/env python3
"""Result 4: the two sides of Silverman's height bound.

dm.py, dm2.py and dm3.py use |h(x(Q)) - hhat(Q)| <= B for every rational point Q, with B = Sage's
silverman_height_bound() (computed by eclib) and hhat = Sage's canonical height (twice Silverman's). eclib
documents B as a bound for h - hhat only. Silverman (Math. Comp. 55 (1990), Theorem 1.1): for a Weierstrass
equation with integral coefficients, discriminant D, invariant j and b2 = a1^2 + 4 a2, put
    mu = log|D|/12 + log+|j|/12 + log+|b2/12|/2 + (log 2)/2 [if b2 != 0];
then -h(j)/24 - mu - 0.973 <= hhat_S(Q) - h(x(Q))/2 <= mu + 1.07, where hhat_S = hhat/2. So
    h - hhat <= h(j)/12 + 2 mu + 1.946   (eclib returns this with 1.922, Bremner's correction: libsrc/htconst.cc),
    hhat - h <= 2 mu + 2.14 = B - h(j)/12 + 0.218.
Hence B bounds both sides when h(j) >= 2.616. For the model of every line of dm_evidence.txt this script
  1. recomputes B from eclib's formula, exactly from the coefficients, and compares it with the recorded B
     (this also shows that B was computed for the model whose points the scripts list);
  2. checks U = 2 mu + 2.14 <= B.
No dependencies. The statement of Silverman's theorem above is quoted from memory and from eclib's source;
it should be compared with the paper."""
import os, sys, math
from fractions import Fraction as Fr

def invariants(a1, a2, a3, a4, a6):
    b2 = a1*a1 + 4*a2; b4 = 2*a4 + a1*a3; b6 = a3*a3 + 4*a6
    b8 = a1*a1*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3*a3 - a4*a4
    c4 = b2*b2 - 24*b4
    D = -b2*b2*b8 - 8*b4**3 - 27*b6*b6 + 9*b2*b4*b6
    return c4, D

def model(curve, m, n, fields):
    S = m*m + n*n
    if curve == 'E04': return (0, 2*S*S, 0, 4*S*S*m*m*n*n, 0)
    if curve == 'D1':
        c3 = 4*S*S*m*m*n*n; c2 = -16*S*m**3*n**3; c1 = 4*S*S*(m*m - n*n)**2; c0 = -16*S*m*n*(m*m - n*n)**2
        return (0, c1, 0, c0*c2, c0*c0*c3)
    for f in fields:
        if f.startswith('model='): return tuple(int(x) for x in f[6:].split(','))
    raise ValueError(curve)

here = os.path.dirname(os.path.abspath(__file__))
worst = 0.0; count = 0; bad = []; dmax = 0.0; hjmin = 1e9
for line in open(os.path.join(here, 'dm_evidence.txt')):
    if line.startswith('#') or not line.strip(): continue
    w = line.split(); run = w[1]; curve = w[2]
    m, n = (int(x) for x in run.split(','))
    B = float([f for f in w if f.startswith('B=')][0][2:])
    a = model(curve, m, n, w)
    assert all(isinstance(x, int) for x in a)
    c4, D = invariants(*a); assert D != 0
    j = Fr(c4**3, D)
    hj = max(math.log(abs(j.numerator)), math.log(j.denominator)) if j != 0 else 0.0
    logplus = lambda x: math.log(x) if x > 1 else 0.0
    b2 = a[0]*a[0] + 4*a[1]
    twomu = math.log(abs(D))/6 + logplus(abs(j.numerator)/j.denominator)/6 + logplus(abs(b2)/12) + (math.log(2) if b2 else 0.0)
    Bf = 1.922 + hj/12 + twomu
    U = twomu + 2.14
    count += 1; worst = max(worst, U / B); dmax = max(dmax, abs(Bf - B)); hjmin = min(hjmin, hj)
    if U > B or abs(Bf - B) > 0.002: bad.append((w[0], curve, round(U, 3), round(Bf, 4), B))
    if '-v' in sys.argv: print(w[0], curve, 'B=%.3f' % B, 'formula=%.4f' % Bf, 'U=%.3f' % U, 'h(j)=%.2f' % hj)
print('lines: %d; largest |formula - recorded B|: %.4f; smallest h(j): %.2f; largest U/B: %.4f; lines failing: %d' % (count, dmax, hjmin, worst, len(bad)))
for b in bad: print('  ', b)
assert count == 121 and not bad
print('height_sides OK')
