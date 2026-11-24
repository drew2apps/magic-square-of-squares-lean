#!/usr/bin/env python3
"""Control: real Heron triangles with two rational medians (Buchholz-Rathbun). Base AB = the side whose median is NOT rational, so that
m_a, m_b (forms z+1, z-2) are rational: forms 0..3 are squares; the 4-form first step must keep the tuple of the real point."""
import sys, os
from math import gcd, isqrt
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'cover_descent'))
import pt, tuple_descent as td
issq = lambda x: x >= 0 and isqrt(x)**2 == x
tris = [(73, 51, 26), (626, 875, 291), (4368, 1241, 3673), (14791, 14384, 11257), (28779, 13816, 15155)]
for tri in tris:
    med = {x: issq(2*sum(y*y for y in tri) - 3*x*x) for x in tri}     # 4 m_x^2 = 2(sum) - 3 x^2
    c = [x for x in tri if not med[x]]
    A16 = (sum(tri))*(sum(tri) - 2*tri[0])*(sum(tri) - 2*tri[1])*(sum(tri) - 2*tri[2])
    print('triangle', tri, 'rational median to side:', med, '16*Area^2 square:', issq(A16))
    if len(c) != 1: continue
    c = c[0]; a, b = [x for x in tri if x != c]                      # a = |BC|, b = |AC|
    # C = z*c with x = (b^2 + c^2 - a^2)/(2c), y = 2*Area/c
    x = F(b*b + c*c - a*a, 2*c*c); y = F(isqrt(A16), 2*c*c)
    r, s = y.numerator, y.denominator; u, v = x.numerator, x.denominator
    forms, sup, known = pt.triangle(r, s)
    vals = [(al[0]*u + be[0]*v, al[1]*u + be[1]*v) for al, be in forms]
    sq = [issq(p*p + q*q) for p, q in vals]
    tup = tuple(td.cls_of(vals[j]) for j in range(4))
    T4, used, sizes = td.descend(forms[:4], sup[:4], [k[:4] for k in known], 150)
    T5, known5 = pt.line(r, s)
    print('  z = %s + %s i; norms square (z, z-1, z+1, z-2, 2z-1): %s' % (x, y, sq))
    print('  tuple of the real point (4 forms): %s; kept by the 4-form first step: %s (of %d left; candidates %d)' % (tup, tup in T4, len(T4), sizes[0]**4))
    print('  squarefree part of d1 d2 d3 d4: %d' % td.sqf(tup[0]*tup[1]*tup[2]*tup[3]))
    print('  5-form first step on this line: %d left: %s' % (len(T5), T5[:8]))
