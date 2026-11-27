#!/usr/bin/env python3
"""z11: symbolic check of the map of Washington (Elliptic Curves, Thm 2.17) used by z11_pair.py and t10_pairct.py.
Written by an AI model (review agent z11); not read by a person.  sympy only, no PARI.
For v^2 = a s^4 + b s^3 + c s^2 + d s + q^2:  x = (2q(v+q) + d s)/s^2,  y = (4q^2(v+q) + 2q(d s + c s^2) - d^2 s^2/(2q))/s^3
lie on y^2 + a1 x y + a3 y = x^3 + a2 x^2 + a4 x + a6 (a1 = d/q, a2 = c - d^2/(4q^2), a3 = 2qb, a4 = -4q^2 a, a6 = a2 a4),
with inverse s = (2q(x + c) - d^2/(2q))/y, v = -q + s(s x - d)/(2q); residues are computed modulo v^2 - g(s).
    python3 -I wash_check.py      (prints wash_check.txt)"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _z11common import Out
import sympy as sp
say = Out('wash_check.txt')
a, b, c, d, q, s, v = sp.symbols('a b c d q s v')
g = a*s**4 + b*s**3 + c*s**2 + d*s + q**2
x = (2*q*(v+q) + d*s)/s**2
y = (4*q**2*(v+q) + 2*q*(d*s + c*s**2) - d**2*s**2/(2*q))/s**3
a1 = d/q; a2 = c - d**2/(4*q**2); a3 = 2*q*b; a4 = -4*q**2*a; a6 = a2*a4
def res(e): return sp.simplify(sp.rem(sp.expand(sp.numer(sp.together(e))), v**2 - g, v))
r1 = res(y**2 + a1*x*y + a3*y - (x**3 + a2*x**2 + a4*x + a6))
r2 = res((2*q*(x + c) - d**2/(2*q))/y - s)
r3 = res(-q + s*(s*x - d)/(2*q) - v)
say('wash_check.py (review z11; written by an AI model): Washington Thm 2.17, symbolically, modulo v^2 = g(s)')
say('  the Weierstrass equation: residue %s' % r1)
say('  inverse, s: residue %s' % r2)
say('  inverse, v: residue %s' % r3)
say('ALL ZERO' if (r1, r2, r3) == (0, 0, 0) else 'NOT ZERO')
say.finish()
