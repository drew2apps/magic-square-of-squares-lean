#!/usr/bin/env python3
"""Proposition 1 of RATIONAL_DISTANCE.md, every computation in its proof, in exact arithmetic over Q(i).

    C:   w3^2 = a w1^2 + b w2^2,   w4^2 = c w1^2 + d w2^2,      a = 1 - i, b = i, c = -i, d = 1 + i.

(1) On C, X = (w1/w2)^2 and Y = w1 w3 w4 / w2^3 satisfy  E':  Y^2 = X (aX + b)(cX + d).
(2) With kappa = -(1 + i), beta = 1 - i:  x = kappa (X - beta/2), y = kappa Y  turns E' into  E: y^2 = x^3 - x.
(3) E has the eight points O, (0,0), (1,0), (-1,0), (i, +-(1 - i)), (-i, +-(1 + i)) over Q(i), and 8 points over F_5 and
    over F_13 (the residue fields of Q(i) at the primes above 5 and 13).
(4) The X-coordinates of the eight points are oo, beta/2, 0, beta, -i, 1.  X = (w1/w2)^2 has to be a square of Q(i):
    of the five finite values only 0 and 1 are squares, and X = 0, X = oo need i or 1 - i to be a square.
No Jacobian and no computer algebra system is used; sympy only expands polynomials.  Prints one line per step and
'PROP1 OK' at the end; any failure is an AssertionError."""
from fractions import Fraction as F
from math import isqrt
from sympy import symbols, I, expand, Rational

# (1) and (2): polynomial identities
a, b, c, d = 1 - I, I, -I, 1 + I
w1, w2, X = symbols('w1 w2 X')
W3, W4 = a*w1**2 + b*w2**2, c*w1**2 + d*w2**2                      # w3^2 and w4^2 on C
lhs = w1**2*W3*W4/w2**6                                            # Y^2
Xv = w1**2/w2**2
assert expand(lhs*w2**6 - (Xv*(a*Xv + b)*(c*Xv + d))*w2**6) == 0
print('(1) Y^2 = X(aX + b)(cX + d) on C: identity checked')
kappa, beta = -(1 + I), 1 - I
x = kappa*(X - beta/2)
assert expand(kappa**2*X*(a*X + b)*(c*X + d) - (x**3 - x)) == 0
print('(2) kappa^2 X(aX + b)(cX + d) = x^3 - x with x = kappa(X - beta/2), kappa = -(1 + i), beta = 1 - i: identity checked')

# exact arithmetic in Q(i): pairs of Fractions
def mul(p, q): return (p[0]*q[0] - p[1]*q[1], p[0]*q[1] + p[1]*q[0])
def sub(p, q): return (p[0] - q[0], p[1] - q[1])
def rat_sqrt(q):
    if q < 0: return None
    n, m = isqrt(q.numerator), isqrt(q.denominator)
    return F(n, m) if n*n == q.numerator and m*m == q.denominator else None
def is_square(g):
    """Is g = (re, im) a square in Q(i)?  (p + qi)^2 = g  iff  p^2 - q^2 = re, 2pq = im; then p^2 = (re +- n)/2, n^2 = norm."""
    if g == (0, 0): return True
    n = rat_sqrt(g[0]*g[0] + g[1]*g[1])
    if n is None: return False
    for s in (n, -n):
        p = rat_sqrt((g[0] + s)/2)
        if p is None: continue
        if p == 0:
            q = rat_sqrt(-g[0])
            if q is not None and g[1] == 0: return True
            continue
        q = g[1]/(2*p)
        if p*p - q*q == g[0]: return True
    return False
assert is_square((F(0), F(2))) and is_square((F(-1), F(0))) and not is_square((F(0), F(1))) and not is_square((F(2), F(0)))

# (3) the eight points
one, i_ = (F(1), F(0)), (F(0), F(1))
pts = [((F(0), F(0)), (F(0), F(0))), (one, (F(0), F(0))), ((F(-1), F(0)), (F(0), F(0))),
       (i_, (F(1), F(-1))), (i_, (F(-1), F(1))), ((F(0), F(-1)), (F(1), F(1))), ((F(0), F(-1)), (F(-1), F(-1)))]
for xx, yy in pts: assert mul(yy, yy) == sub(mul(xx, mul(xx, xx)), xx)
print('(3) seven affine points of y^2 = x^3 - x over Q(i) checked; with O that is 8')
for p in (5, 13):
    n = 1 + sum(1 for xx in range(p) for yy in range(p) if (yy*yy - xx**3 + xx) % p == 0)
    assert n == 8
    print('    points of y^2 = x^3 - x over F_%d: %d' % (p, n))

# (4) pull back: X = beta/2 + x/kappa = (beta/2)(1 - x)
half_beta = (F(1, 2), F(-1, 2))
Xs = {}
for xx in [(F(0), F(0)), one, (F(-1), F(0)), i_, (F(0), F(-1))]:
    Xs[xx] = mul(half_beta, sub(one, xx))
assert Xs[(F(0), F(0))] == half_beta and Xs[one] == (0, 0) and Xs[(F(-1), F(0))] == (F(1), F(-1))
assert Xs[i_] == (F(0), F(-1)) and Xs[(F(0), F(-1))] == one
sq = dict((xx, is_square(v)) for xx, v in Xs.items())
def show(g):
    if g[1] == 0: return str(g[0])
    im = ('' if abs(g[1]) == 1 else str(abs(g[1])) + '*') + 'i'
    return ('-' if g[1] < 0 else '') + im if g[0] == 0 else '%s %s %s' % (g[0], '-' if g[1] < 0 else '+', im)
print('(4) X-coordinates of the affine points: %s' % '; '.join('x = %s -> X = %s (%s)' % (show(xx), show(v), 'a square' if sq[xx] else 'not a square') for xx, v in Xs.items()))
assert [v for xx, v in Xs.items() if sq[xx]] == [(0, 0), one]
assert not is_square(i_) and not is_square((F(1), F(-1)))          # X = 0 needs b = i a square; X = oo needs a = 1 - i a square
print('    X = 0 gives w1 = 0 and w3^2 = i w2^2; X = oo gives w2 = 0 and w3^2 = (1 - i) w1^2: i and 1 - i are not squares of Q(i)')
print('    so every point of C over Q(i) has w1^2 = w2^2: the eight points (1 : +-1 : +-1 : +-1), where w3^2 = w4^2 = w1^2 too')
assert a + b == 1 and c + d == 1
print('PROP1 OK')
