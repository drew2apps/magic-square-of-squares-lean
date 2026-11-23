#!/usr/bin/env python3
"""The model and the three identities of dm4.py, with m and n as indeterminates (exact computer algebra).

C: y^2 = F0 GL = (n^2 t^2 + m^2)(2S t^2 - 8mn t + 2S), S = m^2 + n^2, with the point t = m/n, as an elliptic curve by
Connell's formulas (the curve E0 of dm3.py and dm4.py).  Printed:
  - the three points of order 2, and the model  Y^2 = x (x - 2A)(x - 2B)  with x = X + 6m^4 + 2n^4,
        A = m^4 - 2m^3 n + 2m^2 n^2 + 2m n^3 + n^4,   B = m^4 + 2m^3 n + 2m^2 n^2 - 2m n^3 + n^4,   A B = m^8 + 14 m^4 n^4 + n^8;
  - (I0)    x(2Q - R0)      = 4 Pi(t)^2 / (F0 GL),  Pi = 2mn (n^2 t^2 - m^2) + (m^4 - n^4) t;
  - (I+-)   x(Q1 + Q2 - R0) = ((y1 - y2)/(2t))^2,       x(Q1 - Q2) = ((y1 + y2)/(2t))^2,       Q2 = (-t, y2), y2^2 = F0 GU;
  - (Iinv)  x(Q1 + Q3 - R0) = ((y1 - y3)/(t^2 - 1))^2,  x(Q1 - Q3) = ((y1 + y3)/(t^2 - 1))^2,  Q3 = (1/t, y3/t^2), y3^2 = F8 GL.
The identities are checked in the fields Q(m, n)(t)(y1, y2) and Q(m, n)(t)(y1, y3), so they hold for every pair (m, n)
and every point where the expressions are defined.  dm4.py checks them again for the pair it is run on.
usage: dm4_identities.py     (needs passagemath-schemes with Singular; a few seconds)"""
import sys, time
from sage.all__sagemath_schemes import *
T0 = time.time()
A = PolynomialRing(QQ, 'm,n'); K0 = FractionField(A); m, n = K0.gens()
Rt = PolynomialRing(K0, 't'); Kt = FractionField(Rt); t = Kt.gen(); Tt = Rt.gen()
S = m*m + n*n
def qpoly(x): return (n*n*x**2 + m*m) * (2*S*x**2 - 8*m*n*x + 2*S)
F0 = n*n*t**2 + m*m; F8 = m*m*t**2 + n*n; GL = 2*S*t**2 - 8*m*n*t + 2*S; GU = 2*S*t**2 + 8*m*n*t + 2*S
t1 = m/n; q0 = 2*m*(m*m - n*n)/n
qs = qpoly(t1 + Tt); e0, d, c, b, a = [qs[i] for i in range(5)]
assert e0 == q0**2
a1 = d/q0; a2 = c - d**2/(4*q0**2); a3 = 2*q0*b; a4 = -4*q0**2*a; a6 = a2*a4
print('a1 =', a1, '\na2 =', a2, '\na3 =', a3, '\na4 =', a4)
RX = PolynomialRing(K0, 'X'); X = RX.gen()
two = 4*(X**3 + a2*X**2 + a4*X + a6) + (a1*X + a3)**2
print('2-division polynomial factors:', two.factor())
e1 = -(6*m**4 + 2*n**4)
Aq = m**4 - 2*m**3*n + 2*m**2*n**2 + 2*m*n**3 + n**4; Bq = m**4 + 2*m**3*n + 2*m**2*n**2 - 2*m*n**3 + n**4
assert two == 4*(X - e1)*(X - e1 - 2*Aq)*(X - e1 - 2*Bq) and Aq*Bq == m**8 + 14*m**4*n**4 + n**8
print('model: Y^2 = x (x - 2A)(x - 2B), A B = m^8 + 14 m^4 n^4 + n^8   [%.0f s]' % (time.time() - T0), flush=True)
class F4:
    r1 = None; r2 = None
    def __init__(s, a=0, b=0, c=0, d=0): s.v = (Kt(a), Kt(b), Kt(c), Kt(d))
    @staticmethod
    def co(o): return o if isinstance(o, F4) else F4(o)
    def __add__(s, o): o = F4.co(o); return F4(*[x + y for x, y in zip(s.v, o.v)])
    __radd__ = __add__
    def __neg__(s): return F4(*[-x for x in s.v])
    def __sub__(s, o): return s + (-F4.co(o))
    def __rsub__(s, o): return (-s) + o
    def __mul__(s, o):
        o = F4.co(o); a, b, c, d = s.v; e, f, g, h = o.v; r1, r2 = F4.r1, F4.r2
        return F4(a*e + b*f*r1 + c*g*r2 + d*h*r1*r2, a*f + b*e + (c*h + d*g)*r2, a*g + c*e + (b*h + d*f)*r1, a*h + d*e + b*g + c*f)
    __rmul__ = __mul__
    def conj1(s): a, b, c, d = s.v; return F4(a, -b, c, -d)
    def inv(s):
        p = s * s.conj1(); assert p.v[1] == 0 and p.v[3] == 0
        nn = p.v[0]**2 - p.v[2]**2 * F4.r2
        return s.conj1() * F4(p.v[0]/nn, 0, -p.v[2]/nn, 0)
    def __truediv__(s, o): return s * F4.co(o).inv()
    def __eq__(s, o): o = F4.co(o); return all(x == y for x, y in zip(s.v, o.v))
def toE0(tt, y):
    tau = tt - t1
    return (2*q0*(y + q0) + d*tau) * (1/tau**2), (4*q0**2*(y + q0) + 2*q0*(d*tau + c*tau**2) - d**2*tau**2/(2*q0)) * (1/tau**3)
def neg(P): return P[0], -P[1] - a1*P[0] - a3
def addX(P, R):
    lam = (R[1] - P[1]) / (R[0] - P[0]); return lam*lam + a1*lam - a2 - P[0] - R[0]
Pi = 2*m*n*(n*n*t**2 - m*m) + (m**4 - n**4)*t
F4.r1 = qpoly(t); F4.r2 = Kt(1); y1 = F4(0, 1, 0, 0)
X0 = addX(toE0(t, y1), neg(toE0(t, -y1)))
print('(I0)  ', X0 - e1 == 4*Pi**2/(F0*GL), '[%.0f s]' % (time.time() - T0), flush=True)
F4.r2 = qpoly(-t); y2 = F4(0, 0, 1, 0)
Q1 = toE0(t, y1); Q2 = toE0(-t, y2); iQ2 = toE0(-t, -y2)
print('(I+-) ', addX(Q1, neg(iQ2)) - e1 == ((y1 - y2)/(2*t))*((y1 - y2)/(2*t)), addX(Q1, neg(Q2)) - e1 == ((y1 + y2)/(2*t))*((y1 + y2)/(2*t)), '[%.0f s]' % (time.time() - T0), flush=True)
F4.r2 = F8*GL; y3 = F4(0, 0, 1, 0)
Q3 = toE0(1/t, y3*(1/t**2)); iQ3 = toE0(1/t, -y3*(1/t**2))
print('(Iinv)', addX(Q1, neg(iQ3)) - e1 == ((y1 - y3)/(t*t - 1))*((y1 - y3)/(t*t - 1)), addX(Q1, neg(Q3)) - e1 == ((y1 + y3)/(t*t - 1))*((y1 + y3)/(t*t - 1)), '[%.0f s]' % (time.time() - T0), flush=True)
