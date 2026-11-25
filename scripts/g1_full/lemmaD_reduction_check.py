"""Exact cross-check of lean_rank/LeanRank/LemmaDReduction.lean (written by an AI model; not read by a person).

E1 : Y^2 = X(X^2 + 8dX - 8d^2), d = 2 + sqrt2, over F = Q(sqrt2); elements of F are pairs (a, b) = a + b sqrt2
of Fractions.  Checks, all exact:
  1. the map (x, w) -> (4(x^2+1)/(x-1+sqrt2)^2, 8w/(x-1+sqrt2)^3) sends the four rational points (+-1, +-4) of
     D2 : w^2 = 2(x^2+1)(x^4-2x^3+2x^2+2x+1) to points of E1, with the coordinates of `pullback`;
  2. the nine affine points of `E1_points_on`, with O, are closed under the chord-tangent law and form a cyclic
     group of order 10 (so `E1Points` says exactly that E1(F) is this group);
  3. the abscissas of the nine points are exactly {0, +-4, +-(12 + 8 sqrt2)} (the list in `E1Points`);
  4. negative tests: a forged point (4, 4 sqrt2) is refused by the curve test; a forged list without 12 + 8 sqrt2
     is refused by check 3; the forged map with 4w in place of 8w is refused by check 1.
Usage: python3 lemmaD_reduction_check.py   (prints ALL CHECKS PASSED)
"""
from fractions import Fraction as Fr

def F(a, b=0): return (Fr(a), Fr(b))
def add(u, v): return (u[0] + v[0], u[1] + v[1])
def sub(u, v): return (u[0] - v[0], u[1] - v[1])
def mul(u, v): return (u[0]*v[0] + 2*u[1]*v[1], u[0]*v[1] + u[1]*v[0])
def inv(u):
    n = u[0]*u[0] - 2*u[1]*u[1]
    assert n != 0
    return (u[0]/n, -u[1]/n)
def div(u, v): return mul(u, inv(v))
def pw(u, k):
    r = F(1)
    for _ in range(k): r = mul(r, u)
    return r

D = F(2, 1)
A2 = mul(F(8), D)            # a2 = 8d
A4 = mul(F(-8), mul(D, D))   # a4 = -8d^2

def on_E1(P):
    if P is None: return True
    X, Y = P
    return mul(Y, Y) == add(add(pw(X, 3), mul(A2, mul(X, X))), mul(A4, X))

def neg(P): return None if P is None else (P[0], (-P[1][0], -P[1][1]))

def addP(P, Q):
    if P is None: return Q
    if Q is None: return P
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2 and add(y1, y2) == F(0): return None
    if P == Q:
        lam = div(add(add(mul(F(3), mul(x1, x1)), mul(mul(F(2), A2), x1)), A4), mul(F(2), y1))
    else:
        lam = div(sub(y2, y1), sub(x2, x1))
    x3 = sub(sub(sub(mul(lam, lam), A2), x1), x2)
    y3 = sub(mul(lam, sub(x1, x3)), y1)
    return (x3, y3)

def image(x, w, k=8):
    z = add(F(x - 1), F(0, 1))
    return (div(mul(F(4), F(x*x + 1)), mul(z, z)), div(mul(F(k), F(w)), pw(z, 3)))

def check(points, xs, k=8):
    # 1. the four rational points of D2 and the coordinates of `pullback`
    for x in (1, -1):
        for w in (4, -4):
            assert 2*(x*x + 1)*(x**4 - 2*x**3 + 2*x*x + 2*x + 1) == w*w
            P = image(x, w, k)
            assert on_E1(P), ('image not on E1', x, w)
            t = Fr(x - 1); Dq = t*t - 2
            a = 4*(x*x + 1)*(t*t + 2)/Dq**2; b = -8*t*(x*x + 1)/Dq**2
            c = 8*t*(t*t + 6)*w/Dq**3; e = -8*(3*t*t + 2)*w/Dq**3
            assert P == ((a, b), (c, e))
    # 2. the nine points and O form a cyclic group of order 10
    for P in points: assert on_E1(P), ('not on E1', P)
    G = [None] + list(points)
    for P in G:
        for Q in G:
            assert addP(P, Q) in G
    gen = (F(-12, -8), F(-32, -24))
    mults, P = [], None
    for _ in range(10):
        mults.append(P); P = addP(P, gen)
    assert P is None and len(set(map(repr, mults))) == 10 and set(map(repr, mults)) == set(map(repr, G))
    # 3. the abscissas
    assert sorted(set(P[0] for P in points)) == sorted(xs)
    return True

POINTS = [(F(0), F(0)), (F(4), F(0, 8)), (F(4), F(0, -8)), (F(-4), F(-16, -8)), (F(-4), F(16, 8)),
          (F(12, 8), F(80, 56)), (F(12, 8), F(-80, -56)), (F(-12, -8), F(-32, -24)), (F(-12, -8), F(32, 24))]
XS = [F(0), F(4), F(-4), F(12, 8), F(-12, -8)]

assert check(POINTS, XS)
# negative tests
def refused(f):
    try: f()
    except AssertionError: return True
    return False
assert refused(lambda: check(POINTS[:1] + [(F(4), F(0, 4))] + POINTS[2:], XS))     # forged point
assert refused(lambda: check(POINTS, [x for x in XS if x != F(12, 8)]))            # forged list
assert refused(lambda: check(POINTS, XS, k=4))                                     # forged map (4w)
print('ALL CHECKS PASSED')
