#!/usr/bin/env python3
"""Exact cross-check of the descent behind lean_rank/LeanRank/E1Descent.lean (rank E1(Q(sqrt2)) = 0).

Written by an AI model (Claude, made by Anthropic) on 26 November 2026; not read by a person.
Untrusted helper: the Lean files are the certificate.  Everything here is exact integer arithmetic in
Z[sqrt2] (pairs (a, b) = a + b*sqrt2) and sympy polynomial identities.

Checks
 1. the constants: A1 = sqrt2^7 eps = 16 + 8 sqrt2 = 8d, B1 = -sqrt2^8 eps^2 = -8d^2, A2 = -sqrt2^9 eps = -16d,
    B2 = 3 sqrt2^12 eps^2 = 96 d^2 = A1^2 - 4 B1, 4 B2 - A2^2 = 768 + 512 sqrt2 (d = 2 + sqrt2, eps = 1 + sqrt2);
 2. the polynomial identities used in V1, V2, M1, prime_dvd_B, real_excl, pow4 (sympy, pi2 a free symbol);
 3. the classes of the x-coordinates of the ten torsion points: E1 gives only 1, -1 (= classes of x and b/x),
    E1' (their images) only 1, 3 -- consistent with the descent result;
 4. a search: for every class delta the Lean proof excludes (E1: +-sqrt2, +-eps sqrt2, +-eps;
    E1': the 14 classes other than 1, 3), no (s, t) with coordinates in [-B, B] makes F_delta(s, t) a nonzero
    square in Z[sqrt2]; for the kept classes (E1: 1, -1; E1': 1, 3) such (s, t) exist;
 5. negative tests: a forged exclusion of a kept class, a forged constant and a forged identity are refused.
"""
from itertools import product
from math import isqrt
import sympy as sp

def mul(x, y): return (x[0]*y[0] + 2*x[1]*y[1], x[0]*y[1] + x[1]*y[0])
def add(x, y): return (x[0] + y[0], x[1] + y[1])
def neg(x): return (-x[0], -x[1])
def pw(x, k):
    r = (1, 0)
    for _ in range(k): r = mul(r, x)
    return r
def norm(x): return x[0]*x[0] - 2*x[1]*x[1]
PI, EPS, THREE, D = (0, 1), (1, 1), (3, 0), (2, 1)

def is_square(z):
    """exact: is z a square in Z[sqrt2]?  (z = (a + b w)^2 = a^2 + 2b^2 + 2ab w)"""
    a, b = z
    if (a, b) == (0, 0): return True
    N = norm(z)
    if N < 0: return False
    n = isqrt(N)
    if n*n != N: return False
    for s in (n, -n):          # x^2 = (a + s)/2 with x = re of the root
        if (a + s) % 2: continue
        x2 = (a + s)//2
        if x2 < 0: continue
        x = isqrt(x2)
        if x*x != x2: continue
        if x == 0:
            if a % 2 == 0 and b == 0 and a//2 >= 0:
                y = isqrt(a//2)
                if 2*y*y == a: return True
            continue
        if b % (2*x): continue
        y = b // (2*x)
        if mul((x, y), (x, y)) == z: return True
    return False

def F(A, B, c, s, t):
    s2, t2 = mul(s, s), mul(t, t)
    inner = add(add(mul(mul(c, c), mul(s2, s2)), mul(mul(A, c), mul(s2, t2))), mul(B, mul(t2, t2)))
    return mul(c, inner)

A1, B1 = mul(pw(PI, 7), EPS), neg(mul(pw(PI, 8), pw(EPS, 2)))
A2, B2 = neg(mul(pw(PI, 9), EPS)), mul(mul(THREE, pw(PI, 12)), pw(EPS, 2))

def check_constants(A1, B1, A2, B2):
    ok = A1 == (16, 8) == mul((8, 0), D)
    ok &= B1 == (-48, -32) == mul((-8, 0), mul(D, D))
    ok &= A2 == (-32, -16) == mul((-16, 0), D)
    ok &= B2 == (576, 384) == add(mul(A1, A1), mul((-4, 0), B1)) == mul((96, 0), mul(D, D))
    ok &= add(mul((4, 0), B2), neg(mul(A2, A2))) == (768, 512)
    return ok

def check_identities():
    p, e, w, s, t, r, u, k1, k2, c, A, B, T3 = sp.symbols('p e w s t r u k1 k2 c A B T3')
    s1, s2_, s3, t1, cc, m, U = sp.symbols('s1 s2 s3 t1 cc m U')
    Fq = lambda A, B, c, s, t: c*(c**2*s**4 + A*c*s**2*t**2 + B*t**4)
    a1, b1, a2, b2 = p**7*e, -(p**8*e**2), -(p**9*e), T3*p**12*e**2
    J = lambda s, t: u*(u**2*s**4 + p**3*e*u*s**2*t**2 - e**2*t**4)
    ids = [
     (r**2*Fq(a1,b1,p*w,p*p*s2_,t), (r*p**4)**2*(p*(-(e**2*w*t**4)+p*(p*(w**3*s2_**4+p**2*e*w**2*s2_**2*t**2))))),
     (r**2*Fq(a1,b1,p*w,p*s1,t), (r*p**3)**2*(p*(w**3*s1**4+p*(p*(p**2*e*w**2*s1**2*t**2-e**2*w*t**4))))),
     (r**2*Fq(a1,b1,p*w,s,t), (r*p)**2*(p*(w**3*s**4+p*(p**5*(e*w**2*s**2*t**2-e**2*w*t**4))))),
     (r**2*Fq(a2,b2,p*w,p**3*s3,t), (r*p**6)**2*(p*(T3*e**2*w*t**4+p*(p*(w**3*s3**4-p**2*e*w**2*s3**2*t**2))))),
     (r**2*Fq(a2,b2,p*w,p**2*s2_,t), (r*p**5)**2*(p*(w**3*s2_**4+p*(p*(-(p**2*e*w**2*s2_**2*t**2)+T3*e**2*w*t**4))))),
     (r**2*Fq(a2,b2,p*w,p*s1,t), (r*p**3)**2*(p*(w**3*s1**4+p*(p**5*(-(e*w**2*s1**2*t**2)+T3*e**2*w*t**4))))),
     (r**2*Fq(a2,b2,p*w,s,t), (r*p)**2*(p*(w**3*s**4+p*(p**7*(-(e*w**2*s**2*t**2)+p**2*T3*e**2*w*t**4))))),
     (Fq(a1,b1,p**4*u,s,t), (p**6)**2*J(s,t)),
     (J(p*s1,t) + (u*e**2)*(t**4-(1+p**5*k2)), -(u*e**2)+p**2*(u*(p**2*u**2*s1**4+p**3*e*u*s1**2*t**2-p**3*e**2*k2))),
     (J(s,p*t1) - u**3*(s**4-(1+p**5*k1)) - u*(u**2-e**2), u*e**2+p**2*(p**3*u**3*k1+p**3*e*u**2*s**2*t1**2-p**2*u*e**2*t1**4)),
     (J(s,t) - (u*s**4*(u**2-e**2)+u*e**2*(s**4-(1+p**5*k1))-u*e**2*(t**4-(1+p**5*k2))), p**3*(u*(e*u*s**2*t**2+p*(p*e**2*(k1-k2))))),
     (r**2*Fq(A,B,p*cc,s,p*t1), (r*p)**2*(p*(cc**3*s**4+p*(A*cc**2*s**2*t1**2+p*B*cc*t1**4)))),
     (r**2*Fq(A,B,p*cc,s,t), r**2*(p*(cc*(p*(p*cc**2*s**4+A*cc*s**2*t**2)+B*t**4)))),
     (4*(c**2*s**4+A*c*s**2*t**2+B*t**4), (2*c*s**2+A*t**2)**2+(4*B-A**2)*t**4),
    ]
    # pow4 modulo p^2 = 2: (1 + p c)^2 = 1 + p^2 (p c + c^2) + (-(p c)) (p^2 - 2)
    ids.append(((1+p*c)**2, 1+p**2*(p*c+c**2) - (p*c)*(p**2-2)))
    return all(sp.expand(l - r_) == 0 for l, r_ in ids), len(ids)

def sclass(z):
    """class of z = (a, b) in K*/K*^2 among products of -1, eps, sqrt2, 3, as a tuple; z an S-unit times square"""
    from fractions import Fraction as Fr
    # brute: find c in the 16 representatives with c*z a square in Q(sqrt2): test c*z*(denominator)^2 integral square
    reps = [(sg, e, i, j) for sg in (0, 1) for e in (0, 1) for i in (0, 1) for j in (0, 1)]
    for rep in reps:
        c = (1, 0)
        if rep[0]: c = neg(c)
        if rep[1]: c = mul(c, EPS)
        if rep[2]: c = mul(c, PI)
        if rep[3]: c = mul(c, THREE)
        if is_square(mul(c, z)): return rep
    return None

def search(A, B, c, bound):
    rng = range(-bound, bound + 1)
    elems = [(a, b) for a in rng for b in rng if (a, b) != (0, 0)]
    for s in elems:
        for t in elems:
            f = F(A, B, c, s, t)
            if f != (0, 0) and is_square(f):
                return (s, t)
    return None

def rep(sg, e, i, j):
    c = (1, 0)
    if sg: c = neg(c)
    if e: c = mul(c, EPS)
    if i: c = mul(c, PI)
    if j: c = mul(c, THREE)
    return c

def main():
    out = []
    say = lambda s: (print(s), out.append(s))
    assert check_constants(A1, B1, A2, B2); say('1 constants A1 = 8d, B1 = -8d^2, A2 = -16d, B2 = 96d^2 = A1^2 - 4B1, 4B2 - A2^2 = 768 + 512 sqrt2: OK')
    ok, k = check_identities(); assert ok; say('2 %d polynomial identities of E1Descent.lean (V1, V2, M1, prime_dvd_B, real_excl, pow4): OK' % k)
    # torsion points of E1 (x-coordinates 4, -4, 12 + 8 sqrt2, -12 - 8 sqrt2; 0 is T): classes of x, and of b/x
    E1x = [(4, 0), (-4, 0), (12, 8), (-12, -8)]
    cls = set()
    for x in E1x:
        cls.add(sclass(x))
    assert cls <= {(0, 0, 0, 0), (1, 0, 0, 0)}, cls
    say('3a classes of x on the eight affine torsion points of E1 with x != 0: %s (only 1, -1)' % sorted(cls))
    # images on E1': X = y^2/x^2 = (x^2 + A x + B)/x  (as an element of Q(sqrt2)); class of X*x^2 = x*(x^2+Ax+B)
    cls2 = set()
    for x in E1x:
        val = mul(x, add(add(mul(x, x), mul(A1, x)), B1))   # = y^2, a square; X = y^2/x^2 has class of y^2 -> 1
        assert is_square(val)
        cls2.add(sclass(val))
    # X(T') = 0 is not a point with X != 0; the class 3 is the image of T': alpha'(T') = B2 = 3 * square
    cls2.add(sclass(B2))
    assert cls2 == {(0, 0, 0, 0), (0, 0, 0, 1)}, cls2
    say("3b classes on E1': images of the torsion points -> 1, alpha'(T') = B2 -> 3: %s" % sorted(cls2))
    bound = 4
    kept1 = {(0, 0, 0, 0), (1, 0, 0, 0)}
    kept2 = {(0, 0, 0, 0), (0, 0, 0, 1)}
    for name, A, B, classes, kept in (('E1', A1, B1, [(sg, e, i, 0) for sg in (0, 1) for e in (0, 1) for i in (0, 1)], kept1),
                                      ("E1'", A2, B2, [(sg, e, i, j) for sg in (0, 1) for e in (0, 1) for i in (0, 1) for j in (0, 1)], kept2)):
        for cl in classes:
            c = rep(*cl)
            # for M1 the Lean proof uses delta = sqrt2^4 * u; the class is the same
            sol = search(A, B, c, bound)
            if cl in kept:
                assert sol is not None, (name, cl)
                say('4 %-3s delta = %-9s kept: F(s,t) a square at (s,t) = %s' % (name, c, sol))
            else:
                assert sol is None, (name, cl, sol)
                say('4 %-3s delta = %-9s excluded in Lean; no square F(s,t), coordinates in [-%d, %d]' % (name, c, bound, bound))
    # 5 negative tests
    neg_ok = 0
    if not check_constants((16, 7), B1, A2, B2): neg_ok += 1           # forged constant
    p, s = sp.symbols('p s')
    if sp.expand((1 + p*s)**2 - (1 + p**2*(p*s + s**2))) != 0: neg_ok += 1   # forged identity (needs p^2 = 2)
    if search(A1, B1, (-1, 0), bound) is not None: neg_ok += 1          # 'delta = -1 excluded on E1' refuted
    if search(A2, B2, THREE, bound) is not None: neg_ok += 1            # 'delta = 3 excluded on E1'' refuted
    assert neg_ok == 4
    say('5 negative tests: 4 forged claims refused')
    say('e1_rank0_check OK')
    open(__file__.replace('.py', '.txt'), 'w').write('\n'.join(out) + '\n')

if __name__ == '__main__':
    main()
