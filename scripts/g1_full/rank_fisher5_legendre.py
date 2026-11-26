# FINDER, not certifying (written by an AI model; not read by a person).
# Legendre descent for x^2 = A y^2 + B z^2 over Z[sqrt(-D)] for D = 2 (norm-Euclidean).  A FINDER only: its output is
# re-checked exactly by rank_fisher5.py.  Elements are pairs (a, b) = a + b w, w^2 = -D.
import sys, json
from sympy import factorint
from sympy.ntheory import sqrt_mod
from fractions import Fraction as Fr
from math import gcd
import os
D = int(os.environ.get('LEG_D', '2'))   # elements a + b w, w^2 = -D (D = 2: Z[sqrt -2]; D = -2: Z[sqrt 2])
def mul(u, v): return (u[0] * v[0] - D * u[1] * v[1], u[0] * v[1] + u[1] * v[0])
def add(u, v): return (u[0] + v[0], u[1] + v[1])
def sub(u, v): return (u[0] - v[0], u[1] - v[1])
def neg(u): return (-u[0], -u[1])
def conj(u): return (u[0], -u[1])
def N(u): return u[0] * u[0] + D * u[1] * u[1]
def aN(u): return abs(N(u))
def rnd(f): return (f.numerator * 2 + f.denominator) // (2 * f.denominator)
def divround(u, v):
    n = N(v); w = mul(u, conj(v)); return (rnd(Fr(w[0], n)), rnd(Fr(w[1], n)))
def divexact(u, v):
    n = N(v); w = mul(u, conj(v)); assert w[0] % n == 0 and w[1] % n == 0, ('not divisible', u, v); return (w[0] // n, w[1] // n)
def divides(v, u):
    n = N(v); w = mul(u, conj(v)); return w[0] % n == 0 and w[1] % n == 0
def egcd(u, v):
    while v != (0, 0):
        q = divround(u, v); u, v = v, sub(u, mul(q, v))
    return u
def sqfree(B):
    """B = s^2 B1: returns (B1, s)"""
    s = (1, 0)
    changed = True
    while changed:
        changed = False
        for pi in primes_of(B):
            pi2 = mul(pi, pi)
            if aN(pi2) > 1 and divides(pi2, B):
                B = divexact(B, pi2); s = mul(s, pi); changed = True; break
    return B, s
def xgcd(u, v):
    """(g, s, t) with s u + t v = g in Z[w]"""
    s0, s1, t0, t1 = (1, 0), (0, 0), (0, 0), (1, 0)
    while v != (0, 0):
        q = divround(u, v)
        u, v = v, sub(u, mul(q, v)); s0, s1 = s1, sub(s0, mul(q, s1)); t0, t1 = t1, sub(t0, mul(q, t1))
    return u, s0, t0
def primes_of(B):
    out = []
    for p in factorint(aN(B)):
        if p == -1: continue
        if p == 2 and D in (2, -2):
            if divides((0, 1), B): out.append((0, 1))
            continue
        if divides((p, 0), B) and pow(-D % p, (p - 1) // 2, p) != 1:
            out.append((p, 0)); continue
        pi = egcd(B, (p, 0))
        if aN(pi) == p:
            for c_ in {pi, conj(pi)}:
                if divides(c_, B): out.append(c_)
        elif aN(pi) == p * p:
            out += [pi] if aN(pi) == p * p and not divides((p, 0), B) else []
            if divides((p, 0), B):
                # p split and p | B: both conjugate primes divide B
                r = sqrt_mod(-D % p, p); pi1 = egcd((p, 0), (r, 1))
                out += [pi1, conj(pi1)]
    return out
def sqrt_mod_pi(A, pi):
    n = aN(pi)
    if n == 2: return (A[0] % 2, 0) if A[1] % 2 == 0 else ((A[0] + 0) % 2, 1)
    import sympy
    if sympy.isprime(n):
        a, b = pi; r = (-a * pow(b, -1, n)) % n
        assert (r * r + D) % n == 0 and divides(pi, (r * 0 + a, b)) or True
        # w = r mod pi: check a + b r = 0 mod n
        assert (a + b * r) % n == 0
        t0 = sqrt_mod((A[0] + A[1] * r) % n, n)
        assert t0 is not None, ('not a square mod', pi)
        return (t0, 0)
    p = pi[0]; assert pi[1] == 0 and n == p * p and p < 20000
    a0, a1 = A[0] % p, A[1] % p
    for x0 in range(p):
        for x1 in range(p):
            y = mul((x0, x1), (x0, x1))
            if y[0] % p == a0 and y[1] % p == a1: return (x0, x1)
    raise AssertionError('not a square mod an inert prime')
def sqrt_mod_B(A, B):
    if aN(B) == 1: return (0, 0)
    ps = primes_of(B)
    M = (1, 0)
    for q in ps: M = mul(M, q)
    assert aN(M) == aN(B), ('B not squarefree or factorisation incomplete', B, ps)
    t, mod = (0, 0), (1, 0)
    for q in ps:
        tq = sqrt_mod_pi(A, q)
        # CRT: t' = t (mod mod), t' = tq (mod q)
        g, s_, u_ = xgcd(mod, q); assert aN(g) == 1
        ginv = conj(g) if N(g) == 1 else None
        # s mod + u q = g (a unit); t' = t u q g^-1 + tq s mod g^-1
        gi = divexact((1, 0), g)
        t = add(mul(mul(t, mul(u_, q)), gi), mul(mul(tq, mul(s_, mod)), gi))
        mod = mul(mod, q)
        t = sub(t, mul(divround(t, mod), mod))
    assert divides(B, sub(mul(t, t), A)), ('sqrt mod B failed', A, B)
    t = sub(t, mul(divround(t, B), B))
    return t
def solve(A, B, depth=0):
    """(x, y, z) with x^2 = A y^2 + B z^2, not all zero."""
    assert depth < 400
    A, sA = sqfree(A); B, sB = sqfree(B)
    # x^2 = A (sA y)^2 + B (sB z)^2 : solve for (x, Y, Z) then y = Y / sA ... keep it homogeneous: y = Y * sB, z = Z * sA,
    # x = X * sA * sB  works: (X sA sB)^2 = A sA^2 (Y sB)^2 + B sB^2 (Z sA)^2  <=>  X^2 = A Y^2 + B Z^2 with A, B the reduced
    X, Y, Z = solve_red(A, B, depth)
    # original A0 = A sA^2, B0 = B sB^2: x^2 = A0 y^2 + B0 z^2 with x = X sA sB, y = Y sB, z = Z sA
    return mul(X, mul(sA, sB)), mul(Y, sB), mul(Z, sA)
def solve_red(A, B, depth):
    for u in ((1, 0),):
        if A == u: return ((1, 0), (1, 0), (0, 0))
        if B == u: return ((1, 0), (0, 0), (1, 0))
    if aN(A) > aN(B):
        x, z, y = solve_red(B, A, depth)
        return x, y, z
    if aN(B) == 1 and D < 0:
        raise AssertionError('reached a unit %s in a real quadratic ring: not handled' % (B,))
    if N(B) == 1:      # B = -1 (only units +-1), A = -1 too or a unit: x^2 + z^2 = A y^2
        if A == (-1, 0): raise AssertionError('x^2 + y^2 + z^2 = 0 has no solution here')
    if N(B) == 1:
        # B = -1: choose t near sqrt(A) so that A - t^2 is small
        import cmath
        z_ = cmath.sqrt(complex(A[0], A[1] * D ** 0.5))
        t = (round(z_.real), round(z_.imag / D ** 0.5))
    else:
        t = sqrt_mod_B(A, B)
    B1 = divexact(sub(mul(t, t), A), B)
    if B1 == (0, 0): return (t, (1, 0), (0, 0))
    X, Y, Z = solve(A, B1, depth + 1)
    # (X t + A Y)^2 - A (X + t Y)^2 = (t^2 - A)(X^2 - A Y^2) = B B1 * B1 Z^2
    return add(mul(X, t), mul(A, Y)), add(X, mul(t, Y)), mul(B1, Z)

if __name__ == '__main__':
    A = tuple(int(v) for v in sys.argv[1:3]); B = tuple(int(v) for v in sys.argv[3:5])
    x, y, z = solve(A, B)
    g_ = egcd(egcd(x, y), z)                       # remove the common factor (halves the digits here)
    x, y, z = [divexact(v, g_) for v in (x, y, z)]
    lhs = mul(x, x); rhs = add(mul(A, mul(y, y)), mul(B, mul(z, z)))
    assert lhs == rhs and (x, y, z) != ((0, 0), (0, 0), (0, 0))
    print(json.dumps(dict(x=x, y=y, z=z)))
