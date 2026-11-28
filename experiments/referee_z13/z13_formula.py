"""Reviewer z13's own implementation of formula (2), the set S' and the form Theta (3) of PAIRING_FORMULA_2.md,
written from the note's text only (the note's closed_formula.py / isogeny_part.py are not imported).
Written by an AI model (reviewer z13), not read by a person.

E: y^2 = x(x - A^2)(x - C^2), (A, B, C) a primitive Pythagorean triple, C odd, A > 0.
Classes are (c_r, c_m) at (x - e_r, x - e_m) = (x, x - A^2); c_* = c_r c_m mod squares.
"""
from fractions import Fraction as Fr
from math import isqrt, gcd
import itertools
import cypari2
pari = cypari2.Pari()
pari.allocatemem(4 * 10 ** 8)

_FC = {}
def factor(n):
    """prime -> exponent of |n| (PARI factor; primality certified by isprime)."""
    n = abs(int(n))
    if n in _FC: return _FC[n]
    if n <= 1: _FC[n] = {}; return {}
    f = pari.factor(n)
    out = {int(q): int(e) for q, e in zip(f[0], f[1])}
    for q in out: assert pari.isprime(q)
    _FC[n] = out
    return out

def sqf(x):
    """signed squarefree part of a non-zero rational."""
    x = Fr(x); assert x != 0
    a = abs(x.numerator * x.denominator); r = 1
    for p, e in factor(a).items():
        if e % 2: r *= p
    return r if x > 0 else -r

def vp(n, p):
    n = abs(n); v = 0
    while n % p == 0: n //= p; v += 1
    return v

def is_sq_Qp(x, p):
    """x a non-zero rational: is x a square in Q_p?"""
    x = Fr(x); n, d = x.numerator, x.denominator
    v = vp(n, p) - vp(d, p)
    if v % 2: return False
    u = (n // p ** vp(n, p)) * (d // p ** vp(d, p))
    if p == 2: return u % 8 == 1
    return pow(u % p, (p - 1) // 2, p) == 1

def hilb(x, y, p):
    return int(pari.hilbert(x, y, p))

def padic_sqrt(x, p):
    """a square root in Q_p of the rational x (assumed a square there), as a PARI p-adic."""
    x = Fr(x)
    N = 90 if p == 2 else 40
    return pari('sqrt((%d)/(%d) + O(%d^%d))' % (x.numerator, x.denominator, p, N))

def triple(m, n, cv):
    """(A, B, C) of E_08 / E_LU for coprime m > n > 0, divided by 2 when m, n are both odd."""
    if cv == '08': A, B, C = m * m - n * n, 2 * m * n, m * m + n * n
    else: A, B, C = 2 * m * n, m * m - n * n, m * m + n * n
    if A % 2 == 0 and C % 2 == 0: A //= 2; B //= 2; C //= 2
    assert A > 0 and B > 0 and C % 2 == 1 and A * A + B * B == C * C and gcd(A, C) == 1
    return A, B, C

def H2(A, C):
    """the note's hypothesis H2: if A is odd, 32 | C - A or 32 | C + A."""
    return A % 2 == 0 or (C - A) % 32 == 0 or (C + A) % 32 == 0

def formula2(A, B, C, dl, b_r, b_m):
    """<(1, dl, dl), (b_r, b_m, b_*)> by formula (2) of the note; None where it is undefined
    (some v | AC (or 2 | A) with neither dl nor -dl a square)."""
    dl = sqf(dl)
    b_s = sqf(b_r * b_m)
    val = 1
    # first product: p | B with C = -A (p); for p = 2 (A odd) 'C = -A mod 8'
    for p in factor(B):
        if p == 2:
            if (C + A) % 8 == 0: val *= hilb(dl, b_s, 2)
        elif (C + A) % p == 0:
            val *= hilb(dl, b_s, p)
    # second product: v | AC, and v = 2 if A is even
    places = set(factor(A * C))
    if A % 2 == 0: places.add(2)
    for p in sorted(places):
        if is_sq_Qp(dl, p):
            val *= hilb(2 * padic_sqrt(dl, p), b_r, p)
        elif is_sq_Qp(-dl, p):
            val *= hilb(-A, b_s, p) * hilb(-C, b_m, p) * hilb((A - C) * padic_sqrt(-dl, p), b_r, p)
        else:
            return None
    return val

def Sprime(A, B, C):
    """S' of the note, from factorisations: squarefree beta built from -1 and the primes of B, +- a square in Q_p at
    every p | AC, and in Q_2 if A is even."""
    gens = [-1] + sorted(factor(B))
    places = set(factor(A * C))
    if A % 2 == 0: places.add(2)
    out = set()
    for mask in itertools.product((0, 1), repeat=len(gens)):
        b = 1
        for g, e in zip(gens, mask):
            if e: b *= g
        if all(is_sq_Qp(b, p) or is_sq_Qp(-b, p) for p in places): out.add(b)
    return out

def Theta(A, B, C, dl, be):
    """the form (3) of the note on S' x S'."""
    val = 1
    for p in factor(B):
        if p == 2:
            if (C + A) % 8 == 0: val *= hilb(dl, be, 2)
        elif (C + A) % p == 0:
            val *= hilb(dl, be, p)
    places = set(factor(A * C))
    if A % 2 == 0: places.add(2)
    for p in places:
        if (not is_sq_Qp(dl, p)) and is_sq_Qp(-dl, p):
            val *= hilb(A * C, be, p)
    return val
