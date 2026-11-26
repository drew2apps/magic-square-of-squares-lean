"""x3 referee: own 2-adic solubility test for the first-descent quartic of the open E'-coset of D1 of (12007,1799).
Q(M, e) = y M^4 + A M^2 e^2 + (B/y) e^4 (A = a', B = b' of the isogenous curve), (M, e) primitive in Z_2^2.
Breadth-first over residues mod 2^k of (M, e): a class is dropped when Q is determined mod 2^(v+3) with v = v_2(Q)
and the value is not v even with unit = 1 mod 8; a class survives when it is determined to be a square.
Written by an AI model; not read by a person."""
import sys
from fractions import Fraction as Fr
import cypari2
pari = cypari2.Pari()
def v2(x):
    k = 0
    while x % 2 == 0: x //= 2; k += 1
    return k
def taylor(f, t):
    """coefficients of f(t + z) in z (f low degree first)"""
    n = len(f); out = []
    g = list(f)
    for i in range(n):
        # value of g at t is coefficient i; then g <- (g - g(t)) / (z) via synthetic division
        acc = 0; q = [0] * (len(g) - 1)
        for j in range(len(g) - 1, -1, -1):
            acc = acc * t + g[j]
            if j: q[j - 1] = acc
        out.append(acc); g = q
        if not g: break
    return out
def soluble(y, A, B, kmax=120, cap=200000):
    """Q(M, e) = y M^4 + A M^2 e^2 + (B/y) e^4 a square (0 allowed) at some point of P^1(Q_2)?  Charts (M, 1), M in Z_2,
    and (1, e), e in 2 Z_2.  Branch on t mod 2^k; Q(t + 2^k z) = Q(t) mod 2^delta, delta = min_i v(Taylor_i(t)) + i k."""
    assert B % y == 0; c = B // y
    charts = [([c, 0, A, 0, y], [(0, 0)]),            # f(M) at e = 1, start: M mod 1
              ([y, 0, A, 0, c], [(0, 1)])]            # f(e) at M = 1, e = 0 mod 2
    undecided = False
    for f, start in charts:
        live = [(t, k) for t, k in start]
        while live:
            t, k = live.pop()
            T = taylor(f, t); val = T[0]
            delta = min([v2(T[i]) + i * k for i in range(1, len(T)) if T[i] != 0] + [10 ** 9])
            if val == 0: return True, 'root'
            w = v2(val)
            if w + 3 <= delta:
                if w % 2 == 0 and (val >> w) % 8 == 1: return True, (t, k)
                continue                                   # determined non-square on the whole disc
            if k >= kmax: undecided = True; continue
            live += [(t, k + 1), (t + 2 ** k, k + 1)]
            if len(live) > cap: return None, 'cap'
    return (None, 'depth') if undecided else (False, 'all discs closed')

def d1(m, n):
    S = m * m + n * n; D2 = (m * m - n * n) ** 2
    c3 = 4 * S * S * m * m * n * n; c2 = -16 * S * m ** 3 * n ** 3; c1 = 4 * S * S * D2; c0 = -16 * S * m * n * D2
    a2, a4, a6 = c1, c0 * c2, c0 * c0 * c3
    rts = [int(r) for r in pari.nfroots(None, pari('x^3+(%d)*x^2+(%d)*x+(%d)' % (a2, a4, a6)))]
    assert len(rts) == 1; e = rts[0]; assert e ** 3 + a2 * e * e + a4 * e + a6 == 0
    return 3 * e + a2, 3 * e * e + 2 * a2 * e + a4
if __name__ == '__main__':
    a, b = d1(12007, 1799); A, B = -2 * a, a * a - 4 * b
    print('D1 of (12007,1799): a =', a, ' b =', b)
    for y in (-6177769598, 6177769598, 73917229246, -73917229246):
        print('open class', y, ':', soluble(y, A, B), flush=True)
    for y in (-1, 244673, -244673):          # point classes: must not come out False
        print('point class', y, ':', soluble(y, A, B), flush=True)
