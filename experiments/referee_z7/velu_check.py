#!/usr/bin/env python3
"""Review z7 of roadmap H, H7 claim 1: the shift and Velu's 2-isogeny for F4's g4*K curve, as polynomial identities.
Written by an AI model (review agent z7, 27 November 2026); not read by a person.  Experiment, not counted.

(1) Velu: for E: y^2 = x(x^2 + ax + b), the map (x, y) -> (y^2/x^2, y(b - x^2)/x^2) sends E to
    E': Y^2 = X(X^2 - 2aX + a^2 - 4b): the residual Y^2 - X(X^2 - 2aX + a^2 - 4b), with y^2 = x(x^2 + ax + b), is 0.
(2) C = dm_g4k.g4k_curve's formula (x = c3 s, y^2 = x^3 + c2 x^2 + c1 c3 x + c0 c3^2, P(s) = g4(s) K(s) with
    g4 = 2S(1 + s), K = S^2 + (2S^2 - 16m^2n^2)s + S^2 s^2): x = -c3 is a root; the shift X = x + c3 gives
    a = -32 m^2 n^2 S, b = 64 m^2 n^2 S^4; 16b is a square for all m, n; E' = X(X - 16Smn(m-n)^2)(X + 16Smn(m+n)^2).
    (run_ct_g4k.py checks these against dm_g4k.g4k_curve itself at each ratio.)
sympy only, symbolic in m, n; no PARI.

    python3 velu_check.py            prints velu_check.txt   (a few seconds)
"""
import sympy as sp
print('# velu_check.py (review z7; written by an AI model)')
x, y, a, b = sp.symbols('x y a b')
X = y ** 2 / x ** 2; Y = y * (b - x ** 2) / x ** 2
num = sp.expand(sp.numer(sp.together(Y ** 2 - X * (X ** 2 - 2 * a * X + a ** 2 - 4 * b))))
num = sp.expand(num.subs(y ** 2, x * (x ** 2 + a * x + b)).subs(y ** 4, (x * (x ** 2 + a * x + b)) ** 2)
                .subs(y ** 6, (x * (x ** 2 + a * x + b)) ** 3))
res = sp.expand(sp.Poly(num, y).as_expr().subs(y, sp.sqrt(x * (x ** 2 + a * x + b))))
print('(1) Velu residual:', sp.simplify(res))
m, n, s = sp.symbols('m n s')
S = m ** 2 + n ** 2
P = sp.expand(2 * S * (1 + s) * (S ** 2 + (2 * S ** 2 - 16 * m ** 2 * n ** 2) * s + S ** 2 * s ** 2))
c0, c1, c2, c3 = [P.coeff(s, i) for i in range(4)]
xx, XX = sp.symbols('xx XX')
cub = xx ** 3 + c2 * xx ** 2 + c1 * c3 * xx + c0 * c3 ** 2
print('(2) c0 = %s, c3 = %s; cubic at x = -c3: %s' % (sp.factor(c0), sp.factor(c3), sp.expand(cub.subs(xx, -c3))))
q = sp.Poly(sp.cancel(sp.expand(cub.subs(xx, XX - c3)) / XX), XX)
A, B = q.all_coeffs()[1], q.all_coeffs()[2]
print('    shifted: a = %s, b = %s' % (sp.factor(A), sp.factor(B)))
r = 32 * m * n * S ** 2
print('    16 b - (32 m n S^2)^2 = %s, so 16 b is a square' % sp.expand(16 * B - r ** 2))
assert sp.expand(16 * B - r ** 2) == 0 and sp.expand(A + 32 * m ** 2 * n ** 2 * S) == 0 and sp.expand(B - 64 * m ** 2 * n ** 2 * S ** 4) == 0
Ep = sp.expand(XX ** 2 - 2 * A * XX + A ** 2 - 4 * B)
target = sp.expand((XX - 16 * S * m * n * (m - n) ** 2) * (XX + 16 * S * m * n * (m + n) ** 2))
print("    E' quadratic minus (X - 16Smn(m-n)^2)(X + 16Smn(m+n)^2):", sp.expand(Ep - target))
print('# end')
