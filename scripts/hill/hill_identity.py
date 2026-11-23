#!/usr/bin/env python3
"""Equations (29), (30) of arXiv:2510.08286v2 as identities in six free variables (HILL_NOTE.md).  sympy, exact.

alpha_id and N_i (i = 1, 2, 3) are indeterminates; alpha_in is given by (12); the squares of the beta by (23)-(26).
Result: (30) is the right side of (29) rewritten, and
    LHS(29) - RHS(29) = 4 N_2^2 N_3^2 alpha_3n alpha_3d * ( 2 h(alpha_2) - h(alpha_1) - h(alpha_3) ),   h(a) = a + 1/a,
where the factor in front is a product of quantities that are not zero.  So (29) says exactly that h(alpha_1), h(alpha_2),
h(alpha_3) are in arithmetic progression, and nothing else.  Also the bracket of (30) is alpha_1n alpha_2n - alpha_1d alpha_2d
in disguise: it equals beta_1n^2 alpha_1n^2 - beta_1d^2 alpha_1d^2, with alpha_1n expanded by (12)."""
from sympy import symbols, simplify, factor, cancel, Rational
a1d, a2d, a3d, N1, N2, N3 = symbols('alpha_1d alpha_2d alpha_3d N_1 N_2 N_3', positive=True)
an = lambda ad, N: 4*ad**2/N + 3*ad + N/2
a1n, a2n, a3n = an(a1d, N1), an(a2d, N2), an(a3d, N3)
b1n2, b1d2, b2n2, b2d2 = a2n/a1n, a2d/a1d, a3n/a2n, a3d/a2d
lhs = N3**2*b2n2*b2d2*(8*b1d2**2*a1d**2 - N2**2)**2 - N2**2*(8*b2d2**2*b1d2**2*a1d**2 - N3**2)**2
rhs29 = 4*N2**2*N3**2*b1n2*b1d2*b2n2*b2d2*(a1n - a1d)**2 - 4*N2**2*N3**2*b2n2*b2d2*(b1n2*a1n - b1d2*a1d)**2
rhs30 = 4*N2**2*N3**2*b2n2*b2d2*(b1d2 - b1n2)*(16/N1**2*b1n2*a1d**4 + 24/N1*b1n2*a1d**3 + (13*b1n2 - b1d2)*a1d**2 + 3*N1*b1n2*a1d + N1**2/4*b1n2)
print('RHS(29) - RHS(30) =', cancel(rhs29 - rhs30))
h = lambda a: a + 1/a
A1, A2, A3 = a1n/a1d, a2n/a2d, a3n/a3d
ratio = cancel((lhs - rhs29)/(2*h(A2) - h(A1) - h(A3)))
print('(LHS(29) - RHS(29)) / (2 h(alpha_2) - h(alpha_1) - h(alpha_3)) =', factor(ratio))
print('that quotient minus 4 N_2^2 N_3^2 alpha_3n alpha_3d =', cancel(ratio - 4*N2**2*N3**2*a3n*a3d))
bracket = 16/N1**2*b1n2*a1d**4 + 24/N1*b1n2*a1d**3 + (13*b1n2 - b1d2)*a1d**2 + 3*N1*b1n2*a1d + N1**2/4*b1n2
print('bracket of (30) minus (beta_1n^2 alpha_1n^2 - beta_1d^2 alpha_1d^2) =', cancel(bracket - (b1n2*a1n**2 - b1d2*a1d**2)))
