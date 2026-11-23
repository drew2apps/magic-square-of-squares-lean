#!/usr/bin/env python3
"""A numerical reading of the last step of O. Hill, arXiv:2510.08286v2 (see HILL_NOTE.md).  Exact arithmetic is not
needed: the point is whether equations (29), (30) hold, and whether beta_1 = 1, on magic squares of real numbers.

A 3x3 magic square with entries v, v+D1, v+2D1 | +D2 | three more | +D2 | three more (the form (15) of the paper; every
magic square of positive numbers has it) is a magic square of squares of real numbers.  Each of the three triples
x^2 < y^2 < z^2 with difference D1 gives  n2 = z - y,  alpha = D1 / n2^2,  and, for any normalisation alpha = alpha_n / alpha_d,
N from (11), (12):  alpha_n = 4 alpha_d^2 / N + 3 alpha_d + N / 2."""
import sys, random
from mpmath import mp, mpf, sqrt
mp.dps = 40
def data(v, D1, D2, scale):
    sq = [v, v + D1, v + 2*D1, v + 2*D1 + D2, v + 3*D1 + D2, v + 4*D1 + D2, v + 4*D1 + 2*D2, v + 5*D1 + 2*D2, v + 6*D1 + 2*D2]
    out = []
    for i in range(3):
        y, z = sqrt(sq[3*i + 1]), sqrt(sq[3*i + 2]); n2 = z - y; al = D1/n2**2
        ad = mpf(scale[i]); an = al*ad                                   # any normalisation of alpha = an / ad
        N = (an - 3*ad) - sqrt((an - 3*ad)**2 - 8*ad**2)                 # a root of (12); the other root gives the same check
        assert abs(an - (4*ad**2/N + 3*ad + N/2)) < mpf(10)**-30
        out.append((al, an, ad, N, n2))
    return out
def check(v, D1, D2, scale=(1, 1, 1)):
    (a1, a1n, a1d, N1, m1), (a2, a2n, a2d, N2, m2), (a3, a3n, a3d, N3, m3) = data(mpf(v), mpf(D1), mpf(D2), scale)
    b1n2, b1d2, b2n2, b2d2 = a2n/a1n, a2d/a1d, a3n/a2n, a3d/a2d           # the squares of (23)-(26)
    lhs = N3**2*b2n2*b2d2*(8*b1d2**2*a1d**2 - N2**2)**2 - N2**2*(8*b2d2**2*b1d2**2*a1d**2 - N3**2)**2
    rhs29 = 4*N2**2*N3**2*b1n2*b1d2*b2n2*b2d2*(a1n - a1d)**2 - 4*N2**2*N3**2*b2n2*b2d2*(b1n2*a1n - b1d2*a1d)**2
    rhs30 = 4*N2**2*N3**2*b2n2*b2d2*(b1d2 - b1n2)*(16/N1**2*b1n2*a1d**4 + 24/N1*b1n2*a1d**3 + (13*b1n2 - b1d2)*a1d**2 + 3*N1*b1n2*a1d + N1**2/4*b1n2)
    h = lambda a: a + 1/a
    return dict(eq29=abs(lhs - rhs29)/abs(lhs), eq30=abs(rhs29 - rhs30)/abs(rhs29), middle_terms_in_AP=abs(2*h(a2) - h(a1) - h(a3)),
                beta1=m1/m2, beta1n2_minus_beta1d2=b1n2 - b1d2, alphas=(a1, a2, a3))
if __name__ == '__main__':
    random.seed(1)
    cases = [('the square of 1..9 (roots sqrt 1 .. sqrt 9)', 1, 1, 1, (1, 1, 1))] + [
        ('random real square %d' % k, random.uniform(1, 50), random.uniform(1, 50), random.uniform(1, 50), tuple(random.uniform(1, 9) for _ in range(3))) for k in range(1, 6)]
    for name, v, D1, D2, sc in cases:
        r = check(v, D1, D2, sc)
        ok = lambda x: 'yes' if x < mpf(10)**-30 else 'NO'
        print('%s: (29) holds to 30 digits: %s; (30): %s; h(alpha_i) in arithmetic progression: %s; beta_1 = %.6f; alpha = %.4f, %.4f, %.4f' % (
            name, ok(r['eq29']), ok(r['eq30']), ok(r['middle_terms_in_AP']), float(r['beta1']), *[float(a) for a in r['alphas']]))
