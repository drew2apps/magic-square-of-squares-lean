#!/usr/bin/env python3
"""Does the Jacobian of H: y^2 = (n^2 t^2 + m^2)(t^2 + 1)(S(t^2 + 1) - 4mnt), S = m^2 + n^2, split over Q?
(The quotient of the G1 curve is the twist y^2 = S (...); a quadratic twist changes the Frobenius polynomial by T -> -T
at most, so the answer is the same.)
If it were isogenous over Q to a product of two elliptic curves, the characteristic polynomial of Frobenius
would be a product of two quadratics at every prime of good reduction. This prints its factorisation type at
a few primes for three pairs; a factor of degree 4 shows that the Jacobian does not split over Q.
Needs passagemath (sage.all__sagemath_schemes)."""
from sage.all__sagemath_schemes import *
for m, n in [(72, 65), (724, 665), (3, 2)]:
    S = m*m + n*n; out = []
    for p in [101, 103, 107, 109, 113, 127, 131, 137, 139, 149]:
        if (2*m*n*S*(m*m - n*n)) % p == 0: continue
        R = PolynomialRing(GF(p), 't'); t = R.gen()
        f = (n*n*t**2 + m*m)*(t**2 + 1)*(S*(t**2 + 1) - 4*m*n*t)
        if f.discriminant() == 0: continue
        chi = HyperellipticCurve(f).frobenius_polynomial()
        out.append('%d:%s' % (p, '+'.join(str(g.degree()) for g, e in chi.factor() for _ in range(e))))
    print('(%d,%d)' % (m, n), ' '.join(out))
