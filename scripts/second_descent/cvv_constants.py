#!/usr/bin/env python3
"""The constants of the explicit bound of Checcoli, Veneziano and Viada (Forum Math. Sigma 7 (2019), Theorem 1.2, N = 2) for the
curve E1 of the last pair, and what the sieve of mw_sieve_big.py would need from the curve C (EXPERIMENT, see SECOND_DESCENT.md).

Their statement, as read from the arXiv version (1602.04097) by an automated reader, not checked by a person:  E: y^2 = x^3 + A x + B
with A, B integers, without complex multiplication; C in E^2 irreducible of genus >= 2; if E(Q) has rank 1 then every P in C(Q) has
    h(P) <= C1 h2(C) deg C + C2(E) (deg C)^2 + C3(E),   C1 = 72.251,  C2(E) = C1 (6.019 + 4 c1(E)),  C3(E) = 4 c2(E),
where h is the sum over the two factors of their Neron-Tate height (three times Silverman's, so 1.5 times PARI's ellheight),
deg C and h2(C) are the degree and the normalised height in P^8 (Weierstrass coordinates, Segre embedding), and for E over Q
    c1(E) = min( (log|D| + hinf(j))/4 + hW(j)/8 + log(|A| + |B| + 3)/2 + 2.919,  3 hWW(E) + 4.709 ),
    c2(E) = min( (log|D| + hinf(j))/4 + log(|A| + |B| + 3)/2 + 3.21,  3 hWW(E)/2 + 2.427 ),   hWW(E) = log max(1, |A|^(1/2), |B|^(1/3)).
usage: cvv_constants.py m,n d a_max      (a_max: the sieve leaves no pair (a, b) with max(a, b) <= a_max except (1, 1))"""
import sys
from math import log
from second_descent import *
import cypari2
pari = cypari2.Pari()
m, n = (int(v) for v in sys.argv[1].split(',')); d = int(sys.argv[2]); amax = int(sys.argv[3])
Q1, r1 = family(m, n)['+13']; E1, _ = weierstrass(d*Q1, r1); Em = E1.minimal_model()
c4, c6 = Em.c4(), Em.c6(); A, B = -27*c4, -54*c6; Es = EllipticCurve([A, B]); assert Es.is_isomorphic(E1)
j = Es.j_invariant(); D = -16*(4*A**3 + 27*B**2)
hW = lambda q: log(max(abs(q.numerator()), abs(q.denominator()))); hinf = lambda q: max(log(abs(float(q))), 0.0) if q != 0 else 0.0
hWW = log(max(1.0, float(abs(A))**0.5, float(abs(B))**(1/3)))
c1 = min((log(abs(D)) + hinf(j))/4 + hW(j)/8 + log(abs(A) + abs(B) + 3)/2 + 2.919, 3*hWW + 4.709)
c2 = min((log(abs(D)) + hinf(j))/4 + log(abs(A) + abs(B) + 3)/2 + 3.21, 3*hWW/2 + 2.427)
C1 = 72.251; C2 = C1*(6.019 + 4*c1); C3 = 4*c2
print('E1 in the form y^2 = x^3 + A x + B: log|A| = %.2f, log|B| = %.2f; complex multiplication: %s' % (log(abs(A)), log(abs(B)), Es.has_cm()))
print('c1(E1) = %.3f, c2(E1) = %.3f; C2(E1) = %.1f, C3(E1) = %.1f' % (c1, c2, C2, C3))
hG = 20.007; low = 1.5*hG*(amax + 1)**2
print('sieve: a non-degenerate solution has a or b above %d, so h(P) >= 1.5 * %.3f * %d^2 = %.3e' % (amax, hG, amax + 1, low))
for deg in (96, 200, 1000):
    H = (low - C2*deg*deg - C3)/(C1*deg)
    print('  if deg C <= %d: the bound is below that as soon as h2(C) <= %.3e' % (deg, H))
