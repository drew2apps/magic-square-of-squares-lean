#!/usr/bin/env python3
"""For which ratios of square type does the rank-1 argument with two maps apply?  EXPERIMENT (SECOND_DESCENT.md, section 6).

For a ratio with m^2 + n^2 a square, E1 is the Jacobian of Y^2 = d q+ L1 L3 (s), and G = P + conj(P) at the degenerate solution
t = 1 is a rational point of E1.  The sieve of mw_sieve.py needs E1(Q) of rank exactly 1 (PARI: both bounds 1) and G of infinite
order.  One line per ratio.  usage: two_map_scope.py pairs123.txt"""
import sys, cypari2
from second_descent import *
pari = cypari2.Pari(); pari.allocatemem(4*10**9)
conj = K.automorphisms()[1] if K.automorphisms()[0](I) == I else K.automorphisms()[0]
tot = r1 = ok = 0
for a in open(sys.argv[1]).read().split():
    m, n = (int(v) for v in a.split(',')); S = m*m + n*n
    if isqrt(S)**2 != S: continue
    tot += 1; Q, r = family(m, n)['+13']; cls = [rational_class(K(v), m, n) for v in (Q(-I), Q(I))]
    if cls[0] is None or cls[0] != cls[1]: print('(%d,%d): class not rational' % (m, n), flush=True); continue
    d = cls[0]; E, A = weierstrass(d*Q, r); EK = E.change_ring(K)
    v = d*Q.change_ring(K)(-I); Y = v.sqrt(); w = 1/(-I - r); P = EK(A*w, A*Y*w*w); G = P + EK(conj(P[0]), conj(P[1])); Gq = E(QQ(G[0]), QQ(G[1]))
    b = pari.ellrank(pari.ellinit([pari(str(c)) for c in E.minimal_model().a_invariants()]), 1); lo, hi = int(b[0]), int(b[1])
    inf = not Gq.has_finite_order(); r1 += (lo == 1 and hi == 1); ok += (lo == 1 and hi == 1 and inf)
    print('(%d,%d): d = %s, rank of E1 in [%d, %d], the point of the degenerate solution has infinite order: %s' % (m, n, d, lo, hi, inf), flush=True)
print('%d ratios of square type; E1 of rank exactly 1: %d; of these with the degenerate point of infinite order (the sieve applies): %d' % (tot, r1, ok))
