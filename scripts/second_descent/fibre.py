#!/usr/bin/env python3
"""The solutions of the last pair whose point A (or B) is that of the degenerate solution (EXPERIMENT; with mw_sieve.py).

A = P + conj(P) = c  means  conj(P) = c - P, and conj acts on the curve by s -> -s, that is x -> g(x) = -A0 x/(2 r x + A0).
So x(P - c) = g(x(P)): a polynomial equation for x(P) over K = Q(i), solved for c = +-G, +-G + T (curve 1) and for the points
of E2(K) that correspond to them (curve 2, B = R - conj(R)).  Every rational t that comes out is classified.
usage: fibre.py m,n d"""
import sys
from second_descent import *
import cypari2
pari = cypari2.Pari(); pari.allocatemem(2*10**9)
m, n = (int(v) for v in sys.argv[1].split(',')); d = int(sys.argv[2])
fam = family(m, n); (Q1, r1), (Q2, r2) = fam['+13'], fam['+24']
E1, A1 = weierstrass(d*Q1, r1); E2, A2 = weierstrass(d*Q2, r2)
Em = E1.minimal_model(); em = pari.ellinit([pari(str(a)) for a in Em.a_invariants()]); g = pari.ellrank(em, 2)[3][0]
G = Em.isomorphism_to(E1)(Em(QQ(str(g[0])), QQ(str(g[1])))); T = [P for P in E1.torsion_points() if not P.is_zero()][0]
e1 = E1.change_ring(K); e2 = E2.change_ring(K); iso = e1.isomorphism_to(e2)
def solve(EK, A, r, consts, label):
    R = PolynomialRing(K, 'x'); x = R.gen(); a2 = EK.a2(); f = x**3 + a2*x**2 + EK.a4()*x + EK.a6(); cand = set(); ident = 0
    for c in consts:
        x0, y0 = c[0], c[1]; cand.add(r + A/x0) if x0 != 0 else None
        Dn = 2*r*x + A; hD = (-A*x + (a2 + x + x0)*Dn)*(x - x0)**2
        pol = (f*Dn - hD) if y0 == 0 else (((f + y0**2)*Dn - hD)**2 - 4*y0**2*f*Dn**2)
        if pol == 0: ident += 1; continue
        for rt, _ in pol.roots():
            if rt != 0: cand.add(r + A/rt)
    ts = {}
    for s in cand:
        t = I*s
        if t[1] == 0:
            tq = Fr(int(QQ(t[0]).numerator()), int(QQ(t[0]).denominator())); ts[str(tq)] = classify_t(m, n, tq)
    print('%s: %d constants, equation identically satisfied for %d; values of s found: %d; rational t among them: %s' % (label, len(consts), ident, len(cand), dict(sorted(ts.items()))), flush=True)
c1 = [e1(G), e1(G) + e1(T)]                                  # the equation is squared, so -c is covered with c
solve(e1, A1, r1, c1, 'curve q+ L1 L3, A = +-G or +-G + T')
solve(e2, A2, r2, [iso(c) for c in c1], 'curve q+ L2 L4, B = +-G or +-G + T (as points of E2 over K)')
