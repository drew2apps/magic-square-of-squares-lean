#!/usr/bin/env python3
"""Are the two maps A = P + conj(P) (curve q+ L1 L3) and B = R - conj(R) (curve q+ L2 L4) into E1 independent?  EXPERIMENT.
Modulo a prime p = 1 mod 4, for every t in F_p at which the four values are squares, all choices of sign are tried for a relation
A = +-B + torsion.  Output: p, (values of t tested, values with some relation, by choice of signs).  A handful of hits is what chance gives.
usage: two_maps_independent.py"""
import sys
from second_descent import *
m, n, d = 627, 364, 521266
fam = family(m, n); (Q1, r1), (Q2, r2) = fam['+13'], fam['+24']
E1, A1 = weierstrass(d*Q1, r1); E2, A2 = weierstrass(d*Q2, r2)
def run(p):
    F = GF(p); iota = F(-1).sqrt(); e1 = E1.change_ring(F); e2 = E2.change_ring(F); q1, q2 = Q1.change_ring(F), Q2.change_ring(F)
    isos = [e2.isomorphism_to(e1)]; isos.append(-isos[0])
    def pt(E, A, Q, r, s, Y):
        w = 1/(s - F(r)); return E(F(A)*w, F(A)*Y*w*w)
    tors = [T for T in e1(0).division_points(2)]
    tot = hit = 0; hitsA = {}
    for t in F:
        if t == 0: continue
        s1, s2 = -iota*t, iota*t
        vals = [F(d)*Q(s) for Q in (q1, q2) for s in (s1, s2)]
        if any(v == 0 or not v.is_square() for v in vals) or s1 in (F(r1), F(r2)) or s2 in (F(r1), F(r2)): continue
        y = [v.sqrt() for v in vals]
        P, Pp = pt(e1, A1, Q1, r1, s1, y[0]), pt(e1, A1, Q1, r1, s2, y[1]); R, Rp = pt(e2, A2, Q2, r2, s1, y[2]), pt(e2, A2, Q2, r2, s2, y[3])
        tot += 1; found = set()
        for na, Ac in (('P+sP', P + Pp), ('P-sP', P - Pp)):
            for nb, Bc in (('R-sR', R - Rp), ('R+sR', R + Rp)):
                for iso in isos[:1]:
                    Bi = iso(Bc)
                    if any((Ac - Bi - T).is_zero() or (Ac + Bi - T).is_zero() for T in tors): found.add((na, nb))
        for f in found: hitsA[f] = hitsA.get(f, 0) + 1
        hit += bool(found)
    return tot, hit, hitsA
for p in (10009, 20089):
    if p % 4 == 1: print(p, run(p), flush=True)
