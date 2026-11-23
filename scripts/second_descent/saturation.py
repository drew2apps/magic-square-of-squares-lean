#!/usr/bin/env python3
"""That G generates E1(Q) modulo torsion (EXPERIMENT; LAST_PAIR.md, Lemma 1).

E1: the Jacobian of Y^2 = d q+ L1 L3 (s);  G = P + conj(P) at the degenerate solution t = 1.
lam: Sage's lower bound for the canonical height of the points of infinite order of E1(Q) (Cremona and Siksek's algorithm,
height_function().min, which is min_gr, the bound for points of good reduction, divided by the square of the lcm of the Tamagawa
exponents and of the number of real components, as the source of Sage 10.8 shows; the height is the one of PARI and Sage, twice Silverman's).  If G = k G0 + torsion then k^2 lam <= h(G).
So it is enough that G is not l times a point plus torsion for the primes l <= sqrt(h(G)/lam).  For such an l a prime p is found
with l dividing the order of E1(F_p) and the reduction of G (for l = 2: of G and of G + T) not divisible by l in E1(F_p).
usage: saturation.py m,n d"""
import sys
from second_descent import *
import cypari2
pari = cypari2.Pari(); pari.allocatemem(4*10**9)
m, n = (int(v) for v in sys.argv[1].split(',')); d = int(sys.argv[2])
Q1, r1 = family(m, n)['+13']; E1, A1 = weierstrass(d*Q1, r1); EK = E1.change_ring(K); conj = [a for a in K.automorphisms() if a(I) == -I][0]
vd = d*Q1.change_ring(K)(-I); wd = 1/(-I - r1); Pd = EK(A1*wd, A1*vd.sqrt()*wd*wd); Gs = Pd + EK(conj(Pd[0]), conj(Pd[1])); G = E1(QQ(Gs[0]), QQ(Gs[1]))
Em = E1.minimal_model(); Gm = E1.isomorphism_to(Em)(G); em = pari.ellinit([pari(str(a)) for a in Em.a_invariants()])
rk = pari.ellrank(em, 2); tors = Em.torsion_points(); T = [P for P in tors if not P.is_zero()]
print('minimal model %s; G = (%s, %s)' % (list(Em.a_invariants()), Gm[0], Gm[1]))
print('rank bounds (PARI ellrank): %d, %d; torsion points other than O: %d' % (int(rk[0]), int(rk[1]), len(T)))
hG = float(pari.ellheight(em, [pari(str(Gm[0])), pari(str(Gm[1]))])); hS = float(Gm.height())
lam = float(Em.height_function().min(.01, 30))
hup = QQ(int(floor(hG*1000)) + 2)/1000; lamlow = QQ(int(floor(lam*10**6)))/10**6     # rationals: above PARI's h(G), below Sage's bound
kmax = isqrt(int(floor(hup/lamlow)))
print('h(G) = %.6f (PARI), %.6f (Sage); lower bound for the height of a point of infinite order (Sage): %.6f; so an index k has k^2 <= %.3f/%.6f, k <= %d' % (hG, hS, lam, float(hup), float(lamlow), kmax))
allok = True
for l in prime_range(2, kmax + 1):
    ok = None
    for p in prime_range(50, 200000):
        if Em.discriminant() % p == 0: continue
        e = Em.change_ring(GF(p))
        if e.cardinality() % l: continue
        Gp = e(Gm); cands = [Gp] + ([Gp + e(t) for t in T] if l == 2 else [])
        if all(not X.division_points(l) for X in cands): ok = p; break
    allok = allok and ok is not None
    print('G is not %d times a point plus torsion: %s' % (l, 'shown modulo %d' % ok if ok else 'NOT SHOWN'))
print('G generates E1(Q) modulo torsion, given the rank bound and the lower bound for the height: %s' % allok)
