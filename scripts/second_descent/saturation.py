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
lam = float(Em.height_function().min(.01, 30))            # Sage returns a double; float() does not change it
_n, _d = lam.as_integer_ratio(); lamQ = QQ(int(_n))/QQ(int(_d))   # the exact value of that double
from height_enclosure import enclose
_lo, _hi, _, _ = enclose(Em, QQ(Gm[0])); assert _lo <= hG <= _hi
hup = QQ(int(_hi.numerator))/QQ(int(_hi.denominator))       # a proved upper bound for h(G) (height_enclosure.py)
lamlow = QQ(floor(lamQ*10**6))/10**6/2                     # HALF of the number Sage returns, rounded down, in exact rationals
kmax = isqrt(int(floor(hup/lamlow)))
print('h(G) = %.6f (PARI), %.6f (Sage); lower bound for the height of a point of infinite order (Sage): %.6f; with half of it, an index k has k^2 <= %.4f/%.6f, k <= %d' % (hG, hS, lam, ceil(float(hup)*10**4)/10**4, float(lamlow), kmax))
D = ZZ(Em.discriminant()); den = ZZ(QQ(Gm[0]).denominator())
def not_divisible(l):
    """A prime p at which the reduction of G is not l times a point (l odd: then G is not l times a point plus torsion either,
    because the torsion has order 2), or None.  p is taken with l dividing the order N of E1(F_p) exactly once, so that
    G is in l E1(F_p) exactly when (N/l) G = O."""
    for p in prime_range(50, 400000):
        if D % p == 0 or den % p == 0: continue
        N = p + 1 - int(pari.ellap(em, p))
        if N % l or (N//l) % l == 0: continue
        e = Em.change_ring(GF(p))
        if not ((N//l)*e(Gm)).is_zero(): return p
    return None
allok = True
for l in prime_range(2, 200):
    if l == 2:
        ok = None
        for p in prime_range(50, 200000):
            if D % p == 0 or den % p == 0: continue
            e = Em.change_ring(GF(p))
            if e.cardinality() % 2: continue
            Gp = e(Gm)
            if all(not X.division_points(2) for X in [Gp] + [Gp + e(t) for t in T]): ok = p; break
    else: ok = not_divisible(l)
    allok = allok and ok is not None
    if l <= max(kmax, 19) or ok is None: print('G is not %d times a point plus torsion: %s' % (l, 'shown modulo %d' % ok if ok else 'NOT SHOWN'))
need = hup/211**2
print('the same for every prime up to 199: %s; so G generates as soon as the points of infinite order have height above %.6f, and Sage\'s bound is %d times that' % (allok, ceil(float(need)*10**6)/10**6, int(floor(lamQ/need))))
print('G generates E1(Q) modulo torsion, given the rank bound and the lower bound for the height: %s' % (allok and kmax < 211))
