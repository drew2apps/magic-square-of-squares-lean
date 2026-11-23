#!/usr/bin/env python3
"""A Mordell-Weil sieve for the last pair (EXPERIMENT; see LAST_PAIR.md, Lemma 6).

Ratio (m,n) with S = m^2+n^2 a square.  C1: Y^2 = d q+ L1 L3 (s), C2: Y^2 = d q+ L2 L4 (s), with Jacobians E1 and E2, and
E2 is the twist of E1 by -1.  E1(Q) has rank 1: generator G, torsion {O, T}.  A solution t gives, with s = -i t,
    P in C1(K), R in C2(K)   (K = Q(i)),     A = P + conj(P) in E1(Q),     B = R - conj(R) in E2^(-1)(Q) = E1(Q),
so A = +-(a G + e T), B = +-(b G + e' T) with integers a, b >= 0.
For a prime p = 1 mod 4 of good reduction, K embeds in F_p x F_p; the reductions of P, conj(P) are the points of C1(F_p) above
s = -iota t and s = +iota t, with square roots whose relative sign is not known.  So (a, b) modulo the order of G mod p lies in
the set S_p computed from all t in P^1(F_p), with every choice of sign.  A pair (a, b) that is not in S_p for some p cannot
come from a solution.  usage: mw_sieve.py m,n d N p1 p2 ..."""
import sys, itertools
from second_descent import *
import cypari2
pari = cypari2.Pari(); pari.allocatemem(4*10**9)
m, n = (int(v) for v in sys.argv[1].split(',')); d = int(sys.argv[2]); N = int(sys.argv[3]); primes = [int(v) for v in sys.argv[4:]]
fam = family(m, n); (Q1, r1), (Q2, r2) = fam['+13'], fam['+24']
E1, A1 = weierstrass(d*Q1, r1); E2, A2 = weierstrass(d*Q2, r2)
Em = E1.minimal_model(); em = pari.ellinit([pari(str(a)) for a in Em.a_invariants()]); rk = pari.ellrank(em, 2)
assert int(rk[0]) == 1 and int(rk[1]) == 1, 'E1 does not have rank exactly 1'
# G: the point of the degenerate solution t = 1 (s = -i), computed exactly; it does not depend on which point PARI returns
conj = [a for a in K.automorphisms() if a(I) == -I][0]; EK1 = E1.change_ring(K)
vd = d*Q1.change_ring(K)(-I); assert vd.is_square(); wd = 1/(-I - r1); Pd = EK1(A1*wd, A1*vd.sqrt()*wd*wd); Gs = Pd + EK1(conj(Pd[0]), conj(Pd[1]))
G = E1(QQ(Gs[0]), QQ(Gs[1])); assert not G.has_finite_order()
T = [P for P in E1.torsion_points() if not P.is_zero()]; assert len(T) == 1; T = T[0]
assert E2.quadratic_twist(-1).is_isomorphic(E1), 'E2 is not the twist of E1 by -1'
Gm = E1.isomorphism_to(Em)(G)
print('E1 has rank 1 (PARI: bounds 1 and 1); G = P + conj(P) at the degenerate solution t = 1, of canonical height %.3f; that G generates E1(Q) modulo torsion is the business of saturation.py' % float(pari.ellheight(em, [pari(str(Gm[0])), pari(str(Gm[1]))])), flush=True)
def red(E, P):
    """The reduction of a rational point: the point at infinity when p divides the denominator of x."""
    if P.is_zero() or QQ(P[0]).denominator() % E.base_ring().characteristic() == 0: return E(0)
    return E(P[0], P[1])
def good(p):
    """p = 1 mod 4, good reduction, the maps s = r + A/x defined modulo p, and no automorphism other than +-1 modulo p."""
    nd = lambda q: QQ(q).numerator()*QQ(q).denominator()
    return p % 4 == 1 and all(ZZ(x) % p != 0 for x in (E1.discriminant(), E2.discriminant(), nd(A1), nd(A2), nd(r1), nd(r2), d)) and GF(p)(E1.j_invariant()) not in (0, 1728)
def local(p):
    assert good(p), 'the prime %d is not allowed' % p
    F = GF(p); iota = F(-1).sqrt(); e1 = E1.change_ring(F); e2 = E2.change_ring(F); q1, q2 = (d*Q1).change_ring(F), (d*Q2).change_ring(F)
    tau = e2.isomorphism_to(e1); Gp = red(e1, G); Tp = red(e1, T); o = Gp.order()        # tau: one of the two isomorphisms, they differ by -1
    dlog = {}; X = e1(0)
    for a in range(o): dlog[X] = a; X += Gp
    def pts(E, A, q, r, s):                              # the points of E above s (s = None is infinity), [] if there is none over F_p
        if s is None:
            v = q.leading_coefficient(); return [E(0, F(A)*y) for y in ([v.sqrt(), -v.sqrt()] if v.is_square() else [])]
        if s == F(r): return [E(0)]
        v = q(s)
        if not v.is_square(): return []
        w = 1/(s - F(r)); y = v.sqrt(); return [E(F(A)*w, F(A)*y*w*w), E(F(A)*w, -F(A)*y*w*w)]
    def res(X):                                           # the residues a mod o with +-(a G + e T) = X
        out = set()
        for e in (0, 1):
            Y = X - e*Tp
            if Y in dlog: out.add(dlog[Y]); out.add((-dlog[Y]) % o)
        return out
    Sp = set()
    for t in list(F) + [None]:
        s1, s2 = (None, None) if t is None else (-iota*t, iota*t)
        P1, P2, R1, R2 = pts(e1, A1, q1, r1, s1), pts(e1, A1, q1, r1, s2), pts(e2, A2, q2, r2, s1), pts(e2, A2, q2, r2, s2)
        if not (P1 and P2 and R1 and R2): continue
        ra = set(); rb = set()
        for X in P1:
            for Y in P2: ra |= res(X + Y)
        for X in R1:
            for Y in R2: rb |= res(tau(X - Y))
        for a in ra:
            for b in rb: Sp.add((a, b))
    return o, Sp
loc = {}
for p in primes:
    loc[p] = local(p); print('p = %d: order of G %d, pairs allowed %d of %d' % (p, loc[p][0], len(loc[p][1]), loc[p][0]**2), flush=True)
# saturation of G at small primes l: G is not l times a rational point (nor l times a point plus T) if this fails modulo some prime
for l in (2, 3, 5, 7, 11, 13, 17, 19):
    ok = False
    for p in prime_range(50, 40000):
        if E1.discriminant() % p == 0 or d % p == 0: continue
        e = E1.change_ring(GF(p)); cnt = e.cardinality()
        if cnt % l: continue
        Gp, Tp = red(e, G), red(e, T)
        if all(not X.division_points(l) for X in ((Gp, Gp + Tp) if l == 2 else (Gp,))): ok = True; break
    print('G is not divisible by %d in E1(Q) modulo torsion: %s' % (l, 'shown modulo %d' % p if ok else 'NOT SHOWN'), flush=True)
def lift(loc, N, say=True):
    """The pairs 0 <= a, b <= N whose residues lie in every local set of loc, and the range actually covered."""
    (o1, S1), (o2, S2) = loc[primes[0]], loc[primes[1]]
    if gcd(o1, o2) != 1:                                  # the two orders are not coprime: lift from the first prime only, in a smaller range
        N = min(N, 30000); cand = set((a, b) for (a0, b0) in S1 for a in range(a0, N + 1, o1) for b in range(b0, N + 1, o1))
        if say: print('the orders of G modulo %d and %d are not coprime: lifting from the first prime only (a simpler path than for coprime orders), range %d' % (primes[0], primes[1], N), flush=True)
        if say: print('pairs 0 <= a, b <= %d allowed by p = %d: %d' % (N, primes[0], len(cand)), flush=True)
        rest = primes[1:]
    else:
        M = o1*o2; inv = pow(o1, -1, o2); N = min(N, M - 1); cand = set()
        crt = lambda x1, x2: x1 + o1*(((x2 - x1)*inv) % o2)
        B2 = {}; A1s = {}
        for (a2, b2) in S2: B2.setdefault(a2, []).append(b2)
        for (a1, b1) in S1: A1s.setdefault(a1, []).append(b1)
        for a1, b1s in A1s.items():
            for a2, b2s in B2.items():
                a = crt(a1, a2)
                if a > N: continue
                for b1 in b1s:
                    for b2 in b2s:
                        b = crt(b1, b2)
                        if b <= N: cand.add((a, b))
        if say: print('pairs 0 <= a, b <= %d allowed by p = %d and %d: %d' % (N, primes[0], primes[1], len(cand)), flush=True)
        rest = primes[2:]
    for p in rest:
        o, Sp = loc[p]; cand = set(c for c in cand if (c[0] % o, c[1] % o) in Sp)
        if say: print('  after p = %d: %d' % (p, len(cand)), flush=True)
    return sorted(cand), N
left, N = lift(loc, N)
# A test of the lifting itself: pairs spread over the range are planted by adding their residues to every local set; the
# lifting has to return each of them, and nothing else but them and what was left before.
plant = sorted(set([(2, 1), (N, 3), (N//2 + 7, N - 5), (31337 % N, 4242421 % N), (loc[primes[0]][0] + 1, N - 2)]))
left2, _ = lift({p: (o, Sp | set((a % o, b % o) for (a, b) in plant)) for p, (o, Sp) in loc.items()}, N, say=False)
print('test of the lifting: %d pairs planted up to %d by adding their residues to the local sets; all of them returned: %s; returned beyond them and the pairs left: %s' % (
    len(plant), N, all(x in left2 for x in plant), sorted(set(left2) - set(plant) - set(left))[:10]))
print('left:', left[:40])
