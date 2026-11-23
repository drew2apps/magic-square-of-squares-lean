#!/usr/bin/env python3
"""The Mordell-Weil sieve of mw_sieve.py for a larger range of (a, b), with arrays (EXPERIMENT; see SECOND_DESCENT.md, section 6).
The pairs allowed by the first two primes are lifted to 0 <= a, b < L * (product of the two orders) and filtered by the other primes.
usage: mw_sieve_big.py m,n d L p1 p2 p3 ...

The text below is that of mw_sieve.py.

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
m, n = (int(v) for v in sys.argv[1].split(',')); d = int(sys.argv[2]); L = int(sys.argv[3]); primes = [int(v) for v in sys.argv[4:]]
fam = family(m, n); (Q1, r1), (Q2, r2) = fam['+13'], fam['+24']
E1, A1 = weierstrass(d*Q1, r1); E2, A2 = weierstrass(d*Q2, r2)
Em = E1.minimal_model(); em = pari.ellinit([pari(str(a)) for a in Em.a_invariants()]); rk = pari.ellrank(em, 2)
assert int(rk[0]) == 1 and int(rk[1]) == 1 and len(rk[3]) == 1
g = rk[3][0]; G = Em.isomorphism_to(E1)(Em(QQ(str(g[0])), QQ(str(g[1])))); T = [P for P in E1.torsion_points() if not P.is_zero()]; assert len(T) == 1; T = T[0]
print('E1 has rank 1; generator of canonical height %.3f (PARI; not checked to be saturated: if G = k G0 the sieve below is on multiples of G only)' % float(pari.ellheight(em, g)), flush=True)
def local(p):
    F = GF(p); iota = F(-1).sqrt(); e1 = E1.change_ring(F); e2 = E2.change_ring(F); q1, q2 = (d*Q1).change_ring(F), (d*Q2).change_ring(F)
    tau = e2.isomorphism_to(e1); Gp = e1(G); Tp = e1(T); o = Gp.order()
    dlog = {}; X = e1(0)
    for a in range(o): dlog[X] = a; X += Gp
    def pts(E, A, q, r, s):                              # the points of E above s (s = None is infinity), [] if there is none over F_p
        if s is None:
            v = q.leading_coefficient(); return [E(0, F(A)*y) for y in ([v.sqrt(), -v.sqrt()] if v.is_square() else [])] if True else []
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
import numpy as np, time
loc = {}
for p in primes:
    t0 = time.time(); loc[p] = local(p); print('p = %d: order of G %d, pairs allowed %d of %d' % (p, loc[p][0], len(loc[p][1]), loc[p][0]**2), flush=True)
(o1, S1), (o2, S2) = loc[primes[0]], loc[primes[1]]; assert gcd(o1, o2) == 1; M = o1*o2; inv = pow(o1, -1, o2)
crt = lambda x1, x2: x1 + o1*(((x2 - x1)*inv) % o2)
base = sorted(set((crt(a1, a2), crt(b1, b2)) for (a1, b1) in S1 for (a2, b2) in S2))
A0 = np.array([x[0] for x in base], dtype=np.int64); B0 = np.array([x[1] for x in base], dtype=np.int64)
print('pairs modulo %d allowed by the first two primes: %d; lifting to 0 <= a, b < %d' % (M, len(base), L*M), flush=True)
def table(p):
    o, Sp = loc[p]; T = np.zeros((o, o), dtype=bool)
    for (a, b) in Sp: T[a, b] = True
    return o, T
o3, T3 = table(primes[2]); lift = (M*np.arange(L, dtype=np.int64)); sa = []; sb = []
for c in range(0, len(base), 500):
    a = A0[c:c+500, None] + lift[None, :]; b = B0[c:c+500, None] + lift[None, :]
    hit = T3[(a % o3)[:, :, None], (b % o3)[:, None, :]]
    ii, jj, kk = np.nonzero(hit); sa.append(a[ii, jj]); sb.append(b[ii, kk])
sa = np.concatenate(sa); sb = np.concatenate(sb)
print('  after p = %d: %d' % (primes[2], len(sa)), flush=True)
for p in primes[3:]:
    o, T = table(p); keep = T[sa % o, sb % o]; sa, sb = sa[keep], sb[keep]
    print('  after p = %d: %d' % (p, len(sa)), flush=True)
print('left:', sorted(zip(sa.tolist(), sb.tolist()))[:40])
