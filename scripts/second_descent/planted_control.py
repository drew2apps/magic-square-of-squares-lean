#!/usr/bin/env python3
"""A positive control for the Mordell-Weil sieve with a planted, non-degenerate point (EXPERIMENT; LAST_PAIR.md, section 8).

For every ratio (m, n), t0 = (m + n)/(m - n) is a rational number at which the norms of Q1(-i t0) and Q2(-i t0) are squares, so
Q1(-i t0) = d1 w1^2 and Q2(-i t0) = d2 w2^2 with d1, d2 rational and w1, w2 in K = Q(i)  (Q1 = q+ L1 L3, Q2 = q+ L2 L4).
So the twisted system  Y1^2 = d1 Q1(-i t),  Y2^2 = d2 Q2(-i t)  has the solution t0, which is not degenerate.  It is the system of
mw_sieve.py with other twists; d1 and d2 are different, so the second curve is not the twist of the first, and B = R - conj(R)
is a rational point of F, the Jacobian of Y^2 = -d2 Q2, with its own generator.  Everything else is the same code path:
points over F_p x F_p with every sign, the sum P + P', the difference R - R' carried to F by (x, y) -> (-x, iota y).

Two tests.  (1) Local, for any ratio: the reductions of the true points A in EA(Q) and B in F(Q) lie in the sets that the sieve
computes from t0 modulo p.  (2) Whole sieve, when EA(Q) and F(Q) have rank exactly 1: the pair (a0, b0) of the planted solution,
computed exactly, must survive every prime; the survivors in a box are printed.
usage: planted_control.py m,n [N [p1 p2 ...]]"""
import sys, itertools
from second_descent import *
import cypari2
pari = cypari2.Pari(); pari.allocatemem(2*10**9)
conj = [a for a in K.automorphisms() if a(I) == -I][0]
def setup(m, n):
    t0 = QQ(m + n)/(m - n); fam = family(m, n); (Q1, r1), (Q2, r2) = fam['+13'], fam['+24']; s0 = -I*t0; out = dict(t0=t0, Q=(Q1, Q2), r=(r1, r2))
    ds = []
    for Q in (Q1, Q2):
        z = Q.change_ring(K)(s0); N = QQ(z.norm()); assert N.is_square(), 'norm not a square'; c = N.sqrt()
        q = 2*(c + z[0]) if c + z[0] != 0 else 2*(-c + z[0]); dd = QQ(q).squarefree_part().abs(); assert (z/dd).is_square(); ds.append(dd)
    d1, d2 = ds; out['d'] = (d1, d2)
    EA, A1 = weierstrass(d1*Q1, r1); ER, A2 = weierstrass(d2*Q2, r2); F, A2m = weierstrass(-d2*Q2, r2)
    assert F.a_invariants() == (0, -ER.a2(), 0, ER.a4(), -ER.a6()) and A2m == -A2
    def pt(E, A, Q, r, d):
        EK = E.change_ring(K); v = d*Q.change_ring(K)(s0); Y = v.sqrt(); w = 1/(s0 - r); return EK, EK(A*w, A*Y*w*w)
    EAK, P = pt(EA, A1, Q1, r1, d1); Asum = P + EAK(conj(P[0]), conj(P[1])); Aq = EA(QQ(Asum[0]), QQ(Asum[1]))
    ERK, R = pt(ER, A2, Q2, r2, d2); Bd = R - ERK(conj(R[0]), conj(R[1])); assert Bd[0] in QQ and (Bd[1]/I) in QQ
    Bq = F(-QQ(Bd[0]), QQ(I*Bd[1]))                      # theta(x, y) = (-x, i y): the anti-invariant points of ER(K) are the points of F(Q)
    out.update(EA=EA, ER=ER, F=F, A1=A1, A2=A2, Aq=Aq, Bq=Bq); return out
def fibres(c, p):
    """The sieve's local computation at t0 mod p: all X + Y on EA and all theta(X - Y) on F, with the points X, Y above s = -iota t, +iota t."""
    Fp = GF(p); iota = Fp(-1).sqrt(); (Q1, Q2), (r1, r2), (d1, d2) = c['Q'], c['r'], c['d']
    ea, er, f = c['EA'].change_ring(Fp), c['ER'].change_ring(Fp), c['F'].change_ring(Fp)
    def pts(E, A, q, r, s):
        if s is None:
            v = q.leading_coefficient(); return [E(0, Fp(A)*y) for y in ([v.sqrt(), -v.sqrt()] if v.is_square() else [])]
        if s == Fp(r): return [E(0)]
        v = q(s)
        if not v.is_square(): return []
        w = 1/(s - Fp(r)); y = v.sqrt(); return [E(Fp(A)*w, Fp(A)*y*w*w), E(Fp(A)*w, -Fp(A)*y*w*w)]
    theta = lambda X: f(0) if X.is_zero() else f(-X[0], iota*X[1])
    q1, q2 = (d1*Q1).change_ring(Fp), (d2*Q2).change_ring(Fp)
    def at(t):
        s1, s2 = (None, None) if t is None else (-iota*t, iota*t)
        P1, P2, R1, R2 = pts(ea, c['A1'], q1, r1, s1), pts(ea, c['A1'], q1, r1, s2), pts(er, c['A2'], q2, r2, s1), pts(er, c['A2'], q2, r2, s2)
        if not (P1 and P2 and R1 and R2): return None
        return set(X + Y for X in P1 for Y in P2), set(theta(X - Y) for X in R1 for Y in R2)
    return Fp, ea, f, at
def red(E, P):
    """The reduction of a rational point: the point at infinity when p divides the denominator of x."""
    if P.is_zero() or QQ(P[0]).denominator() % E.base_ring().characteristic() == 0: return E(0)
    return E(P[0], P[1])
def good(c, p):
    num = lambda q: QQ(q).numerator()*QQ(q).denominator()
    return p % 4 == 1 and all(ZZ(x) % p != 0 for x in (c['EA'].discriminant(), c['F'].discriminant(), num(c['A1']), num(c['A2']), num(c['d'][0]), num(c['d'][1]), num(c['r'][0]), num(c['r'][1]), num(c['t0']))) \
        and all(GF(p)(E.j_invariant()) not in (0, 1728) for E in (c['EA'], c['F']))
def gens(E):
    Em = E.minimal_model(); e = pari.ellinit([pari(str(a)) for a in Em.a_invariants()]); r = pari.ellrank(e, 2)
    if not (int(r[0]) == 1 and int(r[1]) == 1 and len(r[3]) == 1): return None
    g = r[3][0]; return Em.isomorphism_to(E)(Em(QQ(str(g[0])), QQ(str(g[1]))))
def coeff(E, G, X):
    """(k, T) with X = k G + T, T torsion; None if X is not of that form with |k| <= 60."""
    tors = E.torsion_points()
    if not any(T.is_zero() for T in tors): tors = tors + [E(0)]
    k0 = int(round((float(X.height())/float(G.height()))**0.5))
    for k in (k0, -k0):
        for T in tors:
            if k*G + T == X: return abs(k), T
    return None
if __name__ == '__main__':
    m, n = (int(v) for v in sys.argv[1].split(',')); N = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    c = setup(m, n); primes = [int(v) for v in sys.argv[3:]] or [p for p in prime_range(13, 400) if good(c, p)][:14]
    print('ratio (%d,%d), planted t0 = %s, twists d1 = %s, d2 = %s' % (m, n, c['t0'], c['d'][0], c['d'][1]))
    ok = 0
    for p in primes:
        assert good(c, p); Fp, ea, f, at = fibres(c, p); t = Fp(c['t0']); sets = at(t)
        assert sets is not None, 'the planted point is not seen modulo %d' % p
        a_ok = red(ea, c['Aq']) in sets[0] or -red(ea, c['Aq']) in sets[0]; b_ok = red(f, c['Bq']) in sets[1] or -red(f, c['Bq']) in sets[1]; ok += a_ok and b_ok
        assert a_ok and b_ok, 'the reduction of the planted point is not in the set of the sieve modulo %d' % p
    print('local test: for the %d primes %s the reductions of A and of B are in the sets computed by the sieve from t0 mod p: %d of %d' % (len(primes), primes, ok, len(primes)))
    GA, GFg = gens(c['EA']), gens(c['F'])
    if GA is None or GFg is None: print('whole sieve not run: EA or F does not have rank exactly 1 with a point found'); sys.exit(0)
    ca, cb = coeff(c['EA'], GA, c['Aq']), coeff(c["F"], GFg, c["Bq"])
    if ca is None or cb is None: print('whole sieve not run: A or B is not a small multiple of the point found (the point found is not a generator, or the multiple is above 60)'); sys.exit(0)
    a0, b0 = ca[0], cb[0]; print('EA and F have rank 1; the planted solution has A = +-(%d G_A + torsion), B = +-(%d G_F + torsion)' % (a0, b0))
    cand = set((a, b) for a in range(N + 1) for b in range(N + 1))
    for p in primes:
        Fp, ea, f, at = fibres(c, p); Sp = set()
        def res(E, G, X):
            tors = [red(E, T) for T in (c['EA'] if E is ea else c['F']).torsion_points()] + [E(0)]; out = set(); Y = E(0); Gp = red(E, G)
            for k in range(Gp.order()):
                if any(Y + T == X or -(Y + T) == X for T in tors): out.add(k)
                Y += Gp
            return out
        oa, ob = red(ea, GA).order(), red(f, GFg).order(); memo_a, memo_b = {}, {}
        for t in list(Fp) + [None]:
            sets = at(t)
            if sets is None: continue
            ra = set(); rb = set()
            for X in sets[0]:
                if X not in memo_a: memo_a[X] = res(ea, GA, X)
                ra |= memo_a[X]
            for X in sets[1]:
                if X not in memo_b: memo_b[X] = res(f, GFg, X)
                rb |= memo_b[X]
            Sp |= set((a, b) for a in ra for b in rb)
        cand = set(x for x in cand if (x[0] % oa, x[1] % ob) in Sp)
        assert (a0, b0) in cand, 'the planted pair is removed by p = %d' % p
    print('whole sieve, pairs 0 <= a, b <= %d, %d primes: the planted pair (%d, %d) survives; survivors: %s' % (N, len(primes), a0, b0, sorted(cand)[:12]))
