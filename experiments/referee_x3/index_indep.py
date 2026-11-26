#!/usr/bin/env python3
"""x3 referee: independent check of saturation (primes <= 199) and of the CPS lower bound hhat >= hi/210^2, own code.
l odd: good q with l | #C(F_q), (N/l) G != O mod q.  l = 2 (different method from the bundle): good q at which the
quadratic factor is irreducible (so C(F_q)[2] = Z/2, the 2-Sylow is cyclic and P in 2C(F_q) iff (N/2) P = O), with
(N/2) G != O and (N/2)(G + T) != O.  Point counts by numpy Legendre tables (exact integers mod q).
Written by an AI model; not read by a person."""
import sys, numpy as np
from fractions import Fraction as Fr
from math import gcd, isqrt
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dm_indep as D
pari = D.pari

def count(A2, A4, A6, q):
    x = np.arange(q, dtype=np.int64)
    v = (((x + A2 % q) % q * x + A4 % q) % q * x + A6 % q) % q
    sqs = np.zeros(q, dtype=np.int8); sqs[(x * x) % q] = 1
    nz = v != 0
    return q + 1 + int(np.sum(np.where(nz, 2 * sqs[v] - 1, 0)))

def red(P, q):
    if P is None: return None
    if P[0].denominator % q == 0: return None
    return (P[0].numerator * pow(P[0].denominator, -1, q) % q, P[1].numerator * pow(P[1].denominator, -1, q) % q)
def addq(C, P, Q, q):
    if P is None: return Q
    if Q is None: return P
    if P[0] == Q[0]:
        if (P[1] + Q[1]) % q == 0: return None
        l = (3 * P[0] ** 2 + 2 * C.A2 * P[0] + C.A4) * pow(2 * P[1], -1, q) % q
    else: l = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, q) % q
    x = (l * l - C.A2 - P[0] - Q[0]) % q; return (x, (l * (P[0] - x) - P[1]) % q)
def mulq(C, k, P, q):
    R = None
    while k:
        if k & 1: R = addq(C, R, P, q)
        P = addq(C, P, P, q); k >>= 1
    return R

def saturation(C, G, T, lmax=199):
    e = int(T[0]); a, b = 3 * e + C.A2, 3 * e * e + 2 * C.A2 * e + C.A4; disc = a * a - 4 * b
    GT = C.add(G, T); den = G[0].denominator * G[1].denominator * GT[0].denominator * GT[1].denominator
    ls = [int(v) for v in pari.primes([2, lmax])]; done = {}
    q = 3
    while len(done) < len(ls):
        q = int(pari.nextprime(q + 1))
        if C.D % q == 0 or den % q == 0: continue
        N = count(C.A2, C.A4, C.A6, q)
        for l in ls:
            if l in done or N % l: continue
            if l == 2:
                # classical criterion (char != 2): P not 2-torsion is in 2C(F_q) iff every x(P) - e_i is a square in the
                # field of definition of e_i; with an irreducible quadratic this reduces (norm argument) to X(P) square.
                chi = lambda v: pow(v % q, (q - 1) // 2, q)
                def in2(P):
                    R = red(P, q)
                    if R is None: return None
                    X = (R[0] - e) % q
                    if X == 0 or R[1] == 0: return None
                    if chi(disc) == q - 1: return chi(X) == 1
                    sd = int(pari.lift(pari.sqrt(pari.Mod(disc, q)))); assert (sd * sd - disc) % q == 0
                    th = [(-a + sd) * pow(2, -1, q) % q, (-a - sd) * pow(2, -1, q) % q]
                    vals = [X] + [(X - z) % q for z in th]
                    if 0 in vals: return None
                    return all(chi(v) == 1 for v in vals)
                r1, r2 = in2(G), in2(GT)
                # cross-check against the order criterion when the 2-Sylow is cyclic
                if chi(disc) == q - 1 and r1 is not None:
                    assert r1 == (mulq(C, N // 2, red(G, q), q) is None)
                if r1 is False and r2 is False: done[l] = q
            else:
                assert (N // l) % 2 == 0
                if mulq(C, N // l, red(G, q), q) is not None: done[l] = q
        assert q < 200000, ('not shown', sorted(set(ls) - set(done)))
    return done

def vp(x, p):
    k = 0
    while x % p == 0: x //= p; k += 1
    return k
def lower_bound(C, need):
    E = pari.ellinit([0, C.A2, 0, C.A4, C.A6])
    Mm = pari('(E) -> my(v, M = ellminimalmodel(E, &v)); [M[1..5], v]')(E)
    u, r, s, t = [Fr(int(pari.numerator(z)), int(pari.denominator(z))) for z in Mm[1]]
    a1, a2, a3, a4, a6 = 0, C.A2, 0, C.A4, C.A6
    # Silverman III.1 Table 3.1, own transcription
    A1 = (a1 + 2 * s) / u; A2 = (a2 - s * a1 + 3 * r - s * s) / u ** 2; A3 = (a3 + r * a1 + 2 * t) / u ** 3
    A4 = (a4 - s * a3 + 2 * r * a2 - (t + r * s) * a1 + 3 * r * r - 2 * s * t) / u ** 4
    A6 = (a6 + r * a4 + r * r * a2 + r ** 3 - t * a3 - t * t - r * t * a1) / u ** 6
    A = [A1, A2, A3, A4, A6]; assert all(c.denominator == 1 for c in A); A1, A2, A3, A4, A6 = [int(c) for c in A]
    assert A == [Fr(int(v)) for v in Mm[0][:5]]
    b2 = A1 * A1 + 4 * A2; b4 = 2 * A4 + A1 * A3; b6 = A3 * A3 + 4 * A6; b8 = A1 * A1 * A6 + 4 * A2 * A6 - A1 * A3 * A4 + A2 * A3 * A3 - A4 * A4
    c4 = b2 * b2 - 24 * b4; c6 = -b2 ** 3 + 36 * b2 * b4 - 216 * b6; Dm = -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6
    assert C.D == Dm * u ** 12
    bad = [int(p) for p in pari.factor(abs(Dm))[0]]
    m = 12; nontriv = []
    for p in bad:
        if not (vp(Dm, p) < 12 or (c4 and vp(c4, p) < 4) or (c6 and vp(c6, p) < 6)): nontriv.append(p)
        if c4 % p: m = m * vp(Dm, p) // gcd(m, vp(Dm, p))
    j = Fr(c4 ** 3, Dm)
    Binf = D.lnplus_bounds(j)[1] / 4 + D.lnplus_bounds(Fr(b2, 12))[1] + (D.LN2[1] if b2 else 0) + Fr(1946, 1000)
    hb = m * m * need + Binf
    X0 = 1
    while D.ln_bounds_int(X0 + 1)[0] <= hb: X0 += 1
    found = []
    for e in range(1, isqrt(X0) + 1):
        e2 = e * e
        for a in range(-X0, X0 + 1):
            if gcd(a, e) != 1: continue
            v = 4 * a ** 3 + b2 * a * a * e2 + 2 * b4 * a * e2 * e2 + b6 * e2 ** 3
            if v > 0 and isqrt(v) ** 2 == v: found.append(Fr(a, e2))
    return dict(m=m, X0=X0, hb=float(hb), found=found, needs_kraus=nontriv, bad=len(bad))

if __name__ == '__main__':
    for arg in [x for x in sys.argv[1:]] or ['7735,4968', '6557,324', '6881,6233']:
        mm, nn = map(int, arg.split(','))
        r = D.run(mm, nn)
        C, G, T, hi = r['C'], r['G'], r['T'], r['hi']
        sat = saturation(C, G, T)
        lb = lower_bound(C, hi / 210 ** 2)
        print(arg, 'saturated at all', len(sat), 'primes <= 199 (largest q used %d)' % max(sat.values()), lb, flush=True)

def kraus2(A):
    """own check: no change [u=2, r, s, t] makes the integral model A = [a1,a2,a3,a4,a6] integral.  By the group law
    [u,r,s,t][1,r',s',t'] = [u, r + u^2 r', s + u s', t + u^2 s r' + u^3 t'] it suffices to take r mod 4, s mod 2,
    t mod 8 (r, s, t are 2-integral when both models are integral, from the b8, a2, a6 relations)."""
    a1, a2, a3, a4, a6 = A; u = 2
    for r in range(4):
        for s in range(2):
            for t in range(8):
                N = [a1 + 2 * s, a2 - s * a1 + 3 * r - s * s, a3 + r * a1 + 2 * t,
                     a4 - s * a3 + 2 * r * a2 - (t + r * s) * a1 + 3 * r * r - 2 * s * t,
                     a6 + r * a4 + r * r * a2 + r ** 3 - t * a3 - t * t - r * t * a1]
                if all(N[i] % u ** k == 0 for i, k in enumerate((1, 2, 3, 4, 6))): return False
    return True
