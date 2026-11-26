#!/usr/bin/env python3
"""x3 referee: independent re-implementation of the F4 Demjanenko-Manin step on C: Y^2 = g4(w) K(w), exact rationals.
Own code throughout (no import of the bundle or of exact_last_pair/sigma_exact); PARI only for the hint G and the
minimal-model hint, both checked here.  Height-difference bounds as quoted in the bundle (Silverman 1990 Thm 1.1 doubled;
CPS for E0 points on a minimal model).  Written by an AI model; not read by a person."""
import sys
from fractions import Fraction as Fr
from math import gcd, isqrt
import cypari2
pari = cypari2.Pari(); pari.allocatemem(2*10**9)

# ---------- exact rational log enclosures
PREC = 2**-160
def _ln_r(r_lo, r_hi):
    """ln on [r_lo, r_hi] subset [1,2], atanh series, rational outward bounds"""
    def lnb(r, upper):
        z = (r - 1) / (r + 1); z2 = z * z; s = Fr(0); term = z; k = 0
        while True:
            s += term / (2 * k + 1); k += 1; term *= z2
            tail = 2 * term / ((2 * k + 1) * (1 - z2))
            if tail < Fr(1, 10**40): break
        return 2 * s + (tail if upper else 0)
    return lnb(r_lo, False), lnb(r_hi, True)
LN2 = _ln_r(Fr(2), Fr(2))
def ln_bounds_int(N):
    """(lo, hi) with lo <= ln N <= hi, N positive integer"""
    assert N >= 1
    k = N.bit_length() - 1
    if k <= 200:
        r = Fr(N, 2**k); lo, hi = _ln_r(r, r)
    else:
        top = N >> (k - 200)                       # top * 2^(k-200) <= N < (top+1) * 2^(k-200)
        lo, _ = _ln_r(Fr(top, 2**200), Fr(top, 2**200)); _, hi = _ln_r(Fr(top + 1, 2**200), Fr(top + 1, 2**200))
    return lo + k * LN2[0], hi + k * LN2[1]
def lnplus_bounds(q):
    q = abs(Fr(q))
    if q <= 1: return Fr(0), Fr(0)
    a = ln_bounds_int(q.numerator); b = ln_bounds_int(q.denominator)
    return max(Fr(0), a[0] - b[1]), a[1] - b[0]
def h_bounds(q):
    q = Fr(q); return ln_bounds_int(max(abs(q.numerator), q.denominator))

# ---------- curve y^2 = x^3 + A2 x^2 + A4 x + A6
class Curve:
    def __init__(s, A2, A4, A6):
        s.A2, s.A4, s.A6 = A2, A4, A6
        b2, b4, b6 = 4 * A2, 2 * A4, 4 * A6; b8 = 4 * A2 * A6 - A4 * A4
        s.b2 = b2; s.D = -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6
        s.c4 = b2 * b2 - 24 * b4; assert s.D != 0
    def f(s, x): return ((x + s.A2) * x + s.A4) * x + s.A6
    def on(s, P): return P is None or P[1] ** 2 == s.f(P[0])
    def add(s, P, Q):
        if P is None: return Q
        if Q is None: return P
        if P[0] == Q[0]:
            if P[1] + Q[1] == 0: return None
            l = (3 * P[0] ** 2 + 2 * s.A2 * P[0] + s.A4) / (2 * P[1])
        else: l = (Q[1] - P[1]) / (Q[0] - P[0])
        x = l * l - s.A2 - P[0] - Q[0]; return (x, l * (P[0] - x) - P[1])
    def silverman(s):
        """(B, U): h(x) - hhat <= B, hhat - h(x) <= U for this integral model (Silverman 1990 Thm 1.1, doubled)"""
        j = Fr(s.c4 ** 3, s.D)
        twomu = ln_bounds_int(abs(s.D))[1] / 6 + lnplus_bounds(j)[1] / 6 + lnplus_bounds(Fr(s.b2, 12))[1] + (LN2[1] if s.b2 else 0)
        hj = h_bounds(j)[1]
        return twomu + hj / 12 + Fr(1946, 1000), twomu + Fr(214, 100)

def sq(q):
    q = Fr(q); return q >= 0 and isqrt(q.numerator) ** 2 == q.numerator and isqrt(q.denominator) ** 2 == q.denominator
def sqrt_q(q): return Fr(isqrt(q.numerator), isqrt(q.denominator))

def classify(m, n, t):
    """global.py classify_t, retyped from its definition"""
    s = Fr(m * m + n * n, 2)
    F0 = m * m + n * n * t * t; F8 = n * n + m * m * t * t; F4 = s * (1 + t * t); L = F4 - 2 * m * n * t; U = F4 + 2 * m * n * t
    if not all(sq(v) for v in (F0, F8, F4, L, U)): return 'no'
    cells = [(m * t + n) ** 2, (m * t - n) ** 2, (n * t + m) ** 2, (n * t - m) ** 2, F0, F8, F4, L, U]
    return 'SOLUTION' if len(set(cells)) == 9 and 0 not in cells else 'degenerate'

def model(m, n):
    S = m * m + n * n
    g4 = [2 * S, 2 * S]; K = [S * S, 2 * S * S - 16 * m * m * n * n, S * S]   # global.py's lists, low degree first
    P = [0] * 4
    for i, u in enumerate(g4):
        for j, v in enumerate(K): P[i + j] += u * v
    c0, c1, c2, c3 = P
    assert c0 == c3 and c1 == c2, 'not palindromic'
    return P, Curve(c2, c1 * c3, c0 * c3 * c3)

def identities(m, n):
    """g4K(t^2) = 16 F4 L U as polynomials in t, and w^3 g4K(1/w) = g4K(w): checked at 20 rational points each
    (degree 6 and 3 respectively, so this proves the identities)."""
    P, _ = model(m, n)
    gk = lambda w: sum(c * w ** i for i, c in enumerate(P))
    s = Fr(m * m + n * n, 2)
    for k in range(1, 21):
        t = Fr(k, 7) - 1; w = t * t
        F4 = s * (1 + t * t); L = F4 - 2 * m * n * t; U = F4 + 2 * m * n * t
        assert gk(w) == 16 * F4 * L * U
        w = Fr(k, 3) + Fr(1, 5)
        assert w ** 3 * gk(1 / w) == gk(w)
    return True

def run(m, n, k=8, extra=1):
    identities(m, n)
    P, C = model(m, n); c3 = P[3]
    E = pari.ellinit([0, C.A2, 0, C.A4, C.A6])
    pts = pari.ellrank(E, 2)[3]
    if len(pts) == 0: pts = pari.ellrank(E, 4)[3]
    G = pari.ellsaturation(E, pts, 300)[0]
    G = (Fr(int(pari.numerator(G[0])), int(pari.denominator(G[0]))), Fr(int(pari.numerator(G[1])), int(pari.denominator(G[1]))))
    assert C.on(G)
    e = -c3; T = (Fr(e), Fr(0)); assert C.on(T)              # w = -1 is the root of g4
    # torsion: odd part by gcd of #C(F_q); 2-part: one rational root, T not in 2C(Q)
    g = 0
    for q in [int(v) for v in pari.primes([3, 600])]:
        if C.D % q == 0: continue
        a2, a4, a6 = C.A2 % q, C.A4 % q, C.A6 % q
        chi = [pow((x ** 3 + a2 * x * x + a4 * x + a6) % q, (q - 1) // 2, q) for x in range(q)]
        N = q + 1 + sum(1 for v in chi if v == 1) - sum(1 for v in chi if v == q - 1)
        g = gcd(g, N)
    g2 = g
    while g % 2 == 0: g //= 2
    assert g == 1, ('odd torsion', g)
    a, b = 3 * e + C.A2, 3 * e * e + 2 * C.A2 * e + C.A4
    assert a * a - 4 * b < 0, 'quadratic factor has real roots'
    if b >= 0 and isqrt(b) ** 2 == b:
        d = isqrt(b)
        assert not any(v >= 0 and isqrt(v) ** 2 == v for v in (a + 2 * d, a - 2 * d)), 'T in 2C(Q)'
    # hhat(G) enclosure
    B, U = C.silverman()
    Q = G
    for _ in range(k): Q = C.add(Q, Q)
    hx = h_bounds(Q[0])
    lo, hi = (hx[0] - B) / 4 ** k, (hx[1] + U) / 4 ** k
    assert lo > 0
    # DM bound: |hhat(P) - hhat(P')| <= B + U + 2 ln|c3| ;  |a| != |b|  =>  2 max - 1 <= that / lo
    lhs = (B + U + 2 * ln_bounds_int(abs(c3))[1]) / lo
    amax = (lhs + 1) // 2
    amax = int(amax) * extra
    cand = set(); Pa = None
    for a in range(amax + 1):
        for R in (Pa, C.add(Pa, T)):
            if R is not None: cand.add(R[0])
        Pa = C.add(Pa, G)
    # |a| = |b|: x(P') in {x(P), x(P + T)}; x(P + T) = e + beta/(x - e), beta = f'(e)
    beta = 3 * e * e + 2 * C.A2 * e + C.A4
    cand.update([Fr(c3), Fr(-c3)])
    qa, qb, qc = e, beta - e * e - c3 * c3, c3 * c3 * e
    dsc = qb * qb - 4 * qa * qc
    if dsc >= 0 and isqrt(dsc) ** 2 == dsc:
        cand.update([Fr(-qb + isqrt(dsc), 2 * qa), Fr(-qb - isqrt(dsc), 2 * qa)])
    # sanity of the |a| = |b| algebra: for a few points P, x(P) x(P+T) relation formula
    for R in (G, C.add(G, G)):
        RT = C.add(R, T); assert RT[0] == e + beta / (R[0] - e)
    ts = {}
    for X in cand:
        w = X / c3
        if sq(w):
            r = sqrt_q(w)
            for t in {r, -r}: ts[t] = classify(m, n, t)
    ts[Fr(0)] = classify(m, n, Fr(0))
    hits = {str(t): v for t, v in ts.items() if v != 'no'}
    return dict(pair=(m, n), hG=(float(lo), float(hi)), BU=(float(B), float(U)), DMrhs=float(B + U + 2 * ln_bounds_int(abs(c3))[1]),
                amax=amax, ncand=len(cand), nt=len(ts), hits=hits,
                gcd2=g2, status='ALERT' if 'SOLUTION' in ts.values() else 'closed', G=G, C=C, T=T, lo=lo, hi=hi)

if __name__ == '__main__':
    extra = 3 if '--x3' in sys.argv else 1
    for a in [x for x in sys.argv[1:] if not x.startswith('--')] or ['7735,4968', '6557,324', '6881,6233']:
        m, n = map(int, a.split(','))
        r = run(m, n, extra=extra)
        print({k: v for k, v in r.items() if k not in ('G', 'C', 'T', 'lo', 'hi')}, flush=True)
