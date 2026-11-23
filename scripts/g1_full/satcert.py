#!/usr/bin/env python3
"""Saturation certificates for points on a curve  E: Y^2 = x^3 + a2 x^2 + a4 x = x (x - al)(x - be)  with integers al, be.
Pure Python (no Sage, no PARI): used by dm4.py to produce a certificate and by check_evidence.py to verify it.

Statement certified, for a prime p, points P_1..P_r of E(Q) and a finite set T of points that contains generators of
the torsion subgroup:
    (*)  no point Q of E(Q) has  p Q = c_1 P_1 + ... + c_r P_r + (a torsion point)  with some c_i not divisible by p.
If the P_i are independent and their index in E(Q)/torsion is at most K, then (*) for all primes p <= K says that the
P_i are a basis of E(Q) modulo torsion.

Method.  A homomorphism f: E(Q) -> Z/p that vanishes on p E(Q) is called a functional.  If the matrix of the values
of some functionals on P_1..P_r and on the points of T has rank  r + (rank of the rows of T)  over F_p, then (*) holds:
a relation p Q = sum c_i P_i + t gives sum c_i f(P_i) + f(t) = 0 for every f, a linear dependence of the rows in
which a row of some P_i has a nonzero coefficient, and the rank condition excludes that.
Functionals used, for a prime q of good reduction (q odd, q not dividing al be (al - be)):
  p odd:  if p divides N = #E(F_q) exactly once, the p-part of E(F_q) is cyclic of order p, and
          P -> the discrete logarithm of (N/p) (P mod q) in it;  N is counted naively here;
  p = 2:  P -> the quadratic character modulo q of x(P) - e, for e = 0 and e = al  (the 2-descent map; at a point
          with x = e the value is the character of (e - e')(e - e'') for the other two roots).
A certificate is, for each prime p, the list of primes q used.  The verifier checks itself that every q is an odd
prime of good reduction (trial division), counts #E(F_q) itself, and recomputes every value.

Why a functional vanishes on p E(Q): reduction modulo q is a group homomorphism E(Q) -> E(F_q), since the model has
good reduction at q.  For p odd, (N/p) kills p E(F_q) because N kills E(F_q), and the image lies in the subgroup of
order p, which is cyclic because p divides N only once; the discrete logarithm identifies it with Z/p.  For p = 2 the
map x - e modulo squares is the classical 2-descent homomorphism E(F_q) -> F_q^*/squares, which kills 2 E(F_q).
The torsion subgroup: check_evidence.py confirms that the listed torsion points are all of it by comparing their
number with the gcd of #E(F_q) over good primes q >= 17.  (The torsion subgroup of E(Q) has order at most 16 by
Mazur's theorem, so its order is prime to q, and torsion of order prime to q injects into E(F_q).)"""
from fractions import Fraction as Fr

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b'\0\0'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]

def leg(a, q):
    a %= q
    return 0 if a == 0 else (1 if pow(a, (q - 1) // 2, q) == 1 else -1)

def is_prime(q):
    """Trial division (the auxiliary primes are small)."""
    if q < 2 or (q % 2 == 0 and q != 2): return False
    d = 3
    while d * d <= q:
        if q % d == 0: return False
        d += 2
    return True

def good(al, be, q):
    """q is an odd prime at which Y^2 = x (x - al)(x - be) has good reduction (the model is integral, so it is then
    minimal at q and its reduction is the reduction of the curve)."""
    return q > 2 and is_prime(q) and (al * be * (al - be)) % q != 0

def count(al, be, q):
    """#E(F_q), by summing the quadratic character."""
    return q + 1 + sum(leg(x * (x - al) * (x - be), q) for x in range(q))

def red(P, q):
    """Reduction of a rational point (x, y) (fractions) modulo q; None is the point at infinity."""
    x, y = Fr(P[0]), Fr(P[1])
    if x.denominator % q == 0: return None
    return (x.numerator * pow(x.denominator, -1, q) % q, y.numerator * pow(y.denominator, -1, q) % q)

def add(P, Q, a2, a4, q):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % q == 0: return None
        l = (3*x1*x1 + 2*a2*x1 + a4) * pow(2*y1, -1, q) % q
    else:
        l = (y2 - y1) * pow(x2 - x1, -1, q) % q
    x3 = (l*l - a2 - x1 - x2) % q
    return (x3, (l*(x1 - x3) - y1) % q)

def mul(k, P, a2, a4, q):
    R = None
    while k:
        if k & 1: R = add(R, P, a2, a4, q)
        P = add(P, P, a2, a4, q); k >>= 1
    return R

def on_curve(P, al, be):
    x, y = Fr(P[0]), Fr(P[1]); return y*y == x*(x - al)*(x - be)

def functionals(p, q, al, be, pts, N=None):
    """The values (in Z/p) of the functionals at q on the points pts, as a list of columns; [] if q is not usable."""
    if not good(al, be, q): return []
    a2, a4 = (-(al + be)) % q, (al * be) % q
    R = [red(P, q) for P in pts]
    for P in R: assert P is None or (P[1]*P[1] - P[0]*(P[0] - al)*(P[0] - be)) % q == 0
    if p == 2:
        cols = []
        for e, o1, o2 in ((0, al, be), (al, 0, be)):
            col = []
            for P in R:
                if P is None: col.append(0); continue
                v = (P[0] - e) % q
                if v == 0: v = (e - o1) * (e - o2) % q
                col.append(0 if leg(v, q) == 1 else 1)
            cols.append(col)
        return cols
    if N is None: N = count(al, be, q)
    if N % p != 0 or (N // p) % p == 0: return []
    S = [mul(N // p, P, a2, a4, q) for P in R]
    g = next((P for P in S if P is not None), None)
    if g is None: return [[0] * len(pts)]
    table = {None: 0}; W = None
    for k in range(1, p):
        W = add(W, g, a2, a4, q); assert W is not None and W not in table; table[W] = k
    assert add(W, g, a2, a4, q) is None                  # g has order p
    return [[table[P] for P in S]]

def rank_mod(rows, p):
    M = [[v % p for v in row] for row in rows]; rk = 0
    ncol = len(M[0]) if M else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(M)) if M[i][c]), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]; inv = pow(M[rk][c], -1, p)
        M[rk] = [v * inv % p for v in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c]:
                f = M[i][c]; M[i] = [(a - f * b) % p for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk

def holds(p, qs, al, be, gens, tors, counts=None):
    """True if the functionals at the primes qs prove (*) for p."""
    pts = list(gens) + list(tors); cols = []
    for q in qs: cols += functionals(p, q, al, be, pts, None if counts is None else counts.get(q))
    if not cols: return False
    rows = [[col[i] for col in cols] for i in range(len(pts))]
    return rank_mod(rows, p) == len(gens) + rank_mod(rows[len(gens):], p)

def certificate(al, be, gens, tors, K, qmax=40000):
    """For every prime p <= K a list of primes q that proves (*); None if none is found with q <= qmax."""
    assert all(on_curve(P, al, be) for P in list(gens) + list(tors))
    qs_all = [q for q in primes_upto(qmax) if good(al, be, q)]; counts = {}; cert = {}
    for p in primes_upto(K):
        used = []; cols = []; pts = list(gens) + list(tors); done = False
        for q in qs_all:
            if p != 2:
                if q not in counts: counts[q] = count(al, be, q)
                if counts[q] % p != 0 or (counts[q] // p) % p == 0: continue
            new = functionals(p, q, al, be, pts, counts.get(q))
            if not new: continue
            rows_old = [[col[i] for col in cols] for i in range(len(pts))] if cols else [[] for _ in pts]
            rows_new = [[col[i] for col in cols + new] for i in range(len(pts))]
            if cols and rank_mod(rows_new, p) == rank_mod(rows_old, p): continue          # nothing gained
            cols += new; used.append(q)
            if rank_mod(rows_new, p) == len(gens) + rank_mod(rows_new[len(gens):], p): done = True; break
        if not done: return None
        cert[str(p)] = used
    return cert

def verify(cert, al, be, gens, tors, K):
    """The certificate proves (*) for every prime p <= K (everything is recomputed here)."""
    if not all(on_curve(P, al, be) for P in list(gens) + list(tors)): return False
    return all(str(p) in cert and holds(p, cert[str(p)], al, be, gens, tors) for p in primes_upto(K))
