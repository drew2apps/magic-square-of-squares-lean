#!/usr/bin/env python3
"""The Demjanenko-Manin closures of dm_sigma.py, rechecked in exact arithmetic with no Sage and no PARI.

    python3 sigma_exact.py            every line of sigma_evidence.txt

Input: sigma_evidence.txt, one line per closing curve: the ratio, the curve (E+ or E-), the generator G, the
index bound K and the rank bounds with the programs that gave them (Sage/PARI/eclib: these are the inputs this
script does not recompute; see sigma_certify.py).  Everything else is recomputed here from m, n and G alone:

 1. The model.  f = F0 F8 GL (coefficients of degree 0..6, palindromic), g(z) = a3 + a4 z + a5 (z^2 - 2) + a6 (z^3 - 3z)
    so that f(t) = t^3 g(t + 1/t); q(z) = (z + s) g(z), s = 2 for E+ and -2 for E-; z0 = -s; A = q'(z0);
    y^2 = X^3 + a2 X^2 + a4 X + a6 with a2 = q''(z0)/2, a4 = A q'''(z0)/6, a6 = A^2 q''''(z0)/24, X = A/(z - z0).
    Checked: f is palindromic, g reproduces f at 40 values of t, and the point (X, A Y/(z - z0)^2) of a point
    (z, Y) of Y^2 = q(z) lies on the model (at 30 values of z with q(z) a square times a known square).
 2. G lies on the model and is not torsion (its canonical height has a positive lower bound, item 5).
 3. Torsion: the only rational torsion points are O and (e, 0) with e the integer root of the cubic.  The gcd of
    #E(F_q) over good primes 17 <= q < 600 is 2 (torsion of order prime to q injects into E(F_q); Mazur bounds the
    order by 16, so no prime above 16 divides it).
 4. Saturation: for every prime p <= K, a functional E(Q) -> Z/p that vanishes on p E(Q) separates G from
    p E(Q) + torsion (as in satcert.py: a good prime q with p || #E(F_q), P -> (N/p) P mod q in the cyclic
    p-part; for p = 2 two such q with different values on G, since T maps to the nonzero element).  With rank 1
    and index at most K this makes G a generator modulo torsion.
 5. Heights, in rational arithmetic (logarithms enclosed as in height_sides.py):
      h - hhat <= B   (eclib's form of Silverman's bound, enclosed from above),
      hhat - h <= U   (Silverman's archimedean bound as quoted by Cremona-Prickett-Siksek; see height_sides.py),
    hhat(G) enclosed from Q = 2^k G (exact duplication): (h(x(Q)) - B)/4^k <= hhat(G) <= (h(x(Q)) + U)/4^k.
    Demjanenko-Manin: |a^2 - a'^2| hhat(G) <= B + U + log(|A| + 4), so when |a| != |a'|,
    max(|a|, |a'|) <= Xb := (B + U + log(|A| + 4)) / (lower end of hhat(G)); the listing goes to floor(Xb) + 2.
    B uses Silverman's constant 0.973 (eclib and Sage use Bremner's correction 0.961, which gives a smaller B).
    The index bound K of the evidence is checked against sqrt(hhat(G)/lam), with the upper end of the enclosure and the
    recorded Cremona-Siksek lower bound lam (the one floating-point input; sigma_certify.py).
 6. The candidates: X(aG + T) for 0 <= a <= amax and T in {O, (e, 0)} (exact group law); X = 0 (z = infinity);
    P = O (z = z0); P' = O (X = -A/(2 z0)); |a| = |a'|, T'' = O: X = M(X), so X = 0 or X = -A/z0;
    |a| = |a'|, T'' = (e, 0): x(P + T'') = e + f'(e)/(X - e) = M(X) = -A X/(2 z0 X + A), a quadratic in X whose
    rational roots are added.  Each candidate z gives the t with t + 1/t = z, and each t (and 1/t) is tested in
    F0, F8, GL, GU.  Checked: every t with the four corners square is degenerate (t = 0, +-1 or infinity).
"""
import os, sys, math
from fractions import Fraction as Fr
from math import gcd, isqrt

here = os.path.dirname(os.path.abspath(__file__))
QMAX = 4000          # auxiliary primes for the saturation functionals

# ---- rational enclosures of logarithms (as in height_sides.py) ----
def atanh_bounds(z, terms=60):
    s = Fr(0); p = z; z2 = z * z
    for i in range(terms):
        s += p / (2 * i + 1); p *= z2
    return s, s + p / ((2 * terms + 1) * (1 - z2))
LN2 = tuple(2 * x for x in atanh_bounds(Fr(1, 3)))
def ln_bounds(q, bits=96):
    q = Fr(q); assert q > 0
    k = q.numerator.bit_length() - q.denominator.bit_length()
    y = q / Fr(2) ** k
    if y < 1: y *= 2; k -= 1
    one = 1 << bits
    ylo = Fr((y.numerator * one) // y.denominator, one); yhi = Fr(-((-y.numerator * one) // y.denominator), one)
    lo = 2 * atanh_bounds((ylo - 1) / (ylo + 1))[0]; hi = 2 * atanh_bounds((yhi - 1) / (yhi + 1))[1]
    return (lo + k * (LN2[0] if k >= 0 else LN2[1]), hi + k * (LN2[1] if k >= 0 else LN2[0]))
def lnplus_hi(q): return ln_bounds(q)[1] if q > 1 else Fr(0)
def hnaive(x):
    """Enclosure of the logarithmic height of a rational number."""
    x = Fr(x)
    H = max(abs(x.numerator), x.denominator)
    return (Fr(0), Fr(0)) if H == 1 else ln_bounds(H)

# ---- polynomials (lists of integer coefficients, low degree first) ----
def pmul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): c[i + j] += x * y
    return c
def peval(a, x): return sum(c * x**i for i, c in enumerate(a))
def pder(a): return [i * c for i, c in enumerate(a)][1:]

def model(m, n, s):
    S = m*m + n*n
    f = pmul(pmul([m*m, 0, n*n], [n*n, 0, m*m]), [2*S, -8*m*n, 2*S])
    assert len(f) == 7 and all(f[k] == f[6 - k] for k in range(7))
    a = f
    g = [a[3] - 2*a[5], a[4] - 3*a[6], a[5], a[6]]
    for t in [Fr(k, 7) for k in range(1, 41)]:
        assert peval(f, t) == t**3 * peval(g, t + 1/t)
    q = pmul([s, 1], g)
    z0 = -s; assert peval(q, z0) == 0
    d = [q]
    for _ in range(4): d.append(pder(d[-1]))
    v = [peval(p, z0) for p in d]
    A = v[1]; assert A != 0
    assert v[2] % 2 == 0 and (A * v[3]) % 6 == 0 and (A * A * v[4]) % 24 == 0
    a2, a4, a6 = v[2] // 2, A * v[3] // 6, A * A * v[4] // 24
    # the map (z, Y) -> (A/(z - z0), A Y/(z - z0)^2) sends Y^2 = q(z) to the model: check q(z) (A/(z-z0)^2)^2 = cubic(X)
    for z in [Fr(k, 3) for k in range(-15, 16) if Fr(k, 3) != z0]:
        X = Fr(A) / (z - z0)
        assert peval(q, z) * (Fr(A) / (z - z0)**2)**2 == X**3 + a2*X*X + a4*X + a6
    return dict(f=f, g=g, q=q, z0=z0, A=A, a2=a2, a4=a4, a6=a6, S=S)

# ---- the curve y^2 = x^3 + a2 x^2 + a4 x + a6 over Q and over F_q ----
class Curve:
    def __init__(s, a2, a4, a6):
        s.a2, s.a4, s.a6 = a2, a4, a6
        s.b2, s.b4, s.b6 = 4*a2, 2*a4, 4*a6
        s.b8 = 4*a2*a6 - a4*a4
        s.D = -s.b2**2*s.b8 - 8*s.b4**3 - 27*s.b6**2 + 9*s.b2*s.b4*s.b6
        s.c4 = s.b2**2 - 24*s.b4
        assert s.D != 0
    def on(s, P): return P is None or P[1]**2 == P[0]**3 + s.a2*P[0]**2 + s.a4*P[0] + s.a6
    def add(s, P, Q):
        if P is None: return Q
        if Q is None: return P
        if P[0] == Q[0]:
            if P[1] + Q[1] == 0: return None
            l = (3*P[0]**2 + 2*s.a2*P[0] + s.a4) / (2*P[1])
        else:
            l = (Q[1] - P[1]) / (Q[0] - P[0])
        x = l*l - s.a2 - P[0] - Q[0]
        return (x, l*(P[0] - x) - P[1])
    def dup_x(s, x):
        num = x**4 - s.b4*x*x - 2*s.b6*x - s.b8; den = 4*x**3 + s.b2*x*x + 2*s.b4*x + s.b6
        assert den != 0
        return num / den

def is_prime(q):
    if q < 2: return False
    d = 2
    while d * d <= q:
        if q % d == 0: return False
        d += 1
    return True
def leg(a, q):
    a %= q
    return 0 if a == 0 else (1 if pow(a, (q - 1)//2, q) == 1 else -1)
def count_mod(C, q):
    return q + 1 + sum(leg(x**3 + C.a2*x*x + C.a4*x + C.a6, q) for x in range(q))
def red(P, q):
    if P is None: return None
    x, y = P
    return (x.numerator * pow(x.denominator, -1, q) % q, y.numerator * pow(y.denominator, -1, q) % q)
def add_mod(C, P, Q, q):
    if P is None: return Q
    if Q is None: return P
    if P[0] == Q[0]:
        if (P[1] + Q[1]) % q == 0: return None
        l = (3*P[0]*P[0] + 2*C.a2*P[0] + C.a4) * pow(2*P[1], -1, q) % q
    else:
        l = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, q) % q
    x = (l*l - C.a2 - P[0] - Q[0]) % q
    return (x, (l*(P[0] - x) - P[1]) % q)
def mul_mod(C, k, P, q):
    R = None
    while k:
        if k & 1: R = add_mod(C, R, P, q)
        P = add_mod(C, P, P, q); k >>= 1
    return R

def issq_int(x): return x >= 0 and isqrt(x)**2 == x
def issq(x): x = Fr(x); return issq_int(x.numerator) and issq_int(x.denominator)
def four_corners(m, n, t):
    F0 = m*m + n*n*t*t; F8 = n*n + m*m*t*t
    GL = (m+n)**2*(t-1)**2 + (m-n)**2*(t+1)**2; GU = (m-n)**2*(t-1)**2 + (m+n)**2*(t+1)**2
    return all(issq(v) for v in (F0, F8, GL, GU))

def check(m, n, cv, Gx, Gy, K, lam):
    s = 2 if cv == 'E+' else -2
    M = model(m, n, s); A, z0 = M['A'], M['z0']
    C = Curve(M['a2'], M['a4'], M['a6'])
    G = (Fr(Gx), Fr(Gy)); assert C.on(G), 'G not on the model'
    # 3. torsion
    roots = [e for e in divisors_cands(C) if e**3 + C.a2*e*e + C.a4*e + C.a6 == 0]
    assert len(roots) == 1, roots
    e = roots[0]; T = (Fr(e), Fr(0))
    g = 0; qs = 0
    for q in range(17, 600):
        if is_prime(q) and C.D % q:
            g = gcd(g, count_mod(C, q)); qs += 1
    assert g == 2, g
    tors = [None, T]
    # 5. heights
    b2, D, j = C.b2, C.D, Fr(C.c4**3, C.D)
    jsize = max(abs(j.numerator), j.denominator); assert jsize >= 32
    common_hi = ln_bounds(abs(D))[1]/6 + lnplus_hi(abs(j))/6 + lnplus_hi(Fr(abs(b2), 12)) + (LN2[1] if b2 else 0)
    Bhi = common_hi + ln_bounds(jsize)[1]/12 + Fr(1946, 1000)     # Silverman's published 0.973, doubled (eclib uses Bremner's 0.961)
    Uhi = common_hi + Fr(214, 100)
    Bfloat = float(Bhi)
    x = G[0]; k = 0
    while True:
        hx = hnaive(x)
        if hx[1] > 60000 or k >= 12: break
        x = C.dup_x(x); k += 1
    hx = hnaive(x)
    lo = (hx[0] - Bhi) / 4**k; hi = (hx[1] + Uhi) / 4**k
    assert lo > 0
    logK = ln_bounds(abs(A) + 4)[1]
    Xb = (Bhi + Uhi + logK) / lo
    amax = int(Xb) + 2                      # floor, and two more multiples for margin
    # the index bound: if G = k G1 + torsion, hhat(G1) = hhat(G)/k^2 >= lam, so k^2 <= hi/lam; K must reach that
    lam_safe = Fr(lam) * Fr(99999, 100000)  # lam is recorded to 6 digits
    Kmin = 1
    while (Kmin + 1)**2 * lam_safe <= hi: Kmin += 1
    assert K >= Kmin, ('index bound K below sqrt(hhat(G)/lam)', K, Kmin)
    # 4. saturation at every prime p <= K
    assert K * K * Fr(1) >= 0
    Gden = G[0].denominator * G[1].denominator
    for p in [p for p in range(2, K + 1) if is_prime(p)]:
        vals = []
        for q in range(3, QMAX):
            if not is_prime(q) or C.D % q == 0 or Gden % q == 0: continue
            N = count_mod(C, q)
            if N % p or (N // p) % p == 0: continue
            v = mul_mod(C, N // p, red(G, q), q)
            if p > 2:
                if v is not None: vals.append(q); break
            else:
                assert mul_mod(C, N // 2, red(T, q), q) is not None      # T maps to the nonzero element
                vals.append((q, v is not None))
                if len(set(b for _, b in vals)) == 2: break
        assert (p > 2 and vals) or (p == 2 and len(set(b for _, b in vals)) == 2), ('saturation', p)
    # 6. candidates
    cand = set([('z0',)])
    def addX(X):
        cand.add(('inf',) if X == 0 else ('z', Fr(z0) + Fr(A) / X))
    addX(Fr(-A, 2*z0)); addX(Fr(0)); addX(Fr(-A, z0)); addX(T[0])
    aG = None
    for a in range(0, amax + 1):
        for Tt in tors:
            P = C.add(aG, Tt)
            if P is not None: addX(P[0])
        aG = C.add(aG, G)
    # x(P + T) = e + f'(e)/(X - e) = -A X/(2 z0 X + A):  -A X (X - e) = (e (X - e) + fp)(2 z0 X + A)
    fp = 3*e*e + 2*C.a2*e + C.a4
    # coefficients of  -A X^2 + A e X - [ (e X - e^2 + fp)(2 z0 X + A) ]
    c2 = -A - 2*z0*e
    c1 = A*e - (e*A + 2*z0*(fp - e*e))
    c0 = -(fp - e*e)*A
    if c2 == 0 and c1 == 0: assert c0 != 0
    elif c2 == 0: addX(Fr(-c0, c1))
    else:
        disc = c1*c1 - 4*c2*c0
        if issq_int(disc):
            r = isqrt(disc)
            for sg in (1, -1): addX(Fr(-c1 + sg*r, 2*c2))
    four = set()
    for c in cand:
        if c[0] == 'inf':
            ts = [Fr(0)]
            if issq(2*(m*m + n*n)): four.add('inf')
        elif c[0] == 'z0': ts = tsol(Fr(z0))
        else: ts = tsol(c[1])
        for t in ts:
            for tt in ([t, 1/t] if t != 0 else [t]):
                if four_corners(m, n, tt): four.add(str(tt))
    bad = [t for t in four if t not in ('0', '1', '-1', 'inf')]
    return dict(k=k, lo=lo, hi=hi, B=Bhi, U=Uhi, amax=amax, ncand=len(cand), four=sorted(four), bad=bad, tors_q=qs, Kmin=Kmin)

def tsol(z):
    D = z*z - 4
    if D < 0 or not issq(D): return []
    r = Fr(isqrt(D.numerator), isqrt(D.denominator))
    return [(z + r)/2, (z - r)/2]

def divisors_cands(C):
    """Integer roots of the monic cubic: divisors of a6 (a6 != 0 here) found through the rational root theorem
    by a numeric root and an exact check."""
    import cmath
    out = []
    a, b, c = C.a2, C.a4, C.a6
    # Newton-free: integer roots satisfy |e| <= 1 + max(|a|, |b|, |c|); find them as rounded real roots
    for r in real_roots(a, b, c):
        for e in (math.floor(r) - 1, math.floor(r), math.floor(r) + 1, math.floor(r) + 2):
            if e**3 + a*e*e + b*e + c == 0 and e not in out: out.append(e)
    return out

def real_roots(a, b, c):
    """Real roots of x^3 + a x^2 + b x + c by bisection on exact integers (sign changes between critical points)."""
    def f(x): return x**3 + a*x*x + b*x + c
    # critical points of f: 3x^2 + 2a x + b = 0
    dsc = a*a - 3*b
    pts = []
    if dsc > 0:
        r = isqrt(dsc)
        pts = [(-a - r)//3 - 1, (-a + r)//3 + 1]
    Mx = 2 + abs(a) + abs(b) + abs(c)
    edges = sorted(set([-Mx] + pts + [Mx]))
    out = []
    for lo, hi in zip(edges, edges[1:]):
        flo, fhi = f(lo), f(hi)
        if flo == 0: out.append(lo)
        if flo * fhi < 0:
            while hi - lo > 1:
                mid = (lo + hi) // 2
                if f(mid) * flo > 0: lo = mid
                else: hi = mid
            out.append(lo)
    if f(edges[-1]) == 0: out.append(edges[-1])
    return out

def negative():
    """Altered evidence must be refused: G replaced by 2G, 3G or 4G (not saturated),
    a point off the curve, the generator of one curve given for the other."""
    lines = [l for l in open(os.path.join(here, 'sigma_evidence.txt')) if l.strip() and not l.startswith('#')]
    refused = total = 0
    for line in lines:
        w = line.split(); m, n = (int(v) for v in w[0].split(',')); cv = w[1]
        f = dict(x.split('=', 1) for x in w[2:]); Gx, Gy = (Fr(v) for v in f['G'].split(':'))
        M = model(m, n, 2 if cv == 'E+' else -2); C = Curve(M['a2'], M['a4'], M['a6'])
        G = (Gx, Gy); e = [r for r in divisors_cands(C)][0]; T = (Fr(e), Fr(0))
        alts = [C.add(G, G), C.add(C.add(G, G), G), (Gx, Gy + 1), (Gx + 1, Gy)]
        if 2 < int(f['K']): alts.append(C.add(C.add(G, G), C.add(G, G)))      # 4G
        K0 = int(f['K']); p = K0 + 1
        while any(p % q == 0 for q in range(2, p)): p += 1
        Pp = None
        for _ in range(p): Pp = C.add(Pp, G)
        alts.append(Pp)                                                       # pG, p the least prime above K
        other = 'E-' if cv == 'E+' else 'E+'
        for P in alts:
            total += 1
            try:
                r = check(m, n, cv, P[0], P[1], int(f['K']), Fr(f['lam'])); print('NOT REFUSED', m, n, cv)
            except AssertionError:
                refused += 1
        total += 1
        try:
            check(m, n, other, Gx, Gy, int(f['K']), Fr(f['lam'])); print('NOT REFUSED (other curve)', m, n, cv)
        except AssertionError:
            refused += 1
        total += 1
        try:                                                                  # K lowered
            r0 = check(m, n, cv, Gx, Gy, int(f['K']), Fr(f['lam']))
            check(m, n, cv, Gx, Gy, r0['Kmin'] - 1, Fr(f['lam'])); print('NOT REFUSED (K lowered)', m, n, cv)
        except AssertionError:
            refused += 1
    print('negative tests: %d altered records, %d refused' % (total, refused))
    return refused == total


def same_group(path):
    """A fresh sigma_certify.py output against sigma_evidence.txt: PARI may return another generator (-G, G + T, -G + T),
    so the generators are compared up to sign and the point of order 2, exactly, and K and amax must agree."""
    def read(f):
        out = {}
        for line in open(f):
            if line.startswith('#') or not line.strip(): continue
            w = line.split(); d = dict(x.split('=', 1) for x in w[2:])
            out[(w[0], w[1])] = (tuple(Fr(v) for v in d['G'].split(':')), d['K'], d['amax'])
        return out
    old, new = read(os.path.join(here, 'sigma_evidence.txt')), read(path)
    assert set(old) == set(new), (sorted(set(old) ^ set(new)))
    for k in old:
        (G, K, a), (G2, K2, a2) = old[k], new[k]
        assert (K, a) == (K2, a2), (k, K, K2, a, a2)
        m, n = (int(v) for v in k[0].split(',')); M = model(m, n, 2 if k[1] == 'E+' else -2)
        C = Curve(M['a2'], M['a4'], M['a6']); e = divisors_cands(C)[0]; T = (Fr(e), Fr(0))
        neg = lambda P: (P[0], -P[1])
        assert any(G2 == Q for Q in (G, neg(G), C.add(G, T), neg(C.add(G, T)))), ('another group', k)
    print('same generators up to sign and torsion, same K and amax: %d curves' % len(old))


if __name__ == '__main__':
    if len(sys.argv) > 2 and sys.argv[1] == '--same-group':
        same_group(sys.argv[2]); sys.exit(0)
    if '--negative' in sys.argv:
        sys.exit(0 if negative() else 1)
    n_ok = 0
    for line in open(os.path.join(here, 'sigma_evidence.txt')):
        if line.startswith('#') or not line.strip(): continue
        w = line.split(); m, n = (int(v) for v in w[0].split(',')); cv = w[1]
        f = dict(x.split('=', 1) for x in w[2:])
        Gx, Gy = f['G'].split(':')
        r = check(m, n, cv, Fr(Gx), Fr(Gy), int(f['K']), Fr(f['lam']))
        assert not r['bad'], (m, n, cv, r['bad'])
        n_ok += 1
        print('%d,%d %s: hhat(G) in [%.6f, %.6f] (k = %d), B <= %.3f, U <= %.3f, |a| <= %d (dm_sigma.py listed %s), '
              '%d candidates, four corners square only at t in {%s}; saturated at p <= %s (needed: %d); torsion Z/2 (%d primes)'
              % (m, n, cv, r['lo'], r['hi'], r['k'], r['B'], r['U'], r['amax'], f['amax'], r['ncand'], ', '.join(r['four']),
                 f['K'], r['Kmin'], r['tors_q']), flush=True)
    print('sigma_exact OK: %d curves' % n_ok)
