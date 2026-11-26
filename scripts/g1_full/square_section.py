#!/usr/bin/env python3
"""Is the class 2 of Sel(E) (SELMER_SQUARE.md, Theorem S1) the image of a rational point at every square-type ratio?
Experiment (task t8), not counted, nothing in the census.  Written by an AI model (Claude); not read by a person.

    python3 square_section.py [--mmax N] [--bound B] [--write] [--negative]
      --mmax N    square-type coprime ratios with m <= N (default 300), both curves E+ and E-
      --bound B   search box |s|, |t| <= B on the covers D_lambda (default 60)
      --write     write square_section.txt and square_section.tsv (otherwise print only)
      --negative  run the forged inputs (output square_section_negative.txt with --write)
      --tsv-to P  also write the .tsv to the path P (for CI; nothing is written into the repository)
      --controls  add the controls C1 (local solver against a box search) and C2 (PARI points, guide)
      --find-hints  PARI's ellrank (a finder only) on the curves left open by the box search: a combination of its
                  points of class 2, written to square_section_hints.json with --write (every hint is checked
                  exactly by the main run, which reads that file and needs no PARI)

Model (SELMER_SQUARE.md): E+ = sigma_exact.model(m, n, 2), E- = model(m, -n, 2); X = x + e, e = -2 S^3, gives
E: y^2 = x(x^2 + a x + b).  Here a, b are recomputed from the model on every curve, and it is checked on every curve
that a = -4 a0 and b = 4 a0^2 + kappa^2 with kappa = 4 m n (m - n)(m + n)^3.

The covering of the class 2:  C_2: N^2 = 2 M^4 + a M^2 e^2 + (b/2) e^4 = 2 (U^2 + V^2) = X1^2 + X2^2 with
U = M^2 - a0 e^2, V = (kappa/2) e^2, X1 = U + V = M^2 - c1 e^2, X2 = U - V = M^2 - c2 e^2, c1 = a0 - kappa/2,
c2 = a0 + kappa/2.  C_2(Q) is non-empty iff the class 2 is alpha of a rational point (x = 2 M^2/e^2).

Second descent (exact).  A rational point has e != 0 (N^2 = 2 M^4 has none); take M, e coprime integers.  Then
g = gcd(X1, X2) divides (c2 - c1) = kappa, and (X1, X2, N) is g times a primitive Pythagorean triple; every rational
point of X1^2 + X2^2 = N^2 is (X1, X2) = lam (s^2 - t^2, 2 s t) with (s : t) rational and lam rational, and the
square class of lam lies in Lambda = { +-d : d squarefree, d | 2 kappa }.  Eliminating M, e:
    F1(s, t) = lam kappa (s^2 - 2 s t - t^2) = (kappa e)^2,
    F2(s, t) = lam kappa (kappa (s^2 - t^2) + c1 (s^2 - 2 s t - t^2)) = (kappa M)^2.
So C_2(Q) != {} implies: for some lam in Lambda there is (s : t) in P^1(Q) with F1, F2 both rational squares
(cover D_lam).  If for every lam in Lambda some completion Q_p (or R) has no (s : t) with F1, F2 both squares in it,
then C_2(Q) = {} and the class 2 is not the image of a rational point: it lies in Sha (C_2 is everywhere locally
soluble, being in the Selmer group).

Local test at p (exact, exhaustive over residue discs of P^1(Z_p)): a disc t = t0 + p^k Z_p (or s = s0 + p^k Z_p,
s0 in pZ) is
  - soluble if at its centre (an integer point) both F1 and F2 are squares in Q_p (0 counts: e = 0 or M = 0 there),
  - insoluble if F1 or F2 has, on the whole disc, one non-square class: F(centre) = u p^r, F(disc) in
    F(centre) + p^k Z_p (integral forms), so for k >= r + 1 (p odd) or k >= r + 3 (p = 2) the class is constant,
  - refined otherwise (cap: depth 60; an undecided disc makes the answer 'undecided', never 'insoluble').
At the real place: the sign of F1 and F2 on the arcs between their real roots, from exact isolating intervals.
Primes of good reduction: if p is odd, p does not divide lam kappa, and F1, F2 are separable mod p without a common root
mod p, then e^2 = F1, M^2 = F2 is a smooth curve of genus 1 over F_p, which has a point (Hasse), and a smooth point
lifts (Hensel); so only p | 2 lam kappa Res(F1, F2) disc(F1) disc(F2) is tested, and these p are found from the
factorisations of m, n, m -+ n, k, H, H- (checked by multiplication), with an exact assert that the product
Res disc disc divides a product of powers of them.

Points: on a cover D_lam that is soluble everywhere, (s, t) with |s|, |t| <= B are tried; a hit gives M, e and the
point X = e + 2 M^2/e^2 of E, checked on the curve exactly, and its class 2 checked exactly.
Hints (since 26 November 2026): square_section_hints.json lists points of class 2 found by PARI's ellrank (effort 1,
R supplied; --find-hints, a finder only) on curves where the box search finds none; each is checked exactly on the
curve and in the class 2, and such a curve must have a cover soluble everywhere (asserted).  Also asserted on every
curve: 2 in Sel(E) (from selmer_square_5000.json), and R on the curve with 2R = P and alpha(R) = 2H.
"""
import os, sys, json, itertools
from fractions import Fraction as Fr
from math import isqrt, gcd
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import sigma_exact as SE
import sympy

WRITE = '--write' in sys.argv
NEG = '--negative' in sys.argv
def arg(name, default):
    return int(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default
MMAX = arg('--mmax', 300)
BOUND = arg('--bound', 60)
out = []
def say(s=''): print(s, flush=True); out.append(s)

def H(m, n): return m**4 - 2*m**3*n + 2*m*m*n*n + 2*m*n**3 + n**4
def issq_int(z): return z >= 0 and isqrt(z)**2 == z
def issq(q):
    q = Fr(q)
    return q >= 0 and issq_int(q.numerator) and issq_int(q.denominator)

def vp(z, p):
    assert z != 0
    v = 0
    while z % p == 0: z //= p; v += 1
    return v, z

def padic_square(z, p):
    """z an integer (or Fraction); True iff z is a square in Q_p (0 is a square)."""
    z = Fr(z)
    if z == 0: return True
    vn, un = vp(z.numerator, p); vd, ud = vp(z.denominator, p)
    if (vn - vd) % 2: return False
    u = un * ud            # same square class as un/ud
    if p == 2: return u % 8 == 1
    return pow(u % p, (p - 1)//2, p) == 1

def form_eval(F, s, t): A, B, C = F; return A*s*s + B*s*t + C*t*t

def content_v(F, p):
    return min(vp(c, p)[0] for c in F if c != 0)

def disc_class(F, s0, t0, p, k):
    """'sq' / 'nsq' if F has one square class on the disc (s0, t0) + p^k (both coordinates moved by p^k Z_p), else None.
    F(disc) lies in F(centre) + p^(k + c) Z_p, c = v_p(content of F)."""
    val = form_eval(F, s0, t0)
    if val == 0: return None
    r, u = vp(val, p)
    if k + content_v(F, p) < r + (3 if p == 2 else 1): return None
    return 'sq' if padic_square(val, p) else 'nsq'

def legendre(z, p): return pow(z % p, (p - 1)//2, p)

def first_discs(F1, F2, p):
    """The residue discs of depth 1 to search.  For p > 200 (a prime of b only, p not dividing 2 lam kappa) two
    exact shortcuts: (i) if some Fi is, mod p, c L^2 with c a non-square (a double root), every disc but the one at the
    root of L has Fi a non-square unit: only that disc is returned; (ii) otherwise a disc centre with F1, F2 both
    non-zero squares mod p (hence squares in Q_p) is looked for by scanning; if found, 'soluble' is returned.
    Neither shortcut is used if p divides the content of F1 or F2 (then all p + 1 discs are searched)."""
    if p <= 200 or any(all(c % p == 0 for c in F) for F in (F1, F2)):
        return [((1, t0), 1, 't') for t0 in range(p)] + [((0, 1), 1, 's')]
    for F in (F1, F2):
        A, B, C = (c % p for c in F)
        assert (A, B, C) != (0, 0, 0)
        if (B*B - 4*A*C) % p == 0:
            if A != 0:
                c = A; rho = (-B * pow(2*A, -1, p)) % p          # F = A (s - rho t)^2 mod p
                root = ((0, 1), 1, 's') if rho == 0 else ((1, pow(rho, -1, p)), 1, 't')
            else:
                c = C; root = ((1, 0), 1, 't')                    # F = C t^2 mod p, root (1 : 0)
            if legendre(c, p) == p - 1:
                return [root]
    for t0 in range(min(p, 5000)):
        v1, v2 = form_eval(F1, 1, t0) % p, form_eval(F2, 1, t0) % p
        if v1 and v2 and legendre(v1, p) == 1 and legendre(v2, p) == 1:
            return 'soluble'
    assert p < 10**6, 'no shortcut at p = %d' % p
    return [((1, t0), 1, 't') for t0 in range(p)] + [((0, 1), 1, 's')]   # exhaustive fallback

def local_soluble(F1, F2, p, maxdepth=60):
    """True / False / None (undecided) : is there (s : t) in P^1(Q_p) with F1, F2 squares in Q_p?"""
    undecided = False
    # divide each form by the largest even power of p in its content (square classes unchanged)
    F1 = tuple(c // p**(2*(content_v(F1, p)//2)) for c in F1)
    F2 = tuple(c // p**(2*(content_v(F2, p)//2)) for c in F2)
    # chart 1: (1, t), t in Z_p ; chart 2: (s, 1), s in pZ_p
    stack = first_discs(F1, F2, p)
    if stack == 'soluble': return True
    while stack:
        (s0, t0), k, ch = stack.pop()
        if padic_square(form_eval(F1, s0, t0), p) and padic_square(form_eval(F2, s0, t0), p):
            return True
        c1 = disc_class(F1, s0, t0, p, k); c2 = disc_class(F2, s0, t0, p, k)
        if c1 == 'nsq' or c2 == 'nsq': continue
        if k >= maxdepth: undecided = True; continue
        for j in range(p):
            if ch == 't': stack.append(((1, t0 + j * p**k), k + 1, 't'))
            else: stack.append(((s0 + j * p**k, 1), k + 1, 's'))
    return None if undecided else False

def real_soluble(F1, F2):
    """Exact: is there (s : t) in P^1(R) with F1 >= 0 and F2 >= 0 ?  F1, F2 have together four distinct simple real
    roots in u = s/t (asserted: P1 P2 is squarefree of degree 4 with 4 real roots), so P^1(R) minus the roots is four
    open arcs on which both signs are constant.  A root of Fi with Fj > 0 there has an arc next to it with Fi > 0 and
    Fj > 0, so the answer is: does one test point per arc give F1 > 0 and F2 > 0.  Test points: rational numbers
    between consecutive exact isolating intervals, and one beyond the last (that arc contains (1 : 0))."""
    u = sympy.symbols('u')
    P1 = sympy.Poly(F1[0]*u**2 + F1[1]*u + F1[2], u); P2 = sympy.Poly(F2[0]*u**2 + F2[1]*u + F2[2], u)
    prod = P1 * P2
    assert prod.degree() == 4 and sympy.gcd(prod, prod.diff(u)).degree() == 0
    ivs = [(Fr(str(lo)), Fr(str(hi))) for (lo, hi), mult in prod.intervals(eps=Fr(1, 10**9))]
    assert len(ivs) == 4
    ivs.sort()
    pts = []
    for i in range(3):
        assert ivs[i][1] < ivs[i + 1][0]
        pts.append((ivs[i][1] + ivs[i + 1][0]) / 2)
    pts.append(ivs[3][1] + 1)
    return any(form_eval(F1, q, 1) > 0 and form_eval(F2, q, 1) > 0 for q in pts)

def curve_data(m, nn):
    """E at (m, nn) (nn signed): a, b recomputed from the model; a0, kappa, c1; the 2-torsion abscissa e."""
    M = SE.model(m, nn, 2)
    S = m*m + nn*nn; k = isqrt(S); assert k*k == S, 'not of square type'
    e = -2 * S**3
    assert e**3 + M['a2']*e*e + M['a4']*e + M['a6'] == 0
    a = 3*e + M['a2']; b = 3*e*e + 2*M['a2']*e + M['a4']
    # X^3 + a2 X^2 + a4 X + a6 = x^3 + a x^2 + b x at X = x + e
    for x in range(-3, 4):
        X = x + e
        assert X**3 + M['a2']*X*X + M['a4']*X + M['a6'] == x**3 + a*x*x + b*x
    assert a % 4 == 0
    a0 = -a // 4
    kappa = 4*m*nn*(m - nn)*(m + nn)**3
    assert b == 4*a0*a0 + kappa*kappa, 'b = 4 a0^2 + kappa^2 fails'
    assert b == 4 * S*S * H(m, nn) * H(m, -nn)
    assert kappa % 2 == 0 and b % 2 == 0
    c1 = a0 - kappa // 2
    return dict(M=M, S=S, k=k, e=e, a=a, b=b, a0=a0, kappa=kappa, c1=c1)

def primes_of(m, nn, k):
    """odd and even primes of 2 m n (m-n)(m+n) k H H-, with an exact check of the factorisations."""
    ps = set([2])
    for z in (m, nn, m - nn, m + nn, k, H(m, nn), H(m, -nn)):
        z = abs(z)
        f = sympy.factorint(z)
        prod = 1
        for p, ex in f.items():
            assert sympy.isprime(p); prod *= p**ex
        assert prod == z
        ps |= set(f)
    return sorted(ps)

def covers(cd):
    kappa, c1 = cd['kappa'], cd['c1']
    Q1 = (1, -2, -1)                                   # s^2 - 2 s t - t^2
    G = (kappa + c1, -2*c1, -kappa - c1)               # kappa (s^2 - t^2) + c1 Q1
    return Q1, G

def lambdas(kappa):
    ps = sorted(set(sympy.factorint(abs(kappa))) | {2})
    res = []
    for bits in itertools.product((0, 1), repeat=len(ps)):
        d = 1
        for p, bt in zip(ps, bits):
            if bt: d *= p
        res += [d, -d]
    return res

def second_descent(m, nn, bound=BOUND, want_points=True):
    cd = curve_data(m, nn)
    kappa, c1 = cd['kappa'], cd['c1']
    Q1, G = covers(cd)
    # Res and discriminants, exact: the bad primes of the covers lie in primes_of
    s_, t_ = sympy.symbols('s t')
    q1 = Q1[0]*s_**2 + Q1[1]*s_*t_ + Q1[2]*t_**2; g = G[0]*s_**2 + G[1]*s_*t_ + G[2]*t_**2
    res = sympy.resultant(q1.subs(t_, 1), g.subs(t_, 1), s_)
    assert res == -4*kappa**2 and G[0] != 0, (res, kappa)
    dG = G[1]**2 - 4*G[0]*G[2]
    assert dG == 2 * cd['b'], 'disc of G is 2b'
    plist = primes_of(m, nn, cd['k'])
    rad = 1
    for p in plist: rad *= p
    for z in (res, dG, 8):
        z = abs(z)
        while gcd(z, rad) > 1: z //= gcd(z, rad)
        assert z == 1                                   # every bad prime of the covers is in plist
    survivors = []
    reasons = {}
    realc = {}
    for lam in lambdas(kappa):
        F1 = tuple(lam*kappa*c for c in Q1); F2 = tuple(lam*kappa*c for c in G)
        why = None
        sg = lam > 0
        if sg not in realc: realc[sg] = real_soluble(F1, F2)   # F1, F2 scale by |lam|: only the sign matters
        if not realc[sg]: why = 'R'
        else:
            for p in plist:
                r = local_soluble(F1, F2, p)
                if r is None: why = 'undecided at %d' % p; break
                if r is False: why = 'Q_%d' % p; break
        if why is None: survivors.append(lam)
        reasons[lam] = why
    point = None
    if want_points and survivors:
        for lam in survivors:
            F1 = tuple(lam*kappa*c for c in Q1); F2 = tuple(lam*kappa*c for c in G)
            for s in range(0, bound + 1):
                for t in range(-bound, bound + 1):
                    if gcd(s, t) != 1: continue
                    v1, v2 = form_eval(F1, s, t), form_eval(F2, s, t)
                    if v1 <= 0 or v2 < 0: continue
                    if issq_int(v1) and issq_int(v2):
                        ke, kM = isqrt(v1), isqrt(v2)
                        P = point_from(cd, Fr(kM, kappa), Fr(ke, kappa))
                        if P is not None:
                            point = (lam, s, t, P); break
                if point: break
            if point: break
    return cd, reasons, survivors, point

def point_from(cd, M, e):
    """the point of E (original model) with x = 2 M^2 / e^2, checked exactly; None if e = 0."""
    if e == 0: return None
    x = 2 * M*M / (e*e)
    X = x + cd['e']
    Mo = cd['M']
    rhs = X**3 + Mo['a2']*X*X + Mo['a4']*X + Mo['a6']
    assert issq(rhs), 'the point is not on the curve'
    assert issq((X - cd['e']) / 2), 'class of the point is not 2'
    return (X, Fr(isqrt(rhs.numerator), isqrt(rhs.denominator)))

def check_class2(cd, X, Y):
    """a point (X, Y) given from outside (a hint): on the original model and of class 2, both exactly."""
    X, Y = Fr(X), Fr(Y); Mo = cd['M']
    assert Y*Y == X**3 + Mo['a2']*X*X + Mo['a4']*X + Mo['a6'], 'the hint is not on the curve'
    assert X != cd['e'] and issq((X - cd['e']) / 2), 'class of the hint is not 2'
    return (X, Y)

def check_R(m, nn, cd):
    """the point R of UNIFORM_RANK.md (identity A5) on this curve: on the model, 2R = P (exact doubling), and
    alpha(R) = 2SH = 2H mod squares (S is a square), which lies in K."""
    Mo = cd['M']; S = cd['S']
    xR = Fr(-4*m*nn*S*(m*m - nn*nn)); yR = Fr(-8*m*nn*(m + nn)*S*H(m, nn))
    assert yR*yR == xR**3 + Mo['a2']*xR*xR + Mo['a4']*xR + Mo['a6'], 'R is not on the curve'
    lam = (3*xR*xR + 2*Mo['a2']*xR + Mo['a4']) / (2*yR)
    x2 = lam*lam - Mo['a2'] - 2*xR; y2 = lam*(xR - x2) - yR
    xP = Fr(-(m + nn)**2 * S*S); yP = Fr(S**3 * (m - nn) * (m + nn)**2)
    assert x2 == xP and y2*y2 == yP*yP, '2R is not P'
    aR = xR - cd['e']; assert aR.denominator == 1
    assert sqfree(int(aR)) == sqfree(2*H(m, nn)), 'alpha(R) is not 2H'
    return True

def load_hints():
    path = os.path.join(here, 'square_section_hints.json')
    if not os.path.exists(path): return {}
    return {(m, nn): (Fr(X), Fr(Y)) for m, nn, X, Y in json.load(open(path))['points']}

def sqfree(z):
    sg = -1 if z < 0 else 1
    f = sympy.factorint(abs(z)); r = 1
    for p, ex in f.items():
        if ex % 2: r *= p
    return sg * r

def two_squares(d):
    """alpha, beta with alpha^2 + beta^2 = d (exact search; d squarefree, a sum of two squares)."""
    for al in range(isqrt(d) + 1):
        r = d - al*al
        if issq_int(r): return al, isqrt(r)
    raise AssertionError('%d is not a sum of two squares' % d)

def descent_class(m, nn, d, bound=BOUND, want_points=True):
    """Second descent on C_d: N^2 = d M^4 + a M^2 e^2 + (b/d) e^4 for a positive squarefree d (a sum of two squares).
    d N^2 = X1^2 + X2^2 with X1 = d M^2 - 2 a0 e^2, X2 = kappa e^2; d = alpha^2 + beta^2; every rational point is
    X1 + i X2 = lam (alpha + i beta)(s + i t)^2 (Hilbert 90), lam in +-(squarefree divisors of 2 kappa d) (the proof:
    for p not dividing 2 d kappa, v_p(lam) odd would force p | gcd(X1, X2), which divides d kappa for coprime M, e).
    Cover: F1 = lam kappa q2 = (kappa e)^2, F2 = d lam kappa (kappa q1 + 2 a0 q2) = (d kappa M)^2 with
    q1 = alpha (s^2 - t^2) - 2 beta s t, q2 = beta (s^2 - t^2) + 2 alpha s t."""
    cd = curve_data(m, nn)
    kappa, a0, b = cd['kappa'], cd['a0'], cd['b']
    assert d > 0 and sqfree(d) == d and b % d == 0
    al, be = two_squares(d)
    q1 = (al, -2*be, -al); q2 = (be, 2*al, -be)
    G = tuple(kappa*x + 2*a0*y for x, y in zip(q1, q2))
    s_ = sympy.symbols('s')
    res = sympy.resultant(sympy.Poly(q2[0]*s_**2 + q2[1]*s_ + q2[2], s_), sympy.Poly(G[0]*s_**2 + G[1]*s_ + G[2], s_))
    plist = primes_of(m, nn, cd['k'])
    rad = 1
    for p in plist: rad *= p
    for z in (res, q2[1]**2 - 4*q2[0]*q2[2], G[1]**2 - 4*G[0]*G[2], d):
        z = abs(int(z)); assert z != 0
        while gcd(z, rad) > 1: z //= gcd(z, rad)
        assert z == 1                                   # every bad prime of the covers is in plist
    ps = sorted(set(sympy.factorint(abs(kappa))) | set(sympy.factorint(d)) | {2})
    lams = []
    for bits in itertools.product((0, 1), repeat=len(ps)):
        x = 1
        for p, bt in zip(ps, bits):
            if bt: x *= p
        lams += [x, -x]
    reasons = {}; survivors = []; realc = {}
    for lam in lams:
        F1 = tuple(lam*kappa*c for c in q2); F2 = tuple(d*lam*kappa*c for c in G)
        why = None
        sg = lam > 0
        if sg not in realc: realc[sg] = real_soluble(F1, F2)
        if not realc[sg]: why = 'R'
        else:
            for p in plist:
                r = local_soluble(F1, F2, p)
                if r is None: why = 'undecided at %d' % p; break
                if r is False: why = 'Q_%d' % p; break
        if why is None: survivors.append(lam)
        reasons[lam] = why
    point = None
    if want_points:
        for lam in survivors:
            F1 = tuple(lam*kappa*c for c in q2); F2 = tuple(d*lam*kappa*c for c in G)
            for s in range(0, bound + 1):
                for t in range(-bound, bound + 1):
                    if gcd(s, t) != 1: continue
                    v1, v2 = form_eval(F1, s, t), form_eval(F2, s, t)
                    if v1 <= 0 or v2 < 0: continue
                    if issq_int(v1) and issq_int(v2):
                        e = Fr(isqrt(v1), kappa); Mq = Fr(isqrt(v2), d*kappa)
                        x = d * Mq*Mq / (e*e); X = x + cd['e']; Mo = cd['M']
                        rhs = X**3 + Mo['a2']*X*X + Mo['a4']*X + Mo['a6']
                        assert issq(rhs) and issq((X - cd['e']) / d)
                        point = (lam, s, t, X); break
                if point: break
            if point: break
    return dict(d=d, reasons=reasons, survivors=survivors, point=point, ncovers=len(lams))

def cosets(m, nn, sel):
    """the cosets of K = <alpha(T), alpha(R)> = {1, H H-, 2H, 2H-} in Sel(E), each as its least positive element."""
    h, hm = H(m, nn), H(m, -nn)
    K = sorted(set(sqfree(c) for c in (1, h*hm, 2*h, 2*hm)))
    assert len(K) == 4 and all(c in sel for c in K)
    seen = set(K); reps = []
    for c in sorted(sel, key=abs):
        if c in seen: continue
        cs = set(sqfree(c*k) for k in K)
        assert cs <= set(sel)
        seen |= cs; reps.append(min(cs, key=abs))
    return K, reps

def rank_one(m, nn, sel, selp, bound=BOUND):
    """alpha(E(Q)) = K exactly if every coset of K in Sel(E) other than K has a representative d with C_d(Q) empty
    (all covers locally insoluble).  With Sel(E') = {1, -1} = alpha'(<T'>) this gives rank E(Q) = 1 (2^r = 4*2/4)."""
    K, reps = cosets(m, nn, sel)
    results = [descent_class(m, nn, d, bound=bound) for d in reps]
    allsha = all(not r['survivors'] for r in results)
    return K, reps, results, allsha

def square_type_ratios(mmax):
    rs = []
    for m in range(2, mmax + 1):
        for n in range(1, m):
            if gcd(m, n) == 1 and (m - n) % 2 == 1:
                S = m*m + n*n
                if isqrt(S)**2 == S: rs.append((m, n))
    return rs

def negative():
    res = []
    def expect_fail(label, fn):
        try:
            fn(); res.append(False); say('  NOT REFUSED: ' + label)
        except AssertionError:
            res.append(True); say('  refused: ' + label)
    say('Negative tests')
    # 1. a ratio not of square type
    expect_fail('ratio (5,2), not of square type', lambda: curve_data(5, 2))
    # 2. a forged point: x = 2 M^2/e^2 with M, e not from a cover
    cd = curve_data(12, 5)
    expect_fail('forged point M = 1, e = 1 on E+ (12,5)', lambda: point_from(cd, Fr(1), Fr(1)))
    # 3. a forged kappa (off by a factor 3): b = 4 a0^2 + kappa^2 must fail
    def forged_kappa():
        d = curve_data(12, 5)
        assert d['b'] == 4*d['a0']**2 + (3*d['kappa'])**2
    expect_fail('kappa multiplied by 3', forged_kappa)
    # 4. the local solver must find the point of a soluble pair, and refuse an insoluble one
    def solver_insoluble():
        # F1 = 3 (s^2 + t^2) is not a square at any point of P^1(Q_3) (3 (s^2 + t^2) has odd valuation or class 3u), so
        # the solver must answer False; claiming True is refused
        assert local_soluble((3, 0, 3), (1, 0, 1), 3) is True
    expect_fail('local solver claims a point on 3(s^2+t^2) = square over Q_3', solver_insoluble)
    def solver_soluble():
        assert local_soluble((1, 0, 0), (1, 0, 1), 5) is False
    expect_fail('local solver claims no point on s^2 = sq, s^2 + t^2 = sq over Q_5', solver_soluble)
    # 5. a class known to be the image of a point must survive somewhere: on E- of (12,5) PARI (guide) gives rank 2
    #    and the class 2 is the image of a point (found below); claiming all covers insoluble there is refused
    def forged_sha():
        cd, reasons, surv, pt = second_descent(12, -5, bound=40)
        assert pt is not None
        assert not surv, 'covers survive'
    expect_fail('E- (12,5): "all covers insoluble" (a point of class 2 exists)', forged_sha)
    # 6. alpha(E(Q)) = K claimed where a point of class 2 exists (E- of (4,3)): refused
    def forged_rank_one():
        rec = {(r['m'], r['n']): r for r in json.load(open(os.path.join(here, 'selmer_square_5000.json')))}[(4, -3)]
        K, reps, res, allsha = rank_one(4, -3, rec['S'], rec['Sp'], bound=10)
        assert any(r['point'] for r in res)
        assert allsha, 'a coset survives'
    expect_fail('E- (4,3): "alpha(E(Q)) = K" (a point of class 2 exists)', forged_rank_one)
    # 7. a class d that is not in Sel(E) as a sum of two squares: 3 is refused by two_squares
    expect_fail('class 3 on E+ (12,5) (not a sum of two squares)', lambda: descent_class(12, 5, 3))
    # 8. the cosets of K must lie in the Selmer group: a forged Selmer list missing H H- is refused
    def forged_sel():
        rec = {(r['m'], r['n']): r for r in json.load(open(os.path.join(here, 'selmer_square_5000.json')))}[(15, 8)]
        h, hm = H(15, 8), H(15, -8)
        sel = [c for c in rec['S'] if c != sqfree(h*hm)]
        cosets(15, 8, sel)
    expect_fail('Sel(E) of E+ (15,8) with the class of T removed', forged_sel)
    # 9. real place: -s^2 + t^2 >= 0 needs |s/t| <= 1, s^2 - 4 t^2 >= 0 needs |s/t| >= 2: no real point
    def real_forged():
        assert real_soluble((-1, 0, 4), (1, 0, -1)) is True          # control: 1 <= |u| <= 2 is a real point
        assert real_soluble((-1, 0, 1), (1, 0, -4)) is True          # forged claim
    expect_fail('real solver claims a point on -s^2 + t^2 >= 0, s^2 - 4t^2 >= 0', real_forged)
    # 11. a hint must be on the curve and of class 2: R of E+ (12,5) (on the curve, class 2H) is refused,
    #     and so is a hint moved off the curve
    def forged_hint_R():
        cd = curve_data(12, 5); S = cd['S']
        check_class2(cd, Fr(-4*12*5*S*(144 - 25)), Fr(-8*12*5*17*S*H(12, 5)))
    expect_fail('hint R on E+ (12,5) offered as a point of class 2 (its class is 2H)', forged_hint_R)
    def forged_hint_off():
        cd = curve_data(4, -3)
        check_class2(cd, Fr(-1968), Fr(1))
    expect_fail('hint (X, Y) = (-1968, 1) on E- (4,3) (X is the class-2 point, Y is not)', forged_hint_off)
    # 12. R must be on the curve with 2R = P: R of (12,5) checked on the curve of (12,-5) is refused
    def forged_R():
        cd = curve_data(12, -5); S = cd['S']; Mo = cd['M']
        xR = Fr(-4*12*5*S*(144 - 25)); yR = Fr(-8*12*5*17*S*H(12, 5))
        assert yR*yR == xR**3 + Mo['a2']*xR*xR + Mo['a4']*xR + Mo['a6'], 'R is not on the curve'
    expect_fail('R of E+ (12,5) on E- (12,5)', forged_R)
    say('negative tests refused: %d of %d' % (sum(res), len(res)))
    assert all(res)

def brute_control(ntests=300, seed=1):
    """Control of the local solver (not part of any proof): random pairs of forms at p in {2, 3, 5, 7}; a point
    (s, t), |s|, |t| <= 60, with both values squares in Q_p must not exist when the solver says 'insoluble', and when
    the solver says 'soluble' such a point is looked for (a miss is only reported)."""
    import random
    rnd = random.Random(seed)
    bad = 0; agree = 0; miss = 0
    for _ in range(ntests):
        p = rnd.choice((2, 3, 5, 7))
        F1 = tuple(rnd.randint(-60, 60) * p**rnd.randint(0, 3) for _ in range(3))
        F2 = tuple(rnd.randint(-60, 60) * p**rnd.randint(0, 3) for _ in range(3))
        if F1 == (0, 0, 0) or F2 == (0, 0, 0): continue
        if (F1[1]**2 - 4*F1[0]*F1[2]) == 0 or (F2[1]**2 - 4*F2[0]*F2[2]) == 0: continue
        r = local_soluble(F1, F2, p, maxdepth=40)
        found = any(padic_square(form_eval(F1, s, t), p) and padic_square(form_eval(F2, s, t), p)
                    for s in range(0, 61) for t in range(-60, 61) if gcd(s, t) == 1)
        if r is False and found: bad += 1
        elif r is True and not found: miss += 1
        else: agree += 1
    return bad, agree, miss

def find_hints(data):
    """Finder only (PARI is not trusted): on each curve left open by the box search, ellrank (effort 1, R supplied);
    a combination of its points and T of class exactly 2, computed by PARI's elladd and checked exactly here."""
    import cypari2
    pari = cypari2.Pari(); pari.allocatemem(2 * 10**9)
    found = []; nopen = 0
    for (m, n) in square_type_ratios(MMAX):
        for nn in (n, -n):
            cd, reasons, surv, pt = second_descent(m, nn)
            if pt is not None or not surv: continue
            nopen += 1
            Mo = cd['M']; S = cd['S']
            E = pari.ellinit([0, Mo['a2'], 0, Mo['a4'], Mo['a6']])
            xR = -4*m*nn*S*(m*m - nn*nn)
            rk = pari.ellrank(E, 1, [[xR, pari.ellordinate(E, xR)[0]]])
            cand = []
            for P in rk[3]:
                X = Fr(str(P[0])); cand.append((P, sqfree(int((X - cd['e']).numerator * (X - cd['e']).denominator))))
            cand.append(([cd['e'], 0], sqfree(cd['b'])))                        # T, alpha(T) = b
            hit = None
            for r_ in range(1, len(cand) + 1):
                for sub in itertools.combinations(range(len(cand)), r_):
                    c = 1
                    for i in sub: c = sqfree(c * cand[i][1])
                    if c == 2:
                        Q = [0]
                        for i in sub: Q = pari.elladd(E, Q, cand[i][0])
                        hit = Q; break
                if hit is not None: break
            if hit is not None and len(hit) == 2:
                X, Y = check_class2(cd, Fr(str(hit[0])), Fr(str(hit[1])))
                found.append([m, nn, str(X), str(Y)])
                say('%4d %5d  hint: X = %s' % (m, nn, X))
    say('find_hints: %d open curves, a point of class 2 for %d (PARI ellrank effort 1, guide; checked exactly)' % (nopen, len(found)))
    if WRITE:
        json.dump({'note': 'points of class 2 on E+- at square-type ratios, found by PARI ellrank (a finder only) by '
                           'square_section.py --find-hints; square_section.py checks each exactly', 'points': found},
                  open(os.path.join(here, 'square_section_hints.json'), 'w'), indent=1)

def pari_control(mmax):
    """Guide only (PARI is not trusted): points of ellrank (effort 1, R supplied) are checked exactly on the curve;
    the class of each (and of each product of two) must not lie in a coset that the descent declares empty."""
    import cypari2
    pari = cypari2.Pari(); pari.allocatemem(8 * 10**8)
    data = json.load(open(os.path.join(here, 'selmer_square_5000.json')))
    nchk = 0; nbad = 0
    for r in data:
        m, nn = r['m'], r['n']
        if m > mmax: continue
        cd = curve_data(m, nn); Mo = cd['M']; S = cd['S']
        E = pari.ellinit([0, Mo['a2'], 0, Mo['a4'], Mo['a6']])
        xR = -4*m*nn*S*(m*m - nn*nn)
        rk = pari.ellrank(E, 1, [[xR, pari.ellordinate(E, xR)[0]]])
        cls = []
        for P in rk[3]:
            X, Y = Fr(str(P[0])), Fr(str(P[1]))
            assert Y*Y == X**3 + Mo['a2']*X*X + Mo['a4']*X + Mo['a6']
            if X != cd['e']: cls.append(sqfree(int((X - cd['e']).numerator * (X - cd['e']).denominator)))
        gens = set(cls)
        prods = set(gens) | set(sqfree(x*y) for x in gens for y in gens)
        K, reps = cosets(m, nn, r['S'])
        for d in reps:
            coset = set(sqfree(d*k) for k in K)
            if coset & prods:
                nchk += 1
                res = descent_class(m, nn, d, want_points=False)
                if not res['survivors']: nbad += 1; say('  CONTRADICTION at %d %d class %d' % (m, nn, d))
    return nchk, nbad

if __name__ == '__main__':
    if NEG:
        negative()
        if WRITE:
            open(os.path.join(here, 'square_section_negative.txt'), 'w').write('\n'.join(out) + '\n')
        sys.exit(0)
    data = {(r['m'], r['n']): r for r in json.load(open(os.path.join(here, 'selmer_square_5000.json')))}
    if '--find-hints' in sys.argv:
        find_hints(data); sys.exit(0)
    hints = load_hints()
    say('square_section.py: second descent on the coverings C_d of E+- at square-type ratios, m <= %d, box %d' % (MMAX, BOUND))
    say('Part 1. The class 2: is C_2(Q) empty?  (status: point = a rational point of class 2 found in the box, or given')
    say('        in square_section_hints.json (found by PARI, a finder only), and checked exactly here;')
    say('        sha = every cover D_lam locally insoluble somewhere, so C_2(Q) is empty; open = otherwise)')
    say('        On every curve: 2 is in Sel(E) (selmer_square_5000.json), R is on the curve with 2R = P and alpha(R) = 2H.')
    rows = []
    tally = {'point': 0, 'sha': 0, 'open': 0}; nhint = 0
    for (m, n) in square_type_ratios(MMAX):
        for nn, lab in ((n, 'E+'), (-n, 'E-')):
            assert 2 in data[(m, nn)]['S'], ('2 not in Sel(E)', m, nn)          # Theorem S1, from the exact descent
            cd, reasons, surv, pt = second_descent(m, nn)
            check_R(m, nn, cd)
            r2 = descent_class(m, nn, 2, want_points=False)
            assert sorted(r2['survivors']) == sorted(surv), 'the two parametrisations disagree'
            nl = len(reasons)
            if pt is None and (m, nn) in hints:
                X, Y = check_class2(cd, *hints[(m, nn)])
                assert surv, ('a point of class 2 where the descent says sha', m, nn)
                nhint += 1
                status = 'point'; detail = 'hint (PARI ellrank, checked exactly here) X=%s' % X
            elif pt is not None:
                lam, s, t, P = pt
                status = 'point'; detail = 'lam=%d (s,t)=(%d,%d) X=%s' % (lam, s, t, P[0])
            elif not surv:
                status = 'sha'; detail = 'all %d covers insoluble (%s)' % (nl, ','.join(sorted(set(reasons.values()))))
            else:
                status = 'open'; detail = 'survivors lam in %s, no point with |s|,|t| <= %d' % (surv, BOUND)
            tally[status] += 1
            rows.append([m, nn, lab, status, len(surv), nl, detail])
            say('%4d %5d %s  %-5s survivors %d/%d  %s' % (m, nn, lab, status, len(surv), nl, detail[:150]))
    say('Part 1 tally: %s (point: %d in the box, %d from hints)' % (tally, tally['point'] - nhint, nhint))
    say('')
    say('Part 2. alpha(E(Q)) = K = {1, HH-, 2H, 2H-} exactly, when every other coset of K in Sel(E) (Sel from')
    say('        selmer_square_5000.json, the exact descent of selmer_square.py) has an empty C_d; then')
    say('        rank = log2 |alpha\'(E\'(Q))| with {1,-1} <= alpha\' <= Sel(E\'); rank = 1 exactly when |Sel(E\')| = 2.')
    r1 = []; nK = 0
    for row in rows:
        m, nn = row[0], row[1]
        rec = data[(m, nn)]
        K, reps, res, allsha = rank_one(m, nn, rec['S'], rec['Sp'], bound=0)
        if allsha:
            nK += 1
            lp = len(rec['Sp']).bit_length() - 1
            if lp == 1:
                r1.append((m, nn)); concl = 'rank E(Q) = 1 exactly'
            else:
                concl = '1 <= rank <= %d' % lp
            say('%4d %5d %s  alpha(E(Q)) = K; %d cosets excluded (reps %s); |Sel(E\')| = %d: %s'
                % (m, nn, row[2], len(reps), reps[:6], len(rec['Sp']), concl))
        row.append('alphaK' if allsha else '-')
    say('Part 2: alpha(E(Q)) = K proved on %d curves; rank exactly 1 on %d: %s' % (nK, len(r1), r1))
    if '--controls' in sys.argv:
        say('')
        bad, agree, miss = brute_control()
        say('Control C1 (local solver vs a box search, random forms, p in 2,3,5,7): contradictions %d, agree %d, '
            'soluble but no point in the box %d' % (bad, agree, miss))
        assert bad == 0
        nchk, nbad = pari_control(min(MMAX, 130))
        say('Control C2 (PARI ellrank points, guide; m <= %d): cosets hit by a PARI point %d, declared empty %d'
            % (min(MMAX, 130), nchk, nbad))
        assert nbad == 0
    def write_tsv(path):
        with open(path, 'w') as f:
            f.write('m\tn_signed\tcurve\tclass2_status\tsurvivors\tcovers\tdetail\talphaE_equals_K\n')
            for r in rows: f.write('\t'.join(map(str, r)) + '\n')
    if '--tsv-to' in sys.argv:                      # for CI: the .tsv to a path named here (under /tmp)
        write_tsv(sys.argv[sys.argv.index('--tsv-to') + 1])
    if WRITE:
        open(os.path.join(here, 'square_section.txt'), 'w').write('\n'.join(out) + '\n')
        write_tsv(os.path.join(here, 'square_section.tsv'))
