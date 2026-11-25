#!/usr/bin/env python3
"""Demjanenko-Manin on the sigma-quotients of the four-corner curve (no centre form).

For a ratio m:n put t = u/v and
    F0 = m^2 + n^2 t^2,  F8 = n^2 + m^2 t^2,
    GL = (m+n)^2 (t-1)^2 + (m-n)^2 (t+1)^2,  GU = GL(-t).
These are the four corners of the square (4a^2 = F0, 4i^2 = F8, 16c^2 = GL, 16g^2 = GU, scaled).
f_L = F0 F8 GL is palindromic of degree 6, so f_L(t) = t^3 g(z) with z = t + 1/t and g a cubic.
Since z + 2 = (t+1)^2/t and z - 2 = (t-1)^2/t, a t with f_L(t) a square gives a point with z = s := t + 1/t on
    E^+ : Y^2 = (z + 2) g(z)   and on   E^- : Y^2 = (z - 2) g(z).
f_U = F0 F8 GU = f_L(-t) gives g_U(z) = -g(-z), so a t with f_U(t) a square gives a point with z = -s on the
same two curves (E^- from GU and E^+, crossed).  So a t with all four corners square gives two points P, P'
on E^+ (and on E^-) with z(P') = -z(P).

Model: z0 = -2 (E^+) or 2 (E^-), a root of q(z) = (z - z0) g(z); A = q'(z0).  With X = A/(z - z0),
    y^2 = X^3 + b2 X^2 + b4 X + b6,  b2 = q''(z0)/2, b4 = A q'''(z0)/6, b6 = A^2 q''''(z0)/24.
z -> -z is the Moebius map X -> M(X) = -A X / (2 z0 X + A), an involution, and
    H(M(X)) <= (|A| + 2|z0|) H(X),   so  |h(X(P')) - h(X(P))| <= log(|A| + 4).
If E(Q) = <G> + torsion (rank exactly 1), P = aG + T, P' = a'G + T', then with Silverman's B
    |a^2 - a'^2| hhat(G) <= 2B + log(|A| + 4),
so max(|a|, |a'|) <= (2B + log(|A| + 4)) / hhat(G) when |a| != |a'|; when |a| = |a'|, X(P + T'') = M(X(P)) for a
torsion point T'', an algebraic equation in X.  Every candidate z gives the t with t + 1/t = z, and every such t
is tested in the four corner forms (and the five forms of the full square).  A superset of candidates is harmless.
"""
import sys, os, json, math, time
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sage.all__sagemath_schemes import *
from dm import classify_t, issq, sqrt_fr, tight_bounds


def quartic(m, n, s):
    """q(z) = (z + s) g(z) for E^{sign s}, s = +2 or -2, with integer coefficients."""
    R = PolynomialRing(ZZ, 't'); t = R.gen()
    F0 = m*m + n*n*t**2; F8 = n*n + m*m*t**2
    GL = (m + n)**2*(t - 1)**2 + (m - n)**2*(t + 1)**2
    f = F0*F8*GL
    a = [f[k] for k in range(7)]
    assert all(a[k] == a[6 - k] for k in range(7))
    Z = PolynomialRing(ZZ, 'z'); z = Z.gen()
    g = a[3] + a[4]*z + a[5]*(z**2 - 2) + a[6]*(z**3 - 3*z)
    return Z, (z + s)*g


def four_square(m, n, t):
    t = Fr(t)
    F0 = m*m + n*n*t*t; F8 = n*n + m*m*t*t
    GL = (m + n)**2*(t - 1)**2 + (m - n)**2*(t + 1)**2
    GU = (m - n)**2*(t - 1)**2 + (m + n)**2*(t + 1)**2
    return all(issq(v) for v in (F0, F8, GL, GU))


def dm_sigma(m, n, s):
    Z, q = quartic(m, n, s)
    z0 = -s
    assert q(z0) == 0
    d = [q.derivative(k)(z0) for k in range(5)]
    A = d[1]
    assert A != 0, 'double root'
    b2 = d[2] // 2; b4 = A*d[3] // 6; b6 = A*A*d[4] // 24
    assert 2*b2 == d[2] and 6*b4 == A*d[3] and 24*b6 == A*A*d[4]
    E = EllipticCurve(QQ, [0, b2, 0, b4, b6])
    out = dict(pair=[m, n], curve='E%s' % ('+' if s > 0 else '-'), A_digits=len(str(abs(A))))
    lo, hi, _ = tight_bounds(E)
    if lo < hi and EFFORT > 1:                     # more effort on the isogeny class
        from cypari2 import Pari
        pari = Pari(); pari.allocatemem(2*10**9)
        isos = [pari.ellinit(x) for x in pari.ellisomat(E.minimal_model().pari_curve(), 0, 1)[0]]
        for eff in range(2, EFFORT + 1):
            for F in isos:
                rk = pari.ellrank(F, eff); lo = max(lo, int(rk[0])); hi = min(hi, int(rk[1]))
            if lo == hi: break
    mg = MAGMA.get(((m, n), out['curve']))
    if mg is not None and mg < hi:
        assert mg >= lo, (lo, hi, mg)
        hi = mg; out['upper'] = 'Magma TwoPowerIsogenyDescentRankBound: %d' % mg
    out['rank'] = [lo, hi]
    if lo == 0 and hi == 0:                        # rank 0: the points are torsion points
        cand = set([QQ(z0)])
        for T in E.torsion_points():
            if not T.is_zero(): cand.add(None if T[0] == 0 else QQ(z0) + QQ(A)/T[0])
        return finish(m, n, cand, out, 'rank 0')
    if not (lo == 1 and hi == 1):
        out['status'] = 'rank not 1'; return out
    Em = E.minimal_model(); iso = Em.isomorphism_to(E); Ep = Em.pari_curve()
    pts = []
    for effort in (0, 1, 2, 3, 4, 5):
        pts = Ep.ellrank(effort)[3]
        if len(pts) >= 1: break
    if len(pts) != 1:
        out['status'] = 'no point found'; return out
    P0 = Em([QQ(pts[0][0]), QQ(pts[0][1])]); h0 = float(P0.height())
    lam = float(Em.height_function().min(0.01, 20))         # lower bound for hhat on points of infinite order
    if lam <= 0:
        out['status'] = 'no height lower bound'; return out
    kmax = int(math.floor(math.sqrt(h0 / lam))) + 1
    sat = Ep.ellsaturation([pts[0]], kmax)
    G = iso(Em([QQ(sat[0][0]), QQ(sat[0][1])]))
    hG = float(G.height()); B = float(E.silverman_height_bound())
    # exact enclosure of hhat(G): |h(x(2^k G)) - hhat(2^k G)| <= B, hhat(2^k G) = 4^k hhat(G)
    def hx(P):
        X = P[0]; return math.log(max(abs(int(X.numerator())), abs(int(X.denominator()))))
    k = 7; Q = G
    for _ in range(k): Q = 2*Q
    hQ = hx(Q); hG_lo = (hQ - B)/4**k; hG_hi = (hQ + B)/4**k
    assert hG_lo <= hG <= hG_hi, (hG_lo, hG, hG_hi)
    j = E.j_invariant(); hj = math.log(max(abs(int(j.numerator())), abs(int(j.denominator()))))
    assert hj >= 2.62          # so the upper side hhat - h <= B - h(j)/12 + 0.218 is below B (see G1_FULL_SYSTEM.md)
    logK = math.log(abs(A) + 4)
    amax = int(math.floor((2*B + logK)/hG_lo)) + 1            # from the exact lower end of the enclosure
    tors = E.torsion_points()
    out.update(G=str(G.xy()[0])[:60], hG=hG, hG_enclosure=[hG_lo, hG_hi], B=B, hj=hj, logK=logK, amax=amax,
               ntors=len(tors), sat='<= %d' % kmax, lam=lam)
    cand = set()           # candidate z-coordinates (None = infinity)
    cand.add(QQ(z0))       # P = O
    def z_of_X(X):
        return None if X == 0 else QQ(z0) + QQ(A)/X
    # P' = O: 2 z0 X + A = 0
    cand.add(z_of_X(QQ(-A)/(2*z0)))
    aG = E(0); worst = 0.0
    for a in range(0, amax + 1):
        for T in tors:
            P = aG + T
            if not P.is_zero():
                cand.add(z_of_X(P[0]))
                worst = max(worst, abs(hx(P) - a*a*hG))     # a check of B on the listed points
        aG = aG + G
    assert worst <= B
    out['max |h - hhat| on the list'] = worst
    # |a| = |a'|: X(P + T'') = M(X(P))
    R = PolynomialRing(QQ, 'x'); x = R.gen()
    f = x**3 + b2*x**2 + b4*x + b6
    D = 2*z0*x + A
    for X in (0, QQ(-A)/z0):                # T'' = O: X = M(X)
        cand.add(z_of_X(QQ(X)))
    for T in tors:
        if T.is_zero(): continue
        xT, yT = QQ(T[0]), QQ(T[1])
        cand.add(z_of_X(xT))
        ND = (-A*x + (b2 + x + xT)*D)*(x - xT)**2
        lhs = (f + yT**2)*D - ND
        eq = lhs if yT == 0 else (4*yT**2*f*D**2 - lhs**2)
        if eq == 0:
            out['status'] = 'identity'; return out
        for r, _ in eq.roots(): cand.add(z_of_X(r))
    return finish(m, n, cand, out, 'DM')


def finish(m, n, cand, out, how):
    ts = set([Fr(0)])                       # z = infinity: t = 0 or infinity (t = 0 stands for both)
    for zz in cand:
        if zz is None: continue
        zf = Fr(int(zz.numerator()), int(zz.denominator()))
        Dz = zf*zf - 4
        if Dz < 0 or not issq(Dz): continue
        r = sqrt_fr(Dz); ts.add((zf + r)/2); ts.add((zf - r)/2)
    four = {}
    for t in ts:
        for tt in ([t, Fr(1, 1)/t] if t != 0 else [t]):
            if four_square(m, n, tt): four[str(tt)] = classify_t(m, n, tt)
    # t = infinity: F0 ~ n^2 t^2, F8 ~ m^2 t^2, GL, GU ~ 2(m^2+n^2) t^2: corners square iff 2(m^2+n^2) is a square
    if issq(2*(m*m + n*n)): four['inf'] = 'degenerate'
    out['ncand'] = len(cand); out['nt'] = len(ts)
    out['four'] = four
    out['status'] = 'ALERT' if 'SOLUTION' in four.values() else 'closed'
    out['how'] = how
    return out


EFFORT = 1
MAGMA = {}
def read_magma(path):
    """Upper rank bounds from Magma's calculator (lines 'BOUNDS m n curve r' of a transcript)."""
    import re
    for m, n, c, r in re.findall(r'<line>BOUNDS (\d+) (\d+) (E[+-]) (\d+)</line>', open(path).read()):
        MAGMA[((int(m), int(n)), c)] = int(r)

if __name__ == '__main__':
    for arg in sys.argv[1:]:
        if arg.startswith('--effort='): EFFORT = int(arg[9:]); continue
        if arg.startswith('--magma='):
            for path in arg[8:].split(','): read_magma(path)
            continue
        m, n = (int(v) for v in arg.split(','))
        for s in (2, -2):
            t0 = time.time()
            try: r = dm_sigma(m, n, s)
            except Exception as ex: r = dict(pair=[m, n], curve='E%s' % ('+' if s > 0 else '-'), status='error ' + repr(ex)[:120])
            r['sec'] = round(time.time() - t0, 1)
            print(json.dumps(r, default=str), flush=True)
