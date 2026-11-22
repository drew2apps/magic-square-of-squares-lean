#!/usr/bin/env python3
"""Demjanenko-Manin for the quotient y^2 = F0(t) GL(t), which occurs four times in the Jacobian.

C: y^2 = q(t) = (n^2 t^2 + m^2)(2S t^2 - 8mn t + 2S), S = m^2 + n^2, has the rational point
t1 = m/n, y1 = 2m(m^2 - n^2)/n.  A solution t0 of the five-form system gives rational points on C with
t-coordinates t0 (from F0 L) and -t0 (from F0 U).  Here t is a degree-2 function that is not even, so
its height is compared with canonical heights through the fibre {Q, R0 - Q} of t:
with tau = t - t1, sigma = 1/tau and X the Weierstrass coordinate (Connell's formulas),
    X(Q) + X(R0 - Q) = 4 q0^2 sigma^2 + 2 d sigma,
    X(Q) X(R0 - Q)   = (d^2 - 4 q0^2 c) sigma^2 - 4 q0^2 b sigma - 4 q0^2 a,
so [1 : sum : product] is the image of sigma under a degree-2 map P^1 -> P^2, and
    | 2 h(sigma) - h(X(Q)) - h(X(R0 - Q)) | <= C1 + log 2,
C1 = log of the larger row-sum norm of the integer matrix of the map and of its adjugate.
With Silverman's bound B:  | 2 h(sigma(Q)) - hhat(Q) - hhat(R0 - Q) | <= K = C1 + log 2 + 2B,
and | h(t) - h(sigma) | <= h(t1) + log 2.
If E(Q) = <G> + torsion (rank 1), Q = aG + T, R0 = rho G + T0:  hhat(Q) + hhat(R0 - Q) = hG (a^2 + (rho - a)^2).
For Q1 = a1 G + T1, Q2 = a2 G + T2 with t(Q2) = -t(Q1):  2 |a1 - a2| |a1 + a2 - rho| hG <= 2K + 4(h(t1) + log 2).
So either |a1| <= W + |rho|/2, with W the right side over 2 hG: these points are listed, for every torsion point T;
or a1 = a2 or a1 + a2 = rho: then Q2 = Q + T'' with T'' torsion and Q = Q1 or R0 - Q1 (which has the same t),
so t(Q + T'') = -t(Q). For T'' = O this forces t = 0; for the other torsion points it is a polynomial equation
in t, solved at the end. The output records C1, K, rho, W and the largest observed value of
|2 h(sigma) - hhat(Q) - hhat(R0 - Q)|, which must be below K."""
import sys, json, math, time
from fractions import Fraction as Fr
from sage.all__sagemath_schemes import *
from dm import classify_t, tight_bounds, issq

def hQ(x):
    x = QQ(x); return math.log(max(abs(int(x.numerator())), abs(int(x.denominator())), 1))

def dm3(m, n):
    S = m*m + n*n
    R = PolynomialRing(QQ, 'T'); T = R.gen()
    q = (n*n*T**2 + m*m) * (2*S*T**2 - 8*m*n*T + 2*S)
    t1 = QQ(m)/n; q0 = QQ(2*m*(m*m - n*n))/n
    assert q(t1) == q0**2
    qs = q(t1 + T); e0, d, c, b, a = [qs[i] for i in range(5)]
    assert e0 == q0**2
    a1 = d/q0; a2 = c - d**2/(4*q0**2); a3 = 2*q0*b; a4 = -4*q0**2*a; a6 = a2*a4
    E0 = EllipticCurve(QQ, [a1, a2, a3, a4, a6])
    Em = E0.minimal_model(); phi = E0.isomorphism_to(Em); psi = Em.isomorphism_to(E0)
    lo, hi, _ = tight_bounds(Em)
    out = dict(pair=[m, n], curve='E0L', rank=[lo, hi], model=[str(v) for v in Em.ainvs()])
    if not (lo == 1 and hi == 1): out['status'] = 'rank not 1'; return out
    tors = Em.torsion_points(); out['ntors'] = len(tors)
    Ep = Em.pari_curve(); pts = []
    for effort in (0, 1, 2, 3):
        pts = Ep.ellrank(effort)[3]
        if len(pts) >= 1: break
    if len(pts) != 1: out['status'] = 'no point found'; return out
    P0 = Em([QQ(pts[0][0]), QQ(pts[0][1])]); h0 = float(P0.height())
    lam = float(Em.height_function().min(0.01, 20))
    if lam <= 0: out['status'] = 'no height lower bound'; return out
    kmax = int(math.floor(math.sqrt(h0 / lam))) + 1
    sat = Ep.ellsaturation([pts[0]], kmax)
    G = Em([QQ(sat[0][0]), QQ(sat[0][1])]); hG = float(G.height())
    out['gens'] = 'pari point, saturated at primes <= %d (height lower bound %.4f)' % (kmax, lam)
    out['G'] = [str(G[0]), str(G[1])]; out['hG'] = hG
    B = float(Em.silverman_height_bound())
    R0 = phi(E0([-a2, a1*a2 - a3]))
    # R0 = rho G + torsion
    rho = None
    for k in range(0, int(math.sqrt(float(R0.height()) / hG)) + 3):
        for sgn in (1, -1):
            if (R0 - sgn*k*G).has_finite_order(): rho = sgn*k
        if rho is not None: break
    if rho is None: out['status'] = 'R0 not a multiple of G'; return out
    # X' = (X - r)/u^2 : find u^2 and r from two points
    A1 = 4*q0**2; B1 = 2*d; A2 = d**2 - 4*q0**2*c; B2 = -4*q0**2*b; C2 = -4*q0**2*a
    Pa = psi(G); Pb = psi(2*G)
    u2 = (Pa[0] - Pb[0]) / (G[0] - (2*G)[0]); r = Pa[0] - u2*G[0]
    assert psi(3*G)[0] == u2*(3*G)[0] + r
    rows = [[QQ(0), QQ(0), QQ(1)],
            [A1/u2, B1/u2, -2*r/u2],
            [(A2 - r*A1)/u2**2, (B2 - r*B1)/u2**2, (C2 + r*r)/u2**2]]
    L = lcm([QQ(v).denominator() for row in rows for v in row])
    M = matrix(ZZ, [[ZZ(v*L) for v in row] for row in rows])
    if M.det() == 0: out['status'] = 'degenerate map'; return out
    Madj = M.adjugate()
    norm = lambda MM: max(sum(abs(int(v)) for v in row) for row in MM.rows())
    C1 = math.log(max(norm(M), norm(Madj)))
    K = C1 + math.log(2) + 2*B
    K2 = hQ(t1) + math.log(2)
    Delta = 2*K + 4*K2
    W = Delta / (2*hG)
    amax = int(math.floor(1.02*(W + abs(rho)))) + 3
    out.update(B=B, C1=C1, K=K, rho=rho, W=W, amax=amax)
    def sigma_of(Q):                       # Q on Em, not O; returns sigma = 1/(t - t1), or None for sigma = infinity
        Qs = R0 - Q
        if Q.is_zero() or Qs.is_zero(): return None
        x1 = psi(Q)[0]; x2 = psi(Qs)[0]; S1 = x1 + x2; S2 = x1*x2
        den = A2*B1 - A1*B2
        assert den != 0
        sg = (A2*S1 - A1*S2 + A1*C2) / den
        assert A1*sg**2 + B1*sg == S1 and A2*sg**2 + B2*sg + C2 == S2
        return sg
    verdicts = {'0': classify_t(m, n, Fr(0))}; worst = 0.0; npts = 0
    def test_t(tt, on_curve=False):
        tf = Fr(int(QQ(tt).numerator()), int(QQ(tt).denominator()))
        if on_curve: assert issq((n*n*tf*tf + m*m) * (2*S*tf*tf - 8*m*n*tf + 2*S))   # the point is on C
        verdicts[str(tf)] = classify_t(m, n, tf)
    Q = Em(0)
    for k in range(0, amax + 1):
        for Qk0 in ([Q] if k == 0 else [Q, -Q]):
            for Tt in tors:
                Qk = Qk0 + Tt
                sg = sigma_of(Qk)
                if sg is None: tt = t1                     # tau = 0
                elif sg == 0: continue                      # t = infinity
                else:
                    tt = t1 + 1/sg
                    N = float(Qk.height()) + float((R0 - Qk).height())
                    worst = max(worst, abs(2*hQ(sg) - N))    # numerical check of the height inequality
                npts += 1; test_t(tt, True)
        Q = Q + G
    # the cases a1 = a2 and a1 + a2 = rho:  t(Q + T'') = -t(Q) for a torsion point T''
    Kt = FractionField(PolynomialRing(QQ, 't')); t = Kt.gen()
    class F2:                                   # alpha + beta*y in Q(t)[y]/(y^2 - q(t))
        def __init__(s, al, be=0): s.al = Kt(al); s.be = Kt(be)
        def __add__(s, o): o = o if isinstance(o, F2) else F2(o); return F2(s.al + o.al, s.be + o.be)
        __radd__ = __add__
        def __neg__(s): return F2(-s.al, -s.be)
        def __sub__(s, o): return s + (-(o if isinstance(o, F2) else F2(o)))
        def __rsub__(s, o): return (-s) + o
        def __mul__(s, o):
            o = o if isinstance(o, F2) else F2(o)
            return F2(s.al*o.al + s.be*o.be*qt, s.al*o.be + s.be*o.al)
        __rmul__ = __mul__
        def inv(s):
            nn = s.al**2 - s.be**2*qt; return F2(s.al/nn, -s.be/nn)
        def __truediv__(s, o): o = o if isinstance(o, F2) else F2(o); return s * o.inv()
    qt = q(t); tau = t - t1; y = F2(0, 1)
    Xg = (2*q0*(y + q0) + d*tau) * F2(1/tau**2)
    Yg = (4*q0**2*(y + q0) + 2*q0*(d*tau + c*tau**2) - d**2*tau**2/(2*q0)) * F2(1/tau**3)
    nalg = 0
    for Tt in tors:
        if Tt.is_zero(): continue
        TT = psi(Tt); xT, yT = QQ(TT[0]), QQ(TT[1])
        lamb = (Yg - yT) / (Xg - xT); nu = yT - lamb*xT
        X3 = lamb*lamb + a1*lamb - a2 - Xg - xT
        Y3 = -(lamb + a1)*X3 - nu - a3
        taup = (2*q0*(X3 + c) - d**2/(2*q0)) / Y3
        ee = taup + (t1 + t)                    # t(Q + T) + t(Q)
        Nn = ee.al**2 - ee.be**2*qt
        num = Nn.numerator()
        if num == 0: out['status'] = 'identity'; return out
        for rr, _ in num.roots(): test_t(rr); nalg += 1
        test_t(t1 + (2*q0*(xT + c) - d**2/(2*q0))/yT if yT != 0 else t1)
    out['nalg'] = nalg
    out['npts'] = npts; out['worst_observed'] = worst
    if worst > K: out['status'] = 'height inequality violated'; return out
    out['t'] = {k: v for k, v in verdicts.items() if v != 'no' and len(k) < 40}
    out['status'] = 'ALERT' if 'SOLUTION' in verdicts.values() else 'closed'
    return out

if __name__ == '__main__':
    for s in sys.argv[1:]:
        m, n = (int(v) for v in s.split(','))
        t0 = time.time()
        try: r = dm3(m, n)
        except Exception as ex: r = dict(pair=[m, n], curve='E0L', status='error ' + repr(ex)[:120])
        r['sec'] = round(time.time() - t0, 1)
        print(json.dumps(r, default=str), flush=True)
