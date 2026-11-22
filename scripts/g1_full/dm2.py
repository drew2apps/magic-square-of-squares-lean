#!/usr/bin/env python3
"""Demjanenko-Manin for the second doubled quotient of the G1 curve.

z = t + 1/t.  F0 F8 F4 L = t^4 Q(z) 2S z (2S z - 8mn)  and  F0 F8 F4 U = t^4 Q(z) 2S z (2S z + 8mn),
Q(z) = m^2 n^2 z^2 + (m^2 - n^2)^2, S = m^2 + n^2.  So a solution t gives rational points with
z-coordinates z0 and -z0 on  D: y^2 = z C(z),  C = 2S (2S z - 8mn) Q(z) = c3 z^3 + c2 z^2 + c1 z + c0.
With x = c0 / z:  E: Y^2 = x^3 + c1 x^2 + c0 c2 x + c0^2 c3, and the two points P, P' have x(P') = -x(P).
If E(Q) has rank exactly 1: P = aG + T, P' = a'G + T', |a^2 - a'^2| hhat(G) <= 2B (Silverman bound B),
so max(|a|, |a'|) <= 2B / hhat(G), or |a| = |a'| and x(P) + x(P + T'') = 0 with T'' torsion."""
import sys, json, math, time
from fractions import Fraction as Fr
from sage.all__sagemath_schemes import *
from dm import classify_t, issq, sqrt_fr, tight_bounds

def dm2(m, n):
    S = m*m + n*n; D2 = (m*m - n*n)**2
    c3 = 4*S*S*m*m*n*n; c2 = -16*S*m**3*n**3; c1 = 4*S*S*D2; c0 = -16*S*m*n*D2
    E = EllipticCurve(QQ, [0, c1, 0, c0*c2, c0*c0*c3])
    lo, hi, _pts = tight_bounds(E)
    out = dict(pair=[m, n], curve='D1', rank=[lo, hi])
    if not (lo == 1 and hi == 1): out['status'] = 'rank not 1'; return out
    Em = E.minimal_model(); iso = Em.isomorphism_to(E)
    Ep = Em.pari_curve(); pts = []
    for effort in (0, 1, 2, 3, 4, 5):
        pts = Ep.ellrank(effort)[3]
        if len(pts) >= 1: break
    if len(pts) != 1: out['status'] = 'no point found'; return out
    P0 = Em([QQ(pts[0][0]), QQ(pts[0][1])]); h0 = float(P0.height())
    # rigorous lower bound for hhat on points of infinite order (Cremona-Siksek, as in Sage)
    lam = float(Em.height_function().min(0.01, 20))
    if lam <= 0: out['status'] = 'no height lower bound'; return out
    kmax = int(math.floor(math.sqrt(h0 / lam))) + 1
    sat = Ep.ellsaturation([pts[0]], kmax)
    G1 = Em([QQ(sat[0][0]), QQ(sat[0][1])])
    gens = [iso(G1)]; out['gens'] = 'pari point, saturated at primes <= %d (height lower bound %.4f)' % (kmax, lam)
    if len(gens) != 1: out['status'] = 'gens %d' % len(gens); return out
    G = gens[0]; hG = float(G.height()); B = float(E.silverman_height_bound())
    amax = int(math.floor(1.02 * 2*B / hG)) + 2                              # 2% and +2: margin for rounding
    tors = E.torsion_points()
    out.update(G=str(G[0]), hG=hG, B=B, amax=amax, ntors=len(tors))
    cand = set([QQ(0)])
    aG = E(0)
    for a in range(0, amax + 1):
        for T in tors:
            P = aG + T
            if not P.is_zero(): cand.add(QQ(P[0]))
        aG = aG + G
    R = PolynomialRing(QQ, 'x'); x = R.gen()
    f = x**3 + c1*x**2 + c0*c2*x + c0*c0*c3
    for T in tors:
        if T.is_zero(): continue
        xT, yT = QQ(T[0]), QQ(T[1])
        cand.add(xT)
        Rx = f + yT**2 - (x - xT)**2 * (c1 + xT)
        eq = Rx if yT == 0 else (4*yT**2 * f - Rx**2)
        if eq == 0: out['status'] = 'identity'; return out
        for r, _ in eq.roots(): cand.add(QQ(r))
    verdicts = {}
    ts = set([Fr(0)])
    for xx in cand:
        if xx == 0: continue                      # z = infinity: t = 0 (or infinity)
        z = Fr(c0) / Fr(int(xx.numerator()), int(xx.denominator()))
        Dz = z*z - 4
        if Dz < 0 or not issq(Dz): continue
        r = sqrt_fr(Dz); ts.add((z + r)/2); ts.add((z - r)/2)
    for t in ts: verdicts[str(t)] = classify_t(m, n, t)
    out['ncand'] = len(cand); out['nt'] = len(ts)
    out['t'] = {k: v for k, v in verdicts.items() if v != 'no' and len(k) < 40}
    out['status'] = 'ALERT' if 'SOLUTION' in verdicts.values() else 'closed'
    return out

if __name__ == '__main__':
    for s in sys.argv[1:]:
        m, n = (int(v) for v in s.split(','))
        t0 = time.time()
        try: r = dm2(m, n)
        except Exception as ex: r = dict(pair=[m, n], curve='D1', status='error ' + str(ex)[:80])
        r['sec'] = round(time.time() - t0, 1)
        print(json.dumps(r, default=str), flush=True)
