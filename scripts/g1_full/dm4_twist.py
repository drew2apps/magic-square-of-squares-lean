#!/usr/bin/env python3
"""The height relations of the four-point argument (FOUR_POINTS.md, Proposition 1) on actual rational points.

No rational t0 is known at which all four points Q1..Q4 of C: y^2 = F0 GL are rational (a search for one, also on
quadratic twists, is recorded in twist_search.txt).  But two of the four are rational on a quadratic twist of C for
infinitely many t0:  if GL GU is a square at t0 = p/q, then F0 GL and F0 GU have the same squarefree part d, and Q1
(above t0) and Q2 (above -t0) are rational points of  d y^2 = F0 GL.  The values t0 are the abscissas of the rational
points of the curve  w^2 = GL GU,  which has genus 1 and positive rank for the pairs used; multiples of a point give
t0 of any height.  All the points of the argument that involve Q1 and Q2,
    2 Q1 - R0,   2 Q2 - R0,   Q1 + Q2 - R0,   Q1 - Q2,
are then rational points of the twist  E_d: Y^2 = X (X - 2dA)(X - 2dB)  with X = d x, and the identities (I0), (I+-)
give their abscissas as explicit fractions.  This script
  1. builds these four points from the fractions and checks with the group law of E_d that
     (Q1 + Q2 - R0) + (Q1 - Q2) = +-(2 Q1 - R0)  and  (Q1 + Q2 - R0) - (Q1 - Q2) = +-(2 Q2 - R0)   (as abscissas);
  2. computes their canonical heights (Sage; the canonical height does not change under a twist) and prints
     N(a_i) - H = hhat(2 Q_i - R0)/4 - h(t0)   and   <a1, a2> = (hhat(Q1 + Q2 - R0) - hhat(Q1 - Q2))/4
     next to the bounds of Proposition 1.  Two sets of constants are printed: [-eL, eU] and kappa of the record of the
     pair (they use bounds for points over Q on E_T), and the same with Silverman's bounds for points over the algebraic
     closure, which is what applies to points that are rational only on a twist.
The same is done for Q1 and Q3 (above 1/t0) when F0 F8 is a square at t0, with (Iinv).
usage: dm4_twist.py m,n [max multiple] [max height]     (needs passagemath-schemes)"""
import sys, math, json, time
from sage.all__sagemath_schemes import *

def quartic_points(g, t1, w1, extra, kmax):
    """Abscissas of multiples of the point 'extra' on w^2 = g(t), with base point (t1, w1) (Connell's formulas)."""
    R = g.parent(); T = R.gen(); gs = g(t1 + T); q0 = w1
    a, b, c, d = gs[4], gs[3], gs[2], gs[1]; assert gs[0] == q0*q0 and q0 != 0
    E = EllipticCurve(QQ, [d/q0, c - d*d/(4*q0*q0), 2*q0*b, -4*q0*q0*a, (c - d*d/(4*q0*q0))*(-4*q0*q0*a)])
    def to(t, w):
        Tv = t - t1; return E([(2*q0*(w + q0) + d*Tv)/Tv**2, (4*q0*q0*(w + q0) + 2*q0*(d*Tv + c*Tv*Tv) - d*d*Tv*Tv/(2*q0))/Tv**3])
    out = []
    Em = E.minimal_model(); back = Em.isomorphism_to(E)
    gens = [back(Em([QQ(P[0]), QQ(P[1])])) for P in Em.pari_curve().ellrank(1)[3]]       # points of infinite order (PARI), if any
    tors = E.torsion_points()
    for Gp in gens[:2]:
        for k in range(1, kmax + 1):
            for Tt in tors[:2]:
                Pk = k*Gp + Tt
                if Pk.is_zero() or Pk[1] == 0: continue
                Tv = (2*q0*(Pk[0] + c) - d*d/(2*q0))/Pk[1]
                if Tv != 0: out.append(t1 + Tv)
    return out

def run(m, n, kmax, hmax, rec):
    S = m*m + n*n; A = m**4 - 2*m**3*n + 2*m*m*n*n + 2*m*n**3 + n**4; B = m**4 + 2*m**3*n + 2*m*m*n*n - 2*m*n**3 + n**4
    al, be = 2*A, 2*B
    ET = EllipticCurve(QQ, [0, -(al + be), 0, al*be, 0]); j = ET.j_invariant(); D = ET.discriminant()
    hj = math.log(max(abs(j.numerator()), abs(j.denominator()))); b2 = ET.b2()
    mu = math.log(abs(D))/12 + max(0.0, math.log(abs(j)))/12 + max(0.0, math.log(abs(b2/12)))/2 + (math.log(2)/2 if b2 != 0 else 0)
    U = 2*mu + 2.14; BS = hj/12 + 2*mu + 1.946        # -U <= h(x) - hhat <= BS for points over any number field: Silverman's bound in the form
    # implemented by Sage's silverman_height_bound (documented there as valid over all number fields), with 0.973 for its 0.961
    octic = m**8 + 14*m**4*n**4 + n**8; r0 = m*m*n*n*octic; r1 = (m*m - n*n)**2*octic
    lo4 = math.log(n*n*(m - n)**2) - math.log(2) - 2*math.log(r0) - 2*math.log(r1); up4 = math.log(2) + 2*math.log(2*m**3*n + m**4 - n**4)
    wpm = math.log(32*m*m*n*n*S) + math.log(2*m*m); winv = math.log(4*(m*m - n*n)**2*(m + n)**2) + math.log(8*(m + n)**2)
    eLS = (BS - lo4)/4; eUS = (up4 + U)/4; kpmS = (BS + U + wpm)/4; kinvS = (BS + U + winv)/4
    print('pair (%d,%d): bounds over the algebraic closure: N(a) - H in [%.1f, %.1f], |<a1,a2>| <= %.1f, |<a1,a3>| <= %.1f' % (m, n, -eLS, eUS, kpmS, kinvS)
          + ('; record: [%.1f, %.1f], kappa %.1f' % (-rec['eL'], rec['eU'], rec['kappa']) if rec else ''))
    R = PolynomialRing(QQ, 't'); t = R.gen()
    gLU = 4*S*S*(t*t + 1)**2 - 64*m*m*n*n*t*t; g08 = (m*m + n*n*t*t)*(n*n + m*m*t*t)
    b2c, b4c, b8c = -4*(al + be), 2*al*be, -(al*be)**2
    def hhat(x):
        """An interval that contains the canonical height of a point with abscissa x on E_T (any quadratic twist): exact
        doubling of the abscissa seven times, then  -U <= h(x) - hhat <= BS  for the last one."""
        x = QQ(x)
        for _ in range(7): x = (x**4 - b4c*x*x - b8c) / (4*x**3 + b2c*x*x + 2*b4c*x)
        hx = math.log(int(max(abs(x.numerator()), x.denominator()))); return ((hx - BS)/4**7, (hx + U)/4**7)
    def cubic(x): return x*(x - al)*(x - be)
    rows = []
    for kind, g, base, extra in (('pm', gLU, (0, 2*S), (1, 4*(m*m - n*n))), ('inv', g08, (0, m*n), (1, S))):
        ts = quartic_points(g, QQ(base[0]), QQ(base[1]), (QQ(extra[0]), QQ(extra[1])), kmax)
        seen = set()
        for t0 in ts:
            p, q = abs(t0.numerator()), abs(t0.denominator()); H = math.log(max(p, q))
            if p*q == 0 or p == q or H > hmax or (p, q) in seen: continue
            seen.add((p, q)); T0 = time.time()
            Psi0 = m*m*q*q + n*n*p*p
            def Pi(pp, qq): return 2*m*n*(n*n*pp*pp - m*m*qq*qq) + (m**4 - n**4)*pp*qq
            def Psi(pp, qq): return S*(pp*pp + qq*qq) - 4*m*n*pp*qq
            x1 = QQ(2*Pi(p, q)**2)/(Psi0*Psi(p, q))                             # x(2 Q1 - R0), by (I0)
            if kind == 'pm':
                x2 = QQ(2*Pi(-p, q)**2)/(Psi0*Psi(-p, q))                       # x(2 Q2 - R0), Q2 above -t0
                G = 2*S*(p*p + q*q); rr = ZZ(G*G - 64*m*m*n*n*p*p*q*q).sqrt(); xs = [QQ(32*m*m*n*n*Psi0)/(G + rr), QQ(32*m*m*n*n*Psi0)/(G - rr)]     # (I+-)
            else:
                x2 = QQ(2*Pi(q, p)**2)/((m*m*p*p + n*n*q*q)*Psi(q, p))          # x(2 Q3 - R0), Q3 above 1/t0
                G = S*(p*p + q*q); rr = ZZ(G*G - (m*m - n*n)**2*(p*p - q*q)**2).sqrt(); xs = [QQ(2*(m*m - n*n)**2*Psi(p, q))/(G + rr), QQ(2*(m*m - n*n)**2*Psi(p, q))/(G - rr)]   # (Iinv)
            xa, xb = xs
            # the four abscissas belong to rational points of one quadratic twist: the values of the cubic are in one square class
            assert all((cubic(x1)*cubic(v)).is_square() for v in (x2, xa, xb)), 'not on one twist'
            # group law on the abscissas: x(Pa + Pb) and x(Pa - Pb) are x1 and x2
            den = (xa - xb)**2
            assert x1 + x2 == (2*xa*xb*(xa + xb) + 4*(-(al + be))*xa*xb + 2*al*be*(xa + xb))/den and x1*x2 == (xa*xb - al*be)**2/den, 'group law check fails'
            h1, h2, ha, hb = hhat(x1), hhat(x2), hhat(xa), hhat(xb)
            N1 = (h1[0]/4 - H, h1[1]/4 - H); N2 = (h2[0]/4 - H, h2[1]/4 - H); pr = ((ha[0] - hb[1])/4, (ha[1] - hb[0])/4)
            if pr[0] > -pr[1]: pass
            kb = kpmS if kind == 'pm' else kinvS
            ok = (-eLS <= N1[0] and N1[1] <= eUS) and (-eLS <= N2[0] and N2[1] <= eUS) and max(abs(pr[0]), abs(pr[1])) <= kb
            mid = lambda iv: (iv[0] + iv[1])/2
            rows.append(dict(pair=[m, n], kind=kind, t0=str(t0) if len(str(t0)) < 50 else 'height %.1f' % H, H=round(H, 3), N1_minus_H=round(mid(N1), 3), N2_minus_H=round(mid(N2), 3),
                             pairing=round(mid(pr), 3), error=round((BS + U)/4**8, 4), bound_pairing=round(kb, 2), inside_bounds=ok))
            print('  %-3s H = %8.3f  N(a1) - H = %7.3f  N(a2) - H = %7.3f  pairing = %7.3f  (each +- %.3f)  %s   [%.0f s]' % (kind, H, mid(N1), mid(N2), mid(pr), (BS + U)/4**8, 'inside' if ok else 'OUTSIDE', time.time() - T0), flush=True)
    return rows

if __name__ == '__main__':
    args = [a for a in sys.argv[1:]]; pairs = [tuple(int(v) for v in a.split(',')) for a in args if ',' in a]; nums = [float(a) for a in args if ',' not in a]
    kmax = int(nums[0]) if nums else 6; hmax = nums[1] if len(nums) > 1 else 45.0
    recs = {}
    try:
        import os
        for f in ('dm4_log.jsonl', 'dm4_cross.jsonl'):
            for l in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f)): jj = json.loads(l); recs[tuple(jj['pair'])] = jj
    except Exception: pass
    allrows = []
    for m, n in pairs: allrows += run(m, n, kmax, hmax, recs.get((m, n)))
    bad = [r for r in allrows if not r['inside_bounds']]
    print('%d configurations on twists, %d outside the bounds; largest |N(a) - H| = %.3f, largest |pairing| = %.3f' % (len(allrows), len(bad),
          max([max(abs(r['N1_minus_H']), abs(r['N2_minus_H'])) for r in allrows] + [0]), max([abs(r['pairing']) for r in allrows] + [0])))
