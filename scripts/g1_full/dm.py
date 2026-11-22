#!/usr/bin/env python3
"""G1 family, pairs where every elliptic quotient has positive rank: the Demjanenko-Manin argument.

E = E(m, n): y^2 = x (x + 2 S m^2)(x + 2 S n^2),  S = m^2 + n^2,  kappa = 2 S n^2.
A solution t of the five-form system gives P, P' in E(Q) with x(P) = kappa t^2 and x(P) x(P') = kappa^2
(P from F0, F4; P' from F8, F4 after w -> 1/w).  If E(Q) = <G> + torsion (rank exactly 1), write
P = aG + T, P' = a'G + T'.  With |h(x(Q)) - hhat(Q)| <= B for all Q (Silverman),
    |a^2 - a'^2| hhat(G) <= 2B + 2 log kappa,
so either max(|a|, |a'|) <= (2B + 2 log kappa)/hhat(G)  (finitely many P, listed), or |a| = |a'|, i.e.
P' = +-P + T'' with T'' torsion, and x(P) x(P + T'') = kappa^2 is an algebraic equation for x(P).
Every candidate t is then tested in the full system, so a superset of candidates is harmless."""
import sys, json, math, time
from fractions import Fraction as Fr
from sage.all__sagemath_schemes import *

def issq_int(x):
    if x < 0: return False
    r = math.isqrt(x); return r * r == x
def issq(q):
    q = Fr(q); return issq_int(q.numerator) and issq_int(q.denominator)
def sqrt_fr(q): return Fr(math.isqrt(q.numerator), math.isqrt(q.denominator))
def classify_t(m, n, t):
    s = Fr(m*m + n*n, 2)
    F0 = m*m + n*n*t*t; F8 = n*n + m*m*t*t; F4 = s * (1 + t*t); L = F4 - 2*m*n*t; U = F4 + 2*m*n*t
    if not all(issq(v) for v in (F0, F4, F8, L, U)): return 'no'
    cells = [(m*t + n)**2, (m*t - n)**2, (n*t + m)**2, (n*t - m)**2, F0, F8, F4, L, U]
    if len(set(cells)) == 9 and all(c != 0 for c in cells): return 'SOLUTION'
    return 'degenerate'

def tight_bounds(E):
    """Rank bounds from PARI ellrank on the curve and the curves isogenous to it, with effort 0 and 1.
    Returns (lo, hi, points on the minimal model found by PARI)."""
    Ep = E.minimal_model().pari_curve()
    from cypari2 import Pari
    pari = Pari()
    isos = [pari.ellinit(x) for x in pari.ellisomat(Ep, 0, 1)[0]]
    lo, hi = 0, 99; pts = []
    for eff in (0, 1):
        for i, F in enumerate(isos):
            rk = pari.ellrank(F, eff); lo = max(lo, int(rk[0])); hi = min(hi, int(rk[1]))
            if i == 0 and len(rk[3]) > len(pts): pts = rk[3]
            if lo == hi: return lo, hi, pts
    return lo, hi, pts

def dm(m, n):
    S = m*m + n*n; kappa = 2*S*n*n
    E = EllipticCurve(QQ, [0, 2*S*S, 0, 4*S*S*m*m*n*n, 0])
    lo, hi, _pts = tight_bounds(E); rk = E.pari_curve().ellrank(1)
    out = dict(pair=[m, n], rank=[lo, hi])
    if not (lo == 1 and hi == 1): out['status'] = 'rank not 1'; return out
    try:
        gens = E.gens(proof=True); out['gens'] = 'mwrank'
    except Exception as ex:
        pts = rk[3]
        if len(pts) == 1: P0 = E([QQ(pts[0][0]), QQ(pts[0][1])])
        else: P0 = E([kappa, 4*S*S*n*n])          # the point with w = 1, on E for every (m, n)
        if P0.has_finite_order(): out['status'] = 'gens error ' + str(ex)[:60]; return out
        sat, idx, reg = E.saturation([P0], max_prime=-1)
        gens = sat; out['gens'] = 'pari+saturation index %s' % idx
    if len(gens) != 1: out['status'] = 'gens %d' % len(gens); return out
    G = gens[0]; hG = float(G.height())
    B = float(E.silverman_height_bound())
    amax = int(math.floor((2*B + 2*math.log(kappa)) / hG)) + 1
    tors = E.torsion_points()
    out.update(hG=hG, B=B, amax=amax, ntors=len(tors))
    cand = set()          # candidate x-coordinates
    aG = E(0)
    for a in range(0, amax + 1):
        for T in tors:
            P = aG + T
            if not P.is_zero(): cand.add(QQ(P[0]))
        aG = aG + G
    R = PolynomialRing(QQ, 'x'); x = R.gen()
    f = x * (x + 2*S*m*m) * (x + 2*S*n*n); A2 = 2*S*S
    cand.add(QQ(kappa)); cand.add(QQ(-kappa))
    for T in tors:
        if T.is_zero(): continue
        xT, yT = QQ(T[0]), QQ(T[1])
        cand.add(xT)
        # x R(x) with R = f + yT^2 - (x - xT)^2 (A2 + x + xT + kappa^2/x)
        xR = x*(f + yT**2) - (x - xT)**2 * (x*(A2 + x + xT) + kappa**2)
        eq = xR if yT == 0 else (4*yT**2 * f * x**2 - xR**2)
        if eq == 0: out['status'] = 'identity'; return out
        for r, _ in eq.roots(): cand.add(QQ(r))
    verdicts = {}
    for xx in cand:
        w = Fr(int(xx.numerator()), int(xx.denominator())) / kappa
        if w < 0 or not issq(w): continue
        t = sqrt_fr(w)
        for tt in (t, -t):
            verdicts[str(tt)] = classify_t(m, n, tt)
    out['ncand'] = len(cand); out['t'] = {k: v for k, v in verdicts.items() if v != 'no' and len(k) < 40}; out['nt'] = len(verdicts)
    out['status'] = 'ALERT' if 'SOLUTION' in verdicts.values() else 'closed'
    return out

if __name__ == '__main__':
    pairs = [tuple(int(v) for v in s.split(',')) for s in sys.argv[2:] if ',' in s] if sys.argv[1] == '-' else \
            [tuple(int(v) for v in s.split(',')) for s in json.load(open(sys.argv[1]))]
    fo = open(sys.argv[1] + '.dm.jsonl', 'a') if sys.argv[1] != '-' else None
    single = '--single' in sys.argv
    for (m, n) in pairs:
        for (mm, nn) in ([(m, n)] if single else [(m, n), ((m + n, m - n) if (m + n) % 2 else ((m + n)//2, (m - n)//2))]):
            t0 = time.time()
            try: r = dm(mm, nn)
            except Exception as ex: r = dict(pair=[mm, nn], status='error ' + str(ex)[:80])
            r['of'] = [m, n]; r['sec'] = round(time.time() - t0, 1)
            print(json.dumps(r, default=str), flush=True)
            if fo: fo.write(json.dumps(r, default=str) + '\n'); fo.flush()
