"""Demjanenko-Manin on rank-1 quotients of the G1 system, ported to PARI (no Sage). EXPERIMENT,
written by an AI model (Claude), not read by a person.  Curves (as in the repo's dm.py, dm2.py, dm_sigma.py):
  E   : y^2 = x(x + 2 S m^2)(x + 2 S n^2), relation x(P) x(P') = kappa^2, kappa = 2 S n^2, t^2 = x / kappa
  D1  : dm2.py's model, relation x(P') = -x(P), z = c0 / x, t + 1/t = z
  Sig+/Sig- : dm_sigma.py's models, relation X(P') = M(X(P)), z = z0 + A / X, t + 1/t = z
Generic step: rank exactly 1 (PARI ellrank on the isogeny class), generator G from PARI, saturated by
ellsaturation at primes <= SATB (NOT a proof that G generates: label 'conditional'), height bounds
BUP, U of Cremona-Prickett-Siksek / Silverman as in the repo's hbounds.py, then every aG + T with
|a| <= amax, and the roots of the polynomial for P' = +-P + T''.  Every candidate t is tested in the
five forms (degenerate / SOLUTION).  usage: dm_pari.py pairs.json out.jsonl"""
import sys, json, math, time
from fractions import Fraction as Fr
import cypari2
pari = cypari2.Pari(); pari.allocatemem(10**9)
SATB = 1000
x = pari('x')

def issq_int(v): return v >= 0 and math.isqrt(v) ** 2 == v
def issq(q): q = Fr(q); return issq_int(q.numerator) and issq_int(q.denominator)
def sqrt_fr(q): return Fr(math.isqrt(q.numerator), math.isqrt(q.denominator))
def fr(z): z = pari(z); return Fr(int(pari.numerator(z)), int(pari.denominator(z)))
def classify_t(m, n, t):
    s = Fr(m*m + n*n, 2)
    F0 = m*m + n*n*t*t; F8 = n*n + m*m*t*t; F4 = s*(1 + t*t); L = F4 - 2*m*n*t; U = F4 + 2*m*n*t
    if not all(issq(v) for v in (F0, F4, F8, L, U)): return 'no'
    cells = [(m*t + n)**2, (m*t - n)**2, (n*t + m)**2, (n*t - m)**2, F0, F8, F4, L, U]
    return 'SOLUTION' if len(set(cells)) == 9 and all(c != 0 for c in cells) else 'degenerate'

def alpha(code, c):
    if c == 1: return Fr(0)
    if code >= 5:
        mm = code - 4
        if mm % 2 == 0 and c in (2, mm): return Fr(mm, 4)
        if mm % 2 == 1 and c == mm: return Fr(mm*mm - 1, 4*mm)
    if code == 3 and c == 2: return Fr(1, 2)
    if code == 4 and c == 3: return Fr(2, 3)
    if code == -1 and c in (2, 4): return Fr(1)
    if code <= -5:
        mm = -code - 4
        if c == 2: return Fr(1)
        if c == 4: return Fr(mm + 4, 4)
    if code == -4 and c == 3: return Fr(4, 3)
    if code == -3 and c == 2: return Fr(3, 2)
    raise ValueError((code, c))
def lp(q): q = abs(Fr(q)); return math.log(q.numerator) - math.log(q.denominator) if q > 1 else 0.0
def height_bounds(E):
    """BUP >= h(x) - hhat and U >= hhat - h(x) for points of this integral model (hbounds.py, in floats + margin)."""
    D = int(E[11]); b2 = int(E[5]); j = fr(E[12])
    na = 0.0; Dmin = int(pari.ellminimalmodel(E)[11])
    for p, e in zip(*[list(map(int, c)) for c in pari.factor(abs(D))]):
        lr = pari.elllocalred(E, p)
        emin = int(pari.valuation(Dmin, p))
        na += float(alpha(int(lr[1]), int(lr[3])) + Fr(e - emin, 6)) * math.log(p)
    tail = lp(Fr(abs(b2), 12)) + (math.log(2) if b2 else 0)
    BUP = na + lp(j) / 4 + tail + 1.946
    U = math.log(abs(D)) / 6 + lp(j) / 6 + tail + 2.14
    return BUP * 1.0001 + 1e-6, U * 1.0001 + 1e-6

def rank_bounds(E):
    isos = [pari.ellinit(c) for c in pari.ellisomat(pari.ellminimalmodel(E), 0, 1)[0]]
    lo, hi = 0, 99
    for eff in (0, 1):
        for F in isos:
            rk = pari.ellrank(F, eff); lo = max(lo, int(rk[0])); hi = min(hi, int(rk[1]))
            if lo == hi: return lo, hi
    return lo, hi

def dm_curve(m, n, ainv, rel, Cextra, to_t, name):
    """rel: (alpha, beta) PARI polynomials in x with relation alpha(x) + beta(x) x' = 0."""
    E = pari.ellinit(ainv)
    out = dict(pair=[m, n], curve=name)
    lo, hi = rank_bounds(E); out['rank'] = [lo, hi]
    if not (lo == 1 and hi == 1): out['status'] = 'rank not 1'; return out
    pts = pari.ellrank(E, 2)[3]
    if len(pts) != 1:
        # isogenous curve found the point; search harder here
        pts = pari.ellrank(E, 4)[3]
    if len(pts) != 1: out['status'] = 'no generator found'; return out
    G = pari.ellsaturation(E, pts, SATB)[0]
    hG = float(pari.ellheight(E, G))
    BUP, U = height_bounds(E); B = max(BUP, U)
    amax = int(math.floor(1.02 * (2*B + Cextra) / hG)) + 2
    out.update(hG=hG, B=B, amax=amax, gens='pari, saturated at primes <= %d (conditional)' % SATB)
    T = pari.elltors(E); tors = [pari('[0]')]
    for g, d in zip(list(T[2]), [int(v) for v in T[1]]):
        tors = [pari.elladd(E, P, pari.ellmul(E, g, k)) for P in tors for k in range(d)]
    out['ntors'] = len(tors)
    cand = set()
    for a in range(0, amax + 1):
        aG = pari.ellmul(E, G, a)
        for Tp in tors:
            P = pari.elladd(E, aG, Tp)
            if len(P) == 2: cand.add(fr(P[0]))
    al, be = rel
    a2 = ainv[1]; f = x**3 + ainv[1]*x**2 + ainv[3]*x + ainv[4]
    for Tp in tors:
        if len(Tp) == 1:
            eq = al + be * x
        else:
            xT, yT = pari(Tp[0]), pari(Tp[1]); cand.add(fr(xT))
            A0 = f + yT**2 - (a2 + x + xT) * (x - xT)**2
            lhs = al * (x - xT)**2 + be * A0
            eq = lhs if yT == 0 else lhs**2 - 4 * yT**2 * be**2 * f
        if eq == 0: out['status'] = 'identity'; return out
        fa = pari.factor(eq)
        for fac in fa[0]:
            if pari.poldegree(fac) == 1:
                cand.add(fr(-pari.polcoef(fac, 0) / pari.polcoef(fac, 1)))
    verdicts = {}
    for X in cand:
        for t in to_t(X):
            verdicts[str(t)] = classify_t(m, n, t)
    out['ncand'] = len(cand); out['nt'] = len(verdicts)
    out['t'] = {k: v for k, v in verdicts.items() if v != 'no' and len(k) < 40}
    out['status'] = 'ALERT' if 'SOLUTION' in verdicts.values() else 'closed'
    return out

def t_from_z(z):
    # t + 1/t = z  ->  t^2 - z t + 1 = 0
    disc = z*z - 4
    if disc < 0 or not issq(disc): return []
    r = sqrt_fr(disc); return list({(z + r) / 2, (z - r) / 2} - {Fr(0)})

def curve_E(m, n):
    S = m*m + n*n; kappa = 2*S*n*n; res = []
    # E
    def tE(X):
        if X == 0: return [Fr(0)]
        w = X / kappa
        if w < 0 or not issq(w): return []
        t = sqrt_fr(w); return [t, -t]
    res.append(([0, 2*S*S, 0, 4*S*S*m*m*n*n, 0], (-kappa**2 + 0*x, x), 2*math.log(kappa), tE, 'E'))
    return res
def all_curves(m, n):
    S = m*m + n*n; res = curve_E(m, n)
    # D1 (dm2.py)
    D2 = (m*m - n*n)**2
    c3 = 4*S*S*m*m*n*n; c2 = -16*S*m**3*n**3; c1 = 4*S*S*D2; c0 = -16*S*m*n*D2
    def tD(X):
        if X == 0: return []
        return t_from_z(Fr(c0) / X)
    res.append(([0, c1, 0, c0*c2, c0*c0*c3], (x, 1 + 0*x), 0.0, tD, 'D1'))
    # sigma curves (dm_sigma.py)
    t = pari('t'); z = pari('z')
    F0 = m*m + n*n*t**2; F8 = n*n + m*m*t**2; GL = (m+n)**2*(t-1)**2 + (m-n)**2*(t+1)**2
    fL = F0 * F8 * GL
    a = [int(pari.polcoef(fL, k, 't')) for k in range(7)]
    assert all(a[k] == a[6-k] for k in range(7))
    gz = a[3] + a[4]*z + a[5]*(z**2 - 2) + a[6]*(z**3 - 3*z)
    for z0, nm in ((-2, 'Sig+'), (2, 'Sig-')):
        q = (z - z0) * gz
        A = int(pari.subst(pari.deriv(q, 'z'), 'z', z0))
        b2 = fr(pari.subst(pari.deriv(pari.deriv(q, 'z'), 'z'), 'z', z0)) / 2
        b4 = A * fr(pari.subst(pari.deriv(pari.deriv(pari.deriv(q, 'z'), 'z'), 'z'), 'z', z0)) / 6
        b6 = A * A * Fr(int(pari.polcoef(q, 4, 'z')) * 24) / 24
        assert b2.denominator == b4.denominator == b6.denominator == 1
        def tS(X, z0=z0, A=A):
            if X == 0: return []
            return t_from_z(z0 + Fr(A) / X)
        # X' = -A X / (2 z0 X + A):  A x + (2 z0 x + A) x' = 0
        res.append(([0, int(b2), 0, int(b4), int(b6)], (A * x, 2*z0*x + A), math.log(abs(A) + 4), tS, nm))
    return res

if __name__ == '__main__':
    pairs = json.load(open(sys.argv[1]))['open']; fo = open(sys.argv[2], 'a')
    for (m, n) in pairs:
        mm, nn = (m + n, m - n) if (m + n) % 2 else ((m + n)//2, (m - n)//2)
        jobs = [(m, n, c) for c in all_curves(m, n)] + [(mm, nn, (c[0], c[1], c[2], c[3], c[4] + '-partner')) for c in all_curves(mm, nn)]
        for (pm, pn, (ainv, rel, C, tt, nm)) in jobs:
            t0 = time.time()
            try: r = dm_curve(pm, pn, ainv, rel, C, tt, nm)
            except Exception as ex: r = dict(pair=[m, n], curve=nm, status='error ' + str(ex)[:100])
            r['sec'] = round(time.time() - t0, 1); r['of'] = [m, n]
            print(json.dumps(r, default=str), flush=True); fo.write(json.dumps(r, default=str) + '\n'); fo.flush()
