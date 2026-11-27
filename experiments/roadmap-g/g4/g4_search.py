"""G4: a search for a solution on the square-type ratios (EXPERIMENT, written by an AI model; not read by a person;
nothing here changes the census).
A solution t of the five-form system (F0, F8, F4, L, U all squares; global.py's classify_t) gives a rational point on
  * E+ and E-: Y^2 = (z +- 2) g(z), z = t + 1/t (UNIFORM_RANK.md section 1; sigma_exact.model(m, n, +-2)), because
    F0 F8 L = t^3 g(z)/4 and z +- 2 = (t +- 1)^2 / t;
  * each quotient y^2 = P(w), w = t^2, of global.py's curves().
For each ratio: PARI's ellrank on each curve; on every curve whose rank bounds meet (lo = hi <= RMAX) the generators
are saturated (ellsaturation at primes <= 100) and, with the torsion, every point of canonical height <= H of the
group they generate is listed (Fincke-Pohst on the height matrix; H is the run-time bound).  Each point is mapped back
to t (z^2 - 4 a square, resp. w a square), and t is tested exactly: the number of the five forms that are squares, and
classify_t.  Reported: every t with at least four of the five forms squares, and any SOLUTION.
Completeness is not claimed: the generators are PARI's, the rank is PARI's (not proved here), saturation only at
primes <= 100, and only points of height <= H are listed.  The test is that the known degenerate t = +-1 (all five
forms squares on a square-type ratio) is found on every curve that carries it (planted test).
usage: g4_search.py LIST i N [H] [--write]"""
import sys, os, json, math
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get('MSS_REPO', os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(REPO, 'scripts', 'g1_full'))
import importlib
G = importlib.import_module('global'); SE = importlib.import_module('sigma_exact')
pari = G.pari
issq, sqrt_fr = G.issq, G.sqrt_fr
RMAX = 3

def forms(m, n, t):
    s = Fr(m*m + n*n, 2)
    F0 = m*m + n*n*t*t; F8 = n*n + m*m*t*t; F4 = s*(1 + t*t); L = F4 - 2*m*n*t; U = F4 + 2*m*n*t
    return [F0, F8, F4, L, U]

def fr(x): return Fr(int(pari.numerator(x)), int(pari.denominator(x)))

def curve_list(m, n):
    """(name, E, map point -> list of t)."""
    out = []
    for s, name in ((2, 'E+'), (-2, 'E-')):
        M = SE.model(m, n, s); E = pari.ellinit([0, M['a2'], 0, M['a4'], M['a6']])
        def to_t(P, M=M):
            z = Fr(M['z0']) if len(P) == 1 else (None if P[0] == 0 else M['z0'] + Fr(M['A']) / fr(P[0]))
            if z is None: return ['inf', Fr(0)]                       # z = infinity: t = 0 or infinity
            d = z*z - 4
            if not issq(d): return []
            r = sqrt_fr(d); return list({(z + r)/2, (z - r)/2})
        out.append((name, E, to_t, [Fr(M['A'], z - M['z0']) for z in (2, -2) if z != M['z0']]))
    for name, poly in G.curves(m, n).items():
        if len(poly) == 5: c = [poly[4], poly[3], poly[2], poly[1]]; inv = True
        else: c = poly; inv = False
        c0, c1, c2, c3 = c
        E = pari.ellinit([0, c2, 0, c1*c3, c0*c3*c3])
        def to_t(P, c3=c3, inv=inv):
            if len(P) == 1: w = Fr(0) if inv else None
            else:
                x = fr(P[0]) / c3
                w = (None if x == 0 else 1/x) if inv else x
            if w is None: return ['inf']
            if w < 0 or not issq(w): return []
            r = sqrt_fr(w); return list({r, -r})
        out.append((name, E, to_t, [Fr(c3)]))                         # w = 1, t = +-1
    return out

def points(E, gens, H):
    """All points T + sum a_i g_i with canonical height <= H (exact list of coefficient vectors by Fincke-Pohst)."""
    tors = pari.elltors(E); T = [pari('[0]')]
    for g, d in zip(list(tors[2]), [int(x) for x in tors[1]]):
        T = [pari.elladd(E, P, pari.ellmul(E, g, k)) for P in T for k in range(d)]
    r = len(gens); out = []
    if r == 0: return T, [[]]
    Hm = pari.ellheightmatrix(E, gens)
    vecs = [[0]*r]
    fp = pari.qfminim(Hm, H + 1e-9, 10**5, 2)          # all nonzero v with v^t Hm v <= H, up to sign
    Mv = fp[2]; ncol = int(pari.matsize(Mv)[1]) if r else 0
    for k in range(ncol):
        v = [int(Mv[i, k]) for i in range(r)]; vecs.append(v); vecs.append([-x for x in v])
    for v in vecs:
        Q = pari('[0]')
        for a, g in zip(v, gens): Q = pari.elladd(E, Q, pari.ellmul(E, g, a))
        for t0 in T: out.append(pari.elladd(E, Q, t0))
    return out, vecs

def run(m, n, K):
    rec = dict(ratio=[m, n], curves={}, hits=[], solutions=[])
    for name, E, to_t, Xs in curve_list(m, n):
        kp = []
        for X in Xs:
            y2 = X**3 + fr(E[1])*X*X + fr(E[3])*X + fr(E[4])
            if issq(y2):
                P = pari('[%s, %s]' % (X, sqrt_fr(y2)))
                if pari.ellorder(E, P) == 0: kp.append(P)
        kv = pari('[' + ','.join(str(P) for P in kp) + ']') if kp else pari('[]')
        for eff in (0, 1, 2):
            rk = pari.ellrank(E, eff, kv); lo, hi = int(rk[0]), int(rk[1])
            if lo == hi: break
        info = dict(rank=[lo, hi], effort=eff)
        gens = list(rk[3])
        for P in kp:                                   # the point of t = +-1, if PARI's list misses its direction
            cand = gens + [P]
            Hm = pari.ellheightmatrix(E, cand); d = float(pari.matdet(Hm)); scale = 1.0
            for i in range(len(cand)): scale *= float(Hm[i, i])
            if d > 1e-9 * scale: gens = cand
        if lo == hi <= RMAX and len(gens) == lo:
            if gens: gens = list(pari.ellsaturation(E, gens, 100))
            hmax = max([float(pari.ellheight(E, P)) for P in gens] or [0.0])
            pts, vecs = points(E, gens, K*K*hmax)
            info.update(searched=True, points=len(pts), H=round(K*K*hmax, 3))
            ts = set()
            for P in pts:
                for t in to_t(P): ts.add(t)
            for t in sorted((x for x in ts if x != 'inf'), key=lambda x: (abs(x.numerator) + x.denominator, x)):
                f = forms(m, n, t); k = sum(issq(v) for v in f)
                if k >= 4:
                    c = G.classify_t(m, n, t)
                    rec['hits'].append(dict(curve=name, t=str(t), squares=k, which=[int(issq(v)) for v in f], verdict=c))
                    if c == 'SOLUTION': rec['solutions'].append(str(t))
            found = [t for t in (Fr(1), Fr(-1)) if t in ts]
            info['planted'] = dict(carried=len(kp) > 0, found=[str(t) for t in found])
            if kp: assert found, ('planted test fails', m, n, name)        # the point of t = +-1 must be listed
        else:
            info['searched'] = False; info['gens'] = len(gens)
        rec['curves'][name] = info
    return rec

if __name__ == '__main__' and '--negative' in sys.argv:
    # the report path is live: with classify_t told that every t with five squares is a solution, the search must
    # report t = +-1 as SOLUTIONS on small square-type ratios
    real = G.classify_t
    G.classify_t = lambda m, n, t: 'SOLUTION' if real(m, n, t) != 'no' else 'no'
    lines = []
    for m, n in ((12, 5), (24, 7), (40, 9)):
        r = run(m, n, 4)
        lines.append('(%d,%d) planted solutions reported: %s; %s' % (m, n, sorted(set(r['solutions'])), 'as it must be' if set(r['solutions']) == {'1', '-1'} else 'FAILURE'))
    lines.append('# %d of %d as they must be' % (sum('as it must be' in l for l in lines), len(lines)))
    out = '\n'.join(lines) + '\n'; f = os.path.join(HERE, 'g4_negative.txt')
    if '--write' in sys.argv: open(f, 'w').write(out)
    print(out, end='')
    if '--write' not in sys.argv and os.path.exists(f): assert open(f).read() == out; print('# matches g4_negative.txt')
    sys.exit()
if __name__ == '__main__':
    a = [x for x in sys.argv[1:] if not x.startswith('--')]
    pairs = json.load(open(a[0]))['open'][int(a[1])::int(a[2])]
    K = float(a[3]) if len(a) > 3 else 8
    out = []
    for m, n in pairs:
        r = run(m, n, K); out.append(r)
        print(json.dumps(r), flush=True)
    if '--write' in sys.argv:
        json.dump(out, open(os.path.join(HERE, 'g4_%s_%s.json' % (a[1], a[2])), 'w'))
