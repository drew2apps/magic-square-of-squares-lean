#!/usr/bin/env python3
"""search.py (conductor-t13): point search on the 2-covering quartics y^2 = G(t1, t2) of covers.json over Z[i].
EXPERIMENT; written by an AI model (build agent t13), not read by a person.  Standard library and numpy only.

    python3 search.py --curve=E --H=100000 [--classes=coset] [--write]   -> search_E_H100000.txt / .json

REGION.  covers.py stores, for each quartic, the covariant point (z, T) of hyperbolic 3-space (reduced into the
fundamental domain of PSL2(Z[i])).  The search takes every pair (t1, t2) of Gaussian integers with
    |t1 - z t2|^2 / T + T |t2|^2 <= H,
t2 = 0 (then only t1 = 1) or t2 in the quadrant re t2 > 0, im t2 >= 0 (one t2 per class modulo the units; (t1, t2)
and (u t1, u t2) give the same point).  This is an ellipsoid adapted to G: for a reduced G it is, up to a bounded
factor, the set of points of P^1(K) of naive height (relative to the covariant) at most H.  Pairs with
gcd(t1, t2) != 1 are not removed (they repeat points; harmless).

SIEVE.  For split primes p = 1 mod 4 (13 <= p < 1000, both primes r and -r above p, those where G mod p is not
identically zero) a table allowed[t2 mod p][t1 mod p] = [G(t1, t2) is a square mod p] (0 counts as a square) is
precomputed; a pair survives if it is allowed at every prime ideal used.  Survivors are tested EXACTLY: G(t1, t2) is
evaluated in Z[i] with Python integers and its square root is extracted exactly (norm a perfect square, then
a^2 = (N + x)/2, b^2 = (N - x)/2).  Every exact hit is written out; certify.py maps it back to E and re-checks it
with Sage (on the curve over K, Kummer images).  The sieve can only remove pairs where G(t1, t2) is a non-square
modulo a prime ideal, i.e. where W^2 = G(t1, t2) has no solution in Z[i]: nothing with a solution is lost.
"""
import sys, os, json, time, math
from fractions import Fraction as Fr
from math import isqrt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ARGS = dict(a.lstrip('-').split('=', 1) if '=' in a else (a.lstrip('-'), '1') for a in sys.argv[1:])
WRITE = 'write' in ARGS
OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(str(s))

def gint(pair):
    a, b = Fr(pair[0]), Fr(pair[1]); assert a.denominator == 1 and b.denominator == 1
    return (int(a), int(b))
def gmul(x, y): return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])
def gadd(x, y): return (x[0] + y[0], x[1] + y[1])

def G_exact(g, t1, t2):
    """G(t1, t2) = sum g_k t1^(4-k) t2^k in Z[i] (Horner in t1)."""
    p2 = [(1, 0)]
    for k in range(4): p2.append(gmul(p2[-1], t2))
    v = g[0]
    for k in range(1, 5):
        v = gadd(gmul(v, t1), gmul(g[k], p2[k]))
    return v

def gauss_sqrt(z):
    x, y = z
    N2 = x * x + y * y; N = isqrt(N2)
    if N * N != N2: return None
    if (N + x) % 2: return None
    a2, b2 = (N + x) // 2, (N - x) // 2
    a, b = isqrt(a2), isqrt(b2)
    if a * a != a2 or b * b != b2: return None
    for w in ((a, b), (a, -b)):
        if gmul(w, w) == (x, y): return w
    return None

def gdivmod(x, y):
    n = y[0] * y[0] + y[1] * y[1]
    num = gmul(x, (y[0], -y[1]))
    q = ((2 * num[0] + n) // (2 * n), (2 * num[1] + n) // (2 * n))
    r = (x[0] - gmul(q, y)[0], x[1] - gmul(q, y)[1])
    return q, r

def ggcd(x, y):
    while y != (0, 0):
        x, y = y, gdivmod(x, y)[1]
    return x

def primitive(t1, t2):
    """(t1, t2) divided by gcd(t1, t2), then by the unit that puts t2 (or t1 if t2 = 0) in re > 0, im >= 0."""
    g = ggcd(t1, t2)
    t1 = gdivmod(t1, g)[0]; t2 = gdivmod(t2, g)[0]
    lead = t2 if t2 != (0, 0) else t1
    for u in ((1, 0), (0, 1), (-1, 0), (0, -1)):
        l = gmul(lead, u)
        if l[0] > 0 and l[1] >= 0:
            return (gmul(t1, u), gmul(t2, u))
    raise AssertionError

def split_primes(lo, hi):
    out = []
    for p in range(lo, hi):
        if p % 4 != 1 or any(p % q == 0 for q in range(2, isqrt(p) + 1)): continue
        r = next(r for r in range(2, p) if (r * r + 1) % p == 0)
        out += [(p, r), (p, p - r)]
    return out

def tables(g, primes, maxn=44):
    tabs = []
    for p, r in primes:
        gk = [(a + b * r) % p for a, b in g]
        if all(c == 0 for c in gk): continue
        sq = np.zeros(p, dtype=bool); sq[(np.arange(p, dtype=np.int64) ** 2) % p] = True
        X = np.arange(p, dtype=np.int64)[None, :]; Y = np.arange(p, dtype=np.int64)[:, None]
        v = np.full((p, p), gk[0], dtype=np.int64)
        Yp = np.ones((p, 1), dtype=np.int64)
        for k in range(1, 5):
            Yp = (Yp * Y) % p
            v = (v * X + gk[k] * Yp) % p
        allowed = sq[v]                       # allowed[t2 mod p, t1 mod p]
        tabs.append((p, r, allowed, float(allowed.mean())))
        if len(tabs) >= maxn: break
    tabs.sort(key=lambda t: t[3])             # the most selective first
    return tabs

def search_one(g, zc, T, H, chunk=1 << 20, H0=0.0):
    tabs = tables(g, split_primes(13, 1000))
    hits = []; npairs = 0; nsurv = 0; nt2 = 0
    # t2 = 0: the point (1 : 0)
    if 1.0 / T > H0:                          # h_p(1, 0) = 1/T
        npairs += 1
        if gauss_sqrt(G_exact(g, (1, 0), (0, 0))) is not None: hits.append(((1, 0), (0, 0)))
    r2max = H / T
    cmax = isqrt(int(r2max)) + 1
    t2s = [(c, d) for c in range(1, cmax + 1) for d in range(0, cmax + 1) if c * c + d * d <= r2max]
    for (c, d) in t2s:
        R2 = H * T - T * T * (c * c + d * d)
        if R2 < 0: continue
        nt2 += 1
        R02 = H0 * T - T * T * (c * c + d * d)    # the inner ellipsoid h_p <= H0 (already searched) is skipped
        cen = zc * complex(c, d); R = math.sqrt(R2)
        a0, a1 = math.ceil(cen.real - R), math.floor(cen.real + R)
        b0, b1 = math.ceil(cen.imag - R), math.floor(cen.imag + R)
        avals = np.arange(a0, a1 + 1, dtype=np.int64); na = len(avals)
        da2 = (avals - cen.real) ** 2
        rows = max(1, chunk // max(1, na))
        rowtabs = []
        for p, r, allowed, fr in tabs:
            rowtabs.append((p, allowed[(c + d * r) % p], (avals % p).astype(np.int32), r))
        for bs in range(b0, b1 + 1, rows):
            bvals = np.arange(bs, min(b1, bs + rows - 1) + 1, dtype=np.int64)
            db2 = (bvals - cen.imag) ** 2
            dd = da2[None, :] + db2[:, None]
            disc = (dd <= R2) & (dd > R02)
            npairs += int(disc.sum())
            p, row, am, r = rowtabs[0]
            bm = ((bvals * r) % p).astype(np.int32)
            m = disc & row[(am[None, :] + bm[:, None]) % p]
            ib, ia = np.nonzero(m)
            for p, row, am, r in rowtabs[1:]:
                if len(ia) == 0: break
                bm = ((bvals * r) % p).astype(np.int32)
                keep = row[(am[ia] + bm[ib]) % p]
                ib, ia = ib[keep], ia[keep]
            nsurv += len(ia)
            for j in range(len(ia)):
                t1 = (int(avals[ia[j]]), int(bvals[ib[j]])); t2 = (c, d)
                if gauss_sqrt(G_exact(g, t1, t2)) is not None:
                    hits.append((t1, t2))
    return dict(npairs=npairs, nt2=nt2, nsurv=nsurv, hits=hits, nprimes=len(tabs))

def main():
    cv = ARGS.get('curve', 'E'); H = float(ARGS.get('H', '10000')); only = ARGS.get('classes', '')
    H0 = float(ARGS.get('H0', '0'))
    src = ARGS.get('models', 'covers.json')         # or altmodels.json
    data = json.load(open(os.path.join(HERE, src)))['curves'][cv]
    if H0:
        say('search.py (conductor-t13): curve %s, H0 < h_p <= H with H0 = %g, H = %g.  Written by an AI model (build agent t13); not read by a person.' % (cv, H0, H))
    else:
        say('search.py (conductor-t13): curve %s, H = %g.  Written by an AI model (build agent t13); not read by a person.' % (cv, H))
    res = {}
    for cl in data['classes']:
        if only and not cl['label'].startswith(only): continue
        if ARGS.get('skipbase') and cl['label'].endswith('m00'): continue    # the start model (searched from covers.json)
        g = [gint(x) for x in cl['G']]
        zc = complex(cl['covariant'][0], cl['covariant'][1]); T = cl['covariant'][2]
        t0 = time.time()
        r = search_one(g, zc, T, H, H0=H0)
        dt = time.time() - t0
        prim = sorted(set(primitive(t1, t2) for t1, t2 in r['hits']))
        say('  %-11s T = %.4g: %d values of t2, %d pairs (t1, t2), %d sieve primes, %d sieve survivors, %d exact hits, %d primitive' % (
            cl['label'], T, r['nt2'], r['npairs'], r['nprimes'], r['nsurv'], len(r['hits']), len(prim)))
        print('     (%.0f s)' % dt, file=sys.stderr, flush=True)
        r['hits'] = prim
        for t1, t2 in r['hits'][:20]:
            say('     hit t1 = %d%+di, t2 = %d%+di' % (t1[0], t1[1], t2[0], t2[1]))
        res[cl['label']] = {'npairs': r['npairs'], 'hits': r['hits']}
    if WRITE:
        tag = 'search_%s%s_H%d%s%s' % (cv, '_alt' if src != 'covers.json' else '', int(H), ('_from%d' % int(H0)) if H0 else '', ('_' + only.replace(' ', '').replace('+', '')) if only else '')
        open(os.path.join(HERE, tag + '.txt'), 'w').write('\n'.join(OUT) + '\n')
        json.dump(res, open(os.path.join(HERE, tag + '.json'), 'w'), indent=1, sort_keys=True)

if __name__ == '__main__':
    main()
