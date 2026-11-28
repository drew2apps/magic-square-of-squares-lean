import os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__))  # repository copy: paths relative to this folder
"""Task 3: the four-case table at odd p, against the local term computed from random local points (zct), for every
a = (1, d, d) in the local image Im_p and every b in Im_p, on E_08 and E_LU for many (m, n) incl. edge cases."""
import sys, random, collections, itertools
sys.path.insert(0, HERE)
import zct
from fractions import Fraction as Fr
from sympy import factorint, isprime
pari = zct.pari
def nonres(p):
    return next(e for e in range(2, p) if pow(e, (p - 1) // 2, p) == p - 1)
def creps(p):
    e = nonres(p); return [1, e, p, e * p]
def cbits(x, p):
    v, u = zct.cls(x, p); return (v % 2, 0 if pow(u, (p - 1) // 2, p) == 1 else 1)
def rep_of(bits, p):
    return {(0, 0): 1, (0, 1): nonres(p), (1, 0): p, (1, 1): nonres(p) * p}[bits]
def local_image(A, C, p, rng):
    """Im_p as a set of (bits(x - 0), bits(x - A^2)) for p odd"""
    es = (0, A * A, C * C)
    def kum(x):
        if x == es[0]: return (cbits((es[0] - es[1]) * (es[0] - es[2]), p), cbits(es[0] - es[1], p))
        if x == es[1]: return (cbits(es[1] - es[0], p), cbits((es[1] - es[0]) * (es[1] - es[2]), p))
        return (cbits(x - es[0], p), cbits(x - es[1], p))
    S = set([((0, 0), (0, 0))])
    def add(c):
        nonlocal S
        S |= {((a[0][0] ^ c[0][0], a[0][1] ^ c[0][1]), (a[1][0] ^ c[1][0], a[1][1] ^ c[1][1])) for a in S}
    for e in es: add(kum(e))
    for it in range(200000):
        if len(S) >= 4: break
        e = rng.choice(es); x = e + Fr(rng.randrange(1, p ** 5)) * Fr(p) ** rng.randrange(-6, 12) * rng.choice((1, -1))
        y2 = x * (x - es[1]) * (x - es[2])
        if y2 != 0 and zct.issq(y2, p): add(kum(x))
    return S
def hil(a, b, p): return zct.hil_vu(zct.cls(a, p), zct.cls(b, p), p)
def psqrt_cls(x, p):
    r = pari('sqrt((%d)+O(%d^40))' % (x, p)); return r
def hil_padic(r, b, p):
    v = int(pari.valuation(r, p)); u = r / pari(p) ** v
    return zct.hil_vu((v, int(pari.lift(pari.Mod(pari.truncate(u), p)))), zct.cls(b, p), p)
def table(A, C, dl, br, bm, bs, p):
    if (C - A) % p == 0: return 'p|C-A', 1
    if (C + A) % p == 0: return 'p|C+A', hil(dl, bs, p)
    if zct.issq(Fr(dl), p): return 'dl square', hil_padic(2 * psqrt_cls(dl, p), br, p)
    if zct.issq(Fr(-dl), p):
        r = psqrt_cls(-dl, p); return '-dl square', hil(-A, bs, p) * hil(-C, bm, p) * hil_padic((A - C) * r, br, p)
    return 'none', None
def local_term(A, C, dl, b, p, rng):
    F, npts = zct.local_forms(A, C, dl, p, rng)
    if not F: return 'nopoint'
    br, bm = b; bstar = br * bm
    bb = (zct.cls(br, p), zct.cls(bm, p), zct.cls(bstar, p))
    vals = set()
    for f in F:
        h = 1
        for i in range(3): h *= zct.hil_vu(f[i], bb[i], p)
        vals.add(h)
    return vals.pop() if len(vals) == 1 else 'nonconst'
def run(pairs, rng):
    cnt = collections.Counter(); ex = []; sol = collections.Counter()
    for m, n, tag in pairs:
        for cv in ('08', 'LU'):
            A, C = (m * m - n * n, m * m + n * n) if cv == '08' else (2 * m * n, m * m + n * n)
            if A == 0: continue
            B2 = C * C - A * A
            ps = sorted(q for q in set(factorint(abs(A * B2 * C))) if q != 2)
            for p in ps:
                Im = local_image(A, C, p, rng)
                if len(Im) != 4: cnt[(tag, 'IMAGE SIZE', len(Im))] += 1; continue
                for a in Im:
                    if a[0] != (0, 0): continue
                    dl = rep_of(a[1], p)
                    # solubility claims
                    if A % p == 0: sol[(tag, 'p|A: p∤δ and ±δ square', a[1][0] == 0 and (zct.issq(Fr(dl), p) or zct.issq(Fr(-dl), p)))] += 1
                    if C % p == 0: sol[(tag, 'p|C: δ square', zct.issq(Fr(dl), p))] += 1
                    for bcl in Im:
                        br, bm = rep_of(bcl[0], p), rep_of(bcl[1], p)
                        bs = br * bm
                        lt = local_term(A, C, dl, (br, bm), p, rng)
                        kind, pr = table(A, C, dl, br, bm, bs, p)
                        ok = (pr == lt)
                        cnt[(tag, cv, kind, 'agree' if ok else ('NOPRED' if pr is None else ('LT ' + str(lt) if isinstance(lt, str) else 'DISAGREE')))] += 1
                        if not ok and len(ex) < 12: ex.append((tag, cv, m, n, p, kind, 'dl', dl, 'b', (br, bm), 'table', pr, 'local', lt))
    return cnt, ex, sol
if __name__ == '__main__':
    rng = random.Random(11)
    pairs = []
    from math import gcd
    for m in range(2, 40):
        for n in range(1, m):
            if gcd(m, n) == 1: pairs.append((m, n, 'coprime'))
    pairs = rng.sample(pairs, 60) + [(465, 416, 'coprime'), (5623, 1017, 'coprime'), (3320, 2303, 'coprime'), (7, 4, 'coprime'), (4, 1, 'coprime'), (5, 2, 'coprime')]
    pairs += [(6, 3, 'noncoprime'), (9, 3, 'noncoprime'), (15, 5, 'noncoprime'), (10, 5, 'noncoprime'), (21, 14, 'noncoprime'), (35, 15, 'noncoprime'), (12, 9, 'noncoprime'), (22, 11, 'noncoprime')]
    pairs += [(n, m, 'm<n') for m, n, t in pairs[:15] if t == 'coprime']
    cnt, ex, sol = run(pairs, rng)
    for k, v in sorted(cnt.items(), key=str): print(k, v)
    print('solubility claims (tag, claim, holds):', dict(sol))
    for e in ex: print('  ', e)
