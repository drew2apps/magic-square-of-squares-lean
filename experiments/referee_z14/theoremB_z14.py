"""z14 (scope extension): independent test of EXPLORATION_5 §1 "Theorem B":
    for a = (1, d, d), b = (1, beta, beta) in Sel^2 of y^2 = (x - e1)(x - e2)(x - e3), K = Q(sqrt d), beta = N(pi):
        <a, b> = prod_w (gamma_w, pi)_{K_w},   gamma = (u1 + sqrt d u3)(u1 + sqrt d u2)
    at local points of the covering u1^2 = x - e1, d u2^2 = x - e2, d u3^2 = x - e3,
against the repository's pairing (rank_second.full_selmer2, t10_pairct.DFullX; DFull on basis x basis extended
bilinearly, and DFull directly on 2 pairs per curve), on E_mix, E_08, E_LU of ratios outside the pilot (ratios.txt),
also with the roots relabelled (first component at e3).  Right side: local points from DFull's own search (its
local_xs), u_i by PARI p-adic square roots lifted to rationals (relative precision 2^60 resp. p^30), Hilbert symbols of
K by PARI's nfhilbert at every prime of K above 2, d, beta, z0, N(x0 + y0 sqrt d) and the e_i - e_j, real places by
signs (100-digit mpmath); each finite place evaluated at two local points (the factor must not change).
Written by an AI model (reviewer z14); not read by a person.   usage: theoremB_z14.py OUT.jsonl i0 i1"""
import sys, os, json, time, random
HERE = os.path.dirname(os.path.abspath(__file__))
ARGS = sys.argv[1:]; sys.argv = sys.argv[:1]
sys.path.insert(0, os.path.join(HERE, '..', 'conductor-t10'))
import t10_pairct as T
RS = T.RS; RI = T.RI
PARI = T.PARI
PARI.allocatemem(2 * 10 ** 8, 2 * 10 ** 9)
from fractions import Fraction as Fr
from math import gcd
import mpmath
mpmath.mp.dps = 100
NPAIRS = 40

def sqf(x): x = Fr(x); return RI.sqfree_part(x.numerator * x.denominator)
def primes(n):
    n = abs(int(n))
    if n < 2: return set()
    f = PARI.factor(n); return set(int(q) for q in f[0])

def rat_sqrt_padic(t, p, N):
    """a rational r with r^2 = t to relative precision p^N (t a non-zero p-adic square)"""
    t = Fr(t)
    s = PARI('sqrt((%d)/(%d) + O(%d^%d))' % (t.numerator, t.denominator, p, RS.vp(t.numerator, p) - RS.vp(t.denominator, p) + N))
    v = int(PARI.valuation(s, p)); u = PARI.lift(s * PARI(p) ** (-v))
    r = Fr(int(u)) * Fr(p) ** v
    assert r * r == t or RS.vp((r * r - t).numerator, p) - RS.vp((r * r - t).denominator, p) >= RS.vp(t.numerator, p) - RS.vp(t.denominator, p) + N - 2
    return r

def nfel(d, a, b):
    a, b = Fr(a), Fr(b)
    return PARI('Mod((%d)/(%d) + (%d)/(%d)*y, y^2 - (%d))' % (a.numerator, a.denominator, b.numerator, b.denominator, d))

def rhs(es, d, beta, D, rng, stats):
    """prod_w (gamma_w, pi)_w, or None if beta is not a norm from Q(sqrt d) (returns ('skip', why))."""
    s = PARI.qfsolve(PARI.matdiagonal([1, -d, -beta]))
    if str(PARI.type(s)) != 't_COL': return 'skip'
    x0, y0, z0 = [int(c) for c in s]
    assert x0 * x0 - d * y0 * y0 == beta * z0 * z0 and z0 != 0
    nf = PARI.nfinit(PARI('y^2 - (%d)' % d))
    pi = nfel(d, Fr(x0, z0), Fr(y0, z0))
    e1, e2, e3 = es; a = (1, d, d)
    S = {2}
    for q in (d, beta, z0, x0 * x0 - d * y0 * y0, e2 - e1, e3 - e1, e3 - e2): S |= primes(q)
    tot = 1
    for p in sorted(S):
        xs = D.local_xs(p, rng)
        facs = []
        for x in xs[:2]:
            N = 60 if p == 2 else 30
            u = [rat_sqrt_padic((Fr(x) - es[i]) / a[i], p, N) for i in range(3)]
            g = nfel(d, u[0], u[2]) * nfel(d, u[0], u[1])
            f = 1
            for pr in PARI.idealprimedec(nf, p):
                f *= int(PARI.nfhilbert(nf, g, pi, pr))
            facs.append(f)
        if len(set(facs)) > 1: stats['local factor depends on the point'] = stats.get('local factor depends on the point', 0) + 1
        tot *= facs[0]
    if d > 0:
        # a real point of the covering: x - e1 > 0, (x - e2)/d > 0, (x - e3)/d > 0
        xr = D.local_xs('inf', rng)[0]
        u = [mpmath.sqrt(mpmath.mpf((Fr(xr) - es[i]).numerator) / (Fr(xr) - es[i]).denominator / a[i]) for i in range(3)]
        sd = mpmath.sqrt(d)
        for sg in (1, -1):
            gv = (u[0] + sg * sd * u[2]) * (u[0] + sg * sd * u[1]); pv = (x0 + sg * sd * y0) / z0
            assert abs(gv) > mpmath.mpf(10) ** -50 and abs(pv) > mpmath.mpf(10) ** -50
            if gv < 0 and pv < 0: tot = -tot
    return tot

def coords_in(v, bv):
    piv = []
    for i, b in enumerate(bv):
        c = 1 << i
        for pr, pc in piv:
            if b ^ pr < b: b ^= pr; c ^= pc
        assert b
        piv.append((b, c)); piv.sort(reverse=True)
    c = 0
    for pr, pc in piv:
        if v ^ pr < v: v ^= pr; c ^= pc
    return c if v == 0 else None

def run_curve(name, es, rng, out, ratio):
    t0 = time.time()
    bad, loc, sel, kum, enc, places = RS.full_selmer2(es, random.Random(1))
    k = len(sel)
    D = [T.DFullX(es, a) for a in sel]
    M = [[None] * k for _ in range(k)]
    for i in range(k):
        for j in range(i, k):
            M[i][j] = M[j][i] = T.safe_pair(D[i], sel[j])[0]
    def elt(mask):
        c1 = c2 = 1
        for i in range(k):
            if mask >> i & 1: c1 *= sel[i][0]; c2 *= sel[i][1]
        return (sqf(c1), sqf(c2))
    def pred(x, y):
        s = 0
        for i in range(k):
            for j in range(k):
                if x >> i & 1 and y >> j & 1 and M[i][j] == -1: s ^= 1
        return -1 if s else 1
    V = [mask for mask in range(1, 2 ** k) if elt(mask)[0] == 1 and elt(mask)[1] != 1]
    pairs = [(x, y) for x in V for y in V]
    if len(pairs) > NPAIRS: pairs = rng.sample(pairs, NPAIRS)
    st = {}; bads = []
    direct = 0
    for x, y in sorted(pairs):
        d, beta = elt(x)[1], elt(y)[1]
        if d == beta: st['d = beta (trivial)'] = st.get('d = beta (trivial)', 0) + 1
        Da = T.DFullX(es, elt(x))
        r = rhs(es, d, beta, Da, rng, st)
        if r == 'skip': st['skipped: beta not a norm from Q(sqrt d)'] = st.get('skipped: beta not a norm from Q(sqrt d)', 0) + 1; continue
        want = pred(x, y)
        if direct < 2:
            h = T.safe_pair(Da, elt(y))[0]; direct += 1
            st['DFull direct = bilinear'] = st.get('DFull direct = bilinear', 0) + (h == want)
        key = 'agree' if r == want else 'DISAGREE'
        st[key] = st.get(key, 0) + 1
        st['pairing %d' % want] = st.get('pairing %d' % want, 0) + 1
        if key == 'DISAGREE' and len(bads) < 3: bads.append((d, beta, r, want))
    rec = dict(m=ratio[0], n=ratio[1], curve=name, es=list(es), dim=k, nV=len(V), npairs=len(pairs), stats=st, bad=bads, time=round(time.time() - t0, 1))
    with open(out, 'a') as f: f.write(json.dumps(rec, sort_keys=True) + '\n')
    sys.stderr.write('%s %s %s\n' % (ratio, name, st))

if __name__ == '__main__':
    out, i0, i1 = ARGS[0], int(ARGS[1]), int(ARGS[2])
    lines = [l.split() for l in open(os.path.join(HERE, 'ratios.txt')) if l.strip() and l[0] != '#']
    done = set()
    if os.path.exists(out):
        for l in open(out): r = json.loads(l); done.add((r['m'], r['n'], r['curve']))
    for l in lines[i0:i1]:
        m, n = map(int, l[0].split(','))
        rng = random.Random(7 * m + 3 * n)
        A, B, C = m * m - n * n, 2 * m * n, m * m + n * n
        es_mix = (-2 * C * C + 2 * A * B, 0, 4 * A * B)        # E_mix translated (middle root 0)
        curves = [('mix', es_mix), ('mix rotated', (es_mix[2], es_mix[0], es_mix[1])), ('08', (0, A * A, C * C)), ('LU', (0, B * B, C * C))]
        for name, es in curves:
            if (m, n, name) in done: continue
            run_curve(name, es, rng, out, (m, n))
