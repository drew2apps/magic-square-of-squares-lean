#!/usr/bin/env python3
"""z11 (review agent, run 9): independent recomputation of t10's pair-curve Cassels-Tate test.
Written by an AI model (review agent z11); not read by a person.  Experiment, not counted.

Own code: forms (from LAST_PAIR.md s.11), base points, Washington's map (checked symbolically in wash_check.py),
identification of the roots of E_jk with the three factorisations of F_j F_k by a FINITE-FIELD test (not t10's
resolvent formula), and the class of a tuple by the divisor form of the Kummer map:
    delta_r(psi P) = F_j(t_P) / F_j(t0),    delta_*(psi P) = +-a_j a_k(t_P) / a_j a_k(t0) (mod Q(i)*^2),
so c(tau) = (F_j(t0), +-d_j d_k D(t0), +-d_j d_k D'(t0)) with D, D' the rational classes of a_j a_k(t0), a_j conj a_k(t0).
Both signs are tested separately.
Reused (read, with care): z7's my2desc.Curve (complete 2-descent) and ct_mine.Cov (Cassels-Tate pairing), with the
primes of the discriminant from a PARI factorisation (primality checked) instead of sympy's factorint, and with the
real local points drawn in (e1, e2) and (e3, e3 + 3(e3 - e1)) (z7's draw misses them when the roots cluster);
G1's exact_local.py for the residual tuples (as G2 and t10 do).
PARI: factor/isprime of integers, factor of the cubic (roots checked exactly), sqrt mod p; ellrank only as a FINDER of
points in analyse(points=True) (every point checked; nothing relies on it); z7's code uses PARI p-adic square roots,
hilbert and qfsolve (points checked), and mpmath at the real place (as z7 says: a cross-check, not a certificate).

    python3 -I z11_pair.py        the degenerate controls (prints z11_pair.txt): every known tuple survives every
                                  pair curve, both signs; (12,5)'s extra tuples are killed on 0L, 0U, 8L, 8U
    python3 -I z11_pair.py m,n    any ratio (prints)
"""
import sys, os, random, json
from fractions import Fraction as Fr
from math import isqrt, gcd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _z11common import REPO, Out
import cypari2
pari = cypari2.Pari()
pari.allocatemem(10 ** 9, silent=True)
import my2desc as MY
import ct_mine as CT

def fac(x):
    x = abs(int(x))
    if x <= 1: return {}
    f = pari.factor(x); out = {int(p): int(e) for p, e in zip(f[0], f[1])}
    for p in out: assert pari.isprime(p)
    return out

def sqf(q, signed=True):
    q = Fr(q); assert q != 0
    a = abs(q.numerator * q.denominator); d = 1
    for p, e in fac(a).items():
        if e % 2: d *= p
    return -d if (signed and q < 0) else d

def issq(q):
    q = Fr(q)
    return q > 0 and isqrt(q.numerator) ** 2 == q.numerator and isqrt(q.denominator) ** 2 == q.denominator

# ---- forms: a_j(t) = A + B t, Gaussian (re, im)
def aform(m, n, j):
    return {'0': ((m, 0), (0, n)), '8': ((n, 0), (0, m)),
            'L': ((-(m + n), m - n), (m + n, m - n)), 'U': ((-(m - n), m + n), (m - n, m + n))}[j]
def aval(m, n, j, t):
    A, B = aform(m, n, j); t = Fr(t)
    return (A[0] + B[0] * t, A[1] + B[1] * t)
def gm(x, y): return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])
def cj(x): return (x[0], -x[1])
def N(x): return x[0] * x[0] + x[1] * x[1]
def Fcoef(m, n, j):
    A, B = aform(m, n, j)
    return [N(A), 2 * (A[0] * B[0] + A[1] * B[1]), N(B)]
def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for k, y in enumerate(b): out[i + k] += x * y
    return out
def rclass(z):
    """rational class (positive squarefree) of z in Q(i) with square norm (Lemma S1)"""
    w = Fr(isqrt((N(z)).numerator), isqrt((N(z)).denominator)); assert w * w == N(z)
    r = 2 * (z[0] + w)
    if r == 0: r = 2 * (z[0] - w)
    return sqf(r, signed=False)

BASE = {'08': Fr(0), 'LU': Fr(0), '0L': None, '0U': None, '8L': None, '8U': None}
def base_point(m, n, jk):
    return {'08': Fr(0), 'LU': Fr(0), '0L': Fr(m, n), '0U': Fr(-m, n), '8L': Fr(n, m), '8U': Fr(-n, m)}[jk]

class PairCurve:
    def __init__(s, m, n, jk):
        s.m, s.n, s.jk = m, n, jk; j, k = jk[0], jk[1]; s.j, s.k = j, k
        G = pmul(Fcoef(m, n, j), Fcoef(m, n, k))            # G(t) = sum G[i] t^i
        s.G = G
        t0 = base_point(m, n, jk); s.t0 = t0
        # g(s) = G(t0 + s)
        g = [Fr(0)] * 5
        for i, c in enumerate(G):
            p = [Fr(c)]
            for _ in range(i): p = pmul(p, [t0, Fr(1)])
            for ii, x in enumerate(p): g[ii] += x
        e, d, c, b, a = g
        assert e != 0 and issq(e), 'base point not on C'
        q = Fr(isqrt(e.numerator), isqrt(e.denominator)); s.g = g; s.q = q
        s.wc = (a, b, c, d, q)
        a1 = d / q; a2 = c - d * d / (4 * q * q); a3 = 2 * q * b; a4 = -4 * q * q * a; a6 = a2 * a4
        s.ai = (a1, a2, a3, a4, a6)
        A2, A4, A6 = a2 + a1 * a1 / 4, a4 + a1 * a3 / 2, a6 + a3 * a3 / 4
        mu = 1
        while not all((cc * mu ** (2 * i)).denominator == 1 for cc, i in ((A2, 1), (A4, 2), (A6, 3))): mu += 1
        s.mu = mu
        C2, C4, C6 = int(A2 * mu ** 2), int(A4 * mu ** 4), int(A6 * mu ** 6)
        s.cub = (C2, C4, C6)
        f = pari('x^3 + (%d)*x^2 + (%d)*x + (%d)' % (C2, C4, C6))
        fa = pari.factor(f); es = []
        for i in range(len(fa[0])):
            pol = fa[0][i]; assert int(pari.poldegree(pol)) == 1
            es += [Fr(str(-pari.polcoef(pol, 0) / pari.polcoef(pol, 1)))] * int(fa[1][i])
        assert len(set(es)) == 3 and all(x.denominator == 1 for x in es)
        s.es = sorted(int(x) for x in es)
        for r in s.es: assert r ** 3 + C2 * r * r + C4 * r + C6 == 0
    # psi on points over a field: x on the short integral model
    def X_of(s, sv, v):
        a, b, c, d, q = s.wc
        return s.mu ** 2 * (2 * q * (v + q) + d * sv) / sv ** 2
    def from_X(s, X, Y):
        """(X, Y) on the integral model -> (t, v) on C (t = t0 + s), exact; None for O"""
        a1, a2, a3, a4, a6 = s.ai; a, b, c, d, q = s.wc
        x = Fr(X) / s.mu ** 2; y = Fr(Y) / s.mu ** 3 - (a1 * x + a3) / 2
        if y == 0: return None
        sv = (2 * q * (x + c) - d * d / (2 * q)) / y
        v = -q + sv * (sv * x - d) / (2 * q)
        assert v * v == sum(cc * sv ** i for i, cc in enumerate(s.g)), 'inverse map wrong'
        return s.t0 + sv, v
    def iota_base_X(s):
        a1, a2, a3, a4, a6 = s.ai
        return s.mu ** 2 * (-a2)
    # the three factorisations, as functions of t (Gaussian values over Q, or mod p with i -> ip)
    def parts(s, t):
        aj, ak = aval(s.m, s.n, s.j, t), aval(s.m, s.n, s.k, t)
        return {'r': N(aj), '*': gm(aj, ak), "'": gm(aj, cj(ak))}
    def ffield_identify(s, primes):
        """for each root e_i and each factorisation P: is chi(X - e_i) == chi(q_P(t) q_P(t0)) at every F_p point?"""
        res = {}
        for p in primes:
            ip = int(pari.lift(pari.sqrt(pari.Mod(-1, p))))
            def red(x): x = Fr(x); return x.numerator * pow(x.denominator, -1, p) % p
            def chi(x): x %= p; return 0 if x == 0 else (1 if pow(x, (p - 1) // 2, p) == 1 else -1)
            gp = [red(c) for c in s.g]; a, b, c, d, q = [red(x) for x in s.wc]; mu = s.mu % p
            t0p = red(s.t0)
            def qp(t, P):
                aj = [red(x) for x in aform(s.m, s.n, s.j)[0]], [red(x) for x in aform(s.m, s.n, s.j)[1]]
                A, B = aform(s.m, s.n, s.j); C, D = aform(s.m, s.n, s.k)
                ajr, aji = (A[0] + B[0] * t) % p, (A[1] + B[1] * t) % p
                akr, aki = (C[0] + D[0] * t) % p, (C[1] + D[1] * t) % p
                if P == 'r': return (ajr * ajr + aji * aji) % p
                zj = (ajr + ip * aji) % p
                zk = (akr + ip * aki) % p if P == '*' else (akr - ip * aki) % p
                return zj * zk % p
            q0 = {P: qp(t0p, P) for P in ('r', '*', "'")}
            ok = {(i, P): set() for i in range(3) for P in q0}
            npts = 0
            for sv in range(1, p):
                gv = sum(cc * pow(sv, i, p) for i, cc in enumerate(gp)) % p
                if gv == 0 or chi(gv) != 1: continue
                r0 = int(pari.lift(pari.sqrt(pari.Mod(gv, p))))
                for v in (r0, p - r0):
                    X = mu * mu * (2 * q * (v + q) + d * sv) * pow(sv * sv, -1, p) % p
                    t = (t0p + sv) % p
                    for i, e in enumerate(s.es):
                        ce = chi(X - e)
                        if ce == 0: continue
                        for P in q0:
                            cq = chi(qp(t, P) * q0[P])
                            if cq == 0: continue
                            ok[(i, P)].add(ce == cq)
                    npts += 1
            for key, v in ok.items(): res.setdefault(key, set()).update(v)
        match = {}
        for i in range(3):
            good = [P for P in ('r', '*', "'") if res[(i, P)] == {True}]
            assert len(good) == 1, ('identification not unique', i, good, res)
            match[good[0]] = i
        assert len(match) == 3
        s.role = match
        return match
    def tuple_class(s, tp, sign):
        """(c_e1, c_e2, c_e3) for a tuple tp = (d0, d8, dL, dU)"""
        IDX = {'0': 0, '8': 1, 'L': 2, 'U': 3}
        dd = tp[IDX[s.j]] * tp[IDX[s.k]]
        P0 = s.parts(s.t0)
        cr = sqf(P0['r']); cs = sign * sqf(dd * rclass(P0['*'])); cp = sqf(cr * cs)
        # independent check of the third coordinate (up to sign): +-dd D'(t0)
        assert abs(cp) == sqf(dd * rclass(P0["'"]), signed=False), 'third coordinate inconsistent'
        out = [None] * 3
        out[s.role['r']] = cr; out[s.role['*']] = cs; out[s.role["'"]] = cp
        return tuple(out)

class Curve2(MY.Curve):
    def __init__(self, e1, e2, e3):
        self.e = [int(e1), int(e2), int(e3)]
        assert len(set(self.e)) == 3
        a, b, c = self.e
        Nn = 2 * abs((a - b) * (a - c) * (b - c))
        self.primes = sorted(fac(Nn))
        self.S = [-1] + self.primes

_orig_lp = CT.Cov.local_points
def _local_points(self, p, rng, want):
    """z7's search, but at the real place x is drawn in (e1, e2) and (e3, e3 + 3(e3 - e1)) (z7's draw in
    [-10H, 10H] misses the real locus when the roots cluster)."""
    if p != 0: return _orig_lp(self, p, rng, want)
    e1, e2, e3 = sorted(self.es); out = []
    for _ in range(200000):
        if len(out) >= want: break
        if rng.randint(0, 1): x = Fr(e1) + Fr(rng.randint(1, 10 ** 6 - 1), 10 ** 6) * (e2 - e1)
        else: x = Fr(e3) + Fr(rng.randint(1, 3 * 10 ** 6), 10 ** 6) * (e3 - e1)
        if self.ok(x, p) and x not in out: out.append(x)
    assert len(out) == want, ('no local points', p)
    return out
CT.Cov.local_points = _local_points

class Pairing:
    def __init__(s, es, seed=1):
        s.es = list(es)
        s.E = Curve2(*es)
        s.dim = s.E.selmer(seed=seed + 100)
        s.sel = [(CT.sqf(a), CT.sqf(b)) for a, b in s.E.sel]
        s.pl = CT.places_of(s.E, (3, 5, 7, 11, 13))
        s.rng = random.Random(seed)
        s.covs = {}
    def insel(s, c): return s.E.sel_contains((Fr(c[0]), Fr(c[1])))
    def cov(s, c):
        if c not in s.covs: s.covs[c] = CT.Cov(s.es, c, s.rng)
        return s.covs[c]
    def pair(s, b, c):
        """<b, D_c> with the local values"""
        D = s.cov(c); tot = 1; per = {}
        for p in s.pl:
            v = D.local_value(b, p, s.rng, 4)
            per[p] = v; tot *= v
            if p not in s.E.primes and p != 0: assert v == 1
        return tot, [p for p, v in per.items() if v == -1]
    def test(s, c):
        """'outside' | ('killed', b, places) | 'kernel'"""
        c = (int(c[0]), int(c[1]))
        if not hasattr(s, 'rc'): s.rc = {}
        if c in s.rc: return s.rc[c]
        s.rc[c] = s.test0(c); return s.rc[c]
    def test0(s, c):
        if not s.insel(c): return ('outside',)
        for b in s.sel:
            v, pls = s.pair(b, c)
            if v == -1: return ('killed', b, pls)
        return ('kernel',)

def residual(m, n, pmax=400):
    import exact_local as X
    TD = X.TD
    forms, sup, known = TD.magic(m, n)
    divs = []
    for S in sup:
        ds = [1]
        for qq in S: ds += [x * qq for x in ds]
        divs.append(ds)
    bad = set()
    for x in (2 * m * n * (m * m - n * n), m * m + n * n, m ** 8 + 14 * m ** 4 * n ** 4 + n ** 8): bad |= set(fac(x))
    loc = []
    for p in sorted(bad):
        L, al = X.allowed_at(forms, p)
        if len(al) < len(L.labels()) ** len(forms): loc.append((p, L, al))
    left = X.join(divs, loc)
    kill = {}
    alive = [t for t in left if t not in known]
    for qq in pari.primes([3, pmax]):
        qq = int(qq)
        if qq in bad or not alive: continue
        L, al = X.allowed_at(forms, qq)
        still = []
        for t in alive:
            if tuple(L.dlabel(d) for d in t) in al: still.append(t)
            else: kill[t] = qq
        alive = still
    return left, kill, known

PAIRS = ['08', '0L', '0U', '8L', '8U', 'LU']
FFP = [101, 109, 113]

def analyse(m, n, tuples=None, pairs=PAIRS, points=True, log=print):
    if tuples is None:
        left, kill, known = residual(m, n)
    else:
        left, known = tuples; kill = {}
    log('ratio (%d,%d): %d tuples (%d known): %s' % (m, n, len(left), len(known), sorted(left)))
    out = {t: {} for t in left}
    for jk in pairs:
        C = PairCurve(m, n, jk)
        primes = [p for p in FFP if all(Fr(x).denominator % p for x in C.g) and C.mu % p and C.q.numerator % p]
        role = C.ffield_identify(primes)
        Xi = C.iota_base_X()
        di = [sqf(Xi - e) if Xi != e else None for e in C.es]
        PR = Pairing(C.es)
        # rational points (finder: PARI ellrank), checked; their classes against the divisor formula
        ptinfo = ''
        if points:
            Ep = pari.ellinit([0, C.cub[0], 0, C.cub[1], C.cub[2]])
            pari.setrand(1); r = pari.ellrank(Ep)
            npt = 0
            for P in r[3]:
                X, Y = Fr(str(P[0])), Fr(str(P[1]))
                assert Y * Y == (X - C.es[0]) * (X - C.es[1]) * (X - C.es[2])
                tv = C.from_X(X, Y)
                if tv is None: continue
                t, v = tv
                if X in C.es: continue
                parts = C.parts(t); P0 = C.parts(C.t0)
                assert sqf(X - C.es[role['r']]) == sqf(parts['r'] * P0['r']), 'divisor formula fails at r'
                assert abs(sqf(X - C.es[role['*']])) == sqf(rclass(parts['*']) * rclass(P0['*']), signed=False), 'fails at *'
                cpt = (sqf(X - C.es[0]), sqf(X - C.es[1]))
                assert PR.insel(cpt), 'point image outside my Sel2'
                # the pairing of a point image with every Selmer basis element is +1
                for b in PR.sel: assert PR.pair(b, cpt)[0] == 1, 'a point image pairs to -1'
                npt += 1
            ptinfo = '; %d point(s) found by PARI checked (divisor formula, in my Sel2, pair +1 with all of Sel2)' % npt
        log('  E_%s: t0 = %s; e = %s; roles %s (F_p test at %s); delta(psi(iota t0)) = %s; my dim Sel2 = %d%s'
            % (jk, C.t0, C.es, {k: 'e%d' % (v + 1) for k, v in sorted(role.items())}, primes,
               {P: di[i] for P, i in role.items()}, PR.dim, ptinfo))
        for t in left:
            res = []
            for sign in (1, -1):
                c3 = C.tuple_class(t, sign)
                c = (c3[0], c3[1])
                res.append((c, PR.test(c)))
            v = [r[1][0] for r in res]
            assert v[0] == v[1] or set(v) <= {'kernel', 'killed'}, v
            out[t][jk] = res
            log('    %s: class(+) %s: %s | class(-) %s: %s' % (t, res[0][0], res[0][1], res[1][0], res[1][1]))
    log('  summary (%d,%d):' % (m, n))
    for t in left:
        ks = [jk for jk in pairs if out[t][jk][0][1][0] in ('killed', 'outside')]
        ks2 = [jk for jk in pairs if out[t][jk][1][1][0] in ('killed', 'outside')]
        log('    %s known=%s goodprime=%s killed(+) by %s; killed(-) by %s' % (t, t in known, kill.get(t), ks, ks2))
    return out, kill, known

CONTROLS = [(4, 3), (12, 5), (15, 8), (23, 7), (21, 20), (31, 17)]

if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    if args:
        for a in args:
            m, n = map(int, a.split(',')); analyse(m, n)
    else:
        say = Out('z11_pair.txt')
        say('z11_pair.py: degenerate controls (written by an AI model, review z11; not read by a person). m^2 + n^2 is a')
        say('square or twice a square; the known tuples are classes of rational points and must survive every pair curve.')
        ok = True
        for m, n in CONTROLS:
            out, kill, known = analyse(m, n, log=say)
            for t in known:
                for jk in PAIRS:
                    ok &= all(r[1][0] == 'kernel' for r in out[t][jk])
        say('ALL KNOWN TUPLES SURVIVE, BOTH SIGNS' if ok else 'A KNOWN TUPLE IS KILLED: THE METHOD IS WRONG')
        say.finish()
