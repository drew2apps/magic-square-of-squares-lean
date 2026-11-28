"""t10: are the good-prime kills of the four-class descent shadowed by a Cassels-Tate obstruction on a pair curve?
(EXPERIMENT, written by an AI model; not read by a person; nothing here changes the census.  Tier D.)

THE SETTING.  For a ratio m : n the four-class descent (tuple_descent.magic; LAST_PAIR.md sections 11-12) attaches
to a solution t = u/v the tuple (d0, d8, dL, dU) of rational classes of the Gaussian linear forms (h1_pairs.py's
convention, t = u/v)
    a0 = m + i n t,  a8 = n + i m t,  aL = (m + n)(t - 1) + i (m - n)(t + 1),  aU = (m - n)(t - 1) + i (m + n)(t + 1),
(tuple_descent's forms are v a_j, or v conj(a_8): the same products d_j d_k).  RESIDUAL tuples: those allowed at every
bad prime (every prime of 2mn(m^2 - n^2), m^2 + n^2 and m^8 + 14 m^4 n^4 + n^8; odd primes by G1's exact_local.py,
2 by the repository's test, as G2 does); the real place imposes nothing (the forms are positive definite).  For each
residual tuple the KILLING PRIME is the least good prime at which its labels are not allowed (exact_local.py).

THE PAIR CURVES.  For two forms j, k the curve C_jk: y^2 = G(t) = F_j(t) F_k(t) (F = N(a)) has genus one; its
Jacobian E_jk has full rational 2-torsion, the three points T_r, T*, T' belonging to the factorisations
{F_j, F_k} (over Q), {a_j a_k, conj} and {a_j conj a_k, conj} (over Q(i)).  A point of the four-form curve with tuple
tau lies on C_jk, with F_j(t) a rational square and a_j a_k(t) of class D = d_j d_k.
  BASE POINT.  The method needs a rational point t0 of C_jk: t0 = 0 for (0,8) and (L,U); for the other pairs a point
  found by PARI's hyperellratpoints (a FINDER; the point is checked).  t = (alpha s + beta)/(gamma s + delta) moves
  t0 to s = 0, and the map of Washington (Elliptic Curves, Thm 2.17) gives psi: C_jk -> E_jk, (0, q) -> O, then an
  integral model y^2 = (x - e1)(x - e2)(x - e3).  Every point used is mapped back and tested on the quartic.
  ROOTS.  For the 2-torsion point T_b = (e_b, 0), P_b = psi^-1(T_b) is the point of C with P_b ~ (zeros of g1) - iota O
  + ..., so that g1(0) g2(s_b) = g2(0) g1(s_b) holds for the factorisation {g1, g2} of T_b (derived in the notes; here
  a test): this identifies e_r (for {F_j, F_k}) and e* (for {a_j a_k, conj}) exactly; then
  (x - e_r) psi = kappa_r F_j(t) and (x - e*) psi = kappa* D(t) (the second modulo -1, since -1 is a square in Q(i))
  modulo squares, with constants kappa read off every point found (2-torsion, generators, sums): all must agree
  (asserted), and the identification of the roots must be the only one consistent with the points (asserted).
  THE CLASS OF A TUPLE.  A tuple point P would have delta(psi(P)) = c(tau) := (x - e_r : kappa_r, x - e* : +-kappa* D,
  x - e' : the product) in H^1(Q, E[2]) = (Q*/Q*^2)^2.  The two signs differ by delta(psi(iota t0)), the class of a
  rational point, so they behave alike below (checked).
  THE TEST (as H1's "Cassels-Tate kernel test"): Sel^2(E_jk) by rank_second.full_selmer2 (complete 2-descent, exact);
  the images of the points found (PARI's ellrank as a finder; checked on the curve) span a subgroup of
  delta(E(Q)); the Cassels-Tate pairing <b, D_c> by rank_second.DFull.pairing (exact, conic points from PARI's
  qfsolve as hints, checked).  The class c(tau) of a tuple point is in delta(E(Q)), hence in the kernel of the
  pairing.  So tau is PAIRING-KILLED at (j, k) if c(tau) is not in Sel^2, or if <b, D_c(tau)> = -1 for some b in
  Sel^2 (b runs over a basis of Sel^2 modulo the point images; the pairing with a point image is 1, checked on one).
  Otherwise c(tau) is in the kernel: "image" if it is the class of a point found, "kernel" if not.
What it rests on: the complete 2-descent and the Cassels-Tate pairing as rank_second.py computes them (Cassels 1998;
Skorobogatov 6.2; reviews/review_t2.md), exact_local.py and the repository's test at 2 for the residual sets (tested,
not proved), and the theory quoted above (Washington's map, the Kummer map).  No PARI output is relied on.
usage: t10_pairct.py RATIO [RATIO ...]           e.g. 465,416 (prints)
       t10_pairct.py --list FILE [--sample N] [--seed S] [--from i] [--to j] [--out NAME] [--write]   (t10_NAME.txt)
       t10_pairct.py --negative [--write]        (t10_negative.txt: controls and forged inputs)
       t10_pairct.py --h1control [--write]       (t10_h1control.txt: (1109503,597671) against h1_pairs.txt)
       t10_pairct.py --deep 'm,n;m,n' --effort k [--write]   (t10_deep.txt: more points for the tuples not killed)
       t10_pairct.py --summarize f1,f2,... [--write]          (t10_summary.txt: the per-tuple lines counted again)
Hints (base points, PARI's points, conic points) are read from every t10_hints*.json here and written, with --write,
to t10_hints_NAME.json.  Timing goes to stderr.
The repository is found through MSS_REPO, or else ../..; G1's exact_local.py in experiments/roadmap-g/g1."""
import sys, os, json, random, time
from fractions import Fraction as Fr
from math import isqrt, gcd
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get('MSS_REPO', os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(REPO, 'scripts', 'g1_full'))
sys.path.insert(0, os.path.join(REPO, 'experiments', 'roadmap-g', 'g1'))
import rank_indep as RI
import rank_second as RS
import exact_local as X
TD = X.TD
import cypari2
PARI = cypari2.Pari(); PARI.allocatemem(5 * 10 ** 8, 4 * 10 ** 9)          # the stack grows as needed
WRITE = '--write' in sys.argv
def _arg(name, default=None):
    return sys.argv[sys.argv.index(name) + 1] if name in sys.argv else default
HFILE = os.path.join(HERE, 't10_hints%s.json' % ('_' + _arg('--out') if _arg('--out') else ''))
HINTS = {}
for _f in sorted(os.listdir(HERE)):
    if _f.startswith('t10_hints') and _f.endswith('.json'):
        for _k, _v in json.load(open(os.path.join(HERE, _f))).items(): HINTS.setdefault(_k, {}).update(_v)
RS.HINTS['conic'].update(HINTS.get('rs_conic', {}))
RI.HINTS.setdefault('conic', {}).update(HINTS.get('ri_conic', {}))
RI.HINTS.setdefault('factor', {})
OUT = []
def say(s=''):
    OUT.append(s); print(s, flush=True)

def factorint_proved(x):
    x = abs(int(x))
    if x <= 1: return {}
    f = PARI.factor(x); out = {int(q): int(e) for q, e in zip(f[0], f[1])}
    for q in out: assert PARI.isprime(q), ('not proved prime', q)
    return out
TD.factorint = factorint_proved

# ------------------------------------------------------------------------------------------- the four-class descent
def badprimes(m, n):
    bad = set()
    for x in (2 * m * n * (m * m - n * n), m * m + n * n, m ** 8 + 14 * m ** 4 * n ** 4 + n ** 8): bad |= set(factorint_proved(x))
    return sorted(bad)

def residual(m, n, pmax=5000):
    """(residual tuples, {tuple: killing good prime or None}, bad primes, known tuples)."""
    forms, sup, known = TD.magic(m, n)
    divs = []
    for S in sup:
        ds = [1]
        for q in S: ds += [x * q for x in ds]
        divs.append(ds)
    bad = badprimes(m, n)
    loc = []
    for p in bad:
        L, al = X.allowed_at(forms, p)
        for kt in known: assert tuple(L.dlabel(d) for d in kt) in al, ('known point rejected', m, n, p)
        if len(al) < len(L.labels()) ** len(forms): loc.append((p, L, al))
    left = X.join(divs, loc)
    kill = {}; alive = [t for t in left if t not in known]
    for kt in known: kill[kt] = None
    for q in PARI.primes([3, pmax]):
        q = int(q)
        if q in bad or not alive: continue
        L, al = X.allowed_at(forms, q)
        for kt in known: assert tuple(L.dlabel(d) for d in kt) in al, ('known point rejected', m, n, q)
        still = []
        for t in alive:
            if tuple(L.dlabel(d) for d in t) in al: still.append(t)
            else: kill[t] = q
        alive = still
    for t in alive: kill[t] = None
    return left, kill, bad, known

# ------------------------------------------------------------------------------------------------ the pair curves
IDX = {'0': 0, '8': 1, 'L': 2, 'U': 3}
PAIRS = [('0', '8'), ('0', 'L'), ('0', 'U'), ('8', 'L'), ('8', 'U'), ('L', 'U')]
def aform(m, n, j):
    """a_j(t) = (r0 + r1 t) + i (i0 + i1 t) as ((r0, r1), (i0, i1))"""
    return {'0': ((m, 0), (0, n)), '8': ((n, 0), (0, m)), 'L': ((-(m + n), m + n), (m - n, m - n)),
            'U': ((-(m - n), m - n), (m + n, m + n))}[j]
def a_hom(m, n, j, u, v):
    (r0, r1), (i0, i1) = aform(m, n, j); return (r0 * v + r1 * u, i0 * v + i1 * u)
def gm(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def gconj(a): return (a[0], -a[1])
def N(a): return a[0] * a[0] + a[1] * a[1]
def isq(q): q = Fr(q); return q >= 0 and isqrt(q.numerator) ** 2 == q.numerator and isqrt(q.denominator) ** 2 == q.denominator
def rsqrt(q): q = Fr(q); assert isq(q); return Fr(isqrt(q.numerator), isqrt(q.denominator))
def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i + j] += x * y
    return out

FB = set()
def sqf(q, signed=True):
    """squarefree part of a non-zero rational (with sign if signed): primes of the factor base FB first, the rest by
    a proved factorisation."""
    q = Fr(q); assert q != 0
    a = abs(q.numerator * q.denominator); d = 1
    for p in sorted(FB):
        if a % p: continue
        e = 0
        while a % p == 0: a //= p; e += 1
        if e % 2: d *= p
    if not (isqrt(a) ** 2 == a):
        for p, e in factorint_proved(a).items():
            if e % 2: d *= p
    return -d if (signed and q < 0) else d

class Pair:
    """C_jk with a rational base point moved to s = 0, Washington's map to E_jk, an integral model with full
    2-torsion (as h1_pairs.Pair, for a base point that need not be t = 0 or an integer)."""
    def __init__(s, m, n, j, k, t0=None):
        s.m, s.n, s.j, s.k = m, n, j, k
        s.Gt = s.quartic_t()
        if t0 is None: t0 = s.find_base()
        s.t0 = t0
        # Moebius t = (al s + be) / (ga s + de) with (be : de) = t0
        if t0 == 'inf': al, be, ga, de = 0, 1, 1, 0
        else:
            t0 = Fr(t0); be, de = t0.numerator, t0.denominator
            g, x, y = egcd(be, de)          # x be + y de = 1
            al, ga = -y, x                  # al de - be ga = -y de - x be = -1
        s.M = (al, be, ga, de)
        co = [0] * 5
        for i, c in enumerate(s.Gt):        # sum c t^i -> sum c (al s + be)^i (ga s + de)^(4 - i)
            p = [c]
            for _ in range(i): p = pmul(p, [be, al])
            for _ in range(4 - i): p = pmul(p, [de, ga])
            for ii, x in enumerate(p): co[ii] += x
        s.co = co
        assert co[0] != 0 and isq(co[0]), 'the base point is not a point of C_jk'
        e, d, c, b, a = [Fr(x) for x in co]
        s.q = q = rsqrt(e)
        s.ai = [d / q, c - d * d / (4 * q * q), 2 * q * b, -4 * q * q * a, (c - d * d / (4 * q * q)) * (-4 * q * q * a)]
        a1, a2, a3, a4, a6 = s.ai
        A2, A4, A6 = a2 + a1 * a1 / 4, a4 + a1 * a3 / 2, a6 + a3 * a3 / 4
        u = 1
        while not all((cc * u ** (2 * i)).denominator == 1 for cc, i in ((A2, 1), (A4, 2), (A6, 3))): u += 1
        s.u = u; s.A = (int(A2 * u ** 2), int(A4 * u ** 4), int(A6 * u ** 6))
        f = PARI('x^3 + (%d)*x^2 + (%d)*x + (%d)' % s.A)
        fac = PARI.factor(f)
        es = []
        for i in range(len(fac[0])):
            pol = fac[0][i]
            assert int(PARI.poldegree(pol)) == 1, 'full rational 2-torsion expected'
            es += [int(-PARI.polcoef(pol, 0) / PARI.polcoef(pol, 1))] * int(fac[1][i])
        es = tuple(sorted(es))
        assert len(set(es)) == 3
        for r in es: assert r ** 3 + s.A[0] * r * r + s.A[1] * r + s.A[2] == 0
        s.es = es
    def quartic_t(s):
        m, n, j, k = s.m, s.n, s.j, s.k
        def Fp(jj):
            (r0, r1), (i0, i1) = aform(m, n, jj)
            return [r0 * r0 + i0 * i0, 2 * (r0 * r1 + i0 * i1), r1 * r1 + i1 * i1]
        return pmul(Fp(j), Fp(k))
    def find_base(s):
        G = s.Gt
        if G[0] != 0 and isq(G[0]): return 0
        if G[4] != 0 and isq(G[4]): return 'inf'
        key = '%d,%d %s%s' % (s.m, s.n, s.j, s.k)
        if key in HINTS.get('base', {}):
            t0 = HINTS['base'][key]
            return 'inf' if t0 == 'inf' else Fr(t0)
        for t in (Fr(s.m, s.n), Fr(-s.m, s.n), Fr(s.n, s.m), Fr(-s.n, s.m)):   # t = m/n is a point of C_0L, etc.
            v = G[0] + G[1] * t + G[2] * t ** 2 + G[3] * t ** 3 + G[4] * t ** 4
            if v != 0 and isq(v): return t
        pol = PARI('Pol(%s)' % str(list(reversed(G))).replace(' ', ''))
        best = None
        for B in (30, 300, 3000):
            pts = PARI.hyperellratpoints(pol, B)
            for P in pts:
                t = Fr(str(P[0]))
                if G[0] + G[1] * t + G[2] * t ** 2 + G[3] * t ** 3 + G[4] * t ** 4 == 0: continue
                h = max(abs(t.numerator), t.denominator)
                if best is None or h < best[0]: best = (h, t)
            if best: break
        if best is None: raise LookupError('no base point found')
        HINTS.setdefault('base', {})[key] = str(best[1])
        return best[1]
    def on(s, P): x, y = P; return y * y == (x - s.es[0]) * (x - s.es[1]) * (x - s.es[2])
    def to_s(s, P):
        """(x, y) on the integral model -> s on the quartic in s (None for the point over the base point), checked."""
        a1, a2, a3, a4, a6 = s.ai; e, d, c, b, a = [Fr(x) for x in s.co]
        x = Fr(P[0]) / s.u ** 2; y = Fr(P[1]) / s.u ** 3 - (a1 * x + a3) / 2
        if y == 0: return None
        sv = (2 * s.q * (x + c) - d * d / (2 * s.q)) / y
        assert isq(sum(cc * sv ** i for i, cc in enumerate(s.co))), 'the map to the quartic is wrong'
        return sv
    def iota_base(s):
        """psi(iota t0) = psi((0, -q)) on the integral model: on Washington's model x = d^2/(4q^2) - c, y = 4 q^2 v3
        with v = -q + v1 s + v2 s^2 + v3 s^3 + ... the branch through (0, -q) (the limit of Washington's formulas);
        checked on the curve."""
        a1, a2, a3, a4, a6 = s.ai; e, d, c, b, a = [Fr(x) for x in s.co]; q = s.q
        v1 = -d / (2 * q); v2 = (v1 * v1 - c) / (2 * q); v3 = (2 * v1 * v2 - b) / (2 * q)
        x = d * d / (4 * q * q) - c; y = 4 * q * q * v3
        assert y * y + a1 * x * y + a3 * y == x ** 3 + a2 * x * x + a4 * x + a6, 'psi(iota t0) is not on the curve'
        P = (x * s.u ** 2, (y + (a1 * x + a3) / 2) * s.u ** 3)
        assert s.on(P) and s.to_s(P) in (0, None)
        return P
    def point_uv(s, P):
        """(x, y) on the integral model -> (u, v) on C_jk (t = u / v), or None for the base point; the points over
        s = infinity (Washington's y = 0) included."""
        a1, a2, a3, a4, a6 = s.ai; e, d, c, b, a = [Fr(x) for x in s.co]
        x = Fr(P[0]) / s.u ** 2; y = Fr(P[1]) / s.u ** 3 - (a1 * x + a3) / 2
        num = 2 * s.q * (x + c) - d * d / (2 * s.q)
        if y == 0:
            if num == 0: return None
            assert isq(a), 'a point at s = infinity, but the leading coefficient is not a square'
            al, be, ga, de = s.M; g = gcd(al, ga); return (al // g, ga // g)
        return s.to_uv(s.to_s(P))
    def phis(s):
        """the roots of the resolvent cubic, a (th1 th2 + th3 th4), for the pairings r = {F_j, F_k}, * = {a_j a_k,
        conj}, ' = {a_j conj a_k, conj} (exact, in the s coordinate)."""
        al, be, ga, de = s.M; out = {}
        def root(jj):
            (r0, r1), (i0, i1) = aform(s.m, s.n, jj)
            A = (r0 * de + r1 * be, i0 * de + i1 * be); B = (r0 * ga + r1 * al, i0 * ga + i1 * al)   # a_jj = A + B s
            nb = Fr(N(B)); assert nb != 0
            z = gm((-A[0], -A[1]), gconj(B)); return (Fr(z[0]) / nb, Fr(z[1]) / nb)
        tj, tk = root(s.j), root(s.k); a = Fr(s.co[4])
        out['r'] = a * (N(tj) + N(tk))
        out['*'] = 2 * a * gm(tj, tk)[0]
        out["'"] = 2 * a * gm(tj, gconj(tk))[0]
        return out
    def to_uv(s, sv):
        """s -> (u, v) integers with t = u / v (v = 0 for t = infinity)."""
        al, be, ga, de = s.M
        u, v = al * sv + be, ga * sv + de
        L = u.denominator * v.denominator // gcd(u.denominator, v.denominator)
        u, v = int(u * L), int(v * L); g = gcd(u, v); return u // g, v // g
    def add(s, P, Q):
        A2, A4, A6 = s.A
        if P is None: return Q
        if Q is None: return P
        (x1, y1), (x2, y2) = P, Q
        if x1 == x2:
            if y1 + y2 == 0: return None
            l = (3 * x1 * x1 + 2 * A2 * x1 + A4) / (2 * y1)
        else: l = (y2 - y1) / (x2 - x1)
        x3 = l * l - A2 - x1 - x2; y3 = -(y1 + l * (x3 - x1))
        return (x3, y3)

def egcd(a, b):
    if b == 0: return (a, 1, 0) if a >= 0 else (-a, -1, 0)
    g, x, y = egcd(b, a % b); return g, y, x - (a // b) * y

def classes_at(Pr, uv):
    """(F_j class with sign, D = class of a_j a_k (mod -1)) at the point (u : v) of C_jk."""
    m, n, j, k = Pr.m, Pr.n, Pr.j, Pr.k
    u, v = uv
    aj, ak = a_hom(m, n, j, u, v), a_hom(m, n, k, u, v)
    Fj = N(aj)
    z = gm(aj, ak); r = rsqrt(N(z)); x = z[0] + r
    if x == 0: x = z[0] - r
    return sqf(Fj), sqf(2 * x, signed=False)

def find_points(Pr):
    key = '%d,%d %s%s %s' % (Pr.m, Pr.n, Pr.j, Pr.k, Pr.t0)
    if key in HINTS.get('points', {}):
        h = HINTS['points'][key]
        return [(Fr(x), Fr(y)) for x, y in h['pts']], h['pari']
    E = PARI.ellinit([0, Pr.A[0], 0, Pr.A[1], Pr.A[2]])
    r = PARI.ellrank(E)
    pts = [(Fr(str(P[0])), Fr(str(P[1]))) for P in r[3]]
    tors = PARI.elltors(E)
    pts += [(Fr(str(P[0])), Fr(str(P[1]))) for P in tors[2] if len(P) == 2]
    pari = [int(r[0]), int(r[1]), int(r[2]), [int(x) for x in tors[1]]]
    HINTS.setdefault('points', {})[key] = {'pts': [[str(x), str(y)] for x, y in pts], 'pari': pari}
    return pts, pari

class DFullX(RS.DFull):
    """rank_second.DFull with a wider search for local points when its own search finds none (the pairing does not
    depend on the local points chosen; rank_second checks that on every place)."""
    def local_xs(s, p, rng, want=RS.NPTS):
        try: return RS.DFull.local_xs(s, p, rng, want)
        except AssertionError: pass
        out = []; lo, hi = min(s.es), max(s.es); w = hi - lo + 1
        for _ in range(400000):
            if len(out) >= want: break
            if p == 'inf': x = Fr(rng.randrange(lo - 4 * w, hi + 4 * w), rng.randrange(1, 60))
            else:
                e = s.es[rng.randrange(3)]
                x = e + Fr(rng.randrange(1, p ** 8)) * Fr(p) ** rng.randrange(-20, 40) * rng.choice((1, -1))
                if x.denominator % p == 0 and RS.vp(x.denominator, p) % 2: continue
            if s.ok_x(x, p) and x not in out: out.append(x)
        assert len(out) == want, ('no local point found', p)
        s.cache[p] = out
        return out

def safe_pair(D, b):
    try: return D.pairing(b)
    except AssertionError as ex:
        if 'no local point' in str(ex): return None, str(ex)
        raise

class PairData:
    """Everything about one pair curve that does not depend on the tuple."""
    def __init__(s, m, n, j, k, forge=None):
        s.m, s.n, s.j, s.k = m, n, j, k
        Pr = s.Pr = Pair(m, n, j, k)
        es = Pr.es
        G, s.pari = find_points(Pr)
        for P in G: assert Pr.on(P), 'a hint point is not on the curve'
        Pi = Pr.iota_base()
        if Pi[1] != 0 and Pi not in G and (Pi[0], -Pi[1]) not in G: G = G + [Pi]
        # rank_second's local point search draws x near the roots and in [-8H, 8H]: shift the model so that the
        # middle root is 0 (the classes x - e_i do not change)
        s.sh = es[1]; s.esh = tuple(e - s.sh for e in es)
        s.bad, loc, sel, kumsh, enc, places = RS.full_selmer2(s.esh, random.Random(1))
        s.kum = lambda x: kumsh(x - s.sh)
        s.sel = sel
        s.gens = [-1] + s.bad
        T2 = [(Fr(e), Fr(0)) for e in es]
        # points (for the constants): the combinations sum c_i P_i, c_i in {-1, 0, 1}, of the points found and the
        # 2-torsion, at most 80
        neg = lambda P: None if P is None else (P[0], -P[1])
        pts = [None]
        for g in list(G) + T2[:2]:
            new = []
            for P in pts:
                for Q in (Pr.add(P, g), Pr.add(P, neg(g))):
                    if Q is not None and Q not in pts and Q not in new: new.append(Q)
            pts = (pts + new)[:80]
        pts = [P for P in pts if P is not None]
        def kum3(P):
            if P[1] == 0:
                b = es.index(int(P[0])); o = [e for e in es if e != es[b]]
                c = [Fr(es[b] - o[0]) * (es[b] - o[1]) if i == b else Fr(int(P[0]) - es[i]) for i in range(3)]
                return c
            return [Fr(P[0]) - e for e in es]
        # the roots e_r, e*: the exact test on the 2-torsion points (see the docstring)
        al, be, ga, de = Pr.M
        def gj(jj, sv):  # a_jj in the s coordinate (homogeneous in (al s + be, ga s + de))
            return a_hom(m, n, jj, al * sv + be, ga * sv + de)
        # primary: the affine map between the roots of the resolvent cubic and the e's (exactly one permutation)
        ph = Pr.phis(); from itertools import permutations
        fits = []
        for pi in permutations(range(3)):
            (p1, p2, p3), (f1, f2, f3) = [Fr(es[pi[i]]) for i in range(3)], (ph['r'], ph['*'], ph["'"])
            if f1 == f2: continue
            lam = (p1 - p2) / (f1 - f2); mu = p1 - lam * f1
            if lam != 0 and p3 == lam * f3 + mu: fits.append(pi)
        assert len(fits) == 1, ('the resolvent does not identify the roots', fits)
        role = {'r': fits[0][0], '*': fits[0][1], "'": fits[0][2]}
        s.role_resolvent = dict(role)
        if forge == 'wrong root': role['r'], role['*'] = role['*'], role['r']
        # check: the 2-torsion points that map to a finite s (the factorisation test of the docstring)
        s.t2checked = 0
        for b, T in enumerate(T2):
            sv = Pr.to_s(T)
            if sv is None: continue
            aj0, ak0, ajs, aks = gj(j, 0), gj(k, 0), gj(j, sv), gj(k, sv)
            tests = {'r': N(aj0) * N(aks) == N(ak0) * N(ajs),
                     '*': gm(gm(aj0, ak0), gconj(gm(ajs, aks)))[1] == 0,
                     "'": gm(gm(aj0, gconj(ak0)), gconj(gm(ajs, gconj(aks))))[1] == 0}
            hit = [kk for kk, v in tests.items() if v]
            assert len(hit) == 1 and s.role_resolvent[hit[0]] == b, ('the 2-torsion point disagrees with the resolvent', b, hit)
            s.t2checked += 1
        # the constants, from every point (all must agree); roots not identified above: by agreement
        data = []
        FB.clear(); FB.update(s.bad); FB.update(q for x in (m, n, m - n, m + n, m * m + n * n) for q in factorint_proved(x))
        for P in pts:
            uv = Pr.point_uv(P)
            if uv is None: continue
            fj, D = classes_at(Pr, uv)
            data.append((kum3(P), fj, D))
        def consts(b, which):
            vals = set()
            for c, fj, D in data:
                vals.add(sqf(c[b] * fj) if which == 'r' else sqf(c[b] * D, signed=False))
            return vals
        for which in ('r', '*'):
            ok = [b for b in range(3) if len(consts(b, which)) == 1]
            assert role[which] in ok, ('the identified root disagrees with the points', which)
        s.npts = len(data)
        assert role['r'] != role['*'] and s.npts > 0, 'no point to read the constants from'
        s.role = role
        s.kr = consts(role['r'], 'r').pop(); s.ks = consts(role['*'], '*').pop()
        # uniqueness of the identification on the points (the other roots must not be consistent too)
        s.unique = (len([b for b in range(3) if len(consts(b, 'r')) == 1]) == 1, len([b for b in range(3) if len(consts(b, '*')) == 1]) == 1)
        # point images
        s.imgs = [s.kum(P[0]) for P in G] + [s.kum(e) for e in es[:2]]
        s.imgv = [s.vec(c) for c in s.imgs]
        s.selv = [s.vec(c) for c in sel]
        assert len(RS.f2_span(s.selv + s.imgv)) == len(sel), 'a point image outside Sel^2: the descent is wrong'
        s.ri = len(RS.f2_span(s.imgv))
        s.pcache = {}; s.kill_places = []; s.dcache = {}
        s.bad4 = set(badprimes(m, n)) | {2}
        s.badE_in_bad4 = set(s.bad) <= s.bad4
    def dfull(s, c):
        key = tuple(c)
        if key not in s.dcache: s.dcache[key] = DFullX(s.esh, c)
        return s.dcache[key]
    def vec(s, c):
        """(c1, c2) -> bit vector over the generators [-1] + bad (None if a prime outside)."""
        v = 0; L = len(s.gens)
        for jj, d in enumerate(c):
            d = Fr(d); d = RI.sqfree_part(d.numerator * d.denominator); rest = abs(d)
            for i, g in enumerate(s.gens):
                if (g == -1 and d < 0) or (g != -1 and d % g == 0): v |= 1 << (i + jj * L)
            for g in s.bad:
                if rest % g == 0: rest //= g
            if rest != 1: return None
        return v
    def cls_of_tuple(s, tp, sign=1):
        """c(tau) as (x - e1, x - e2) classes."""
        D = sqf(Fr(tp[IDX[s.j]] * tp[IDX[s.k]]), signed=False)
        comp = [None] * 3
        comp[s.role['r']] = s.kr
        comp[s.role['*']] = sign * s.ks * D
        o = 3 - s.role['r'] - s.role['*']
        comp[o] = comp[s.role['r']] * comp[s.role['*']]
        return (RI.sqfree_part(comp[0]), RI.sqfree_part(comp[1]))
    def test(s, c, forge=None):
        """('outside Sel2' | 'killed' | 'image' | 'kernel', detail)"""
        v = s.vec(c)
        if v is None or len(RS.f2_span(s.selv + [v])) != len(s.sel): return 'outside Sel2', None
        if len(RS.f2_span(s.imgv + [v])) == s.ri: return 'image', None
        key = tuple(c)
        if key in s.pcache: return s.pcache[key]
        # basis of Sel^2 modulo the images
        comp = []
        base = list(RS.f2_span(s.imgv))
        for b, bv in zip(s.sel, s.selv):
            if len(RS.f2_span(base + [bv])) > len(base): base = RS.f2_span(base + [bv]); comp.append(b)
        Dc = s.dfull(c)
        vals = []
        for b in comp:
            h, per = safe_pair(Dc, b)
            vals.append((b, h, per))
        # checks: alternating on c, trivial on a point image
        h_cc, _ = safe_pair(Dc, c)
        assert h_cc in (1, None), 'the pairing is not alternating'
        nz = [x for x in s.imgs if s.vec(x) != 0]
        if nz:
            cc = (RI.sqfree_part(Fr(nz[0][0]).numerator * Fr(nz[0][0]).denominator), RI.sqfree_part(Fr(nz[0][1]).numerator * Fr(nz[0][1]).denominator))
            assert safe_pair(Dc, cc)[0] in (1, None), 'a point image pairs non-trivially'
        bad = [(b, per) for b, h, per in vals if h == -1]
        fail = [per for b, h, per in vals if h is None]
        if forge == 'no pairing': bad = []
        if bad:
            b, per = bad[0]
            # symmetry check on the killing pair
            hb, _ = safe_pair(s.dfull(b), c)
            assert hb in (-1, None), 'the pairing is not symmetric'
            places = [p for p, x in per.items() if x == -1]
            s.kill_places.append(places)
            res = ('killed', 'b = %s, local -1 at %s%s' % (b, places, '' if hb == -1 else ' (symmetry not checked: %s)' % _))
        elif fail: res = ('undecided', 'local points not found: %s' % fail)
        else: res = ('kernel', 'paired with %d classes' % len(comp))
        s.pcache[key] = res
        return res

def analyse_ratio(m, n, pairs=PAIRS, forge=None, tuples=None):
    t1 = time.time()
    if tuples is not None: left, known = tuples; kill = {t: 'not computed' for t in left}; left = list(known) + [t for t in left if t not in known]
    else: left, kill, bad, known = residual(m, n)
    say('ratio (%d,%d): %d residual tuple(s) after the bad primes; known tuples %d' % (m, n, len(left), len(known)))
    res = {t: {} for t in left}
    for j, k in pairs:
        try: Pd = PairData(m, n, j, k, forge=forge)
        except LookupError as ex:
            say('  E_%s%s: %s' % (j, k, ex))
            for t in left: res[t][j + k] = 'no base'
            continue
        say('  E_%s%s: base t0 = %s; e = %s; e_r = e%d, e* = e%d (resolvent; 2-torsion checks %d; constants agree on %d points, alone unique: %s); bad primes of E in those of the ratio: %s; dim Sel2 = %d, point images %d; PARI ellrank [%d,%d] s=%d tors %s'
            % (j, k, Pd.Pr.t0, Pd.Pr.es, Pd.role['r'] + 1, Pd.role['*'] + 1, Pd.t2checked, Pd.npts, Pd.unique, Pd.badE_in_bad4, len(Pd.sel), Pd.ri, *Pd.pari))
        STATS['curves'] += 1; STATS['badE_outside'] += (not Pd.badE_in_bad4)
        for t in left:
            Pd.kill_places = []
            c = Pd.cls_of_tuple(t); c2 = Pd.cls_of_tuple(t, -1)
            # the two signs differ by delta(psi(iota t0)), the class of a rational point: checked here, so that the
            # pairing of the one is the pairing of the other
            v1, v2 = Pd.vec(c), Pd.vec(c2)
            assert (v1 is None) == (v2 is None)
            if v1 is not None: assert len(RS.f2_span(Pd.imgv + [v1 ^ v2])) == Pd.ri, 'the two signs do not differ by a point image'
            r1 = Pd.test(c, forge=forge)
            v = r1[0]
            res[t][j + k] = v
            if v in ('killed', 'outside Sel2'):
                say('    tuple %s: class %s: %s %s' % (t, c, v, r1[1] or ''))
                for pl in Pd.kill_places[:1]:
                    for p in pl: PLACES[p if p == 'inf' else ('bad' if p in Pd.bad4 else 'GOOD')] = PLACES.get(p if p == 'inf' else ('bad' if p in Pd.bad4 else 'GOOD'), 0) + 1
            if t in known: assert v in ('image', 'kernel'), 'a tuple of a rational point is killed: the method is wrong'
    say('  per tuple (killing good prime; pairs):')
    for t in left:
        ks = [jk for jk, v in res[t].items() if v in ('killed', 'outside Sel2')]
        say('    %s  good prime %s  %s  => %s' % (t, kill.get(t), ' '.join('%s:%s' % (jk, v) for jk, v in res[t].items()),
                                                 ('PAIRING-KILLED by ' + ','.join(ks)) if ks else 'not pairing-killed'))
    sys.stderr.write('  time %.1f s\n' % (time.time() - t1))
    return left, kill, res, known

STATS = {'curves': 0, 'badE_outside': 0}
PLACES = {}

def read_list(path):
    out = []
    for l in open(path):
        l = l.strip()
        if not l or l[0] == '#' or 'left' not in l: continue
        m, n = map(int, l.split()[0].split(',')); out.append((m, n))
    return out

def run_list(pairs, name):
    say('t10_pairct.py: %s, %d ratios. Written by an AI model; not read by a person. Experiment, not counted; tier D.' % (name, len(pairs)))
    tot = {'ratios': 0, 'tuples': 0, 'pk': 0, 'goodkill': 0, 'pk_and_good': 0, 'pk_not_good': 0, 'good_not_pk': 0, 'neither': 0}
    perpair = {}; survivors = []; ratios_all = 0
    for m, n in pairs:
        left, kill, res, known = analyse_ratio(m, n)
        tot['ratios'] += 1; allk = True
        for t in left:
            if t in known: continue
            tot['tuples'] += 1
            ks = [jk for jk, v in res[t].items() if v in ('killed', 'outside Sel2')]
            for jk in ks: perpair[jk] = perpair.get(jk, 0) + 1
            g = kill.get(t) is not None
            tot['pk'] += bool(ks); tot['goodkill'] += g
            tot['pk_and_good'] += bool(ks) and g; tot['pk_not_good'] += bool(ks) and not g
            tot['good_not_pk'] += g and not ks; tot['neither'] += (not g) and not ks
            if not ks: survivors.append((m, n, t, kill.get(t), res[t])); allk = False
        ratios_all += allk
    say('')
    say('SUMMARY (%s): %d ratios, %d residual tuples (known tuples of degenerate ratios left out)' % (name, tot['ratios'], tot['tuples']))
    say('  killed by a good prime <= 5000: %d; pairing-killed (some pair curve): %d' % (tot['goodkill'], tot['pk']))
    say('  good-prime-killed and pairing-killed: %d; good-prime-killed, NOT pairing-killed: %d' % (tot['pk_and_good'], tot['good_not_pk']))
    say('  pairing-killed, not killed by a good prime <= 5000: %d; neither: %d' % (tot['pk_not_good'], tot['neither']))
    say('  ratios all of whose residual tuples are pairing-killed: %d of %d' % (ratios_all, tot['ratios']))
    say('  kills per pair curve: %s' % ' '.join('%s:%d' % (jk, perpair.get(jk, 0)) for jk in [a + b for a, b in PAIRS]))
    say('  places carrying a -1 (first killing pair of each kill): %s' % PLACES)
    say('  pair curves: %d; with a bad prime of E outside those of the ratio: %d' % (STATS['curves'], STATS['badE_outside']))
    for m, n, t, g, r in survivors: say('  NOT pairing-killed: (%d,%d) %s good prime %s %s' % (m, n, t, g, r))
    return tot

CONTROLS = [(4, 3), (7, 1), (12, 5), (17, 7), (627, 364), (991, 263)]

def deep(m, n, effort):
    """For the tuples of (m, n) that no pair curve kills: more points of each E_jk (PARI's ellrank with a higher
    effort, a FINDER; checked on the curve); a class in the span of their images is the class of a rational point of
    E_jk, so no obstruction on that pair curve (at any level) can remove the tuple."""
    left, kill, res, known = analyse_ratio(m, n)
    surv = [t for t in left if not any(v in ('killed', 'outside Sel2') for v in res[t].values())]
    say('deep (%d,%d), ellrank effort %d: %d tuple(s) not pairing-killed' % (m, n, effort, len(surv)))
    for j, k in PAIRS:
        Pd = PairData(m, n, j, k)
        E = PARI.ellinit([0, Pd.Pr.A[0], 0, Pd.Pr.A[1], Pd.Pr.A[2]])
        r = PARI.ellrank(E, effort)
        pts = [(Fr(str(P[0])), Fr(str(P[1]))) for P in r[3]]
        for P in pts: assert Pd.Pr.on(P), 'a point is not on the curve'
        imgv = Pd.imgv + [Pd.vec(Pd.kum(P[0])) for P in pts]
        ri = len(RS.f2_span(imgv))
        say('  E_%s%s: PARI ellrank [%d,%d] s=%d, %d point(s); images span %d of dim Sel2 = %d' % (j, k, int(r[0]), int(r[1]), int(r[2]), len(pts), ri, len(Pd.sel)))
        for t in surv:
            v = Pd.vec(Pd.cls_of_tuple(t))
            img = len(RS.f2_span(imgv + [v])) == ri
            say('    %s: %s' % (t, 'the class of a rational point of E_%s%s' % (j, k) if img else 'not in the span of the points found (%s)' % res[t][j + k]))

def h1_control():
    """Cross-check with H1 (h1_pairs.txt): (1109503,597671), its two known tuples and four extra tuples (G1), with
    h1_pairs's points of E_0L, E_0U as hints (the same quartics and the same integral models: base point t = 0)."""
    say('t10_pairct.py --h1control: (1109503,597671), the tuples of h1_pairs.py. Written by an AI model; not read by a person.')
    H1 = os.path.join(REPO, 'experiments', 'roadmap-h', 'h1', 'h1_hints.json')
    h = json.load(open(H1))
    RS.HINTS['conic'].update(h.get('rs_conic', {}))
    for key, pts in h.get('points', {}).items():
        mm, jk = key.split()
        HINTS.setdefault('points', {}).setdefault('%s %s 0' % (mm, jk), {'pts': pts, 'pari': [-1, -1, -1, []]})
    M_, N_ = 1109503, 597671
    KNOWN = [(1109503, 597671, 2, 1), (1195342, 2219006, 2, 1)]
    EXTRA = [(1109503, 597671, 2, 959), (1109503, 597671, 1918, 1), (1195342, 2219006, 2, 959), (1195342, 2219006, 1918, 1)]
    left, kill, res, known = analyse_ratio(M_, N_, pairs=[('0', 'L'), ('0', 'U')], tuples=(EXTRA, KNOWN))
    # h1_pairs.txt: E_0U kills (.., 2, 959), (1195342, .., 2, 959); E_0L kills (.., 1918, 1) and (1195342, .., 1918, 1)
    want = {EXTRA[0]: {'0U'}, EXTRA[1]: {'0L'}, EXTRA[2]: {'0U'}, EXTRA[3]: {'0L'}}
    agree = all({jk for jk, v in res[t].items() if v in ('killed', 'outside Sel2')} == want[t] for t in EXTRA)
    agree &= all(res[t][jk] in ('image', 'kernel') for t in KNOWN for jk in res[t])
    say('AGREES WITH h1_pairs.txt' if agree else 'DIFFERS FROM h1_pairs.txt')

def negative():
    say('t10_pairct.py --negative: forged inputs that must be refused, and controls that must survive.')
    say('Written by an AI model; not read by a person. Experiment, not counted; tier D.')
    ok_all = True
    # 1. degenerate ratios: their known tuples (classes of rational points) must survive every pair curve
    for m, n in CONTROLS:
        left, kill, res, known = analyse_ratio(m, n)
        surv = all(res[t][jk] in ('image', 'kernel') for t in known for jk in res[t])
        ok_all &= surv
        say('  1. (%d,%d): %d known tuple(s): %s' % (m, n, len(known), 'all survive, as they must' if surv else 'ONE IS KILLED'))
    # 2. a forged class: a residual class of (465,416) times a prime outside the bad primes of E_0L must be refused
    m, n = 465, 416
    left, kill, bad, known = residual(m, n)
    Pd = PairData(m, n, '0', 'L')
    q = next(int(p) for p in PARI.primes(100) if int(p) not in Pd.bad)
    c = Pd.cls_of_tuple(left[0]); cf = (c[0] * q, c[1])
    r = Pd.test(cf)
    ok_all &= r[0] == 'outside Sel2'
    say('  2. forged class %s (the class of %s times %d in x - e1): %s' % (cf, left[0], q, 'refused (outside Sel^2), as it must be' if r[0] == 'outside Sel2' else 'NOT REFUSED: ' + r[0]))
    # 3. the roots e_r and e* exchanged: the constants must disagree on the points
    try:
        PairData(m, n, '0', 'L', forge='wrong root'); say('  3. roots exchanged: NOT REFUSED'); ok_all = False
    except AssertionError as ex: say('  3. roots exchanged: refused (%s), as it must be' % str(ex).split(',')[0].strip("('"))
    # 4. a forged base point: t = 1 is not a point of C_0L for (465,416)
    try:
        Pair(m, n, '0', 'L', t0=1); say('  4. forged base point t0 = 1: NOT REFUSED'); ok_all = False
    except AssertionError as ex: say('  4. forged base point t0 = 1: refused (%s), as it must be' % ex)
    # 5. with the pairing left out nothing may be killed (the classes are in Sel^2 and not images of points found)
    left2, kill2, res2, _ = analyse_ratio(m, n, pairs=[('0', '8')], forge='no pairing')
    nk = all(res2[t]['08'] in ('image', 'kernel', 'undecided') for t in left2)
    ok_all &= nk
    say('  5. pairing left out, (465,416) on E_08: %s' % ('nothing killed, as it must be' if nk else 'SOMETHING KILLED'))
    say('ALL NEGATIVE TESTS AND CONTROLS AS THEY MUST BE' if ok_all else 'A TEST FAILED')

BADC = {}
def summarize(files):
    """the per-tuple lines of one or more outputs of --list, counted again"""
    import re
    tot = dict(ratios=0, tuples=0, pk=0, good=0, both=0, good_not_pk=0, pk_not_good=0, neither=0)
    perpair = {}; surv = []; ratio = None; rall = {}; places = {}
    for f in files:
        for l in open(f):
            if l.startswith('ratio ('):
                ratio = l.split()[1].rstrip(':'); tot['ratios'] += 1; rall[ratio] = True
            mm = re.match(r'    (\([0-9, ]+\))  good prime (\S+)  (.*)  => (.*)$', l.rstrip('\n'))
            if mm:
                t, g, prs, verdict = mm.groups()
                tot['tuples'] += 1
                ks = [x.split(':')[0] for x in prs.split() if x.split(':', 1)[1] in ('killed', 'outside')]
                for jk in ks: perpair[jk] = perpair.get(jk, 0) + 1
                gk = g != 'None'
                tot['pk'] += bool(ks); tot['good'] += gk; tot['both'] += bool(ks) and gk
                tot['good_not_pk'] += gk and not ks; tot['pk_not_good'] += bool(ks) and not gk; tot['neither'] += not gk and not ks
                if not ks: surv.append((ratio, t, g, prs)); rall[ratio] = False
            mm = re.match(r'    tuple .* local -1 at \[(.*)\]', l)
            if mm:
                m_, n_ = map(int, ratio.strip('()').split(','))
                if ratio not in BADC: BADC[ratio] = set(badprimes(m_, n_)) | {2}
                for p in mm.group(1).split(','):
                    p = p.strip().strip("'")
                    kind = 'inf' if p == 'inf' else ('bad prime of the ratio' if int(p) in BADC[ratio] else 'OTHER')
                    places[kind] = places.get(kind, 0) + 1
    say('combined summary of %s' % ', '.join(os.path.basename(f) for f in files))
    say('  %d ratios, %d residual tuples; killed by a good prime: %d; pairing-killed: %d' % (tot['ratios'], tot['tuples'], tot['good'], tot['pk']))
    say('  good-prime-killed and pairing-killed: %d; good-prime-killed, NOT pairing-killed: %d; pairing-killed only: %d; neither: %d'
        % (tot['both'], tot['good_not_pk'], tot['pk_not_good'], tot['neither']))
    say('  ratios all of whose residual tuples are pairing-killed: %d of %d' % (sum(rall.values()), len(rall)))
    say('  kills per pair curve: %s' % ' '.join('%s:%d' % (a + b, perpair.get(a + b, 0)) for a, b in PAIRS))
    say('  places carrying a -1 in the kill lines (every pair, one b each): %s' % places)
    say('  ratios with a tuple not pairing-killed: %s' % sorted(set(r for r, *_ in surv)))
    for r, t, g, prs in surv: say('    %s %s good prime %s: %s' % (r, t, g, prs))
    return tot, places

def save_hints():
    HINTS['rs_conic'] = RS.HINTS['conic']; HINTS['ri_conic'] = RI.HINTS.get('conic', {})
    json.dump(HINTS, open(HFILE, 'w'), indent=0, sort_keys=True)

if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    if '--summarize' in sys.argv:
        files = [os.path.join(HERE, x) for x in _arg('--summarize').split(',')]
        summarize(files)
        if WRITE: open(os.path.join(HERE, 't10_summary.txt'), 'w').write('\n'.join(OUT) + '\n')
    elif '--deep' in sys.argv:
        for a in _arg('--deep').split(';'):
            m, n = map(int, a.split(',')); deep(m, n, int(_arg('--effort', '2')))
        if WRITE: open(os.path.join(HERE, 't10_deep.txt'), 'w').write('\n'.join(OUT) + '\n')
    elif '--h1control' in sys.argv:
        h1_control()
        if WRITE: open(os.path.join(HERE, 't10_h1control.txt'), 'w').write('\n'.join(OUT) + '\n')
    elif '--negative' in sys.argv:
        negative()
        f = os.path.join(HERE, 't10_negative.txt')
        if WRITE: open(f, 'w').write('\n'.join(OUT) + '\n'); save_hints()
    elif '--list' in sys.argv:
        path = _arg('--list'); pairs = read_list(path)
        if '--sample' in sys.argv:
            rnd = random.Random(int(_arg('--seed', '1'))); pairs = sorted(rnd.sample(pairs, int(_arg('--sample'))))
        lo, hi = int(_arg('--from', '0')), int(_arg('--to', str(len(pairs))))
        pairs = pairs[lo:hi]
        name = _arg('--out', os.path.basename(path).split('.')[0])
        run_list(pairs, name)
        if WRITE:
            open(os.path.join(HERE, 't10_%s.txt' % name), 'w').write('\n'.join(OUT) + '\n'); save_hints()
    else:
        for a in args:
            m, n = map(int, a.split(','))
            analyse_ratio(m, n)
