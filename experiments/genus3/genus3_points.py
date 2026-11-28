#!/usr/bin/env python3
"""Claim (B) of the separate session's note (27 November 2026, section 2): none of the 27 curves E1, E2, E3 of the nine
open ratios has rank 0.  Here: for each curve a rational point of infinite order, checked exactly; for the curves where
PARI finds two points, independence in E(Q)/2E(Q), exactly; the curves pairwise non-isogenous and non-CM, exactly.

Written by an AI model (build agent t12); not read by a person.  Experiment, tier D.  Nothing here changes a count.

    python3 -I genus3_points.py [--write]       genus3_points.txt
    python3 -I genus3_points.py --negative      forged points and models that must be refused

Models (own, checked symbolically in sympy for general m, n): with A = m^2 n^2, B = m^4 + n^4, C = S^2, D = 2S^2 - 16m^2n^2,
  E2: X = ACz, W = ACw:  W^2 = (X + 2AC)(X + BC)(X + AD);   E3: W^2 = (X - 2AC)(X + BC)(X + AD);
  E1: the quartic w^2 = (z^2 - 4)(Az + B)(Cz + D) with z = 2 sent to infinity, s = 1/(z - 2), v = w s^2, X = c s,
      W = c v, c = 16 S^2 (m^2 - n^2)^2:  W^2 = (X + 4S^2(m^2-n^2)^2)(X + 16A(m^2-n^2)^2)(X + 4S^4).
Then each model y^2 = (x - e1)(x - e2)(x - e3) is moved to e1 = 0 and divided by the largest square u^2 dividing e2, e3
(x -> (x - e1)/u^2, y -> y/u^3), an isomorphism over Q.
What rests on PARI (cypari2): only the search for the points (ellrank, used as a finder; its bounds are printed as
PARI's, not used).  Exact, in integers and fractions: each point is on the curve; n P != O for n = 1, ..., 12, so P has
infinite order (Mazur, quoted: a rational torsion point has order <= 12); also some multiple k P (k <= 8) has a
non-integral x, so P is not torsion by Nagell-Lutz (quoted) as a second check.  Two points P1, P2: the 2-descent map
delta(P) = (x - e1, x - e2) mod squares (the standard fix at e1, e2) is injective on E(Q)/2E(Q) for full rational
2-torsion (quoted); if delta(T1), delta(T2), delta(P1), delta(P2) are independent over F_2, E(Q)/2E(Q) has dimension
>= 4 = rank + 2, so rank >= 2.  Non-isogeny of E1, E2, E3 at a ratio: a_p differs at one good prime (isogenous curves
over Q have equal a_p; quoted); a_p by counting points in pure Python.  Non-CM: j is not an integer (CM j-invariants
are algebraic integers; quoted).
"""
import sys, os, math, json
from fractions import Fraction as Fr
from math import isqrt, gcd
import sympy as sp

WRITE = '--write' in sys.argv; NEG = '--negative' in sys.argv
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = []
def say(s=''):
    OUT.append(s); print(s, flush=True)

OPEN = [(195, 28), (223, 167), (440, 279), (609, 200), (627, 364), (719, 161), (809, 409), (861, 620), (991, 263)]

# ---- models, checked symbolically -------------------------------------------------------------------------
def check_models():
    m, n, z, w = sp.symbols('m n z w')
    S = m**2 + n**2; A, B, C, D = m**2 * n**2, m**4 + n**4, S**2, 2 * S**2 - 16 * m**2 * n**2
    Pz = {1: (z**2 - 4) * (A * z + B) * (C * z + D), 2: (z + 2) * (A * z + B) * (C * z + D), 3: (z - 2) * (A * z + B) * (C * z + D)}
    X = A * C * z
    ok2 = sp.expand((X + 2 * A * C) * (X + B * C) * (X + A * D) - (A * C)**2 * Pz[2]) == 0
    ok3 = sp.expand((X - 2 * A * C) * (X + B * C) * (X + A * D) - (A * C)**2 * Pz[3]) == 0
    c = 16 * S**2 * (m**2 - n**2)**2
    s_ = 1 / (z - 2); X1 = c * s_
    # W = c v = c w s^2: W^2 - prod(X - e_i) must be c^2 s^4 (w^2 - P1(z))
    e = [-4 * S**2 * (m**2 - n**2)**2, -16 * A * (m**2 - n**2)**2, -4 * S**4]
    lhs = (X1 - e[0]) * (X1 - e[1]) * (X1 - e[2])
    ok1 = sp.cancel(sp.together(lhs - c**2 * s_**4 * Pz[1])) == 0
    return ok1 and ok2 and ok3

def raw_roots(M, N, i):
    S = M * M + N * N; A, B, C, D = M * M * N * N, M**4 + N**4, S * S, 2 * S * S - 16 * M * M * N * N
    if i == 1: return (-4 * S * S * (M * M - N * N)**2, -16 * A * (M * M - N * N)**2, -4 * S**4)
    if i == 2: return (-2 * A * C, -B * C, -A * D)
    return (2 * A * C, -B * C, -A * D)

def square_part(g):
    """largest u with u^2 | g (g > 0), by sympy's factorint (the result is checked: u^2 divides g)"""
    u = 1
    for p, k in sp.factorint(g).items(): u *= p ** (k // 2)
    assert g % (u * u) == 0
    return u

def model(M, N, i):
    e = sorted(raw_roots(M, N, i))
    e1 = e[0]; d2, d3 = e[1] - e1, e[2] - e1
    u = square_part(gcd(d2, d3))
    return (0, d2 // (u * u), d3 // (u * u))

# ---- exact group law on y^2 = (x - e1)(x - e2)(x - e3) ------------------------------------------------------
def coeffs(es):
    e1, e2, e3 = es
    return (-(e1 + e2 + e3), e1 * e2 + e1 * e3 + e2 * e3, -e1 * e2 * e3)   # a2, a4, a6
def on_curve(P, es):
    if P is None: return True
    x, y = P
    return y * y == (x - es[0]) * (x - es[1]) * (x - es[2])
def add(P, Q, es):
    a2, a4, a6 = coeffs(es)
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if y1 + y2 == 0: return None
        lam = (3 * x1 * x1 + 2 * a2 * x1 + a4) / (2 * y1)
    else:
        lam = (y2 - y1) / (x2 - x1)
    x3 = lam * lam - a2 - x1 - x2
    return (x3, -(y1 + lam * (x3 - x1)))
def mult(k, P, es):
    R = None
    for _ in range(k): R = add(R, P, es)
    return R

def infinite_order(P, es):
    """exact: P on the curve, n P != O for n = 1..12 (Mazur); and a multiple k P, k <= 8, with non-integral x
    (Nagell-Lutz on this integral model).  Returns (ok, k)."""
    assert on_curve(P, es), 'the point is not on the curve'
    R = None
    for k in range(1, 13):
        R = add(R, P, es)
        assert R is not None, 'n P = O for some n <= 12: the point is torsion'
        assert on_curve(R, es)
    R = None; kk = None
    for k in range(1, 9):
        R = add(R, P, es)
        if R[0].denominator != 1: kk = k; break
    assert kk is not None, 'no multiple k P (k <= 8) with non-integral x'
    return kk

def issq(q):
    q = Fr(q)
    return q > 0 and isqrt(q.numerator) ** 2 == q.numerator and isqrt(q.denominator) ** 2 == q.denominator
def delta(P, es):
    e1, e2, e3 = es
    x = P[0]
    if x == e1: return (Fr((e1 - e2) * (e1 - e3)), Fr(e1 - e2))
    if x == e2: return (Fr(e2 - e1), Fr((e2 - e1) * (e2 - e3)))
    return (Fr(x - e1), Fr(x - e2))
def f2_rank(vecs):
    """rank over F_2 of classes in (Q*/Q*^2)^2: k minus log2 of the number of subsets with square product"""
    k = len(vecs); sq = 0
    for msk in range(1 << k):
        a = Fr(1); b = Fr(1)
        for j in range(k):
            if msk >> j & 1: a *= vecs[j][0]; b *= vecs[j][1]
        if issq(a) and issq(b): sq += 1
    l = sq.bit_length() - 1; assert 1 << l == sq
    return k - l

# ---- a_p by counting, j ----------------------------------------------------------------------------------
def ap(es, p):
    e1, e2, e3 = [e % p for e in es]
    N = 1
    for x in range(p):
        v = (x - e1) * (x - e2) * (x - e3) % p
        N += 1 if v == 0 else (2 if pow(v, (p - 1) // 2, p) == 1 else 0)
    return p + 1 - N
def jinv(es):
    a2, a4, a6 = coeffs(es)
    b2, b4, b6 = 4 * a2, 2 * a4, 4 * a6
    b8 = 4 * a2 * a6 - a4 * a4          # a1 = a3 = 0
    c4 = b2 * b2 - 24 * b4
    Dl = -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6
    return Fr(c4 ** 3, Dl), Dl

_PARI = []
def gp():
    if not _PARI:
        import cypari2
        p = cypari2.Pari(); p.allocatemem(2 * 10**9); _PARI.append(p)
    return _PARI[0]
HINTS_FILE = os.path.join(HERE, 'genus3_points_hints.json')
HINTS = json.load(open(HINTS_FILE)) if os.path.exists(HINTS_FILE) else {}
def frp(P): return (Fr(int(gp().numerator(P[0])), int(gp().denominator(P[0]))), Fr(int(gp().numerator(P[1])), int(gp().denominator(P[1]))))
def pari_E(es):
    a2, a4, a6 = coeffs(es)
    return gp().ellinit([0, a2, 0, a4, a6])
def pari_interval(es):
    """PARI's ellrank at effort 0 (deterministic): a guide, printed as PARI's"""
    E = pari_E(es); r = gp().ellrank(E, 0)
    return (int(r[0]), int(r[1]), int(r[2])), int(gp().ellrootno(E))
def find_points(key, es):
    """points of infinite order: from the hints file if present (each is checked exactly by the caller), else
    PARI's ellrank with increasing effort (randomised, so the points found are kept as hints)"""
    if key in HINTS: return [(Fr(x), Fr(y)) for x, y in HINTS[key]['points']]
    E = pari_E(es); gp().setrand(1)
    for eff in (0, 1, 2, 4, 6, 8, 10):
        r = gp().ellrank(E, eff)
        if len(r[3]): break
    pts = [frp(P) for P in r[3]]
    HINTS.setdefault(key, {})['points'] = [[str(P[0]), str(P[1])] for P in pts]
    return pts
def find_more_points(key, es):
    """for root number +1: more points, from PARI's ellrank on every curve of the isogeny class (ellisomat), mapped
    back to this curve by PARI's dual isogeny (ellisogenyapply); each image is checked exactly by the caller (on the
    curve, infinite order); PARI is only a finder here"""
    if key in HINTS and 'more' in HINTS[key]: return [(Fr(x), Fr(y)) for x, y in HINTS[key]['more']]
    pari = gp(); pari.setrand(7)
    E = pari_E(es); L = pari.ellisomat(E, 0, 0)[0]
    out = []
    for k in range(len(L)):
        Ej = pari.ellinit(L[k][0]); dual = L[k][2]
        for eff in (0, 2, 4):
            for P in pari.ellrank(Ej, eff)[3]:
                Q = pari.ellisogenyapply(dual, P)
                if len(Q) < 2: continue
                Q = frp(Q)
                if on_curve(Q, es) and Q not in out and order(Q, es) is None: out.append(Q)
    HINTS.setdefault(key, {})['more'] = [[str(P[0]), str(P[1])] for P in out]
    return out
def indep_mod3(pts, Ts, es, nq=60):
    """exact F_3-rank of the images of the points modulo the torsion, through the maps E(Q) -> E(F_q)[3], P -> (N/3) P
    (q a good prime with 3 || N = #E(F_q), N counted; a homomorphism that kills 3E(Q)); a primitive relation
    a P1 + b P2 + ... = T (T torsion) would give the same F_3-relation in every row, so rank(points, torsion) -
    rank(torsion) is a lower bound for the rank of the subgroup generated by the points modulo torsion"""
    disc = (es[0] - es[1]) * (es[0] - es[2]) * (es[1] - es[2])
    a2, a4, a6 = coeffs(es)
    def addq(P, Q, q):
        if P is None: return Q
        if Q is None: return P
        x1, y1 = P; x2, y2 = Q
        if x1 == x2:
            if (y1 + y2) % q == 0: return None
            lam = (3 * x1 * x1 + 2 * a2 * x1 + a4) * pow(2 * y1, -1, q) % q
        else: lam = (y2 - y1) * pow(x2 - x1, -1, q) % q
        x3 = (lam * lam - a2 - x1 - x2) % q
        return (x3, (-(y1 + lam * (x3 - x1))) % q)
    def mulq(k, P, q):
        R = None
        while k:
            if k & 1: R = addq(R, P, q)
            P = addq(P, P, q); k >>= 1
        return R
    cols = list(pts) + list(Ts)
    rows = []
    for q in sp.primerange(5, 10**5):
        if len(rows) >= nq: break
        if disc % q == 0 or any(P[0].denominator % q == 0 or P[1].denominator % q == 0 for P in cols): continue
        N = q + 1 - ap(es, q)
        if N % 3 or N % 9 == 0: continue
        red = [(P[0].numerator * pow(P[0].denominator, -1, q) % q, P[1].numerator * pow(P[1].denominator, -1, q) % q) for P in cols]
        img = [mulq(N // 3, R, q) for R in red]
        g = next((R for R in img if R is not None), None)
        if g is None: rows.append([0] * len(cols)); continue
        g2 = addq(g, g, q)
        rows.append([0 if R is None else (1 if R == g else (2 if R == g2 else None)) for R in img])
        assert None not in rows[-1], 'an image outside the group of order 3'
    def rk3(M):
        M = [r[:] for r in M]; r_ = 0; ncol = len(M[0]) if M else 0
        for c in range(ncol):
            piv = next((i for i in range(r_, len(M)) if M[i][c] % 3), None)
            if piv is None: continue
            M[r_], M[piv] = M[piv], M[r_]
            inv = 1 if M[r_][c] % 3 == 1 else 2
            M[r_] = [(v * inv) % 3 for v in M[r_]]
            for i in range(len(M)):
                if i != r_ and M[i][c] % 3:
                    f_ = M[i][c]; M[i] = [(M[i][j] - f_ * M[r_][j]) % 3 for j in range(ncol)]
            r_ += 1
        return r_
    full = rk3(rows); tors = rk3([r[len(pts):] for r in rows]) if Ts else 0
    return full - tors, len(rows)

def find_torsion(key, es):
    """the torsion points: PARI's elltors as a finder; checked exactly by the caller"""
    if key in HINTS and 'torsion' in HINTS[key]: return [(Fr(x), Fr(y)) for x, y in HINTS[key]['torsion']]
    T = gp().elltors(pari_E(es))
    # all points of the group generated by PARI's generators
    gens = [frp(P) for P in T[2]]
    pts = [None]
    for g in gens:
        new = []
        for P in pts:
            R = P
            while True:
                if R not in new: new.append(R)
                R = add(R, g, es)
                if R == P: break
        pts = new
    out = [P for P in pts if P is not None]
    HINTS.setdefault(key, {})['torsion'] = [[str(P[0]), str(P[1])] for P in out]
    return out
def order(P, es):
    R = P
    for k in range(1, 13):
        if R is None: return k
        R = add(R, P, es)
    return None
def torsion_exact(Ts, es):
    """exact: every listed point has finite order (k P = O, k <= 12), the list with O is closed under addition (a
    subgroup), and its order equals gcd #E(F_q) over good primes 17 <= q < 400 (torsion injects into E(F_q) for good q
    prime to its order; Mazur: order <= 16): so the list is the whole torsion subgroup.  Returns its structure."""
    G_ = [None] + list(Ts)
    for P in Ts:
        assert on_curve(P, es), 'a torsion point is not on the curve'
        assert order(P, es) is not None, 'a listed torsion point has infinite order'
    for P in G_:
        for Q in G_: assert add(P, Q, es) in G_, 'the torsion list is not closed under addition'
    disc = (es[0] - es[1]) * (es[0] - es[2]) * (es[1] - es[2])
    g = 0
    for q in sp.primerange(17, 400):
        if disc % q: g = gcd(g, q + 1 - ap(es, q))
    assert g == len(G_), 'the torsion list is incomplete (gcd #E(F_q) = %d, listed %d)' % (g, len(G_))
    o = max(order(P, es) for P in Ts)
    return 'Z/2 x Z/%d' % o if len(G_) == 2 * o else 'order %d' % len(G_)

def run():
    say('genus3_points.py: a point of infinite order on each of the 27 curves E1, E2, E3 of the nine open ratios')
    ok = check_models()
    say('the three Weierstrass models agree with the note\'s E1, E2, E3 (sympy, general m, n): %s' % ok)
    assert ok
    say()
    summary = []
    for (M, N) in OPEN:
        say('(%d,%d)' % (M, N))
        es_all = {}
        for i in (1, 2, 3):
            es = model(M, N, i); es_all[i] = es
            key = '%d,%d,E%d' % (M, N, i)
            pts = find_points(key, es)
            Ts = find_torsion(key, es)
            tstr = torsion_exact(Ts, es)
            (lo, hi, s), w = pari_interval(es)
            good = []
            for P in pts:
                kk = infinite_order(P, es); good.append((P, kk))
            assert good, 'no point'
            r2 = f2_rank([delta(P, es) for P in Ts] + [delta(P, es) for P, _ in good])
            exact_lo = r2 - 2
            assert exact_lo >= 1
            extra = ''
            if w == 1:
                more = find_more_points(key, es)
                allp = [P for P, _ in good] + [P for P in more if P not in [Q for Q, _ in good]]
                for P in allp: infinite_order(P, es)
                r3, nrow = indep_mod3(allp, Ts, es)
                extra = ('      %d point(s) in all (more from the isogeny class through PARI\'s dual isogenies, each checked'
                         ' exactly); images in E(F_q)[3] at %d primes q: F_3-rank %d modulo torsion: rank >= %d (exact)'
                         % (len(allp), nrow, r3, r3))
                exact_lo = max(exact_lo, r3)
            hP = max(len(str(good[0][0][0].numerator)), len(str(good[0][0][0].denominator)))
            say('  E%d: y^2 = x(x - %d)(x - %d); torsion %s (exact)' % (i, es[1], es[2], tstr))
            say('      %d point(s) (hints from PARI), each checked: on the curve, n P != O for n <= 12, non-integral x at %s;'
                % (len(good), ', '.join('%dP' % kk for _, kk in good)))
            say('      first point x = %s' % (good[0][0][0] if hP <= 60 else '(a fraction with %d digits)' % hP))
            say('      delta of the torsion and of the points spans dimension %d of E(Q)/2E(Q): rank >= %d (exact)' % (r2, exact_lo))
            if extra: say(extra)
            say('      PARI\'s ellrank at effort 0 (a guide): [%d, %d], s = %d; root number %+d (PARI)' % (lo, hi, s, w))
            summary.append(((M, N), i, exact_lo, (lo, hi), w))
        # non-isogeny and non-CM
        pr = [p for p in sp.primerange(5, 400) if all(math.prod(((x - y) for x in es for y in es if x < y)) % p for es in es_all.values())]
        diffs = []
        for i, j in ((1, 2), (1, 3), (2, 3)):
            q = next(p for p in pr if ap(es_all[i], p) != ap(es_all[j], p))
            diffs.append('E%d/E%d at p = %d (%d vs %d)' % (i, j, q, ap(es_all[i], q), ap(es_all[j], q)))
        say('  pairwise non-isogenous over Q, a_p differ: ' + '; '.join(diffs))
        jn = []
        for i in (1, 2, 3):
            j, Dl = jinv(es_all[i]); assert Dl != 0
            assert j.denominator != 1; jn.append('E%d: denominator of j has %d digits' % (i, len(str(j.denominator))))
        say('  non-CM (j not an integer): ' + '; '.join(jn))
    say()
    say('Summary (exact lower bounds; PARI\'s interval in brackets as PARI\'s only):')
    for (MN, i, lo, (a, b), w) in summary:
        say('  %-11s E%d  rank >= %d   [PARI %d..%d, w = %+d]' % ('(%d,%d)' % MN, i, lo, a, b, w))
    say('All 27 curves have a rational point of infinite order: no E_i has rank 0 (exact, given Mazur).')
    two = ['%s E%d' % ('(%d,%d)' % MN, i) for (MN, i, lo, _, _) in summary if lo >= 2]
    say('Rank >= 2 proved exactly for: %s' % (', '.join(two) if two else 'none (on the nine curves of root number +1 every point'
        ' found, also through the isogeny class, is dependent on the first modulo torsion: F_3-rank 1)'))
    return summary

def negative():
    res = []
    def t(name, fn):
        try: fn(); res.append(False); say('  NOT refused: ' + name)
        except AssertionError as ex: res.append(True); say('  refused: %s (%s)' % (name, ex.args[0] if ex.args else ''))
    es = model(627, 364, 3)
    P = find_points('627,364,E3', es)[0]
    t('a point moved off the curve (y + 1)', lambda: infinite_order((P[0], P[1] + 1), es))
    t('a point of order 2 offered as a point of infinite order', lambda: infinite_order((Fr(es[1]), Fr(0)), es))
    def tors4():
        # a point of order 4 on the curve, if any: halve a 2-torsion point; else a point of finite order from PARI
        Ts = find_torsion('627,364,E3', es)
        Q = max(Ts, key=lambda R: order(R, es))     # a point of maximal finite order
        infinite_order(Q, es)
    t('a torsion generator (from PARI elltors) offered as a point of infinite order', tors4)
    t('the same point on the wrong curve (E2 of the same ratio)', lambda: infinite_order(P, model(627, 364, 2)))
    def dep():
        # two dependent points (P and 3P): the F_2 rank must not reach 4
        Q = mult(3, P, es); T = [(Fr(es[0]), Fr(0)), (Fr(es[1]), Fr(0))]
        r2 = f2_rank([delta(R, es) for R in T] + [delta(P, es), delta(Q, es)])
        assert r2 >= 4, 'dependent points give F_2 rank %d < 4' % r2
    t('P and 3P offered as independent points', dep)
    def dep3():
        Ts = find_torsion('627,364,E3', es)
        r3, _ = indep_mod3([P, mult(2, P, es), add(P, Ts[0], es)], Ts, es)
        assert r3 >= 2, 'P, 2P, P + T give F_3-rank %d < 2' % r3
    t('P, 2P and P + T offered to the mod-3 independence test', dep3)
    def badmodel():
        es2 = (es[0], es[1] + 1, es[2])
        infinite_order(P, es2)
    t('an altered model (e2 + 1)', badmodel)
    def tors_incomplete():
        Ts = find_torsion('627,364,E3', es)
        torsion_exact([Q for Q in Ts if Q[1] == 0], es)
    t('the torsion list cut to the points of order 2', tors_incomplete)
    def tors_forged():
        Ts = find_torsion('627,364,E3', es)
        torsion_exact(Ts + [P], es)
    t('a point of infinite order added to the torsion list', tors_forged)
    say('negative tests: %d of %d refused' % (sum(res), len(res)))
    return all(res)

if __name__ == '__main__':
    if NEG:
        say('genus3_points.py --negative')
        ok = negative()
        if WRITE: open(os.path.join(HERE, 'genus3_points_negative.txt'), 'w').write('\n'.join(OUT) + '\n')
        sys.exit(0 if ok else 1)
    run()
    if WRITE:
        open(os.path.join(HERE, 'genus3_points.txt'), 'w').write('\n'.join(OUT) + '\n')
        json.dump(HINTS, open(HINTS_FILE, 'w'), indent=0, sort_keys=True)
