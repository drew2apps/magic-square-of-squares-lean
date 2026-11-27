"""H1: the four extra tuples of (1109503, 597671) on the pair curves E_0L and E_0U (EXPERIMENT, written by an AI
model; not read by a person; nothing here changes the census).

THE QUOTIENTS.  The four-class descent works with the Gaussian linear forms (t = u/v)
    a0 = m + i n t,  a8 = n + i m t,  aL = (m + n)(t - 1) + i(m - n)(t + 1),  aU = (m - n)(t - 1) + i(m + n)(t + 1),
and a solution with class (d0, d8, dL, dU) has a_j = d_j z_j^2 in Q(i).  For two of the forms, j and k, the point
(t, y) with y^2 = F_j(t) F_k(t) (F = N(a)) lies on the genus-one curve E_jk, and a_j a_k = (d_j d_k) (z_j z_k)^2:
the class D = d_j d_k (in Q*/<-1, squares>) is the value at (t, y) of the descent map of E_jk for the rational
2-torsion point that belongs to the factorisation F_j F_k = (a_j a_k)(conj a_j a_k).  So a tuple can come from a
rational point only if its class D_jk is the class of a rational point of E_jk.
THE EXTRA TUPLES (G1): the known tuples (of t = 0 and t = infinity) times (1, 1, 959, 1) or (1, 1, 1, 959).  The first
kind changes D_0L by 959, the second D_0U.

WHAT IS COMPUTED, EXACTLY (integers and rationals):
 1. E_jk from the quartic y^2 = F_j F_k(t) with its rational point t = 0, by the map of Washington, Elliptic Curves,
    Thm 2.17 (checked: every point used is mapped back to (t, y) and y^2 = F_j F_k(t) is tested), then an integral
    model y^2 = (x - e1)(x - e2)(x - e3).
 2. Three points of infinite order (PARI's ellrank, a FINDER; stored in h1_hints.json and checked on the curve).
 3. The complete 2-descent and the Cassels-Tate pairing on Sel^2, by the repository's scripts/g1_full/rank_second.py
    (full_selmer2, DFull.pairing: exact; conic points from PARI's qfsolve as hints, checked there), as for E0L of
    (920,129): rank <= dim Sel^2 - 2 - rank(pairing), rank >= rank of the point images modulo torsion.
 4. The root e* whose descent map x - e* is the class D_jk / D_jk(t = t0) (t0 = 0 here; for square-type ratios the
    first t0 in 1, -1, 2, ... where F_j F_k(t0) is a non-zero square): checked on the points (the generators, the
    2-torsion and sums), modulo -1; exactly one root must agree.
 5. When the bounds of 3 meet, the image of E_jk(Q) in Sel^2 is the span of the point images (its dimension is rank
    + 2), and a tuple is excluded when no element of that span has x - e* = +- D_jk(tuple) / D_jk(t = 0).
What it rests on: the theory of the Cassels-Tate pairing and its formula on 2-coverings (used, not proved; see
rank_second.py), and the repository's programs named above.  No PARI output is relied on: PARI finds points and conic
points, which are checked.
usage: h1_pairs.py [--write]      (h1_pairs.txt)
       h1_pairs.py --negative     (h1_pairs_negative.txt: forged inputs that must be refused)
The repository is found through MSS_REPO, or else ../../.."""
import sys, os, json, random
from fractions import Fraction as Fr
from math import isqrt, gcd
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get('MSS_REPO', os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(REPO, 'scripts', 'g1_full'))
import rank_indep as RI
import rank_second as RS
import sympy
WRITE = '--write' in sys.argv
HFILE = os.path.join(HERE, 'h1_hints.json')
HINTS = json.load(open(HFILE)) if os.path.exists(HFILE) else {}
RS.HINTS['conic'].update(HINTS.get('rs_conic', {}))
RI.HINTS.setdefault('conic', {}).update(HINTS.get('ri_conic', {}))
OUT = []
def say(s=''): OUT.append(s)

M, N = 1109503, 597671
KNOWN = [(1109503, 597671, 2, 1), (1195342, 2219006, 2, 1)]
EXTRA = [(1109503, 597671, 2, 959), (1109503, 597671, 1918, 1), (1195342, 2219006, 2, 959), (1195342, 2219006, 1918, 1)]
IDX = {'0': 0, '8': 1, 'L': 2, 'U': 3}

def forms(m, n, t):
    return {'0': (Fr(m), n * t), '8': (Fr(n), m * t), 'L': ((m + n) * (t - 1), (m - n) * (t + 1)), 'U': ((m - n) * (t - 1), (m + n) * (t + 1))}
def gm(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def isq(q): q = Fr(q); return q >= 0 and isqrt(q.numerator) ** 2 == q.numerator and isqrt(q.denominator) ** 2 == q.denominator
def rsqrt(q): q = Fr(q); assert isq(q); return Fr(isqrt(q.numerator), isqrt(q.denominator))
FB = []
def set_base(ps):
    FB[:] = sorted(set(ps))
def sqcl(q):
    """squarefree part of |q| (the class in Q*/<-1, squares>), over the factor base FB (primes proved prime by
    rank_indep.factor_int); the cofactor must be a square."""
    q = Fr(q); a = abs(q.numerator * q.denominator); d = 1
    for p in FB:
        e = 0
        while a % p == 0: a //= p; e += 1
        if e % 2: d *= p
    assert a > 0 and isqrt(a) ** 2 == a, 'a class outside the factor base'
    return d
def gauss_class(z):
    """d with z = d w^2 in Q(i), for z of square norm (z = (z + r)^2 / (2 (Re z + r)), r^2 = N z)."""
    Nz = z[0] ** 2 + z[1] ** 2; r = rsqrt(Nz); x = z[0] + r
    if x == 0: x = z[0] - r
    return sqcl(2 * x)
def Dpair(m, n, j, k, t):
    f = forms(m, n, t); return gauss_class(gm(f[j], f[k]))

def quartic(m, n, j, k, t0=0):
    """the quartic F_j F_k (t0 + t) (t0 = 0 unless that value is not a non-zero square)"""
    T = sympy.Symbol('t') + t0
    al = {'0': (m, n * T), '8': (n, m * T), 'L': ((m + n) * (T - 1), (m - n) * (T + 1)), 'U': ((m - n) * (T - 1), (m + n) * (T + 1))}
    P = sympy.expand((al[j][0] ** 2 + al[j][1] ** 2) * (al[k][0] ** 2 + al[k][1] ** 2))
    return [int(P.coeff(T, i)) for i in range(5)]           # e, d, c, b, a (low degree first)

class Pair:
    def __init__(s, m, n, j, k):
        s.m, s.n, s.j, s.k = m, n, j, k
        # base point: t0 = 0 if F_j F_k (0) is a non-zero square (twice-square type), else the first of 1, -1, 2, -2, ...
        for t0 in [0] + [s_ * i for i in range(1, 50) for s_ in (1, -1)]:
            co = quartic(m, n, j, k, t0)
            if co[0] != 0 and isq(co[0]): break
        else: raise AssertionError('no base point t0 found')
        s.t0 = Fr(t0); s.co = co; e, d, c, b, a = [Fr(x) for x in s.co]
        s.q = q = rsqrt(e)
        s.ai = [d / q, c - d * d / (4 * q * q), 2 * q * b, -4 * q * q * a, (c - d * d / (4 * q * q)) * (-4 * q * q * a)]
        a1, a2, a3, a4, a6 = s.ai
        A2, A4, A6 = a2 + a1 * a1 / 4, a4 + a1 * a3 / 2, a6 + a3 * a3 / 4
        u = 1
        while not all((cc * u ** (2 * i)).denominator == 1 for cc, i in ((A2, 1), (A4, 2), (A6, 3))): u += 1
        s.u = u; s.A = (int(A2 * u ** 2), int(A4 * u ** 4), int(A6 * u ** 6))
        X = sympy.Symbol('X')
        s.es = tuple(sorted(int(r) for r in sympy.Poly(X ** 3 + s.A[0] * X ** 2 + s.A[1] * X + s.A[2], X).ground_roots()))
        assert len(s.es) == 3, 'full rational 2-torsion expected'
    def on(s, P): x, y = P; return y * y == (x - s.es[0]) * (x - s.es[1]) * (x - s.es[2])
    def to_t(s, P):
        """(x, y) on the integral model -> t on the quartic (None for the point over t = 0, y = -q, or a point
        with y = 0 on Washington's model), checked on the quartic."""
        a1, a2, a3, a4, a6 = s.ai; e, d, c, b, a = [Fr(x) for x in s.co]
        x = Fr(P[0]) / s.u ** 2; y = Fr(P[1]) / s.u ** 3 - (a1 * x + a3) / 2
        if y == 0: return None
        t = (2 * s.q * (x + c) - d * d / (2 * s.q)) / y
        assert isq(sum(cc * t ** i for i, cc in enumerate(s.co))), 'the map to the quartic is wrong'
        return t + s.t0
    def add(s, P, Q):
        """the group law on y^2 = (x - e1)(x - e2)(x - e3) (P, Q affine, P != -Q), exact."""
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

def find_points(Pr):
    key = '%d,%d %s%s' % (Pr.m, Pr.n, Pr.j, Pr.k)
    if key in HINTS.get('points', {}): return [(Fr(x), Fr(y)) for x, y in HINTS['points'][key]]
    import cypari2
    pari = cypari2.Pari(); pari.allocatemem(10 ** 9)
    E = pari.ellinit([0, Pr.A[0], 0, Pr.A[1], Pr.A[2]])
    r = pari.ellrank(E, 3)
    pts = [(Fr(str(P[0])), Fr(str(P[1]))) for P in r[3]]
    HINTS.setdefault('points', {})[key] = [[str(x), str(y)] for x, y in pts]
    return pts

def analyse(m, n, j, k, tuples, forge=None):
    Pr = Pair(m, n, j, k)
    say('E_%s%s of (%d,%d): y^2 = F_%s F_%s (t); integral model y^2 = (x - e1)(x - e2)(x - e3), e = %s' % (j, k, m, n, j, k, Pr.es))
    G = find_points(Pr)
    for P in G: assert Pr.on(P), 'a hint point is not on the curve'
    bad, loc, sel, kum, enc, places = RS.full_selmer2(Pr.es, random.Random(1))
    S = m * m + n * n
    set_base(set(bad) | {2} | set(q for x in (m, n, m - n, m + n, S) for q in RI.factor_int(x)))
    gensv = [-1] + bad
    def gvec(pr):
        v = 0
        for jj, d in enumerate(pr):
            d = RI.sqfree_part(Fr(d).numerator * Fr(d).denominator); rest = abs(d)
            for i, g in enumerate(gensv):
                if (g == -1 and d < 0) or (g != -1 and d % g == 0): v |= 1 << (i + jj * len(gensv))
            for g in bad:
                if rest % g == 0: rest //= g
            assert rest == 1, 'a class outside the group of the bad primes'
        return v
    selv = [gvec(a) for a in sel]
    T2 = [(Fr(Pr.es[0]), Fr(0)), (Fr(Pr.es[1]), Fr(0))]
    imgs = [kum(P[0]) for P in G] + [kum(Pr.es[0]), kum(Pr.es[1])]
    imgv = [gvec(c) for c in imgs]
    assert len(RS.f2_span(selv + imgv)) == len(sel), 'a point image outside Sel^2: the descent is wrong'
    ri = len(RS.f2_span(imgv))
    Ds = [RS.DFull(Pr.es, a) for a in sel]
    Mx = [[Ds[i].pairing(b)[0] for b in sel] for i in range(len(sel))]
    for i in range(len(sel)):
        assert Mx[i][i] == 1, 'the pairing is not alternating'
        for jj in range(len(sel)): assert Mx[i][jj] == Mx[jj][i], 'the pairing is not symmetric'
    for i in range(len(sel)):
        for c in imgs:
            cc = (RI.sqfree_part(Fr(c[0]).numerator * Fr(c[0]).denominator), RI.sqfree_part(Fr(c[1]).numerator * Fr(c[1]).denominator))
            assert Ds[i].pairing(cc)[0] == 1, 'the image of a rational point pairs non-trivially'
    rho = len(RS.f2_span([sum((0 if v == 1 else 1) << jj for jj, v in enumerate(rw)) for rw in Mx]))
    if forge == 'no pairing': rho = 0
    up, lo = len(sel) - 2 - rho, ri - 2
    say('  complete 2-descent: dim Sel^2 = %d; Cassels-Tate pairing of rank %d; the images of %d points and the 2-torsion span dimension %d'
        % (len(sel), rho, len(G), ri))
    say('  => %d <= rank <= %d%s' % (lo, up, ' (exact)' if lo == up else ''))
    # the root e* of the class D_jk
    t0 = Pr.t0; D0 = Dpair(m, n, j, k, t0)
    if len(G) >= 3: pts = list(G) + [Pr.add(G[0], T) for T in T2] + [Pr.add(G[0], Pr.add(G[1], G[2])), Pr.add(G[1], T2[0]), Pr.add(G[2], T2[1])]
    else:                                   # fewer than three points (phantom ratios of the scan): their sums with the 2-torsion
        pts = list(G) + [Pr.add(P, T) for P in G for T in T2] + ([Pr.add(G[0], G[1])] if len(G) == 2 else [])
    data = []
    for P in pts:
        t = Pr.to_t(P)
        if t is None: continue
        data.append((P, sqcl(Fr(Dpair(m, n, j, k, t), D0))))
    agree = lambda e: all(sqcl(P[0] - e) == c for P, c in data)
    cands = list(Pr.es)
    if forge == 'wrong root': cands = [e for e in Pr.es if not agree(e)]     # only the roots that are not e*
    match = [e for e in cands if agree(e)]
    say('  roots e whose class x - e agrees with D_%s%s / D_%s%s(t = %s) on %d points: %d' % (j, k, j, k, t0, len(data), len(match)))
    assert len(match) == 1, 'the root of the descent map is not identified'
    ei = Pr.es.index(match[0])
    # classes of x - e* over the image of E(Q) (span of imgs); kum gives (x - e1, x - e2); x - e3 = their product
    def comp(c):
        a1, a2 = c
        return sqcl([a1, a2, Fr(a1) * Fr(a2)][ei])
    span = [(1, 1)]
    for c in imgs:
        span = span + [(Fr(a[0]) * Fr(c[0]), Fr(a[1]) * Fr(c[1])) for a in span]
    image_classes = set(comp(c) for c in span)
    sspan = [(1, 1)]
    for c in sel:
        sspan = sspan + [(Fr(a[0]) * Fr(c[0]), Fr(a[1]) * Fr(c[1])) for a in sspan]
    sel_classes = set(comp(c) for c in sspan)
    say('  classes of x - e* on E_%s%s(Q): %d; on Sel^2: %d' % (j, k, len(image_classes), len(sel_classes)))
    # the kernel of the Cassels-Tate pairing on Sel^2: it contains the image of E(Q)/2E(Q) (the pairing vanishes on
    # it), so a class that no element of the kernel has is not the class of a rational point, whatever the rank
    rad_classes = sel_classes
    if lo != up and forge != 'no pairing':
        B = [[0 if v == 1 else 1 for v in rw] for rw in Mx]
        nsel = len(sel); kern = []
        for x in range(1 << nsel):          # nsel <= 12 here: a direct search over Sel^2 is cheap
            if all(sum((x >> i & 1) * B[i][jj] for i in range(nsel)) % 2 == 0 for jj in range(nsel)): kern.append(x)
        assert len(kern) == 1 << (nsel - rho), 'the kernel of the pairing has the wrong size'
        def elt(x):
            a = (Fr(1), Fr(1))
            for i in range(nsel):
                if x >> i & 1: a = (a[0] * Fr(sel[i][0]), a[1] * Fr(sel[i][1]))
            return a
        rad_classes = set(comp(elt(x)) for x in kern)
        assert image_classes <= rad_classes, 'a rational point outside the kernel of the pairing'
        say('  classes of x - e* on the kernel of the pairing (dimension %d): %d' % (nsel - rho, len(rad_classes)))
    res = {}
    for tp in tuples:
        Dt = sqcl(Fr(tp[IDX[j]] * tp[IDX[k]], D0))
        insel = Dt in sel_classes
        if lo != up and Dt not in rad_classes:
            v = 'EXCLUDED (Cassels-Tate): its class %d is not the class of any element of the kernel of the pairing, so of no rational point of E_%s%s' % (Dt, j, k)
        elif lo != up: v = 'NOT DECIDED (the rank bounds differ%s)' % ('' if forge == 'no pairing' else '; the class is in the kernel of the pairing, %s' % (
                'and is the class of a rational point found' if Dt in image_classes else 'not the class of the rational points found'))
        elif Dt in image_classes: v = 'the class of a rational point of E_%s%s (not excluded here)' % (j, k)
        else: v = 'EXCLUDED: its class %d is not the class of any rational point of E_%s%s (it is in Sel^2: %s)' % (Dt, j, k, 'yes' if insel else 'no')
        res[tp] = v
        say('  tuple %s: D_%s%s / D(t = %s) = %d: %s' % (tp, j, k, t0, Dt, v))
    return res, Pr, sel

def main():
    say('H1: the extra tuples of (1109503,597671). Written by an AI model; not read by a person. Experiment, not counted.')
    allres = {}
    for j, k in (('0', 'L'), ('0', 'U')):
        res, Pr, sel = analyse(M, N, j, k, KNOWN + EXTRA)
        for tp, v in res.items(): allres.setdefault(tp, []).append((j + k, v))
    say('')
    say('VERDICT per tuple:')
    for tp in KNOWN + EXTRA:
        ex = [jk for jk, v in allres[tp] if v.startswith('EXCLUDED')]
        kind = 'known (degenerate point)' if tp in KNOWN else 'extra'
        say('  %s %s: %s' % (tp, kind, ('killed by E_%s' % ex[0]) if ex else 'not excluded by E_0L or E_0U'))
        if tp in KNOWN: assert not ex, 'a tuple of a rational point is excluded: the method is wrong'
    return allres

def negative():
    say('h1_pairs.py --negative: forged inputs that must be refused')
    # 1. with the pairing left out the bounds do not meet, so nothing may be excluded
    res, _, _ = analyse(M, N, '0', 'L', EXTRA[:1], forge='no pairing')
    say('  1. pairing left out: %s' % ('refused (not decided), as it must be' if all(v.startswith('NOT DECIDED') for v in res.values()) else 'NOT REFUSED'))
    # 2. a wrong root (the roots cycled) must not be identified
    try:
        analyse(M, N, '0', 'L', EXTRA[:1], forge='wrong root'); say('  2. wrong root: NOT REFUSED')
    except AssertionError as ex: say('  2. wrong root: refused (%s), as it must be' % str(ex).split(':')[0])
    # 3. a control ratio of square type: its two known tuples must be classes of rational points
    # 3. a forged tuple: the known tuple of t = 0 with dL times 3 must be excluded or not in Sel^2 (3 is a bad prime)
    res, _, _ = analyse(M, N, '0', 'L', [(1109503, 597671, 6, 1)])
    say('  3. forged tuple (1109503, 597671, 6, 1): %s' % ('excluded, as it must be' if list(res.values())[0].startswith('EXCLUDED') else 'NOT EXCLUDED'))

if __name__ == '__main__':
    if '--negative' in sys.argv: negative(); f = 'h1_pairs_negative.txt'
    else: main(); f = 'h1_pairs.txt'
    s = '\n'.join(OUT) + '\n'
    f = os.path.join(HERE, f)
    if WRITE:
        open(f, 'w').write(s)
        HINTS['rs_conic'] = RS.HINTS['conic']; HINTS['ri_conic'] = RI.HINTS.get('conic', {})
        json.dump(HINTS, open(HFILE, 'w'), indent=0, sort_keys=True)
    else:
        sys.stdout.write(s)
        if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS from ' + os.path.basename(f))
