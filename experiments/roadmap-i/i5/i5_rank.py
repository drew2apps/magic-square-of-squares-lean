"""I5: H7's step, or else rank_indep.py's exact descent, on the census tables' rows that are conditional on a rank being 1
(EXPERIMENT, written by an AI model; not read by a person; nothing here changes the census).
For each curve: the model is rebuilt from m, n (and d) by the census's own code (census-wide/sd_pari.py: family(m,n)
['+13'], weierstrass(d Q, r) for the route curve E1; the curve is the one route_pari.py and rank_hard.py use); PARI
gives a minimal model and the isogeny class (finders), and each is checked exactly: the minimal model by c4, c6 up to
u^4, u^6, each isogenous curve by equal point counts mod every good prime < 1000.
  H7's step needs a curve with full rational 2-torsion in the class (then the repository's complete 2-descent with the
  Cassels-Tate pairing, rank_second.py through H1's bound_full).  Otherwise: rank_indep.py's rank_bound (an exact
  descent via 2-isogeny with an exact second descent, unchanged), on the model shifted to y^2 = x (x^2 + a x + b).
  A point of infinite order (PARI's ellrank as a hint; checked on the curve, n G != O for n <= 12) gives rank >= 1.
usage: i5_rank.py [label ...] [--negative] [--write]     (i5_rank.txt, i5_rank_negative.txt)
Folders: MSS_REPO (the repository), H1_DIR (roadmap H's h1, for bound_full)."""
import sys, os, io, contextlib, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get('MSS_REPO', os.path.join(HERE, '..', '..', '..'))
H1 = os.environ.get('H1_DIR', next((d for d in (os.path.join(HERE, '..', '..', 'roadmap-h', 'h1'), os.path.join(HERE, '..', 'h1')) if os.path.isdir(d)), os.path.join(HERE, '..', '..', 'roadmap-h', 'h1')))
sys.path.insert(0, os.path.join(REPO, 'experiments', 'census-wide')); sys.path.insert(0, os.path.join(REPO, 'scripts', 'g1_full')); sys.path.insert(0, H1)
_argv = sys.argv; sys.argv = [sys.argv[0]]
with contextlib.redirect_stdout(io.StringIO()):
    import sd_pari as SD
    import rank_indep as RI
    import h1_rank0_exact as R
sys.argv = _argv
pari = SD.pari
OUT = []
def say(s=''): OUT.append(s)
# (label, source, m, n, d, rows of the tables that rest on rank 1 of this curve); source E1 = the route curve
# (census-wide/sd_pari.py), otherwise a curve of census2000/dm_pari.py's all_curves(m, n) by its name
CURVES = [('E1(3540,2581)', 'E1', 3540, 2581, -11740078, '(3540,2581) [CENSUS5000], (6121,959) [CENSUS7000, its partner]'),
          ('E1(4440,2231)', 'E1', 4440, 2231, -13342, '(4440,2231) [CENSUS5000], (6671,2209) [CENSUS7000, its partner]'),
          ('E1(6480,4961)', 'E1', 6480, 4961, 709342, '(6480,4961) [CENSUS7000]'),
          ('E1(8249,3720)', 'E1', 8249, 3720, -108415202, '(8249,3720) [CENSUS9000]'),
          ('E1(10512,145)+', 'E1', 10512, 145, 220962238, '(10512,145), (10657,10367) [CENSUS11000] (route log_10512_220962238)'),
          ('E1(10512,145)-', 'E1', 10512, 145, -220962238, '(10512,145), (10657,10367) [CENSUS11000] (route log_10512_-220962238)'),
          ('D1(12007,1799)', 'D1', 12007, 1799, None, '(6903,5104) [CENSUS7000] (DM on D1 of its partner, assuming rank 1); (12007,1799) [G3]'),
          ('Sig+(12204,9853)', 'Sig+', 12204, 9853, None, 'G3 open (12204,9853): its only curve with PARI bounds [1,3]'),
          ('Sig+(22057,2351)', 'Sig+', 22057, 2351, None, 'G3 open (12204,9853), through its partner (22057,2351): PARI bounds [1,3]')]
with contextlib.redirect_stdout(io.StringIO()):
    sys.path.insert(0, os.path.join(REPO, 'experiments', 'census2000')); import dm_pari as DM

def c46(ai):
    E = pari.ellinit(ai); return Fr(str(E[9])), Fr(str(E[10]))

def model(src, m, n, d):
    if src == 'E1':
        Q, r = SD.family(m, n)['+13']; ain, A = SD.weierstrass(d * Q, r)
        ai = [pari(x.numerator) / x.denominator for x in ain]
    else:
        ai = [c[0] for c in DM.all_curves(m, n) if c[4] == src][0]; ai = [pari(int(v)) for v in ai]
    Em = pari.ellminimalmodel(pari.ellinit(ai))
    mi = [int(v) for v in Em[:5]]
    (c4, c6), (d4, d6) = c46(ai), c46(mi)
    u2 = (c6 / d6) / (c4 / d4)
    assert u2 * u2 == c4 / d4 and u2 ** 3 == c6 / d6, 'the minimal model is not the curve'
    return mi

def npts(ai, p):
    a1, a2, a3, a4, a6 = [x % p for x in ai]; s = 1
    for x in range(p):
        # y^2 + a1 x y + a3 y = f(x): discriminant in y
        D = ((a1 * x + a3) ** 2 + 4 * (x ** 3 + a2 * x * x + a4 * x + a6)) % p
        s += 1 if D == 0 else (2 if pow(D, (p - 1) // 2, p) == 1 else 0)
    return s

def isog_ok(ai, bi):
    Ea, Eb = pari.ellinit(ai), pari.ellinit(bi)
    N = 6 * int(Ea[11]) * int(Eb[11])          # the discriminants of the two models
    return all(npts(ai, p) == npts(bi, p) for p in range(3, 1000) if RI.is_prime_td(p) and N % p)

def short(ai):
    """y^2 = x^3 + A2 x^2 + A4 x + A6 (integral, a1 = a3 = 0) isomorphic over Q to ai (complete the square, scale by 4)"""
    a1, a2, a3, a4, a6 = ai
    b2, b4, b6 = a1 * a1 + 4 * a2, 2 * a4 + a1 * a3, a3 * a3 + 4 * a6
    return [0, b2, 0, 8 * b4, 16 * b6]          # y^2 = x^3 + b2 x^2 + 8 b4 x + 16 b6  (x -> 4x, y -> 8y + ...)

def int_roots(s):
    import sympy
    X = sympy.Symbol('X')
    return sorted(int(r) for r in sympy.Poly(X**3 + s[1]*X**2 + s[3]*X + s[4], X).ground_roots() if r.is_Integer)

def run(lab, src, m, n, d, rows):
    t0 = time.time()
    say('%s: curve %s of (%d,%d)%s; the rows resting on its rank: %s' % (lab, src, m, n, '' if d is None else ' with d = %d' % d, rows))
    mi = model(src, m, n, d)
    say('  minimal model %s (checked: same c4, c6 up to u^4, u^6)' % mi)
    cls = sorted([int(v) for v in pari.ellminimalmodel(pari.ellinit(c))[:5]] for c in pari.ellisomat(pari.ellinit(mi), 0, 1)[0])
    for c in cls: assert isog_ok(mi, c), 'an isogenous curve with other point counts'
    roots = [int_roots(short(c)) for c in cls]
    say('  isogeny class (PARI, each checked by point counts mod p < 1000): %d curves, rational 2-torsion points %s' % (len(cls), [len(r) for r in roots]))
    best = None
    for c, rs in zip(cls, roots):
        if len(rs) == 3:
            sc = short(c); e = rs[0]; a, b = 3*e + sc[1], 3*e*e + 2*sc[1]*e + sc[3]
            es = R.roots_of(a, b); r, info = R.bound_full(es)
            say('  full rational 2-torsion on %s: H7\'s step (complete 2-descent and Cassels-Tate pairing): %s, rank <= %d' % (c, info, r))
        elif len(rs) == 1:
            sc = short(c)
            e, a, b = RI.shift_to_two_torsion(sc[1], sc[3], sc[4]); ap, bp = -2 * a, a * a - 4 * b
            pE = [(P[0] - e, P[1]) for P in RI.find_points(sc)]
            say('  on %s (one rational 2-torsion point): rank_indep.py\'s descent via 2-isogeny and second descent' % c)
            r, rec = RI.rank_bound(a, b, pE, RI.find_points([0, ap, 0, bp, 0]), say=say)
        else: continue
        best = r if best is None else min(best, r)
    sc = short(mi); pts = RI.find_points(sc); inf = False
    if pts:
        rs = int_roots(sc)
        if rs:
            e = rs[0]; a, b = 3*e + sc[1], 3*e*e + 2*sc[1]*e + sc[3]; G = [Fr(pts[0][0]) - e, Fr(pts[0][1])]
            inf = all(R.mul(a, b, k, G) is not None for k in range(1, 13))
    say('  exact rank bound (least over the class) %s; a point of infinite order %s' % (best, 'known, so rank >= 1' if inf else 'not checked'))
    say('  root number %d (PARI ellrootno; under the parity conjecture the rank is %s)' % (int(pari.ellrootno(pari.ellinit(mi))), 'odd' if int(pari.ellrootno(pari.ellinit(mi))) == -1 else 'even'))
    if best == 1 and inf: say('  so rank = 1 (computed, exact)')
    return best

def negative():
    say('I5 negative and positive controls. Written by an AI model; not read by a person.')
    mi = model('E1', 4440, 2231, -13342)
    bad = list(mi); bad[3] += 4 * 10 ** 6
    say('forged isogenous curve (a4 + 4*10^6): %s' % ('refused' if not isog_ok(mi, bad) else 'ACCEPTED')); assert not isog_ok(mi, bad)
    good = sorted([int(v) for v in pari.ellminimalmodel(pari.ellinit(c))[:5]] for c in pari.ellisomat(pari.ellinit(mi), 0, 1)[0])
    say('positive control, the 2-isogenous curve PARI finds: %s' % ('accepted' if all(isog_ok(mi, c) for c in good) else 'REFUSED')); assert all(isog_ok(mi, c) for c in good)
    sc = short(mi); P = RI.find_points(sc)[0]; Q = (P[0], P[1] + 1)
    on = lambda T: Fr(T[1]) ** 2 == Fr(T[0]) ** 3 + sc[1] * Fr(T[0]) ** 2 + sc[3] * Fr(T[0]) + sc[4]
    say('forged point (y + 1): %s' % ('refused (not on the curve)' if not on(Q) else 'ACCEPTED')); assert not on(Q)
    say('positive control, the hint point: %s' % ('on the curve' if on(P) else 'NOT on the curve')); assert on(P)

if __name__ == '__main__' and '--negative' in sys.argv:
    negative(); s = '\n'.join(OUT) + '\n'; f = os.path.join(HERE, 'i5_rank_negative.txt')
    if '--write' in sys.argv: open(f, 'w').write(s)
    else: sys.stdout.write(s); print('# matches' if os.path.exists(f) and open(f).read() == s else '# DIFFERS')
elif __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    say('I5: exact rank bounds for the curves of the rows conditional on a rank being 1. Written by an AI model; not read by a person. Experiment, not counted.')
    for c in CURVES:
        if not args or c[0] in args: run(*c)
    s = '\n'.join(OUT) + '\n'; f = os.path.join(HERE, 'i5_rank.txt')
    if '--write' in sys.argv and not args: open(f, 'w').write(s)
    else:
        sys.stdout.write(s)
        if not args and os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS')
