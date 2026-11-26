#!/usr/bin/env python3
"""F4: an exact rank bound (no PARI in the bound) for C: y^2 = g4 K, by the repository's rank_indep.rank_bound
(2-isogeny descent and second descent, commit 24d6ef6, unchanged, with sympy).
Written by an AI model (Claude); not read by a person.  Experiment, not counted.   usage: g4k_rank_exact.py [--write] [--negative] [m,n ...]  (no ratio: the three closures, or with --negative two ratios of rank >= 2)"""
import sys, os, io, contextlib, json
HERE = os.path.dirname(os.path.abspath(__file__))
_IN_REPO = os.path.join(HERE, '..', '..', '..')                   # experiments/census-wide/g4k -> repository root
REPO = os.environ.get('MSS_REPO', _IN_REPO)
sys.path.insert(0, HERE)
import dm_g4k
sys.path.insert(0, os.path.join(REPO, 'scripts', 'g1_full')); os.chdir(os.path.join(REPO, 'scripts', 'g1_full'))
import rank_indep as RI
from fractions import Fraction as Fr
import cypari2
pari = cypari2.Pari()

_orig_fcp = RI.find_conic_point
def find_conic_point(d, a, c):
    """rank_indep.find_conic_point, but when qfsolve's point has N = 0, move to the second intersection of the conic
    with a line through it (still only a hint: rank_indep checks the point it is given)"""
    try: return _orig_fcp(d, a, c)
    except AssertionError as ex:
        if 'N = 0' not in str(ex): raise
    from math import gcd
    M = pari.matrix(3, 3, [2 * d, a, 0, a, 2 * c, 0, 0, 0, -2]); P0 = [int(v) for v in pari.qfsolve(M)]
    F = lambda X: d * X[0] ** 2 + a * X[0] * X[1] + c * X[1] ** 2 - X[2] ** 2
    Bf = lambda X, Y: (F([X[i] + Y[i] for i in range(3)]) - F(X) - F(Y))
    for D in ([1, 0, 1], [0, 1, 1], [1, 1, 1], [1, 2, 1], [2, 1, 3]):
        if F(D) == 0: continue
        lam = Fr(-Bf(P0, D), F(D)); X = [P0[i] + lam * D[i] for i in range(3)]
        den = 1
        for v in X: den = den * v.denominator // gcd(den, v.denominator)
        X = [int(v * den) for v in X]; g = gcd(gcd(X[0], X[1]), X[2]); X = [v // g for v in X]
        if X[2] != 0 and F(X) == 0: return tuple(X)
    raise AssertionError('no conic point with N != 0 found')
RI.find_conic_point = find_conic_point

NEG = '--negative' in sys.argv; WRITE = '--write' in sys.argv
ARGS = [a for a in sys.argv[1:] if not a.startswith('--')] or (['3540,2581', '6671,2209'] if NEG else ['7735,4968', '6557,324', '6881,6233'])
LINES = []
_print = print
def print(*a): LINES.append(' '.join(str(v) for v in a)); _print(*a)
for arg in ARGS:
    m, n = (int(v) for v in arg.split(','))
    ainv = dm_g4k.g4k_curve(m, n)[0]; _, a2, _, a4, a6 = ainv
    rts = [int(r) for r in pari.nfroots(None, pari('x^3 + (%d)*x^2 + (%d)*x + (%d)' % (a2, a4, a6)))]
    assert len(rts) >= 1
    e = rts[0]; assert e ** 3 + a2 * e * e + a4 * e + a6 == 0
    a, b = 3 * e + a2, 3 * e * e + 2 * a2 * e + a4            # y^2 = X(X^2 + a X + b), X = x - e
    E = pari.ellinit(ainv); pts = pari.ellrank(E, 2)[3]
    if len(pts) == 0: pts = pari.ellrank(E, 4)[3]
    G = [Fr(int(pari.numerator(pts[0][0])), int(pari.denominator(pts[0][0]))) - e,
         Fr(int(pari.numerator(pts[0][1])), int(pari.denominator(pts[0][1])))]
    assert RI.on_curve([0, a, 0, b, 0], G), 'hint point not on the shifted curve'
    say = []
    with contextlib.redirect_stdout(io.StringIO()):
        bound, rec = RI.rank_bound(a, b, [G], RI.find_points([0, -2 * a, 0, a * a - 4 * b, 0]), say=say.append)
    print('(%d,%d) g4*K: 2-torsion roots %d; shifted model a = %d, b = %d' % (m, n, len(rts), a, b))
    for s in say: print(s)
    print('exact rank bound (rank_indep, no PARI in the bound): %d; a point of infinite order is known, so rank >= 1' % bound)
    if NEG:
        lo = int(pari.ellrank(E, 2)[0])
        print('  negative control: PARI finds %d independent points, so the exact bound must be >= %d: %s' % (lo, lo, 'yes' if bound >= lo else 'NO'))
        assert bound >= lo
txt = '\n'.join(LINES) + '\n'
name = os.path.join(HERE, 'g4k_rank_exact_negative.txt' if NEG else 'g4k_rank_exact.txt')
if WRITE: open(name, 'w').write(txt)
elif not sys.argv[1:] or sys.argv[1:] == ['--negative']:
    assert txt == open(name).read(), 'output differs from ' + os.path.basename(name)
    _print('# matches ' + os.path.basename(name))
