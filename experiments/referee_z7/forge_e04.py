#!/usr/bin/env python3
"""Review z7 of roadmap H, H7 claim 4: forgeries against h7_e04_rank.py (imported, unchanged).
Written by an AI model (review agent z7, 27 November 2026); not read by a person.  Experiment, not counted.

(1) A known rank >= 1 never gets a lower bound: on 24 random curves E04(m,n) (m <= 400, seed 2026) and 9 congruent-
    number curves y^2 = x^3 - n^2 x, the bundle's least_bound (rank_second's descent and pairing over the walk) is
    compared with an exact lower bound for the rank: the F2-dimension, modulo the 2-torsion images, of the images in
    Sel^2 (my2desc.py's Kummer map) of the points stored in z7_points.json (found by PARI's ellrank with --find; each
    checked on the curve exactly here).  A violation (bound < lower bound) would be an error of the bundle.
(2) A point of order 4 as the "generator": x = -186116 on E04(23,7) (torsion Z/4 x Z/2); checked here to have order 4
    by the exact group law (H1's mul); h7_e04_rank.run_curve must refuse it by its infinite-order test (the committed
    --negative only tries x = 0, which is refused earlier, as y = 0).
Rests on: the bundle's code for what is tested; nothing quoted beyond ct_mine/my2desc.
PARI: inside the bundle's rank_second.py (conic points); with --find, ellrank to find points (a finder only).

    python3 forge_e04.py           prints forge_e04.txt   (about 4 minutes)
    python3 forge_e04.py --find    rewrite the forge part of z7_points.json (PARI ellrank), then exit
"""
import sys, os, io, json, random, contextlib
from fractions import Fraction as Fr
from math import gcd, isqrt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
REPO = os.environ.get('MSS_REPO', os.path.join(HERE, '..', '..'))
H7 = os.environ.get('H7_DIR', os.path.join(HERE, '..', 'roadmap-h', 'h7'))
H1 = os.environ.get('H1_DIR', os.path.join(HERE, '..', 'roadmap-h', 'h1'))
os.environ.setdefault('MSS_REPO', REPO); os.environ.setdefault('H1_DIR', H1)
sys.path.insert(0, H7)
_a = sys.argv; sys.argv = [sys.argv[0]]
with contextlib.redirect_stdout(io.StringIO()):
    import h7_e04_rank as H
sys.argv = _a
from my2desc import Curve, span
PTS = os.path.join(HERE, 'z7_points.json')

def tests():
    rng = random.Random(2026); out = []
    while len(out) < 24:
        m = rng.randint(3, 400); n = rng.randint(1, m - 1)
        if gcd(m, n) == 1: out.append(('E04', m, n))
    return out + [('cong', k, 0) for k in (5, 6, 7, 34, 41, 65, 137, 210, 1254)]

def model(kind, m, n):
    if kind == 'E04':
        S = m * m + n * n; return (0, -2 * S * m * m, -2 * S * n * n)
    return (-m, 0, m)

def find():
    import cypari2
    pari = cypari2.Pari(); pari.setrand(1)
    data = json.load(open(PTS)) if os.path.exists(PTS) else {}
    out = {}
    for kind, m, n in tests():
        e1, e2, e3 = model(kind, m, n)
        E = pari.ellinit([0, -(e1 + e2 + e3), 0, e1 * e2 + e1 * e3 + e2 * e3, -e1 * e2 * e3])
        out['%s %d,%d' % (kind, m, n)] = [[str(P[0]), str(P[1])] for P in pari.ellrank(E, 1)[3]]
    data['forge_e04'] = out
    json.dump(data, open(PTS, 'w'), indent=1, sort_keys=True)
    print('wrote', PTS)

def main():
    print('# forge_e04.py (review z7; written by an AI model): forgeries against h7_e04_rank.py')
    pts = json.load(open(PTS))['forge_e04']; bad = 0
    for kind, m, n in tests():
        es = model(kind, m, n); E = Curve(*es)
        f = lambda x: (x - es[0]) * (x - es[1]) * (x - es[2])
        P = [(Fr(x), Fr(y)) for x, y in pts['%s %d,%d' % (kind, m, n)]]
        for x, y in P: assert y * y == f(x), 'a stored point is not on the curve'
        tors = [E.classvec(E.kummer(Fr(es[0]))), E.classvec(E.kummer(Fr(es[1])))]
        lo = len(span(tors + [E.classvec(E.kummer(x)) for x, _ in P])) - len(span(tors))
        rb, info = H.least_bound(tuple(sorted(es)))
        ok = rb is not None and rb >= lo
        bad += not ok
        nm = 'E04 (%d,%d)' % (m, n) if kind == 'E04' else 'y^2 = x^3 - %d^2 x' % m
        print('%s: exact lower bound %d (images of %d stored point(s)); bundle least_bound %s: %s' %
              (nm, lo, len(P), rb, 'consistent' if ok else 'VIOLATION'))
    print('violations: %d' % bad)
    m, n, x = 23, 7, Fr(-186116)
    S = m * m + n * n; a, b = 2 * S * S, 4 * S * S * m * m * n * n
    y2 = x * (x * x + a * x + b); y = Fr(int(y2.numerator ** 0.5 + 0.5));
    from math import isqrt
    y = Fr(isqrt(y2.numerator), isqrt(y2.denominator)); assert y * y == y2
    Q = (x, y)
    assert H.R.mul(a, b, 2, Q) is not None and H.R.mul(a, b, 4, Q) is None
    try:
        H.run_curve(m, n, x); print('order-4 point x = %s on E04(23,7) as the generator: ACCEPTED' % x)
    except AssertionError as ex:
        print('order-4 point x = %s on E04(23,7) (2P != O, 4P = O, exact) as the generator: refused (%s)' % (x, ex))
    print('# end')

if __name__ == '__main__':
    find() if '--find' in sys.argv else main()
