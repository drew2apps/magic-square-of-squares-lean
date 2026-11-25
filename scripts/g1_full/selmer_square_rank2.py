#!/usr/bin/env python3
"""Rank exactly 2 at the square-type curves with Selmer bound 2 (selmer_square_5000.json), when a point is known whose
class under alpha is outside <alpha(T), alpha(R)>.  Experiment, not counted.

    python3 selmer_square_rank2.py [--find] [--write]
      --find   ask PARI (ellrank) for points and record them as hints (selmer_square_rank2_hints.json, with --write)
      without --find the hints file is read and every point is checked here exactly; PARI is not trusted.

Written by an AI model (Claude); not read by a person.

The argument, for one curve E = E+ or E- with s = 2 (exact upper bound, selmer_square.py): 2^rank =
|alpha(E(Q))| |alpha'(E'(Q))| / 4 (Silverman-Tate III.6).  alpha'(T') = -1 is not trivial, so |alpha'(E'(Q))| >= 2.
alpha(T) = H H- and alpha(R) = 2 S H (UNIFORM_RANK.md) generate a subgroup of order 4 (Theorem 1 there).  If a
rational point G has x(G) - e outside the four classes 1, H H-, 2H, 2H- (S is a square; tested by: no product
(x(G) - e) c, c in {1, H H-, 2H, 2H-}, is a rational square), then |alpha(E(Q))| >= 8, so rank >= 2, and with s = 2,
rank = 2.
"""
import os, sys, json
from fractions import Fraction as Fr
from math import isqrt
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import sigma_exact as SE

WRITE = '--write' in sys.argv
HFILE = os.path.join(here, 'selmer_square_rank2_hints.json')
out = []
def say(s): print(s, flush=True); out.append(s)

def H(m, n): return m**4 - 2*m**3*n + 2*m*m*n*n + 2*m*n**3 + n**4
def issq(q):
    q = Fr(q)
    return q >= 0 and isqrt(q.numerator)**2 == q.numerator and isqrt(q.denominator)**2 == q.denominator

def check_curve(m, nn, pts):
    """Return the first hint point (as strings) that proves rank >= 2, or None."""
    M = SE.model(m, nn, 2)
    C = SE.Curve(M['a2'], M['a4'], M['a6'])
    S = m*m + nn*nn; k = isqrt(S); assert k*k == S
    e = -2 * S**3
    assert e**3 + C.a2*e*e + C.a4*e + C.a6 == 0
    xR = Fr(-4*m*nn*S*(m*m - nn*nn))
    assert issq(xR**3 + C.a2*xR*xR + C.a4*xR + C.a6)                   # R is a rational point
    h, hm = H(m, nn), H(m, -nn)
    known = [1, h*hm, 2*h, 2*hm]
    assert issq((xR - e) * 2 * h)                                      # alpha(R) = 2 S H = 2H mod squares
    for (xs, ys) in pts:
        P = (Fr(xs), Fr(ys))
        assert C.on(P), 'a hint point is not on the curve'
        if P[0] == e: continue
        if not any(issq((P[0] - e) * c) for c in known):
            return (xs, ys)
    return None

def find_points(m, nn):
    import cypari2
    pari = cypari2.Pari(); pari.allocatemem(4 * 10**8)
    M = SE.model(m, nn, 2)
    E = pari.ellinit([0, M['a2'], 0, M['a4'], M['a6']])
    S = m*m + nn*nn; xR = -4*m*nn*S*(m*m - nn*nn)
    R = [xR, pari.ellordinate(E, xR)[0]]
    for eff in (1,):
        r = pari.ellrank(E, eff, [R])
        pts = [(str(P[0]), str(P[1])) for P in r[3]]
        if check_curve(m, nn, pts): return pts
    return pts

if __name__ == '__main__':
    DATA = sys.argv[sys.argv.index('--data') + 1] if '--data' in sys.argv else os.path.join(here, 'selmer_square_5000.json')
    data = json.load(open(DATA))
    curves = [(r['m'], r['n']) for r in data if r['bound'] == 2]
    hints = json.load(open(HFILE)) if os.path.exists(HFILE) else {}
    if '--negative' in sys.argv:
        # a point off the curve, and a point whose class is in <alpha(T), alpha(R)> (R itself), must not prove rank 2
        m, nn = curves[0]
        res = []
        key = '%d,%d' % (m, nn)
        xs, ys = hints[key][0] if key in hints and hints[key] else ('0', '0')
        try:
            check_curve(m, nn, [(xs, str(Fr(ys) + 1))]); res.append(False); say('  NOT REFUSED: point off the curve')
        except AssertionError:
            res.append(True); say('  refused: a point off the curve (%d,%d)' % (m, nn))
        M = SE.model(m, nn, 2); S = m*m + nn*nn
        C = SE.Curve(M['a2'], M['a4'], M['a6'])
        xR = Fr(-4*m*nn*S*(m*m - nn*nn)); y2 = xR**3 + C.a2*xR*xR + C.a4*xR + C.a6
        yR = Fr(isqrt(y2.numerator), isqrt(y2.denominator))
        P3 = C.add(C.add((xR, yR), (xR, yR)), (xR, yR))                  # 3R: class of R, must not count
        ok = check_curve(m, nn, [(str(xR), str(yR)), (str(P3[0]), str(P3[1]))]) is None
        res.append(ok); say('  %s: R and 3R (class of R) given as the second point (%d,%d)' % ('refused' if ok else 'NOT REFUSED', m, nn))
        say('negative tests: %d of %d refused' % (sum(res), len(res)))
        sys.exit(0 if all(res) else 1)
    say('selmer_square_rank2.py: the %d curves with Selmer bound 2 in selmer_square_5000.json' % len(curves))
    proved, notyet = [], []
    for (m, nn) in curves:
        key = '%d,%d' % (m, nn)
        if '--find' in sys.argv and key not in hints:
            hints[key] = find_points(m, nn)
            if WRITE: json.dump(hints, open(HFILE, 'w'), indent=0, sort_keys=True)
        G = check_curve(m, nn, hints.get(key, []))
        (proved if G else notyet).append((m, nn))
    say('rank exactly 2 (exact check of the recorded points; upper bound s = 2 from selmer_square.py): %d curves: %s'
        % (len(proved), ' '.join('(%d,%d)' % x for x in proved)))
    say('no point outside <alpha(T), alpha(R)> recorded: %d curves: %s' % (len(notyet), ' '.join('(%d,%d)' % x for x in notyet)))
    if WRITE:
        json.dump(hints, open(HFILE, 'w'), indent=0, sort_keys=True)
        open(os.path.join(here, 'selmer_square_rank2.txt'), 'w').write('\n'.join(out) + '\n')
