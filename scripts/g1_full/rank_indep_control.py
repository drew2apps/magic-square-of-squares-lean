#!/usr/bin/env python3
"""A control of the second descent of rank_indep.py on curves where the answer is known from rational points.

    python3 rank_indep_control.py [count]      output: rank_indep_control.txt

Written by an AI model (Claude); not read by a person.

For curves y^2 = x(x^2 + a x + b) with small a, b: PARI's ellrank gives points (a hint; each is checked on the curve).
The class of x(P) of every point P found (the points themselves, not their sums) is in the image of alpha, so the curve C_d
has a rational point and the second descent must NOT exclude it; the same on E' with the points PARI finds there.
The curves are random (a, b) with |a| <= 200, |b| <= 3000 for which the first descent is above PARI's lower bound or
that lower bound is at least 2.  The program also records the rank bound of rank_indep.py with no known points given (only the
torsion), and checks that it is at least the number of independent points PARI found (if it were smaller, the descent
would be wrong).  It also counts the curves where the second descent lowers the bound of the first descent.
"""
import os, sys, itertools
WRITE = '--write' in sys.argv   # outputs and hint files are written only with --write, so CI can compare stdout with them
from fractions import Fraction as Fr
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import rank_indep as RI
import cypari2
pari = cypari2.Pari()

def pts_of(ainvs):
    E = pari.ellinit(ainvs)
    r = pari.ellrank(E, 1)
    pts = [(Fr(str(P[0])), Fr(str(P[1]))) for P in r[3]]
    return int(r[0]), int(r[1]), pts, E

def main(count):
    out = []
    def say(s): print(s, flush=True); out.append(s)
    n = 0; lowered = 0; checked_classes = 0
    import random
    rng = random.Random(2026)
    tried = 0
    while tried < 20000:
        tried += 1
        a, b = rng.randint(-200, 200), rng.randint(-3000, 3000)
        if n >= count: break
        if b == 0 or a * a - 4 * b == 0: continue
        bp = a * a - 4 * b
        if bp > 0 and RI.isqrt(bp) ** 2 == bp: continue
        if b > 0 and RI.isqrt(b) ** 2 == b: continue                 # E' would have more 2-torsion
        R = RI.first_descent(a, b)
        fb = RI.pow2_floor(len(R['S'])) + RI.pow2_floor(len(R['Sp'])) - 2
        lo, hi, pts, Ep = pts_of([0, a, 0, b, 0])
        if not (fb > lo or lo >= 2): continue                           # keep curves with Sha[phi] or with rank >= 2
        loP, hiP, ptsP, _ = pts_of([0, -2 * a, 0, bp, 0])
        n += 1
        for P in pts: assert RI.on_curve([0, a, 0, b, 0], P)
        for P in ptsP: assert RI.on_curve([0, -2 * a, 0, bp, 0], P)
        # every class of a known point (and of sums) must survive the second descent
        for (A, B, ps, S) in ((a, b, pts, R['S']), (-2 * a, bp, ptsP, R['Sp'])):
            cls = set()
            Ell = RI
            for P in ps:
                if P[0] != 0: cls.add(RI.sq_class_key(P[0]))
            for d in sorted(cls):
                if d in (1, RI.sqfree_part(B)): continue
                assert d in S, ('class of a point outside the Selmer group', a, b, d)
                surv, _ = RI.second_descent(d, A, B, RI.find_conic_point(d, A, B // d))
                assert surv, ('the second descent excludes the class of a rational point', a, b, d)
                checked_classes += 1
        r, rec = RI.rank_bound(a, b, [], [], say=lambda s: None)
        assert r >= lo, ('rank bound below the number of independent points', a, b, r, lo)
        if r < fb: lowered += 1
        say('a=%d b=%d: PARI [%d,%d]; first descent <= %d; with the second descent <= %d' % (a, b, lo, hi, fb, r))
    say('rank_indep_control: %d curves; %d classes of rational points tested, none excluded; the bound never below '
        'PARI\'s lower bound; the second descent lowered the bound on %d curves' % (n, checked_classes, lowered))
    if WRITE: open(os.path.join(here, 'rank_indep_control.txt'), 'w').write('\n'.join(out) + '\n')

def census(keys):
    """Curves of the census (magma_sigma*.json) where Magma's bound is 2 or 3: the bound here, with no known point
    given, must be at least PARI's lower bound."""
    import re
    A = RI.magma_ainvs(); mag = {}
    for f in ('magma_sigma2_raw.txt', 'magma_sigma3_raw.txt'):
        for mm, nn, c, r in re.findall(r'BOUNDS (\d+) (\d+) (E[+-]) (\d+)', open(os.path.join(here, f)).read()):
            mag[((int(mm), int(nn)), c)] = int(r)
    out = []
    for k in keys:
        ai = A[k]
        r = pari.ellrank(pari.ellminimalmodel(pari.ellinit(ai)), 1)
        e, a, b = RI.shift_to_two_torsion(ai[1], ai[3], ai[4])
        rb, _ = RI.rank_bound(a, b, [], [], say=lambda s: None)
        assert rb >= int(r[0]), ('below PARI\'s lower bound', k)
        out.append('%d,%d %s: Magma %d, PARI [%d,%d], rank_indep (no known points given) %d' % (k[0][0], k[0][1], k[1], mag[k], int(r[0]), int(r[1]), rb))
        print(out[-1], flush=True)
    if WRITE: open(os.path.join(here, 'rank_indep_control_census.txt'), 'w').write('\n'.join(out) + '\n')

if __name__ == '__main__':
    if '--census' in sys.argv:
        census([((224, 207), 'E-'), ((260, 69), 'E-'), ((144, 17), 'E-'), ((224, 207), 'E+'), ((273, 136), 'E+'), ((208, 105), 'E-')])
    else:
        main(int(sys.argv[1]) if len(sys.argv) > 1 else 40)
