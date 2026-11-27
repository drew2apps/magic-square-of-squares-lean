"""Review z9 of roadmap H (run 7), written by an AI model; not read by a person.  Experiment, not counted.
(a) 20 random ratios of 15000 < m <= 30000 that have a local point everywhere (rejection sampling, seed 919): class or
no class by the repository's criterion with sympy's factorint, against whether h2_ex_15000_*of3.txt lists them (and the
dimension); (b) the whole count line of m = 25789 (the largest P) and m = 29999, recomputed for every n: coprime,
degenerate, no point, no class, class, and the list of class-keeping ratios with their dimensions.
Rests on: the repository's criterion (scripts/cover_descent/structure).  usage: python3 -I sample_h2b.py [--write]
(sample_h2b.txt; about 2 minutes)"""
import os, sys, random
from math import gcd, isqrt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import *
O = Out('sample_h2b.txt')
sys.path.insert(0, STRUCT)
import class_selmer as C, class_selmer_linalg as CL
from sympy import factorint
C.fac = lambda x: sorted(factorint(abs(int(x))))
O('Review z9: class status of 20 ratios with points everywhere, and two whole values of m, against H2. Written by an AI model; not read by a person.')
rows = {}; mline = {}
for i in range(3):
    for l in open(os.path.join(H2, 'h2_ex_15000_%dof3.txt' % i)):
        if l.startswith('m='):
            w = l.split(); mline[int(w[0][2:])] = dict((x.split('=')[0], int(x.split('=')[1])) for x in w[1:])
        elif l.strip():
            a, dim, ls, P = l.split(); m, n = map(int, a.split(',')); rows[(m, n)] = int(dim)
def deg(m, n): N = m * m + n * n; return isqrt(N) ** 2 == N or isqrt(2 * N) ** 2 == 2 * N
rnd = random.Random(919); got = 0; bad = 0
while got < 20:
    m = rnd.randrange(15001, 30001); n = rnd.randrange(1, m)
    if gcd(m, n) != 1 or deg(m, n): continue
    sy = CL.system(m, n)
    if sy is None: continue
    ok, dim = CL.solve(*sy); mine = dim if ok else None; h = rows.get((m, n)); got += 1; bad += mine != h
    O('points everywhere (%d,%d): here %s; H2 %s; %s' % (m, n, 'class dim %d' % mine if ok else 'no class', 'class dim %d' % h if h is not None else 'not listed', 'AGREE' if mine == h else 'DISAGREE'))
for m in (25789, 29999):
    c = {'coprime': 0, 'degenerate': 0, 'no_point': 0, 'no_class': 0, 'class': 0}; cl = set()
    for n in range(1, m):
        if gcd(m, n) != 1: continue
        c['coprime'] += 1
        if deg(m, n): c['degenerate'] += 1; continue
        sy = CL.system(m, n)
        if sy is None: c['no_point'] += 1; continue
        ok, dim = CL.solve(*sy)
        if not ok: c['no_class'] += 1
        else: c['class'] += 1; cl.add((m, n, dim))
    hcl = set((mm, nn, d) for (mm, nn), d in rows.items() if mm == m); a = c == mline[m] and cl == hcl; bad += not a
    O('m=%d here %s | H2 %s | class lists equal: %s | %s' % (m, ' '.join('%s=%d' % kv for kv in c.items()), ' '.join('%s=%d' % kv for kv in mline[m].items()), cl == hcl, 'AGREE' if a else 'DISAGREE'))
O('disagreements: %d' % bad)
O.finish()
