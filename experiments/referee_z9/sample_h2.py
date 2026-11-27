"""Review z9 of roadmap H (run 7), written by an AI model; not read by a person.  Experiment, not counted.
Reruns a sample of H2's exhaustive range 15000 < m <= 30000 without H2's code: 20 random class-keeping ratios of
h2_ex_15000_*of3.txt (seed 909), the 3 with the largest P, and 20 uniformly random coprime ratios of the range.
Class status: the repository's criterion (class_selmer_linalg.system / solve) with sympy's factorint in place of PARI.
P: G2's g2_killsets.run (roadmap-g/g2, run 6), the irredundant kill set; P = its largest good prime (0 if none).
Compared: class or not, the dimension, log2 of the tuple space, P.  Rests on: the repository's criterion, G2's code and
G1's exact_local.py (run 6; tested, not proved), PARI's factorisations inside G2 (each factor proved prime).
usage: python3 -I sample_h2.py [--write]   (sample_h2.txt; about 1 minute)"""
import os, sys, random
from math import gcd, isqrt, log2
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import *
O = Out('sample_h2.txt')
os.environ.setdefault('MSS_REPO', REPO)
sys.path.insert(0, G2); sys.path.insert(0, G1)
argv = sys.argv; sys.argv = [argv[0]]
import g2_killsets as K
sys.argv = argv
sys.path.insert(0, STRUCT)
import class_selmer as C, class_selmer_linalg as CL
from sympy import factorint
C.fac = lambda x: sorted(factorint(abs(int(x))))
O('Review z9: sample of H2 15000 < m <= 30000 against the repository criterion (sympy factorisations) and G2 kill sets. Written by an AI model; not read by a person.')
rows = {}
for i in range(3):
    for l in open(os.path.join(H2, 'h2_ex_15000_%dof3.txt' % i)):
        if not l.startswith('m=') and l.strip():
            a, dim, ls, P = l.split(); m, n = map(int, a.split(',')); rows[(m, n)] = (int(dim), int(ls), P)
rnd = random.Random(909); keys = sorted(rows)
samp = rnd.sample(keys, 20) + sorted(keys, key=lambda r: (-int(rows[r][2]), r))[:3]
rand = []
while len(rand) < 20:
    m = rnd.randrange(15001, 30001); n = rnd.randrange(1, m)
    if gcd(m, n) == 1: rand.append((m, n))
def deg(m, n): N = m * m + n * n; return isqrt(N) ** 2 == N or isqrt(2 * N) ** 2 == 2 * N
bad_n = 0
for tag, lst in (('class-keeping', samp), ('random', rand)):
    for m, n in lst:
        st, dim = 'degenerate', None
        if not deg(m, n):
            sy = CL.system(m, n)
            if sy is None: st = 'no point'
            else:
                ok, dim = CL.solve(*sy); st = 'class' if ok else 'no class'
        h = rows.get((m, n)); agree = (st == 'class') == (h is not None); line = '%s (%d,%d): here %s' % (tag, m, n, st)
        if st == 'class':
            space, lb, kill, known = K.run(m, n)
            badset = set()
            for x in (2*m*n*(m*m - n*n), m*m + n*n, m**8 + 14*m**4*n**4 + n**8): badset |= set(factorint(x))
            g = [p for p in kill if p not in badset] if kill is not None else None
            P = 'SURVIVOR' if kill is None else str(max(g)) if g else '0'
            line += ' dim %d log2 %.0f P %s (bad primes alone: %s; kill set %s)' % (dim, log2(space), P, 'EMPTY' if not lb else '%d left' % len(lb), kill)
            if h is not None: agree = agree and h == (dim, round(log2(space)), P)
        line += ' | H2 %s | %s' % ('dim %d log2 %d P %s' % h if h else 'not listed', 'AGREE' if agree else 'DISAGREE')
        bad_n += not agree; O(line)
O('disagreements: %d of %d' % (bad_n, len(samp) + len(rand)))
O.finish()
