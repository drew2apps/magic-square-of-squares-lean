#!/usr/bin/env python3
"""Rank bounds of E1 and of its twist E2 for the ratios of a scope file, with more effort (EXPERIMENT; LAST_PAIR.md, section 10).

For each ratio of the scope file: PARI's ellrank with effort 6 on E1 and on E2 (the Jacobian of Y^2 = d q+ L2 L4, the twist
of E1 by -1).  Where E2 has rank exactly 1 the route of LAST_PAIR.md applies with -d in place of d (the two conditions are the
same, because -1 is a square in Q(i), and the two curves change places).  One line per ratio.  Slow; not rerun in CI.
usage: rank_pairs.py scope.txt"""
import sys, re
from second_descent import *
import cypari2
pari = cypari2.Pari(); pari.allocatemem(4*10**9)
tally = {}
for l in open(sys.argv[1]):
    mm = re.match(r'\((\d+),(\d+)\): d = (\d+), rank of E1 in \[(\d+), (\d+)\]', l)
    if not mm: continue
    m, n, d, lo, hi = (int(x) for x in mm.groups()); fam = family(m, n); out = []
    for key in ('+13', '+24'):
        Q, r = fam[key]; E, A = weierstrass(d*Q, r)
        b = pari.ellrank(pari.ellinit([pari(str(c)) for c in E.minimal_model().a_invariants()]), 6); out.append((int(b[0]), int(b[1])))
    kind = 'sum of the upper bounds %d' % (out[0][1] + out[1][1])
    tally[(out[0], out[1])] = tally.get((out[0], out[1]), 0) + 1
    print('(%d,%d): E1 in [%d, %d], E2 in [%d, %d]; %s' % (m, n, out[0][0], out[0][1], out[1][0], out[1][1], kind), flush=True)
print('bounds (E1, E2) and number of ratios:', sorted(tally.items()))
