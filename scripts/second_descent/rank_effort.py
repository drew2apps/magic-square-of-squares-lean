#!/usr/bin/env python3
"""PARI's rank bounds for E1 with more effort, for the ratios of a scope file whose first bounds are 1 and 3 (EXPERIMENT).

For each such ratio: ellrank with effort 3, and with effort 6 if the upper bound is still above the lower one; each call in a
process of its own with a limit of 150 seconds.  One line per call: RES m n effort lower upper.  Not rerun in CI.
usage: rank_effort.py scope.txt"""
import sys, re, subprocess
code = r'''
import sys
sys.argv = ['x']
from second_descent import *
import cypari2
pari = cypari2.Pari(); pari.allocatemem(4*10**9)
m, n, d, eff = %d, %d, %d, %d
Q, r = family(m, n)['+13']; E, A = weierstrass(d*Q, r)
b = pari.ellrank(pari.ellinit([pari(str(c)) for c in E.minimal_model().a_invariants()]), eff)
print('RES', m, n, eff, int(b[0]), int(b[1]))
'''
for l in open(sys.argv[1]):
    if '[1, 3]' not in l: continue
    m, n, d = (int(x) for x in re.match(r'\((\d+),(\d+)\): d = (\d+)', l).groups())
    for eff in (3, 6):
        try:
            r = subprocess.run([sys.executable, '-c', code % (m, n, d, eff)], capture_output=True, text=True, timeout=150)
            out = [x for x in r.stdout.splitlines() if x.startswith('RES')]
            print(out[0] if out else 'FAIL %d %d %d' % (m, n, eff), flush=True)
            if out and out[0].split()[-1] == out[0].split()[-2]: break
        except subprocess.TimeoutExpired:
            print('TIMEOUT', m, n, eff, flush=True)
