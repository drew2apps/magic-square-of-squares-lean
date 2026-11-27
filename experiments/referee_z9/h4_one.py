"""Review z9 of roadmap H (run 7), written by an AI model; not read by a person.  Experiment, not counted.
Reruns H4 (experiments/roadmap-h/h4) on G1's three out-of-memory ratios with h4_split.line (G2's kill-set computation,
the join split by the class d0; read in review z9: the split join is the full join, since exact_local.join is a
filter on each tuple), and compares each line with the committed h4_split.txt.
Rests on: G2's g2_killsets.py and G1's exact_local.py (run 6; tested, not proved), PARI's factorisations (each factor
proved prime).  usage: ulimit -v 6500000; python3 -I h4_one.py [--write]   (h4_one.txt; about 2 minutes of CPU per ratio
on PARI 2.17.4, 350 s for the three; not for CI)"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import *
O = Out('h4_one.txt')
os.environ.setdefault('MSS_REPO', REPO); os.environ.setdefault('G1_DIR', G1); os.environ.setdefault('G2_DIR', G2)
sys.path.insert(0, H4)
argv = sys.argv; sys.argv = [argv[0]]
import h4_split as H
sys.argv = argv
O('Review z9: H4 rerun (h4_split.line) against h4_split.txt. Written by an AI model; not read by a person.')
ref = open(os.path.join(H4, 'h4_split.txt')).read().splitlines()
for m, n in H.RATIOS:
    l = H.line(m, n); r = [x for x in ref if x.startswith('%d,%d ' % (m, n))]
    O(l); O('    identical to the line of h4_split.txt: %s' % (r == [l]))
O.finish()
