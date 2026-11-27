#!/usr/bin/env python3
"""z11: H1's control with my code: (1109503,597671) (twice-square type), its 2 known and 4 extra tuples on E_0L, E_0U,
with my base points t0 = m/n, -m/n (H1 and t10 use t = 0).  The extra tuples are locally soluble everywhere (H1, review
z8); the expected kills are read from the VERDICT lines of ../roadmap-h/h1/h1_pairs.txt.  Written by an AI model
(review agent z11); not read by a person.  PARI as in z11_pair.py (no ellrank here).
    python3 -I h1_mine.py      (prints h1_mine.txt)"""
import sys, os, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _z11common import H1, Out
from z11_pair import analyse
say = Out('h1_mine.txt')
say('h1_mine.py (review z11; written by an AI model): (1109503,597671) against h1_pairs.txt')
want = {}
for l in open(os.path.join(H1, 'h1_pairs.txt')):
    m = re.match(r'  (\(\d+, \d+, \d+, \d+\)) (known|extra).*?(killed by E_(\w\w)|not excluded)', l)
    if m: want[eval(m.group(1))] = ({m.group(4)} if m.group(4) else set(), m.group(2) == 'known')
KNOWN = sorted(t for t, (k, kn) in want.items() if kn); EXTRA = sorted(t for t, (k, kn) in want.items() if not kn)
say('h1_pairs.txt VERDICT: %s' % '; '.join('%s %s' % (t, ','.join(sorted(want[t][0])) or 'not excluded') for t in sorted(want)))
out, kill, known = analyse(1109503, 597671, tuples=(KNOWN + EXTRA, KNOWN), pairs=['0L', '0U'], points=False, log=say)
ok = len(want) == 6 and all({jk for jk in ('0L', '0U') if out[t][jk][i][1][0] in ('killed', 'outside')} == want[t][0]
                            for t in want for i in (0, 1))
say('MY CODE AGREES WITH h1_pairs.txt (both signs)' if ok else 'MY CODE DIFFERS FROM h1_pairs.txt')
say.finish()
