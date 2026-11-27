"""One table from the per-range outputs g1_k13.txt ... g1_k29.txt (bookkeeping only).  usage: g1_summary.py [--write]"""
import os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__)); rows = []; T = [0]*5; bad = []
rows.append('k   range                         drawn  degen  no point  points everywhere  class kept     survivors  controls')
for k in range(13, 30):
    s = open(os.path.join(HERE, 'g1_k%d.txt' % k)).read()
    rng = re.search(r'range k = \d+: (\d+) <= m < (\d+); (\d+) random', s)
    dg = int(re.search(r'degenerate among the random ratios: (\d+)', s).group(1))
    npnt, pe = map(int, re.search(r'no local point somewhere: (\d+); a local point everywhere: (\d+)', s).groups())
    ck = int(re.search(r'a class kept by the criterion: (\d+)', s).group(1))
    sv = re.search(r'SURVIVORS[^:]*: (.*)', s).group(1)
    c1, c2 = map(int, re.search(r'controls \([^)]*\): (\d+), (\d+) keep', s).groups())
    bad += [l.split(' degenerate')[0] for l in s.splitlines() if 'NOT the degenerate points' in l or 'FAILED' in l]
    for i, v in enumerate((int(rng.group(3)), dg, npnt, pe, ck)): T[i] += v
    rows.append('%-3d %-29s %6s %6d %9d %18d %5d (%5.2f%%) %9s  %d/%d' % (k, '[%s, %s)' % rng.group(1, 2), rng.group(3), dg, npnt, pe, ck, 100.0*ck/pe, sv, c2, c1))
rows.append('all %-29s %6d %6d %9d %18d %5d (%5.2f%%)' % ('', T[0], T[1], T[2], T[3], T[4], 100.0*T[4]/T[3]))
rows.append('controls not keeping exactly their degenerate points: %s' % (', '.join(bad) or 'none'))
out = '\n'.join(rows) + '\n'
if '--write' in sys.argv: open(os.path.join(HERE, 'g1_summary.txt'), 'w').write(out)
print(out, end='')
