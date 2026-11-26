"""Harder PARI rank bounds for dm_pari's D1 curve on a few ratios. usage: d1rank.py m,n effort.  EXPERIMENT."""
import sys, time
sys.path.insert(0, '../../census2000')
from dm_pari import *
pari.allocatemem(4*10**9)
m, n = map(int, sys.argv[1].split(',')); eff = int(sys.argv[2])
c = [r for r in all_curves(m, n) if r[4] == 'D1'][0]
Em = pari.ellminimalmodel(pari.ellinit(c[0]))
print(m, n, 'D1 conductor', pari.ellglobalred(Em)[0], 'root number', pari.ellrootno(Em), flush=True)
for k, a in enumerate(pari.ellisomat(Em, 0, 1)[0]):
    t = time.time(); b = pari.ellrank(pari.ellinit(a), eff)
    print(m, n, 'D1 curve', k, 'effort', eff, 'bounds', int(b[0]), int(b[1]), 'points', len(b[3]), '%.0f s' % (time.time()-t), flush=True)
