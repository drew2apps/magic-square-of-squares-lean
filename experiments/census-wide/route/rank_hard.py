import sys, time
sys.path.insert(0, '..')
from sd_pari import *
pari.allocatemem(4*10**9)
m, n, d = map(int, sys.argv[1].split(',')); eff = int(sys.argv[2])
Q, r = family(m, n)['+13']; ain, A = weierstrass(d*Q, r)
Em = pari.ellminimalmodel(pari.ellinit([pari(x.numerator)/x.denominator for x in ain]))
for k, c in enumerate(pari.ellisomat(Em, 0, 1)[0]):
    t = time.time(); b = pari.ellrank(pari.ellinit(c), eff); print(m, n, d, 'curve', k, 'effort', eff, 'bounds', int(b[0]), int(b[1]), 'points', len(b[3]), '%.0f s' % (time.time() - t), flush=True)
