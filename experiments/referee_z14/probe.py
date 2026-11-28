import sys, os, time, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'conductor-t10'))
sys.argv = sys.argv[:1]
import t10_pairct as T
RS = T.RS
from fractions import Fraction as Fr
for (m, n) in [(465, 416), (881, 49), (1234, 567), (7, 4)]:
    A, B, C = m*m-n*n, 2*m*n, m*m+n*n
    Rp, Rm = C*C+A*B, C*C-A*B
    for j, k in [('0','L'),('0','U'),('8','L'),('8','U')]:
        t0=time.time()
        Pr = T.Pair(m, n, j, k)
        es = Pr.es
        print((m,n), j+k, 'base', Pr.t0, 'diffs ok', (es[1]-es[0], es[2]-es[0]) == (2*Rm, 2*Rp), 'u', Pr.u, '%.2fs'%(time.time()-t0))
es = (-2*C*C, -2*A*B, 2*A*B)
t0=time.time()
sh = es[1]; esh = tuple(e - sh for e in es)
bad, loc, sel, kum, enc, places = RS.full_selmer2(esh, random.Random(1))
print('sel', len(sel), sel, '%.2fs'%(time.time()-t0))
for a in sel[:3]:
    D = T.DFullX(esh, a)
    for b in sel[:3]:
        t0=time.time(); h, per = T.safe_pair(D, b); print(a, b, h, '%.2fs'%(time.time()-t0))
