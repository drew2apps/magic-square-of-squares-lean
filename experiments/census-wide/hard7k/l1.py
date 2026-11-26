"""L(E,1) for a quotient curve by PARI's lfun.  If L(E,1) != 0 then rank E(Q) = 0 (Gross-Zagier-Kolyvagin; this is a theorem,
but the computation is numerical).  usage: l1.py m,n,name.  EXPERIMENT."""
import sys, time
R = '/tmp/claude-0/-home-user/9ec7a8d3-42dc-53c8-8518-1a663ac809f1/scratchpad/repo/scripts/g1_full'
import importlib.util
spec = importlib.util.spec_from_file_location('g', R + '/global.py'); g = importlib.util.module_from_spec(spec); spec.loader.exec_module(g)
pari = g.pari; pari.allocatemem(24*10**9); pari.set_real_precision(19)
m, n, name = sys.argv[1].split(','); m, n = int(m), int(n)
poly = g.curves(m, n)[name]
c = [poly[4], poly[3], poly[2], poly[1]] if len(poly) == 5 else poly
c0, c1, c2, c3 = c
E = pari.ellminimalmodel(pari.ellinit([0, c2, 0, c1*c3, c0*c3*c3]))
t = time.time()
L = pari.ellL1(E, 0)
print(m, n, name, 'L(E,1) =', L, 'seconds', round(time.time()-t), flush=True)
t = time.time(); print('BSD-predicted Sha (if rank 0):', pari.ellbsd(E), 'omega', flush=True)
