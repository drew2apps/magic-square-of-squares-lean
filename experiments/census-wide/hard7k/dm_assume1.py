"""dm_pari's Demjanenko-Manin step on one curve, ASSUMING rank 1 (PARI bounds only [1,3]).  Conditional result.
usage: dm_assume1.py m,n curve_name.  EXPERIMENT."""
import sys, json
sys.path.insert(0, '../../census2000')
import dm_pari
from dm_pari import *
pari.allocatemem(4*10**9)
m, n = map(int, sys.argv[1].split(',')); name = sys.argv[2]
real = dm_pari.rank_bounds
dm_pari.rank_bounds = lambda E: (1, 1)
c = [r for r in all_curves(m, n) if r[4] == name][0]
r = dm_curve(m, n, *c)
r['rank'] = list(real(pari.ellinit(c[0]))); r['assumed'] = 'rank 1'
print(json.dumps(r, default=str))
