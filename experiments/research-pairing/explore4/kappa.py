# fit the class constants of the mixed curves on the pilot's lines: class = (2, kappa * d_j d_k) at (e1, e2)
import sys, collections
sys.path.insert(0,'/tmp/t10'); import os; os.chdir('/tmp/t10')
from parse import load, sqf
R=load(); idx={'0':0,'8':1,'L':2,'U':3}; tab=collections.Counter(); c1s=collections.Counter()
for c in R:
    cv=c['cv']
    if cv in ('08','LU'): continue
    m,n=c['ratio']
    for l in c['lines']:
        D=l['tup'][idx[cv[0]]]*l['tup'][idx[cv[1]]]
        k=sqf(l['cls'][1]*D); base=sqf(m*n*(m*m-n*n))
        q=sqf(k*base)        # kappa / sqf(mn(m^2-n^2)) mod squares
        par='m even' if m%2==0 else 'n even' if n%2==0 else 'both odd'
        tab[(cv,par,q)]+=1; c1s[(cv,l['cls'][0])]+=1
for k,v in sorted(tab.items()): print(v,k)
print(dict(c1s))
