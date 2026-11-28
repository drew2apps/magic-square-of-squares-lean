import sys,time,json
sys.argv=['x']; sys.path.insert(0,'/tmp/rn')
from pipeline import *
for m,n in ((24599,5207),(14903,9696)):
    t0=time.time()
    out,dims,known=run(m,n)
    print((m,n),'C=',fact(m*m+n*n),'ntuples',len(out),'dims',dims,'time',round(time.time()-t0,1))
    for t,d in out.items():
        a=sqf(m*n*t[0]*t[1]); b=sqf((m*m-n*n)*t[2]*t[3])
        gh=(a,b)==(sqf(2*m*n),1) or (a,b)==(1,sqf(2*(m*m-n*n)))
        print('  ',t,'ghost' if gh else 'other', 'SURVIVES' if not any(d.values()) else 'killed on '+','.join(k for k in d if d[k]))
