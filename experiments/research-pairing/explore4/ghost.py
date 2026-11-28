import os,sys
os.environ['MSS_REPO']='/tmp/rn/repo'; sys.path.insert(0,'/tmp/rn/repo/scripts/cover_descent')
import tuple_descent as TD
from sympy import factorint
def sqf(x):
    r=1
    for q,e in factorint(abs(x)).items():
        if e%2: r*=q
    return r
for m,n in ((4,3),(12,5),(15,8),(21,20),(35,12),(45,28),(63,16),(77,36),(7,1),(17,7),(23,7),(31,17),(41,1),(49,23)):
    forms,sup,known=TD.magic(m,n)
    out=[]
    for t in known:
        a=sqf(m*n*t[0]*t[1]); b=sqf((m*m-n*n)*t[2]*t[3])
        s='S1' if (a,b)==(sqf(2*m*n),1) else 'S2' if (a,b)==(1,sqf(2*(m*m-n*n))) else 'other (a,b)=%s'%((a,b),)
        out.append((t,s))
    print((m,n),'S=',m*m+n*n,out)
