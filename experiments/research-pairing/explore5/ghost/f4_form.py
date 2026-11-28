"""3(b): add the centre as a fifth Gaussian linear form gamma_4 = (1+i)(m+in)(u+iv), of norm F4 = 2(m^2+n^2)(u^2+v^2),
to the repository's local test (which uses the four forms F0, F8, GL, GU only), and recompute the residual set at the
bad primes.  Support of d4: primes of 2(m^2+n^2) (det of gamma_4 and its conjugate is -4i(m^2+n^2))."""
import sys,json
sys.argv=['x']; sys.path.insert(0,'/tmp/rn')
from pipeline import *
gm=TD.gmul
def residual5(m,n,extra=()):
    forms,sup,known=TD.magic(m,n)
    c=gm((1,1),(m,n))
    f4=(c, gm(c,(0,1)))      # alpha u + beta v with alpha=c, beta=c*i
    forms5=forms+[f4]
    C=m*m+n*n
    sup5=sup+[sorted({2}|set(fact(C)))]
    divs=[]
    for S in sup5:
        ds=[1]
        for q in S: ds+=[x*q for x in ds]
        divs.append(ds)
    bad=set(extra)
    for x in (2*m*n*(m*m-n*n), C, m**8+14*m**4*n**4+n**8): bad|=set(fact(x))
    loc=[]; 
    for p in sorted(bad):
        L,al=X.allowed_at(forms5,p)
        if len(al)<len(L.labels())**5: loc.append((p,L,al))
    return sorted(X.join(divs,loc)), forms5
if __name__=='__main__':
    # sanity: degenerate ratios keep the known solutions (with their d4)
    for m,n in ((4,3),(12,5),(7,1),(17,7),(21,20),(41,1)):
        T5,forms5=residual5(m,n)
        kn=[]
        C=m*m+n*n
        import math
        pts=[(1,1),(-1,1)] if math.isqrt(C)**2==C else [(0,1),(1,0)]
        for u,v in pts:
            kn.append(tuple(TD.cls_of((a[0]*u+b[0]*v,a[1]*u+b[1]*v)) for a,b in forms5))
        print('degenerate',(m,n),'5-tuples left',len(T5),'known',kn,'all known kept:',all(k in T5 for k in kn))
    for m,n in ((3320,2303),(5623,1017),(9113,3200),(12313,5913),(14903,9696),(24599,5207)):
        T4,_=residual(m,n); T5,_=residual5(m,n)
        proj=sorted(set(t[:4] for t in T5))
        print((m,n),'4-form residual',len(T4),'| 5-form residual',len(T5),'| 4-tuples that extend',len(proj))
        for t in T5: print('    ',t)
