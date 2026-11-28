# Soundness check of the five-form local test: at random coprime (u,v), whenever no form is dead at p,
# the tuple of labels must be in the allowed set (the repo's own selftest, applied to five forms).
import sys; sys.argv=['x']; sys.path.insert(0,'/tmp/rn'); sys.path.insert(0,'/tmp/ghost')
from pipeline import *
from f4_form import residual5
import random, math
tot=[0,0]
for m,n in ((14903,9696),(9113,3200),(3320,2303),(465,416),(627,364),(4,3),(7,1)):
    forms,sup,known=TD.magic(m,n); c=TD.gmul((1,1),(m,n)); forms5=forms+[(c,TD.gmul(c,(0,1)))]
    bad=set()
    for x in (2*m*n*(m*m-n*n), m*m+n*n): bad|=set(fact(x))
    small=[p for p in sorted(bad) if p<200]
    h,f=TD.selftest(forms5,small,samples=4000,seed=5)
    # the exact allowed sets used by residual5, tested the same way
    h2=f2=0; rnd=random.Random(9)
    for p in small:
        L,al=X.allowed_at(forms5,p)
        for _ in range(3000):
            u,v=rnd.randrange(-10**6,10**6),rnd.randrange(1,10**6)
            if math.gcd(u,v)!=1: continue
            try: res=[L.cls((a[0]*u+b[0]*v,a[1]*u+b[1]*v),10**6) for a,b in forms5]
            except Exception as ex: h2=-1; break
            if any(r[0] or not r[1] for r in res): continue
            h2+=1; f2+=tuple(r[2] for r in res) not in al
    print((m,n),'primes',small,'repo Local: points',h,'failures',f,'| exact_local: points',h2,'failures',f2)
