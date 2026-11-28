import random, itertools, collections
from fractions import Fraction as Fr
from sympy import factorint
from parse import load, sqf
from ct08 import issq
def soluble(e,a,p,rng,tries=4000):
    # a_i u_i^2 = x - e_i ; search x
    for it in range(tries):
        j=rng.randint(0,10); w=rng.randint(1,p**3 if p<60 else p*p)
        base=rng.choice(list(e)+[None])
        x=Fr(w,p**(2*j)) if base is None else Fr(base)+w*p**j
        if all(x!=ei and issq((x-ei)/Fr(ai),p) for ei,ai in zip(e,a)): return True
    return False
R=load(); rng=random.Random(3); cnt=collections.Counter()
for which in ('LU','08','0L'):
    k=0
    for c in R:
        if c['cv']!=which or not c['lines']: continue
        k+=1
        if k>6: break
        l=c['lines'][0]; e=c['e']; c1,c2=l['cls']
        m,n=c['ratio']
        ps=sorted(set(factorint(abs((e[0]-e[1])*(e[0]-e[2])*(e[1]-e[2])*c1*c2))))
        for (i,j) in itertools.permutations(range(3),2):
            a=[None]*3; a[i]=c1; a[j]=c2; a[3-i-j]=sqf(c1*c2)
            ok=all(soluble(e,a,p,rng) for p in ps)
            cnt[(which,'c1 at e%d, c2 at e%d'%(i+1,j+1),ok)]+=1
for k,v in sorted(cnt.items()): print(k,v)
