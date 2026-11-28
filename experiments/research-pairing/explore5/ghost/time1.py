import sys,time,json,random
sys.argv=['x']; sys.path.insert(0,'/tmp/rn')
from pipeline import *
from math import gcd
rng=random.Random(5); ts=[]; hits=0
for _ in range(40):
    while True:
        m=rng.randrange(15001,40000); n=rng.randrange(1,m)
        if gcd(m,n)==1: break
    t0=time.time(); T,known=residual(m,n); ts.append(time.time()-t0); hits+=bool(T)
    print(m,n,len(T),round(ts[-1],2),flush=True)
print(sum(ts)/len(ts),hits)
