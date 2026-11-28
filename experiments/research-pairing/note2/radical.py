"""How big is the radical of Theta on S' ?  Random coprime ratios with a 2-adic point; factorisations only."""
import random, collections, sys
from math import gcd, isqrt
from isogeny_part import Sprime, Theta
from closed_formula import triple
def is_sq(x): r=isqrt(x); return r*r==x
def basis(SP):
    B=[]; span={1}
    for b in sorted(SP,key=abs):
        if b in span: continue
        B.append(b); span|={sqf_mul(b,s) for s in span}
    return B
from closed_formula import sqf
def sqf_mul(a,b): return sqf(a*b)
def rank2(M):
    M=[r[:] for r in M]; rk=0; n=len(M)
    for c in range(n):
        piv=next((i for i in range(rk,n) if M[i][c]),None)
        if piv is None: continue
        M[rk],M[piv]=M[piv],M[rk]
        for i in range(n):
            if i!=rk and M[i][c]: M[i]=[x^y for x,y in zip(M[i],M[rk])]
        rk+=1
    return rk
def stats(pairs,label):
    out=collections.Counter()
    for m,n in pairs:
        for cv in ('08','LU'):
            A,C=triple(m,n,cv); SP=Sprime(A,C); B=basis(SP)
            G=[[0 if Theta(A,C,x,y)==1 else 1 for y in B] for x in B]
            out[(label,cv,'dim radical of Theta = %d'%(len(B)-rank2(G)))]+=1
    return out
rng=random.Random(11); nd=[]; H=int(sys.argv[1]) if len(sys.argv)>1 else 20000
while len(nd)<1500:
    m=rng.randint(2,H); n=rng.randint(1,m-1)
    if gcd(m,n)!=1: continue
    if not ((m*n)%4==0 or (m-n)%8==0 or (m+n)%8==0): continue
    S=m*m+n*n
    if is_sq(S) or (S%2==0 and is_sq(S//2)): continue
    nd.append((m,n))
dg=[]
for a in range(2,60):
    for b in range(1,a):
        if gcd(a,b)==1 and (a-b)%2==1:
            m,n=max(a*a-b*b,2*a*b),min(a*a-b*b,2*a*b)
            if (m*n)%4==0: dg.append((m,n))
dg=dg[:400]
res=stats(nd,'non-degenerate, random, m <= %d'%H); res.update(stats(dg,'square type'))
for k,v in sorted(res.items()): print(v,k)
