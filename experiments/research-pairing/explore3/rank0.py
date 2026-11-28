"""How often does a 2-descent with the Cassels-Tate pairing prove rank 0 for E_A: y^2=x(x-A^2)(x-C^2) or its twist by -1, y^2=x(x-B^2)(x-C^2)?
Rank 0 on either curve excludes the ratio: the 8 torsion points lie over t in {0, +-1, infinity}."""
import sys, random, collections, time
sys.path.insert(0,'/tmp/mix')
from math import gcd, isqrt
from cassels import *
from sel import *
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
def bound(e,seed):
    rng=random.Random(seed)
    S=sorted({2}|set().union(*[primes_of(abs(e[i]-e[j])) for i in range(3) for j in range(i)]))
    LI={p:local_image(e,p,rng) for p in S}
    g2=[-1]+S; gen_pairs=[(g,1) for g in g2]+[(1,g) for g in g2]
    def t_full(gp):
        out=[(LI[p],cls2(gp[0],gp[1],p)) for p in S]; out.append(([],1 if gp[0]<0 else 0)); return out
    ker=nullspace(gen_pairs,t_full)
    basis=[]
    for k in ker:
        c1=c2=1
        for i,gp in enumerate(gen_pairs):
            if k>>i&1: c1*=gp[0]; c2*=gp[1]
        basis.append((sqf(c1),sqf(c2)))
    d=len(basis)
    G=[[0]*d for _ in range(d)]
    for i in range(d):
        for j in range(i+1,d):
            v=cassels(e,triple(basis[i]),triple(basis[j]),seed=seed)[0]
            G[i][j]=G[j][i]=0 if v==1 else 1
    return d, d-rank2(G)-2
def is_sq(x): r=isqrt(x); return r*r==x
rng=random.Random(5); H=int(sys.argv[1]); N=int(sys.argv[2]); kind=sys.argv[3]
pairs=[]
if kind=='nd':
    while len(pairs)<N:
        m=rng.randint(2,H); n=rng.randint(1,m-1)
        if gcd(m,n)!=1 or not ((m*n)%4==0 or (m-n)%8==0 or (m+n)%8==0): continue
        S=m*m+n*n
        if is_sq(S) or (S%2==0 and is_sq(S//2)): continue
        pairs.append((m,n))
else:
    for a in range(2,80):
        for b in range(1,a):
            if gcd(a,b)==1 and (a-b)%2==1:
                m,n=max(a*a-b*b,2*a*b),min(a*a-b*b,2*a*b)
                if (m*n)%4==0 and m<=H: pairs.append((m,n))
    rng.shuffle(pairs); pairs=pairs[:N]
out=collections.Counter(); t0=time.time(); fails=0
for (m,n) in pairs:
    A,B,C=m*m-n*n,2*m*n,m*m+n*n
    try:
        dA,bA=bound((0,A*A,C*C),7); dB,bB=bound((0,B*B,C*C),7)
    except Exception as ex:
        fails+=1; continue
    out[('bound on E_A',bA)]+=1; out[('bound on twist',bB)]+=1
    out[('min of the two bounds',min(bA,bB))]+=1
    out[('parity of the two bounds',(bA%2,bB%2))]+=1
print(kind,'ratios',len(pairs),'failures',fails,'time %.0f s'%(time.time()-t0))
for k,v in sorted(out.items()): print('  ',v,k)
