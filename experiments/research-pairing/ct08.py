"""Cassels' pairing on E: y^2 = x(x-A^2)(x-C^2), C^2-A^2=B^2, for a class a=(1,delta,delta), with the explicit conic points
(A,0,1), (C,0,1), (1,1,0).  <a,b> = prod_v (u-A, b_*)_v (u-C, b_m)_v (u_m-u_*, b_r)_v  at a local point  u^2=x, delta u_m^2 = x-A^2, delta u_*^2 = x-C^2."""
import cypari2, random, sys, collections
from fractions import Fraction as Fr
from sympy import factorint
from parse import load, sqf
pari=cypari2.Pari(); pari.allocatemem(10**9)
def val(x,p):
    x=Fr(x); n,d=x.numerator,x.denominator; v=0
    while n%p==0: n//=p; v+=1
    while d%p==0: d//=p; v-=1
    return v,n,d
def issq(x,p):
    if x==0: return False
    v,n,d=val(x,p)
    if v%2: return False
    u=n*pow(d,-1,p**3 if p==2 else p)%(8 if p==2 else p)
    return u%8==1 if p==2 else pow(u,(p-1)//2,p)==1
def hil(x,y,p):
    return int(pari.hilbert(x,y,p))
def local_point(A,C,dl,p,rng):
    ok=lambda u: issq((u*u-A*A)/Fr(dl),p) and issq((u*u-C*C)/Fr(dl),p)
    cands=[Fr(0)]+[Fr(1,p**(2*k)) *w for k in (3,4) for w in (1,3,5,7)]
    for u in cands:
        if u!=0 and not issq(Fr(dl),p): continue
        if ok(u): return u
    for it in range(20000):
        j=rng.randint(0,12); w=rng.randint(1,p**3 if p<50 else p*p); base=rng.choice([A,-A,C,-C,0])
        u=Fr(base)+Fr(w)*Fr(p)**rng.choice([j,-j]) if base==0 else Fr(base)+w*p**j
        if u!=0 and ok(u): return u
    return None
def psqrt(x,p,N=60):
    return pari(f'sqrt(({x.numerator})/({x.denominator})+O({p}^{N}))')
def pairing(A,C,dl,b_r,b_s,seed=1,conv=0,detail=False):
    b_m=sqf(b_r*b_s); rng=random.Random(seed)
    B2=C*C-A*A
    primes=sorted(set([2])|set(factorint(abs(A)))|set(factorint(B2))|set(factorint(abs(C)))|set(factorint(abs(dl)))|set(factorint(abs(b_r)))|set(factorint(abs(b_s))))
    tot=1; loc={}
    for p in primes:
        u=local_point(A,C,dl,p,rng)
        if u is None: return None,{}
        t=1
        if u-A==0 or u-C==0: u=u+Fr(p)**40  # stay off the special points (they remain local points to this precision? avoid)
        x1,x2=(b_s,b_m) if conv==0 else (b_m,b_s)
        t*=hil(pari(f'({(u-A).numerator})/({(u-A).denominator})'),x1,p)*hil(pari(f'({(u-C).numerator})/({(u-C).denominator})'),x2,p)
        if b_r!=1:
            um=psqrt((u*u-A*A)/Fr(dl),p); us=psqrt((u*u-C*C)/Fr(dl),p); d=um-us
            if pari.valuation(d,p)>40: d=um+us
            t*=hil(d,b_r,p)
        loc[p]=t; tot*=t
    return tot,loc
if __name__=='__main__':
    R=load(); idx={'0':0,'8':1,'L':2,'U':3}; res=collections.Counter(); bad=[]
    which=sys.argv[1] if len(sys.argv)>1 else '08'
    for c in R:
        if c['cv']!=which: continue
        m,n=c['ratio']; S=m*m+n*n
        for l in c['lines']:
            D=l['tup'][idx[which[0]]]*l['tup'][idx[which[1]]]
            if which=='08': A,C=m*m-n*n,S; dl=sqf(m*n*D); assert l['cls'][0]==1 and sqf(l['cls'][1])==dl
            else:
                A,C=8*m*n,4*S; dl=sqf((m*m-n*n)*D)     # c = (2S, 2D) at (e_r, e_m); times the image (2S, 2(m^2-n^2)) of a point of order 4
                assert sqf(l['cls'][0])==sqf(2*S) and sqf(l['cls'][1])==sqf(2*D)
            b_r,b_m=l['b']; b_s=sqf(b_r*b_m)           # printed components are at (e_r, e_m)
            for conv in (0,):
                v,loc=pairing(A,C,dl,b_r,b_s,conv=conv)
                res[(conv,v,'b_r=1' if b_r==1 else 'b_r!=1')]+=1
                if conv==0 and v!=-1 and len(bad)<5: bad.append((m,n,l['tup'],dl,l['b'],v,loc))
    print(which,dict(res)); print(bad[:3])
