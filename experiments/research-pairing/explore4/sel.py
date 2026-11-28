"""local images of E(Q_p)/2 and Selmer linear algebra over F_2 (bit masks)."""
import sys,random
sys.path.insert(0,'/tmp/mix')
from cassels import *
def rank_reduce(rows):
    """row-reduce list of int bitmasks; returns basis (echelon, distinct leading bits)"""
    basis=[]
    for r in rows:
        for b in basis:
            if r ^ b < r: r^=b
        if r: basis.append(r); basis.sort(reverse=True)
    return basis
def inspan(basis,v):
    for b in basis:
        if v ^ b < v: v^=b
    return v==0
def bits(t):
    r=0
    for x in t: r=(r<<1)|x
    return r
def cls2(c1,c2,p): return bits(class_vec(c1,p)+class_vec(c2,p))
def local_image(e,p,rng,tries=200000):
    """basis of the image of E(Q_p) in (Q_p^*/sq)^2 (components at e1,e2); proves dimension 2 (p odd) / 3 (p=2)"""
    target=3 if p==2 else 2
    rows=[]
    for i in range(3):
        t=tors_image(e,i); rows.append(cls2(t[0],t[1],p))
    basis=rank_reduce(rows)
    vmax=max(vp(e[i]-e[j],p) for i in range(3) for j in range(i)); n=0
    while len(basis)<target and n<tries:
        n+=1
        x=sample_x(e,p,rng,vmax); nn,d=x.numerator,x.denominator
        f=[(nn-e[i]*d)*d for i in range(3)]
        if 0 in f or not issq_int(f[0]*f[1]*f[2],p): continue
        v=cls2(f[0],f[1],p)
        if not inspan(basis,v): basis=rank_reduce(basis+[v])
    if len(basis)!=target: raise RuntimeError('local image rank %d at %d'%(len(basis),p))
    return basis
def nullspace(vecs,test):
    """vecs: list of generators g; test(g)-> list of (basis,vector) local data.  Returns the subgroup (as list of
    index-masks over generators) of products whose local vector lies in the local image at every place."""
    # eliminate: maintain list of (combo mask, residual per place)
    n=len(vecs); cur=[(1<<i, None) for i in range(n)]
    data=[test(g) for g in vecs]          # data[i][place] = (basis, vec)
    nplaces=len(data[0])
    combos=[1<<i for i in range(n)]; res=[[d[k][1] for k in range(nplaces)] for d in data]
    for k in range(nplaces):
        basis=data[0][k][0]
        # reduce residual modulo local image
        def red(v):
            for b in basis:
                if v ^ b < v: v^=b
            return v
        r=[red(x[k]) for x in res]
        # gaussian elimination on r to find combos with zero residual at place k
        piv={}; newc=[]; newres=[]
        for i in range(len(combos)):
            v=r[i]; c=combos[i]; rr=list(res[i])
            while v:
                h=v.bit_length()
                if h in piv:
                    pv,pc,pr=piv[h]; v^=pv; c^=pc; rr=[x^y for x,y in zip(rr,pr)]
                else:
                    piv[h]=(v,c,rr); break
            if v==0: newc.append(c); newres.append(rr)
        combos,res=newc,newres
    return combos
