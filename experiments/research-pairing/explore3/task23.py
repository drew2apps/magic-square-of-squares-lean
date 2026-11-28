import sys,time,random,collections,json,itertools
sys.path.insert(0,'/tmp/mix'); sys.path.insert(0,'/tmp/t10')
from parse import load
from cassels import *
from sel import *
import mpmath
R=load(); T0=time.time()
SEED=int(sys.argv[1]) if len(sys.argv)>1 else 11
NR=int(sys.argv[2]) if len(sys.argv)>2 else 60
PER=int(sys.argv[3]) if len(sys.argv)>3 else 50
EXTRA=int(sys.argv[4]) if len(sys.argv)>4 else 150
rng=random.Random(SEED)
curves={}
for c in R:
    if c['cv']=='0L': curves[c['ratio']]=c
ratios=sorted(curves); sample=sorted(rng.sample(ratios,NR))
def conic_pairing(m,n,d,beta,seed=1,detail=False):
    """prod_v (ell(t^2,z^2,w), beta)_v, ell = tangent form of w^2 = d T^2 + 8C^2 T Z + (16A^2B^2/d) Z^2 at a qfsolve point,
    (t,z=1,w) a local point of w^2 = d t^4 + 8 C^2 t^2 + 16A^2B^2/d found by an independent search."""
    rg=random.Random(seed)
    A=m*m-n*n;B=2*m*n;C=m*m+n*n; W=(C*C+A*B)*(C*C-A*B)
    c0=(16*A*A*B*B)//d; assert c0*d==16*A*A*B*B
    G=pari.matrix(3,3,[d,4*C*C,0,4*C*C,c0,0,0,0,-1])
    s=pari.qfsolve(G)
    if s.type()!='t_COL': return None,{'fail':'conic insoluble'}
    den=1
    for x in s: den=den*int(x.denominator())//gcd(den,int(x.denominator()))
    pt=[int(x*den) for x in s]; g=gcd(gcd(pt[0],pt[1]),pt[2]); T0_,Z0_,w0_=[x//g for x in pt]
    assert d*T0_**2+8*C*C*T0_*Z0_+c0*Z0_**2==w0_**2
    L=[d*T0_+4*C*C*Z0_, 4*C*C*T0_+c0*Z0_, -w0_]; g=gcd(gcd(L[0],L[1]),L[2]); L=[x//g for x in L]
    P={2}|primes_of(A)|primes_of(B)|primes_of(C*C+A*B)|primes_of(C*C-A*B)|primes_of(d)|primes_of(beta)
    for x in L: P|=primes_of(x)
    q=lambda t: d*t**4+8*C*C*t*t+c0
    tot=1; loc={}
    # real place
    if d>0: t=Fr(rg.randrange(1,10**6),rg.randrange(1,10**3))
    else:
        t0=mpmath.sqrt(mpmath.mpf(4*C*C)/abs(d)); t=Fr(int(t0*10**6),10**6)
    qv=q(t)
    if qv<=0: return None,{'fail':'real'}
    w=mpmath.sqrt(mpmath.mpf(qv.numerator)/qv.denominator)*rg.choice([1,-1])
    val=L[0]*mpmath.mpf(t.numerator)**2/t.denominator**2+L[1]+L[2]*w
    if abs(val)<mpmath.mpf(10)**(-60)*(abs(L[1])+abs(L[2]*w)): return None,{'fail':'real cancel'}
    tv=-1 if (val<0 and beta<0) else 1
    loc[-1]=tv; tot*=tv
    for p in sorted(P):
        N=prec_for(p); vmax=vp(16*A*A*B*B,p); done=False
        for attempt in range(60000):
            if p==2: u=rg.randrange(1,2**12)
            else: u=rg.randrange(1,min(p*p,10**12))
            j=rg.randint(-(vmax+4),vmax+4)
            t=Fr(u)*Fr(p)**j; qv=q(t)
            if qv==0 or not issq_int(qv.numerator*qv.denominator,p): continue
            w=psqrt(qv,p,N)*rg.choice([1,-1])
            val=L[0]*pari('(%d)/(%d)'%((t*t).numerator,(t*t).denominator))+L[1]+L[2]*w
            if val==0 or int(val.padicprec(p))-int(pari.valuation(val,p))<4: continue
            tv=int(pari.hilbert(val,beta,p)); done=True; break
        if not done: return None,{'fail':('no local point on quartic',p)}
        loc[p]=tv; tot*=tv
    return tot,(loc if detail else {})
def kron_all(x,plist):
    r=1
    for p in plist: r*=int(pari.kronecker(x,p))
    return r
def redei_guess(d,beta,Wodd):
    """solve x^2-d y^2=beta z^2 (qfsolve); for p | W to odd power: Legendre(x0+y0 s) with s^2=d mod p.
    returns dict or None if insoluble"""
    s=pari.qfsolve(pari.matdiagonal([1,-d,-beta]))
    if s.type()!='t_COL': return None
    den=1
    for x in s: den=den*int(x.denominator())//gcd(den,int(x.denominator()))
    pt=[int(x*den) for x in s]; g=gcd(gcd(pt[0],pt[1]),pt[2]); x0,y0,z0=[x//g for x in pt]
    assert x0*x0-d*y0*y0==beta*z0*z0
    out=dict(sol=(x0,y0,z0),dsq=True,bsq=True,zero=False,signdep=False)
    prod1=1; prod2=1
    for p in Wodd:
        if pari.kronecker(d,p)!=1: out['dsq']=False; return out
        if pari.kronecker(beta,p)!=1: out['bsq']=False
        sp=int(pari('lift(sqrt(Mod(%d,%d)))'%(d%p,p)))
        v1=int(pari.kronecker((x0+y0*sp)%p,p)); v2=int(pari.kronecker((x0-y0*sp)%p,p))
        if v1==0 or v2==0:
            out['zero']=True
            v1=v2=(v1 or v2)
        if v1!=v2: out['signdep']=True
        prod1*=v1; prod2*=v2
    out['S']=prod1; out['S2']=prod2
    # other primitive solutions through qfparam: is S independent of the solution?
    G=pari.matdiagonal([1,-d,-beta]); M=pari.qfparam(G,s)
    vals=set(); nsol=0
    for (u,v) in [(1,0),(0,1),(1,1),(1,-1),(2,1),(1,2),(3,1),(1,3),(2,3),(3,2),(5,2),(2,5),(7,3),(3,-7)]:
        q=[int(M[i,0])*u*u+int(M[i,1])*u*v+int(M[i,2])*v*v for i in range(3)]
        g=gcd(gcd(q[0],q[1]),q[2])
        if g==0: continue
        q=[x//g for x in q]; assert q[0]**2-d*q[1]**2==beta*q[2]**2
        pr=1
        for p in Wodd:
            sp=int(pari('lift(sqrt(Mod(%d,%d)))'%(d%p,p)))
            v1=int(pari.kronecker((q[0]+q[1]*sp)%p,p))
            if v1==0: v1=int(pari.kronecker((q[0]-q[1]*sp)%p,p))
            pr*=v1
        if pr!=0: vals.add(pr); nsol+=1
    out['Svals']=sorted(vals); out['nsol']=nsol
    return out
# symbolic identity check
import sympy
A_,B_,C_,d_,T_,Z_=sympy.symbols('A B C d T Z')
Rp_=C_**2+A_*B_; Rm_=C_**2-A_*B_; s_=d_*T_+4*A_*B_*Z_; D_=d_*T_-4*A_*B_*Z_
ident2=sympy.expand(Rp_*s_**2-Rm_*D_**2-2*A_*B_*d_*(d_*T_**2+8*C_**2*T_*Z_+16*A_**2*B_**2/d_*Z_**2))
print('symbolic identity R+ s^2 - R- D^2 - 2ABd*(dT^2+8C^2TZ+(16A^2B^2/d)Z^2) == 0 :',ident2==0)
res=dict(normdim=collections.Counter(),dims=collections.Counter(),sel2ok=collections.Counter(),t2=collections.Counter(),t3=collections.Counter(),pairs=0)
rows=[]   # per pair records for task 3 analysis
bad=[]
npairs=0
for (m,n) in sample:
    c=curves[(m,n)]; e=c['e']
    A=m*m-n*n;B=2*m*n;C=m*m+n*n; Rp=C*C+A*B; Rm=C*C-A*B; W=Rp*Rm
    assert e[2]-e[0]==2*Rp and e[1]-e[0]==2*Rm and e[2]-e[1]==4*A*B
    fW=pari.factor(W); Wodd=[int(p) for p,k in zip(fW[0],fW[1]) if int(k)%2 and int(p)!=2]
    S=sorted({2}|primes_of(A)|primes_of(B)|primes_of(C)|primes_of(Rp)|primes_of(Rm))
    lrng=random.Random(1000+m)
    LI={p:local_image(e,p,lrng) for p in S}
    gens=[-1]+sorted(primes_of(2*A*B))
    # (1,d,d) at (e1,e2,e3) -> components (c1,c2)=(1,d); real place: always soluble (sign pattern (+,+,+) or (+,-,-))
    ker=nullspace(gens,lambda g:[(LI[p],cls2(1,g,p)) for p in S])
    def elt(mask):
        r=1
        for i,g in enumerate(gens):
            if mask>>i&1: r*=g
        return r
    basis=[elt(k) for k in ker]; k=len(basis); res['dims'][k]+=1
    # validation of the local-image machinery: full Sel^2 dimension vs the data file
    g2=[-1]+sorted({2}|primes_of(A)|primes_of(B)|primes_of(Rp)|primes_of(Rm))
    gen_pairs=[(g,1) for g in g2]+[(1,g) for g in g2]
    def t_full(gp):
        out=[(LI[p],cls2(gp[0],gp[1],p)) for p in S]
        # real place: c1>0 required (image of E(R): (+,+,+),(+,-,-)): encode as residual bit
        out.append(([],1 if gp[0]<0 else 0)); return out
    kfull=nullspace(gen_pairs,t_full)
    res['sel2ok'][len(kfull)==c['dim']]+=1
    # all elements
    elts=[]
    for mask in range(1,1<<k):
        r=1
        for i in range(k):
            if mask>>i&1: r=sqf(r*basis[i])
        elts.append(r)
    pairs=[(x,y) for x in basis for y in basis]
    allp=[(x,y) for x in elts for y in elts]
    lrng.shuffle(allp)
    for pr in allp:
        if len(pairs)>=PER: break
        if pr not in pairs: pairs.append(pr)
    pairs=pairs[:PER]
    print('ratio',(m,n),'gens',len(gens),'dim of {d: (1,d,d) ELS} =',k,'basis',basis,'| Sel2 dim (mine,file)',len(kfull),c['dim'],'| pairs',len(pairs),'t=%.0f'%(time.time()-T0),flush=True)
    for (d,beta) in pairs:
        npairs+=1
        a=(1,d,d); b=(1,beta,beta)
        v1,i1=cassels(e,a,b,seed=3+SEED,detail=True)
        v2,i2=conic_pairing(m,n,d,beta,seed=5+7*SEED,detail=True)
        key=('agree' if (v1 is not None and v1==v2) else 'DISAGREE' if (v1 is not None and v2 is not None) else 'fail')
        res['t2'][key]+=1; res['t2'][('cassels value',v1)]+=1
        if key!='agree' and len(bad)<20: bad.append(((m,n),d,beta,v1,v2,i1 if v1 is None else '',i2 if v2 is None else ''))
        # symmetry / alternation
        if d==beta: res['t2'][('diag <d,d>',v1)]+=1
        # Task 3
        rg=redei_guess(d,beta,Wodd)
        if rg is None: res['t3']['x^2-dy^2=beta z^2 insoluble']+=1; res['t3'][('insoluble, pairing',v1)]+=1; continue
        if not rg['dsq']: res['t3']['d not a square mod some p|W']+=1; continue
        res['t3']['soluble']+=1
        if not rg['bsq']: res['t3']['beta nonsquare mod some p|W']+=1
        if rg['signdep']: res['t3']['depends on sign of s']+=1
        if rg['zero']: res['t3']['x0+y0 s = 0 mod p for some p']+=1
        x0,y0,z0=rg['sol']
        Jm1=kron_all(-1,Wodd); J2=kron_all(2,Wodd)
        rows.append(dict(ratio=(m,n),d=d,beta=beta,P=v1,S=rg['S'],S2=rg['S2'],Svals=rg['Svals'],nsol=rg['nsol'],main=True,Jm1=Jm1,J2=J2,Jz=kron_all(z0,Wodd),Jx=kron_all(x0,Wodd),Jy=kron_all(y0,Wodd),
                         sol=(x0,y0,z0),nW=len(Wodd),Wmod8=W%8,par=(m+n)%2,trivial=(d==1 or beta==1 or d==beta)))
    # norm group: d over gens with (d,W)_v=1 for all v  (d a norm from Q(sqrt W))
    PW=sorted({2}|primes_of(A)|primes_of(B)|set(int(p) for p in fW[0]))
    kn=nullspace(gens,lambda g:[([],bits([0 if int(pari.hilbert(g,W,p))==1 else 1 for p in PW]))])
    res['normdim'][(k,len(kn))]+=1
    assert all(inspan(rank_reduce(kn),x) for x in ker)
    # extra Task-3 pairs: random nontrivial pairs with x^2-dy^2=beta z^2 soluble (pairing from Task 1 code only)
    extra=0
    for (d,beta) in allp:
        if extra>=EXTRA: break
        if (d,beta) in pairs or d==beta or d==1 or beta==1: continue
        rg=redei_guess(d,beta,Wodd)
        if rg is None or not rg['dsq']: continue
        v1,_=cassels(e,(1,d,d),(1,beta,beta),seed=3)
        extra+=1; x0,y0,z0=rg['sol']
        rows.append(dict(ratio=(m,n),d=d,beta=beta,P=v1,S=rg['S'],S2=rg['S2'],Svals=rg['Svals'],nsol=rg['nsol'],main=False,Jm1=kron_all(-1,Wodd),J2=kron_all(2,Wodd),Jz=kron_all(z0,Wodd),Jx=kron_all(x0,Wodd),Jy=kron_all(y0,Wodd),
                         sol=(x0,y0,z0),nW=len(Wodd),Wmod8=W%8,par=(m+n)%2,trivial=False,signdep=rg['signdep'],bsq=rg['bsq']))
res['pairs']=npairs
print('\n=== TASK 2 ===')
print('ratios',len(sample),'dimension histogram of {d | (1,d,d) in Sel^2}:',dict(sorted(res['dims'].items())))
print('full Sel^2 dimension from my local images equals dim in data file:',dict(res['sel2ok']))
print('pairs',npairs,dict((str(k),v) for k,v in res['t2'].items())); print('bad',bad)
print('\n=== TASK 3 ===')
print(dict(res['t3']))
print('(dim of ELS set, dim of {d: norm from Q(sqrt W)}) histogram:',dict(sorted(res['normdim'].items())))
ALL=rows
for label,rows in [('pairs of the Task-2 sample',[r for r in ALL if r['main']]),('Task-2 sample, nontrivial (d!=1,beta!=1,d!=beta)',[r for r in ALL if r['main'] and not r['trivial']]),('extra random nontrivial soluble pairs',[r for r in ALL if not r['main']]),('ALL',ALL)]:
  print('---',label)
  n=len(rows)
  def cntm(f): return sum(1 for r in rows if f(r))
  print(' S takes more than one value over the tested primitive solutions:',cntm(lambda r:len(r['Svals'])>1),'/',n,'; min #solutions tested',min([r['nsol'] for r in rows] or [0]))
  print(' pairing = S for every tested solution:',cntm(lambda r:r['Svals']==[r['P']]),'/',n)
  print('pairs with a solution and d square mod all p|W (odd mult):',n)
  print(' pairing = S (plain):',cntm(lambda r:r['P']==r['S']),'/',n)
  print(' pairing = S2 (other sign of s):',cntm(lambda r:r['P']==r['S2']),'/',n)
  print(' pairing = -S :',cntm(lambda r:r['P']==-r['S']))
  for nm,f in [('S*J(-1)',lambda r:r['S']*r['Jm1']),('S*J(2)',lambda r:r['S']*r['J2']),('S*J(z0)',lambda r:r['S']*r['Jz']),('S*J(x0)',lambda r:r['S']*r['Jx']),('S*J(y0)',lambda r:r['S']*r['Jy']),
               ('S*sign: -1 if d<0 and beta<0',lambda r:r['S']*(-1 if r['d']<0 and r['beta']<0 else 1))]:
      print(' pairing =',nm,':',cntm(lambda r:r['P']==f(r)),'/',n, '| J factor = -1 on',cntm(lambda r:f(r)!=r['S']))
  nt=[r for r in rows if r['P']==-1]
  print(' restricted to pairing=-1 (%d pairs): S=-1 on'%len(nt),sum(1 for r in nt if r['S']==-1))
  pt=[r for r in rows if r['P']==1]
  print(' restricted to pairing=+1 (%d pairs): S=+1 on'%len(pt),sum(1 for r in pt if r['S']==1))
  print(' value table (P,S):',dict(collections.Counter((r['P'],r['S']) for r in rows)))
  print(' by parity of m+n:',dict(collections.Counter((r['par'],r['P'],r['S']) for r in rows)))
json.dump(dict(dims=dict(res['dims']),normdim={str(k):v for k,v in res['normdim'].items()},t2={str(k):v for k,v in res['t2'].items()},t3={str(k):v for k,v in res['t3'].items()},rows=ALL,bad=[str(x) for x in bad],sample=sample),open('/tmp/mix/task23_seed%d.json'%SEED,'w'))
print('partial factorizations',STATS,'time %.0f'%(time.time()-T0))
