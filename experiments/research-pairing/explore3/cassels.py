"""General Cassels pairing on Sel^2 of y^2=(x-e1)(x-e2)(x-e3) (full rational 2-torsion).
D_a: a_i u_i^2 = x - e_i.  H_i: a_j u_j^2 - a_k u_k^2 = (e_k-e_j) u_0^2.  L_i = tangent form at a qfsolve point.
<a,b> = prod_v prod_i (L_i(P_v), b_i)_v.   Exact: rationals + PARI p-adics; real place with 150-digit mpmath."""
import cypari2, random, sys, signal, math
from fractions import Fraction as Fr
from math import gcd
import mpmath
mpmath.mp.dps=150
pari=cypari2.Pari(); pari.allocatemem(10**9)
sys.path.insert(0,'/tmp/t10')

class Timeout(Exception): pass
def _alarm(*a): raise Timeout()
signal.signal(signal.SIGALRM,_alarm)
_fcache={}
STATS={'partial_factor':0}
def primes_of(n, tmo=20):
    n=abs(int(n))
    if n<2: return set()
    if n in _fcache: return _fcache[n]
    try:
        signal.alarm(tmo)
        f=pari.factor(n); signal.alarm(0)
        s=set(int(x) for x in f[0])
    except Timeout:
        STATS['partial_factor']+=1
        f=pari.factor(n,10**7); s=set(int(x) for x in f[0] if int(x)<10**7)
    except Exception:
        signal.alarm(0); raise
    _fcache[n]=s; return s
def sqf(n):
    n=int(n); s=1 if n>0 else -1; r=1
    f=pari.factor(abs(n))
    for q,e in zip(f[0],f[1]):
        if int(e)%2: r*=int(q)
    return s*r
def issq_int(n,p):
    """nonzero integer n a square in Q_p ?"""
    if n==0: return False
    v=0
    while n%p==0: n//=p; v+=1
    if v%2: return False
    if p==2: return n%8==1
    return pow(n%p,(p-1)//2,p)==1
def vp(n,p):
    v=0
    while n%p==0: n//=p; v+=1
    return v
def class_vec(n,p):
    """class of nonzero integer n in Q_p^*/squares as a tuple of bits"""
    v=0
    while n%p==0: n//=p; v+=1
    if p==2:
        u=n%8
        return (v%2, 1 if u in (3,7) else 0, 1 if u in (5,7) else 0)   # u = (-1)^a 5^b mod 8
    return (v%2, 0 if pow(n%p,(p-1)//2,p)==1 else 1)
def sample_x(e,p,rng,vmax):
    """random rational x aimed at covering E(Q_p)/2: near a root, generic, or near infinity"""
    t=rng.random()
    if p==2: w=rng.randrange(1,2**12)
    else: w=rng.randrange(1,min(p*p,10**12))
    if t<0.6:
        j=rng.randint(0,vmax+5); return Fr(rng.choice(e)+w*p**j)
    if t<0.8:
        j=rng.randint(1,vmax//2+4); return Fr(w,p**(2*j))
    return Fr(rng.randrange(-10**6,10**6))+Fr(rng.choice(e))
def local_point(e,a,p,rng,tries=40000):
    """rational x with (x-e_i)/a_i a square in Q_p for i=1,2,3 (none zero)."""
    vmax=max(vp(e[i]-e[j],p) for i in range(3) for j in range(i))
    for _ in range(tries):
        x=sample_x(e,p,rng,vmax); n,d=x.numerator,x.denominator
        ok=True
        for i in range(3):
            t=(n-e[i]*d)*d*a[i]
            if t==0 or not issq_int(t,p): ok=False; break
        if ok: return x
    return None
def prec_for(p):
    return 400 if p==2 else max(10,int(500/math.log2(p))+2)
def psqrt(x,p,N):
    return pari('sqrt((%d)/(%d)+O(%d^%d))'%(x.numerator,x.denominator,p,N))
def tangent_forms(e,a):
    """for i: (j,k,(alpha,beta,gamma)) with L_i = alpha*u_j+beta*u_k+gamma*u_0, primitive integral"""
    out=[]
    for i in range(3):
        j,k=[t for t in range(3) if t!=i]
        c=(a[j],-a[k],-(e[k]-e[j]))
        G=pari.matdiagonal([c[0],c[1],c[2]])
        s=pari.qfsolve(G)
        if s.type()!='t_COL': raise ValueError('conic H_%d not soluble: %r'%(i+1,s))
        s=[x for x in s]; den=1
        for x in s: den=den*int(x.denominator())//gcd(den,int(x.denominator()))
        pt=[int(x*den) for x in s]; g=gcd(gcd(pt[0],pt[1]),pt[2]); pt=[x//g for x in pt]
        assert sum(c[t]*pt[t]**2 for t in range(3))==0
        L=[c[t]*pt[t] for t in range(3)]; g=gcd(gcd(L[0],L[1]),L[2]); L=[x//g for x in L]
        out.append((j,k,tuple(L),tuple(pt)))
    return out
def real_point(e,a):
    sg=tuple(1 if x>0 else -1 for x in a)
    if sg==(1,1,1): return Fr(e[2])+abs(e[2]-e[0])+Fr(7,3)
    if sg==(1,-1,-1): return Fr(2*e[0]+e[1],3)
    return None
def cassels(e,a,b,seed=1,detail=False,forms=None,extra_primes=()):
    """returns (value, info).  a,b triples of nonzero integers (squarefree reps).  value None on failure."""
    rng=random.Random(seed)
    assert all(x!=0 for x in a) and all(x!=0 for x in b)
    if forms is None: forms=tangent_forms(e,a)
    P={2}|set(extra_primes)
    for i in range(3):
        for j in range(i): P|=primes_of(e[i]-e[j])
        P|=primes_of(a[i]); P|=primes_of(b[i])
    for (j,k,L,pt) in forms:
        for c in L: P|=primes_of(c)
    tot=1; loc={}
    # real place
    x=real_point(e,a)
    if x is None: return None,{'fail':'real'}
    u=[mpmath.sqrt(mpmath.mpf(x.numerator-e[i]*x.denominator)/(x.denominator*a[i]))*rng.choice([1,-1]) for i in range(3)]
    t=1
    for i,(j,k,L,pt) in enumerate(forms):
        val=L[0]*u[j]+L[1]*u[k]+L[2]
        if abs(val)<mpmath.mpf(10)**(-60)*(abs(L[0]*u[j])+abs(L[1]*u[k])+abs(L[2])): return None,{'fail':'real cancel'}
        if val<0 and b[i]<0: t=-t
    loc[-1]=t; tot*=t
    for p in sorted(P):
        N=prec_for(p); done=False
        for attempt in range(30):
            x=local_point(e,a,p,rng)
            if x is None: return None,{'fail':('no local point',p)}
            u=[psqrt((x-e[i])/a[i],p,N) for i in range(3)]
            u=[ui*rng.choice([1,-1]) for ui in u]
            t=1; good=True
            for i,(j,k,L,pt) in enumerate(forms):
                if b[i]==1: continue
                val=L[0]*u[j]+L[1]*u[k]+L[2]
                if val==0 or int(val.padicprec(p))-int(pari.valuation(val,p))<4: good=False; break
                t*=int(pari.hilbert(val,b[i],p))
            if good: done=True; break
        if not done: return None,{'fail':('precision',p)}
        loc[p]=t; tot*=t
    return tot,(loc if detail else {})
def tors_image(e,i):
    """image in Sel^2 of the 2-torsion point (e_i,0)"""
    r=[0,0,0]; j,k=[t for t in range(3) if t!=i]
    r[j]=e[i]-e[j]; r[k]=e[i]-e[k]; r[i]=r[j]*r[k]
    return tuple(sqf(x) for x in r)
def triple(c): return (c[0],c[1],sqf(c[0]*c[1]))
