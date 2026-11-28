"""Residual tuples (the repository's exact local test at the bad primes) against the Cassels-Tate pairing on the six pair
curves (general code of /tmp/mix).  For each ratio: which residual tuples are killed on which curve."""
import sys, os, random, collections, json, time
os.environ['MSS_REPO']='/tmp/rn/repo'
sys.path.insert(0,'/tmp/rn/repo/experiments/roadmap-g/g1'); sys.path.insert(0,'/tmp/mix')
import exact_local as X
TD=X.TD
from cassels import cassels, triple, sqf, primes_of
from sel import local_image, nullspace, cls2, inspan
import cypari2
P=cypari2.Pari(); P.allocatemem(2*10**9)
def fact(x):
    x=abs(int(x))
    if x<=1: return {}
    f=P.factor(x); return {int(q):int(e) for q,e in zip(f[0],f[1])}
TD.factorint=fact
def residual(m,n):
    forms,sup,known=TD.magic(m,n)
    divs=[]
    for S in sup:
        ds=[1]
        for q in S: ds+=[x*q for x in ds]
        divs.append(ds)
    bad=set()
    for x in (2*m*n*(m*m-n*n), m*m+n*n, m**8+14*m**4*n**4+n**8): bad|=set(fact(x))
    loc=[]
    for p in sorted(bad):
        L,al=X.allowed_at(forms,p)
        if len(al)<len(L.labels())**len(forms): loc.append((p,L,al))
    return sorted(X.join(divs,loc)), known
class Curve:
    def __init__(s,e,seed=7):
        s.e=e; rng=random.Random(seed)
        s.S=sorted({2}|set().union(*[primes_of(abs(e[i]-e[j])) for i in range(3) for j in range(i)]))
        s.LI={p:local_image(e,p,rng) for p in s.S}
        g2=[-1]+s.S; gp=[(g,1) for g in g2]+[(1,g) for g in g2]
        def t(x): return [(s.LI[p],cls2(x[0],x[1],p)) for p in s.S]+[([],1 if x[0]<0 else 0)]
        s.basis=[]
        for k in nullspace(gp,t):
            c1=c2=1
            for i,x in enumerate(gp):
                if k>>i&1: c1*=x[0]; c2*=x[1]
            s.basis.append((sqf(c1),sqf(c2)))
        s.cache={}
    def insel(s,c): return c[0]>0 and all(inspan(s.LI[p],cls2(c[0],c[1],p)) for p in s.S)
    def killed(s,c):
        c=(sqf(c[0]),sqf(c[1]))
        if c not in s.cache:
            if not s.insel(c): s.cache[c]='outside'
            else: s.cache[c]=any(cassels(s.e,triple(c),triple(b),seed=3)[0]==-1 for b in s.basis)
        return s.cache[c]
def curves(m,n):
    A,B,C=m*m-n*n,2*m*n,m*m+n*n; Rp,Rm=C*C+A*B,C*C-A*B; k=sqf(m*n*(m*m-n*n))
    EA=Curve((0,A*A,C*C)); EB=Curve((0,B*B,C*C)); EM=Curve((0,2*Rm,2*Rp))
    return {'08':(EA,lambda t:(1,m*n*t[0]*t[1])),'LU':(EB,lambda t:(2*C,2*t[2]*t[3])),
            '0L':(EM,lambda t:(2,k*t[0]*t[2])),'8U':(EM,lambda t:(2,k*t[1]*t[3])),
            '0U':(EM,lambda t:(2,2*k*t[0]*t[3])),'8L':(EM,lambda t:(2,2*k*t[1]*t[2]))}
def run(m,n):
    T,known=residual(m,n); cs=curves(m,n); out={}
    for t in T:
        out[t]={cv:E.killed(f(t)) for cv,(E,f) in cs.items()}
    dims={cv:len(E.basis) for cv,(E,f) in cs.items()}
    return out,dims,known
if __name__=='__main__':
    lst=sys.argv[1]; i,N=int(sys.argv[2]),int(sys.argv[3]); outf=sys.argv[4]
    pairs=[tuple(int(x) for x in l.split()[0].split(',')) for l in open(lst) if 'left' in l][i::N]
    with open(outf,'w') as fo:
        for m,n in pairs:
            t0=time.time()
            try:
                out,dims,known=run(m,n)
                rec=dict(m=m,n=n,dims=dims,tuples=[[list(t),{cv:(v if v=='outside' else int(v)) for cv,v in d.items()}] for t,d in out.items()],time=round(time.time()-t0,1))
            except Exception as ex:
                rec=dict(m=m,n=n,error=repr(ex)[:200])
            fo.write(json.dumps(rec)+'\n'); fo.flush()
