import sys,time,random,collections,re,ast,itertools
sys.path.insert(0,'/tmp/mix'); sys.path.insert(0,'/tmp/t10')
from parse import load
from cassels import *
from sel import *
R=load(); T0=time.time(); rng=random.Random(99)
# (1) all 08 and LU lines
for cv in ('08','LU'):
    cnt=collections.Counter()
    for c in R:
        if c['cv']!=cv: continue
        for l in c['lines']:
            cnt[cassels(c['e'],triple(l['cls']),triple(l['b']))[0]]+=1
    print('all lines',cv,dict(cnt),'t=%.0f'%(time.time()-T0),flush=True)
# (2) seeds really give different local points
M=[(c,l) for c in R if c['cv'] in ('0L','0U','8L','8U') for l in c['lines']]
diff=0; tot=0
for c,l in rng.sample(M,50):
    a=triple(l['cls']); e=c['e']
    for p in sorted({2}|primes_of(e[1]-e[0])|primes_of(e[2]-e[0])|primes_of(e[2]-e[1])):
        x1=local_point(e,a,p,random.Random(1)); x2=local_point(e,a,p,random.Random(777)); tot+=1; diff+=(x1!=x2)
print('seed sanity: local points differ between seeds 1 and 777 at',diff,'of',tot,'(line,prime) pairs')
# (3) brute-force cross-check of the ELS set {d:(1,d,d)} by direct local point search, 3 ratios with fewest generators
curves={c['ratio']:c for c in R if c['cv']=='0L'}
sample=ast.literal_eval(open('/tmp/mix/task23_seed11.json').read().split('"sample": ')[1].rstrip('}').replace('[','(').replace(']',')').join(['[',']'])) if False else None
import json
J=json.load(open('/tmp/mix/task23_seed11.json')); sample=[tuple(x) for x in J['sample']]
def setup(m,n):
    c=curves[(m,n)]; e=c['e']; A=m*m-n*n;B=2*m*n;C=m*m+n*n; Rp=C*C+A*B; Rm=C*C-A*B
    S=sorted({2}|primes_of(A)|primes_of(B)|primes_of(C)|primes_of(Rp)|primes_of(Rm))
    lr=random.Random(1000+m); LI={p:local_image(e,p,lr) for p in S}
    gens=[-1]+sorted(primes_of(2*A*B))
    ker=nullspace(gens,lambda g:[(LI[p],cls2(1,g,p)) for p in S])
    return e,S,LI,gens,ker
small=sorted(sample,key=lambda r:len(primes_of(2*(r[0]**2-r[1]**2)*2*r[0]*r[1])))[:3]
for (m,n) in small:
    e,S,LI,gens,ker=setup(m,n); kb=rank_reduce(ker); ok=collections.Counter()
    for mask in range(1<<len(gens)):
        d=1
        for i,g in enumerate(gens):
            if mask>>i&1: d*=g
        member=inspan(kb,mask)
        found=all(local_point(e,(1,d,d),p,rng,tries=3000) is not None for p in S)
        ok[(member,found)]+=1
    print('brute force ratio',(m,n),'gens',len(gens),'(in linear-algebra set, local points found at all tested primes):',dict(ok),'t=%.0f'%(time.time()-T0),flush=True)
# (4) how much of the pairing on V x V is forced by the Redei values on soluble pairs + symmetry + <d,d>=1 ?
logl=[l for l in open('/tmp/mix/task23_seed11.log') if l.startswith('ratio (')]
full=collections.Counter(); frac=[]
for l in logl:
    mm=re.match(r'ratio \((\d+), (\d+)\).*basis (\[[-\d, ]+\])',l); m,n=int(mm[1]),int(mm[2]); basis=ast.literal_eval(mm[3]); k=len(basis)
    P=sorted({2}|set().union(*[primes_of(b) for b in basis]))
    H=[[[0 if int(pari.hilbert(basis[i],basis[j],p))==1 else 1 for j in range(k)] for i in range(k)] for p in P]
    H.append([[1 if (basis[i]<0 and basis[j]<0) else 0 for j in range(k)] for i in range(k)])
    # row masks
    HM=[[sum(h[i][j]<<j for j in range(k)) for i in range(k)] for h in H]
    def sol(x,y):
        for hm in HM:
            t=0
            for i in range(k):
                if x>>i&1: t^=hm[i]
            if bin(t&y).count('1')&1: return False
        return True
    def tens(x,y):
        r=0
        for i in range(k):
            if x>>i&1: r|=y<<(i*k)
        return r
    rows=[]; nsol=0
    for x in range(1,1<<k):
        rows.append(tens(x,x))
        for y in range(1,1<<k):
            if x<y: rows.append(tens(x,y)^tens(y,x))
            if sol(x,y): rows.append(tens(x,y)); nsol+=1
    rk=len(rank_reduce(list(set(rows))))
    # symmetric+alternating relations alone have rank k(k+1)/2 ; unknown part of the form has dim k(k-1)/2
    forced=rk-k*(k+1)//2; full[(forced==k*(k-1)//2)]+=1; frac.append((k,forced,k*(k-1)//2,nsol,(2**k-1)**2))
print('pairing determined by soluble-pair values + symmetry + alternation (ratios):',dict(full))
print(' (k, #independent entries forced, k(k-1)/2, #soluble ordered pairs, #pairs) per ratio:',frac)
print('time %.0f'%(time.time()-T0))
