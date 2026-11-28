"""Task 1 sampler beyond m = 15000.  usage: sampler.py MODE SEED SECONDS OUT.jsonl
MODE biased: m+n, m-n built from primes = +-1 mod 8 (necessary for the S1 ghost shape on the 499 ratios), half of them
             replaced by the partner (m+n, m-n) (S2 candidates).
MODE random: random coprime m > n, 15000 < m <= 60000, with the mod-8 pattern seen on all 9,736 stage-1 ratios.
MODE fewp  : m^2+n^2 = p, 2p, pq, 2pq or k^2 p (few prime factors), random otherwise (same mod-8 pattern)."""
import sys,time,json,random
mode,seed,secs,outf=sys.argv[1],int(sys.argv[2]),float(sys.argv[3]),sys.argv[4]
sys.argv=['x']; sys.path.insert(0,'/tmp/rn')
from pipeline import *
from math import gcd
rng=random.Random(seed)
def rprime(lo,hi):
    while True:
        p=int(P.nextprime(rng.randrange(lo,hi)))
        if p%8 in (1,7): return p
def mod8ok(m,n):
    if (m+n)%2: return (m*n)%4==0
    return (m-n)%8==0 or (m+n)%8==0
def cand():
    if mode=='biased':
        while True:
            p=rprime(20000,150000) if rng.random()<0.7 else rprime(100,400)*rprime(100,400)
            r=rng.random()
            q=1 if r<0.1 else rprime(3,p) if r<0.7 else rprime(3,300)*rprime(3,300)
            if q>=p or gcd(p,q)>1: continue
            m,n=(p+q)//2,(p-q)//2
            if m<=15000: continue
            if rng.random()<0.5: return (m,n)
            return (p,q)
    if mode=='random':
        while True:
            m=rng.randrange(15001,60000); n=rng.randrange(1,m)
            if gcd(m,n)==1 and mod8ok(m,n): return (m,n)
    if mode=='fewp':
        while True:
            m=rng.randrange(15001,60000); n=rng.randrange(1,m)
            if gcd(m,n)!=1 or not mod8ok(m,n): continue
            C=m*m+n*n; f=fact(C)
            odd=[q for q,e in f.items() if e%2 and q>2]
            if len(odd)<=2 and not (len(odd)==0): return (m,n)
def isdeg(m,n):
    C=m*m+n*n; return int(P.issquare(C)) or int(P.issquare(2*C))
t0=time.time(); tried=0; seen=set()
with open(outf,'w') as fo:
    while time.time()-t0<secs:
        m,n=cand()
        if (m,n) in seen or isdeg(m,n): continue
        seen.add((m,n)); tried+=1
        try:
            T,known=residual(m,n)
            if not T: continue
            if len(T)>400:
                fo.write(json.dumps(dict(m=m,n=n,skipped=len(T)))+'\n'); fo.flush(); continue
            cs=curves(m,n); tl=[]
            for t in T:
                k={cv:E.killed(f(t)) for cv,(E,f) in cs.items()}
                tl.append([list(t),{cv:(v if v=='outside' else int(v)) for cv,v in k.items()}])
            fo.write(json.dumps(dict(m=m,n=n,tuples=tl,tried=tried))+'\n'); fo.flush()
        except Exception as ex:
            fo.write(json.dumps(dict(m=m,n=n,error=repr(ex)[:200]))+'\n'); fo.flush()
    fo.write(json.dumps(dict(done=1,tried=tried,secs=round(time.time()-t0)))+'\n')
