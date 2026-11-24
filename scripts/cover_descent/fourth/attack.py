import sys, os, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare import *
from mytuple import *
from math import gcd
random.seed(2026)
PR = [2,3,5,7,11,13,17,29,37,41]
st = dict(cases=0, L_in_R=0, R_in_U=0, U_in_R=0, all3=0)
bugs=[]; loose=[]; tight=[]
def one(name, forms, primes=PR, extra=()):
    for p in primes:
        U, L, N, how = local_best(forms, p, extra=extra); Ls=set(L); R = repo_allowed(forms, p)
        st['cases']+=1; a=Ls<=R; b=R<=U; c=U<=R
        st['L_in_R']+=a; st['R_in_U']+=b; st['U_in_R']+=c; st['all3']+=(Ls==R==U)
        if not a: bugs.append((name,forms,p,sorted(Ls-R)[0],L[sorted(Ls-R)[0]])); print("!!! BUG",name,forms,p,sorted(Ls-R)[0],L[sorted(Ls-R)[0]],flush=True)
        if not b: loose.append((name,forms,p,len(Ls),len(R),len(U),N))
        if not c: tight.append((name,forms,p,len(Ls),len(R),len(U),N,how))
def rnd(c=5):
    while True:
        a=(random.randint(-c,c),random.randint(-c,c)); b=(random.randint(-c,c),random.randint(-c,c))
        if a!=(0,0) or b!=(0,0): return (a,b)
t0=time.time()
# 1. random small systems
for i in range(120): one('rand%d'%i, [rnd() for _ in range(4)])
print("random:", st, len(bugs), "%.0fs"%(time.time()-t0), flush=True)
# 2. designed: two forms vanishing at the same point mod p / exactly; rational zeros; p | coefficients; (1+i) powers
for p in (2,3,5,7,13,17):
    for i in range(12):
        u0=random.randint(0,p*p); 
        f1=((1,0),(-u0,0))                        # u - u0 v : exact rational zero at (u0,1)
        f2=((1,1),(-u0+p*random.randint(0,3),-u0))    # zero mod p at the same point
        f3=((p*random.randint(1,3),random.randint(-3,3)),(random.randint(-3,3),p*random.randint(0,3)))   # p | coefficients
        k=random.randint(0,5); z=(1,0)
        for _ in range(k): z=(z[0]-z[1],z[0]+z[1])                                 # (1+i)^k
        f4=((z[0]*random.randint(1,3),z[1]),(z[1]-z[0]*random.randint(0,2),z[0]+z[1]))
        fs=[f1,f2,f3,f4]; random.shuffle(fs)
        if any(f==((0,0),(0,0)) for f in fs): continue
        one('designed p=%d #%d'%(p,i), fs, primes=sorted({2,p,5,13}))
print("designed:", st, len(bugs), "%.0fs"%(time.time()-t0), flush=True)
# 3. positive controls: forms through a chosen point (u0,v0) with values c_j*z_j^2
pc=0; pcfail=0
for i in range(60):
    while True:
        u0=random.randint(-9,9); v0=random.choice([1,1,1,2,3,5])
        if gcd(u0,v0)==1: break
    forms=[]; cls=[]
    for j in range(4):
        c=random.choice([1,2,3,5,6,7,10,13,15,17,21,34,65]); z=(random.randint(-6,6),random.randint(-6,6))
        if z==(0,0): z=(1,2)
        val=(c*(z[0]*z[0]-z[1]*z[1]), c*2*z[0]*z[1])
        # a*u0 + b*v0 = val : choose a random with a*u0 = val mod v0
        while True:
            a=(random.randint(-8,8),random.randint(-8,8))
            if (val[0]-a[0]*u0)%v0==0 and (val[1]-a[1]*u0)%v0==0: break
        b=((val[0]-a[0]*u0)//v0,(val[1]-a[1]*u0)//v0)
        forms.append((a,b)); cls.append(c)
    for p in [q for q in primes_upto(100)]:
        R = repo_allowed(forms,p); tp=tuple(hclass(c,p) for c in cls)
        et = exact_tuple(forms,u0,v0,p); assert et==tp, ("my exact class disagrees with construction",forms,p,et,tp)
        pc+=1
        if tp not in R: pcfail+=1; print("!!! planted point rejected by repo", forms,(u0,v0),cls,p,flush=True)
    one('planted%d'%i, forms, primes=[2,3,5,13,17], extra=[(u0,v0)])
print("planted points: %d (system,prime) checks at primes <100, rejected by repo: %d"%(pc,pcfail))
print("ALL:", st); print("BUGS:", bugs)
print("R not in U:", len(loose), loose[:6]); print("U not in R:", len(tight)); 
for x in tight[:12]: print("   ", x)
