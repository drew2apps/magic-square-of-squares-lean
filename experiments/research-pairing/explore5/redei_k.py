# Test of the identity  <(1,d,d),(1,beta,beta)> = prod_w (gamma_w, pi)_{K_w},  K = Q(sqrt d),
# gamma_w = (u1 + sqrt(d) u3)(u1 + sqrt(d) u2) at a local point of the 2-covering, N(pi) = beta.
import sys, random, collections, itertools
sys.argv=['x']; sys.path.insert(0,'/tmp/rn'); sys.path.insert(0,'/tmp/mix')
from pipeline import Curve, cassels, triple as tr3, sqf, P
from cassels import local_point, primes_of, vp
from fractions import Fraction as Fr
from math import gcd
gp=P
def padic_sqrt_rat(t,p,N):
    """rational r with r^2 = t up to relative precision p^N (t a rational that is a p-adic square)"""
    s=gp('(t,p,N)->my(r=sqrt(t+O(p^(valuation(t,p)+N))));lift(r*p^(-valuation(r,p)))*p^valuation(r,p)')(gp(f'{t.numerator}/{t.denominator}'),p,N)
    return Fr(int(gp.numerator(s)),int(gp.denominator(s)))
def rhs(e,d,beta,rng):
    sol=gp.qfsolve(gp.matdiagonal([1,-d,-beta]))
    if str(gp.type(sol))!='t_COL': return None
    x0,y0,z0=[int(c) for c in sol]
    if z0==0: return None
    nf=gp.nfinit(gp(f'y^2-({d})'))
    def el(a,b): return gp(f'Mod(({a.numerator}/{a.denominator})+({b.numerator}/{b.denominator})*y,y^2-({d}))')
    pi=el(Fr(x0,z0),Fr(y0,z0))
    S=set()
    for n in (2,d,beta,e[1]-e[0],e[2]-e[0],e[2]-e[1],z0,x0*x0-d*y0*y0):
        S|=set(primes_of(abs(n))) if abs(n)>1 else set()
    S.add(2); val=1; a=(1,d,d)
    for p in sorted(S):
        x=local_point(e,a,p,rng)
        if x is None: return None
        N=60 if p==2 else 24
        u=[padic_sqrt_rat((x-e[i])/a[i],p,N) for i in range(3)]
        for pr in gp.idealprimedec(nf,p):
            # gamma as an nf element: (u1 + y u3)(u1 + y u2)
            g=el(u[0],u[2])*el(u[0],u[1])
            val*=int(gp.nfhilbert(nf,g,pi,pr))
    # real places
    if d>0:
        xr=max(e)+1+rng.random()*abs(e[2]-e[0]); import math
        u=[math.sqrt((xr-e[i])/a[i]) for i in range(3)]; sd=math.sqrt(d)
        for s in (1,-1):
            g=(u[0]+s*sd*u[2])*(u[0]+s*sd*u[1]); pv=(x0+s*sd*y0)/z0
            if g<0 and pv<0: val=-val
    return val
def run(e,rng,cnt,label):
    E=Curve(e); els=[(1,1)]
    def mul(x,y): return tuple(sqf(x[i]*y[i]) for i in (0,1))
    for b in E.basis: els+=[mul(x,b) for x in els]
    one=[c for c in els if c[0]==1 and c[1]!=1]
    for a in one:
        for b in one:
            r=rhs(e,a[1],b[1],rng)
            if r is None: cnt[(label,'no pi or no point')]+=1; continue
            ref=cassels(e,tr3(a),tr3(b),seed=7)[0]
            cnt[(label,r==ref)]+=1
rng=random.Random(11); cnt=collections.Counter()
ratios=[(465,416),(511,260),(627,364),(3320,2303),(9113,3200),(23,10),(41,12),(88,35),(121,56),(200,77)]
for m,n in ratios:
    A,B,C=m*m-n*n,2*m*n,m*m+n*n
    run((0,2*(C*C-A*B),2*(C*C+A*B)),rng,cnt,'E_mix'); run((0,A*A,C*C),rng,cnt,'E_08'); run((0,B*B,C*C),rng,cnt,'E_LU')
for _ in range(60):
    e=sorted(rng.sample(range(-300,300),3)); run(tuple(e),rng,cnt,'random e')
for k,v in sorted(cnt.items(),key=str): print(k,v)
