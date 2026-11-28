"""Closed forms for the local terms at odd primes, against the numerical terms.  Curve y^2=x(x-A^2)(x-C^2), class a=(1,dl,dl), b=(b_r,b_m,b_s)."""
import sys, collections
from fractions import Fraction as Fr
from parse import load, sqf
from ct08 import pairing, issq, hil, pari
def psq(x,p,N=40): return pari(f'sqrt(({x})+O({p}^{N}))')
def predict(A,C,dl,b_r,b_m,b_s,p):
    if (C-A)%p==0: return 'p|C-A',1
    if (C+A)%p==0: return 'p|C+A',hil(dl,b_s,p)
    if issq(Fr(dl),p): return 'dl square',hil(2*psq(dl,p),b_r,p)
    if issq(Fr(-dl),p):
        r=psq(-dl,p)
        return '-dl square',hil(-A,b_s,p)*hil(-C,b_m,p)*hil((A-C)*r,b_r,p)
    return 'none',None
R=load(); idx={'0':0,'8':1,'L':2,'U':3}; agree=collections.Counter(); tot=collections.Counter(); two=collections.Counter(); ex=[]
for which in sys.argv[1:]:
    for c in R:
        if c['cv']!=which: continue
        m,n=c['ratio']; S=m*m+n*n
        for l in c['lines']:
            D=l['tup'][idx[which[0]]]*l['tup'][idx[which[1]]]
            if which=='08': A,C,dl=m*m-n*n,S,sqf(m*n*D)
            else: A,C,dl=2*m*n,S,sqf((m*m-n*n)*D)
            b_r,b_m=l['b']; b_s=sqf(b_r*b_m)
            v,loc=pairing(A,C,dl,b_r,b_s)
            tot[(which,v)]+=1
            for p,t in loc.items():
                if p==2: two[(which,t)]+=1; continue
                kind,pr=predict(A,C,dl,b_r,b_m,b_s,p)
                agree[(which,kind,pr==t)]+=1
                if pr!=t and len(ex)<6: ex.append((which,m,n,l['tup'],dl,(b_r,b_m),p,kind,pr,t))
print('totals (curve, pairing value):',dict(tot))
print('odd local terms (curve, case, prediction agrees):'); 
for k,v in sorted(agree.items()): print('  ',k,v)
print('2-adic terms:',dict(two)); print('examples of disagreement:',ex)
