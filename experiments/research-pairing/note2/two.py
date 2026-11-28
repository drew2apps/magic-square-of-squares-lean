import collections, random
from fractions import Fraction as Fr
from parse import load, sqf
from ct08 import issq, hil, pari, local_point, psqrt
def cls2(x):
    x=Fr(x); n,d=x.numerator,x.denominator; v=0
    while n%2==0: n//=2; v+=1
    while d%2==0: d//=2; v-=1
    u=(n*pow(d,-1,8))%8
    return (v%2,u)
NAME={(0,1):'1',(0,7):'-1',(0,5):'5',(0,3):'-5',(1,1):'2',(1,7):'-2',(1,5):'10',(1,3):'-10'}
def t2_numeric(A,C,dl,b_r,b_m,b_s,rng):
    u=local_point(A,C,dl,2,rng)
    if u is None: return None
    if u-A==0 or u-C==0: u=u+Fr(2)**40
    t=hil(pari(f'({(u-A).numerator})/({(u-A).denominator})'),b_s,2)*hil(pari(f'({(u-C).numerator})/({(u-C).denominator})'),b_m,2)
    um=psqrt((u*u-A*A)/Fr(dl),2,80); us=psqrt((u*u-C*C)/Fr(dl),2,80); d=um-us
    if pari.valuation(d,2)>60: d=um+us
    return t*hil(d,b_r,2)
def t2_formula(A,C,dl,b_r,b_m,b_s):
    if issq(Fr(dl),2): return hil(2*psqrt(Fr(dl),2,80),b_r,2)
    if issq(Fr(-dl),2): return hil(-A,b_s,2)*hil(-C,b_m,2)*hil((A-C)*psqrt(Fr(-dl),2,80),b_r,2)
    return None
if __name__=='__main__':
    R=load(); idx={'0':0,'8':1,'L':2,'U':3}; rng=random.Random(2); st=collections.Counter(); brc=collections.Counter()
    for c in R:
        which=c['cv']
        if which not in ('08','LU'): continue
        m,n=c['ratio']; S=m*m+n*n
        for l in c['lines']:
            D=l['tup'][idx[which[0]]]*l['tup'][idx[which[1]]]
            if which=='08': A,C,dl=m*m-n*n,S,sqf(m*n*D)
            else: A,C,dl=2*m*n,S,sqf((m*m-n*n)*D)
            while A%2==0 and C%2==0: A//=2; C//=2
            shape='E (A even)' if A%2==0 else 'O (A odd)'
            b_r,b_m=l['b']; b_s=sqf(b_r*b_m)
            num=t2_numeric(A,C,dl,b_r,b_m,b_s,rng); f=t2_formula(A,C,dl,b_r,b_m,b_s)
            st[(which,shape,'delta ~ '+NAME[cls2(dl)],'formula '+('n/a' if f is None else 'agrees' if f==num else 'DISAGREES'))]+=1
            if shape.startswith('O'): brc[('b_r ~ '+NAME[cls2(b_r)],'b_m/b_* ~ '+NAME[cls2(Fr(b_m,b_s))])]+=1
    for k,v in sorted(st.items()): print(v,k)
    print('shape O: classes of b at 2:',dict(brc))
