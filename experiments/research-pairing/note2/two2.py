import collections, random
from fractions import Fraction as Fr
from parse import load, sqf
from ct08 import hil
from two import t2_numeric, cls2, NAME
def v2(x):
    v=0
    while x%2==0: x//=2; v+=1
    return v
R=load(); idx={'0':0,'8':1,'L':2,'U':3}; rng=random.Random(4); st=collections.Counter(); bad=collections.Counter()
for c in R:
    which=c['cv']
    if which not in ('08','LU'): continue
    m,n=c['ratio']; S=m*m+n*n
    for l in c['lines']:
        D=l['tup'][idx[which[0]]]*l['tup'][idx[which[1]]]
        if which=='08': A,C,dl=m*m-n*n,S,sqf(m*n*D)
        else: A,C,dl=2*m*n,S,sqf((m*m-n*n)*D)
        while A%2==0 and C%2==0: A//=2; C//=2
        if A%2==0: continue
        b_r,b_m=l['b']; b_s=sqf(b_r*b_m)
        num=t2_numeric(A,C,dl,b_r,b_m,b_s,rng)
        minus=(C-A)%8==0; k=v2(C-A) if minus else v2(C+A)
        pred=1 if minus else hil(dl,b_m,2)
        st[('O-' if minus else 'O+','k=%d'%min(k,7),'agree' if pred==num else 'DIFF')]+=1
        if pred!=num: bad[('O-' if minus else 'O+','k=%d'%k,'delta~'+NAME[cls2(dl)],'beta~'+NAME[cls2(b_m)])]+=1
for k,v in sorted(st.items()): print(v,k)
print('disagreements by class:'); 
for k,v in sorted(bad.items()): print('  ',v,k)
