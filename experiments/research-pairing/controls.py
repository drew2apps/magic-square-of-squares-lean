import collections, random
from parse import load, sqf
from ct08 import pairing
R=load(); idx={'0':0,'8':1,'L':2,'U':3}; res=collections.Counter(); rng=random.Random(5)
L=[(c,l) for c in R if c['cv'] in ('08','LU') for l in c['lines']]
for c,l in rng.sample(L,150):
    which=c['cv']; m,n=c['ratio']; S=m*m+n*n
    D=l['tup'][idx[which[0]]]*l['tup'][idx[which[1]]]
    if which=='08': A,B,C,dl=m*m-n*n,2*m*n,S,sqf(m*n*D)
    else: A,B,C,dl=2*m*n,m*m-n*n,S,sqf((m*m-n*n)*D)
    for name,(br,bm) in {'a itself':(1,dl),'2-torsion (1,-1,-1)':(1,-1),'point of order 4':(sqf(C*(C+B)),sqf(B*(B+C))),'other point of order 4':(sqf(C*(C-B)),sqf(-B*(C-B)))}.items():
        v,_=pairing(A,C,dl,br,sqf(br*bm))
        res[(name,v)]+=1
    # bilinearity in b: <a,b>*<a,b'> = <a,bb'> with b the killing class and b' the point of order 4
    br,bm=l['b']; pr,pm=sqf(C*(C+B)),sqf(B*(B+C))
    v1,_=pairing(A,C,dl,br,sqf(br*bm)); v2,_=pairing(A,C,dl,sqf(br*pr),sqf(br*pr*bm*pm))
    res[('killing b times the point of order 4',v2)]+=1
for k,v in sorted(res.items()): print(k,v)
