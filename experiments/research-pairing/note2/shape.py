import json, collections, sys
from parse import sqf
from sympy import factorint
def odd_sqf(n):
    r=1
    for p,e in factorint(n).items():
        if p!=2 and e%2: r*=p
    return r
st=collections.Counter(); detail=collections.Counter(); dims=collections.Counter(); ex=[]
for fn,cv in (('full0.jsonl','08'),('full1.jsonl','LU')):
    for line in open(fn):
        d=json.loads(line); m,n=d['m'],d['n']; S=m*m+n*n; c=odd_sqf(S)
        A=(m*m-n*n) if cv=='08' else 2*m*n
        for row in d['rows']:
            if not any(w[0]=='tuple' for w in row['why']): continue
            ca=row['ct_all']; bs=ca['b']; vs=ca['v']
            killers=[b for b,v in zip(bs,vs) if v==-1]
            if not killers: st[(cv,'not killed')]+=1; continue
            k1=[b for b in killers if abs(sqf(b[0]))==1]
            kS=[b for b in killers if c % abs(odd_sqf(abs(b[0])))==0 and abs(sqf(b[0])) not in (1,)]
            kind='killed by some b with b_r = ±1' if k1 else 'only by b with b_r != ±1'
            st[(cv,kind)]+=1
            if not k1:
                # what divides b_r in the killers: primes of c (odd squarefree part of S), of A, 2?
                tags=set()
                for b in killers:
                    br=abs(sqf(b[0])); t=[]
                    for p in factorint(br):
                        t.append('2' if p==2 else 'c' if c%p==0 else 'S' if S%p==0 else 'A' if A%p==0 else 'B')
                    tags.add(''.join(sorted(set(t))))
                detail[(cv,tuple(sorted(tags)))]+=1
                if len(ex)<4: ex.append((m,n,cv,row['a'],d['sel'],killers[:3]))
for k,v in sorted(st.items()): print(v,k)
print('when only b_r != ±1 kills: prime types in b_r of the killers (c = odd squarefree part of S):')
for k,v in sorted(detail.items(),key=lambda x:-x[1])[:12]: print('  ',v,k)
for e in ex: print(e)
