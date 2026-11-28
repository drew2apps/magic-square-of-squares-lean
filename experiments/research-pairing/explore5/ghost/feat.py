import json,collections
from sympy import factorint
rows=[json.loads(l) for f in ('out_997.jsonl','out_15000_0.jsonl','out_15000_1.jsonl') for l in open('/tmp/rn/'+f)]
hit=set((r['m'],r['n']) for r in rows)
allr=[]
for f in ('g2_997_all.txt','g2_15000_all.txt'):
    for l in open('/tmp/rn/'+f):
        if l.strip() and l[0].isdigit():
            m,n=map(int,l.split()[0].split(',')); allr.append((m,n))
print(len(allr),len(hit))
def feats(m,n):
    C=m*m+n*n; fc=factorint(C)
    return dict(par=(m%2,n%2), m16=(m%16,n%16), wC=len(fc), Csq=tuple(sorted(e for e in fc.values())), w=len(factorint(m*n*(m*m-n*n))))
import random
for name,key in (('par',lambda f:f['par']),('wC',lambda f:f['wC']),('Cexp',lambda f:f['Csq']),('w',lambda f:f['w'])):
    c1=collections.Counter(); c2=collections.Counter()
    for (m,n) in allr:
        f=feats(m,n); (c1 if (m,n) in hit else c2)[key(f)]+=1
    print(name,'HIT',sorted(c1.items())); print(name,'stage1 only',sorted(c2.items()))
# m16 classes
c1=collections.Counter(); c2=collections.Counter()
for (m,n) in allr:
    (c1 if (m,n) in hit else c2)[(m%8,n%8)]+=1
print('mod8 hit',sorted(c1.items())); print('mod8 stage1',sorted(c2.items()))
