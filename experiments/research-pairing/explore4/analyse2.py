import json, itertools, collections
from sympy import factorint
def sqf(x):
    r=1
    for q,e in factorint(abs(x)).items():
        if e%2: r*=q
    return r
rows=[json.loads(l) for f in ('out_997.jsonl','out_15000_0.jsonl','out_15000_1.jsonl') for l in open(f)]
shape=collections.Counter(); aff=collections.Counter(); sizes=collections.Counter()
special=[]; nsurv_by_shape=collections.Counter()
for r in rows:
    m,n=r['m'],r['n']; T=[tuple(t) for t,_ in r['tuples']]
    for t,k in r['tuples']:
        a=sqf(m*n*t[0]*t[1]); b=sqf((m*m-n*n)*t[2]*t[3])
        s='other'
        if (a,b)==(sqf(2*m*n),1): s='S1: d0d8~2, dLdU~m2-n2'
        elif (a,b)==(1,sqf(2*(m*m-n*n))): s='S2: d0d8~mn, dLdU~2'
        elif a==1 and b==1: s='a=b=1'
        elif a==1: s='a=1 only'
        elif b==1: s='b=1 only'
        killed=any(k.values())
        shape[(s,killed)]+=1
        if s[0]=='S': special.append((m,n,t,k))
    # affine structure: vectors over F2
    def vec(t):
        v=set()
        for i,d in enumerate(t):
            for q,e in factorint(d).items():
                if e%2: v.add((i,q))
        return frozenset(v)
    V=set(vec(t) for t in T); L=list(V); ok=True
    if len(L)<=64:
        v0=L[0]
        for x in L:
            for y in L:
                if (x^y^v0) not in V: ok=False;break
            if not ok:break
        aff[ok]+=1
    else: aff['big']+=1
    sizes[len(L)]+=1
print('ratios',len(rows))
for k,v in sorted(shape.items(),key=str): print('  shape,killed',k,v)
print('affine (coset of an F2-subspace)?',dict(aff))
print('sizes of residual sets',sorted(sizes.items()))
print('special-shape tuples:',len(special),' ratios:',len(set((m,n) for m,n,_,_ in special)))
c=collections.Counter()
for m,n,t,k in special: c[tuple(sorted(x for x in k if k[x]))]+=1
for kk,v in c.most_common(): print('   killed on',kk,v)
