import json,collections
from sympy import factorint
def sqf(x):
    r=1
    for q,e in factorint(abs(x)).items():
        if e%2: r*=q
    return r
rows=[json.loads(l) for f in ('out_997.jsonl','out_15000_0.jsonl','out_15000_1.jsonl') for l in open('/tmp/rn/'+f)]
S1={};S2={}
for r in rows:
    m,n=r['m'],r['n']
    for t,k in r['tuples']:
        a=sqf(m*n*t[0]*t[1]); b=sqf((m*m-n*n)*t[2]*t[3])
        if (a,b)==(sqf(2*m*n),1): S1.setdefault((m,n),[]).append((tuple(t),k))
        if (a,b)==(1,sqf(2*(m*m-n*n))): S2.setdefault((m,n),[]).append((tuple(t),k))
allr=[]
for f in ('g2_997_all.txt','g2_15000_all.txt'):
    for l in open('/tmp/rn/'+f):
        if l.strip() and l[0].isdigit():
            m,n=map(int,l.split()[0].split(',')); allr.append((m,n))
print('S1 ratios',len(S1),'S2 ratios',len(S2),'both',len(set(S1)&set(S2)))
json.dump({'S1':[[k,v] for k,v in S1.items()],'S2':[[k,v] for k,v in S2.items()],'all':allr},open('ghostdata.json','w'))
# features
def pm8(x,odd_only=True):
    return collections.Counter(q%8 for q,e in factorint(x).items() if q>2 and (e%2 or not odd_only))
c=collections.Counter()
for (m,n) in allr:
    h=(m,n) in S1
    f=pm8((m-n)*(m+n)) ; ok=all(k in (1,7) for k in f)
    f2=pm8((m-n)*(m+n),False); ok2=all(k in (1,7) for k in f2)
    c[('S1' if h else 'no',ok,ok2)]+=1
print('S1 vs [odd-exp primes of m^2-n^2 all +-1 mod 8], [all primes]:',sorted(c.items()))
c=collections.Counter()
for (m,n) in allr:
    h=(m,n) in S2
    f=pm8(m*n) ; ok=all(k in (1,7) for k in f)
    f2=pm8(m*n,False); ok2=all(k in (1,7) for k in f2)
    c[('S2' if h else 'no',ok,ok2)]+=1
print('S2 vs [odd-exp primes of mn all +-1 mod 8], [all primes]:',sorted(c.items()))
