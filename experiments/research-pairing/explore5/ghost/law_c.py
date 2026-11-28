"""Task 2(c): S1-shape tuples at the 333 ratios: E_08 class (1,sqf(2mn)) killed or not, against quartic-residue predicates."""
import json,collections
from sympy import factorint
d=json.load(open('ghostdata.json'))
def is4th(a,p):
    a%=p
    if a==0: return None
    if pow(a,(p-1)//2,p)!=1: return False
    if p%4==3: return True
    return pow(a,(p-1)//4,p)==1
rows=[]
for (m,n),tl in d['S1']:
    k08=set(k['08'] for t,k in tl); kLU=set(k['LU'] for t,k in tl)
    assert len(k08)==1
    A,C=m*m-n*n,m*m+n*n
    fa=factorint(A); fb=factorint(m*n); fc=factorint(C)
    q1=all(is4th(C,p) for p in fa if p>2)
    q1o=all(is4th(C,p) for p,e in fa.items() if p>2 and e%2)
    q2=all(is4th(C,p) for p in fb if p>2)
    # primes of C: is A*B (or 2mn, m^2-n^2) a QR / 4th power mod q
    q3=all(is4th(2*m*n,q) for q,e in fc.items() if e%2)
    q4=all(is4th(A,q) for q,e in fc.items() if e%2)
    q3b=all(pow(2*m*n%q,(q-1)//2,q)==1 for q,e in fc.items() if e%2)
    q4b=all(pow(A%q,(q-1)//2,q)==1 for q,e in fc.items() if e%2)
    rows.append(dict(m=m,n=n,k08=k08.pop(),kLU=sorted(kLU),q1=q1,q1o=q1o,q2=q2,q3=q3,q4=q4,q3b=q3b,q4b=q4b,par=(m%2,n%2),nC=len([1 for e in fc.values() if e%2]),
                     kmix=sorted(set(tuple(k[c] for c in ('0L','8U','0U','8L')) for t,k in tl))))
print('S1 ratios',len(rows),'parities',collections.Counter(r['par'] for r in rows),'killed on E_08:',sum(r['k08'] for r in rows),'killed on E_LU ever:',sum(1 in r['kLU'] for r in rows))
for name in ('q1','q1o','q2','q3','q4','q3b','q4b'):
    c=collections.Counter((r['k08'],r[name]) for r in rows); print(name,sorted(c.items()))
c=collections.Counter((r['k08'],r['q1'],r['q2'],r['q3b']) for r in rows); print('k08,q1,q2,q3b',sorted(c.items()))
c=collections.Counter((r['k08'],r['q1'] and r['q2']) for r in rows); print('k08 vs q1&q2',sorted(c.items()))
c=collections.Counter((r['k08'],r['q1'] and r['q2'] and r['q3']) for r in rows); print('k08 vs q1&q2&q3',sorted(c.items()))
c=collections.Counter((r['k08'],r['nC']) for r in rows); print('k08 vs number of odd-exponent primes of C',sorted(c.items()))
json.dump(rows,open('law_c_rows.json','w'))
