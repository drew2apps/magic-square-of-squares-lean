"""Task 2(b): when is the ghost class (1,2mn,2mn) (= class of the 4-torsion point times (C,1,C)) in Sel^2(E_08)?
Random coprime ratios m > n (not only stage-1 ratios)."""
import sys,json,random,collections
sys.argv=['x']; sys.path.insert(0,'/tmp/rn')
from pipeline import *
from math import gcd
rng=random.Random(11); N=3000; rows=[]
while len(rows)<N:
    m=rng.randrange(2,4000); n=rng.randrange(1,m)
    if gcd(m,n)!=1: continue
    A,B,C=m*m-n*n,2*m*n,m*m+n*n
    if int(P.issquare(C)): continue
    E=Curve((0,A*A,C*C)); ins=E.insel((1,sqf(2*m*n)))
    fa=fact(A)
    odd5=any(q%8==5 for q in fa)                 # a prime = 5 mod 8 divides m^2-n^2
    rows.append((m,n,ins,odd5))
json.dump(rows,open('law_b.json','w'))
c=collections.Counter((r[2],r[3]) for r in rows)
print('N =',N,' (in Sel2, some p=5 mod 8 divides m^2-n^2):',sorted(c.items()))
# 2-adic part among those with no p = 5 mod 8
tab=collections.defaultdict(collections.Counter)
for m,n,ins,o in rows:
    if not o: tab[(m%16,n%16)][ins]+=1
mixed=[(k,dict(v)) for k,v in tab.items() if len(v)>1]
print('classes (m,n) mod 16 with both outcomes, given the odd condition:',len(mixed),mixed[:10])
# guess a 2-adic rule
def rule2(m,n):
    if (m+n)%2==0: return (m-n)%8==0 or (m+n)%8==0       # both odd: 16 | m^2-n^2
    return (m*n)%4==0                                     # opposite parity: 4 | mn  (then m^2+n^2 = 1 mod 8)
c=collections.Counter()
for m,n,ins,o in rows: c[(ins,(not o) and rule2(m,n))]+=1
print('test: in Sel2  <=>  [no prime = 5 mod 8 divides m^2-n^2] and [4 | mn, or m,n odd with 8 | m-n or 8 | m+n]:',sorted(c.items()))
tab2=collections.Counter()
for m,n,ins,o in rows:
    if not o: tab2[((m%8,n%8),ins)]+=1
print(sorted(tab2.items()))
