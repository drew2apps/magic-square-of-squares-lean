#!/usr/bin/env python3
"""The test of class_complete.py on random ratios.  A test, not a proof.  Python only.

Random coprime pairs m > n >= 1 with m up to a bound (fixed seed).  A pair with no point at all modulo some bad prime
(the four norms are squares at no point of the projective line modulo p^k: the ordinary local obstruction, no class
involved) is left out, and so is a pair with more than 13 primes in 2mn(m^2 - n^2).  For the others, every squarefree
product d of the primes of 2mn(m^2 - n^2) is tested as in class_complete.py: the five laws against the local tests at
the bad primes below the bound.
Usage: class_complete_random.py <bound on the primes> <number of pairs> <bound on m>
  200 1500 3000   (about 25 s; output class_complete_random.txt, rerun by CI)
  400 6000 20000  (about 2 minutes; output class_complete_random_20000.txt, not rerun by CI)
"""
import sys, random, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import lean_cert as L
from math import gcd
from collections import Counter
B0=int(sys.argv[1]); NP=int(sys.argv[2]); MAXM=int(sys.argv[3]); random.seed(20261124)
SP=[q for q in range(3,B0,2) if all(q%r for r in range(3,int(q**.5)+1,2))]
def K(p): return {2:6,3:4,5:3,7:2,11:2,13:2,17:2,19:2}.get(p,1 if p>40 else 2)
pairs=set()
while len(pairs)<NP:
    m=random.randint(2,MAXM); n=random.randint(1,m-1)
    if gcd(m,n)==1: pairs.add((m,n))
new=[]; obstructed=0; weak=Counter(); nl=nb=ncand=0
for (m,n) in sorted(pairs):
    P=sorted(L.factor(2*m*n*(m*m-n*n)))
    if len(P)>13: continue
    C=L.prods(P); ncand+=len(C)
    N1=(m*m+n*n)*(m**8+14*m**4*n**4+n**8); N=2*m*n*(m*m-n*n)*N1
    def law(d):
        if d%2: return 'F'
        for q in P:
            if q>2 and d%q==0 and q%8 not in (1,7): return 'E'
        if L.exact5(m,n) and d%5 in (1,4): return 'L5'
        for q in SP:
            if N1%q==0 and not L.eul(d,q): return 'GO'
        return ''
    A=set(d for d in C if not law(d)); Bs=set(C)
    bad=[p for p in [2]+SP if N%p==0]; SV={p:L.survivors(m,n,p**K(p),p) for p in bad}
    if any(not SV[p] for p in bad): obstructed+=1; ncand-=len(C); continue
    for p in bad:
        M=p**K(p); surv=SV[p]
        for d in list(Bs):
            if L.killed(d,surv,M):
                Bs.discard(d)
                if d in A: new.append(((m,n),p,M,d))
    for d in Bs-A: weak[law(d)]+=1
    nl+=len(C)-len(A); nb+=len(C)-len(Bs)
print('random coprime pairs %d with m <= %d, candidates %d, bound %d'%(len(pairs),MAXM,ncand,B0))
print('pairs with no point at all modulo some bad prime (left out): %d'%obstructed)
print('excluded by the laws: %d; by the local tests at the bad primes: %d'%(nl,nb))
print('excluded locally at a bad prime but allowed by the laws: %d' % len(new))
if new: print('  by prime: %s; first: %s' % (sorted(Counter(e[1] for e in new).items()), new[:6]))
print('excluded by a law but not by the local test at the precision used:',dict(weak))
