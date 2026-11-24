#!/usr/bin/env python3
"""Are the five laws on the class everything that the bad primes say?  A test, not a proof.  Python only; about 10 s.

For each of the 425 ratios of `exceptions1000` and each squarefree product d of the primes of 2mn(m^2 - n^2) (all of them:
odd ones too, and primes not = ±1 mod 8), compare
  the laws:   F (d even), E (odd primes of d are ±1 mod 8), the law at 5, G and O at the primes below the bound;
  the points: the local test of lean_cert.py (`survivors`, `killed`) at every prime below the bound that divides
              2mn(m^2 - n^2)(m^2 + n^2)(m^8 + 14 m^4 n^4 + n^8), modulo 2^6, 3^4, 5^3, p^2 for p < 40 and p above.
The laws are theorems, so a class that they exclude has no point; the test asks for the converse.
Usage: class_complete.py [bound, default 100]
"""
import sys, re, time, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import lean_cert as L
from collections import Counter
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..')
src = open(os.path.join(ROOT, 'MagicSquareOfSquares/G1Census1000.lean')).read()
body = src[src.index('def exceptions1000'):]; body = body[:body.index(']')]
pairs = [tuple(int(x) for x in t) for t in re.findall(r'\((\d+), (\d+)\)', body)]
B0=int(sys.argv[1]) if len(sys.argv)>1 else 100
SP=[q for q in range(3,B0,2) if all(q%r for r in range(3,int(q**.5)+1,2))]
def K(p): return {2:6,3:4,5:3,7:2,11:2,13:2,17:2,19:2}.get(p,1 if p>40 else 2)
t0=time.time(); new=[]; weak=Counter(); nl=nb=0; ncand=0
for (m,n) in pairs:
    P=sorted(L.factor(2*m*n*(m*m-n*n))); C=L.prods(P); ncand+=len(C)
    N1=(m*m+n*n)*(m**8+14*m**4*n**4+n**8); N=2*m*n*(m*m-n*n)*N1
    def law(d):
        if d%2: return 'F'
        for q in P:
            if q>2 and d%q==0 and q%8 not in (1,7): return 'E'
        if L.exact5(m,n) and d%5 in (1,4): return 'L5'
        for q in SP:
            if N1%q==0 and not L.eul(d,q): return 'GO'
        return ''
    A=set(d for d in C if not law(d))
    bad=[p for p in [2]+SP if N%p==0]; Bs=set(C)
    for p in bad:
        M=p**K(p); surv=L.survivors(m,n,M,p)
        for d in list(Bs):
            if L.killed(d,surv,M):
                Bs.discard(d)
                if d in A: new.append(((m,n),p,M,d,'supp' if p in P else 'other', (m*n*(m-n)*(m+n))%(p*p)==0, d%p==0))
    for d in A-Bs: pass
    for d in Bs-A: weak[law(d)]+=1
    nl+=len(C)-len(A); nb+=len(C)-len(Bs)
print('ratios %d, candidates %d, bound %d'%(len(pairs),ncand,B0))
print('excluded by the laws (restricted to primes < %d): %d; by the local tests at the bad primes < %d: %d'%(B0,nl,B0,nb))
print('excluded locally at a bad prime but allowed by the laws:',len(new)); print('  by prime: %s' % sorted(Counter(e[1] for e in new).items()))
print('excluded by a law but not by the local test at the precision used:',dict(weak))
