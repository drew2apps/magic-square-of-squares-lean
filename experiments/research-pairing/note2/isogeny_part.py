"""The part of Sel^2 with first component 1, and the pairing on it, from factorisations alone.
E: y^2 = x(x-A^2)(x-C^2), (A,B,C) primitive Pythagorean, C odd.
S' = { beta squarefree, supported on -1 and the primes of B : beta is +- a square in Q_p for every prime p | AC, and in Q_2 if A is even }.
Theta(delta,beta) = prod_{p | B, C = -A mod p (mod 8 for p = 2)} (delta,beta)_p * prod_{v | AC or v = 2 | A : delta not a square, -delta a square in Q_v} (AC,beta)_v."""
import sys, json, collections, itertools
from sympy import factorint
from closed_formula import issq, hil, sqf, triple
def Sprime(A,C):
    B2=C*C-A*A; gens=[-1]+sorted(factorint(B2)); places=sorted(set(factorint(abs(A*C)))|({2} if A%2==0 else set()))
    out=set()
    for mask in itertools.product((0,1),repeat=len(gens)):
        b=1
        for g,e in zip(gens,mask):
            if e: b*=g
        if all(issq(b,p) or issq(-b,p) for p in places): out.add(b)
    return out
def Theta(A,C,dl,be):
    B2=C*C-A*A; val=1
    for p in factorint(B2):
        if p==2:
            if A%2==1 and (C+A)%8==0: val*=hil(dl,be,2)
        elif (C+A)%p==0: val*=hil(dl,be,p)
    for p in set(factorint(abs(A*C)))|({2} if A%2==0 else set()):
        if (not issq(dl,p)) and issq(-dl,p): val*=hil(A*C,be,p)
    return val
if __name__=='__main__':
    res=collections.Counter(); dims=collections.Counter()
    for fn,cv in (('full0.jsonl','08'),('full1.jsonl','LU')):
        for line in open(fn):
            d=json.loads(line); A,C=triple(d['m'],d['n'],cv); SP=Sprime(A,C); seen=False
            for row in d['rows']:
                ca=row.get('ct_all') or {}
                if 'b' not in ca or row['a'][0]!=1: continue
                if not seen:
                    seen=True
                    data={sqf(b[1]) for b in ca['b'] if sqf(b[0])==1}
                    res[(cv,'S\' from factorisations = {b in Sel2 : b_r = 1}',data==SP)]+=1
                    dims[(cv,len(SP).bit_length()-1)]+=1
                dl=sqf(row['a'][1])
                res[(cv,'a in S\'',dl in SP)]+=1
                for b,v in zip(ca['b'],ca['v']):
                    if sqf(b[0])==1: res[(cv,'Theta = pairing',Theta(A,C,dl,sqf(b[1]))==v)]+=1
    for k,v in sorted(res.items()): print(v,k)
    print('dim S\':',dict(sorted(dims.items())))
