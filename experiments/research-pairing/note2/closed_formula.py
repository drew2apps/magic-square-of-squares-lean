"""The pairing on E: y^2 = x(x-A^2)(x-C^2), (A,B,C) a primitive Pythagorean triple with C odd, by the closed formula alone:
   <(1,dl,dl),(b_r,b_m,b_s)> = prod_{p | B, C = -A mod p (mod 8 if p = 2)} (dl,b_s)_p  *  prod_{v | AC, and v = 2 if A is even} t_v,
   t_v = (2 sqrt(dl), b_r)_v if dl is a square in Q_v, else (-A,b_s)_v (-C,b_m)_v ((A-C) sqrt(-dl), b_r)_v.
No local point is searched for.  Returns None if at some v | AC neither dl nor -dl is a square (then (1,dl,dl) is not locally soluble)."""
import sys, json, collections
from fractions import Fraction as Fr
from sympy import factorint
import cypari2
pari=cypari2.Pari(); pari.allocatemem(10**9)
def issq(x,p):
    x=Fr(x); n,d=x.numerator,x.denominator; v=0
    while n%p==0: n//=p; v+=1
    while d%p==0: d//=p; v-=1
    if v%2: return False
    u=n*pow(d,-1,8 if p==2 else p)%(8 if p==2 else p)
    return u%8==1 if p==2 else pow(u,(p-1)//2,p)==1
def hil(x,y,p): return int(pari.hilbert(x,y,p))
def psq(x,p): return pari(f'sqrt(({x})+O({p}^{80 if p==2 else 30}))')
def sqf(n):
    s=1 if n>0 else -1; r=1
    for p,e in factorint(abs(int(n))).items():
        if e%2: r*=p
    return s*r
def triple(m,n,cv):
    A,C=(m*m-n*n,m*m+n*n) if cv=='08' else (2*m*n,m*m+n*n)
    while A%2==0 and C%2==0: A//=2; C//=2
    return A,C
def closed(A,C,dl,b_r,b_m):
    b_s=sqf(b_r*b_m); B2=C*C-A*A; val=1
    for p in factorint(B2):
        if p==2:
            if A%2==0: continue
            if (C+A)%8==0: val*=hil(dl,b_s,2)
            elif (C-A)%8!=0: return None
        elif (C+A)%p==0: val*=hil(dl,b_s,p)
    places=set(factorint(abs(A*C)))|({2} if A%2==0 else set())
    for p in places:
        if issq(dl,p): val*=hil(2*psq(dl,p),b_r,p)
        elif issq(-dl,p): val*=hil(-A,b_s,p)*hil(-C,b_m,p)*hil((A-C)*psq(-dl,p),b_r,p)
        else: return None
    return val
if __name__=='__main__':
    res=collections.Counter()
    for fn,cv in (('full0.jsonl','08'),('full1.jsonl','LU')):
        for line in open(fn):
            d=json.loads(line); A,C=triple(d['m'],d['n'],cv)
            for row in d['rows']:
                if row['a'][0]!=1: res[(cv,'a_r != 1 skipped')]+=1; continue
                dl=row['a'][1]
                if 'b' not in (row.get('ct_all') or {}): res[(cv,'row without the full table')]+=1; continue
                for b,v in zip(row['ct_all']['b'],row['ct_all']['v']):
                    c=closed(A,C,dl,b[0],b[1])
                    res[(cv,'agree' if c==v else 'DISAGREE' if c is not None else 'formula undefined')]+=1
    for k,v in sorted(res.items()): print(v,k)
