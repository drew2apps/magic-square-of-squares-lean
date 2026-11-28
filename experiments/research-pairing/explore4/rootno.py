import cypari2, random, collections, re
from math import gcd, isqrt
pari=cypari2.Pari(); pari.allocatemem(2*10**9)
def w(m,n):
    A,B,C=m*m-n*n,2*m*n,m*m+n*n
    EA=pari.ellinit([0,-(A*A+C*C),0,A*A*C*C,0]); EB=pari.ellinit([0,-(B*B+C*C),0,B*B*C*C,0])
    Em=pari.ellinit(pari.ellfromeqn(pari(f'y^2-(x+{2*C*C})*(x-{2*A*B})*(x+{2*A*B})')))
    return int(pari.ellrootno(EA)),int(pari.ellrootno(EB)),int(pari.ellrootno(Em))
def load(fn): return [tuple(int(x) for x in l.split()[0].split(',')) for l in open(fn) if re.match(r'\d+,\d+',l)]
def is_sq(x): r=isqrt(x); return r*r==x
def tab(pairs,label):
    c=collections.Counter(); cm=collections.Counter()
    for m,n in pairs:
        a,b,x=w(m,n); c[(a,b)]+=1; cm[x]+=1
    N=len(pairs)
    print(label,N,'| (w_A,w_B):',{k:'%d (%.1f%%)'%(v,100*v/N) for k,v in sorted(c.items())},'| w_mix:',dict(cm))
K1=load('g2_997_all.txt'); K2=load('g2_15000_all.txt')
tab(K1,'class-keeping non-degenerate, m<=5000'); 
rng=random.Random(1); tab(rng.sample(K2,2500),'class-keeping non-degenerate, 5000<m<=15000 (sample)')
need=[tuple(int(x) for x in l.split()[0].split(',')) for l in open('g2_997_all.txt') if 'left' in l]
tab(need,'  of which bad primes leave tuples (m<=5000)')
def rand(H,N,cond):
    out=[]
    while len(out)<N:
        m=rng.randint(2,H); n=rng.randint(1,m-1)
        S=m*m+n*n
        if gcd(m,n)==1 and not is_sq(S) and not (S%2==0 and is_sq(S//2)) and cond(m,n): out.append((m,n))
    return out
tab(rand(15000,2500,lambda m,n: True),'random non-degenerate, m<=15000')
tab(rand(15000,2500,lambda m,n: (m*n)%4==0 or (m-n)%8==0 or (m+n)%8==0),'random non-degenerate with a 2-adic point')
