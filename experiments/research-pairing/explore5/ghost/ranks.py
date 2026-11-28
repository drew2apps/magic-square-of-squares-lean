"""3(c): ranks of E_08, E_LU, E_mix at the ghost ratios (PARI ellrank, effort 2) and whether the ghost class of each
surviving tuple is the image of a rational point (span of the images of the points PARI found and of the 2-torsion)."""
import sys,json
sys.argv=['x']; sys.path.insert(0,'/tmp/rn')
from pipeline import *
from cassels import tors_image
from fractions import Fraction as Fr
def vec(c):
    s=set()
    for i,x in enumerate(c):
        x=int(x)
        if x<0: s.add((i,-1))
        for q,e in fact(x).items():
            if e%2: s.add((i,q))
    return frozenset(s)
def inspanF2(gens,v):
    basis=[]  # list of (pivot, vec)
    def red(x):
        for piv,b in basis:
            if piv in x: x=x^b
        return x
    for g in gens:
        g=red(g)
        if g: basis.append((sorted(g,key=str)[0],g))
    # re-reduce fully
    v=red(v)
    while True:
        w=red(v)
        if w==v: break
        v=w
    return len(v)==0
def ptclass(e,x):
    x=Fr(int(x.numerator()),int(x.denominator()))
    out=[]
    for i in (0,1):
        t=(x-e[i]); out.append(t.numerator*t.denominator)
    return out
S=json.load(open('/tmp/rn/survivors.json'))+[[24599,5207,[[5207,24599,2,1],[49198,10414,2,1]],{}]]
for m,n,alive,_ in S:
    A,B,C=m*m-n*n,2*m*n,m*m+n*n
    cs=curves(m,n)
    print((m,n),flush=True)
    for name,e in (('E_08',(0,A*A,C*C)),('E_LU',(0,B*B,C*C)),('E_mix',(0,2*(C*C-A*B),2*(C*C+A*B)))):
        EP=P.ellinit([0,-(e[0]+e[1]+e[2]),0,e[0]*e[1]+e[0]*e[2]+e[1]*e[2],-e[0]*e[1]*e[2]])
        r=P.ellrank(EP,2); lo,hi=int(r[0]),int(r[1]); pts=[p for p in r[3]]
        gens=[vec(tors_image(e,i)[:2]) for i in range(3)]
        ok=True
        for p in pts:
            x=p[0]
            if any(x==ei for ei in e): continue
            gens.append(vec(ptclass(e,x)))
        tors=str(P.elltors(EP)[1])
        line='   %s rank in [%d,%d], %d points found, torsion %s, root no %d'%(name,lo,hi,len(pts),tors,int(P.ellrootno(EP)))
        keys=['08'] if name=='E_08' else ['LU'] if name=='E_LU' else ['0L','8U','0U','8L']
        for t in alive:
            for k in keys:
                E,f=cs[k]; c=f(tuple(t)); c=(sqf(c[0]),sqf(c[1]))
                # torsion-only?
                tonly=inspanF2(gens[:3],vec(c))
                # 4-torsion images too
                line+='\n      tuple %s on %s: class %s; in span of 2-torsion images: %s; in span of images of found points + 2-torsion: %s'%(t,k,c,tonly,inspanF2(gens,vec(c)))
        print(line,flush=True)
