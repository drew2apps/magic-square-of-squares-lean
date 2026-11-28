"""3(c), second pass.  For each ghost ratio and E in {E_08, E_LU}: Sel^2, the kernel K of the Cassels-Tate pairing on it,
the span T of the images of the rational torsion (Z/4 x Z/2: dimension 2), whether the ghost class lies in T, and a point
search (PARI hyperellratpoints) on the 2-covering of the ghost class
      u1^2 - a u2^2 = (e2-e1),  u1^2 - a u3^2 = (e3-e1)      (class (1,a,a) at (e1,e2,e3), x - e1 = u1^2)
written as the quartic  w^2 = a [ E2 (1 + a s^2)^2 - E3 (1 - a s^2)^2 ],  E2 = e2-e1 (a square), E3 = e3-e1."""
import sys,json
sys.argv=['x']; sys.path.insert(0,'/tmp/rn')
from pipeline import *
from cassels import tors_image
H=int(sys.argv[1]) if len(sys.argv)>1 else 100000
def vec(c):
    s=set()
    for i,x in enumerate(c):
        x=int(x)
        if x<0: s.add((i,-1))
        for q,e in fact(x).items():
            if e%2: s.add((i,q))
    return frozenset(s)
def span(gens):
    S={frozenset()}
    for g in gens: S|={s^g for s in S}
    return S
S=json.load(open('/tmp/rn/survivors.json'))+[[24599,5207,[[5207,24599,2,1],[49198,10414,2,1]],{}]]
for m,n,alive,_ in S:
    A,B,C=m*m-n*n,2*m*n,m*m+n*n
    cs=curves(m,n); print((m,n),'C =',fact(C),flush=True)
    for name,e,key in (('E_08',(0,A*A,C*C),'08'),('E_LU',(0,B*B,C*C),'LU')):
        E,f=cs[key]; d=len(E.basis)
        G=[[0 if (i==j or cassels(e,triple(E.basis[i]),triple(E.basis[j]),seed=3)[0]==1) else 1 for j in range(d)] for i in range(d)]
        # kernel of the pairing
        ker=[]
        for mask in range(1<<d):
            if all(sum(G[i][j] for j in range(d) if mask>>j&1)%2==0 for i in range(d)):
                c1=c2=1
                for j in range(d):
                    if mask>>j&1: c1*=E.basis[j][0]; c2*=E.basis[j][1]
                ker.append(vec((sqf(c1),sqf(c2))))
        ker=set(ker)
        EP=P.ellinit([0,-(e[0]+e[1]+e[2]),0,e[0]*e[1]+e[0]*e[2]+e[1]*e[2],-e[0]*e[1]*e[2]])
        tg=[vec(tors_image(e,i)[:2]) for i in range(3)]
        for T in P.elltors(EP)[2]:
            x=T[0]
            if not any(x==ei for ei in e): tg.append(vec((int(x-e[0]),int(x-e[1]))))
        Ts=span(tg)
        import math
        dk=int(math.log2(len(ker))); dt=int(math.log2(len(Ts)))
        g=sorted(set(tuple(sorted(vec(tuple(sqf(x) for x in f(tuple(t)))),key=str)) for t in alive))
        gv=[frozenset(x) for x in g]
        r=P.ellrank(EP,1)
        print('   %s: dim Sel2 %d, dim kernel of pairing %d, dim torsion image %d => rank <= %d; ellrank [%d,%d]; root number %d'%(name,d,dk,dt,dk-dt,int(r[0]),int(r[1]),int(P.ellrootno(EP))))
        for v in gv:
            print('      ghost class in kernel: %s; in torsion image: %s'%(v in ker, v in Ts))
            if v in Ts: continue
            # reduce the ghost class by torsion to the form (1,a) if possible and search
            cands=[v^t for t in Ts]
            best=[c for c in cands if all(i==1 for i,q in c)]
            for c in best:
                a=1
                for i,q in c: a*=q
                E2=e[1]-e[0]; E3=e[2]-e[0]
                pol='(%d)*((%d)*(1+(%d)*x^2)^2-(%d)*(1-(%d)*x^2)^2)'%(a,E2,a,E3,a)
                pts=P.hyperellratpoints(P(pol),H)
                # second model: parametrise the conic u1^2 - a u3^2 = E3 from (C,0)
                pol2='(%d)*((%d)*(1+(%d)*x^2)^2-(%d)*(1-(%d)*x^2)^2)'%(a,E3,a,E2,a)
                pts2=P.hyperellratpoints(P(pol2),H)
                print('      covering of class (1,%d,%d): rational points of height <= %d on the two quartic models: %d, %d  %s'%(a,a,H,len(pts),len(pts2),str(pts)[:200] if len(pts) else ''),flush=True)
