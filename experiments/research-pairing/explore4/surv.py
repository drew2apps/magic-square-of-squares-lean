import sys, json
sys.argv=['x']
from pipeline import *
def rank2(M):
    M=[r[:] for r in M]; rk=0; n=len(M)
    for c in range(n):
        piv=next((i for i in range(rk,n) if M[i][c]),None)
        if piv is None: continue
        M[rk],M[piv]=M[piv],M[rk]
        for i in range(n):
            if i!=rk and M[i][c]: M[i]=[x^y for x,y in zip(M[i],M[rk])]
        rk+=1
    return rk
for m,n,alive,dims in json.load(open('survivors.json')):
    A,B,C=m*m-n*n,2*m*n,m*m+n*n
    print((m,n),'S =',fact(C),'mn sqf',sqf(m*n),'m^2-n^2 sqf',sqf(A))
    for name,e,tors in (('E_08',(0,A*A,C*C),2),('E_LU',(0,B*B,C*C),2),('E_mix',(0,2*(C*C-A*B),2*(C*C+A*B)),2)):
        E=Curve(e); d=len(E.basis)
        G=[[0 if (i==j or cassels(e,triple(E.basis[i]),triple(E.basis[j]),seed=3)[0]==1) else 1 for j in range(d)] for i in range(d)]
        bound=d-rank2(G)-tors
        EP=P.ellinit([0,-(e[0]+e[1]+e[2]),0,e[0]*e[1]+e[0]*e[2]+e[1]*e[2],-e[0]*e[1]*e[2]])
        r=P.ellrank(EP,1)
        print('   ',name,'dim Sel2',d,'pairing rank',rank2(G),'=> rank <=',bound,'| PARI ellrank',[int(r[0]),int(r[1])],'root number',int(P.ellrootno(EP)))
    for t in alive:
        print('    tuple',t,'d0d8 ~',sqf(t[0]*t[1]),' dLdU ~',sqf(t[2]*t[3]),' classes: E_08 delta =',sqf(m*n*t[0]*t[1]),' E_LU (c_r,c_m) =',(sqf(2*C),sqf(2*t[2]*t[3])))
