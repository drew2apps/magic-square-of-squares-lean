# [build agent t12] Committed unchanged apart from this line, from a separate AI session (27 Nov 2026); untrusted; rerun under PARI 2.17.4: output identical (../NOTES.md).
import cypari2, sys
pari = cypari2.Pari(); pari.allocatemem(3*10**9)
OPEN=[(195,28),(223,167),(440,279),(719,161),(609,200),(809,409),(627,364),(991,263),(861,620)]
def E_of(polystr):
    return pari.ellminimalmodel(pari.ellinit(pari.ellfromeqn(pari('y^2-('+polystr+')'))))
for (m,n) in OPEN:
    S=m*m+n*n; A=(m*n)**2; B=m**4+n**4; C=S*S; D=2*S*S-16*m*m*n*n
    cs={'E1':f'(x^2-4)*({A}*x+{B})*({C}*x+{D})','E2':f'(x+2)*({A}*x+{B})*({C}*x+{D})','E3':f'(x-2)*({A}*x+{B})*({C}*x+{D})'}
    for name,ps in cs.items():
        E=E_of(ps)
        r0=pari.ellrank(E,1); lo0,hi0=int(r0[0]),int(r0[1])
        if lo0==hi0: print((m,n),name,'rank',lo0); continue
        M=pari.ellisomat(E,2,1)[0]
        best=(lo0,hi0); det=[]
        for k in range(len(M)):
            Ek=pari.ellinit(M[k])
            try:
                r=pari.ellrank(Ek,1); lo,hi=int(r[0]),int(r[1])
            except Exception as ex:
                lo,hi=0,99
            det.append((lo,hi,str(pari.elltors(Ek)[1])))
            best=(max(best[0],lo),min(best[1],hi))
        print((m,n),name,'class size',len(M),'bounds per curve',det,'=> rank in',list(best),'root',pari.ellrootno(E), flush=True)
