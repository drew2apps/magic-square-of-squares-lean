# [build agent t12] Committed unchanged apart from this line, from a separate AI session (27 Nov 2026); untrusted; rerun under PARI 2.17.4: output identical (../NOTES.md).
import cypari2, sys
pari = cypari2.Pari(); pari.allocatemem(2*10**9)
m,n = 627,364
pari(f"m={m};n={n};")
pari("nf=nfinit(w^2+1); ii=Mod(w,w^2+1); p=m+ii*n; pb=m-ii*n;")
pari("Q(s)=(m-n*s)*(n+m*s)*(p*s-pb)*(pb*s+p);")
pari("Qc(s)=(m-n*s)*(n+m*s)*(pb*s-p)*(p*s+pb);")
pari("E=ellinit(ellfromeqn(y^2-Q(x)),nf); Ec=ellinit(ellfromeqn(y^2-Qc(x)),nf);")
print("j =", pari("lift(E.j)"))
print("jc=", pari("lift(Ec.j)"))
print("conj check:", pari("lift(Ec.j)==subst(lift(E.j),w,-w)"))
print("tors", pari("elltors(E)"))
print("cond norm", pari("idealnorm(nf,ellglobalred(E)[1])"), pari("factor(idealnorm(nf,ellglobalred(E)[1]))"))
# isogeny class
r = pari("M=ellisomat(E,2,1); [vector(#M[1],k,lift(ellinit(M[1][k],nf).j)), M[2]]")
js = r[0]
print("class size", len(js))
jc = pari("lift(Ec.j)")
for k,j in enumerate(js):
    print(k, "is conj j" if j==jc else "", str(j)[:80])
print(r[1])
