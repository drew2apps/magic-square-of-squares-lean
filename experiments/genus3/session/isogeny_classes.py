# [build agent t12] Committed unchanged apart from this line, from a separate AI session (27 Nov 2026); untrusted; rerun under PARI 2.17.4: output identical (../NOTES.md).
import cypari2, itertools, sys
pari = cypari2.Pari(); pari.allocatemem(2*10**9)
m,n = (int(sys.argv[1]),int(sys.argv[2])) if len(sys.argv)>2 else (627,364)
S=m*m+n*n
x=pari('x')
forms={'F0':f'({m*m}+{n*n}*x^2)','F8':f'({n*n}+{m*m}*x^2)','F4':f'({2*S})*(x^2+1)',
 'L':f'({(m+n)**2}*(x-1)^2+{(m-n)**2}*(x+1)^2)','U':f'({(m-n)**2}*(x-1)^2+{(m+n)**2}*(x+1)^2)'}
def E_of(polystr):
    return pari.ellminimalmodel(pari.ellinit(pari.ellfromeqn(pari('y^2-('+polystr+')'))))
C={}
for a,b in itertools.combinations(forms,2):
    C[a+'*'+b]=E_of(forms[a]+'*'+forms[b])
A=(m*n)**2; B=m**4+n**4; Cc=S*S; D=2*S*S-16*m*m*n*n
C['E1']=E_of(f'(x^2-4)*({A}*x+{B})*({Cc}*x+{D})')
C['E2']=E_of(f'(x+2)*({A}*x+{B})*({Cc}*x+{D})')
C['E3']=E_of(f'(x-2)*({A}*x+{B})*({Cc}*x+{D})')
def Epm(m,n):
    S=m*m+n*n
    return E_of(f'(x+{2*S**3})*(x^2+{16*m*m*n*n*(m+n)**2}*x+{16*m*m*n*n*(m+n)**4*S*S})')
C['E+']=Epm(m,n); C['E-']=Epm(m,-n)
# with F4: genus-3 curves from 4 forms incl F4 that are even in t: F4*F0*F8 is a triple; quartets: F4 F0 L U, F4 F8 L U (even in t)
# y^2 = h(x), x=t^2 (genus 1) and y^2 = x h(x) (genus 2)
LU=f'({4*S*S}*(x+1)^2-{64*m*m*n*n}*x)'
C['[F4F0LU|x=t^2]']=E_of(f'{2*S}*(x+1)*({m*m}+{n*n}*x)*'+LU)
C['[F4F8LU|x=t^2]']=E_of(f'{2*S}*(x+1)*({n*n}+{m*m}*x)*'+LU)
C['[F0F8LU|x=t^2]']=E_of(f'({m*m}+{n*n}*x)*({n*n}+{m*m}*x)*'+LU)
C['[F4F0F8|x(..)]']=E_of(f'{2*S}*x*(x+1)*({m*m}+{n*n}*x)*({n*n}+{m*m}*x)')  # y^2 = x F4F0F8 in x: genus 1 (deg 4)
C['[F4LU|x(..)]']=E_of(f'{2*S}*x*(x+1)*'+LU)
C['[F4LU|x=t^2]q']=E_of(f'{2*S}*(x+1)*'+LU+'*1')
C['[F0LU|x(..)]']=E_of(f'x*({m*m}+{n*n}*x)*'+LU)
C['[F8LU|x(..)]']=E_of(f'x*({n*n}+{m*m}*x)*'+LU)
bad=2*3*m*n*(m*m-n*n)*S*(m**8+14*m**4*n**4+n**8)
ps=[p for p in pari.primes(120) if bad%int(p)!=0][:40]
vec={k:tuple(int(pari.ellap(E,p)) for p in ps) for k,E in C.items()}
names=list(C)
print('ratio',(m,n),'primes used',len(ps))
grp={}
for k in names: grp.setdefault(vec[k],[]).append(k)
for v,ks in grp.items():
    E=C[ks[0]]; 
    try: r=pari.ellrank(E,1); rk=[int(r[0]),int(r[1])]
    except Exception: rk='?'
    print(ks,'rank',rk,'tors',pari.elltors(E)[1])
# twists
for a,b in itertools.combinations(names,2):
    if vec[a]!=vec[b] and all(abs(u)==abs(v) for u,v in zip(vec[a],vec[b])): print('twist pair',a,b)
