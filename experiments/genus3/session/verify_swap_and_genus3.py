# [build agent t12] Committed unchanged apart from this line, from a separate AI session (27 Nov 2026); untrusted; rerun under PARI 2.17.4: output identical (../NOTES.md).
import sympy as sp, cypari2
m,n,u,v=sp.symbols('m n u v')
nrm=lambda p: sp.expand(p[0]**2+p[1]**2)
P=lambda m,n,u,v: [(m*v,n*u),(n*v,m*u),((m+n)*(u-v),(m-n)*(u+v)),((m-n)*(u-v),(m+n)*(u+v))]
A=[nrm(p) for p in P(m,n,u,v)]; B=[nrm(p) for p in P(u,v,m,n)]
print('swap (m,n)<->(u,v) fixes each corner norm:',[sp.simplify(a-b)==0 for a,b in zip(A,B)])
F4=2*(m*m+n*n)*(u*u+v*v); print('sum L+U == 2*F4:',sp.expand(A[2]+A[3]-2*F4)==0,'; F0+F8 == F4/2:',sp.expand(A[0]+A[1]-F4/2)==0)
# Jacobian of y^2=F0F8LU: traces
pari=cypari2.Pari(); pari.allocatemem(10**9)
M,N=627,364; S=M*M+N*N; a=(M*N)**2; b=M**4+N**4; c=S*S; d=2*S*S-16*M*M*N*N
E=lambda s: pari.ellinit(pari.ellfromeqn(pari('y^2-('+s+')')))
Es=[E(f'(x^2-4)*({a}*x+{b})*({c}*x+{d})'),E(f'(x+2)*({a}*x+{b})*({c}*x+{d})'),E(f'(x-2)*({a}*x+{b})*({c}*x+{d})')]
def f(t,p):
    return ((M*M+N*N*t*t)*(N*N+M*M*t*t)*(((M+N)*(t-1))**2+((M-N)*(t+1))**2)*(((M-N)*(t-1))**2+((M+N)*(t+1))**2))%p
ok=[]
for p in [17,23,31,41,43,47,53,59,61,67,71,73]:
    leg=lambda x: 0 if x%p==0 else (1 if pow(x,(p-1)//2,p)==1 else -1)
    lead=(4*M*M*N*N*S*S)%p
    tr=-sum(leg(f(t,p)) for t in range(p)) - leg(lead)     # degree 8: two points at infinity iff lead is a square
    ok.append(tr==sum(int(pari.ellap(e,p)) for e in Es))
print('trace of Frobenius on y^2=F0F8LU equals a_p(E1)+a_p(E2)+a_p(E3) at 12 primes:',ok)
