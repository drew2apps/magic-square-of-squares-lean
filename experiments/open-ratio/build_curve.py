# Rebuild the (627,364) bound-2 curve over Q(i): delta*Y^2 = Q(s), Q = L1*L3*M1*M2 (labels of kclose.py)
import cypari2; pari = cypari2.Pari(); pari.allocatemem(2*10**9)
print(pari('default(version)') if False else pari.version())
gp = pari
gp('m=627; n=364; p=m+I*n; pb=m-I*n;')
gp('L1=m-n*s; L3=n+m*s; M1=p*s-pb; M2=pb*s+p; Q=L1*L3*M1*M2;')
print('Q =', gp('Q'))
print('Q(-I)=', gp('subst(Q,s,-I)'), ' Q(I)=', gp('subst(Q,s,I)'))
print('factor Q(-I) in Z[i]:', gp('K=nfinit(y^2+1); idealfactor(K, subst(subst(Q,s,-I),I,y))'))
