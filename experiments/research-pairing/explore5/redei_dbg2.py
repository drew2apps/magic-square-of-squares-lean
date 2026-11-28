import sys, random, collections
exec(open('redei_k.py').read().split("rng=random.Random(11); cnt=")[0])
rng=random.Random(3)
def place_by_place(e,d,beta):
    a=(1,d,d); nf=gp.nfinit(gp(f'y^2-({d})'))
    def el(a_,b_): return gp(f'Mod(({a_.numerator}/{a_.denominator})+({b_.numerator}/{b_.denominator})*y,y^2-({d}))')
    x0,y0,z0=[int(c) for c in gp.qfsolve(gp.matdiagonal([1,-d,-beta]))]; pi=el(Fr(x0,z0),Fr(y0,z0))
    c2,c3=e[2]-e[0],e[1]-e[0]       # u1^2-d u3^2=c2 ; u1^2-d u2^2=c3
    q=[]
    for c in (c2,c3):
        s=[int(t) for t in gp.qfsolve(gp.matdiagonal([1,-d,-c]))]; q.append((Fr(s[0],s[2]),Fr(s[1],s[2])))
    S=set()
    for n in (2,d,beta,c2,c3,e[2]-e[1],z0,x0*x0-d*y0*y0)+tuple(t.denominator*t.numerator for qq in q for t in qq):
        if abs(n)>1: S|=set(primes_of(abs(n)))
    tot_l=tot_r=tot_g=1
    for p in sorted(S):
        x=local_point(e,a,p,rng); N=80 if p==2 else 30
        u=[padic_sqrt_rat((x-e[i])/a[i],p,N) for i in range(3)]
        L2=q[0][0]*u[0]-d*q[0][1]*u[2]-c2; L3=q[1][0]*u[0]-d*q[1][1]*u[1]-c3
        lhs=int(gp.hilbert(gp(f'{(L2*L3).numerator}/{(L2*L3).denominator}'),beta,p))
        g=el(u[0],u[2])*el(u[0],u[1]); gq=el(q[0][0],q[0][1])*el(q[1][0],q[1][1])
        r=1; rg=1
        for pr in gp.idealprimedec(nf,p):
            r*=int(gp.nfhilbert(nf,g/gq,pi,pr)); rg*=int(gp.nfhilbert(nf,g,pi,pr))
        const=int(gp.hilbert(4*c2*c3,beta,p))   # (2c2d*2c3d, beta) = (4 c2 c3, beta)
        print('  p=%d  (L2L3,beta)=%d  const=%d  prod_w(gamma_w/gamma,pi)=%d   prod_w(gamma_w,pi)=%d'%(p,lhs,const,r,rg), ' agree' if lhs*const==r else ' DIFFER')
        tot_l*=lhs; tot_r*=r; tot_g*=rg
    print('  totals over finite places: Cassels',tot_l,' ratio form',tot_r,' gamma_w form',tot_g, ' ref engine',cassels(e,tr3((1,d)),tr3((1,beta)),seed=7)[0], 'q=',q,'pi=',(x0,y0,z0))
for e,a,b in (((-131,-130,-4),(1,2),(1,14)),((-131,-130,-4),(1,42),(1,7))):
    print(e,a,b); place_by_place(e,a[1],b[1])
