from fractions import Fraction as Fr
from math import isqrt
import itertools
class G:
    __slots__=('a','b')
    def __init__(s,a,b=0): s.a=Fr(a); s.b=Fr(b)
    def __add__(s,o): o=g(o); return G(s.a+o.a,s.b+o.b)
    def __sub__(s,o): o=g(o); return G(s.a-o.a,s.b-o.b)
    def __mul__(s,o): o=g(o); return G(s.a*o.a-s.b*o.b,s.a*o.b+s.b*o.a)
    def inv(s): d=s.a*s.a+s.b*s.b; return G(s.a/d,-s.b/d)
    def __truediv__(s,o): return s*g(o).inv()
    def __neg__(s): return G(-s.a,-s.b)
    def key(s): return (s.a,s.b)
    def zero(s): return s.a==0 and s.b==0
def g(x): return x if isinstance(x,G) else G(x)
def rsq(q):
    if q<0: return None
    n,d=q.numerator,q.denominator; rn,rd=isqrt(n),isqrt(d)
    return Fr(rn,rd) if rn*rn==n and rd*rd==d else None
def is_sq(z):
    if z.zero(): return True
    N=rsq(z.a*z.a+z.b*z.b)
    if N is None: return False
    return any(rsq((z.a+e*N)/2) is not None and (z.b==0 or rsq((z.a+e*N)/2)!=0 or True) and rsq((e*N-z.a)/2) is not None for e in (1,-1))
m,n=627,364; I=G(0,1); p=G(m,n); pb=G(m,-n); r=G(Fr(m,n)); w=p/pb
F={'L1':lambda s: G(m)-s*n,'L2':lambda s: G(m)+s*n,'L3':lambda s: G(n)+s*m,'L4':lambda s: G(n)-s*m,
   'M1':lambda s: p*s-pb,'M2':lambda s: pb*s+p,'M3':lambda s: p*s+pb,'M4':lambda s: pb*s-p}
quartics=[('L1','L3','M1','M2'),('L1','L3','M3','M4'),('L2','L4','M1','M2'),('L2','L4','M3','M4')]
seeds=[G(0),G(1),I,r,w,r*w,r*I,w*I,G(2),G(Fr(1,2)),G(1,1),G(3),r*r,w*w,r+I,r+1,w+1,w+I,G(Fr(m+n,m-n)),p,pb,p/n,p/m]
maps=[lambda s:-s,lambda s:s.inv() if not s.zero() else None,lambda s:(G(1)-s)/(G(1)+s) if not (G(1)+s).zero() else None,lambda s:s*I,
      lambda s:(s-I)/(s+I) if not (s+I).zero() else None, lambda s:s*w, lambda s:s*r]
pts={s.key():s for s in seeds}
for it in range(3):
    new={}
    for s in list(pts.values()):
        for f in maps:
            try: t=f(s)
            except ZeroDivisionError: t=None
            if t is not None and t.key() not in pts and t.a.denominator<10**40: new[t.key()]=t
    pts.update(new)
    if len(pts)>40000: break
print('candidates',len(pts))
hits=[]
for s in pts.values():
    for q in quartics:
        v=G(1)
        for f in q: v=v*F[f](s)
        if not v.zero() and is_sq(v): hits.append((q,s.a,s.b))
print('square values (non-root):',len(hits))
for h in hits[:20]: print(h)
