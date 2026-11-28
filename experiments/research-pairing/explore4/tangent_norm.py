# Identity behind the Redei form: on the conic u1^2 - d u3^2 = c T^2 with a point Q=(q1,q3,qT),
# the tangent form at Q satisfies  L(P) = -N(qT*(u1+sqrt(d)u3) - T*(q1+sqrt(d)q3)) / (2 qT T)  for P on the conic.
from sympy import symbols, expand, simplify
u1,u3,T,q1,q3,qT,d,c=symbols('u1 u3 T q1 q3 qT d c')
L=q1*u1-d*q3*u3-c*qT*T
N=(qT*u1-T*q1)**2-d*(qT*u3-T*q3)**2
diff=expand(-2*qT*T*L-N)
# reduce with the two conic relations u1^2 - d u3^2 = c T^2 and q1^2 - d q3^2 = c qT^2
diff=expand(diff.subs(u1**2,c*T**2+d*u3**2).subs(q1**2,c*qT**2+d*q3**2))
print('(-2 qT T) L(P) - N(difference) on the conic =',simplify(diff))
