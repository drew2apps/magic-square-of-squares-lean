#!/usr/bin/env python3
"""Lemma 2 of LAST_PAIR.md for every ratio at once: identities in the variables m, n, d, t over Q(i) (EXPERIMENT).

With s = -i t and A = q+(-i t), A' = q+(i t), B = (m + i n t)(n - i m t), B' = (m - i n t)(n + i m t):
    f1 = d Q1(-i t) = d A B,   f2 = d Q1(i t) = d A' B',   f3 = d Q2(-i t) = d A B',   f4 = d Q2(i t) = d A' B,
    A A' = GL GU / 4,          B B' = F0 F8,                f1 f2 f3 f4 = d^4 (F0 F8 GL GU / 4)^2,
    Q2(s) = -s^4 Q1(1/s)       (so Y^2 = d Q2 is the twist by -1 of Y^2 = d Q1, for every m, n, d),
where F0, F8, GL, GU are four of the five forms of the family in t.  The discriminant of F0 F8 GL GU in t is printed in
factored form: where it is not zero the eight roots are distinct, A, A', B, B' are pairwise coprime and squarefree, a product
of some of the f is a constant times a square only for none or all of them, and the curve of solutions has two components,
each a cover of the line with group (Z/2)^3 branched at 8 points, of genus 9, and an unramified double cover of the curve
u^2 = F0 F8, v^2 = GL GU of genus 5.  usage: lemma2_general.py"""
from second_descent import *                              # K = Q(i), I, and the Sage names
from sage.arith.misc import factor
R = PolynomialRing(K, 'm,n,d,t,s'); m, n, d, t, s = R.gens()
S, D = m**2 + n**2, m**2 - n**2
qp = lambda x: S*x**2 - 2*D*x + S
Q1 = lambda x: qp(x)*(m - n*x)*(n + m*x); Q2 = lambda x: qp(x)*(m + n*x)*(n - m*x)
A, Ac = qp(-I*t), qp(I*t); B, Bc = (m + I*n*t)*(n - I*m*t), (m - I*n*t)*(n + I*m*t)
f1, f2, f3, f4 = d*Q1(-I*t), d*Q1(I*t), d*Q2(-I*t), d*Q2(I*t)
F0, F8, G4 = m**2 + n**2*t**2, n**2 + m**2*t**2, 2*S*(t**2 + 1); GL, GU = G4 - 8*m*n*t, G4 + 8*m*n*t
checks = [('f1 = d A B', f1 == d*A*B), ("f2 = d A' B'", f2 == d*Ac*Bc), ("f3 = d A B'", f3 == d*A*Bc), ("f4 = d A' B", f4 == d*Ac*B),
          ("A A' = GL GU / 4", 4*A*Ac == GL*GU), ("B B' = F0 F8", B*Bc == F0*F8),
          ('f1 f2 f3 f4 = d^4 (F0 F8 GL GU / 4)^2', 16*f1*f2*f3*f4 == d**4*(F0*F8*GL*GU)**2),
          ('Q2(s) = -s^4 Q1(1/s)', all(Q2(s).coefficient({s: k}) == -Q1(s).coefficient({s: 4 - k}) for k in range(5)))]
for name, ok in checks: print('%-40s %s' % (name, ok))
Rq = PolynomialRing(QQ, 'm,n'); mq, nq = Rq.gens(); T = PolynomialRing(Rq, 't'); tt = T.gen()
Sq, Dq = mq**2 + nq**2, mq**2 - nq**2
g = (mq**2 + nq**2*tt**2)*(nq**2 + mq**2*tt**2)*((Sq*(tt**2 + 1))**2 - 16*mq**2*nq**2*tt**2)
print('discriminant in t of F0 F8 GL GU / 4:', factor(g.discriminant()))
print('all identities hold:', all(ok for _, ok in checks))
