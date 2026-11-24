#!/usr/bin/env python3
"""Second step for the perfect triangle: for each 4 of the 5 points a_j in {0,1,-1,2,1/2}, the curve Y^2 = D * prod (z - a_j),
D = +1 and -1 (the class of the product of four forms at a degenerate tuple g*(1,1,1,1,2) is 1 up to squares, the form 2z-1 contributing 2*2).
Rank over Q(i) = rank(D=1) + rank(D=-1) over Q.  PARI ellrank gives [lower, upper, ...]."""
from itertools import combinations
from fractions import Fraction as F
import cypari2
pari = cypari2.Pari(); pari.allocatemem(10**9)
A = [F(0), F(1), F(-1), F(2), F(1, 2)]; names = ['z', 'z-1', 'z+1', 'z-2', '2z-1']
for Q in combinations(range(5), 4):
    a = [A[j] for j in Q]; a4 = a[3]
    # z = a4 + 1/t :  t^4 Y^2 = D * c * prod_{j<3} (t - t_j),  t_j = -1/(a4 - a_j),  c = prod (a4 - a_j)
    tj = [-1/(a4 - x) for x in a[:3]]; c = (a4 - a[0])*(a4 - a[1])*(a4 - a[2])
    out = []
    for D in (1, -1):
        tw = D*c
        e = [tw*t for t in tj]                                   # y^2 = prod (x - tw t_j)
        a2 = -(e[0] + e[1] + e[2]); a4c = e[0]*e[1] + e[0]*e[2] + e[1]*e[2]; a6 = -e[0]*e[1]*e[2]
        E = pari.ellinit([0, pari(str(a2)), 0, pari(str(a4c)), pari(str(a6))])
        Em = pari.ellminimalmodel(E)
        rk = pari.ellrank(Em)
        tors = pari.elltors(Em)
        out.append('D=%+d: cond %s, rank bounds [%s, %s], torsion %s' % (D, pari.ellglobalred(Em)[0], rk[0], rk[1], tors[1]))
    print('{%s}: %s' % (', '.join(names[j] for j in Q), ' | '.join(out)))
