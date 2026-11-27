#!/usr/bin/env python3
"""Review x4 (run 6): (10512,145): root numbers of E1 for d = +-220962238 (sd_pari's construction, as rank_hard.py), and PARI ellrank of
the g4*K quotient for (10512,145) and its partner (10657,10367) (the F4 curve).  PARI's ellrank is randomised: its lower bounds may differ between runs.  A check, not a certificate.
Written by an AI model; not read by a person.   python3 cond10512.py  (prints only)"""
import sys, os
R = os.environ.get('MSS_REPO', os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
sys.path.insert(0, R + '/experiments/census-wide'); sys.path.insert(0, R + '/scripts/gauss_descent')
sys.path.insert(0, R + '/experiments/census-wide/lvalue')
from sd_pari import pari, family, weierstrass
pari.allocatemem(2 * 10**9)
m, n = 10512, 145
for d in (220962238, -220962238):
    Q, r = family(m, n)['+13']; ain, A = weierstrass(d * Q, r)
    Em = pari.ellminimalmodel(pari.ellinit([pari(x.numerator) / x.denominator for x in ain]))
    w = int(pari.ellrootno(Em)); print('E1, d = %d: root number %d, conductor %s' % (d, w, pari.ellglobalred(Em)[0])); assert w == -1
import lvalue_big as LB
for (a, b) in ((10512, 145), (10657, 10367)):
    E = pari.ellinit(list(LB.model(LB.poly(a, b, 'g4*K'))))
    rk = pari.ellrank(E, 1)
    print('g4*K of (%d,%d): PARI ellrank (effort 1) bounds [%d, %d], root number %d' % (a, b, int(rk[0]), int(rk[1]), int(pari.ellrootno(pari.ellminimalmodel(E)))))
