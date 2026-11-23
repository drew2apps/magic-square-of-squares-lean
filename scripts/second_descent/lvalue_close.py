#!/usr/bin/env python3
"""Second descent with rank 0 read from an L-value.  EXPERIMENT; weaker evidence than second_descent.py.

For a ratio that second_descent.py leaves open: on each curve of the family, a summand E(Q) or E^(-1)(Q) whose 2-descent
bound is positive can still have rank 0 (the bound counts the 2-part of the Tate-Shafarevich group).  If the root number is
+1 and the conductor is small enough, lproof.py computes L(E,1) with an explicit error bound.  When L(E,1) > 0, the rank is
0 (Kolyvagin, with modularity).  Then the equation step of second_descent.py runs with that summand.  A closure found here
is marked "by L-value": it rests on PARI's conductor, root number and coefficients a_n and on the two theorems, not on a descent.

usage: lvalue_close.py LIMIT m,n [m,n ...]      LIMIT: largest conductor for which L(E,1) is computed"""
import sys, json
from math import isqrt
from sage.all__sagemath_schemes import *
import second_descent as sd
import lproof
pari = lproof.pari

def run(m, n, lim):
    Ssq = isqrt(m*m + n*n)**2 == m*m + n*n; res = dict(pair=[m, n], status='open', tried=[])
    for name, (Q, r) in sd.family(m, n).items():
        if Q(sd.sq) == Q(-sd.sq): continue                             # even quartic: the equation is an identity
        cls = [sd.rational_class(sd.K(v), m, n) for v in ((Q(-sd.I), Q(sd.I)) if Ssq else (Q(0), Q.leading_coefficient()))]
        if cls[0] is None or cls[0] != cls[1]: continue
        d = cls[0]; bds = []; info = []
        for tw in (1, -1):
            E, A = sd.weierstrass(tw*d*Q, r); Em = E.minimal_model().pari_curve(); hi = sd.bounds(E)
            N = int(pari.ellglobalred(Em)[0]); rn = int(pari.ellrootno(Em)); rec = dict(curve=name, twist=tw, pari_bound=hi, root_number=rn, conductor='%.3g' % N)
            if hi > 0 and rn == 1 and N <= lim:
                o = lproof.prove([pari(str(c)) for c in E.minimal_model().a_invariants()])      # L(E,1) with an explicit error bound
                rec['L'] = o['L']; rec['tail_bound'] = o['tail_bound']; rec['L_over_Omega'] = o['L_over_Omega']
                if o['positive'] and float(o['L']) > 1e-6: hi = 0; rec['rank0_by'] = 'L-value'
            bds.append(hi); info.append(rec)
        res['tried'] += [i for i in info if 'L' in i]
        if 0 in bds and any(i.get('rank0_by') for i in info):
            o = sd.try_curve(m, n, name, Q, r, bounds_given=tuple(bds))
            if o['status'] == 'closed':
                res['status'] = 'closed by L-value'; res['by'] = name; res['rational_t'] = o['rational_t']; res['closing'] = [i for i in info if i.get('rank0_by')]
                return res
    return res

if __name__ == '__main__':
    lim = float(sys.argv[1])
    for a in sys.argv[2:]:
        m, n = (int(v) for v in a.split(','))
        print(json.dumps(run(m, n, lim)), flush=True)
