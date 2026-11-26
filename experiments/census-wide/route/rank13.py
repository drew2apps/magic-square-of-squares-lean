import sys, json
sys.path.insert(0, '..')
from sd_pari import *
for a in sys.argv[1:]:
    m, n, sg = map(int, a.split(','))
    Q, r = family(m, n)['+13']
    d = sg*rational_class(pari.subst(Q, 's', -I), m, n)
    ain, A = weierstrass(d*Q, r); a = [pari(x.numerator)/x.denominator for x in ain]
    Em = pari.ellminimalmodel(pari.ellinit(a))
    iso = pari.ellisomat(Em, 0, 1)
    out = []
    for c in iso[0]:
        F = pari.ellminimalmodel(pari.ellinit(c))
        for eff in (1, 3):
            b = pari.ellrank(F, eff); out.append((eff, int(b[0]), int(b[1])))
    print(m, n, sg, 'cond', pari.factor(pari.ellglobalred(Em)[0])[0], out, flush=True)
