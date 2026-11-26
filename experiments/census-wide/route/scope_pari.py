#!/usr/bin/env python3
"""PARI port of the repo's two_map_scope.py (and its -d twist variant): for a square-type ratio, the class d of the '+13'
quartic, the rank bounds of E1 (Jacobian of Y^2 = d Q(s)) and of its -1 twist, and whether the degenerate point G = P + conj(P)
(P at s = -i, t = 1) is rational of infinite order.  usage: scope_pari.py m,n ...  EXPERIMENT."""
import sys, json
sys.path.insert(0, '..')
from sd_pari import *
def scope(m, n, sign=1):
    Q, r = family(m, n)['+13']
    cls = [rational_class(pari.subst(Q, 's', v), m, n) for v in (-I, I)]
    if cls[0] is None or cls[0] != cls[1]: return dict(pair=[m, n], status='class not rational')
    d = sign*cls[0]
    ain, A = weierstrass(d*Q, r); a = [pari(x.numerator)/x.denominator for x in ain]
    E = pari.ellinit(a); Em = pari.ellminimalmodel(E)
    rk = pari.ellrank(Em, 1); lo, hi = int(rk[0]), int(rk[1])
    EK = pari.ellinit(a, nf)
    v = pari.lift(d*pari.subst(Q, 's', -I))
    rts = pari.nfroots(nf, X**2 - v)
    res = dict(pair=[m, n], d=int(d), rank=[lo, hi])
    if len(rts) == 0: res['G'] = 'P not over K'; return res
    Y = rts[0]; rr = pari(r.numerator)/r.denominator; Ap = pari(A.numerator)/A.denominator
    w = 1/(-I - rr); P = [pari.lift(Ap*w), pari.lift(Ap*Y*w*w)]
    cP = [pari.subst(c, 'y', -pari('y')) for c in P]
    G = pari.elladd(EK, P, cP)
    if len(G) == 1: res['G'] = 'O'; return res
    g = [pari.simplify(pari.lift(c)) for c in G]
    if any(pari.poldegree(c, 'y') > 0 for c in g): res['G'] = 'not rational'; return res
    res['G'] = [str(c) for c in g]; res['G_order'] = int(pari.ellorder(E, g))
    res['hG'] = float(pari.ellheight(E, g)) if res['G_order'] == 0 else None
    return res
if __name__ == '__main__':
    for a in sys.argv[1:]:
        m, n = map(int, a.split(','))
        for sg in (1, -1): print(json.dumps(scope(m, n, sg)), flush=True)
