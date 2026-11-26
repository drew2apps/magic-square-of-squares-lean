"""x3 referee: controls for dm_indep.py / index_indep.py on (7735,4968): each must come out as stated.
(1) a classifier forged to call t = 1 a solution must give ALERT; (2) 3G in place of G must fail saturation at 3;
(3) an x-relation forged to x x' = 2 c3^2 must not be used silently (the |a| = |b| quadratic changes; checked by
recomputing x(P) x(P + T) for P = G, which must NOT equal 2 c3^2).  Written by an AI model; not read by a person."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
import dm_indep as D, index_indep as I
ok = []
keep = D.classify
D.classify = lambda m, n, t: 'SOLUTION' if t == 1 else keep(m, n, t)
r = D.run(7735, 4968); print('(1) forged classifier: status', r['status'], '(must be ALERT)'); ok.append(r['status'] == 'ALERT')
D.classify = keep
r = D.run(7735, 4968); C, G, T = r['C'], r['G'], r['T']
G3 = C.add(G, C.add(G, G))
try:
    I.saturation(C, G3, T, lmax=7); print('(2) 3G: NOT refused'); ok.append(False)
except AssertionError as ex:
    print('(2) 3G refused:', str(ex)[:60]); ok.append('3' in str(ex))
P, _ = D.model(7735, 4968); c3 = P[3]; GT = C.add(G, T)
print('(3) x(G) x(G+T) == c3^2 ?', G[0] * GT[0] == c3 * c3, '(must be False: G is not a solution point)'); ok.append(G[0] * GT[0] != c3 * c3)
print('(4) unchanged run: status', r['status'], 'hits', r['hits']); ok.append(r['status'] == 'closed' and set(r['hits']) == {'1', '-1'})
print('controls:', 'all as required (%d of %d)' % (sum(ok), len(ok)) if all(ok) else 'FAILED')
