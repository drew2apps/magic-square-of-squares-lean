"""z15: as z15_findpoint.py (same calls, so PARI's fixed random seed gives the same point), then print the point, check
it, and carry it to E_LU of the partner ratio (5623,1017) through the two minimal models (PARI's ellminimalmodel; the
image point is checked on the partner's curve); report whether the partner's surviving class is the class of a
rational point.  Written by an AI model (reviewer z15), not read by a person."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('MSS_REPO', os.path.abspath(os.path.join(HERE, '..', '..')))
sys.argv = [sys.argv[0]]; sys.path.insert(0, HERE)
import io, contextlib
import t10_pairct_copy as T
from fractions import Fraction as Fr
RS = T.RS; P = T.PARI
with contextlib.redirect_stdout(io.StringIO()):
    Pd = T.PairData(3320, 2303, '0', '8'); Pd2 = T.PairData(5623, 1017, 'L', 'U')
E = P.ellinit([0, Pd.Pr.A[0], 0, Pd.Pr.A[1], Pd.Pr.A[2]])
pts = []
for eff in (4, 5, 6):
    r = P.ellrank(E, eff); pts = list(r[3])
    if pts: break
assert pts, 'no point this time'
Q = pts[0]; x, y = Fr(str(Q[0])), Fr(str(Q[1])); assert Pd.Pr.on((x, y))
print('(3320,2303) E_08 (t10 model, A =', Pd.Pr.A, '): point found with ellrank effort', eff, ': x =', x)
print('   canonical height', P.ellheight(E, Q), '; order', P.ellorder(E, Q))
v = Pd.vec(Pd.cls_of_tuple((2, 1, 1, 635399))); im = Pd.imgv + [Pd.vec(Pd.kum(x))]
print('   ghost class (1, 19505) in the span of torsion and this point:', len(RS.f2_span(im + [v])) == len(RS.f2_span(im)))
E2 = P.ellinit([0, Pd2.Pr.A[0], 0, Pd2.Pr.A[1], Pd2.Pr.A[2]])
m1 = P("ellminimalmodel(ellinit(%s), &z15v1)" % str([int(a) for a in E[:5]])); v1 = P('z15v1')
m2 = P("ellminimalmodel(ellinit(%s), &z15v2)" % str([int(a) for a in E2[:5]])); v2 = P('z15v2')
same = list(m1[:5]) == list(m2[:5])
print('(5623,1017) E_LU and (3320,2303) E_08 have the same minimal model:', same)
if same:
    Qm = P.ellchangepoint(Q, v1); Q2 = P.ellchangepointinv(Qm, v2)
    x2, y2 = Fr(str(Q2[0])), Fr(str(Q2[1])); assert Pd2.Pr.on((x2, y2))
    v = Pd2.vec(Pd2.cls_of_tuple((2, 1270798, 1, 2))); im = Pd2.imgv + [Pd2.vec(Pd2.kum(x2))]
    print('   image point on (5623,1017) E_LU checked on the curve; surviving class (16326209, 1) in the span of torsion and it:',
          len(RS.f2_span(im + [v])) == len(RS.f2_span(im)))
