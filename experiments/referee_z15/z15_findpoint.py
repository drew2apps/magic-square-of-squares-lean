"""z15: try to find a rational point of E_08 at (3320,2303) (t10's integral model) with the ghost class, by PARI's
ellrank at higher effort (a finder; the point is checked on the curve and its class computed by t10's Kummer map).
E_LU at (5623,1017) is the same curve up to the scaling x -> 4x (partner ratios).  Written by an AI model (reviewer
z15), not read by a person."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('MSS_REPO', os.path.abspath(os.path.join(HERE, '..', '..')))
sys.argv = [sys.argv[0]]; sys.path.insert(0, HERE)
import io, contextlib
import t10_pairct_copy as T
from fractions import Fraction as Fr
RS = T.RS; P = T.PARI
with contextlib.redirect_stdout(io.StringIO()): Pd = T.PairData(3320, 2303, '0', '8')
E = P.ellinit([0, Pd.Pr.A[0], 0, Pd.Pr.A[1], Pd.Pr.A[2]])
v = Pd.vec(Pd.cls_of_tuple((2, 1, 1, 635399)))
for eff in (4, 5, 6):
    r = P.ellrank(E, eff)
    pts = [(Fr(str(Q[0])), Fr(str(Q[1]))) for Q in r[3]]
    for Q in pts: assert Pd.Pr.on(Q)
    imgv = Pd.imgv + [Pd.vec(Pd.kum(Q[0])) for Q in pts]
    s1 = len(RS.f2_span(imgv)); s2 = len(RS.f2_span(imgv + [v]))
    print('effort %d: %s, %d point(s), images span %d; ghost class %s' % (eff, [int(r[0]), int(r[1])], len(pts), s1, 'IS the class of a point found' if s1 == s2 else 'not in the span'), flush=True)
    if pts: break
