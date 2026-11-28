"""z15: the curves carrying the non-torsion ghost class at the other three survivor ratios (E_08 at (9113,3200) and
(14903,9696), both root number +1 and rank bound 2; E_LU at (12313,5913) is the partner of the first): PARI's ellrank at
rising effort (a finder; points checked on t10's model), torsion orders (elltors), and whether the surviving class is
the class of a rational point found.  Also the torsion order of E_LU at (14903,9696) (rank 0 by the pairing).
Written by an AI model (reviewer z15), not read by a person."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('MSS_REPO', os.path.abspath(os.path.join(HERE, '..', '..')))
sys.argv = [sys.argv[0]]; sys.path.insert(0, HERE)
import io, contextlib
import t10_pairct_copy as T
from fractions import Fraction as Fr
RS = T.RS; P = T.PARI
for (m, n), j, k in (((14903, 9696), 'L', 'U'), ((3320, 2303), 'L', 'U'), ((5623, 1017), '0', '8'), ((9113, 3200), 'L', 'U'), ((12313, 5913), '0', '8')):
    with contextlib.redirect_stdout(io.StringIO()): Pd = T.PairData(m, n, j, k)
    E = P.ellinit([0, Pd.Pr.A[0], 0, Pd.Pr.A[1], Pd.Pr.A[2]])
    print('(%d,%d) E_%s%s: torsion %s' % (m, n, j, k, P.elltors(E)[1]), flush=True)
for (m, n), j, k, t in (((9113, 3200), '0', '8', (1, 2, 511, 1759)), ((14903, 9696), '0', '8', (1, 2, 10414, 49198))):
    with contextlib.redirect_stdout(io.StringIO()): Pd = T.PairData(m, n, j, k)
    E = P.ellinit([0, Pd.Pr.A[0], 0, Pd.Pr.A[1], Pd.Pr.A[2]])
    print('(%d,%d) E_%s%s: torsion %s' % (m, n, j, k, P.elltors(E)[1]), flush=True)
    v = Pd.vec(Pd.cls_of_tuple(t))
    for eff in (3, 4, 5, 6):
        r = P.ellrank(E, eff)
        pts = [(Fr(str(Q[0])), Fr(str(Q[1]))) for Q in r[3]]
        for Q in pts: assert Pd.Pr.on(Q)
        imgv = Pd.imgv + [Pd.vec(Pd.kum(Q[0])) for Q in pts]
        s1 = len(RS.f2_span(imgv)); s2 = len(RS.f2_span(imgv + [v]))
        print('   effort %d: %s, %d point(s), images span %d; surviving class %s' % (eff, [int(r[0]), int(r[1])], len(pts), s1,
              'IS the class of a point found' if s1 == s2 else 'not in the span'), flush=True)
        if len(pts) >= 2: break
