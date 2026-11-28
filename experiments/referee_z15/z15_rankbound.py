"""z15: rank bounds on E_08 and E_LU at the five survivor ratios of EXPLORATION_4, with the repository's code:
Sel^2 by rank_second.full_selmer2 (through t10's PairData, the same integral model t10 uses), the Cassels-Tate
pairing matrix on a basis of Sel^2 by rank_second.DFull.pairing (t10's DFullX wrapper), F_2-rank of the matrix.
rank E(Q) <= dim Sel^2 - 2 - rank(pairing matrix).  PARI's ellrank (a finder; bounds quoted, not relied on).
Also: is the class of each surviving tuple in the span of the images of points PARI finds with ellrank effort 2?
Written by an AI model (reviewer z15), not read by a person.  usage: python3 z15_rankbound.py > z15_rankbound.txt"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('MSS_REPO', os.path.abspath(os.path.join(HERE, '..', '..')))
sys.argv = [sys.argv[0]]
sys.path.insert(0, HERE)
import io, contextlib
import t10_pairct_copy as T
from fractions import Fraction as Fr
RS = T.RS; RI = T.RI
SURV = {(3320, 2303): [(2, 1, 1, 635399), (2, 1, 1270798, 2)], (5623, 1017): [(2, 1270798, 1, 2), (635399, 1, 1, 2)],
        (9113, 3200): [(1, 2, 511, 1759), (1, 2, 3518, 1022)], (12313, 5913): [(1022, 3518, 2, 1), (1759, 511, 2, 1)],
        (14903, 9696): [(1, 2, 10414, 49198), (1, 2, 24599, 5207)]}
def f2rank(M):
    rows = [int(''.join(str(x) for x in r), 2) if r else 0 for r in M]
    return len(RS.f2_span([r for r in rows if r]))
print('z15_rankbound.py: written by an AI model (reviewer z15), not read by a person.')
pairs = [('0', '8'), ('L', 'U')]
if '--mixed' in sys.argv: pairs.append(('0', 'L'))
for (m, n), surv in SURV.items():
    for j, k in pairs:
        t0 = time.process_time()
        with contextlib.redirect_stdout(io.StringIO()):
            Pd = T.PairData(m, n, j, k)
        sel = Pd.sel; d = len(sel)
        M = []
        for a in sel:
            Da = Pd.dfull(a); row = []
            for b in sel:
                h, per = T.safe_pair(Da, b)
                assert h in (1, -1), ('undecided pairing', a, b, per)
                row.append(0 if h == 1 else 1)
            M.append(row)
        sym = all(M[i][jj] == M[jj][i] for i in range(d) for jj in range(d))
        alt = all(M[i][i] == 0 for i in range(d))
        rk = f2rank(M)
        E = T.PARI.ellinit([0, Pd.Pr.A[0], 0, Pd.Pr.A[1], Pd.Pr.A[2]])
        r = T.PARI.ellrank(E, 2)
        pts = [(Fr(str(P[0])), Fr(str(P[1]))) for P in r[3]]
        for P in pts: assert Pd.Pr.on(P)
        imgv = Pd.imgv + [Pd.vec(Pd.kum(P[0])) for P in pts]
        ri = len(RS.f2_span(imgv))
        w = int(T.PARI.ellrootno(E))
        print('(%d,%d) E_%s%s: dim Sel2 %d, pairing matrix rank %d (symmetric %s, alternating %s) => rank <= %d; '
              'PARI ellrank effort 2 [%d,%d]; root number %+d; images of torsion + points found span %d'
              % (m, n, j, k, d, rk, sym, alt, d - 2 - rk, int(r[0]), int(r[1]), w, ri))
        for t in surv:
            c = Pd.cls_of_tuple(t); v = Pd.vec(c)
            inimg = len(RS.f2_span(imgv + [v])) == ri
            # in the kernel of the pairing on all of Sel^2?
            Dc = Pd.dfull(c); ker = all(T.safe_pair(Dc, b)[0] == 1 for b in sel)
            print('    tuple %s: class %s; pairs trivially with all of Sel2: %s; class of a point found: %s' % (t, c, ker, inimg))
        sys.stdout.flush()
        sys.stderr.write('(%d,%d) %s%s %.1f s\n' % (m, n, j, k, time.process_time() - t0))
