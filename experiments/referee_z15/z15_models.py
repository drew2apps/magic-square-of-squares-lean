"""z15: are the session's closed-form models y^2 = x(x - A^2)(x - C^2) (E_08), x(x - B^2)(x - C^2) (E_LU),
x(x - 2(C^2 - AB))(x - 2(C^2 + AB)) (E_mix) isomorphic over Q to t10's E_08, E_LU, E_0L, E_0U, E_8L, E_8U (Washington's
map from the quartic)?  Test: equal minimal models (PARI ellminimalmodel).  Ratios: the five survivor ratios and 15
drawn with seed 5 from the 937.  Written by an AI model (reviewer z15), not read by a person."""
import os, sys, json, random
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('MSS_REPO', os.path.abspath(os.path.join(HERE, '..', '..')))
sys.argv = [sys.argv[0]]; sys.path.insert(0, HERE)
import io, contextlib
import t10_pairct_copy as T
P = T.PARI
def minmodel(a):
    return [int(x) for x in P.ellminimalmodel(P.ellinit(a))[:5]]
def cubic(e): return [0, -(e[0] + e[1] + e[2]), 0, e[0] * e[1] + e[0] * e[2] + e[1] * e[2], -e[0] * e[1] * e[2]]
allr = [tuple(json.loads(l)[k] for k in ('m', 'n')) for l in open(os.path.join(HERE, 'z15_residual_all.jsonl'))]
R = [(3320, 2303), (5623, 1017), (9113, 3200), (12313, 5913), (14903, 9696)] + sorted(random.Random(5).sample(allr, 15))
ok = 0; tot = 0
print('z15_models.py: written by an AI model (reviewer z15), not read by a person.')
for m, n in R:
    A, B, C = m * m - n * n, 2 * m * n, m * m + n * n
    ses = {'08': cubic((0, A * A, C * C)), 'LU': cubic((0, B * B, C * C)), 'mix': cubic((0, 2 * (C * C - A * B), 2 * (C * C + A * B)))}
    res = []
    for j, k in T.PAIRS:
        with contextlib.redirect_stdout(io.StringIO()): Pr = T.Pair(m, n, j, k)
        mine = minmodel([0, Pr.A[0], 0, Pr.A[1], Pr.A[2]])
        key = j + k if j + k in ses else 'mix'
        same = mine == minmodel(ses[key]); ok += same; tot += 1
        res.append('%s%s:%s' % (j, k, 'iso' if same else 'NOT ISO'))
    print((m, n), ' '.join(res))
print('isomorphic: %d of %d' % (ok, tot))
