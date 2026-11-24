import sys, itertools, time
import numpy as np
from mytuple import *

def primes_upto(n):
    return [q for q in range(2, n + 1) if all(q % r for r in range(2, int(q ** .5) + 1))]

def global_filter(supports, p, S):
    """boolean 4-d array over candidate tuples: local image in S"""
    cl = [cands(s) for s in supports]
    hc = [np.array([hclass(d, p) for d in c]) for c in cl]
    nr = len(reps(p))
    T = np.zeros((nr,) * 4, bool)
    for t in S:
        T[t] = True
    return T[np.ix_(hc[0], hc[1], hc[2], hc[3])]

def run(forms, supports, bound=61, verbose=True):
    ps = sorted(set(sum(supports, [])) | set(primes_upto(bound)))
    keepU = None; keepL = None; info = {}
    for p in ps:
        N = None if p <= 70 else 1
        U, L, N = local_sets(forms, p, N)
        Ls = set(L)
        assert Ls <= U, ("LOWER NOT IN UPPER", p, Ls - U)
        info[p] = (U, Ls, N, L)
        fU = global_filter(supports, p, U); fL = global_filter(supports, p, Ls)
        keepU = fU if keepU is None else keepU & fU
        keepL = fL if keepL is None else keepL & fL
        if verbose:
            print("  p=%4d N=%2d |H^4|=%3d |U_p|=%3d |L_p|=%3d %s  global after: U %d  L %d" % (
                p, N, len(reps(p)) ** 4, len(U), len(Ls), "exact" if U == Ls else "GAP", keepU.sum(), keepL.sum()))
    cl = [cands(s) for s in supports]
    left = [tuple(cl[j][i[j]] for j in range(4)) for i in np.argwhere(keepU)]
    return info, left, keepL.sum()

if __name__ == '__main__':
    kind = sys.argv[1]; a, b = int(sys.argv[2]), int(sys.argv[3])
    forms, sup = magic(a, b) if kind == 'magic' else ratdist(a, b)
    t0 = time.time()
    info, left, nl = run(forms, sup, int(sys.argv[4]) if len(sys.argv) > 4 else 61)
    print("left (upper):", left[:20], len(left), " lower-count:", nl, " time %.1f" % (time.time() - t0))
