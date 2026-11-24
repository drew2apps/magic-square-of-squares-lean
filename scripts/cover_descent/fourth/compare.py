import sys, os, time, random, importlib.util
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mytuple import *
from drive import primes_upto, global_filter
spec = importlib.util.spec_from_file_location('repo_td', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tuple_descent.py'))
TD = importlib.util.module_from_spec(spec); spec.loader.exec_module(TD)

def repoK(p): return 12 if p == 2 else (6 if p < 50 else 3)

def repo_allowed(forms, p, K=None):
    L = TD.Local(forms, p, K or repoK(p)); al = L.allowed()
    lab = [L.dlabel(r) for r in reps(p)]
    assert len(set(lab)) == len(lab)
    inv = {l: k for k, l in enumerate(lab)}
    out = set()
    for t in al:
        assert all(l in inv for l in t), ("repo label that is not a rational class", p, t)
        out.add(tuple(inv[l] for l in t))
    return out

STAT = dict(cases=0, L_in_R=0, R_in_U=0, U_in_R=0, all3=0, exactLU=0)
BUGS = []; WEAK = []; LOOSE = []

def compare(name, forms, supports=None, bound=61, primes=None, verbose=False, extra=()):
    ps = primes or sorted(set(sum(supports, [])) | set(primes_upto(bound)))
    keepU = keepR = keepL = None
    for p in ps:
        U, L, N, how = local_best(forms, p, extra=extra)
        Ls = set(L); R = repo_allowed(forms, p)
        STAT['cases'] += 1
        a = Ls <= R; b = R <= U; c = U <= R
        STAT['L_in_R'] += a; STAT['R_in_U'] += b; STAT['U_in_R'] += c; STAT['all3'] += (Ls == R == U); STAT['exactLU'] += (Ls == U)
        if not a:
            t = sorted(Ls - R)[0]
            BUGS.append((name, p, t, L[t], [reps(p)[k] for k in t])); print("  !!! BUG: certified tuple rejected by repo", name, p, t, L[t])
        if not b:
            WEAK.append((name, p, len(R - U), len(R), len(U), N))
        if not c:
            LOOSE.append((name, p, 'N=%d' % N, 'exhaustion keeps %d, descent allows %d, with witness %d' % (len(U), len(R), len(Ls))))
        if verbose or not a:
            print("  %s p=%d N=%d(%s) |L|=%d |R|=%d |U|=%d  L<=R:%s R<=U:%s U<=R:%s" % (name, p, N, how, len(Ls), len(R), len(U), a, b, c))
        if supports:
            fU = global_filter(supports, p, U); fR = global_filter(supports, p, R); fL = global_filter(supports, p, Ls)
            keepU = fU if keepU is None else keepU & fU; keepR = fR if keepR is None else keepR & fR
            keepL = fL if keepL is None else keepL & fL
    if supports:
        cl = [cands(s) for s in supports]
        f = lambda K: sorted(tuple(cl[j][i[j]] for j in range(4)) for i in np.argwhere(K))
        return f(keepL), f(keepR), f(keepU)

if __name__ == '__main__':
    t0 = time.time()
    mag = [(627,364),(820,799),(889,320),(15,8),(7,1),(4,3),(12,5)]
    more = [l.split()[0] for l in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'gauss_descent', 'pairs157.txt')) if ',' in l.split()[0:1][0] if l.strip()] if False else []
    words = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'gauss_descent', 'pairs157.txt')).read().replace('\n', ' ').split()
    prs = [tuple(int(x) for x in w.strip('()').split(',')) for w in words if ',' in w and w.strip('()').replace(',', '').isdigit()]
    random.seed(11); mag += random.sample([q for q in prs if q not in mag], 10)
    for (m, n) in mag:
        forms, sup = magic(m, n)
        fr, sr, kn = TD.magic(m, n); assert [tuple(map(tuple, f)) for f in fr] == [tuple(map(tuple, f)) for f in forms], "forms differ"
        l, r, u = compare('magic(%d,%d)' % (m, n), forms, sup)
        S, used, sizes = TD.descend(fr, sr, kn, 61)
        print("magic (%d,%d): final  L %d  repo(local sets) %d  U %d ; repo descend(61): %d %s ; U==repo: %s ; known %s  [%.0fs]" % (
            m, n, len(l), len(r), len(u), len(S), S[:3], u == [tuple(s) for s in S], kn, time.time() - t0), flush=True)
    rd = [(1,3),(2,5),(3,7),(1,8),(7,16),(26,225),(7,30),(17,40),(23,48),(2,9),(3,10),(4,11),(7,15),(7,17),(1,18),(6,23),(14,31),(7,24),(31,72)]
    for (r_, s_) in rd:
        forms, sup = ratdist(r_, s_)
        l, r, u = compare('rd(%d/%d)' % (r_, s_), forms, sup)
        S, used, sizes = TD.descend(forms, sup, [], 61)
        const = lambda X: sum(1 for t in X if len(set(t)) == 1)
        print("rd %d/%d: final  L %d  repo %d  U %d (non-constant in U: %d) ; repo descend(61) on my forms: %d ; U==repo: %s  [%.0fs]" % (
            r_, s_, len(l), len(r), len(u), len(u) - const(u), len(S), u == [tuple(s) for s in S], time.time() - t0), flush=True)
    for (m, n) in [(15,8),(3,2),(2,1)]:
        f3, s3, kn = TD.cuboid(m, n)
        forms = [tuple(map(tuple, f)) for f in f3] + [tuple(map(tuple, f3[2]))]
        sup = s3 + [s3[2]]
        l, r, u = compare('cuboid(%d,%d)' % (m, n), forms, sup)
        S, used, sizes = TD.descend(f3, s3, kn, 61)
        pj = lambda X: sorted(set(t[:3] for t in X if t[2] == t[3]))
        print("cuboid (%d,%d): final  L %d  repo %d  U %d ; repo descend(61), 3 forms: %d ; U==repo: %s  [%.0fs]" % (
            m, n, len(pj(l)), len(pj(r)), len(pj(u)), len(S), pj(u) == [tuple(s) for s in S], time.time() - t0), flush=True)
    print("STAT", STAT); print("BUGS", BUGS); print("R not inside U (repo weaker than exhaustion):", WEAK); print("U not inside R (the exhaustion keeps a tuple that the descent removes):", LOOSE)
