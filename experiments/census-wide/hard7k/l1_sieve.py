"""L(E,1) for root number +1 by the plain series L(E,1) = 2 * sum a_n/n * exp(-2 pi n/sqrt(N)),
a_p from PARI's ellap, a_n assembled by a numpy sieve.  For conductors where PARI's ellL1 runs out of
memory/time.  usage: l1_sieve.py <curve: m,n,name | gp vector> X.  EXPERIMENT."""
import sys, time, math
import numpy as np
import cypari2
pari = cypari2.Pari(); pari.allocatemem(2*10**9)
arg, X = sys.argv[1], int(float(sys.argv[2]))
if arg.startswith('['):
    E = pari.ellinit(pari(arg))
else:
    import importlib.util
    R = '/tmp/claude-0/-home-user/9ec7a8d3-42dc-53c8-8518-1a663ac809f1/scratchpad/repo/scripts/g1_full'
    spec = importlib.util.spec_from_file_location('g', R + '/global.py'); g = importlib.util.module_from_spec(spec); spec.loader.exec_module(g)
    m, n, name = arg.split(','); poly = g.curves(int(m), int(n))[name]
    c = [poly[4], poly[3], poly[2], poly[1]] if len(poly) == 5 else poly
    c0, c1, c2, c3 = c
    E = pari.ellinit([0, c2, 0, c1*c3, c0*c3*c3])
E = pari.ellminimalmodel(E)
N = int(pari.ellglobalred(E)[0]); w = int(pari.ellrootno(E))
print('N', N, 'root number', w, 'X', X, 'exp(-2piX/sqrtN) =', math.exp(-2*math.pi*X/math.sqrt(N)), flush=True)
assert w == 1
t0 = time.time()
# a_p for all primes <= X, in chunks
f = pari('(E,a,b)->my(v=List());forprime(p=a,b,listput(v,ellap(E,p)));Vec(v)')
P = np.array([], dtype=np.int64); A = np.array([], dtype=np.int64)
step = 10**7
for lo in range(2, X + 1, step):
    hi = min(lo + step - 1, X)
    ps = pari.primes([lo, hi]); aps = f(E, lo, hi)
    if len(ps):
        P = np.concatenate([P, np.array(str(ps)[1:-1].split(','), dtype=np.int64)])
        A = np.concatenate([A, np.array(str(aps)[1:-1].split(','), dtype=np.int64)])
    print('a_p to', hi, round(time.time() - t0), 's', flush=True)
bad = set(int(p) for p in pari.factor(N)[0])
# multiplicative sieve
a = np.ones(X + 1, dtype=np.int64); a[0] = 0
r = math.isqrt(X)
for p, ap in zip(P.tolist(), A.tolist()):
    if p <= r:
        L = X // p                      # multiples p*1..p*L
        tmp = np.full(L, ap, dtype=np.int64)
        pk, prev, cur, k = p * p, 1, ap, 1
        while pk <= X:
            nxt = cur * ap if p in bad else ap * cur - p * prev
            prev, cur = cur, nxt
            j = pk // p                 # index j-1 in tmp corresponds to n = p*j
            tmp[j - 1::j] = cur
            pk *= p
        a[p::p] *= tmp
    else:
        if ap != 1:
            a[p::p] *= ap
print('sieve done', round(time.time() - t0), 's', flush=True)
s = math.sqrt(N); cst = 2 * math.pi / s
total = 0.0; marks = {}
ch = 10**7
for lo in range(1, X + 1, ch):
    hi = min(lo + ch, X + 1)
    nn = np.arange(lo, hi, dtype=np.float64)
    total += float(np.sum(a[lo:hi] / nn * np.exp(-cst * nn)))
    marks[hi - 1] = 2 * total
print('partial sums 2*sum_{n<=x}:', {k: round(v, 12) for k, v in list(marks.items())[-6:]})
print('L(E,1) ~', 2 * total, 'seconds', round(time.time() - t0), flush=True)
