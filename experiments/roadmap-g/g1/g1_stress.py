"""G1: a stress test of the reduction (R) on random ratios far beyond m <= 15000 (EXPERIMENT, written by an AI model;
not read by a person; nothing here changes the census).

For each dyadic range [2^k, 2^(k+1)) cut to 10^4 <= m <= 10^9 (k = 13 .. 29), random coprime ratios m > n >= 1 with m
uniform in the range and n uniform below m (seed: 1000 + k), until COUNT ratios are drawn (default 20000).  For each:
  1. the repository's class criterion: class_selmer.empty (no local point somewhere) and class_selmer_linalg.system /
     solve (the linear system of Theorems C and C' over F_2), unchanged, with the factorisations it needs made by PARI
     and each prime factor proved prime by PARI's isprime (asserted);
  2. for every ratio that keeps a class, the four-class descent of tuple_descent.py (LAST_PAIR.md, section 12) with the
     local test made exact at every odd prime (exact_local.py; at p = 2 the repository's test), at every prime of the
     supports and at the good primes up to 150, then up to 1500, 5000 and 20000 while something other than the known tuples is left.
Controls in each range: degenerate ratios (square and twice-square type, from Pythagorean parameters, seed 2000 + k)
must keep exactly the tuples of their degenerate points.
A non-degenerate ratio with a tuple left is a SURVIVOR and is written first.
usage: g1_stress.py k [COUNT] [CONTROLS] [--write]   (writes g1_k<k>.txt; without --write prints and compares)"""
import sys, os, random, json, time
from math import gcd, isqrt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import exact_local as X
TD = X.TD
sys.path.insert(0, os.path.join(X.REPO, 'scripts', 'cover_descent', 'structure'))
import class_selmer as C
import class_selmer_linalg as CL
import cypari2
P = cypari2.Pari(); P.allocatemem(2 * 10**9)

def fac_proved(x):
    x = abs(int(x)); f = P.factor(x); ps = [int(q) for q in f[0]]
    for q in ps: assert P.isprime(q), ('not proved prime', q)
    return sorted(ps)
C.fac = fac_proved                                    # the criterion's factorisations, with every factor proved prime
def factorint_proved(x):
    x = abs(int(x))
    if x <= 1: return {}
    f = P.factor(x); out = {int(q): int(e) for q, e in zip(f[0], f[1])}
    for q in out: assert P.isprime(q), ('not proved prime', q)
    return out
TD.factorint = factorint_proved                       # tuple_descent's factorisations (supports, squarefree parts), proved

def is_deg(m, n): N = m*m + n*n; return isqrt(N)**2 == N or isqrt(2*N)**2 == 2*N

def descent(m, n):
    forms, sup, known = TD.magic(m, n)
    for pmax in (150, 1500, 5000, 20000):
        primes = sorted(set(q for S in sup for q in S) | set(int(q) for q in P.primes([2, pmax])))
        S, used = X.descend_exact(forms, sup, known, primes)
        if S == known or not S: break
    return S, known, pmax

def one(m, n):
    sy = CL.system(m, n)
    if sy is None: return 'no point', None
    ok, dim = CL.solve(*sy)
    if not ok: return 'no class', None
    S, known, pmax = descent(m, n)
    return 'class dim %d' % dim, (S, known, pmax)

def deg_ratios(lo, hi, count, seed):
    rnd = random.Random(seed); out = set(); tries = 0
    amax = isqrt(2*hi) + 2
    while len(out) < count and tries < 10**6:
        tries += 1
        a = rnd.randrange(2, amax); b = rnd.randrange(1, a)
        if (a - b) % 2 == 0 or gcd(a, b) != 1: continue
        x, y = a*a - b*b, 2*a*b
        for m, n in ((max(x, y), min(x, y)), (max(x + y, abs(x - y)), min(x + y, abs(x - y)))):
            if lo <= m < hi and gcd(m, n) == 1: out.add((m, n))
    return sorted(out)[:count]

def run(k, count, controls):
    lo, hi = max(2**k, 10**4), min(2**(k + 1), 10**9 + 1)
    rnd = random.Random(1000 + k); seen = set(); stats = {}; surv = []; kept = []
    t0 = time.time()
    while len(seen) < count:
        m = rnd.randrange(lo, hi); n = rnd.randrange(1, m)
        if gcd(m, n) != 1 or (m, n) in seen: continue
        seen.add((m, n))
        if is_deg(m, n): stats['degenerate drawn'] = stats.get('degenerate drawn', 0) + 1; continue
        st, d = one(m, n)
        key = st.split(' dim')[0]; stats[key] = stats.get(key, 0) + 1
        if d is not None:
            S, known, pmax = d
            kept.append('(%d,%d) %s; descent (good primes to %d): %s' % (m, n, st, pmax, 'nothing is left' if not S else 'LEFT %s' % S))
            if S: surv.append((m, n, S))
    t1 = time.time()
    ctl = []
    for m, n in deg_ratios(lo, hi, controls, 2000 + k):
        st, d = one(m, n)
        if d is None: ctl.append('(%d,%d) degenerate: %s; FAILED (must keep its points)' % (m, n, st)); continue
        S, known, pmax = d
        ctl.append('(%d,%d) degenerate: %s; descent leaves %s; %s' % (m, n, st, S, 'exactly the degenerate points' if S == known else 'NOT the degenerate points'))
    t2 = time.time()
    pe = stats.get('no class', 0) + stats.get('class', 0)
    lines = ['range k = %d: %d <= m < %d; %d random coprime ratios (seed %d)' % (k, lo, hi, count, 1000 + k)]
    lines.append('SURVIVORS (non-degenerate ratios with a tuple left by the four-class descent): %s' % (surv if surv else 'none'))
    lines.append('degenerate among the random ratios: %d' % stats.get('degenerate drawn', 0))
    lines.append('no local point somewhere: %d; a local point everywhere: %d' % (stats.get('no point', 0), pe))
    lines.append('of those with points everywhere, a class kept by the criterion: %d (%s)' % (stats.get('class', 0), '%.2f%%' % (100.0*stats.get('class', 0)/pe) if pe else '-'))
    lines += kept
    lines.append('controls (degenerate ratios of the range, seed %d): %d, %d keep exactly their degenerate points' % (2000 + k, len(ctl), sum('exactly the degenerate points' in c for c in ctl)))
    lines += ctl
    sys.stderr.write('k=%d random %.0f s, controls %.0f s\n' % (k, t1 - t0, t2 - t1))
    return lines

def negative():
    """The survivor test is not vacuous: degenerate ratios, passed through the same code with the degeneracy check
    switched off, must each be reported as SURVIVORS (their degenerate points leave tuples), and a degenerate ratio whose
    known points are removed from the descent's list must still leave those tuples."""
    lines = []
    for k in (14, 20, 26):
        lo, hi = max(2**k, 10**4), min(2**(k + 1), 10**9 + 1)
        for m, n in deg_ratios(lo, hi, 3, 3000 + k):
            st, d = one(m, n)
            S = d[0] if d else []
            lines.append('(%d,%d) degenerate, check off: %s; %s' % (m, n, st, 'SURVIVOR (as it must be)' if S else 'NOT REPORTED (failure)'))
    for m, n in deg_ratios(10**4, 2**14, 3, 3999):
        forms, sup, known = TD.magic(m, n)
        primes = sorted(set(q for S in sup for q in S) | set(int(q) for q in P.primes([2, 150])))
        S, used = X.descend_exact(forms, sup, [], primes)
        lines.append('(%d,%d) known points withheld: descent leaves %d tuples; %s' % (m, n, len(S), 'they include the known ones' if set(known) <= set(S) and S else 'FAILURE'))
    lines.append('# %d of %d refused as they must be' % (sum(('as it must be' in l) or ('include the known' in l) for l in lines), len(lines)))
    return lines

if __name__ == '__main__' and '--negative' in sys.argv:
    out = '\n'.join(negative()) + '\n'; f = os.path.join(HERE, 'g1_negative.txt')
    if '--write' in sys.argv: open(f, 'w').write(out)
    print(out, end='')
    if '--write' not in sys.argv and os.path.exists(f): assert open(f).read() == out; print('# matches g1_negative.txt')
    sys.exit()
if __name__ == '__main__':
    a = [x for x in sys.argv[1:] if not x.startswith('--')]
    k = int(a[0]); count = int(a[1]) if len(a) > 1 else 20000; controls = int(a[2]) if len(a) > 2 else 40
    out = '\n'.join(run(k, count, controls)) + '\n'; f = os.path.join(HERE, 'g1_k%d.txt' % k)
    if '--write' in sys.argv: open(f, 'w').write(out)
    print(out, end='')
    if '--write' not in sys.argv and os.path.exists(f): assert open(f).read() == out, 'differs from %s' % f; print('# matches %s' % os.path.basename(f))
