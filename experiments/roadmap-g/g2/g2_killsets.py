"""G2: which primes empty the four-class set?  (EXPERIMENT, written by an AI model; not read by a person; nothing here
changes the census.)
For each ratio of a list (the 997 of tuple_selmer_997.txt, the 8,739 of class_selmer_15000_list.txt):
  * the bad primes B: every prime of 2mn(m^2-n^2)(m^2+n^2)(m^8+14m^4n^4+n^8), factorisations by PARI with each factor
    proved prime; the local test at every odd p in B is exact (exact_local.py: a node argument for large p, in the
    manner of Proposition B, so no prime of B is left out); at p = 2 the repository's test;
  * test 6 of COMPLETE.md section 7: does B alone empty the four-class set?
  * an irredundant killing set: the primes of B in increasing order, then good primes in increasing order (to 1500,
    then 5000), until the set is empty; then each prime, largest first, is dropped if the rest still empties it.
    Irredundant (no prime of it can be dropped), not proved of least size.
Output line: m,n  log2|space|  bad-alone: EMPTY | k left  kill: p:tag ...   (tag b = bad, g = good)
usage: g2_killsets.py LIST i N [--write]     (worker i of N; writes g2_<list>_<i>.txt)"""
import sys, os
from math import log2
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'g1'))
import exact_local as X
TD = X.TD
import cypari2
P = cypari2.Pari(); P.allocatemem(2 * 10**9)
def factorint_proved(x):
    x = abs(int(x))
    if x <= 1: return {}
    f = P.factor(x); out = {int(q): int(e) for q, e in zip(f[0], f[1])}
    for q in out: assert P.isprime(q), ('not proved prime', q)
    return out
TD.factorint = factorint_proved

def read_list(path):
    out = []
    for l in open(path):
        l = l.strip()
        if not l or l[0] == '#': continue
        t = l.split()[0].strip('():').split(':')[0].strip('()')
        m, n = map(int, t.split(',')); out.append((m, n))
    return out

def divisors(supports):
    divs = []
    for S in supports:
        ds = [1]
        for q in S: ds += [x*q for x in ds]
        divs.append(ds)
    return divs

def run(m, n):
    forms, sup, known = TD.magic(m, n)
    divs = divisors(sup); space = 1
    for d in divs: space *= len(d)
    bad = set()
    for x in (2*m*n*(m*m - n*n), m*m + n*n, m**8 + 14*m**4*n**4 + n**8): bad |= set(factorint_proved(x))
    bad = sorted(bad)
    cache = {}
    def loc(p):
        if p not in cache:
            L, al = X.allowed_at(forms, p)
            for kt in known: assert tuple(L.dlabel(d) for d in kt) in al, ('known point rejected', m, n, p)
            cache[p] = (p, L, al) if len(al) < len(L.labels())**len(forms) else None
        return cache[p]
    def S_of(ps): return X.join(divs, [c for c in (loc(p) for p in ps) if c is not None])
    left_bad = S_of(bad)
    order = list(bad); good = []
    if left_bad:
        for pmax in (1500, 5000):
            good = [int(q) for q in P.primes([3, pmax]) if int(q) not in bad]
            if not S_of(bad + good): break
        # shortest prefix of the good primes
        lo, hi = 0, len(good)
        if S_of(bad + good): return space, left_bad, None, known
        while lo < hi:
            mid = (lo + hi)//2
            if S_of(bad + good[:mid]): lo = mid + 1
            else: hi = mid
        order = bad + good[:lo]
    kill = [p for p in order if loc(p) is not None]
    for p in sorted(kill, reverse=True):
        rest = [q for q in kill if q != p]
        if not S_of(rest): kill = rest
    return space, left_bad, sorted(kill), known

if __name__ == '__main__':
    a = [x for x in sys.argv[1:] if not x.startswith('--')]
    path, i, N = a[0], int(a[1]), int(a[2])
    pairs = read_list(path)[i::N]
    name = os.path.basename(path).split('.')[0]
    out = []
    for m, n in pairs:
        space, lb, kill, known = run(m, n)
        bad = 'EMPTY' if not lb else '%d left' % len(lb)
        badset = set()
        for x in (2*m*n*(m*m - n*n), m*m + n*n, m**8 + 14*m**4*n**4 + n**8): badset |= set(factorint_proved(x))
        ks = 'NOT EMPTIED' if kill is None else ' '.join('%d:%s' % (p, 'b' if p in badset else 'g') for p in kill)
        out.append('%d,%d  %.0f  bad-alone: %s  kill: %s' % (m, n, log2(space), bad, ks))
        print(out[-1], flush=True)
    if '--write' in sys.argv: open(os.path.join(HERE, 'g2_%s_%d.txt' % (name, i)), 'w').write('\n'.join(out) + '\n')
