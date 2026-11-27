"""H4: G2's kill sets for the three ratios of G1 that ran out of memory at 8 GB (EXPERIMENT, written by an AI model;
not read by a person; nothing here changes the census).
The same computation as G2's g2_killsets.run (bad primes exact, then good primes in increasing order, then the
irredundant kill set), with one change: the join of tuple_descent is run on the tuple space split by the class d0
(one chunk for each divisor d0 of the support of the first form), and the chunks are joined one at a time.  The union
of the chunks' results is the join of the whole space, so the sets are the same; only the memory differs.  An
emptiness test stops at the first chunk with a tuple left.
usage: h4_split.py [--write]   (h4_split.txt)      h4_split.py --negative   (h4_negative.txt)
G2's folder is found through G2_DIR, or else ../../roadmap-g/g2; G1's through G1_DIR, or else ../../roadmap-g/g1."""
import sys, os
from math import log2
HERE = os.path.dirname(os.path.abspath(__file__))
def find(env, name):
    d = os.environ.get(env, os.path.join(HERE, '..', '..', 'roadmap-g', name))
    return d if os.path.isdir(d) else os.path.join(HERE, '..', name)
sys.path.insert(0, find('G1_DIR', 'g1')); sys.path.insert(0, find('G2_DIR', 'g2'))
argv = sys.argv; sys.argv = [argv[0]]
import g2_killsets as K
sys.argv = argv
X = K.X; TD = K.TD; P = K.P
RATIOS = [(88373348, 70002375), (120756965, 14608109), (711343101, 552021980)]

def chunked(divs, loc, stop_first):
    out = []
    for d0 in divs[0]:
        S = X.join([[d0]] + divs[1:], loc)
        if S:
            out += S
            if stop_first: return out
    return sorted(out)

def run(m, n, chunks=True):
    forms, sup, known = TD.magic(m, n)
    divs = K.divisors(sup); space = 1
    for d in divs: space *= len(d)
    bad = set()
    for x in (2*m*n*(m*m - n*n), m*m + n*n, m**8 + 14*m**4*n**4 + n**8): bad |= set(K.factorint_proved(x))
    bad = sorted(bad); cache = {}
    def loc(p):
        if p not in cache:
            L, al = X.allowed_at(forms, p)
            for kt in known: assert tuple(L.dlabel(d) for d in kt) in al, ('known point rejected', m, n, p)
            cache[p] = (p, L, al) if len(al) < len(L.labels())**len(forms) else None
        return cache[p]
    def S_of(ps, first=True):
        lc = [c for c in (loc(p) for p in ps) if c is not None]
        return chunked(divs, lc, first) if chunks else X.join(divs, lc)
    left_bad = S_of(bad, first=False)
    order = list(bad); good = []
    if left_bad:
        for pmax in (1500, 5000):
            good = [int(q) for q in P.primes([3, pmax]) if int(q) not in bad]
            if not S_of(bad + good): break
        if S_of(bad + good): return space, left_bad, None, bad
        lo, hi = 0, len(good)
        while lo < hi:
            mid = (lo + hi)//2
            if S_of(bad + good[:mid]): lo = mid + 1
            else: hi = mid
        order = bad + good[:lo]
    kill = [p for p in order if loc(p) is not None]
    for p in sorted(kill, reverse=True):
        rest = [q for q in kill if q != p]
        if not S_of(rest): kill = rest
    return space, left_bad, sorted(kill), bad

def line(m, n, chunks=True):
    space, lb, kill, bad = run(m, n, chunks)
    b = 'EMPTY' if not lb else '%d left' % len(lb)
    ks = 'NOT EMPTIED' if kill is None else ' '.join('%d:%s' % (p, 'b' if p in bad else 'g') for p in kill)
    return '%d,%d  %.0f  bad-alone: %s  kill: %s' % (m, n, log2(space), b, ks)

if __name__ == '__main__':
    if '--negative' in argv:
        # the chunked join must agree with G2's unsplit join on ratios small enough for it (two need good primes), and on a degenerate ratio
        # (its known tuples must be left, so it must come out NOT EMPTIED or keep tuples)
        out = []
        for m, n in ((291, 224), (308, 185), (465, 416), (511, 260)):
            a, b = line(m, n, True), line(m, n, False)
            out.append('%s  (split: %s)' % (a, 'agrees with the unsplit join' if a == b else 'DIFFERS'))
        forms, sup, known = TD.magic(12, 5); divs = K.divisors(sup)
        lc = [c for c in ([(p,) + X.allowed_at(forms, p) for p in (2, 3, 5, 13)]) if len(c[2]) < len(c[1].labels())**4]
        S = chunked(divs, lc, False)
        out.append('(12,5) degenerate, primes 2, 3, 5, 13: split join keeps %d tuples, the known %s: %s' % (len(S), known, 'kept, as they must be' if set(known) <= set(S) else 'NOT KEPT'))
        f = 'h4_negative.txt'
    else:
        out = ['H4: G2 kill sets for the three ratios of G1 that ran out of memory at 8 GB (split by d0).',
               'Written by an AI model; not read by a person. Experiment, not counted.',
               'm,n  log2|space|  bad-alone: EMPTY | k left  kill: p:tag ...   (tag b = bad, g = good)']
        for m, n in RATIOS: out.append(line(m, n)); sys.stderr.write(out[-1] + '\n')
        f = 'h4_split.txt'
    s = '\n'.join(out) + '\n'; f = os.path.join(HERE, f)
    if '--write' in argv: open(f, 'w').write(s)
    else:
        sys.stdout.write(s)
        if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS from ' + os.path.basename(f))
