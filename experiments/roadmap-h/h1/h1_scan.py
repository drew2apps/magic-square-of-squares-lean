"""H1, part 2: the four-class descent on random degenerate ratios up to m = 10^9 (EXPERIMENT, written by an AI model;
not read by a person; nothing here changes the census).

For each dyadic range k = 13 .. 29 (10^4 <= m < 2^(k+1), cut at 10^9), random coprime degenerate ratios (m^2 + n^2 a
square or twice a square), drawn from Pythagorean parameters by G1's deg_ratios with seed 3000 + k, COUNT per range
(default 1177, so 20,009 in all when every range has enough).  For each: G1's one() (the class criterion, then the
four-class descent with the exact local test at every odd prime, good primes to 150, then 1500, 5000, 20000 while
more than the known tuples is left).  A ratio whose descent leaves more than the tuples of its degenerate points is
listed with the extra tuples (a PHANTOM).
usage: h1_scan.py k [COUNT] [--write]    (h1_scan_k<k>.txt; without --write prints and compares)
       h1_scan.py --negative             a forged run: known tuples withheld on three ratios, which must be listed
G1's folder is found through G1_DIR, or else ../../roadmap-g/g1."""
import sys, os, time
HERE = os.path.dirname(os.path.abspath(__file__))
G1 = os.environ.get('G1_DIR', os.path.join(HERE, '..', '..', 'roadmap-g', 'g1'))
if not os.path.isdir(G1): G1 = os.path.join(HERE, '..', 'g1')
sys.path.insert(0, G1)
argv = sys.argv; sys.argv = [argv[0]]
import g1_stress as G
sys.argv = argv
WRITE = '--write' in argv

def scan(k, count):
    lo, hi = max(2**k, 10**4), min(2**(k + 1), 10**9 + 1)
    rs = G.deg_ratios(lo, hi, count, 3000 + k)
    lines = []; nph = 0; nnoclass = 0
    for m, n in rs:
        st, d = G.one(m, n)
        if d is None: lines.append('(%d,%d) FAILED: %s (a degenerate ratio must keep its points)' % (m, n, st)); continue
        S, known, pmax = d
        extra = sorted(set(S) - set(known))
        if set(known) - set(S): lines.append('(%d,%d) FAILED: a known tuple was removed' % (m, n))
        if extra:
            nph += 1
            lines.append('(%d,%d) PHANTOM: %s; good primes to %d; known %s; extra %s' % (m, n, st, pmax, known, extra))
    head = ['range k = %d: %d <= m < %d; %d degenerate ratios (seed %d)' % (k, lo, hi, len(rs), 3000 + k),
            'ratios with tuples beyond those of their degenerate points: %d' % nph]
    return head + lines

def negative():
    out = []
    for m, n in ((1298876, 247707), (1900421, 472260), (623461796, 439041597)):
        forms, sup, known = G.TD.magic(m, n)
        import exact_local as X
        primes = sorted(set(q for S_ in sup for q in S_) | set(int(q) for q in G.P.primes([2, 150])))
        S, used = X.descend_exact(forms, sup, [], primes)       # known tuples withheld from the comparison
        extra = sorted(set(S) - set())
        out.append('(%d,%d) with the known tuples withheld: %d tuples listed (%s)' % (m, n, len(extra), 'listed, as it must be' if extra and set(known) <= set(S) else 'NOT LISTED'))
    st, d = G.one(1109503, 597671); S, known, pmax = d; extra = sorted(set(S) - set(known))
    out.append('positive control (1109503,597671) (G1): %d extra tuples, good primes to %d (%s)' % (len(extra), pmax, 'listed, as it must be' if len(extra) == 4 else 'WRONG'))
    return out

if __name__ == '__main__':
    if '--negative' in argv:
        lines = negative(); f = os.path.join(HERE, 'h1_scan_negative.txt')
    else:
        k = int(argv[1]); count = int(argv[2]) if len(argv) > 2 and argv[2][0] != '-' else 1177
        t0 = time.time(); lines = scan(k, count); sys.stderr.write('k=%d %.0f s\n' % (k, time.time() - t0))
        f = os.path.join(HERE, 'h1_scan_k%d.txt' % k)
    s = '\n'.join(lines) + '\n'
    if WRITE: open(f, 'w').write(s)
    else:
        sys.stdout.write(s)
        if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS from ' + os.path.basename(f))
