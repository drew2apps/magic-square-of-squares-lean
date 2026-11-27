"""Review x6 (run 6), written by an AI model; not read by a person.  For every ratio of G2's lists that the bad primes alone do not empty, every tuple they leave gets explicit local
points (my own search and label code) at 2 and at every odd bad prime.  Prints only.  usage: python3 -I badalone_witness_all.py   (both lists of roadmap-g/g2; about 25 minutes)"""
import sys, os, random, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import mylocal as ML
X, TD, P = ML.X, ML.TD, ML.P
tot = [0, 0, 0, 0]; bad_rows = []
for f in [os.path.join(HERE, '..', 'roadmap-g', 'g2', x) for x in ('g2_997_all.txt', 'g2_15000_all.txt')]:
    for l in open(f):
        a = l.split()
        if len(a) < 4 or 'SKIPPED' in l or a[3] == 'EMPTY': continue
        m, n = map(int, a[0].split(',')); claimed = int(a[3])
        rnd = random.Random(m * 100003 + n); ML.WIT.clear()
        forms, sup, known = TD.magic(m, n); divs = ML.divisors(sup)
        bad = sorted(set(int(q) for x in (2*m*n*(m*m-n*n), m*m+n*n, m**8+14*m**4*n**4+n**8) for q in P.factor(x)[0]))
        locs = [(p,) + X.allowed_at(forms, p) for p in bad]
        T = X.join(divs, [x for x in locs if len(x[2]) < len(x[1].labels())**4])
        wit = [t for t in T if all(ML.witness(forms, (m, n), q, t, rnd) is not None for q in bad)]
        tot[0] += 1; tot[1] += len(T); tot[2] += len(wit); tot[3] += len(T) == claimed
        if len(wit) != len(T) or len(T) != claimed or not T: bad_rows.append((m, n, claimed, len(T), len(wit)))
        print('%d,%d  left %d (G2 says %d)  with explicit points at every bad prime: %d' % (m, n, len(T), claimed, len(wit)), flush=True)
print('# %d ratios; tuples left %d; with explicit points at 2 and at every odd bad prime %d; count as in G2 for %d; exceptions: %s' % (tot[0], tot[1], tot[2], tot[3], bad_rows or 'none'))
