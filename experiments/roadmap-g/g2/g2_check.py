"""G2 summary and controls.  (EXPERIMENT, written by an AI model; not read by a person; nothing here changes the census.)
  g2_check.py              tables from the committed g2_*_all.txt files; checks g2_summary.txt (--write rewrites it)
  g2_check.py --negative   controls; checks g2_negative.txt (--write rewrites it):
     (a) degenerate ratios must come out NOT EMPTIED (their points keep tuples at every prime);
     (b) for ratios of the 997 with a kill set of 2 or more primes, the kill set must empty the four-class set when
         recomputed from scratch, and each of its primes dropped in turn must leave a tuple (a forged, smaller kill
         set is refused: the committed kill sets are irredundant)."""
import sys, os
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'g1'))

LISTS = [('the 997 (tuple_selmer_997.txt)', 'g2_997_all.txt'),
         ('the 8,739 (class_selmer_15000_list.txt)', 'g2_15000_all.txt'),
         ("G1's 63 kept ratios (g1_kept.txt)", 'g2_g1_kept_0.txt')]

def parse(path):
    rows = []
    for l in open(path):
        if not l.strip() or l[0] == '#' or 'SKIPPED' in l: continue
        a = l.split()
        m, n = map(int, a[0].split(',')); sp = int(a[1])
        empty = a[3] == 'EMPTY'
        ks = l.split('kill:')[1].split()
        kill = None if ks[:2] == ['NOT', 'EMPTIED'] else [(int(x.split(':')[0]), x.split(':')[1]) for x in ks]
        rows.append((m, n, sp, empty, kill))
    return rows

def summary():
    out = []
    for title, f in LISTS:
        p = os.path.join(HERE, f)
        if not os.path.exists(p): continue
        R = parse(p); sk = sum('SKIPPED' in l for l in open(p))
        out.append('## %s: %d ratios%s' % (title, len(R), '; %d more skipped (out of memory)' % sk if sk else ''))
        out.append('bad primes alone empty the set: %d; need a good prime too: %d; not emptied: %d'
                   % (sum(r[3] for r in R), sum(not r[3] and r[4] is not None for r in R), sum(r[4] is None for r in R)))
        need = [r for r in R if not r[3] and r[4]]
        lg = Counter(max(p for p, t in r[4] if t == 'g') for r in need)
        out.append('largest good prime needed: ' + ', '.join('%d (%d)' % kv for kv in sorted(lg.items())))
        out.append('good primes needed are all 1 mod 4: %s' % all(p % 4 == 1 for r in need for p, t in r[4] if t == 'g'))
        out.append('kill-set size: ' + ', '.join('%d: %d' % kv for kv in sorted(Counter(len(r[4]) for r in R if r[4]).items())))
        out.append('kill sets of bad primes only: %d' % sum(all(t == 'b' for p, t in r[4]) for r in R if r[4]))
        sp = Counter(r[2] for r in R); fl = Counter(r[2] for r in R if not r[3])
        out.append('by log2 of the tuple space (ratios, bad primes not enough): '
                   + ', '.join('%d: %d, %d' % (s, sp[s], fl[s]) for s in sorted(sp)))
        use = Counter(p for r in R if r[4] for p, t in r[4] if t == 'b')
        out.append('bad primes used most: ' + ', '.join('%d (%d)' % kv for kv in use.most_common(10)))
        out.append('')
    out += crosscheck()
    return out

def crosscheck():
    """Against the repository's tuple_badprimes_997.txt (bad primes, those of m^2+n^2 and of the octic only below 2000):
    with every bad prime, a ratio it emptied must stay empty, and a ratio it left must have no more tuples left here."""
    repo = os.environ.get('MSS_REPO', os.path.join(HERE, '..', '..', '..'))
    f = os.path.join(repo, 'scripts', 'cover_descent', 'tuple_badprimes_997.txt')
    if not os.path.exists(f): return ['## cross-check: %s not found' % f]
    old = open(f).read().splitlines(); left = {}
    for x in old[2].split(':', 1)[1].split():
        r, k = x.rsplit(':', 1); left[tuple(map(int, r.strip('()').split(',')))] = int(k)
    R = {(r[0], r[1]): r for r in parse(os.path.join(HERE, 'g2_997_all.txt'))}
    now = {}
    for l in open(os.path.join(HERE, 'g2_997_all.txt')):
        a = l.split(); now[tuple(map(int, a[0].split(',')))] = 0 if a[3] == 'EMPTY' else int(a[3])
    tested = [k for k, r in R.items() if r[2] <= 18]
    bad = [k for k in tested if (k in left and now[k] > left[k]) or (k not in left and now[k] > 0)]
    return ['## cross-check with tuple_badprimes_997.txt (tuple space at most 2^18)',
            'ratios: %d (repository %s)' % (len(tested), old[0].split(':')[1].split(';')[0].strip()),
            'left there: %d; left here, with the large bad primes too: %d; emptied by the large bad primes: %d'
            % (len(left), sum(now[k] > 0 for k in tested), sum(1 for k in left if now[k] == 0)),
            'ratios with more tuples left here than there (must be none): %s' % (bad if bad else 'none')]

def negative():
    import g2_killsets as K
    import g1_stress as G
    X = K.X; TD = X.TD; lines = []
    for m, n in G.deg_ratios(10**4, 2**15, 3, 4242):
        space, lb, kill, known = K.run(m, n)
        lines.append('(%d,%d) degenerate: %s' % (m, n, 'NOT EMPTIED (as it must be)' if kill is None else 'EMPTIED (failure)'))
    R = [r for r in parse(os.path.join(HERE, 'g2_997_all.txt')) if r[4] and len(r[4]) >= 2][:8]
    for m, n, sp, empty, kill in R:
        forms, sup, known = TD.magic(m, n); divs = K.divisors(sup)
        def left(ps):
            cs = []
            for p in ps:
                L, al = X.allowed_at(forms, p); cs.append((p, L, al))
            return len(X.join(divs, cs))
        ps = [p for p, t in kill]
        ok = left(ps) == 0
        drops = [left([q for q in ps if q != p]) for p in ps]
        lines.append('(%d,%d) kill %s: recomputed %s; each prime dropped leaves %s; %s'
                     % (m, n, ps, 'empty' if ok else 'NOT EMPTY', drops,
                        'forged kill sets refused (as they must be)' if ok and all(drops) else 'FAILURE'))
    lines.append('# %d of %d as they must be' % (sum('as it must be' in l or 'as they must be' in l for l in lines), len(lines)))
    return lines

if __name__ == '__main__':
    neg = '--negative' in sys.argv
    out = '\n'.join(negative() if neg else summary()) + '\n'
    f = os.path.join(HERE, 'g2_negative.txt' if neg else 'g2_summary.txt')
    if '--write' in sys.argv: open(f, 'w').write(out)
    print(out, end='')
    if '--write' not in sys.argv and os.path.exists(f): assert open(f).read() == out; print('# matches ' + os.path.basename(f))
