"""I8 checkpoint (EXPERIMENT, written by an AI model; not read by a person; nothing here changes the census; tier D).
Takes the records the I8 workers have written (i8_<i>of<N>.jsonl, or the .part while running) and puts the first K
ratios of i2_needgood_list.txt, in list order, into i8_first<K>.jsonl, with a summary i8_first<K>_summary.txt.
It refuses if any of the first K ratios is not done.  No PARI, no recomputation: it only selects and counts.
Per tuple: verdict 'pairing' (removed by the Cassels-Tate pairing on Sel^2 of the pair curves in killed_by),
'rank 0 of E_jk' (the ratio is excluded by a pair curve of rank 0), 'undecided', or 'neither'.
The four mixed curves are 0L, 0U, 8L, 8U; a tuple counted as 'mixed only' is pairing-killed on some of them and on
neither 08 nor LU.
usage: i8_checkpoint.py K N [--write]     N: the number of workers (3)"""
import sys, os, json, re
HERE = os.path.dirname(os.path.abspath(__file__))
I2 = os.environ.get('I2_DIR', os.path.join(HERE, '..', 'i2'))
MIXED = ('0L', '0U', '8L', '8U')

def read_list(path):
    return [tuple(map(int, re.findall(r'\d+', l)[:2])) for l in open(path) if l.strip() and not l.startswith('#')]

def records(N):
    out = {}
    for i in range(N):
        f = os.path.join(HERE, 'i8_%dof%d.jsonl' % (i, N))
        for p in (f, f + '.part'):
            if os.path.exists(p):
                for l in open(p):
                    if l.strip(): o = json.loads(l); out[(o['m'], o['n'])] = o
    return out

def summary(rows, K):
    out = ['I8 checkpoint: the Cassels-Tate pairing test of conductor-t10 on the first %d ratios of i2_needgood_list.txt '
           '(lines 1 to %d, list order). Written by an AI model; not read by a person. Experiment, not counted; tier D.' % (K, K)]
    c = {}; per = {}; mixed_only = 0; allk = 0; err = []; neither = []; r0 = []; src = {}
    for o in rows:
        if 'tuples' not in o: err.append((o['m'], o['n'], o.get('status'))); continue
        src[o['source']] = src.get(o['source'], 0) + 1
        allk += all(t['verdict'] == 'pairing' for t in o['tuples'])
        for t in o['tuples']:
            v = 'rank 0' if t['verdict'].startswith('rank 0') else t['verdict']
            c[v] = c.get(v, 0) + 1
            for jk in t['killed_by']: per[jk] = per.get(jk, 0) + 1
            if t['killed_by'] and set(t['killed_by']) <= set(MIXED): mixed_only += 1
            if v in ('neither', 'undecided'): neither.append((o, t, v))
        ex = sorted(jk for jk, d in (o.get('rank0') or {}).items() if d and d.get('excludes'))
        if ex: r0.append('(%d,%d) by %s' % (o['m'], o['n'], ', '.join('E_' + jk for jk in ex)))
    nr = len(rows) - len(err)
    out.append('ratios: %d done (residual tuples from %s), %d not done' % (nr, ', '.join('%s %d' % kv for kv in sorted(src.items())), len(err)))
    out.append('residual tuples: %d; pairing-killed %d; ratio excluded by a rank-0 pair curve %d; neither %d; undecided %d'
               % (sum(c.values()), c.get('pairing', 0), c.get('rank 0', 0), c.get('neither', 0), c.get('undecided', 0)))
    out.append('ratios all of whose tuples are pairing-killed: %d of %d' % (allk, nr))
    out.append('pairing kills per pair curve (a tuple counts once on each curve that kills it): %s' % ' '.join('%s:%d' % (jk, per.get(jk, 0)) for jk in ('08', '0L', '0U', '8L', '8U', 'LU')))
    out.append('tuples pairing-killed only on the mixed curves 0L, 0U, 8L, 8U (not on 08 or LU): %d' % mixed_only)
    for o, t, v in neither:
        out.append('  %s: (%d,%d) %s, killing good prime %s; per curve: %s' % (v.upper(), o['m'], o['n'], t['tau'], t['kill'],
                   ' '.join('%s:%s' % (jk, t['pairs'][jk]) for jk in ('08', '0L', '0U', '8L', '8U', 'LU'))))
    out.append('ratios that a pair curve of rank 0 excludes on its own (dim Sel^2 = 2, torsion complete, every point over t in {0, +-1, inf}): %d' % len(r0))
    out += ['  ' + s for s in r0]
    for e in err: out.append('  not done: (%d,%d) %s' % e)
    return out

if __name__ == '__main__':
    a = [x for x in sys.argv[1:] if not x.startswith('--')]
    K, N = int(a[0]), int(a[1])
    keys = read_list(os.path.join(I2, 'i2_needgood_list.txt'))[:K]
    recs = records(N)
    missing = [k for k in keys if k not in recs]
    if missing: sys.exit('not done yet: %d of the first %d ratios, the first (%d,%d)' % (len(missing), K, *missing[0]))
    rows = [recs[k] for k in keys]
    s1 = ''.join(json.dumps(o, separators=(', ', ': ')) + '\n' for o in rows)
    s2 = '\n'.join(summary(rows, K)) + '\n'
    f1 = os.path.join(HERE, 'i8_first%d.jsonl' % K); f2 = os.path.join(HERE, 'i8_first%d_summary.txt' % K)
    if '--write' in sys.argv:
        open(f1, 'w').write(s1); open(f2, 'w').write(s2)
    else:
        sys.stdout.write(s2)
        for f, s in ((f1, s1), (f2, s2)):
            if os.path.exists(f): print('# %s %s' % (os.path.basename(f), 'matches' if open(f).read() == s else 'DIFFERS'))
