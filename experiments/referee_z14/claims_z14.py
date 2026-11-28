"""z14 task 5: distinct-case counts behind the numbers of EXPLORATION_3 §1 and EXPLORATION_4 §5 that come from the pilot
(experiments/conductor-t10/t10_*_[ab].txt, read by experiments/research-pairing/parse.py) and from explore3's logs.
A mixed curve is identified by its ratio up to (m, n) ~ (m+n, m-n) (made primitive): same curve up to 2^2 (task 3).
Written by an AI model (reviewer z14); not read by a person.   usage: python3 -I claims_z14.py | diff - claims_z14.txt"""
import os, sys, re, collections
from math import gcd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'research-pairing'))
os.chdir(os.path.join(HERE, '..', 'conductor-t10'))
from parse import load, sqf
R = load()
def partner(m, n):
    a, b = m + n, m - n; g = gcd(a, b); return (a // g, b // g)
def key(m, n): return (m, n) if (m - n) % 2 else partner(m, n)
say = print
say('# z14 claims check. Written by an AI model (reviewer z14); not read by a person.')
files = sorted(set(re.sub(r'.*/', '', f) for f in __import__('glob').glob('t10_*_[ab].txt')))
say('pilot files read by parse.load():', files)
rat = collections.Counter(c['ratio'] for c in R if c['cv'] == '08')
say('ratios: %d curve blocks of E_08 for %d distinct ratios (max repetitions %d); distinct mixed curves %d'
    % (sum(rat.values()), len(rat), max(rat.values()), len(set(key(*r) for r in rat))))
say('pair-curve blocks: %d (%s)' % (len(R), dict(collections.Counter(c['cv'] for c in R))))
# are the four mixed models identical per ratio?
same = collections.Counter()
by = collections.defaultdict(dict)
for c in R: by[c['ratio']][c['cv']] = c['e']
for r, d in by.items():
    mix = [d[k] for k in ('0L', '0U', '8L', '8U') if k in d]
    same[len(set(mix)) == 1 and len(mix) == 4] += 1
say('ratios whose four mixed blocks print the same model e:', dict(same))
# kill lines
kl = [(c, l) for c in R for l in c['lines'] if l['b']]
for cvs, nm in ((('0L', '0U', '8L', '8U'), 'mixed'), (('08',), '08'), (('LU',), 'LU')):
    L = [(c, l) for c, l in kl if c['cv'] in cvs]
    d1 = set((c['ratio'], c['cv'], l['tup']) for c, l in L)
    d2 = set((c['ratio'], c['cv'], l['cls'], l['b']) for c, l in L)
    d3 = set((c['ratio'], l['cls'], l['b']) for c, l in L)
    d4 = set((key(*c['ratio']), l['cls'], l['b']) for c, l in L)
    d5 = set((key(*c['ratio']), l['cls']) for c, l in L)
    say('%s kill lines: %d; distinct (ratio, curve, tuple) %d; distinct (ratio, curve, a, b) %d; distinct (ratio, a, b) %d; distinct (curve up to iso, a, b) %d; distinct (curve up to iso, a) %d'
        % (nm, len(L), len(d1), len(d2), len(d3), len(d4), len(d5)))
    if nm == 'mixed':
        say('  per curve name:', dict(collections.Counter(c['cv'] for c, l in L)))
        # first component of b (at e1)
        fc = collections.Counter(); fcd = collections.Counter(); seen = set()
        for c, l in L:
            b1 = l['b'][0]; cat = 'b1 = 1' if b1 == 1 else ('b1 = 2' if b1 == 2 else 'other')
            fc[cat] += 1
            k = (key(*c['ratio']), l['cls'], l['b'])
            if k not in seen: seen.add(k); fcd[cat] += 1
        say('  first component of b: lines %s; distinct (curve up to iso, a, b) %s' % (dict(fc), dict(fcd)))
        a1 = collections.Counter(l['cls'][0] for c, l in L)
        say('  first component of the killed class a over the lines:', dict(a1), '(so no kill line has a = (1, d, d))')
# EXPLORATION_4 §5 kappa: class = (2, kappa d_j d_k) with kappa = sqf(mn(m^2 - n^2)) * (1 for 0L, 8U; 2 for 0U, 8L)
idx = {'0': 0, '8': 1, 'L': 2, 'U': 3}
tab = collections.Counter(); tabd = collections.Counter(); seen = set()
for c in R:
    cv = c['cv']
    if cv in ('08', 'LU'): continue
    m, n = c['ratio']
    for l in c['lines']:
        D = l['tup'][idx[cv[0]]] * l['tup'][idx[cv[1]]]
        q = sqf(sqf(l['cls'][1] * D) * sqf(m * n * (m * m - n * n)))
        want = 1 if cv in ('0L', '8U') else 2
        ok = (l['cls'][0] == 2 and q == want)
        tab[(cv, ok)] += 1
        k = (c['ratio'], cv, l['tup'])
        if k not in seen: seen.add(k); tabd[(cv, ok)] += 1
say('kappa formula (class = (2, kappa D), kappa = sqf(mn(m^2-n^2)) x (1: 0L, 8U; 2: 0U, 8L)), all tuple lines of the mixed blocks:', dict(sorted(tab.items())))
say('  distinct (ratio, curve, tuple):', dict(sorted(tabd.items())))
blocks = collections.Counter()
for c in R:
    if c['cv'] in ('08', 'LU') or not c['lines']: continue
    m, n = c['ratio']
    blocks[(c['cv'], 'm even' if m % 2 == 0 else 'n even' if n % 2 == 0 else 'both odd')] += 1
say('  the constant kappa is one number per (ratio, curve): blocks with tuple lines', dict(sorted(blocks.items())),
    '; distinct ratios', len(set(c['ratio'] for c in R if c['cv'] == '0L' and c['lines'])))
# explore3 task23 samples: ratios and distinct curves
for f, nm in (('task23_seed11.log', 'seed 11'), ('task23_seed12_all201.log', 'seed 12')):
    rs = []
    for l in open(os.path.join('/home/claude/untrusted/explore3/explore3', f)):
        mm = re.match(r'ratio \((\d+), (\d+)\).*\| pairs (\d+)', l)
        if mm: rs.append((int(mm[1]), int(mm[2]), int(mm[3])))
    ks = collections.Counter(key(m, n) for m, n, p in rs)
    say('task23 %s: %d ratios, %d distinct mixed curves (%d curves sampled twice), Task-2 pairs %d' % (nm, len(rs), len(ks), sum(1 for v in ks.values() if v > 1), sum(p for m, n, p in rs)))
