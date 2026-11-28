"""Referee z17 (written by an AI model; not read by a person): checks of roadmap I's I2 residual file.
(a) the headline counts of i2/i2_residual.jsonl, recounted from the file itself;
(b) each killing prime is good: it divides none of 2mn(m^2-n^2), m^2+n^2, m^8+14m^4n^4+n^8 (exact, no factoring);
(c) a seeded sample (seed 17; 30 ratios from lines 1-2939, the H2 part, and 20 from lines 2940-3876, G2's 937, which
    I8's --compare50 never draws from) recomputed with conductor-t10's own residual() (exact_local of roadmap-g g1 and
    X.join; I2 used h4_split.chunked): same residual set, same killing prime per tuple.
usage: python3 -I z17_i2check.py [--write]     (z17_i2check.txt)"""
import sys, os, json, random, collections
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
I2 = os.path.join(REPO, 'experiments', 'roadmap-i', 'i2')
os.environ['MSS_REPO'] = REPO
argv = sys.argv; sys.argv = [argv[0]]
sys.path.insert(0, os.path.join(REPO, 'experiments', 'conductor-t10'))
import t10_pairct as T
sys.argv = argv
T.PARI.default('nbthreads', 1)
T.say = lambda s='': None
OUT = []
def say(s=''): OUT.append(s)

lst = []
for l in open(os.path.join(I2, 'i2_needgood_list.txt')):
    l = l.strip()
    if l and l[0] != '#':
        a, src = l.split()[:2]; m, n = map(int, a.split(',')); lst.append((m, n, src))
recs = [json.loads(l) for l in open(os.path.join(I2, 'i2_residual.jsonl'))]
say('z17: checks of roadmap-i/i2/i2_residual.jsonl. Written by an AI model (reviewer z17); not read by a person.')
say('(a) counts recomputed from the file')
say('  list lines %d; records %d; same ratios in the same order: %s' % (len(lst), len(recs), [(m, n) for m, n, _ in lst] == [(o['m'], o['n']) for o in recs]))
say('  sources in the list: %s' % dict(sorted(collections.Counter(s for _, _, s in lst).items())))
g2first = next(i for i, (_, _, s) in enumerate(lst) if s.startswith('g2'))
say('  first G2 line (1-based): %d; G2 lines: %d' % (g2first + 1, sum(1 for _, _, s in lst if s.startswith('g2'))))
keys = [sorted(o.keys()) for o in recs]
say('  every record has exactly the fields m, n, residual, kill, killing_prime: %s; field order as documented: %s'
    % (all(k == ['kill', 'killing_prime', 'm', 'n', 'residual'] for k in keys), all(list(o.keys()) == ['m', 'n', 'residual', 'kill', 'killing_prime'] for o in recs)))
sizes = [len(o['residual']) for o in recs]
say('  residual tuples: %d; |S| from %d to %d; every |S| even: %s' % (sum(sizes), min(sizes), max(sizes), all(s % 2 == 0 for s in sizes)))
say('  |S| most common: %s' % ', '.join('%d: %d' % kv for kv in collections.Counter(sizes).most_common(6)))
allk = [k for o in recs for k in o['kill']]
say('  kill entries null: %d; largest killing prime per tuple: %s' % (sum(k is None for k in allk), max(k for k in allk if k is not None)))
say('  killing_prime = max(kill) in every record: %s; residual sorted in every record: %s; no repeated tuple: %s'
    % (all(o['killing_prime'] == max(o['kill']) for o in recs), all(o['residual'] == sorted(o['residual']) for o in recs),
       all(len(set(map(tuple, o['residual']))) == len(o['residual']) for o in recs)))
say('  killing_prime most common: %s' % ', '.join('%d: %d' % kv for kv in collections.Counter(o['killing_prime'] for o in recs).most_common(6)))
bad = 0
for o in recs:
    m, n = o['m'], o['n']; N = 2 * m * n * (m * m - n * n) * (m * m + n * n) * (m ** 8 + 14 * m ** 4 * n ** 4 + n ** 8)
    for k in set(o['kill']):
        if k is None or N % k == 0 or not T.PARI.isprime(k): bad += 1
say('(b) killing primes that are bad, not prime, or null: %d' % bad)
say('(c) recomputed with t10_pairct.residual() (seed 17: 30 of lines 1-%d, 20 of lines %d-%d)' % (g2first, g2first + 1, len(lst)))
rnd = random.Random(17)
idx = sorted(rnd.sample(range(g2first), 30) + rnd.sample(range(g2first, len(lst)), 20))
agree = 0; ntup = 0
for i in idx:
    m, n, src = lst[i]; o = recs[i]
    left, kill, badp, known = T.residual(m, n)
    mine = {tuple(t): k for t, k in zip(o['residual'], o['kill'])}
    theirs = {t: kill.get(t) for t in left if t not in known}
    ok = mine == theirs and not known; agree += ok; ntup += len(mine)
    say('  line %d (%d,%d) %s: |S| = %d, P = %s, known tuples %d: %s' % (i + 1, m, n, src, len(mine), o['killing_prime'], len(known), 'same' if ok else 'DIFFER'))
say('  agree: %d of %d (%d tuples)' % (agree, len(idx), ntup))
s = '\n'.join(OUT) + '\n'; f = os.path.join(HERE, 'z17_i2check.txt')
if '--write' in argv: open(f, 'w').write(s)
else:
    sys.stdout.write(s)
    if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS')
