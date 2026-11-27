"""Review z9 of roadmap H (run 7), written by an AI model; not read by a person.  Experiment, not counted.
Checks G3's table for 11000 < m <= 14000 (experiments/roadmap-g/g3, block 13000 added in run 7):
  * the family by brute force over every pair 11000 < m <= 14000, 0 < n < m (not via Pythagorean parameters): coprime,
    m^2 + n^2 a square or twice a square, no prime = 5 mod 8 to an odd power in m^2 + n^2, not in one of the 484
    classes of scripts/g1_full/modp_classes.txt (p not dividing n) -- the conditions of census2000/gen_pairs.py, read
    from its prose -- against fam_11000_30000.json and the (m, n, type) of the rows of CENSUS_TABLE_G3.tsv;
  * the counts by step and type, and block 13000 < m <= 14000 against stage_13000.json;
  * for the six OPEN ratios and their partners, every record of sd_*, dm_*, r0_* (nothing closes them) and the row of
    the 1000 < m <= 11000 table (census-wide/table) where the partner lies there;
  * the route row (13103,4271): its partner, of square type, and that partner's row in the 1000 < m <= 11000 table.
Rests on: the committed records (it reruns no curve computation).  usage: python3 -I g3_family_z9.py [--write]
(g3_family_z9.txt; about 30 s)"""
import os, sys, json, glob, math, collections
from math import isqrt, gcd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import *
O = Out('g3_family_z9.txt')
O('Review z9: G3 family and table, 11000 < m <= 14000. Written by an AI model; not read by a person.')
cls = collections.defaultdict(set)
for line in open(os.path.join(REPO, 'scripts', 'g1_full', 'modp_classes.txt')):
    p, r = line.split()[:2]; cls[int(p)].add(int(r))
def odd5(N):
    x = N; d = 2
    while d * d <= x:
        if x % d == 0:
            k = 0
            while x % d == 0: x //= d; k += 1
            if d % 8 == 5 and k % 2: return True
        d += 1 if d == 2 else 2
    return x > 1 and x % 8 == 5
def kind(m, n): s = m * m + n * n; return 'square' if isqrt(s) ** 2 == s else 'twice'
def partner(m, n):
    a, b = m + n, abs(m - n)
    if a % 2 == 0 and b % 2 == 0: a, b = a // 2, b // 2
    g = gcd(a, b); return (max(a, b) // g, min(a, b) // g)
mine = []; cnt = collections.Counter()
for m in range(11001, 14001):
    for n in range(1, m):
        s = m * m + n * n; r = isqrt(s)
        if not (r * r == s or (s % 2 == 0 and isqrt(s // 2) ** 2 == s // 2)) or gcd(m, n) != 1: continue
        cnt['degenerate coprime'] += 1
        if odd5(s): cnt['removed: 5 mod 8'] += 1; continue
        if any(n % p and (m * pow(n, -1, p)) % p in rs for p, rs in cls.items()): cnt['removed: a class'] += 1; continue
        mine.append((m, n, kind(m, n)))
fam = sorted(tuple(p) for p in json.load(open(os.path.join(G3, 'fam_11000_30000.json')))['open'] if 11000 < p[0] <= 14000)
tab = [l.rstrip('\n').split('\t') for l in open(os.path.join(G3, 'CENSUS_TABLE_G3.tsv'))][1:]
O('brute force: %s; kept %d (%s); fam_11000_30000.json in range: %d, identical: %s' % (
    dict(cnt), len(mine), dict(collections.Counter(t for _, _, t in mine)), len(fam), sorted((m, n) for m, n, _ in mine) == fam))
O('CENSUS_TABLE_G3.tsv: %d rows; (m, n, type) identical to the brute force: %s' % (len(tab), sorted((int(r[0]), int(r[1]), r[2]) for r in tab) == sorted(mine)))
c = collections.Counter((r[2], r[3]) for r in tab)
for t in ('square', 'twice'): O('  %s: %s' % (t, ', '.join('%s %d' % (s, c[(t, s)]) for s in ('second descent', 'rank 0', 'Demjanenko-Manin', 'special', 'open') if c[(t, s)])))
b13 = [r for r in tab if int(r[0]) > 13000]
O('block 13000-14000: %d rows; %s' % (len(b13), dict(collections.Counter(r[3] for r in b13))))
st = json.load(open(os.path.join(G3, 'stage_13000.json')))
O('  stage_13000.json: ratios %d, after second descent %d, after rank 0 %d, after DM %d, open %s' % (st['ratios'], st['after_second_descent'], st['after_rank0'], st['after_dm'], st['open']))
row12007 = next(r for r in tab if r[:2] == ['12007', '1799'])
O('row (12007,1799) keeps run 6\'s label: %s' % row12007[5])
recs = []
for f in sorted(glob.glob(os.path.join(G3, 'sd_*_*.jsonl'))) + sorted(glob.glob(os.path.join(G3, 'dm_*_*.jsonl'))):
    for l in open(f):
        if l.startswith('{'):
            o = json.loads(l); recs.append(o)
            if o.get('partner'): recs.append(o['partner'])
r0 = collections.defaultdict(list)
for f in sorted(glob.glob(os.path.join(G3, 'r0_*_out*.json'))):
    for k, v in json.load(open(f)).items(): r0[tuple(map(int, k.split(',')))] += v
f0 = {}
for l in open(os.path.join(EXP, 'census-wide', 'table', 'CENSUS_TABLE_1000_11000.tsv')).read().splitlines()[1:]:
    a = l.split('\t'); f0[(int(a[0]), int(a[1]))] = a[3:]
OPEN = [(int(r[0]), int(r[1])) for r in tab if r[3] == 'open']
O('OPEN in the table: %s' % OPEN)
for p in OPEN:
    q = partner(*p); parts = []
    for x in (p, q):
        sd = sorted(set(o['status'] for o in recs if tuple(o['pair']) == x and 'curve' not in o))
        dm = sorted(set(o['status'] for o in recs if tuple(o['pair']) == x and 'curve' in o))
        rz = [v['curve'] for v in r0.get(x, []) if v.get('closes')]
        parts.append('%s: second descent %s, DM %s, rank 0 closes %s%s' % (x, sd, dm, rz or 'none', ', 1000-11000 table: %s %s' % (f0[x][0], f0[x][2]) if x in f0 else ''))
    O('  %s (%s), partner %s: %s' % (p, kind(*p), q, '; '.join(parts)))
q = partner(13103, 4271)
O('route row (13103,4271), %s: partner %s, %s (%d^2 + %d^2 = %d^2), its 1000-11000 row: %s | %s' % (
    kind(13103, 4271), q, kind(*q), q[0], q[1], isqrt(q[0] ** 2 + q[1] ** 2), ' '.join(f0[q][:2]), f0[q][2]))
O('  G3 row: %s' % ' | '.join(next(r for r in tab if r[:2] == ['13103', '4271'])[3:]))
O.finish()
