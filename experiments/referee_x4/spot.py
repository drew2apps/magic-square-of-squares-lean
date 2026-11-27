#!/usr/bin/env python3
"""Review x4 (run 6): spot-check of experiments/census-wide/table/CENSUS_TABLE_1000_11000.tsv.  Samples rows of every
step kind (seeded) and prints the record behind each with a test of what "closes" means for that kind; then checks the
26 hand-entered SPECIAL rows against their route records and partner relations.  Reads committed files; prints only.
Written by an AI model; not read by a person.   python3 spot.py"""
import json, glob, os, sys, random, collections
R = os.environ.get('MSS_REPO', os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
E = R + '/experiments'; CW = E + '/census-wide'; C2 = E + '/census2000'; SD = R + '/scripts/second_descent'
cd = R + '/scripts/cover_descent'
from math import gcd
def part(m, n):
    a, b = m+n, abs(m-n); g = gcd(a, b); a //= g; b //= g; return (max(a, b), min(a, b))
def jl(p):
    for l in open(p):
        if l.startswith('{'): yield json.loads(l)
rows = {}
for l in open(CW + '/table/CENSUS_TABLE_1000_11000.tsv').read().splitlines()[1:]:
    m, n, t, s, d, lab = l.split('\t'); rows[(int(m), int(n))] = (t, s, d, lab)
sd = collections.defaultdict(list)
for p in sorted(glob.glob(SD + '/wide_*.jsonl')) + sorted(glob.glob(CW + '/s*_*.jsonl')):
    if 'lvalue' in p: continue
    for r in jl(p): sd[tuple(r['pair'])].append((os.path.relpath(p, R), r))
lv = collections.defaultdict(list)
for p in sorted(glob.glob(SD + '/wide*lvalue*.jsonl')):
    for r in jl(p): lv[tuple(r['pair'])].append((os.path.relpath(p, R), r))
g = collections.defaultdict(list)
for p in [C2 + '/global_all.json'] + sorted(glob.glob(C2 + '/w3kg*.json')) + sorted(glob.glob(CW + '/g*k_out*.json')):
    for key, cs in json.load(open(p)).items(): g[tuple(map(int, key.split(',')))].append((os.path.relpath(p, R), cs))
dm = collections.defaultdict(list)
for p in sorted(glob.glob(C2 + '/dm[0-9].jsonl')) + sorted(glob.glob(C2 + '/w3kdm*.jsonl')) + sorted(glob.glob(CW + '/dm*k[0-9].jsonl')):
    for r in jl(p): dm[tuple(r['of'])].append((os.path.relpath(p, R), r))
ts = {}
for f in (cd + '/tuple_scan_5000.txt', cd + '/tuple_selmer_15000.txt'):
    for l in open(f): ts[tuple(int(x) for x in l.split(':')[0].strip('()').split(','))] = (os.path.relpath(f, R), l.strip())
listed = set(tuple(int(x) for x in l.split()[0].split(',')) for l in open(cd + '/structure/class_selmer_15000_list.txt'))
rng = random.Random(26)
def pick(f, k):
    c = sorted(x for x in rows if f(x, rows[x])); return rng.sample(c, min(k, len(c)))
S = []
S += [('class<=5000', x) for x in pick(lambda x, r: r[1] == 'class' and x[0] <= 5000, 3)]
S += [('class no class', x) for x in pick(lambda x, r: r[1] == 'class' and 'no class' in r[2], 3)]
S += [('class kept+descent', x) for x in pick(lambda x, r: r[1] == 'class' and 'keeps' in r[2], 2)]
S += [('sd repo direct', x) for x in pick(lambda x, r: r[1] == 'second descent' and r[2].startswith('closed by') and 'wide.py' in r[3], 2)]
S += [('sd here direct', x) for x in pick(lambda x, r: r[1] == 'second descent' and r[2].startswith('closed by') and 'wide.py' not in r[3], 3)]
S += [('sd partner rec', x) for x in pick(lambda x, r: r[1] == 'second descent' and r[2].startswith('partner'), 2)]
S += [('sd through partner', x) for x in pick(lambda x, r: r[1] == 'second descent' and r[2].startswith('closed through'), 2)]
S += [('L-value', x) for x in pick(lambda x, r: r[1] == 'L-value' and not r[2].startswith('partner'), 1)]
S += [('L-value partner', x) for x in pick(lambda x, r: r[1] == 'L-value' and r[2].startswith('partner'), 1)]
S += [('rank0', x) for x in pick(lambda x, r: r[1] == 'rank 0' and not r[2].startswith('partner'), 2)]
S += [('rank0 partner', x) for x in pick(lambda x, r: r[1] == 'rank 0' and r[2].startswith('partner'), 2)]
S += [('DM', x) for x in pick(lambda x, r: r[1] == 'Demjanenko-Manin', 3)]
S += [('repo route', x) for x in pick(lambda x, r: r[1] == 'repository route', 2)]
ok_all = True
for kind, k in S:
    t, s, d, lab = rows[k]; p = part(*k); ok = False; ev = ''
    if kind.startswith('class'):
        if kind == 'class no class': ok = k not in listed; ev = 'not in class_selmer_15000_list.txt: %s' % ok
        else: f, l = ts[k]; ok = l.endswith('nothing is left'); ev = '%s: %s' % (f, l[:110])
    elif kind.startswith('sd'):
        kk = p if kind == 'sd partner rec' else k
        for f, r in sd[kk]:
            if r['status'] == 'closed' and min(r['rank_bounds']) == 0 and set(r['rational_t'].values()) <= {'no', 'degenerate'}:
                ok = True; ev = '%s: %s closed by %s rank_bounds %s t %s' % (f, kk, r['by'], r['rank_bounds'], r['rational_t'])
            if r['status'] == 'closed by partner':
                q = r['partner']
                ok = q['status'] == 'closed' and min(q['rank_bounds']) == 0 and set(q['rational_t'].values()) <= {'no', 'degenerate'} and tuple(q['pair']) == p
                ev = '%s: %s closed by partner %s (%s, bounds %s, t %s)' % (f, k, q['pair'], q['by'], q['rank_bounds'], q['rational_t'])
    elif kind.startswith('L-value'):
        kk = p if 'partner' in kind else k
        for f, r in lv[kk]:
            if r['status'] == 'closed by L-value':
                c = r['closing'][0]; ok = float(c['L']) > 10 * float(c['tail_bound']) and set(r['rational_t'].values()) <= {'no', 'degenerate'}
                ev = '%s: %s curve %s twist %s L %s tail %s t %s' % (f, kk, c['curve'], c['twist'], c['L'], c['tail_bound'], r['rational_t'])
    elif kind.startswith('rank0'):
        kk = p if 'partner' in kind else k
        for f, cs in g[kk]:
            for c in cs:
                if c.get('closes') is True: ok = c['rank'][1] == 0; ev = '%s: %s %s' % (f, kk, json.dumps(c)[:160]); break
    elif kind == 'DM':
        for f, r in dm[k]:
            if r.get('status') == 'closed':
                ok = set(r['t'].values()) <= {'degenerate', 'no'} and r['rank'] == [1, 1]
                ev = '%s: of %s pair %s curve %s rank %s hG %.3f ncand %s t %s gens %s' % (f, k, r['pair'], r['curve'], r['rank'], r['hG'], r['ncand'], r['t'], r['gens']); break
    elif kind == 'repo route':
        L = [l for l in open(SD + '/route_wide_summary.txt') if ('%d,%d' % k) in l]; ok = bool(L) and not L[0].startswith('left'); ev = 'route_wide_summary.txt: ' + L[0][:100].strip()
    ok_all &= ok
    print('%-19s %s %-15s ok=%s | %s' % (kind, k, s, ok, ev))
print('spot-checked %d rows; all close: %s' % (len(S), ok_all))

# ---- the SPECIAL rows (hand-entered in census_table.py): each against its route record, and every partner relation
RT = CW + '/route'
def route_rec(m, n, d=None):
    fs = glob.glob(RT + '/log_%d.txt' % m) + glob.glob(RT + '/log_%d_*.txt' % m) + glob.glob(RT + '/out_%d_%d_*.json' % (m, n))
    out = []
    for f in fs:
        txt = open(f).read().strip().splitlines()
        js = [l for l in txt if l.startswith('{') and '"closed"' in l] if f.endswith('.txt') else [open(f).read()]
        for l in js:
            o = json.loads(l)
            if tuple(o['pair']) == (m, n): out.append((os.path.basename(f), o))
    return out
want_closed = [(3321, 2440), (5073, 2336), (6588, 805), (7811, 3060), (8364, 4123), (8687, 4416)]
want_cond = [(3540, 2581), (4440, 2231), (6480, 4961), (8249, 3720), (10512, 145)]
ok2 = True
for k in want_closed + want_cond:
    rs = route_rec(*k)
    cl = [o['closed'] for _, o in rs]; ru = [o.get('rank_upper_isogeny_class') for _, o in rs]
    good = bool(rs) and (all(cl) if k in want_closed else not any(cl)) and all(set(v for F in o['fibre'].values() for v in F['rational_t'].values()) <= {'degenerate'} for _, o in rs)
    ok2 &= good
    print('special %-14s %-9s records %s closed %s rank upper (isogeny class) %s; fibre t degenerate only -> %s' % (k, 'route' if k in want_closed else 'cond', [f for f, _ in rs], cl, ru, good))
pairs = {(5761, 881): (3321, 2440), (7393, 5783): (6588, 805), (7409, 2737): (5073, 2336), (10871, 4751): (7811, 3060),
         (6121, 959): (3540, 2581), (6671, 2209): (4440, 2231), (6881, 6233): (6557, 324), (6903, 5104): (12007, 1799),
         (10657, 10367): (10512, 145), (1487, 641): (1064, 423), (3833, 977): (2405, 1428), (3631, 1031): (2331, 1300)}
pr = all(part(*a) == b for a, b in pairs.items()); ok2 &= pr
print('partner relations of the 12 SPECIAL partner rows: %s' % pr)
for k in ((2405, 1428), (2331, 1300), (1064, 423)):
    L = [l for l in open(SD + '/route_wide_summary.txt') if ('%d,%d' % k) in l and l.startswith('closed')]; ok2 &= bool(L)
    print('special partner %s in route_wide_summary.txt: %s' % (k, L[0][:60].strip() if L else 'MISSING'))
print('special rows checked: all as labelled: %s' % ok2)
assert ok_all and ok2
