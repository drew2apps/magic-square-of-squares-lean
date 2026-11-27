#!/usr/bin/env python3
"""Review x4 (run 6): an independent recount of experiments/census-wide/table/CENSUS_TABLE_1000_11000.tsv from the
committed records, with code written separately from census_table.py.  Written by an AI model; not read by a person.
The class step is read from the repository's own committed outputs (scripts/cover_descent/...).  For every ratio it
collects every machine-readable closing record, takes the first in census_table.py's order, and compares type and step
with the TSV.  Prints only.   python3 recount.py [--negative]   (--negative: forged inputs that must be detected)"""
import json, glob, os, sys, collections
from math import gcd, isqrt
NEG = '--negative' in sys.argv
R = os.environ.get('MSS_REPO', os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
E = R + '/experiments'; CW = E + '/census-wide'; C2 = E + '/census2000'; SD = R + '/scripts/second_descent'
def sq(x): r = isqrt(x); return r*r == x
def typ(m, n):
    s = m*m + n*n
    if sq(s): return 'square'
    if s % 2 == 0 and sq(s//2): return 'twice'
    return 'other'
def part(m, n):
    a, b = m+n, abs(m-n)
    while a % 2 == 0 and b % 2 == 0: a //= 2; b //= 2
    g = gcd(a, b); a //= g; b //= g
    return (max(a, b), min(a, b))
def jl(p):
    for l in open(p):
        if l.startswith('{'): yield json.loads(l)
# congruence lists
L = set()
srcs = [C2 + '/p2000.json', CW + '/table/pairs_2000_3000.json'] + [CW + '/pairs_%s.json' % r for r in ('3000_5000','5000_7000','7000_9000','9000_11000')]
for f in srcs:
    for m, n in json.load(open(f))['open']: L.add((m, n))
P = sorted(k for k in L if 1000 < k[0] <= 11000)
bad = [k for k in P if gcd(*k) != 1 or not (k[0] > k[1] > 0)]
print('ratios 1000<m<=11000 in the congruence lists:', len(P), 'non-coprime or misordered:', len(bad))
T = collections.Counter(typ(*k) for k in P); print('types', dict(T))
# class step from the REPOSITORY's committed files
cd = R + '/scripts/cover_descent'
nothing = set()
for f in (cd + '/tuple_scan_5000.txt', cd + '/tuple_selmer_15000.txt'):
    for l in open(f):
        k = tuple(int(x) for x in l.split(':')[0].strip('()').split(','))
        if l.rstrip().endswith('nothing is left'): nothing.add(k)
listed = set(tuple(int(x) for x in l.split()[0].split(',')) for l in open(cd + '/structure/class_selmer_15000_list.txt'))
scanned = set(tuple(int(x) for x in l.split()[0].split(',')) for l in open(R + '/scripts/gauss_descent/scan_certificates.txt'))
# evidence per ratio
ev = collections.defaultdict(list)
sd = {}
for p in sorted(glob.glob(SD + '/wide_*.jsonl')):
    if 'lvalue' in p: continue
    for r in jl(p): sd.setdefault(tuple(r['pair']), []).append((r, 'repo'))
for p in sorted(glob.glob(CW + '/s*_*.jsonl')):
    b = os.path.basename(p)
    if not any(b.startswith(x) for x in ('sw_', 's7_', 's9_', 's9r_', 's11_')): continue
    for r in jl(p): sd.setdefault(tuple(r['pair']), []).append((r, b))
sdclosed = {k for k, v in sd.items() if any(r['status'] == 'closed' for r, _ in v)}
sdbyp = {}
for k, v in sd.items():
    for r, _ in v:
        if r['status'] == 'closed by partner':
            pp = tuple(r['partner']['pair']); assert pp == part(*k), (k, pp)
            assert r['partner']['status'] == 'closed', k
            sdbyp[k] = pp
lv = set()
for p in sorted(glob.glob(SD + '/wide*lvalue*.jsonl')):
    for r in jl(p):
        if r['status'] == 'closed by L-value': lv.add(tuple(r['pair']))
        else: assert r['status'] == 'open'
r0 = set()
for p in [C2 + '/global_all.json'] + sorted(glob.glob(C2 + '/w3kg*.json')) + sorted(glob.glob(CW + '/g*k_out*.json')):
    for key, cs in json.load(open(p)).items():
        if any(c.get('closes') is True for c in cs): r0.add(tuple(int(x) for x in key.split(',')))
dmc = collections.defaultdict(set)
for p in sorted(glob.glob(C2 + '/dm[0-9].jsonl')) + sorted(glob.glob(C2 + '/w3kdm*.jsonl')) + sorted(glob.glob(CW + '/dm*k[0-9].jsonl')):
    for r in jl(p):
        if r.get('status') == 'closed': dmc[tuple(r['of'])].add(tuple(r['pair']))
# compare with the TSV
rows = {}
for l in open(CW + '/table/CENSUS_TABLE_1000_11000.tsv').read().splitlines()[1:]:
    m, n, t, s, d, lab = l.split('\t'); rows[(int(m), int(n))] = (t, s, d, lab)
def compare(P, rows, sdclosed, sdbyp, lv, r0, dmc, nothing, listed):
  mine = collections.Counter(); mism = []
  for k in P:
    t = typ(*k); pk = part(*k)
    if t == 'other':
        if k[0] <= 5000: st = 'class' if k in nothing else 'UNCLOSED'
        else: st = 'class' if (k not in listed or k in nothing) else 'UNCLOSED'
        if k[0] <= 5000: assert k in scanned, k
    else:
        st = None
        if k in sdclosed or pk in sdclosed or k in sdbyp: st = 'second descent'
        elif k in lv or pk in lv: st = 'L-value'
        elif k in r0 or pk in r0: st = 'rank 0'
        elif k in dmc: st = 'Demjanenko-Manin'
        else: st = 'NONE'
    mine[st] += 1
    if rows[k][0] != t: mism.append((k, 'type', t, rows[k][0]))
    tst = rows[k][1]
    if st != tst and not (st == 'NONE' and tst in ('special', 'repository route')): mism.append((k, st, tst))
  return mine, mism
print('TSV rows', len(rows), 'same key set as lists:', set(rows) == set(P))
if NEG:
    out = []
    r2 = dict(rows); k = next(x for x in P if rows[x][1] == 'rank 0'); r2[k] = (rows[k][0], 'Demjanenko-Manin', '', ''); out.append(('a rank-0 row relabelled DM', k, r2, {}))
    k = next(x for x in P if rows[x][1] == 'second descent' and x in sdclosed and part(*x) not in sdclosed and x not in lv and x not in r0 and part(*x) not in r0 and x not in dmc and part(*x) not in lv)
    out.append(('the second-descent record of a ratio removed', k, rows, dict(sdclosed=sdclosed - {k})))
    k = next(x for x in P if typ(*x) == 'other' and x[0] > 5000 and x in listed); out.append(('a listed other-type ratio without its four-class record', k, rows, dict(nothing=nothing - {k})))
    k = next(x for x in P if typ(*x) == 'other' and x[0] > 5000 and x not in listed); out.append(('an unlisted other-type ratio added to the list', k, rows, dict(listed=listed | {k})))
    allok = True
    for name, k, rr, forged in out:
        args = dict(sdclosed=sdclosed, sdbyp=sdbyp, lv=lv, r0=r0, dmc=dmc, nothing=nothing, listed=listed); args.update(forged)
        _, mm = compare(P, rr, **args); good = any(x[0] == k for x in mm); allok &= good
        print('forged: %s %s -> %s' % (name, k, 'detected' if good else 'NOT DETECTED'))
    print('negative checks: %s' % ('all detected' if allok else 'FAILED')); assert allok
    sys.exit(0)
mine, mism = compare(P, rows, sdclosed, sdbyp, lv, r0, dmc, nothing, listed)
print('my first-step counts (NONE = no machine record; must be the special/route rows):', dict(mine))
print('TSV step counts:', dict(collections.Counter(v[1] for v in rows.values())))
print('mismatches:', len(mism), mism[:10])
none = [k for k in P if typ(*k) != 'other' and not (k in sdclosed or part(*k) in sdclosed or k in sdbyp or k in lv or part(*k) in lv or k in r0 or part(*k) in r0 or k in dmc)]
print('degenerate ratios with no machine closing record:', len(none))
print(' '.join('%d,%d:%s' % (k[0], k[1], rows[k][1][:7]) for k in none))
# label of each TSV row is one of the defined ones and conditional steps are labelled so
labs = collections.Counter(v[3] for v in rows.values())
for lab, c in sorted(labs.items(), key=lambda x: -x[1]): print('%6d %s' % (c, lab))
# DM via partner's 'of' (not used by the table): extra robustness info
print('DM closed for partner only (not used):', sum(1 for k in P if typ(*k) != 'other' and k not in dmc and part(*k) in dmc))
assert len(P) == 20682 and not bad and set(rows) == set(P) and not mism and len(none) == 35
assert dict(mine) == {'class': 17636, 'second descent': 2582, 'rank 0': 231, 'Demjanenko-Manin': 175, 'L-value': 23, 'NONE': 35}
print('recount: as the table')
