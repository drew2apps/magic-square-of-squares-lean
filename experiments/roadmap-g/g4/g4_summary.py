"""Summary of the G4 search (bookkeeping over g4_<i>_<N>.json, the first pass, and pass2_<i>_<N>.jsonl, the second pass).
usage: g4_summary.py [--write]   (checks g4_summary.txt; --write rewrites it)"""
import json, glob, os, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__))
recs = []
for f in sorted(glob.glob(os.path.join(HERE, 'g4_*_*.json'))): recs += json.load(open(f))
recs.sort(key=lambda r: tuple(r['ratio']))
p2 = {}
for f in sorted(glob.glob(os.path.join(HERE, 'pass2_*_*.jsonl'))):
    for l in open(f):
        if l.startswith('{'): r = json.loads(l); p2[tuple(r['ratio'])] = r
RMAX = 3
ORDER = ['none', 'partial', 'full']
def state(v, w):
    if v.get('searched'): return 'full'
    if w and w.get('searched') == 'partial': return 'partial'
    if w and w.get('searched') is True: return 'full'
    return 'none'
cov = collections.defaultdict(collections.Counter); none = []; other = []; sols = []; planted = [0, 0]; Epm = collections.Counter()
for r in recs:
    q = p2.get(tuple(r['ratio']), {'curves': {}, 'hits': [], 'solutions': []}); k = 0
    for c, v in r['curves'].items():
        lo, hi = v['rank']; w = q['curves'].get(c); s = state(v, w)
        if s == 'full': st = 'full' if v.get('searched') else 'full (pass 2)'
        elif s == 'partial': st = 'partial (pass 2)'
        elif lo != hi: st = 'skipped: rank bounds differ'
        elif lo > RMAX: st = 'skipped: rank above %d' % RMAX
        else: st = 'skipped: no point found' if w else 'skipped: too few points, not in pass 2'
        cov[c][st] += 1; k += s != 'none'
        for vv in [v] + ([w] if w else []):
            p = vv.get('planted', {})
            if vv.get('searched') and p.get('carried'): planted[0] += 1; planted[1] += bool(p.get('found'))
    Epm[max((state(r['curves'][c], q['curves'].get(c)) for c in ('E+', 'E-')), key=ORDER.index)] += 1
    if k == 0: none.append(tuple(r['ratio']))
    for h in r['hits'] + q['hits']:
        if h['t'] not in ('1', '-1'): other.append((tuple(r['ratio']), h['curve'], h['t'], h['squares'], h['verdict']))
    sols += [(tuple(r['ratio']), t) for t in r['solutions'] + q['solutions']]
L = ['square-type ratios: %d (m <= %d); second pass (E+, E-) run on %d of them' % (len(recs), max(r['ratio'][0] for r in recs), len(p2)),
     'SOLUTIONS: %s' % (sols or 'none'),
     'non-degenerate t with at least four of the five forms squares: %s' % (sorted(set(other)) or 'none'),
     'ratios by their best search on E+ or E-: ' + ', '.join('%s %d' % (s, Epm[s]) for s in reversed(ORDER)),
     'ratios with no curve searched at all: %d %s' % (len(none), none[:40]),
     'planted test: searches on a curve carrying the point of t = +-1: %d; in which it is listed: %d' % tuple(planted),
     'per curve (full: rank bounds equal and a full set of PARI points, every point of height <= 64 max h(gens) listed;',
     '  partial: rank bounds equal but fewer points, the subgroup they generate searched the same way):']
for c in sorted(cov):
    L.append('  %-11s ' % c + ', '.join('%s %d' % kv for kv in sorted(cov[c].items(), key=lambda kv: (-kv[1], kv[0]))))
out = '\n'.join(L) + '\n'; f = os.path.join(HERE, 'g4_summary.txt')
if '--write' in sys.argv: open(f, 'w').write(out)
print(out, end='')
if '--write' not in sys.argv and os.path.exists(f): assert open(f).read() == out; print('# matches g4_summary.txt')
