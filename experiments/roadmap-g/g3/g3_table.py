"""G3 census table: every ratio of the two degenerate types (m^2 + n^2 a square, or twice a square) with 11000 < m <= M
that passes the congruences (fam_11000_30000.json), for each 1000-block whose stages are finished, with the first step
that closes it, in the format and with the labels of run5/f0's CENSUS_TABLE_1000_11000.tsv.  Bookkeeping only: it runs
no curve computation.  (EXPERIMENT, written by an AI model; not read by a person; nothing here changes the census.)
Steps, in the order tried (a ratio also closes when its partner closes at that step, as in CENSUS9000.md):
  second descent (sd_<lo>_<k>.jsonl), rank 0 (r0_<lo>_out<k>.json), Demjanenko-Manin (dm_<lo>_<k>.jsonl),
  special: the rank-1 route (route/log_<m>.txt), counted for square-type ratios only (CENSUS9000.md: the route is
  built for square type; on a twice-square ratio its fibre step reports "no" rather than "degenerate");
  partner below 11000: a twice-square ratio's partner ((m+n)/2, (m-n)/2) has m below the ratio's; when it is in
  experiments/census-wide/table/CENSUS_TABLE_1000_11000.tsv (F0_TABLE), its closing row is taken, with its label.
  g3_table.py              checks g3_table.txt and CENSUS_TABLE_G3.tsv (--write rewrites them)
  g3_table.py --negative   for one ratio of each step, removes every record of it and of its partner: it must come out
                           OPEN; checks g3_table_negative.txt (--write rewrites it)"""
import sys, os, json, math, glob, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = {'second descent': 'computed; PARI rank upper bounds',
       'rank 0': 'computed (repository global.py); PARI rank upper bounds',
       'Demjanenko-Manin': 'computed (dm_pari.py); conditional on PARI generators (saturation at primes <= 1000)',
       'special': 'computed (route_pari.py); conditional on (a) G generating E1(Q) mod torsion, (b) the reading of CVV constants',
       'open': 'OPEN'}

def partner(m, n):
    a, b = m + n, abs(m - n)
    if a % 2 == 0 and b % 2 == 0: a, b = a // 2, b // 2
    g = math.gcd(a, b); a, b = a // g, b // g
    return (max(a, b), min(a, b))

def kind(m, n):
    s = m*m + n*n; r = math.isqrt(s)
    return 'square' if r*r == s else 'twice'

def jl(path):
    if not os.path.exists(path): return []
    return [json.loads(l) for l in open(path) if l.startswith('{')]

def records(drop=()):
    """closing records: pair -> (step, detail), the first step that closes the pair itself (no partner yet)."""
    sd, r0, dm, rt = {}, {}, {}, {}
    for f in sorted(glob.glob(os.path.join(HERE, 'sd_*_*.jsonl'))):
        for o in jl(f):
            p = tuple(o['pair']); sd[p] = o
            if o.get('partner'): sd.setdefault(tuple(o['partner']['pair']), o['partner'])
    for f in sorted(glob.glob(os.path.join(HERE, 'r0_*_out*.json'))):
        for key, lst in json.load(open(f)).items():
            c = [x['curve'] for x in lst if x.get('closes')]
            if c: r0[tuple(int(t) for t in key.split(','))] = c[0]
    for f in sorted(glob.glob(os.path.join(HERE, 'dm_*_*.jsonl'))):
        for o in jl(f):
            if o.get('status') == 'closed': dm.setdefault(tuple(o['pair']), o['curve'])
    for f in sorted(glob.glob(os.path.join(HERE, 'route', 'log_*.txt'))):
        o = jl(f)
        if o and o[-1].get('closed') and kind(*o[-1]['pair']) == 'square':
            rt[tuple(o[-1]['pair'])] = 'rank-1 route, class d = %d, E1 rank %s, sieve leaves %s (%s)' % (
                o[-1]['d'], o[-1]['rank_E1'], o[-1]['left'], os.path.relpath(f, HERE))
    f0 = {}
    path = os.environ.get('F0_TABLE', os.path.join(HERE, '..', '..', 'census-wide', 'table', 'CENSUS_TABLE_1000_11000.tsv'))
    for l in open(path).read().splitlines()[1:]:
        a = l.split('\t'); f0[(int(a[0]), int(a[1]))] = a[3:]
    for D in (sd, r0, dm, rt, f0):
        for p in drop: D.pop(p, None)
    return sd, r0, dm, rt, f0

def table(drop=()):
    fam = json.load(open(os.path.join(HERE, 'fam_11000_30000.json')))['open']
    blocks = sorted(int(os.path.basename(f)[6:-5]) for f in glob.glob(os.path.join(HERE, 'stage_*.json')))
    sd, r0, dm, rt, f0 = records(drop)
    rows = []
    for lo in blocks:
        for m, n in sorted(tuple(p) for p in fam if lo < p[0] <= lo + 1000):
            p = (m, n); q = partner(m, n); step = detail = None
            o = sd.get(p)
            if o and o['status'] == 'closed': step, detail = 'second descent', 'class %s' % o['by']
            elif sd.get(q) and sd[q]['status'] == 'closed': step, detail = 'second descent', 'partner (%d, %d): class %s' % (q + (sd[q]['by'],))
            for name, D, fmt in (('rank 0', r0, '%s'), ('Demjanenko-Manin', dm, 'curve %s'), ('special', rt, '%s')):
                if step: break
                if p in D: step, detail = name, fmt % D[p]
                elif q in D: step, detail = name, 'partner (%d, %d): ' % q + fmt % D[q]
            label = None
            if not step and q in f0 and f0[q][0] != 'unclosed' and 'OPEN' not in f0[q][2]:
                step, detail, label = f0[q][0], 'partner (%d, %d), 1000 < m <= 11000 table: %s' % (q + (f0[q][1],)), f0[q][2]
            if not step: step, detail = 'open', 'left after every step'
            rows.append((m, n, kind(m, n), step, detail, label or LAB[step]))
    return blocks, rows

def summary(blocks, rows):
    c = collections.Counter((r[2], r[3]) for r in rows)
    L = ['blocks: %s (11000 < m <= %d); ratios: %d' % (', '.join('%d-%d' % (b, b + 1000) for b in blocks), blocks[-1] + 1000 if blocks else 11000, len(rows))]
    steps = list(LAB) + sorted(set(r[3] for r in rows) - set(LAB))
    for t in ('square', 'twice'):
        L.append('%s: ' % t + ', '.join('%s %d' % (s, c[(t, s)]) for s in steps if c[(t, s)]))
    L.append('OPEN: %s' % ([r[:2] for r in rows if r[3] == 'open'] or 'none'))
    L.append('CONDITIONAL on a rank being 1: %s' % ([r[:2] for r in rows if r[5].startswith('CONDITIONAL')] or 'none'))
    L.append('closed by the route (conditions (a), (b)): %s' % ([r[:2] for r in rows if 'route_pari' in r[5]] or 'none'))
    return L

if __name__ == '__main__':
    neg = '--negative' in sys.argv
    if neg:
        blocks, rows = table(); L = []
        for s in sorted(set(r[3] for r in rows) - {'open'}):
            r = next((r for r in rows if r[3] == s), None)
            if not r: continue
            p = r[:2]; _, rows2 = table(drop=(p, partner(*p)))
            r2 = next(x for x in rows2 if x[:2] == p)
            L.append('(%d,%d) %s, records of it and of its partner removed: %s' % (p + (s, 'OPEN (as it must be)' if r2[3] == 'open' else 'still closed by %s (FAILURE)' % r2[3])))
        L.append('# %d of %d as they must be' % (sum('as it must be' in l for l in L), len(L)))
        out = {'g3_table_negative.txt': '\n'.join(L) + '\n'}
    else:
        blocks, rows = table()
        out = {'g3_table.txt': '\n'.join(summary(blocks, rows)) + '\n',
               'CENSUS_TABLE_G3.tsv': 'm\tn\ttype\tstep\tdetail\tlabel\n' + ''.join('\t'.join(map(str, r)) + '\n' for r in rows)}
    for f, s in out.items():
        if '--write' in sys.argv: open(os.path.join(HERE, f), 'w').write(s)
        elif os.path.exists(os.path.join(HERE, f)): assert open(os.path.join(HERE, f)).read() == s, f
    print(out[list(out)[0]], end='')
    if '--write' not in sys.argv: print('# matches ' + ', '.join(out))
