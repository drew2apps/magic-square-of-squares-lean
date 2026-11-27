#!/usr/bin/env python3
"""F0: one table of every ratio m : n with 1000 < m <= 11000 that passes the congruences (Result 1b of G1_FULL_SYSTEM.md),
with the step that closes it and that step's label.  Every other coprime ratio in the range closes at the congruence
step itself.  Written by an AI model (Claude); not read by a person.  Experiment, not counted: nothing here changes
the census, which stops at m <= 1000.

    python3 census_table.py              print the summary (must equal census_table.txt) and check the table
    python3 census_table.py --write      also rewrite census_table.txt and CENSUS_TABLE_1000_11000.tsv
    python3 census_table.py --negative   altered inputs that must be refused (print; --write: census_table_negative.txt)

It only reads records written earlier; it runs no curve computation.  Inputs:
  congruence lists: ../../census2000/p2000.json (1000 < m <= 2000), pairs_2000_3000.json (this folder, gen_pairs.py
    2000 3000 at 24d6ef6), ../../census-wide/pairs_{3000_5000,5000_7000,7000_9000,9000_11000}.json;
  other type: the repository's class criterion and four-class descent, as committed in scripts/cover_descent/
    (tuple_scan_5000.txt for m <= 5000; structure/class_selmer_15000_list.txt and tuple_selmer_15000.txt above).  The
    session reran those programs at 24d6ef6 (class_step.sh) and every output equalled the committed file (class_step/
    log.txt); with MSS_CLASS_DIR=DIR the table is built from such a rerun instead (DIR/tuple_scan_5000.txt,
    list_5000_10000.txt, list_10000_11000.txt, tuple_selmer_11000.txt).  [Packaging edit, run 6: the default used to be
    the rerun's copies in class_step/, byte-identical to the committed files on this range and not committed again.]
  degenerate type, in this order: second descent (repository wide_*.jsonl for m <= 3000; ../../census-wide
    sw_/s7_/s9_/s9r_/s11_*.jsonl above), directly or through the partner; the repository's L-values (wide_*lvalue*.jsonl);
    rank-0 quotients (global.py records: census2000 global*.json, w3kg*.json, census-wide g5k/g7k/g9k/g11k_out*.json), on
    the ratio or its partner; Demjanenko-Manin (dm*.jsonl, 'closed', for the ratio it serves, key 'of'); then the
    closures listed in SPECIAL from the write-ups (route, conditional, numerical L-value, open).
The table assigns the FIRST step in that order that closes the ratio; a later step may close it too.
"""
import os, sys, json, glob
from math import gcd, isqrt
WRITE = '--write' in sys.argv
HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.abspath(os.path.join(HERE, '..', '..'))
CW = os.path.join(EXP, 'census-wide'); C2 = os.path.join(EXP, 'census2000')
REPO = os.environ.get('MSS_REPO', os.path.abspath(os.path.join(HERE, '..', '..', '..')))
SD_REPO = os.path.join(REPO, 'scripts', 'second_descent')
CS = os.environ.get('MSS_CLASS_DIR')          # a rerun of class_step.sh; default: the repository's committed files
CD = os.path.join(REPO, 'scripts', 'cover_descent')
if CS: CLASS_TUPLE = [os.path.join(CS, f) for f in ('tuple_scan_5000.txt', 'tuple_selmer_11000.txt')]; CLASS_LIST = [os.path.join(CS, f) for f in ('list_5000_10000.txt', 'list_10000_11000.txt')]
else: CLASS_TUPLE = [os.path.join(CD, f) for f in ('tuple_scan_5000.txt', 'tuple_selmer_15000.txt')]; CLASS_LIST = [os.path.join(CD, 'structure', 'class_selmer_15000_list.txt')]

def issq(x): r = isqrt(x); return r * r == x
def typ(m, n):
    s = m * m + n * n
    return 'square' if issq(s) else 'twice' if s % 2 == 0 and issq(s // 2) else 'other'
def partner(m, n):
    a, b = m + n, abs(m - n); g = gcd(a, b); a //= g; b //= g
    return (max(a, b), min(a, b))

L_CONG = 'computed (exact)'
L_CLASS = 'repository programs rerun here at 24d6ef6 (class_selmer_linalg.py, tuple_descent.py); tier D, unreviewed'
L_SD = 'computed; PARI rank upper bounds'
L_SD_REPO = 'repository (wide.py); PARI rank upper bounds'
L_LV = 'repository (wide L-values); numerical L-value with a tail bound'
L_R0 = 'computed (repository global.py); PARI rank upper bounds'
L_DM = 'computed (dm_pari.py); conditional on PARI generators (saturation at primes <= 1000)'
L_ROUTE = 'computed (route_pari.py); conditional on (a) G generating E1(Q) mod torsion, (b) the reading of CVV constants'
L_RROUTE = 'repository route (route_wide.txt / route_twist.txt); tier C/D'
L_COND = 'CONDITIONAL on a rank being 1 (root number -1, PARI bounds [1,3]); plus (a), (b)'
L_COND_DM = 'CONDITIONAL on a rank being 1 (D1 of the partner: root number -1, PARI bounds [1,3]); plus the DM generator condition'   # run 6 (review x4): was L_COND
L_NUM = 'rank 0 from L(E,1) != 0 computed (numerical), tail bound estimated'
L_OPEN = 'OPEN'
# closures written up in census2000/CENSUS3000.md and census-wide/CENSUS5000/7000/9000.md, and the 9000-11000 leftovers
SPECIAL = {
    (1064, 423): ('repository Theorem A route (route_wide.txt), class d = 1906334', L_RROUTE),
    (1487, 641): ('partner (1064,423) closed by the repository route', L_RROUTE),
    (3833, 977): ('partner (2405,1428) closed by the repository route on the twist (route_twist.txt)', L_RROUTE),
    (3631, 1031): ('partner (2331,1300) closed by the repository route (route_wide.txt)', L_RROUTE),
    (3321, 2440): ('rank-1 route, route_pari.py', L_ROUTE),
    (5761, 881): ('partner (3321,2440) closed by the rank-1 route', L_ROUTE),
    (5073, 2336): ('rank-1 route, route_pari.py, d = 40556866', L_ROUTE),
    (6588, 805): ('rank-1 route, route_pari.py, d = -85507438', L_ROUTE),
    (7393, 5783): ('partner (6588,805) closed by the rank-1 route', L_ROUTE),
    (7409, 2737): ('partner (5073,2336) closed by the rank-1 route', L_ROUTE),
    (7811, 3060): ('rank-1 route, route_pari.py, d = 103296242', L_ROUTE),
    (8364, 4123): ('rank-1 route, route_pari.py, d = -105914734', L_ROUTE),
    (8687, 4416): ('rank-1 route, route_pari.py, d = 111925826', L_ROUTE),
    (3540, 2581): ('rank-1 route IF rank E1(Q) = 1', L_COND),
    (6121, 959): ('partner (3540,2581): rank-1 route IF rank E1(Q) = 1', L_COND),
    (4440, 2231): ('rank-1 route IF rank E1(Q) = 1', L_COND),
    (6671, 2209): ('partner (4440,2231): rank-1 route IF rank E1(Q) = 1', L_COND),
    (6480, 4961): ('rank-1 route IF rank E1(Q) = 1', L_COND),
    (6903, 5104): ('Demjanenko-Manin on D1 of the partner (12007,1799) IF its rank is 1', L_COND_DM),
    (8249, 3720): ('rank-1 route IF rank E1(Q) = 1, d = -108415202', L_COND),
    (6557, 324): ('rank-0 quotient with L(E,1) ~ 1.44409 (hard7k/l1_6557.txt)', L_NUM),
    (6881, 6233): ('partner (6557,324): rank-0 quotient with L(E,1) ~ 1.44409', L_NUM),
    (7735, 4968): ('none', L_OPEN),
    # 9000 < m <= 11000 (this run: census-wide/pipe11k.log, route/log_10512_*.txt, route/rh_10512_*.txt)
    (10871, 4751): ('partner (7811,3060) closed by the rank-1 route', L_ROUTE),
    (10512, 145): ('rank-1 route IF rank E1(Q) = 1, d = 220962238 or -220962238 (route/log_10512_*.txt)', L_COND),
    (10657, 10367): ('partner (10512,145): rank-1 route IF rank E1(Q) = 1', L_COND),
}
# the repository route's own closures for m <= 3000 (route_wide_summary.txt), used only if nothing earlier closes
RROUTE_LINE = ('closed by the route (route_wide.txt)', 'closed by the route on the twist', 'their partners in the list')

def jsonl(paths):
    out = []
    for p in paths:
        for l in open(p):
            l = l.strip()
            if l.startswith('{'): out.append(json.loads(l))
    return out

def load(drop=None):
    """all inputs; drop = (kind, pair) removes that pair's records of one kind (for --negative)."""
    D = {}
    lists = [json.load(open(os.path.join(C2, 'p2000.json')))['open'],
             json.load(open(os.path.join(HERE, 'pairs_2000_3000.json')))['open']]
    for r in ('3000_5000', '5000_7000', '7000_9000', '9000_11000'):
        lists.append(json.load(open(os.path.join(CW, 'pairs_%s.json' % r)))['open'])
    D['pairs'] = sorted(set((m, n) for L in lists for m, n in L))
    # class step: ratios with nothing left (m <= 5000: every ratio avoiding the 484 classes; above: the listed ones)
    nothing = set(); listed = set()
    for f in CLASS_TUPLE:
        for l in open(f):
            k = tuple(int(x) for x in l.split(':')[0].strip('()').split(','))
            if l.rstrip().endswith('nothing is left'): nothing.add(k)
    for f in CLASS_LIST:
        for l in open(f):
            listed.add(tuple(int(x) for x in l.split()[0].split(',')))
    scanned = set(tuple(int(x) for x in l.split()[0].split(',')) for l in open(os.path.join(REPO, 'scripts', 'gauss_descent', 'scan_certificates.txt')))
    D['nothing'], D['listed'], D['scanned'] = nothing, listed, scanned
    # second descent
    sd = {}
    for p in sorted(glob.glob(os.path.join(SD_REPO, 'wide_[12]*_[23]*.jsonl'))):
        if 'lvalue' in p: continue
        for r in jsonl([p]): sd.setdefault(tuple(r['pair']), (r, 'repo'))
    for pat in ('sw_*.jsonl', 's7_*.jsonl', 's9_*.jsonl', 's9r_*.jsonl', 's11_*.jsonl'):
        for r in jsonl(sorted(glob.glob(os.path.join(CW, pat)))):
            k = tuple(r['pair'])
            if k not in sd or sd[k][0]['status'] == 'open': sd[k] = (r, 'here')
    lv = {}
    for r in jsonl(sorted(glob.glob(os.path.join(SD_REPO, 'wide*lvalue*.jsonl')))):
        if r['status'] != 'open': lv[tuple(r['pair'])] = r
    r0 = {}
    for p in [os.path.join(C2, 'global_all.json')] + sorted(glob.glob(os.path.join(C2, 'w3kg*.json'))) + \
             sorted(glob.glob(os.path.join(CW, 'g[579]k_out*.json'))) + sorted(glob.glob(os.path.join(CW, 'g11k_out*.json'))):
        for key, curves in json.load(open(p)).items():
            k = tuple(int(x) for x in key.split(','))
            cl = [c['curve'] for c in curves if c.get('closes')]
            if cl: r0.setdefault(k, cl[0])
    dm = {}
    for r in jsonl(sorted(glob.glob(os.path.join(C2, 'dm[0-9].jsonl'))) + sorted(glob.glob(os.path.join(C2, 'w3kdm*.jsonl'))) +
                   sorted(glob.glob(os.path.join(CW, 'dm[579]k[0-9].jsonl'))) + sorted(glob.glob(os.path.join(CW, 'dm11k[0-9].jsonl')))):
        if r.get('status') == 'closed': dm.setdefault(tuple(r['of']), '%s of %s' % (r['curve'], tuple(r['pair'])))
    rr = {}
    for l in open(os.path.join(SD_REPO, 'route_wide_summary.txt')):
        for tag in RROUTE_LINE:
            if l.startswith(tag):
                for tok in l.split(':', 2)[-1].split():
                    for q in tok.split('~'): rr[tuple(int(x) for x in q.split(','))] = tag
    D.update(sd=sd, lv=lv, r0=r0, dm=dm, rr=rr)
    for kind, k in (drop or []):
        if kind == 'class': D['nothing'].discard(k); D['listed'].add(k)
        else: D[kind].pop(k, None)
    return D

def close(D, k):
    """(step, detail, label) for a ratio that passed the congruences."""
    m, n = k; t = typ(m, n)
    if t == 'other':
        if m <= 5000:
            if k in D['nothing']: return ('class', 'four-class descent leaves nothing (tuple_scan_5000)', L_CLASS)
            return ('UNCLOSED', 'other type, not removed by the four-class descent', L_OPEN)
        if k not in D['listed']: return ('class', 'criterion (class_selmer_linalg): no class', L_CLASS)
        if k in D['nothing']: return ('class', 'criterion keeps a class; four-class descent leaves nothing', L_CLASS)
        return ('UNCLOSED', 'other type, listed, not removed', L_OPEN)
    if k in SPECIAL and SPECIAL[k][1] in (L_COND, L_COND_DM, L_OPEN): return ('special', SPECIAL[k][0], SPECIAL[k][1])
    p = partner(m, n)
    for kk, via in ((k, ''), (p, 'partner %s: ' % (p,))):
        r = D['sd'].get(kk)
        if r and r[0]['status'] == 'closed':
            return ('second descent', via + 'closed by %s' % r[0]['by'], L_SD_REPO if r[1] == 'repo' else L_SD)
        if not via and r and r[0]['status'] == 'closed by partner':
            return ('second descent', 'closed through partner %s' % (tuple(r[0]['partner']['pair']),), L_SD_REPO if r[1] == 'repo' else L_SD)
    for kk, via in ((k, ''), (p, 'partner %s: ' % (p,))):
        if kk in D['lv']: return ('L-value', via + 'repository L-value', L_LV)
    for kk, via in ((k, ''), (p, 'partner %s: ' % (p,))):
        if kk in D['r0']: return ('rank 0', via + D['r0'][kk], L_R0)
    if k in D['dm']: return ('Demjanenko-Manin', D['dm'][k], L_DM)
    if k in SPECIAL: return ('special', SPECIAL[k][0], SPECIAL[k][1])
    if k in D['rr']: return ('repository route', D['rr'][k], L_RROUTE)
    return ('UNCLOSED', 'degenerate type, no closing record', L_OPEN)

def build(D):
    rows = []
    for k in D['pairs']:
        if not (1000 < k[0] <= 11000): continue
        rows.append((k, typ(*k)) + close(D, k))
    return rows

def summary(rows):
    from collections import Counter, OrderedDict
    out = []
    rng = [(1000, 2000), (2000, 3000), (3000, 5000), (5000, 7000), (7000, 9000), (9000, 11000)]
    out.append('F0 census table, 1000 < m <= 11000 (experiment; written by an AI model, not read by a person)')
    out.append('ratios passing the congruences: %d' % len(rows))
    for a, b in rng:
        R = [r for r in rows if a < r[0][0] <= b]
        c = Counter(r[1] for r in R)
        out.append('  %d < m <= %d: %d (other %d, square %d, twice-square %d)' % (a, b, len(R), c['other'], c['square'], c['twice']))
    out.append('closing step (first that closes, in the order of the header):')
    c = Counter(r[2] for r in rows)
    for s in ('class', 'second descent', 'L-value', 'rank 0', 'Demjanenko-Manin', 'repository route', 'special', 'UNCLOSED'):
        out.append('  %-17s %d' % (s, c[s]))
    out.append('by label:')
    for lab, cnt in sorted(Counter(r[4] for r in rows).items()):
        out.append('  %5d  %s' % (cnt, lab))
    out.append('conditional on a rank being 1: ' + ' '.join('(%d,%d)' % r[0] for r in rows if r[4] in (L_COND, L_COND_DM)))
    out.append('numerical L-value: ' + ' '.join('(%d,%d)' % r[0] for r in rows if r[4] == L_NUM))
    out.append('open: ' + ' '.join('(%d,%d)' % r[0] for r in rows if r[4] == L_OPEN and r[2] != 'UNCLOSED'))
    un = [r for r in rows if r[2] == 'UNCLOSED']
    out.append('UNCLOSED (a record is missing): %d %s' % (len(un), ' '.join('(%d,%d)' % r[0] for r in un[:30])))
    return out

def tsv(rows):
    L = ['m\tn\ttype\tstep\tdetail\tlabel']
    for (m, n), t, s, d, lab in rows: L.append('%d\t%d\t%s\t%s\t%s\t%s' % (m, n, t, s, d, lab))
    return '\n'.join(L) + '\n'

def main():
    if '--negative' in sys.argv:
        D0 = load(); base = {r[0]: r for r in build(D0)}
        out = ['F0 census table, negative mode: with every record of a ratio and of its partner removed, the ratio must come out UNCLOSED']
        want = [('class', lambda r: r[2] == 'class' and r[0][0] <= 5000),
                ('class', lambda r: r[2] == 'class' and 'keeps a class' in r[3]),
                ('sd', lambda r: r[2] == 'second descent' and r[3].startswith('closed by')),
                ('lv', lambda r: r[2] == 'L-value' and not r[3].startswith('partner')),
                ('r0', lambda r: r[2] == 'rank 0' and not r[3].startswith('partner')),
                ('dm', lambda r: r[2] == 'Demjanenko-Manin'),
                ('rr', lambda r: r[2] == 'repository route')]
        tests = [(kind, next(r[0] for r in sorted(base.values()) if f(r))) for kind, f in want]
        ok = True
        for kind, k in tests:
            p = partner(*k)
            D = load(drop=[(kd, q) for kd in ('class', 'sd', 'lv', 'r0', 'dm', 'rr') for q in (k, p)]); new = close(D, k)
            good = new[0] == 'UNCLOSED'
            ok &= good
            out.append('  drop every record of %s and its partner (%s step): before "%s", after "%s" -> %s' % (k, kind, base[k][2], new[0], 'refused' if good else 'NOT REFUSED'))
        out.append('all refused: %s' % ok)
        print('\n'.join(out))
        if WRITE: open(os.path.join(HERE, 'census_table_negative.txt'), 'w').write('\n'.join(out) + '\n')
        return
    rows = build(load())
    S = summary(rows)
    print('\n'.join(S))
    if WRITE:
        open(os.path.join(HERE, 'census_table.txt'), 'w').write('\n'.join(S) + '\n')
        open(os.path.join(HERE, 'CENSUS_TABLE_1000_11000.tsv'), 'w').write(tsv(rows))
    else:
        f = os.path.join(HERE, 'CENSUS_TABLE_1000_11000.tsv')
        if os.path.exists(f): assert open(f).read() == tsv(rows), 'the table differs from the committed CENSUS_TABLE_1000_11000.tsv'
        g = os.path.join(HERE, 'census_table.txt')
        if os.path.exists(g): assert open(g).read() == '\n'.join(S) + '\n', 'the summary differs from census_table.txt'
        print('# matches census_table.txt and CENSUS_TABLE_1000_11000.tsv')

if __name__ == '__main__':
    main()
