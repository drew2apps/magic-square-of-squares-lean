"""Referee z17 (written by an AI model; not read by a person): I8's checkpoint numbers recounted from
roadmap-i/i8/i8_first360.jsonl with my own code, and its tuples compared with lines 1-360 of i2_residual.jsonl.
Also: the per-tuple verdict recomputed from the per-curve labels (killed / outside Sel2 on some curve -> pairing).
usage: python3 -I z17_i8count.py [--write]    (z17_i8count.txt)"""
import sys, os, json, collections
HERE = os.path.dirname(os.path.abspath(__file__))
RI = os.path.join(HERE, '..', 'roadmap-i')
OUT = []
def say(s=''): OUT.append(s)
say('z17: I8 checkpoint recounted. Written by an AI model (reviewer z17); not read by a person.')
lst = [tuple(map(int, l.split()[0].split(','))) for l in open(os.path.join(RI, 'i2', 'i2_needgood_list.txt')) if l.strip()]
src = [l.split()[1] for l in open(os.path.join(RI, 'i2', 'i2_needgood_list.txt')) if l.strip()]
res = [json.loads(l) for l in open(os.path.join(RI, 'i2', 'i2_residual.jsonl'))]
rows = [json.loads(l) for l in open(os.path.join(RI, 'i8', 'i8_first360.jsonl'))]
say('records %d; ratios = lines 1-360 of the list in order: %s; sources of those lines: %s'
    % (len(rows), [(o['m'], o['n']) for o in rows] == lst[:360], dict(collections.Counter(src[:360]))))
same = all({tuple(t['tau']): t['kill'] for t in o['tuples']} == {tuple(t): k for t, k in zip(r['residual'], r['kill'])} for o, r in zip(rows, res[:360]))
say('tuples and killing primes = those of i2_residual.jsonl lines 1-360: %s; tuples there: %d' % (same, sum(len(r['residual']) for r in res[:360])))
MIX = {'0L', '0U', '8L', '8U'}
ver = collections.Counter(); per = collections.Counter(); mixonly = 0; incons = 0; labels = collections.Counter(); neither = []
for o in rows:
    for t in o['tuples']:
        kb = sorted(jk for jk, v in t['pairs'].items() if v in ('killed', 'outside Sel2'))
        if kb != sorted(t['killed_by']): incons += 1
        v = 'pairing' if kb else ('undecided' if 'undecided' in t['pairs'].values() else 'neither')
        if v != 'pairing' and any(x.get('excludes') for x in o['rank0'].values()): v = 'rank0'
        if (v == t['verdict']) is False and not (v == 'rank0' and t['verdict'].startswith('rank 0')): incons += 1
        ver[v] += 1
        for jk in kb: per[jk] += 1
        if kb and set(kb) <= MIX: mixonly += 1
        for jk, x in t['pairs'].items(): labels[x] += 1
        if v == 'neither': neither.append((o['m'], o['n'], t['tau'], t['kill'], ' '.join('%s:%s' % kv for kv in sorted(t['pairs'].items()))))
say('verdicts recomputed from the labels: %s; records whose killed_by or verdict disagree with the labels: %d' % (dict(sorted(ver.items())), incons))
say('kills per curve: %s; killed only on mixed curves: %d' % (' '.join('%s:%d' % (jk, per[jk]) for jk in ('08', '0L', '0U', '8L', '8U', 'LU')), mixonly))
say('per-curve labels over all tuples: %s' % dict(sorted(labels.items())))
for x in neither: say('  neither: (%d,%d) %s kill %s; %s' % x)
r0 = collections.Counter(); r0rat = []
for o in rows:
    ex = sorted(jk for jk, x in o['rank0'].items() if x.get('excludes'))
    for jk, x in o['rank0'].items(): r0[(jk, x['excludes'], x['complete'])] += 1
    if ex: r0rat.append('(%d,%d) by %s (torsion %s)' % (o['m'], o['n'], ','.join(ex), ','.join('%d/%d' % (o['rank0'][j]['torsion'], o['rank0'][j]['bound']) for j in ex)))
say('pair curves with dim Sel2 = 2 (rank0 entries), by (curve, excludes, complete): %s' % dict(sorted(r0.items())))
say('ratios excluded by a rank-0 pair curve: %d' % len(r0rat))
for x in r0rat: say('  ' + x)
nrank0_not_ex = [(o['m'], o['n'], jk, x['t']) for o in rows for jk, x in o['rank0'].items() if not x['excludes']]
for x in nrank0_not_ex: say('  dim Sel2 = 2 but not excluding: (%d,%d) E_%s points over t in %s' % x)
s = '\n'.join(OUT) + '\n'; f = os.path.join(HERE, 'z17_i8count.txt')
if '--write' in sys.argv: open(f, 'w').write(s)
else:
    sys.stdout.write(s)
    if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS')
