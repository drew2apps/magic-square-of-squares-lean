"""Referee z17 (written by an AI model; not read by a person): conductor-t10's analysis from scratch
(z17_i8rerun.jsonl, z17_i8rerun.py, PARI of this container) against I8's committed records (i8_first360.jsonl).
Per ratio: the residual tuples and killing primes (t10's own residual(), not I2's file), and per tuple the set of
pair curves that kill it (basis-independent) and the per-curve labels ('image' vs 'kernel' depends on the points PARI
finds; reported separately); per curve dim Sel^2 (so whether I8's rank-0 test could apply).
usage: python3 -I z17_i8compare.py [--write]   (z17_i8compare.txt)"""
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__))
ref = {(o['m'], o['n']): o for o in map(json.loads, open(os.path.join(HERE, '..', 'roadmap-i', 'i8', 'i8_first360.jsonl')))}
mine = [json.loads(l) for l in open(os.path.join(HERE, 'z17_i8rerun.jsonl')) if l.startswith('{')]
order = [tuple(map(int, l.split(','))) for l in open(os.path.join(HERE, 'z17_i8_sample.txt')) if l.strip()]
mine = sorted(mine, key=lambda o: order.index((o['m'], o['n'])))
OUT = ['z17: t10 from scratch on %d of I8\'s 360 ratios (the 3 "neither" ratios, then 15 drawn with seed 17) against i8_first360.jsonl. Written by an AI model (reviewer z17); not read by a person.' % len(mine)]
ok_t = ok_k = ok_l = 0; nt = 0
for o in mine:
    r = ref[(o['m'], o['n'])]
    a = {tuple(t['tau']): t for t in o['tuples']}; b = {tuple(t['tau']): t for t in r['tuples']}
    same_t = {k: v['kill'] for k, v in a.items()} == {k: v['kill'] for k, v in b.items()}
    same_k = same_t and all(sorted(a[k]['killed_by']) == sorted(b[k]['killed_by']) for k in a)
    same_l = same_t and all(a[k]['pairs'] == b[k]['pairs'] for k in a)
    ok_t += same_t; ok_k += same_k; ok_l += same_l; nt += len(a)
    neither = [k for k in a if not a[k]['killed_by']]
    dims = ' '.join('%s:%d' % (jk, o['curves'][jk]['dimsel']) for jk in ('08', '0L', '0U', '8L', '8U', 'LU'))
    OUT.append('(%d,%d): %d tuples (known %d); tuples and killing primes %s; killing curves %s; labels %s; not killed: %d; dim Sel2 %s'
               % (o['m'], o['n'], len(a), o['known'], 'same' if same_t else 'DIFFER', 'same' if same_k else 'DIFFER',
                  'same' if same_l else 'differ', len(neither), dims))
    for k in neither:
        OUT.append('    not killed: %s kill %s; %s' % (list(k), a[k]['kill'], ' '.join('%s:%s' % kv for kv in sorted(a[k]['pairs'].items()))))
    if same_t and not same_l:
        for k in a:
            if a[k]['pairs'] != b[k]['pairs']:
                OUT.append('    labels differ at %s: here %s; committed %s' % (list(k), a[k]['pairs'], b[k]['pairs']))
OUT.append('ratios: %d; tuples %d; same tuples and killing primes %d; same killing curves (the verdicts) %d; same labels %d' % (len(mine), nt, ok_t, ok_k, ok_l))
s = '\n'.join(OUT) + '\n'; f = os.path.join(HERE, 'z17_i8compare.txt')
if '--write' in sys.argv: open(f, 'w').write(s)
else:
    sys.stdout.write(s)
    if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS')
