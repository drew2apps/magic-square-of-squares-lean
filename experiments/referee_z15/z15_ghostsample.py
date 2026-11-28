"""z15: ghost-shape ratios for the t10 rerun: the ratios (not the pilot's 201, not in sample80) at which the one E_08
class of the S1 tuples or the one E_LU class of the S2 tuples survives that curve (so the mixed curves decide), in an
order fixed by seed 7; then the other ghost ratios, seed 8.  Written by an AI model (reviewer z15), not read by a person."""
import os, json, random, re
import cypari2
P = cypari2.Pari(); HERE = os.path.dirname(os.path.abspath(__file__))
def sqf(x): return int(P.core(abs(int(x))))
CV = ['08', 'LU', '0L', '0U', '8L', '8U']
rows = [json.loads(l) for f in ('out_997.jsonl', 'out_15000_0.jsonl', 'out_15000_1.jsonl') for l in open(os.path.join(HERE, 'session_data', f))]
pilot = set()
for f in ('t10_997_a.txt', 't10_997_b.txt', 't10_15000_s100_a.txt', 't10_15000_s100_b.txt'):
    for l in open(os.path.join(HERE, '..', 'conductor-t10', f)):
        mm = re.match(r'ratio \((\d+),(\d+)\)', l)
        if mm: pilot.add((int(mm.group(1)), int(mm.group(2))))
s80 = set(tuple(map(int, l.split(','))) for l in open(os.path.join(HERE, 'sample80_seed20261127.txt')))
hot, cold = [], []
for r in rows:
    m, n = r['m'], r['n']
    if (m, n) in pilot or (m, n) in s80: continue
    A = m * m - n * n; isg = False; hotr = False
    for t, d in r['tuples']:
        a = sqf(m * n * t[0] * t[1]); b = sqf(A * t[2] * t[3])
        if (a, b) == (sqf(2 * m * n), 1): isg = True; hotr |= d['08'] == 0
        if (a, b) == (1, sqf(2 * A)): isg = True; hotr |= d['LU'] == 0
    if isg: (hot if hotr else cold).append((m, n))
hot.sort(); cold.sort(); random.Random(7).shuffle(hot); random.Random(8).shuffle(cold)
open(os.path.join(HERE, 'ghost_order.txt'), 'w').write(''.join('%d,%d\n' % r for r in hot + cold))
print('ghost ratios outside the pilot and sample80: %d with a surviving ghost class on E_08/E_LU (first), %d others' % (len(hot), len(cold)))
