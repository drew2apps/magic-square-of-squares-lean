#!/usr/bin/env python3
"""z11: t10's counts recounted from its committed per-run outputs with an own parser, against G2's lists.
Written by an AI model (review agent z11); not read by a person.  No PARI.
    python3 -I recount.py      (prints recount.txt)"""
import sys, os, re, random, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _z11common import T10, G2, Out
say = Out('recount.txt')
FILES = ['t10_997_a.txt', 't10_997_b.txt', 't10_15000_s100_a.txt', 't10_15000_s100_b.txt']
IDX = {'0': 0, '8': 1, 'L': 2, 'U': 3}; PAIRS = ['08', '0L', '0U', '8L', '8U', 'LU']
ratios = {}; heads = {}; per_file = collections.Counter(); perpair = collections.Counter(); firstp = collections.Counter()
kill_lines = 0; nterms = collections.Counter(); pk = 0; mism = 0; kills = set(); prods = set(); fp_ratio = collections.defaultdict(list)
dim2 = collections.defaultdict(list); verdicts = collections.Counter()
for f in FILES:
    ratio = cur = None
    for l in open(os.path.join(T10, f)):
        m = re.match(r'ratio \((\d+),(\d+)\): (\d+) residual tuple\(s\) after the bad primes; known tuples (\d+)', l)
        if m:
            ratio = (int(m.group(1)), int(m.group(2))); ratios[ratio] = f; heads[ratio] = int(m.group(3)); assert m.group(4) == '0'; continue
        m = re.match(r'  E_(\w\w): .*dim Sel2 = (\d+), point images (\d+)', l)
        if m:
            cur = m.group(1)
            if int(m.group(2)) == 2: dim2[ratio].append(cur)
            continue
        m = re.match(r'    tuple (\(.*?\)): class (\(.*?\)): (killed|outside Sel2) .*local -1 at \[(.*)\]', l)
        if m:
            kill_lines += 1; nterms[len(m.group(4).split(','))] += 1; kills.add((ratio, cur, m.group(2))); continue
        m = re.match(r'    (\([0-9, ]+\))  good prime (\S+)  (.*)  => (.*)$', l.rstrip('\n'))
        if m:
            t, g, prs, v = m.groups(); per_file[f] += 1
            ks = [x.split(':')[0] for x in prs.split() if x.split(':')[1] in ('killed', 'outside')]
            for k in ks: perpair[k] += 1
            for x in prs.split(): verdicts[x.split(':')[1]] += 1
            firstp[int(g)] += 1; fp_ratio[ratio].append(int(g)); pk += bool(ks)
            mism += (('PAIRING-KILLED by ' + ','.join(ks)) if ks else 'not pairing-killed') != v
            tt = eval(t)
            for jk in PAIRS: prods.add((ratio, jk, tt[IDX[jk[0]]] * tt[IDX[jk[1]]]))
tot = sum(per_file.values())
say('recount.py (review z11; written by an AI model): t10\'s outputs %s, own parser' % ', '.join(FILES))
say('  ratios %d (%s); residual tuples %d = %s; sum of the header counts %d'
    % (len(ratios), ', '.join('%s %d' % (f, list(ratios.values()).count(f)) for f in FILES), tot, ' + '.join(str(per_file[f]) for f in FILES), sum(heads.values())))
say('  the 101 of g2_997_all.txt (997_a + 997_b): %d tuples; the 15000 sample (_a + _b): %d tuples' % (per_file[FILES[0]] + per_file[FILES[1]], per_file[FILES[2]] + per_file[FILES[3]]))
say('  pairing-killed %d, not %d; verdict text consistent with the pairs on every line: %s' % (pk, tot - pk, mism == 0))
say('  kills per curve: %s' % ' '.join('%s:%d' % (jk, perpair[jk]) for jk in PAIRS))
say('  per (tuple, curve) verdicts: %s' % ' '.join('%s:%d' % kv for kv in sorted(verdicts.items())))
say('  first killing good prime: %s' % ' '.join('%d:%d' % kv for kv in sorted(firstp.items())))
say('  kill lines (one per killed tuple and curve): %d = the sum of the per-curve counts: %s' % (kill_lines, kill_lines == sum(perpair.values())))
say('  local -1 terms on them: %d (terms per line: %s; all odd: %s)' % (sum(k * v for k, v in nterms.items()),
    ' '.join('%d x %d' % (v, k) for k, v in sorted(nterms.items())), all(k % 2 for k in nterms)))
say('  distinct (ratio, curve, class) kills: %d; distinct (ratio, curve, d_j d_k) tested: %d of %d tuple-curve tests' % (len(kills), len(prods), 6 * tot))
# G2's lists
def g2(path):
    out = {}; order = []
    for l in open(path):
        if not l.strip() or l[0] == '#' or 'left' not in l: continue
        r = tuple(map(int, l.split()[0].split(','))); order.append(r)
        out[r] = (int(re.search(r'bad-alone: (\d+) left', l).group(1)),
                  max(int(x.split(':')[0]) for x in re.search(r'kill: (.*)$', l).group(1).split() if x.endswith(':g')))
    return out, order
a, oa = g2(os.path.join(G2, 'g2_997_all.txt')); b, ob = g2(os.path.join(G2, 'g2_15000_need_good.txt'))
say('  G2: %d ratios of g2_997_all.txt need a good prime (m from %d to %d), %d of g2_15000_need_good.txt (m from %d to %d)'
    % (len(a), min(r[0] for r in a), max(r[0] for r in a), len(b), min(r[0] for r in b), max(r[0] for r in b)))
smp = sorted(random.Random(1).sample(ob, 100))
said = {r for r, f in ratios.items()}
say('  t10\'s 997 runs = all %d in file order, 0..79 in 997_a: %s; the 15000 runs = random.Random(1).sample of 100, first 50 in _a: %s'
    % (len(a), set(oa) == {r for r, f in ratios.items() if '997' in f} and set(oa[:79]) == {r for r, f in ratios.items() if f == FILES[0]},
       set(smp) == {r for r, f in ratios.items() if '15000' in f} and set(smp[:50]) == {r for r, f in ratios.items() if f == FILES[2]}))
say('  residual counts equal G2\'s "N left" for every ratio: %s' % all(heads[r] == {**a, **b}[r][0] for r in ratios))
say('  largest first killing prime equals the largest good prime of G2\'s kill set: %d of %d ratios (above it: %d)'
    % (sum(max(v) == {**a, **b}[r][1] for r, v in fp_ratio.items()), len(fp_ratio), sum(max(v) > {**a, **b}[r][1] for r, v in fp_ratio.items())))
say('  ratios with a pair curve of dim Sel2 = 2 in t10\'s lines: %s' % ', '.join('(%d,%d) %s' % (r + (','.join(v),)) for r, v in sorted(dim2.items())))
say.finish()
