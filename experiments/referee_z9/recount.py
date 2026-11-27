"""Review z9 of roadmap H (run 7), written by an AI model; not read by a person.  Experiment, not counted.
Checks H2's counts from its committed outputs, with a parser written here (no H2 code):
  * the exhaustive run h2_ex_15000_{0,1,2}of3.txt: every m with 15000 < m <= 30000 exactly once, the coprime count of
    each m equal to phi(m), the per-m counts summing up, the per-block table of REPORT_H2.md, the totals, the largest P;
  * the biased hunt h2_hunt_k*.txt and the excess-guided hunt h2_design_k*.txt: class-keeping LINES against DISTINCT
    ratios, the ratios with m <= 30000 (inside G2's range or the exhaustive one), the count with m > 30000 (2,131),
    the timeouts, the ranges reached; the hunt's rows inside 15000 < m <= 30000 against the exhaustive output, and its
    rows with m <= 15000 against G2's class_selmer_15000_list.txt;
  * the core-seconds of the exhaustive run in H2's times.txt.
Rests on: the committed outputs only (it reruns nothing).  usage: python3 -I recount.py [--write]   (recount.txt)"""
import os, sys, re, glob, collections
from math import gcd, log2
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import *
O = Out('recount.txt')
O('Review z9: recount of H2 (experiments/roadmap-h/h2). Written by an AI model; not read by a person.')

def phi(m):
    r = m; x = m; p = 2
    while p * p <= x:
        if x % p == 0:
            while x % p == 0: x //= p
            r -= r // p
        p += 1
    if x > 1: r -= r // x
    return r

ms, rows = {}, {}
for i in range(3):
    for l in open(os.path.join(H2, 'h2_ex_15000_%dof3.txt' % i)):
        l = l.strip()
        if not l: continue
        if l.startswith('m='):
            w = l.split(); m = int(w[0][2:]); assert m not in ms, ('m twice', m)
            ms[m] = {k: int(v) for k, v in (x.split('=') for x in w[1:])}
        else:
            a, dim, ls, P = l.split(); m, n = map(int, a.split(',')); assert (m, n) not in rows
            rows[(m, n)] = (int(dim), int(ls), P)
O('EXHAUSTIVE: values of m %d (15001..30000 missing: %d); coprime = phi(m) for every m: %s; per-m sums consistent: %s' % (
    len(ms), sum(1 for m in range(15001, 30001) if m not in ms), all(ms[m]['coprime'] == phi(m) for m in ms),
    all(d['coprime'] == d['degenerate'] + d['no_point'] + d['no_class'] + d['class'] for d in ms.values())))
bym = collections.Counter(m for m, n in rows)
O('class lines = class counts for every m: %s; rows coprime with 0 < n < m: %s' % (all(ms[m]['class'] == bym[m] for m in ms), all(gcd(m, n) == 1 and 0 < n < m for m, n in rows)))
tot = collections.Counter()
for d in ms.values(): tot.update(d)
O('totals: coprime %d, degenerate %d, not degenerate %d, no point %d, no class %d, class %d' % (
    tot['coprime'], tot['degenerate'], tot['coprime'] - tot['degenerate'], tot['no_point'], tot['no_class'], tot['class']))
Ps = [r[2] for r in rows.values()]; num = [int(p) for p in Ps if p != 'SURVIVOR']
O('survivors %d; need a good prime %d; largest P %d at %s' % (Ps.count('SURVIVOR'), sum(1 for p in num if p), max(num), [r for r, v in rows.items() if v[2] == str(max(num))]))
O('P values (count): %s' % ' '.join('%d:%d' % kv for kv in sorted(collections.Counter(num).items())))
O('block       coprime  degenerate  no point  no class  class  need a good prime  largest P')
blk = collections.defaultdict(collections.Counter); bP = collections.defaultdict(list)
for m, d in ms.items(): blk[(m - 1) // 1000 * 1000].update(d)
for (m, n), v in rows.items(): bP[(m - 1) // 1000 * 1000].append(int(v[2]))
for b in sorted(blk):
    c = blk[b]
    O('%5d-%5d %9d %6d %10d %8d %6d %6d %6d' % (b, b + 1000, c['coprime'], c['degenerate'], c['no_point'], c['no_class'], c['class'], sum(1 for p in bP[b] if p), max(bP[b])))
O('largest P per block: %s (range %d to %d; not monotone)' % (', '.join(str(max(bP[b])) for b in sorted(bP)), min(max(v) for v in bP.values()), max(max(v) for v in bP.values())))

# hunts
lines, heads = [], []
for f in sorted(glob.glob(os.path.join(H2, 'h2_hunt_k*.txt'))):
    L = open(f).read().splitlines()
    mm = re.match(r'sampled (\d+) \(draws \d+\); no local point (\d+); no class (\d+); class kept (\d+); timed out: class (\d+), descent (\d+)', L[1])
    s, npt, nc, ck, tc, td = map(int, mm.groups()); k = int(re.search(r'k(\d+)', os.path.basename(f)).group(1))
    rs = [l for l in L[3:] if re.match(r'\d+,\d+ class dim', l)]; assert len(rs) == ck and s == npt + nc + ck + tc
    heads.append((os.path.basename(f), k, s, ck, tc, td))
    for l in rs:
        w = l.split(); lines.append((tuple(map(int, w[0].split(','))), int(w[3]), int(w[4]), w[5], os.path.basename(f)))
hunt = {}
for r, d, ls, P, f in lines: hunt.setdefault(r, (d, ls, P))
dups = collections.Counter(r for r, *_ in lines)
O('HUNT: files %d, sampled %d, class-keeping lines %d, DISTINCT ratios %d (%d ratios occur more than once, %d repeated lines: sampling with replacement)' % (
    len(heads), sum(h[2] for h in heads), len(lines), len(hunt), sum(1 for v in dups.values() if v > 1), sum(v - 1 for v in dups.values())))
O('  repeated: %s' % ' '.join('(%d,%d)x%d' % (r + (v,)) for r, v in sorted(dups.items()) if v > 1))
O('  every repeated line has the same P: %s' % all(len(set(P for rr, d, ls, P, f in lines if rr == r)) == 1 for r in dups))
O('  timeouts: %s' % ', '.join('%s class %d descent %d' % (h[0], h[4], h[5]) for h in heads if h[4] or h[5]))
O('  ranges run: k = %s; largest m of a class-keeping ratio %d (k = %d); survivors %d' % (
    ', '.join(sorted(set(str(h[1]) for h in heads), key=int)), max(r[0] for r in hunt), int(log2(max(r[0] for r in hunt))), sum(1 for v in hunt.values() if v[2] == 'SURVIVOR')))
low = [r for r in hunt if r[0] <= 15000]; mid = [r for r in hunt if 15000 < r[0] <= 30000]; high = [r for r in hunt if r[0] > 30000]
g2l = set(tuple(map(int, l.split()[0].split(','))) for l in open(os.path.join(STRUCT, 'class_selmer_15000_list.txt')) if l.strip())
O('  distinct with m <= 15000: %d (all in class_selmer_15000_list.txt: %s); 15000 < m <= 30000: %d (all equal to the exhaustive row, dim, log2, P: %s); m > 30000: %d' % (
    len(low), all(r in g2l for r in low), len(mid), all(rows.get(r) == hunt[r] for r in mid), len(high)))
O('  m > 30000: largest P %d at %s; need a good prime %d' % (max(int(hunt[r][2]) for r in high), [r for r in high if int(hunt[r][2]) == max(int(hunt[q][2]) for q in high)], sum(1 for r in high if int(hunt[r][2]))))
dz = []; dt = []; taken = 0
for f in sorted(glob.glob(os.path.join(H2, 'h2_design_k*.txt'))):
    L = open(f).read().splitlines()
    for l in L[3:]:
        if not re.match(r'\d+,\d+ ', l): continue
        taken += 1; w = l.split()
        if ' class dim ' in l:
            (dz if w[-1].isdigit() else dt).append((tuple(map(int, w[0].split(','))), w[-1], l))
O('DESIGN: taken %d, class-keeping %d (%.2f%% against the hunt %.2f%%), descent finished %d, unfinished %s' % (
    taken, len(dz) + len(dt), 100.0 * (len(dz) + len(dt)) / taken, 100.0 * len(lines) / sum(h[2] for h in heads), len(dz), [x[2] for x in dt]))
O('  finished: m from %d to %d; need a good prime %d (largest P %d); in the hunt too: %d' % (
    min(r[0] for r, _, _ in dz), max(r[0] for r, _, _ in dz), sum(1 for _, P, _ in dz if int(P)), max(int(P) for _, P, _ in dz), sum(1 for r, _, _ in dz if r in hunt)))
O('DISTINCT class-keeping ratios of the two hunts: %d; with m > 30000: %d + %d = %d (not 2,421)' % (len(hunt) + len(dz), len(high), len(dz), len(high) + len(dz)))
cs = [int(re.search(r'exit=0 (\d+) s', l).group(1)) for l in open(os.path.join(H2, 'times.txt')) if l.startswith('h2_exhaust 15000 30000') and 'exit=0' in l]
O('exhaustive run, times.txt: %d workers, %d core-seconds, %.2f core-hours per block of 1000 on average' % (len(cs), sum(cs), sum(cs) / 3600.0 / 15))
O.finish()
