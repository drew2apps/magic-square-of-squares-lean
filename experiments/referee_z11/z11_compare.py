#!/usr/bin/env python3
"""z11: my code (z11_pair.py) against t10's committed per-run outputs, (tuple, pair curve) by (tuple, pair curve):
the verdict (killed or not, both signs of the class), and, where t10 prints them, the class and the killing b.
Written by an AI model (review agent z11); not read by a person.  PARI as in z11_pair.py (no ellrank here).
    python3 -I z11_compare.py                 4 ratios: (465,416), (5623,1017), (3320,2303), (14885,6076)  (z11_compare.txt)
    python3 -I z11_compare.py --set mid       (2999,1081), (3073,412), (3143,2567), from the middle of t10_997_a.txt
    python3 -I z11_compare.py --set batch30   30 ratios with <= 40 tuples, random.Random(11).sample  (z11_compare_batch30.txt)
    python3 -I z11_compare.py --set big      (9751,5084), 216 tuples  (z11_compare_big.txt)"""
import sys, os, re, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _z11common import T10, Out
from z11_pair import analyse, PAIRS
FILES = ['t10_997_a.txt', 't10_997_b.txt', 't10_15000_s100_a.txt', 't10_15000_s100_b.txt']
def t10_blocks():
    out = {}; heads = []
    for f in FILES:
        cur = None
        for l in open(os.path.join(T10, f)):
            l = l.rstrip('\n')
            m = re.match(r'ratio \((\d+),(\d+)\): (\d+) residual', l)
            if m: cur = (int(m.group(1)), int(m.group(2))); out[cur] = (f, []); heads.append((cur, int(m.group(3)))); continue
            if l.startswith('SUMMARY') or l == '': cur = None
            if cur: out[cur][1].append(l)
    return out, heads
def parse(blk):
    cls = {}; ver = {}; cur = None
    for l in blk:
        m = re.match(r'  E_(\w\w): ', l)
        if m: cur = m.group(1); continue
        m = re.match(r'    tuple (\(.*?\)): class (\(.*?\)): (killed|outside Sel2) b = (\(.*?\))', l)
        if m: cls[(eval(m.group(1)), cur)] = (eval(m.group(2)), eval(m.group(4))); continue
        m = re.match(r'    (\([0-9, ]+\))  good prime (\S+)  (.*)  => (.*)$', l)
        if m:
            for x in m.group(3).split(): jk, v = x.split(':'); ver[(eval(m.group(1)), jk)] = v
    return cls, ver
B, HEADS = t10_blocks()
SETS = {'small': [(465, 416), (5623, 1017), (3320, 2303), (14885, 6076)], 'mid': [(2999, 1081), (3073, 412), (3143, 2567)],
        'big': [(9751, 5084)]}
done = set(SETS['small'] + SETS['mid'] + SETS['big'])
SETS['batch30'] = [r for r, c in sorted(random.Random(11).sample([(r, c) for r, c in HEADS if r not in done and c <= 40], 30))]
name = sys.argv[sys.argv.index('--set') + 1] if '--set' in sys.argv else 'small'
say = Out('z11_compare.txt' if name == 'small' else 'z11_compare_%s.txt' % name)
say('z11_compare.py --set %s (review z11; written by an AI model): my code against t10\'s committed lines' % name)
T = dict(agree=0, dis=0, cls=0, clsdiff=0, bdiff=0, signs=0)
for m, n in SETS[name]:
    f, blk = B[(m, n)]; cls, ver = parse(blk)
    out, kill, known = analyse(m, n, points=False, log=lambda s: None)
    a = d = cd = bd = sd = nc = 0
    for t in out:
        for jk in PAIRS:
            (cp, rp), (cm, rm) = out[t][jk]
            mk = rp[0] in ('killed', 'outside'); tk = ver[(t, jk)] in ('killed', 'outside')
            sd += (rp[0] in ('killed', 'outside')) != (rm[0] in ('killed', 'outside'))
            if mk == tk: a += 1
            else: d += 1
            if (t, jk) in cls:
                nc += 1; tc, tb = cls[(t, jk)]
                cd += tuple(cp) != tc; bd += (rp[0] == 'killed' and tuple(rp[1]) != tb)
    say('(%d,%d) [%s]: %d tuples; verdicts agree %d, disagree %d; t10 classes compared %d, differ %d; killing b differ %d; signs disagree %d'
        % (m, n, f, len(out), a, d, nc, cd, bd, sd))
    T['agree'] += a; T['dis'] += d; T['cls'] += nc; T['clsdiff'] += cd; T['bdiff'] += bd; T['signs'] += sd
say('TOTAL: %d (tuple, curve) verdicts agree, %d disagree; %d classes compared, %d differ; %d killing b differ; both signs disagree on %d'
    % (T['agree'], T['dis'], T['cls'], T['clsdiff'], T['bdiff'], T['signs']))
say.finish()
