#!/usr/bin/env python3
"""Review z7 of roadmap H, H7 claim 4: the curves E04 of tier B (scripts/g1_full/dm_evidence.txt), independently.
Written by an AI model (review agent z7, 27 November 2026); not read by a person.  Experiment, not counted.

For each E04 line of dm_evidence.txt: the curve is dm.py's [0, 2S^2, 0, 4S^2m^2n^2, 0] at the pair the line was run on
(field 2), = y^2 = x(x + 2Sm^2)(x + 2Sn^2); the line's generator abscissa G is checked to be the x of a point (exact);
then my own Sel^2 (my2desc.py) and Cassels-Tate pairing (ct_mine.py) on E04 itself, with the images of T1, T2 and G
checked to lie in the kernel.  Each curve is compared with the first curve of its line in H7's h7_e04_rank.txt (the
bundle's result on E04 itself; for (320,289) the bundle's pairing search gave up there).  Counts at the end.
Rests on: as ct_mine.py.  PARI: inside ct_mine.py only.

    python3 run_ct_e04.py            prints run_ct_e04.txt   (about 3 minutes)
"""
import sys, os, re
from fractions import Fraction as Fr
from math import isqrt
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
REPO = os.environ.get('MSS_REPO', os.path.join(HERE, '..', '..'))
H7 = os.environ.get('H7_DIR', os.path.join(HERE, '..', 'roadmap-h', 'h7'))
import ct_mine as C

def evidence():
    ev = []
    for L in open(os.path.join(REPO, 'scripts', 'g1_full', 'dm_evidence.txt')):
        if L.startswith('#') or not L.strip(): continue
        w = L.split()
        if w[2] != 'E04': continue
        f = dict(x.split('=', 1) for x in w[3:] if '=' in x)
        ev.append((tuple(map(int, w[0].split(','))), tuple(map(int, w[1].split(','))), Fr(f['G'])))
    return ev

def bundle_first():
    b = {}
    for l in open(os.path.join(H7, 'h7_e04_rank.txt')):
        mm = re.match(r'(\d+),(\d+) E04: .*?curve\(s\): (.*?); so', l)
        if mm:
            first = mm.group(3).split(';')[0]
            s = re.search(r'dim Sel\^2 = (\d), Cassels-Tate pairing of rank (\d)', first)
            b[(int(mm.group(1)), int(mm.group(2)))] = (int(s.group(1)), int(s.group(2))) if s else None
    return b

def main():
    print('# run_ct_e04.py (review z7; written by an AI model): my complete 2-descent and Cassels-Tate pairing on every E04')
    print('# curve of dm_evidence.txt, compared with the first curve of each line of H7\'s h7_e04_rank.txt')
    ev = evidence(); runs = []
    for _, rp, gx in ev:
        if (rp, gx) not in runs: runs.append((rp, gx))
    print('E04 lines %d, curves (run pair, G) %d, distinct run pairs %d' % (len(ev), len(runs), len(set(r for r, _ in runs))))
    b = bundle_first(); tally = Counter(); differ = []
    for (m, n), gx in runs:
        S = m * m + n * n
        es = tuple(sorted((0, -2 * S * m * m, -2 * S * n * n)))
        v = (gx - es[0]) * (gx - es[1]) * (gx - es[2])
        assert v > 0 and isqrt(v.numerator) ** 2 == v.numerator and isqrt(v.denominator) ** 2 == v.denominator, 'G not on E04'
        res = C.analyse(es, seed=5, xs=[gx], label='%d,%d:' % (m, n), show_matrix=False)
        tally[(res['dim'], res['rho'], res['bound'])] += 1
        bb = b.get((m, n))
        if bb != (res['dim'], res['rho']):
            differ.append((m, n)); print('   H7 first curve: %s' % ('pairing not computed' if bb is None else 'dim %d, rank %d' % bb))
    print('tally on E04 itself, (dim Sel^2, pairing rank, bound): %s' % ', '.join('%s x%d' % (k, v) for k, v in sorted(tally.items())))
    print('curves where H7\'s first curve differs from mine: %s' % differ)
    print('# end')

if __name__ == '__main__':
    main()
