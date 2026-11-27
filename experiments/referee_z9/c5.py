"""Review z9 of roadmap H (run 7), written by an AI model; not read by a person.  Experiment, not counted.
The genus-1 quotient y^2 = L(t) U(t) of (627,364) (H6's class c5; h6_split's forms, so 4L, 4U of G1_FULL_SYSTEM.md):
  * its Jacobian (PARI's ellfromeqn, then ellminimalmodel) equals that of y^2 = F0 F8 of the partner (991,263)
    (H6's c1 there), and that of y^2 = F0 F8 of (627,364) is printed beside it;
  * torsion subgroup, root number, conductor (PARI: elltors, ellrootno, ellglobalred; mathematical invariants);
  * the rational points of the quartic at t = 0, 1, -1 and at infinity (exact squares), 8 in all;
  * from H6's committed h6_ranks_<m>_<n>.txt: the classes whose PARI lower bound exceeds the number of points PARI
    returned on every curve of the class (so the lower bound rests on parity, i.e. on Sha being finite).
What it shows: with torsion of order 8 and the 8 obvious points, rank 0 is not excluded by counting points; no point of
infinite order is known on c5 of (627,364) (= c1 of (991,263)); root number -1 makes rank 1 the parity prediction.
Rests on: PARI 2.17.4 for the models, torsion, root number and conductor.  No ellrank is called (its output varies).
usage: python3 -I c5.py [--write]   (c5.txt; a few seconds)"""
import os, sys, re
from math import isqrt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import *
O = Out('c5.txt')
import cypari2
pari = cypari2.Pari(); pari.allocatemem(10 ** 9)
sys.path.insert(0, H6)
import h6_split
O('Review z9: the quotient y^2 = L*U of (627,364) and y^2 = F0*F8 of (991,263). Written by an AI model; not read by a person.')
def jac(m, n, A, B):
    F = h6_split.forms(m, n); a, b = F[A], F[B]
    q = '(%d + %d*x + %d*x^2)*(%d + %d*x + %d*x^2)' % tuple(a + b)
    E = pari.ellminimalmodel(pari.ellinit(pari('ellfromeqn(y^2 - %s)' % q)))
    pts = [t for t in (0, 1, -1) if isqrt(v := (a[0] + a[1]*t + a[2]*t*t) * (b[0] + b[1]*t + b[2]*t*t)) ** 2 == v and v > 0]
    inf = isqrt(a[2] * b[2]) ** 2 == a[2] * b[2]
    return E, pts, inf
res = {}
for (m, n, A, B) in ((627, 364, 'L', 'U'), (991, 263, 'F0', 'F8'), (627, 364, 'F0', 'F8')):
    E, pts, inf = jac(m, n, A, B); res[(m, n, A)] = [int(x) for x in E[:5]]
    tors = pari.elltors(E)
    O('(%d,%d) y^2 = %s*%s: minimal model of the Jacobian %s' % (m, n, A, B, res[(m, n, A)]))
    O('    torsion %s (order %d); root number %d; conductor %d' % ([int(x) for x in tors[1]], int(tors[0]), int(pari.ellrootno(E)), int(pari.ellglobalred(E)[0])))
    O('    rational points of the quartic: t = %s with y = +-, and %s at infinity: %d points' % (pts, 2 if inf else 0, 2 * len(pts) + (2 if inf else 0)))
O('Jacobian of L*U (627,364) = Jacobian of F0*F8 (991,263): %s' % (res[(627, 364, 'L')] == res[(991, 263, 'F0')]))
for mn in ('627_364', '991_263'):
    cur = None; info = {}; maps = {}
    for l in open(os.path.join(H6, 'h6_ranks_%s.txt' % mn)):
        mm = re.match(r'(\S+)\s+.*?genus 1\s+(c\d+)\s+rank', l)
        if mm: maps.setdefault(mm.group(2), []).append(mm.group(1)); continue
        mm = re.match(r'(c\d+)\s+(.*?)\s+minimal model', l)
        if mm: cur = mm.group(1); info[cur] = [mm.group(2), []]; continue
        mm = re.match(r'\s+\[.*\]\s+\[(\d+), (\S+), (\d+), (\d+)', l)
        if mm and cur: info[cur][1].append((int(mm.group(1)), int(mm.group(3))))
    par = [(c, v[0], len(maps[c]), max(p for _, p in v[1])) for c, v in sorted(info.items(), key=lambda kv: int(kv[0][1:])) if max(p for _, p in v[1]) < max(lo for lo, _ in v[1])]
    O('(%s) classes whose PARI lower bound exceeds the points returned: %s' % (mn.replace('_', ','), '; '.join('%s %s, %d maps, at most %d points on one curve' % x for x in par)))
O.finish()
