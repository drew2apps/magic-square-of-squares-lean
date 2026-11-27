#!/usr/bin/env python3
"""z11: try to make t10's committed method (../conductor-t10/t10_pairct.py, imported, not changed) kill the class of
a rational point (it must never).  t10's PairData is built WITHOUT PARI's points (only the 2-torsion and psi(iota t0)),
so that the classes of PARI's points (all subset sums) and the known tuples of degenerate ratios are no longer
"image" and go through t10's pairing.  Written by an AI model (review agent z11); not read by a person.
PARI: t10's code as it is; its point finder (ellrank, elltors) with setrand(1) before each call; points checked.
    python3 -I forge_t10.py      (prints forge_t10.txt)"""
import sys, os, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _z11common import T10, Out
from fractions import Fraction as Fr
argv = sys.argv; sys.argv = [argv[0]]
sys.path.insert(0, T10)
fd = os.dup(1); nul = os.open(os.devnull, os.O_WRONLY); os.dup2(nul, 1)
try: import t10_pairct as T                     # its PARI stack message goes to /dev/null
finally: sys.stdout.flush(); os.dup2(fd, 1); os.close(nul); os.close(fd)
sys.argv = argv
say = Out('forge_t10.txt')
say('forge_t10.py (review z11; written by an AI model): t10_pairct.PairData without PARI\'s points')
real_find = T.find_points
def find(Pr):
    T.PARI.setrand(1); return real_find(Pr)
allok = True
for m, n, pairs in [(7, 1, T.PAIRS), (4, 3, T.PAIRS), (12, 5, T.PAIRS), (627, 364, T.PAIRS),
                    (465, 416, [('0', 'L'), ('0', '8'), ('L', 'U')]), (14885, 6076, [('0', 'L'), ('0', '8')])]:
    left, kill, bad, known = T.residual(m, n)
    for j, k in pairs:
        pts, pari_ = find(T.Pair(m, n, j, k))
        T.find_points = lambda Pr: ([], [-1, -1, -1, []])
        try: Pd = T.PairData(m, n, j, k)
        finally: T.find_points = real_find
        for P in pts: assert Pd.Pr.on(P)
        G = [P for P in pts if P[1] != 0]
        combos = [None]
        for g in G: combos = combos + [Pd.Pr.add(P, g) for P in combos]
        cls = set()
        for P in combos:
            if P is None or P[1] == 0: continue
            c = Pd.kum(P[0])
            cls.add(tuple(T.RI.sqfree_part(Fr(x).numerator * Fr(x).denominator) for x in c[:2]))
        res = {c: Pd.test(c)[0] for c in sorted(cls)}
        kres = {t: (Pd.test(Pd.cls_of_tuple(t))[0], Pd.test(Pd.cls_of_tuple(t, -1))[0]) for t in known}
        bad_ = [c for c, v in res.items() if v in ('killed', 'outside Sel2')] + [t for t, v in kres.items() if set(v) & {'killed', 'outside Sel2'}]
        allok &= not bad_
        say('(%d,%d) E_%s%s: dim Sel2 %d, torsion and psi(iota t0) images span %d; %d point classes: %s; known tuples: %s; %s'
            % (m, n, j, k, len(Pd.sel), Pd.ri, len(res), ' '.join('%s:%d' % kv for kv in sorted(collections.Counter(res.values()).items())),
               ' '.join('%s:%s/%s' % (t, *v) for t, v in sorted(kres.items())) or 'none', 'none killed' if not bad_ else 'KILLED: %s' % bad_))
say('t10\'s METHOD KILLS NO CLASS OF A RATIONAL POINT' if allok else 'A CLASS OF A RATIONAL POINT IS KILLED')
say.finish()
