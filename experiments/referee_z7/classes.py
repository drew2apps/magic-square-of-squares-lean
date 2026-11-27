#!/usr/bin/env python3
"""Review z7 of roadmap H, H7 claim 4 and "What failed": isogeny classes, full 2-torsion, and the walk of h7_e04_rank.py.
Written by an AI model (review agent z7, 27 November 2026); not read by a person.  Experiment, not counted.

(1) For the E04 curves of (320,289), (1015,961), (1015,713) (dm.py's model): the isogeny class (PARI's ellisomat, a
    finder: the curves are not proved to be all of the class), which members have full rational 2-torsion (three
    rational roots, PARI's nfroots, checked by substitution), and my own Sel^2 and pairing (ct_mine.py) on each.
(2) The walk of h7_e04_rank.least_bound, replayed with H1's roots_of / isog / shift (imported, unchanged): how many
    distinct curves (by PARI's ellminimalmodel, i.e. up to Q-isomorphism) its "6 curve(s)" are.
(3) D1 (dm2.py's model) and sigma (sigma_exact.model) curves of tier B: no member of the isogeny class (ellisomat) has
    full rational 2-torsion, as REPORT_H7 says.
Rests on: as ct_mine.py; PARI's ellisomat for the list of the class (a finder).
PARI: ellisomat, nfroots, ellminimalmodel; and inside ct_mine.py.

    python3 classes.py            prints classes.txt   (about 1 minute)
"""
import sys, os, io, contextlib
from fractions import Fraction as Fr
from math import lcm
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
REPO = os.environ.get('MSS_REPO', os.path.join(HERE, '..', '..'))
H1 = os.environ.get('H1_DIR', os.path.join(HERE, '..', 'roadmap-h', 'h1'))
import ct_mine as C
pari = C.pari
sys.path.insert(0, H1); sys.path.insert(0, os.path.join(REPO, 'scripts', 'g1_full'))
_a = sys.argv; sys.argv = [sys.argv[0]]
with contextlib.redirect_stdout(io.StringIO()):
    import h1_rank0_exact as R
    import sigma_exact as SE
sys.argv = _a

def int_roots(a4, a6):
    r = [Fr(str(v)) for v in pari.nfroots(None, pari('x^3+(%s)*x+(%s)' % (a4, a6)))]
    if len(r) < 3: return None
    for v in r: assert v ** 3 + Fr(str(a4)) * v + Fr(str(a6)) == 0
    u = 1
    for v in r: u = lcm(u, v.denominator)
    return tuple(sorted(int(v * u * u) for v in r))

def canon(es):
    e1, e2, e3 = es
    E = pari.ellinit([0, -(e1 + e2 + e3), 0, e1 * e2 + e1 * e3 + e2 * e3, -e1 * e2 * e3])
    return str(pari.ellminimalmodel(E)[:5])

def e04(m, n):
    S = m * m + n * n; return [0, 2 * S * S, 0, 4 * S * S * m * m * n * n, 0], (0, -2 * S * m * m, -2 * S * n * n)

def main():
    print('# classes.py (review z7; written by an AI model)')
    print('## (1) isogeny classes of three E04 curves; my descent and pairing on each member with full 2-torsion')
    for m, n in ((320, 289), (1015, 961), (1015, 713)):
        ainv, _ = e04(m, n)
        L, Mx = pari.ellisomat(pari.ellinit(ainv))[:2]
        print('(%d,%d) E04: class of %d curves (ellisomat), degrees from E04 %s' % (m, n, len(L), [int(Mx[0][j]) for j in range(len(L))]))
        for j, c in enumerate(L):
            es = int_roots(c[0][0], c[0][1])
            if es is None: print('  curve %d: one rational 2-torsion point' % j); continue
            C.analyse(es, seed=9, label='  curve %d, full 2-torsion:' % j, show_matrix=False)
    print('## (2) the walk of h7_e04_rank.least_bound (6 visits), up to isomorphism')
    for m, n in ((1015, 713), (1015, 961), (320, 289), (15, 8)):
        _, es = e04(m, n)
        todo = [tuple(sorted(es))]; seen = []; order = []
        while todo and len(seen) < 6:
            es = todo.pop(0); key = tuple(x - es[0] for x in es)
            if key in seen: continue
            seen.append(key); order.append(canon(es))
            for e in es:
                nx = R.roots_of(*R.isog(*R.shift(es, e)))
                if nx: todo.append(tuple(sorted(nx)))
        idx = {}
        for c in order: idx.setdefault(c, len(idx))
        print('(%d,%d): the walk visits %d curves, %d up to isomorphism; sequence %s' % (m, n, len(order), len(idx), [idx[c] for c in order]))
    print('## (3) D1 and sigma curves of tier B: members of the class with full 2-torsion')
    def nfull(ainv):
        L = pari.ellisomat(pari.ellinit(ainv))[0]
        return len(L), sum(1 for c in L if int_roots(c[0][0], c[0][1]) is not None)
    seen = set(); d1 = []
    for Ln in open(os.path.join(REPO, 'scripts', 'g1_full', 'dm_evidence.txt')):
        if Ln.startswith('#') or not Ln.strip(): continue
        w = Ln.split()
        if w[2] != 'D1': continue
        m, n = map(int, w[1].split(','))
        if (m, n) in seen: continue
        seen.add((m, n))
        S = m * m + n * n; D2 = (m * m - n * n) ** 2
        c3 = 4 * S * S * m * m * n * n; c2 = -16 * S * m ** 3 * n ** 3; c1 = 4 * S * S * D2; c0 = -16 * S * m * n * D2
        d1.append(((m, n), nfull([0, c1, 0, c0 * c2, c0 * c0 * c3])))
    print('D1: %d curves; class sizes %s; curves with a full-2-torsion member: %s' % (len(d1), sorted(set(x[1][0] for x in d1)), [x[0] for x in d1 if x[1][1]]))
    sg = []
    for Ln in open(os.path.join(REPO, 'scripts', 'g1_full', 'sigma_evidence.txt')):
        if Ln.startswith('#') or not Ln.strip(): continue
        w = Ln.split(); m, n = map(int, w[0].split(',')); cv = w[1]
        M = SE.model(m, n, 2 if cv == 'E+' else -2)
        sg.append(((m, n), cv, nfull([0, M['a2'], 0, M['a4'], M['a6']])))
    print('sigma: %d curves; class sizes %s; curves with a full-2-torsion member: %s' % (len(sg), sorted(set(x[2][0] for x in sg)), [x[:2] for x in sg if x[2][1]]))
    print('# end')

if __name__ == '__main__':
    main()
