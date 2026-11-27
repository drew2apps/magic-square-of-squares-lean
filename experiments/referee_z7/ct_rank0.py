#!/usr/bin/env python3
"""Review z7 of scripts/second_descent/rank0_ct.py: rank 0 for curve +14 of (609,200) and (809,409), independently.
Written by an AI model (review agent z7, 27 November 2026); not read by a person.  Experiment, not counted.

(1) Same curve: magma_second_ranks.json lists the same a-invariants for (609,200) and (809,409); the rank-0 curves of
    second_123.jsonl for both pairs (twist +1 and -1 of their E) are isomorphic over Q to it (exact: c4' = u^4 c4,
    c6' = u^6 c6 with u^2 a rational square).
(2) The 2-torsion root x0 (sympy's integer roots), the shift x = X + x0 (A, B equal to indep_rank.txt's record), B a
    square, E' = Velu with roots 0, A +- 2 sqrt(B) (velu_check.py), equal to the first full-2-torsion member of the
    record's class_full_2descent.
(3) My own Sel^2 and Cassels-Tate pairing (my2desc.py, ct_mine.py), two seeds, on E' and on the other two
    full-2-torsion members of the class, with the kernel compared with the images of the 2-torsion; the controls of
    rank0_ct.py: (861,620) curve +23, and the L-value curve of (440,279) (lball_exact.CURVE).  Compared with
    rank0_ct.txt.
Rests on: as ct_mine.py.  PARI: inside ct_mine.py only.

    python3 ct_rank0.py            prints ct_rank0.txt   (about 30 s)
"""
import sys, os, re, json
from fractions import Fraction as Fr
from math import isqrt
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
REPO = os.environ.get('MSS_REPO', os.path.join(HERE, '..', '..'))
SD = os.path.join(REPO, 'scripts', 'second_descent')
import ct_mine as C

def c46(a):
    a1, a2, a3, a4, a6 = [Fr(x) for x in a]
    b2 = a1 * a1 + 4 * a2; b4 = 2 * a4 + a1 * a3; b6 = a3 * a3 + 4 * a6
    return b2 * b2 - 24 * b4, -b2 ** 3 + 36 * b2 * b4 - 216 * b6

def isomorphic(a, b):
    (c4, c6), (d4, d6) = c46(a), c46(b)
    u2 = (d6 / c6) / (d4 / c4)
    return u2 > 0 and isqrt(u2.numerator) ** 2 == u2.numerator and isqrt(u2.denominator) ** 2 == u2.denominator \
        and d4 == u2 * u2 * c4 and d6 == u2 ** 3 * c6, u2

def record(pair, curve):
    for l in open(os.path.join(SD, 'indep_rank.txt')):
        if l.startswith('{'):
            d = json.loads(l)
            if d['pair'] == list(pair) and d['curve'] == curve: return d

def eprime(ainvs):
    a1, a2, a3, a4, a6 = ainvs; assert a1 == 0 and a3 == 0
    x = sp.symbols('x')
    rts = [int(r) for r in sp.Poly(x ** 3 + a2 * x ** 2 + a4 * x + a6, x).ground_roots() if r.is_Integer]
    assert len(rts) == 1
    x0 = rts[0]; A = 3 * x0 + a2; B = 3 * x0 * x0 + 2 * a2 * x0 + a4
    s = isqrt(B); assert s * s == B
    return x0, A, B, tuple(sorted((0, A + 2 * s, A - 2 * s)))

def main():
    print('# ct_rank0.py (review z7; written by an AI model): rank0_ct.py\'s curves, independently')
    mj = {tuple(c['pair']): c for c in json.load(open(os.path.join(SD, 'magma_second_ranks.json')))}
    a609, a809 = mj[(609, 200)]['ainvs'], mj[(809, 409)]['ainvs']
    print('(1) magma_second_ranks.json: curve %s of (609,200), curve %s of (809,409); same a-invariants: %s'
          % (mj[(609, 200)]['curve'], mj[(809, 409)]['curve'], a609 == a809))
    for l in open(os.path.join(SD, 'second_123.jsonl')):
        d = json.loads(l)
        if tuple(d.get('pair', [])) in ((609, 200), (809, 409)):
            cs = d['curves']
            if isinstance(cs, str):
                import ast; cs = ast.literal_eval(cs)
            c = [c for c in cs if c.get('curve') == '+14'][0]
            ok, u2 = isomorphic(a609, c['rank0_curve'])
            print('    second_123.jsonl %s: rank-0 curve (twist %s) isomorphic over Q to it: %s (u^2 = %s)' % (d['pair'], c['rank0_twist'], ok, u2))
            assert ok
    d = record((609, 200), '+14'); assert d['ainvs'] == a609
    x0, A, B, es = eprime(a609); iso = d['isogenies'][0]
    print('(2) x0 = %d (record %d), A = %d, B = %d (record: %s); E\' roots %s; record class_full_2descent[0] roots %s'
          % (x0, iso['x0'], A, B, 'same' if (A, B) == (iso['A'], iso['B']) else 'DIFFERENT', es, d['class_full_2descent'][0]['roots']))
    assert (x0, A, B) == (iso['x0'], iso['A'], iso['B']) and tuple(sorted(d['class_full_2descent'][0]['roots'])) == es
    theirs = {}
    for l in open(os.path.join(SD, 'rank0_ct.txt')):
        mm = re.match(r"(.*?): roots \[([-\d, ]+)\]; dim Sel\^2 = (\d+), pairing of rank (\d+)", l)
        if mm: theirs[tuple(int(v) for v in mm.group(2).split(','))] = (int(mm.group(3)), int(mm.group(4)))
    print('(3) my descent and pairing')
    for seed in (3, 4):
        r = C.analyse(es, seed=seed, xs=[], label='  E\' of +14, seed %d:' % seed)
        assert theirs[es] == (r['dim'], r['rho'])
    print('    rank0_ct.txt: dim Sel^2 = %d, pairing rank %d: agrees' % theirs[es])
    for c in d['class_full_2descent'][1:]:
        if 'roots' in c: C.analyse(tuple(sorted(c['roots'])), seed=5, xs=[], label='  class member %s:' % (tuple(c['roots']),))
    d2 = record((861, 620), '+23'); _, _, _, es2 = eprime(d2['ainvs'])
    r = C.analyse(es2, seed=3, xs=[], label='  control (861,620) +23, E\' %s:' % (es2,))
    assert theirs[es2] == (r['dim'], r['rho']); print('    rank0_ct.txt: agrees')
    es3 = (-25921, 0, 516961)       # lball_exact.CURVE = (0, -491040, 0, -13400146081, 0) = x (x - 719^2)(x + 161^2)
    assert -(es3[0] + es3[1] + es3[2]) == -491040 and es3[0] * es3[2] == -13400146081
    r = C.analyse(es3, seed=3, xs=[], label='  control (440,279), lball_exact\'s curve:')
    assert theirs[es3] == (r['dim'], r['rho']); print('    rank0_ct.txt: agrees')
    print('# end')

if __name__ == '__main__':
    main()
