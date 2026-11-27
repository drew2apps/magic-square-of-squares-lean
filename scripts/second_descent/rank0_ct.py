#!/usr/bin/env python3
"""Rank 0 for the curve of (609,200) and (809,409) by an exact complete 2-descent with the Cassels-Tate pairing, on a
2-isogenous curve with full rational 2-torsion.

    python3 rank0_ct.py              prints the result and compares it with rank0_ct.txt
    python3 rank0_ct.py --negative   altered inputs that must be refused, and controls; compares with rank0_ct_negative.txt
    (add --write to rewrite the output instead of comparing)

Written by an AI model (Claude, Opus 5.5) on 27 November 2026 (UTC); not read by a person.  Experiment, not counted.
The method is the one of experiments/roadmap-h/h7 (a separate AI session): when a curve has one rational 2-torsion
point, the 2-isogenous curve may have full rational 2-torsion, and there the complete 2-descent and the Cassels-Tate
pairing of scripts/g1_full/rank_second.py (full_selmer2, DFull.pairing; reviewed in reviews/review_t2.md, which
re-implemented the pairing and found the same matrix) give rank <= dim Sel^2 - 2 - (rank of the pairing matrix).

WHY THIS MATTERS.  SECOND_DESCENT_CHAIN.md, step 7: for (609,200) and (809,409) (the same curve E, 'curve +14' of
indep_rank.txt) rank 0 is "not by any 2-descent": every 2-isogeny descent and the full 2-descents give 4, so if the
rank is 0 then Sha[2] has dimension 4; rank 0 rested on PARI's ellrank (its Cassels pairing, s = 4) and Magma's
TwoPowerIsogenyDescentRankBound, two programs that nobody here re-derived.  This script gets rank 0 from code in this
repository, in exact arithmetic: on E' (full 2-torsion) dim Sel^2 = 6 and the pairing has rank 4, so rank E'(Q) = 0,
and E and E' are isogenous, so rank E(Q) = 0.

WHAT IS CHECKED HERE (exact): E is read from indep_rank.txt (a-invariants and the 2-torsion abscissa x0); the shift
x = X + x0 gives y^2 = X(X^2 + A X + B) with A, B as recorded; B is a square, so E': Y^2 = X(X^2 - 2A X + A^2 - 4B)
(Velu's 2-isogeny with kernel (0,0)) has the three rational roots 0, A + 2 sqrt(B), A - 2 sqrt(B); #E(F_p) = #E'(F_p)
for every good prime p < 1000 (isogenous curves have the same point counts; a necessary check of the formula); the
complete 2-descent and the pairing on E' with two seeds for the local points; the pairing matrix (values +-1) symmetric, with +1 on
the diagonal.  QUOTED: Velu's formula, isogeny invariance of the rank, and the soundness of the complete 2-descent and
of the Cassels-Tate pairing as rank_second.py computes them (its header; Cassels 1998; Skorobogatov, Torsors and
rational points, 6.2).  PARI: only as a finder inside rank_second.py (qfsolve for the conic points, factor for
factorisations), each result checked; no PARI rank, Selmer or pairing output is used.  Reviewed by a separate AI
agent (z7, 27 November 2026: holds; its own complete 2-descent and pairing give 6 and 4 on all three full-2-torsion
curves of the class).  Tier D.
CONTROLS: the curve of (861,620) ('curve +23'), whose rank 0 is already exact (indep_rank.py), gives 0 here too; the
curve of (440,279), whose rank 0 rests on L(E,1) (lball_exact.py), gives only 2 here (the pairing is 0 on Sel^2: this
method cannot decide it, consistent with Sha containing (Z/4)^2); the congruent-number curves y^2 = x^3 - n^2 x for
n = 5, 6, 7 (rank 1, Tunnell / Heegner, quoted) must not get bound 0, and for n = 1, 2, 3 (rank 0) the bound is 0.
"""
import os, sys, json, random
from math import isqrt, gcd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'g1_full'))
import rank_second as RS

WRITE = '--write' in sys.argv
NEG = '--negative' in sys.argv
OUT = os.path.join(HERE, 'rank0_ct_negative.txt' if NEG else 'rank0_ct.txt')
lines = []
def out(s=''):
    print(s, flush=True); lines.append(s)

def record(pair, curve):
    for l in open(os.path.join(HERE, 'indep_rank.txt')):
        l = l.strip()
        if l.startswith('{'):
            d = json.loads(l)
            if d['pair'] == list(pair) and d['curve'] == curve: return d
    raise KeyError((pair, curve))

def count(f, p):
    """#E(F_p) for y^2 = f(x), f a monic cubic given by its coefficient function, p odd: p + 1 + sum of Legendre symbols"""
    s = 0
    for x in range(p):
        v = f(x) % p
        if v: s += 1 if pow(v, (p - 1) // 2, p) == 1 else -1
    return p + 1 + s

def odd_primes_below(n):
    return [p for p in range(3, n) if all(p % q for q in range(2, isqrt(p) + 1))]

def isogenous_full(d, iso):
    """E' with full rational 2-torsion for the 2-torsion point iso['x0'] of E, after the checks of the header"""
    a1, a2, a3, a4, a6 = d['ainvs']; assert a1 == 0 and a3 == 0
    x0, A, B = iso['x0'], iso['A'], iso['B']
    assert x0 ** 3 + a2 * x0 ** 2 + a4 * x0 + a6 == 0, 'x0 is not a root'
    assert A == 3 * x0 + a2 and B == 3 * x0 ** 2 + 2 * a2 * x0 + a4, 'the shift does not give (A, B)'
    s = isqrt(B) if B >= 0 else -1
    assert s >= 0 and s * s == B, 'B is not a square: no full 2-torsion on the isogenous curve'
    es = sorted([0, A + 2 * s, A - 2 * s])
    assert len(set(es)) == 3
    disc = (A * A - 4 * B) * B
    for p in odd_primes_below(1000):
        if disc % p == 0 or (A * A - 4 * B) % p == 0: continue
        nE = count(lambda x: x ** 3 + A * x * x + B * x, p)
        nEp = count(lambda x: (x - es[0]) * (x - es[1]) * (x - es[2]), p)
        assert nE == nEp, ('point counts differ: not isogenous', p, nE, nEp)
    return es

def bound(es, seeds=(1, 2), pairing=True):
    res = set()
    for seed in seeds:
        bad, loc, sel, kum, enc, places = RS.full_selmer2(es, random.Random(seed))
        if not pairing:
            res.add((len(sel), 0)); continue
        Ds = [RS.DFull(es, s_) for s_ in sel]
        Mx = [[Ds[i].pairing(b_)[0] for b_ in sel] for i in range(len(sel))]
        for i in range(len(sel)):
            assert Mx[i][i] == 1, 'the pairing is not alternating'
            for j in range(len(sel)): assert Mx[i][j] == Mx[j][i], 'the pairing is not symmetric'
        rho = len(RS.f2_span([sum((0 if v == 1 else 1) << j for j, v in enumerate(r)) for r in Mx]))
        res.add((len(sel), rho))
    assert len(res) == 1, ('the result depends on the seed', res)
    (dim, rho), = res
    return dim, rho, dim - 2 - rho

def show(label, es, **kw):
    dim, rho, b = bound(es, **kw)
    out('%s: roots %s; dim Sel^2 = %d, pairing of rank %d, rank <= %d' % (label, es, dim, rho, b))
    return b

if not NEG:
    out('# rank0_ct.py: rank 0 by the complete 2-descent and the Cassels-Tate pairing on a 2-isogenous curve with full 2-torsion')
    d = record((609, 200), '+14')
    out('curve +14 of (609,200) and (809,409): a-invariants %s, 2-torsion x0 = %s' % (d['ainvs'], d['rational_2_torsion_x']))
    es = isogenous_full(d, d['isogenies'][0])
    out("E': y^2 = x (x - e2)(x - e3) by Velu from (x0, 0); #E(F_p) = #E'(F_p) for every good odd p < 1000")
    b = show("E' of (609,200)", es)
    assert b == 0
    out('so rank E(Q) = 0 for the curve of (609,200) and (809,409) (isogeny invariance), from this repository\'s exact code')
    out('control, (861,620) curve +23 (rank 0 already exact by indep_rank.py):')
    d2 = record((861, 620), '+23')
    b2 = show("E' of (861,620)", isogenous_full(d2, d2['isogenies'][0]))
    assert b2 == 0
    out('control, (440,279) (rank 0 by L(E,1) in lball_exact.py; this method cannot decide it):')
    b3 = show('E of (440,279)', [-25921, 0, 516961])
    assert b3 == 2
else:
    out('# rank0_ct.py --negative: forged inputs refused; controls of known rank')
    d = record((609, 200), '+14'); iso = dict(d['isogenies'][0])
    for name, change in (('x0 + 1', ('x0', 1)), ('A + 2', ('A', 2))):
        f = dict(iso); f[change[0]] += change[1]
        try: isogenous_full(d, f); out('NOT REFUSED: %s' % name)
        except AssertionError as ex: out('refused, %s: %s' % (name, str(ex.args[0] if isinstance(ex.args[0], str) else ex.args[0][0])))
    d3 = record((223, 167), '1356'); i3 = [i for i in d3['isogenies'] if i['x0'] == -462297360][0]
    try: isogenous_full(d3, i3); out('NOT REFUSED: (223,167), x0 = -462297360')
    except AssertionError as ex: out('refused, (223,167) with x0 = -462297360, whose isogenous curve has no full 2-torsion: %s' % ex.args[0])
    es = isogenous_full(d, iso)
    forged = [es[0], es[1], es[2] + 2]
    A, B = iso['A'], iso['B']; s = isqrt(B)
    try:
        for p in odd_primes_below(200):
            if ((A * A - 4 * B) * B) % p == 0: continue
            assert count(lambda x: x ** 3 + A * x * x + B * x, p) == count(lambda x: (x - forged[0]) * (x - forged[1]) * (x - forged[2]), p), 'point counts differ'
        out('NOT REFUSED: a forged E\' (e3 + 2)')
    except AssertionError as ex: out('refused, a forged E\' (e3 + 2): %s' % ex.args[0])
    dim, rho, b = bound(es, pairing=False)
    out("without the pairing, E' of (609,200): dim Sel^2 = %d, rank <= %d (2-descent alone does not decide, as SECOND_DESCENT_CHAIN.md says)" % (dim, b))
    assert b == 4
    for n, r in ((5, 1), (6, 1), (7, 1), (1, 0), (2, 0), (3, 0)):
        b = show('congruent-number curve n = %d (rank %d)' % (n, r), [-n, 0, n])
        assert b >= r and (r > 0 or b == 0), 'a control of known rank is contradicted'
out('# end')

text = '\n'.join(lines) + '\n'
if WRITE:
    open(OUT, 'w').write(text); print('wrote', OUT)
else:
    old = open(OUT).read() if os.path.exists(OUT) else None
    assert old == text, 'output differs from %s' % os.path.basename(OUT)
    print('# matches', os.path.basename(OUT))
