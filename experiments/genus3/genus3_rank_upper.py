#!/usr/bin/env python3
"""Exact upper rank bounds for the 27 curves E1, E2, E3 of the nine open ratios (models of genus3_points.py, full
rational 2-torsion): the complete 2-descent and the Cassels-Tate pairing on Sel^2 of the repository's
scripts/g1_full/rank_second.py (full_selmer2 and DFull.pairing, imported unchanged), on the curve and on the curves
with full 2-torsion reached from it by rational 2-isogenies (Velu, quoted; at most 6 curves); the least bound is kept
(isogenous curves have the same rank).  rank <= dim Sel^2 - 2 - rank of the pairing.

Written by an AI model (build agent t12); not read by a person.  Experiment, tier D.  Nothing here changes a count.
What it rests on: rank_second.py as reviewed (reviews/review_t2.md; the Cassels-Tate pairing in its Brauer-Manin
form, Cassels 1998, Skorobogatov 6.2, quoted), its factorisations (PARI or sympy as finders, primes proved by
rank_indep.factor_int).  Lower bounds: genus3_points.py (exact).

    python3 -I genus3_rank_upper.py [--write]     genus3_rank_upper.txt   (MSS_REPO: the repository, default ../..)
"""
import sys, os, random, io, contextlib, time
from math import isqrt, gcd
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get('MSS_REPO', os.path.join(HERE, '..', '..'))
WRITE = '--write' in sys.argv
sys.path.insert(0, os.path.join(REPO, 'scripts', 'g1_full')); sys.path.insert(0, HERE)
_argv = sys.argv; sys.argv = [sys.argv[0]]
with contextlib.redirect_stdout(io.StringIO()):
    import rank_second as RS
    import genus3_points as G
sys.argv = _argv
OUT = []
def say(s=''):
    OUT.append(s); print(s, flush=True)

def bound_full(es):
    bad, loc, sel, kum, enc, places = RS.full_selmer2(es, random.Random(1))
    Ds = [RS.DFull(es, s_) for s_ in sel]
    Mx = [[Ds[i].pairing(b_)[0] for b_ in sel] for i in range(len(sel))]
    for i in range(len(sel)):
        assert Mx[i][i] == 1, 'the pairing is not alternating'
        for j in range(len(sel)): assert Mx[i][j] == Mx[j][i], 'the pairing is not symmetric'
    rho = len(RS.f2_span([sum((0 if v == 1 else 1) << j for j, v in enumerate(r)) for r in Mx]))
    return len(sel) - 2 - rho, len(sel), rho

def normal(es):
    e = sorted(es); d2, d3 = e[1] - e[0], e[2] - e[0]
    u = G.square_part(gcd(d2, d3))
    return (0, d2 // (u * u), d3 // (u * u))
def neighbours(es):
    """curves with full 2-torsion 2-isogenous to y^2 = (x - e1)(x - e2)(x - e3): for the kernel (e, 0), with u, v the
    other roots minus e, y^2 = x(x^2 + a x + b), a = -(u + v), b = u v, the image is y^2 = x(x^2 - 2a x + a^2 - 4b),
    whose quadratic has discriminant 16 u v: full 2-torsion iff u v is a square"""
    out = []
    for e in es:
        u, v = [x - e for x in es if x != e]
        if u * v > 0 and isqrt(u * v) ** 2 == u * v:
            a, b = -(u + v), u * v
            A2, B2 = -2 * a, a * a - 4 * b
            r = isqrt(A2 * A2 - 4 * B2); assert r * r == A2 * A2 - 4 * B2 and (-A2 + r) % 2 == 0
            out.append(normal((0, (-A2 + r) // 2, (-A2 - r) // 2)))
    return out

def bound_class(es0, limit=6):
    todo = [normal(es0)]; seen = []; best = None; info = []
    while todo and len(seen) < limit:
        es = todo.pop(0)
        if es in seen: continue
        seen.append(es)
        try:
            r, s, rho = bound_full(es)
            info.append('%d (Sel %d, pairing %d)' % (r, s, rho))
            best = r if best is None else min(best, r)
        except AssertionError as ex:
            if 'local' not in str(ex) and 'conic' not in str(ex): raise
            info.append('not computed (%s)' % str(ex)[:40])
        todo += [x for x in neighbours(es) if x not in seen]
    return best, len(seen), info

if __name__ == '__main__':
    say('genus3_rank_upper.py: exact upper rank bounds (complete 2-descent + Cassels-Tate pairing, rank_second.py)')
    only = [a for a in sys.argv[1:] if not a.startswith('--')]
    rows = []
    for (M, N) in G.OPEN:
        for i in (1, 2, 3):
            if only and '%d,%d,%d' % (M, N, i) not in only: continue
            t0 = time.time()
            best, k, info = bound_class(G.model(M, N, i))
            say('  (%d,%d) E%d: %d curve(s) with full 2-torsion in the 2-isogeny class; bounds %s => rank <= %s'
                % (M, N, i, k, ', '.join(info), best))
            rows.append(((M, N), i, best))
    say('Summary with the exact lower bound 1 of genus3_points.py:')
    for (MN, i, b) in rows:
        say('  %-11s E%d: %s' % ('(%d,%d)' % MN, i, 'rank = 1 (exact)' if b == 1 else '1 <= rank <= %s' % b))
    if WRITE: open(os.path.join(HERE, 'genus3_rank_upper.txt'), 'w').write('\n'.join(OUT) + '\n')
