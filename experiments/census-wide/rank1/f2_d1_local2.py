#!/usr/bin/env python3
"""F2: the E'-side class left open by rank_indep.rank_bound on D1 of (12007,1799) (the Demjanenko-Manin curve behind
(6903,5104)) has no Q_2-point on its first-descent quartic, so it is not in the Selmer group of the dual isogeny and the
exact bound drops from 4 to 3.  Uses only the repository's rank_indep.py at 24d6ef6 (unchanged).  Written by an AI model (Claude); not read by
a person.  Experiment, not counted.

    python3 f2_d1_local2.py              print the record (must equal f2_d1_local2.txt)
    python3 f2_d1_local2.py --write      also write f2_d1_local2.txt
    python3 f2_d1_local2.py --negative   controls that must come out the other way (--write: f2_d1_local2_negative.txt)

What it rests on: rank_indep.first_descent and rank_bound (the repository's exact 2-isogeny descent); rank_indep.solvable,
whose False is a certificate (every residue branch of the Hensel search dies before the depth); the record's arithmetic
(left classes on each side, bound = log2|left_E| + log2|left_E'| - 2)."""
import sys, os, io, contextlib
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get('MSS_REPO', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))   # experiments/census-wide/rank1 -> repository root
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, 'scripts', 'g1_full')); os.chdir(os.path.join(REPO, 'scripts', 'g1_full'))
import rank_indep as RI
import f2_curves
WRITE = '--write' in sys.argv; DEPTH = 40

def quartic(y, A, B):
    """the first-descent quartic of the class y for y^2 = x(x^2 + A x + B): y w^2 = y^2 + A y t^2 + B t^4, written as the
    form y (B t^4 + A y t^2 + y^2) that must be a square (low degree first, as rank_indep.solvable takes it)"""
    return [B * y, 0, A * y * y, 0, y ** 3]

def run():
    name, a, b, G = f2_curves.curves()[4]
    ap, bp = -2 * a, a * a - 4 * b
    with contextlib.redirect_stdout(io.StringIO()):
        bound, rec = RI.rank_bound(a, b, [G], RI.find_points([0, ap, 0, bp, 0]), say=lambda *x: None)
    out = ['F2, %s: the open class on the E\' side and the prime 2 (experiment; written by an AI model, not read by a person)' % name]
    out.append('rank_indep.rank_bound (exact): bound %d; E: Selmer %d, known %d, open cosets %d; E\': Selmer %d, known %d, open cosets %d' % (
        bound, len(rec['E']['selmer']), len(rec['E']['known']), len(rec['E']['open']),
        len(rec["E'"]['selmer']), len(rec["E'"]['known']), len(rec["E'"]['open'])))
    assert bound == 4 and len(rec["E'"]['open']) == 1
    coset = rec["E'"]['open'][0]
    dead = []
    for y in coset:
        s = RI.solvable([quartic(y, ap, bp)], 2, depth=DEPTH)
        out.append('  class %d: Q_2-point on its quartic: %s' % (y, {True: 'yes', False: 'no (certified)', None: 'undecided'}[s]))
        dead.append(s is False)
    assert all(dead), 'a class of the open coset has (or may have) a Q_2-point'
    leftp = set(rec["E'"]['known'])            # the other E' cosets were excluded by the second descent in rank_bound
    k = RI.pow2_floor(len(leftp)); kE = rec['E']['log2']
    out.append('so the whole coset lies outside Sel(phi^): E\' image <= 2^%d, E image <= 2^%d, rank <= %d + %d - 2 = %d' % (k, kE, kE, k, kE + k - 2))
    assert kE + k - 2 == 3
    out.append('result: exact rank bound 3 for %s (was 4); PARI ellrank gives [1,3]; rank 1 is still open' % name)
    return out, rec, ap, bp

def negative():
    out, rec, ap, bp = run()
    res = ['F2, D1 of (12007,1799): controls (each must come out the other way from the open class)']
    ok = True; found = 0
    for y in rec["E'"]['known']:
        if y == 1: continue
        s = RI.solvable([quartic(y, ap, bp)], 2, depth=12)
        good = s is not False; ok &= good; found += (s is True)
        res.append('  known class %d (image of a rational point, so a Q_2-point exists): solvable at depth 12 says %s -> %s' % (
            y, {True: 'yes', False: 'no', None: 'undecided'}[s], 'never "no", as it must' if good else 'WRONG: a false certificate'))
    ok &= found >= 1
    res.append('  (undecided is harmless: only "no" is used as a certificate; at least one known class is found soluble: %s)' % (found >= 1))
    y = rec["E'"]['open'][0][0]
    s = RI.solvable([quartic(y, ap, bp)], 2, depth=2)
    res.append('  open class %d at depth 2: %s -> %s' % (y, s, 'not decided at a shallow depth, as expected' if s is None else 'decided already'))
    res.append('all controls as they must be: %s' % ok)
    return res

if __name__ == '__main__':
    if '--negative' in sys.argv:
        res = negative(); print('\n'.join(res))
        if WRITE: open(os.path.join(HERE, 'f2_d1_local2_negative.txt'), 'w').write('\n'.join(res) + '\n')
    else:
        out = run()[0]; print('\n'.join(out))
        f = os.path.join(HERE, 'f2_d1_local2.txt')
        if WRITE: open(f, 'w').write('\n'.join(out) + '\n')
        elif os.path.exists(f): assert open(f).read() == '\n'.join(out) + '\n', 'differs from f2_d1_local2.txt'
