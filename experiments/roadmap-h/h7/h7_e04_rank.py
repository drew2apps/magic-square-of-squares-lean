"""Beyond H7: an exact rank bound for the curves E04 of the census's Demjanenko-Manin tier (EXPERIMENT, written by an
AI model; not read by a person; nothing here changes the census).
Tier B of LEDGER.md counts 122 ratios closed by Demjanenko-Manin on curves of rank 1, 91 of them on
    E04: y^2 = x (x + 2 S m^2)(x + 2 S n^2),  S = m^2 + n^2   (G1_FULL_SYSTEM.md, Result 4; dm.py's model [0, 2S^2, 0, 4S^2m^2n^2, 0]),
and their upper rank bound 1 rests on PARI, Sage and Magma.  Here, for every curve E04 in scripts/g1_full/dm_evidence.txt:
  (1) the generator's abscissa of that line is checked to be the x of a point of E04 (x (x^2 + 2S^2 x + 4S^2m^2n^2) a
      non-zero square) and the point of infinite order (n G != O for n = 1 .. 12, Mazur, quoted), so rank >= 1;
  (2) E04 has full rational 2-torsion; the repository's complete 2-descent and Cassels-Tate pairing (rank_second.py's
      full_selmer2 and DFull.pairing, through h1/h1_rank0_exact.bound_full, unchanged) give rank <= dim Sel^2 - 2 - rank
      of the pairing; when that is not 1 (or the pairing's local search gives up), the same is tried on the curves with
      full 2-torsion reached from E04 by 2-isogenies (Velu, quoted; at most 6 curves), and the least bound is kept.
Nothing else of the closure (generator, heights, candidates) is touched: this replaces only the upper rank bound.
Output h7_e04_rank.txt.  usage: h7_e04_rank.py [--negative] [--write]   (h7_e04_rank_negative.txt)
Finds the repository through MSS_REPO, or else ../../..; H1's folder through H1_DIR, or else ../h1."""
import sys, os, io, contextlib
from fractions import Fraction as Fr
from math import isqrt
HERE = os.path.dirname(os.path.abspath(__file__))
WRITE = '--write' in sys.argv; NEG = '--negative' in sys.argv
REPO = os.environ.get('MSS_REPO', os.path.join(HERE, '..', '..', '..'))
H1 = os.environ.get('H1_DIR', os.path.join(HERE, '..', 'h1'))
sys.path.insert(0, H1)
_argv = sys.argv; sys.argv = [sys.argv[0]]
with contextlib.redirect_stdout(io.StringIO()):
    import h1_rank0_exact as R
sys.argv = _argv
EVID = os.path.join(REPO, 'scripts', 'g1_full', 'dm_evidence.txt')
OUT = []
def say(s=''): OUT.append(s)

def issq(q):
    q = Fr(q); return q >= 0 and isqrt(q.numerator) ** 2 == q.numerator and isqrt(q.denominator) ** 2 == q.denominator
def fsqrt(q):
    q = Fr(q); return Fr(isqrt(q.numerator), isqrt(q.denominator))

def evidence():
    lines = []
    for L in open(EVID):
        if L.startswith('#') or not L.strip(): continue
        w = L.split()
        if w[2] != 'E04': continue
        f = dict(x.split('=', 1) for x in w[3:] if '=' in x)
        lines.append((tuple(int(v) for v in w[0].split(',')), tuple(int(v) for v in w[1].split(',')), Fr(f['G'])))
    return lines

def point(a, b, x):
    v = x * (x * x + a * x + b)
    assert v != 0 and issq(v), 'the generator abscissa is not the x of a point of infinite order on E04'
    return (x, fsqrt(v))

def infinite_order(a, b, P):
    for k in range(1, 13):
        if R.mul(a, b, k, P) is None: return False
    return True

def least_bound(es, limit=6):
    todo = [tuple(sorted(es))]; seen = []; best = None; info = []
    while todo and len(seen) < limit:
        es = todo.pop(0); key = tuple(x - es[0] for x in es)
        if key in seen: continue
        seen.append(key)
        try: r, inf = R.bound_full(es)
        except AssertionError as ex:
            if 'no local point found' not in str(ex): raise
            info.append('pairing not computed'); r = None
        if r is not None:
            info.append('rank <= %d (%s)' % (r, inf))
            if best is None or r < best: best = r
            if best <= 1: break
        for e in es:
            nx = R.roots_of(*R.isog(*R.shift(es, e)))
            if nx: todo.append(tuple(sorted(nx)))
    return best, '%d curve(s): %s' % (len(seen), '; '.join(info))

def run_curve(m, n, gx, forge=None):
    S = m * m + n * n
    if forge == 'model': S += 1
    a, b = 2 * S * S, 4 * S * S * m * m * n * n
    es = (0, -2 * S * m * m, -2 * S * n * n)
    assert all(e * (e * e + a * e + b) == 0 for e in es), 'the stated roots are not the 2-torsion of the model'
    P = point(a, b, gx)
    if not infinite_order(a, b, P): raise AssertionError('the point has finite order')
    return least_bound(es)

def main():
    say('E04 of the census Demjanenko-Manin tier: exact upper rank bounds. Written by an AI model; not read by a person. Experiment, not counted.')
    ev = evidence(); runs = []
    for _, rp, gx in ev:
        if (rp, gx) not in runs: runs.append((rp, gx))
    rank1 = set()
    for (m, n), gx in runs:
        r, info = run_curve(m, n, gx)
        say('%d,%d E04: G on the curve, of infinite order; %s; so %s' % (m, n, info, 'rank = 1' if r == 1 else 'rank <= %s' % r))
        if r == 1: rank1.add((m, n))
    nr = sum(1 for _, rp, _ in ev if rp in rank1)
    say('curves: %d; rank exactly 1 (computed, exact): %d; ratios of the 91 whose E04 curve this covers: %d' % (len(runs), len(rank1), nr))
    say('not covered: %s' % ', '.join('%d,%d (curve %d,%d)' % (p + rp) for p, rp, _ in ev if rp not in rank1))

def negative():
    say('E04 exact rank bounds, negative and positive controls. Written by an AI model; not read by a person.')
    ev = evidence(); p, (m, n), gx = ev[0]
    for nm, fn in (('a forged model (S + 1) for %d,%d' % (m, n), lambda: run_curve(m, n, gx, forge='model')),
                   ('a forged generator x + 1 for %d,%d' % (m, n), lambda: run_curve(m, n, gx + 1)),
                   ('the 2-torsion point x = 0 as the generator', lambda: run_curve(m, n, Fr(0)))):
        try: fn(); say('%s: ACCEPTED' % nm); raise SystemExit(1)
        except AssertionError as ex: say('%s: refused (%s)' % (nm, str(ex)))
    # positive control: the full-2-torsion curve 2-isogenous to g4*K of (11165,7548), where PARI finds two independent points
    es = (-85759826449103309871318720, 0, 3204011570888134178800320)
    r, info = least_bound(es)
    say('positive control (rank >= 2 by PARI): %s; bound %s: %s' % (info, r, 'consistent' if r is not None and r >= 2 else 'NOT consistent'))
    assert r is not None and r >= 2

if __name__ == '__main__':
    negative() if NEG else main()
    s = '\n'.join(OUT) + '\n'; f = os.path.join(HERE, 'h7_e04_rank_negative.txt' if NEG else 'h7_e04_rank.txt')
    if WRITE: open(f, 'w').write(s)
    else:
        sys.stdout.write(s)
        if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS from ' + os.path.basename(f))
