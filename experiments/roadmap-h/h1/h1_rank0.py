"""H1, part 6: rank-0 quotients of the four-form curve C at the phantom ratios (EXPERIMENT, written by an AI model; not
read by a person; nothing here changes the census).
C is the curve of h1: F_0(t), F_8(t), F_L(t), F_U(t) all squares (the four forms of tuple_descent's magic system, in
the t of h1_pairs.forms).  With w = t^2, f0 = m^2 + n^2 w, f8 = n^2 + m^2 w and K = (m^2 + n^2)^2 (w + 1)^2 - 16 m^2 n^2 w
(= F_L F_U), the repository's scripts/g1_full/global.py lists six genus-1 quotients of C that do not use the fifth form
g4: y^2 = w f0 f8, f0 K, f8 K, w K, w f0 K, w f8 K (a point t of C gives a point with first coordinate w = t^2).  If
one of them has rank 0, every point of C has w among the first coordinates of its torsion points; then every tuple of a
point of C is among the tuples of those t = +-sqrt(w) that lie on C, and an extra tuple that is not one of them is the
class of no rational point.  global.torsion_w (PARI's ellrank, with the isogenous curves for the upper bound, and
elltors) is used unchanged: the rank bound is PARI's (a 2-descent, not re-checked here), the torsion points and the
tuples are exact.
Output h1_rank0.txt: per phantom, the six curves' rank bounds, and for a rank-0 curve the tuples of C over its torsion
points and the verdict on each extra tuple.
usage: h1_rank0.py [--write]      The repository is found through MSS_REPO, or else ../../.."""
import sys, os, re, ast, glob
from fractions import Fraction as Fr
from math import isqrt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
WRITE = '--write' in sys.argv
_argv = sys.argv; sys.argv = [sys.argv[0]]
import h1_pairs as HP
sys.argv = _argv
REPO = os.environ.get('MSS_REPO', os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(REPO, 'scripts', 'g1_full'))
import importlib
GLB = importlib.import_module('global')
NAMES = ['w*f0*f8', 'f0*K', 'f8*K', 'w*K', 'w*f0*K', 'w*f8*K']

def phantoms():
    out = [(1109503, 597671, [(1109503, 597671, 2, 1), (1195342, 2219006, 2, 1)],
            [(1109503, 597671, 2, 959), (1109503, 597671, 1918, 1), (1195342, 2219006, 2, 959), (1195342, 2219006, 1918, 1)])]
    for f in sorted(glob.glob(os.path.join(HERE, 'h1_scan_k*.txt'))):
        for l in open(f):
            mm = re.match(r'\((\d+),(\d+)\) PHANTOM: .*known (\[.*\]); extra (\[.*\])', l.strip())
            if mm: out.append((int(mm.group(1)), int(mm.group(2)), ast.literal_eval(mm.group(3)), ast.literal_eval(mm.group(4))))
    return sorted(out)

def tuple_at(m, n, t):
    """the tuple of t if t (a rational number, or None for infinity) is on C, else None"""
    if t is None: f = {'0': (Fr(0), Fr(n)), '8': (Fr(0), Fr(m)), 'L': (Fr(m + n), Fr(m - n)), 'U': (Fr(m - n), Fr(m + n))}
    else: f = HP.forms(m, n, Fr(t))
    tp = []
    for j in '08LU':
        z = f[j]; Nz = z[0] ** 2 + z[1] ** 2
        if Nz == 0 or not HP.isq(Nz): return None
        tp.append(HP.gauss_class(z))
    return tuple(tp)

def main():
    say = HP.say
    say('H1: rank-0 quotients of C at the phantom ratios. Written by an AI model; not read by a person. Experiment, not counted.')
    say('rank bounds: PARI ellrank (and isogenous curves), as in global.py; torsion and tuples exact')
    tally = []
    for m, n, known, extra in phantoms():
        HP.set_base(set(q for x in (2, m, n, m - n, m + n, m * m + n * n) for q in HP.RI.factor_int(x)) | {2})
        P = GLB.curves(m, n); say(''); say('(%d,%d): known %s; extra %s' % (m, n, known, extra))
        dead = None
        for nm in NAMES:
            lo, hi, ws, T = GLB.torsion_w(P[nm])
            line = '  y^2 = %s: rank %d..%d, torsion %d' % (nm, lo, hi, T)
            if hi == 0:
                ts = set()
                for w in ws:
                    if w is None: ts.add(None); continue
                    if w >= 0 and HP.isq(w): r = HP.rsqrt(w); ts |= {r, -r}
                tps = sorted(set(x for x in (tuple_at(m, n, t) for t in ts) if x is not None))
                line += '; t over its torsion points: %s; their tuples on C: %s' % (sorted(ts, key=lambda x: (x is None, x if x is not None else 0)), tps)
                if dead is None: dead = (nm, tps)
            say(line)
        if dead:
            nm, tps = dead
            for k in known: assert k in tps, 'a known tuple is not the tuple of a torsion point: the method is wrong'
            v = [tp for tp in extra if tp not in tps]
            say('  => y^2 = %s has rank 0 (PARI): every point of C is over its torsion points; extra tuples that are no tuple of those: %d of %d' % (nm, len(v), len(extra)))
            tally.append(((m, n), len(v) == len(extra), nm))
        else:
            say('  => no rank-0 curve among the six')
            tally.append(((m, n), False, None))
    say(''); say('SUMMARY:')
    for mn, ok, nm in tally: say('  %s: %s' % (mn, ('all extra tuples excluded by y^2 = %s (rank 0 by PARI)' % nm) if ok else 'not excluded here'))
    say('phantom ratios: %d; all extras excluded by a rank-0 quotient: %d' % (len(tally), sum(1 for _, ok, _ in tally if ok)))

if __name__ == '__main__':
    main()
    s = '\n'.join(HP.OUT) + '\n'; f = os.path.join(HERE, 'h1_rank0.txt')
    if WRITE: open(f, 'w').write(s)
    else:
        sys.stdout.write(s)
        if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS from ' + os.path.basename(f))
