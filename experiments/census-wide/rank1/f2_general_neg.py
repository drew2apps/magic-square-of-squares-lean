"""F2 controls for f2_general.py (Theta_2 matrix of rank_second.analyse_general on E1 of (8249,3720)).
Written by an AI model (Claude); not read by a person.  Experiment, not counted.
    python3 f2_general_neg.py [--write]     (--write: f2_general_negative.txt)
1. The repository's own rank_second.negative() (altered inputs that must be refused), run with sympy.
2. Positive control on E1 of (8249,3720): a class of the E' side that is the image of a known rational point must pair
   to +1 with the lift of every open coset of the E side (a point pairs trivially), so the matrix rows that lower the
   bound cannot come from point classes."""
import sys, os, io, contextlib
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get('MSS_REPO', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))   # experiments/census-wide/rank1 -> repository root
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, 'scripts', 'g1_full')); os.chdir(os.path.join(REPO, 'scripts', 'g1_full'))
import rank_indep as RI, rank_second as R2
import f2_curves
OUT = []
def say(s): OUT.append(s); print(s, flush=True)
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    ok, res = R2.negative()
lines = [l for l in buf.getvalue().splitlines() if l.startswith('  refused') or l.startswith('  NOT REFUSED')]
say('1. repository rank_second.negative(): %d of %d refused' % (sum(res), len(res)))
for l in lines: say('   ' + l.split(' (')[0].strip())
name, a, b, G = f2_curves.curves()[3]
with contextlib.redirect_stdout(io.StringIO()):
    r, rec = RI.rank_bound(a, b, [G], RI.find_points([0, -2 * a, 0, a * a - 4 * b, 0]), say=lambda s: None)
X = [c[0] for c in rec['E']['open']]
known = [k for k in rec["E'"]['known'] if k != 1]
say('2. %s: open E-side cosets %s; known E\'-side classes %s' % (name, X, known))
good = True
for y in known[:2]:
    for x in X:
        with contextlib.redirect_stdout(io.StringIO()):
            v, lams = R2.lift_pairing(x, a, b, y)
        say('   <%d (a point class), lift of C_%d> = %+d (must be +1)' % (y, x, v)); good &= (v == 1)
say('negative checks: %s' % ('all as required' if all(res) and good else 'FAILED'))
if '--write' in sys.argv: open(os.path.join(HERE, 'f2_general_negative.txt'), 'w').write('\n'.join(OUT) + '\n')
