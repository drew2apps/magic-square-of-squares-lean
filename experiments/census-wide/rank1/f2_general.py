"""F2: rank_second.analyse_general (the matrix of Theta_2 = <y, lift(x)> between the open cosets of the two sides of the
2-isogeny descent; its F_2-rank rho lowers the bound by 2 rho), unchanged from commit 24d6ef6, on the curves of
f2_curves.py with more than one open coset.  Written by an AI model (Claude); not read by a person.  Experiment, not counted.
    python3 f2_general.py k     curve k of f2_curves.CURVES"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get('MSS_REPO', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))   # experiments/census-wide/rank1 -> repository root
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, 'scripts', 'g1_full'))
os.chdir(os.path.join(REPO, 'scripts', 'g1_full'))
import rank_indep as RI, rank_second as R2
import f2_curves
which = [int(x) for x in sys.argv[1:] if x.isdigit()]
for k, (name, a, b, G) in enumerate(f2_curves.curves()):
    if which and k not in which: continue
    R2.say('== ' + name)
    res = R2.analyse_general(name, a, b, [G], RI.find_points([0, -2 * a, 0, a * a - 4 * b, 0]))
    R2.say('  result: %s' % (res,))
