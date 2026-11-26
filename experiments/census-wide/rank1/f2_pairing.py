"""F2: the exact 2-isogeny descent and the Cassels-Tate pairing <y0, y1> on Sel^2 (the repository's rank_indep.py and
rank_second.py, unchanged, at commit 24d6ef6) on the five curves of f2_curves.py.

Written by an AI model (Claude); not read by a person.  Experiment, not counted.
    python3 f2_pairing.py [k]      curve k of f2_curves.CURVES only
"""
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
    try:
        res = R2.analyse(name, a, b, G)
        R2.say('  result: %s' % (res,))
    except AssertionError as ex:
        import traceback; R2.say("  stopped: " + traceback.format_exc().strip().splitlines()[-3].strip())
