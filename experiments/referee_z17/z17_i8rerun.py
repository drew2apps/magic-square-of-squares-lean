"""Referee z17 (written by an AI model; not read by a person): conductor-t10's own analysis, from scratch, on one ratio.
t10_pairct.analyse_ratio(m, n) with no tuples given, so the residual tuples and killing primes come from t10's own
residual() (not from I2's file), and the six pair curves are built by t10's PairData, unchanged.  Prints one JSON line:
  {"m", "n", "tuples": [{"tau", "kill", "pairs": {jk: label}, "killed_by"}], "curves": {jk: {"dimsel", "ri"}}}
dimsel = dim Sel^2 (t10's complete 2-descent), ri = rank of the span of the point images found (PARI points, checked).
usage: python3 -I z17_i8rerun.py m,n"""
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
os.environ['MSS_REPO'] = REPO
argv = sys.argv; sys.argv = [argv[0]]
sys.path.insert(0, os.path.join(REPO, 'experiments', 'conductor-t10'))
import t10_pairct as T
sys.argv = argv
T.PARI.default('nbthreads', 1)
T.say = lambda s='': None
PDS = []
_P = T.PairData
class _PD(_P):
    def __init__(s, *a, **k):
        _P.__init__(s, *a, **k); PDS.append(s)
T.PairData = _PD
m, n = map(int, argv[1].split(','))
left, kill, res, known = T.analyse_ratio(m, n)
tups = []
for t in left:
    if t in known: continue
    ks = [jk for jk, v in res[t].items() if v in ('killed', 'outside Sel2')]
    tups.append({'tau': list(t), 'kill': kill.get(t), 'pairs': res[t], 'killed_by': ks})
cur = {P.j + P.k: {'dimsel': len(P.sel), 'ri': P.ri, 'pari': list(P.pari[:2])} for P in PDS}
print(json.dumps({'m': m, 'n': n, 'known': len(known), 'tuples': tups, 'curves': cur}, sort_keys=True, separators=(', ', ': ')))
