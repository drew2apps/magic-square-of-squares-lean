"""If the named quotient had rank 0, would global.py's torsion check close the ratio?  usage: tors_check.py m,n,name ...  EXPERIMENT."""
import sys, json
R = '/tmp/claude-0/-home-user/9ec7a8d3-42dc-53c8-8518-1a663ac809f1/scratchpad/repo/scripts/g1_full'
import importlib.util
spec = importlib.util.spec_from_file_location('g', R + '/global.py'); g = importlib.util.module_from_spec(spec); spec.loader.exec_module(g)
for a in sys.argv[1:]:
    m, n, name = a.split(','); m, n = int(m), int(n)
    _, _, ws, T = g.torsion_w(g.curves(m, n)[name]); verdicts = set()
    for w in ws:
        if w is None: verdicts.add('degenerate'); continue
        if w < 0 or not g.issq(w): continue
        tt = g.sqrt_fr(w)
        for s in {tt, -tt}: verdicts.add(g.classify_t(m, n, s))
    print(json.dumps(dict(pair=[m, n], curve=name, torsion=T, w=[str(x) for x in ws], verdicts=sorted(verdicts), closes_if_rank0='SOLUTION' not in verdicts)))
