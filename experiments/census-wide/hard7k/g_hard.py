"""The repository's global.py rank-0 test, with PARI ellrank at a higher effort over the 2-isogeny class.
Curves, torsion and t-classification are global.py's own functions.  usage: g_hard.py effort m,n ...  EXPERIMENT."""
import sys, json, time
R = '/tmp/claude-0/-home-user/9ec7a8d3-42dc-53c8-8518-1a663ac809f1/scratchpad/repo/scripts/g1_full'
import importlib.util
spec = importlib.util.spec_from_file_location('g', R + '/global.py'); g = importlib.util.module_from_spec(spec); spec.loader.exec_module(g)
pari = g.pari; pari.allocatemem(4*10**9); eff = int(sys.argv[1])
for a in sys.argv[2:]:
    m, n = map(int, a.split(','))
    for name, poly in g.curves(m, n).items():
        c = [poly[4], poly[3], poly[2], poly[1]] if len(poly) == 5 else poly
        c0, c1, c2, c3 = c
        E = pari.ellinit([0, c2, 0, c1*c3, c0*c3*c3]); t = time.time()
        lo, hi = 0, 99
        for item in pari.ellisomat(pari.ellminimalmodel(E), 0, 1)[0]:
            rk = pari.ellrank(pari.ellinit(item), eff); lo = max(lo, int(rk[0])); hi = min(hi, int(rk[1]))
        res = dict(pair=[m, n], curve=name, rank=[lo, hi], effort=eff, sec=round(time.time()-t, 1))
        if hi == 0:
            # reuse global.py's torsion enumeration and t-classification (its own rank call is ignored)
            _, _, ws, T = g.torsion_w(poly); verdicts = set()
            for w in ws:
                if w is None: verdicts.add('degenerate'); continue
                if w < 0 or not g.issq(w): continue
                tt = g.sqrt_fr(w)
                for s in {tt, -tt}: verdicts.add(g.classify_t(m, n, s))
            res.update(torsion=T, verdicts=sorted(verdicts), closes='SOLUTION' not in verdicts)
        print(json.dumps(res), flush=True)
