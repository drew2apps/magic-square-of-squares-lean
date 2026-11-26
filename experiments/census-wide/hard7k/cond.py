import sys
R = '/tmp/claude-0/-home-user/9ec7a8d3-42dc-53c8-8518-1a663ac809f1/scratchpad/repo/scripts/g1_full'
import importlib.util
spec = importlib.util.spec_from_file_location('g', R + '/global.py'); g = importlib.util.module_from_spec(spec); spec.loader.exec_module(g)
pari = g.pari
for a in sys.argv[1:]:
    m, n, name = a.split(','); m, n = int(m), int(n)
    poly = g.curves(m, n)[name]
    c = [poly[4], poly[3], poly[2], poly[1]] if len(poly) == 5 else poly
    c0, c1, c2, c3 = c
    E = pari.ellminimalmodel(pari.ellinit([0, c2, 0, c1*c3, c0*c3*c3]))
    N = pari.ellglobalred(E)[0]
    print(m, n, name, 'conductor', N, 'digits', len(str(N)), 'factor', pari.factor(N), 'ainv', E[:5])
