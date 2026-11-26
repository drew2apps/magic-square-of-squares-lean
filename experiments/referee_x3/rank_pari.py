"""x3: PARI 2.17.4 ranks (hints only) for g4*K of the closures and of the conditional ratios"""
import sys; import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dm_indep as D
pari = D.pari
for a in sys.argv[1:]:
    m, n = map(int, a.split(','))
    P, C = D.model(m, n)
    E = pari.ellinit([0, C.A2, 0, C.A4, C.A6])
    out = []
    for Fc in pari.ellisomat(pari.ellminimalmodel(E), 0, 1)[0]:
        F = pari.ellinit(Fc); rk = pari.ellrank(F)
        out.append([int(rk[0]), int(rk[1]), int(rk[2])])
    print(a, 'isogeny class ellrank [r,R,s]:', out, 'rootno', int(pari.ellrootno(E)), flush=True)
