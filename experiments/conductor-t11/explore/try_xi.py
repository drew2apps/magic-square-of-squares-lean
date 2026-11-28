import sys, time, json, os; sys.path.insert(0, '.')
import t11_lib as L
fn, lab, name, side, other, dd, y, dS, dSp = sys.argv[1:10]
C = L.curve(lab); C3 = json.load(open(fn))
G = L.theta4_gen_module(C, name, side, other, int(dd), int(y), C3, (int(dS), int(dSp)), os.path.abspath('t11_hints.json'))
P = G.prepare(name)
D, Cc, RD = P['D'], P['C'], P['RD']
nd, ps = G.pencil_primes(D, known=sorted(set([2] + RD['bad'])))
B0 = sorted(set([2] + RD['bad'] + ps))
M = G.T.Model(D['QA'], D['QB'])
places = sorted(set([2, 3, 5, 7, 11, 13] + B0))
sols, K, Rs = G.choose_Xi_lin(M, D, Cc, places)
print('sols', sols[:5]); print('real signs', sorted(Rs))
for p in places: print(p, sorted(K[p]))
R = L.RI.first_descent(C['a'], C['b'])
for p in places:
    for nm, loc in (('E', R['locS']), ("E'", R['locSp'])):
        if p in loc: print(p, nm, 'poss', len(loc[p][0]), 'cert', len(loc[p][1]), end='; ')
    print('sampled', len(K[p]))
