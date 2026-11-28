import sys, time, json; sys.path.insert(0, '.')
import t11_lib as L
fn, lab, name, side, dd, pair = sys.argv[1:7]
C = L.curve(lab); Rd = L.RI.first_descent(C['a'], C['b'])
grp = Rd['Sp'] if side == "E'" else Rd['S']
basis = L.RF.rref_ints(sorted(grp, key=abs))
print('basis of this side S_1:', basis)
C3 = json.load(open(fn)); H = json.load(open('t11_hints.json'))
R = L.analyse3(name, side, C3, int(dd), hints=H)
out, S, cl = L.theta3(C, Rd, C3, R, pair, basis, side=side, say=lambda s: None)
print('Ob:', out)
