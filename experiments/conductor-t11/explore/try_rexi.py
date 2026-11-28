import sys, json; sys.path.insert(0, '.')
import t11_lib as L
fn, lab, side, s, out = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]
C = L.curve(lab); R = L.RI.first_descent(C['a'], C['b'])
C3 = json.load(open(fn))
C3b = L.recertify_xi(C3, R, side, s)
json.dump(C3b, open(out, 'w'), indent=1)
