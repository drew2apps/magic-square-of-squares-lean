import sys, time, json; sys.path.insert(0, '.')
import t11_lib as L
fn, lab, name, side, dd, pair = sys.argv[1:7]
ys = [int(x) for x in sys.argv[7].split(',')]
C = L.curve(lab); Rd = L.RI.first_descent(C['a'], C['b'])
C3 = json.load(open(fn)); H = json.load(open('t11_hints.json'))
t0 = time.time()
R = L.analyse3(name, side, C3, int(dd), hints=H)
out, S, cl = L.theta3(C, Rd, C3, R, pair, ys)
print(out, time.time() - t0)
