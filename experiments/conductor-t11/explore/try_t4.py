import sys, time, json, os; sys.path.insert(0, '.')
import t11_lib as L
fn, lab, name, side, other, dd, y, dS, dSp = sys.argv[1:10]
C = L.curve(lab)
C3 = json.load(open(fn))
G = L.theta4_gen_module(C, name, side, other, int(dd), int(y), C3, (int(dS), int(dSp)), os.path.abspath('t11_hints.json'))
if '--images' in sys.argv:
    G.choose_Xi_lin = L.choose_Xi_lin_images(G, L.RI.first_descent(C['a'], C['b']) if side == "E'" else L.RI.first_descent(C['ap'], C['bp']))
t0 = time.time()
if '--negative' in sys.argv: G.negative(name)
else: th = G.main(name)
print('time', time.time() - t0)
