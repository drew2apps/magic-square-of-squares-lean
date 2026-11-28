import sys, time, json, os; sys.path.insert(0, '.')
import t11_lib as L
fn, name, side, dd, PAIR, big = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5], int(sys.argv[6])
C3 = json.load(open(fn))
R = L.analyse3(name, side, C3, dd)
key = '%s|%s' % (name, side) + ('|V' if PAIR == 'V' else '')
h = L.find_conic(R, PAIR, key, big=big)
print(h)
if h:
    H = json.load(open('t11_hints.json')); H[key] = h; json.dump(H, open('t11_hints.json', 'w'), indent=1, sort_keys=True)
    R = L.analyse3(name, side, C3, dd, hints=H)
    print('rechecked:', R['pushout'] is not None, R['pushout3'] is not None)
