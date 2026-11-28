import sys, time, json; sys.path.insert(0, '.')
import t11_lib as L
fn, name, side, dd, PAIR = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]
C3 = json.load(open(fn))
t0 = time.time()
R = L.analyse3(name, side, C3, dd)
import rank_fisher5 as F
print('dU core', F.core_small(R['dU']) if len(str(abs(R['dU']))) < 60 else 'big', 'dV core', F.core_small(R['dV']), 'Xi digits', len(str(R['Xi'])), time.time()-t0, flush=True)
key = '%s|%s' % (name, side) + ('|V' if PAIR == 'V' else '')
h = L.find_conic_norm(R, PAIR) if __import__("os").environ.get("NORM") else L.find_conic(R, PAIR, key)
print(h)
if h:
    H = json.load(open('t11_hints.json')) if __import__('os').path.exists('t11_hints.json') else {}
    H[key] = h; json.dump(H, open('t11_hints.json', 'w'), indent=1, sort_keys=True)
    R = L.analyse3(name, side, C3, dd, hints=H)
    print('rechecked:', R['pushout'] is not None, R['pushout3'] is not None, time.time()-t0)
