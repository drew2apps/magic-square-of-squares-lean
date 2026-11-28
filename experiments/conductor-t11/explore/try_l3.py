import sys, time, json; sys.path.insert(0, '.')
import t11_lib as L
lab, side, dd = sys.argv[1], sys.argv[2], int(sys.argv[3])
t0 = time.time()
C = L.curve(lab); R = L.RI.first_descent(C['a'], C['b'])
print('undecided local images:', [p for p in R['bad'] if R['locS'][p][0] != R['locS'][p][1] or R['locSp'][p][0] != R['locSp'][p][1]])
import os
sk = int(os.environ.get('SKIP', '0'))
C3, pv = L.level3_lift(C, R, side, dd, skip=sk)
print(time.time() - t0)
json.dump(C3, open('l3_%s_%s_%d%s.json' % (lab.replace('(','_').replace(')','').replace(',','_').replace('+','p').replace('-','m'), side.replace("'", 'p'), dd, ('_s%d' % sk) if sk else ''), 'w'), indent=1)
