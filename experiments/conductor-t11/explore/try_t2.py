import sys, time, json; sys.path.insert(0, '.')
import t11_lib as L
RI, R2 = L.RI, L.R2
lab = sys.argv[1]; side = sys.argv[2]; dds = [int(x) for x in sys.argv[3].split(',')]; etas = [int(x) for x in sys.argv[4].split(',')]
C = L.curve(lab); R = RI.first_descent(C['a'], C['b'])
AA, BB = (C['a'], C['b']) if side == 'E' else (C['ap'], C['bp'])
for dd in dds:
    P0 = RI.find_conic_point(dd, AA, BB // dd)
    surv, info = RI.second_descent(dd, AA, BB, P0, stop_at_first=False)
    R2.check_conic_forms(dd, AA, BB // dd, info['U'], info['V'], P0[2])
    n = 0
    for lam in sorted(surv, key=abs):
        try:
            D = R2.DLam(lam, info['U'], info['V'], bad=info['primes'])
            pv = {eta: D.pairing(eta)[0] for eta in etas}
        except AssertionError as ex:
            print('  lam', lam, 'skipped:', str(ex)[:60]); continue
        print(side, 'class', dd, 'lam', lam, 'pairings', pv, flush=True)
        n += 1
        if n >= 2: break
