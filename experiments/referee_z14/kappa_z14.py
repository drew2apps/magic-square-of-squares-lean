"""z14: EXPLORATION_4 §5's class constants on the mixed curves, on ratios outside the pilot (ratios.txt): the tuple class
of t10_pairct is (x - e_r : kappa_r, x - e* : +-kappa* D), and kappa.py's fit says kappa_r = 2 and kappa* = sqf(mn(m^2-n^2))
times 1 (0L, 8U) or 2 (0U, 8L), e_r = e1, e* = e2.  t10_pairct.PairData reads the constants off the rational points it
knows; here it is given no PARI ellrank points (an empty hint), so it uses psi(iota t0), the 2-torsion and their sums
(all constants must agree: asserted by PairData).
Written by an AI model (reviewer z14); not read by a person.   usage: python3 -I kappa_z14.py | grep -v '^PARI' | diff - kappa_z14.txt"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.argv = sys.argv[:1]
sys.path.insert(0, os.path.join(HERE, '..', 'conductor-t10'))
import t10_pairct as T
T.PARI.allocatemem(2 * 10 ** 8, 8 * 10 ** 8)
RI = T.RI
from fractions import Fraction as Fr
def sqf(x): x = Fr(x); return RI.sqfree_part(x.numerator * x.denominator)
lines = [l.split() for l in open(os.path.join(HERE, 'ratios.txt')) if l.strip() and l[0] != '#']
pick = lines[0:12] + lines[120:128]
print('# z14 kappa check. Written by an AI model (reviewer z14); not read by a person.')
tot = {}
for l in pick:
    m, n = map(int, l[0].split(','))
    out = []
    for j, k in (('0', 'L'), ('0', 'U'), ('8', 'L'), ('8', 'U')):
        Pr = T.Pair(m, n, j, k)
        key = '%d,%d %s%s %s' % (m, n, j, k, Pr.t0)
        T.HINTS.setdefault('points', {})[key] = {'pts': [], 'pari': [-1, -1, -1, []]}
        Pd = T.PairData(m, n, j, k)
        want = sqf(m * n * (m * m - n * n)) * (1 if j + k in ('0L', '8U') else 2)
        ok = (Pd.role['r'], Pd.role['*']) == (0, 1) and Pd.kr == 2 and sqf(Pd.ks * want) == 1
        okpm = (Pd.role['r'], Pd.role['*']) == (0, 1) and Pd.kr == 2 and abs(sqf(Pd.ks * want)) == 1
        out.append('%s%s:%s' % (j, k, 'ok' if ok else ('ok up to sign' if okpm else 'DIFFERS kr=%d ks=%d want=%d roles=%s' % (Pd.kr, Pd.ks, want, Pd.role))))
        tot[out[-1].split(':')[1][:7]] = tot.get(out[-1].split(':')[1][:7], 0) + 1
    print('%d,%d %s %s' % (m, n, ' '.join(l[1:3]), ' '.join(out)), flush=True)
print('total', dict(sorted(tot.items())))
