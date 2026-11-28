#!/usr/bin/env python3
"""reach.py (conductor-t13): what the searched regions mean on E (a guide, not a bound).  EXPERIMENT; written by an
AI model (build agent t13), not read by a person.

    python3 reach.py --H=50000 [--write]   -> reach_H50000.txt

For each coset covering of covers.json, Gaussian integer pairs t near the boundary h_p(t) = H of the searched
ellipsoid (32 directions, fixed seed) are mapped by t -> x = e1 + d1 (u/z)^2, (u : v : z) = X(M t), and the absolute
logarithmic height h(x) (Sage global_height) is printed (min / median / max).  A point of E(K) found on that boundary
would have an x-coordinate of about that height, and canonical height (PARI/Sage normalisation, h(x) + O(1)) of
about that size, up to the height difference of E, which is not bounded here.  For E3 the point would be mapped to
E by the dual isogeny, which doubles the canonical height.  Not a certified statement about E(K).
"""
from sage.all__sagemath_schemes import *
import sys, os, json, random, math
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'scripts', 'second_descent'))
import kctp_recheck as R

K = R.K; I = R.I
ARGS = dict(a.lstrip('-').split('=', 1) if '=' in a else (a.lstrip('-'), '1') for a in sys.argv[1:])
WRITE = 'write' in ARGS
OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(str(s))
def kq(p): return K(QQ(p[0])) + K(QQ(p[1])) * I

def main():
    H = float(ARGS.get('H', '50000'))
    rng = random.Random(20261127)
    cov = json.load(open(os.path.join(HERE, 'covers.json')))['curves']
    say('reach.py (conductor-t13), H = %g.  Written by an AI model (build agent t13); not read by a person.  A guide.' % H)
    for cv in ('E', 'E3'):
        es = [kq(x) for x in cov[cv]['e']]
        for cl in cov[cv]['classes']:
            if not cl['label'].startswith('coset'): continue
            zc = complex(cl['covariant'][0], cl['covariant'][1]); T = cl['covariant'][2]
            (a, b), (c, d) = [[kq(x) for x in row] for row in cl['M']]
            d1 = kq(cl['d'][0])
            hs = []
            for k in range(32):
                # a direction: t2 = r2 e^{i phi2}, t1 = z t2 + r1 e^{i phi1}, with r1^2/T + T r2^2 = H
                f = rng.random(); r2 = math.sqrt(f * H / T); r1 = math.sqrt((1 - f) * H * T)
                p1, p2 = rng.random() * 2 * math.pi, rng.random() * 2 * math.pi
                w2 = complex(r2 * math.cos(p2), r2 * math.sin(p2)); w1 = zc * w2 + complex(r1 * math.cos(p1), r1 * math.sin(p1))
                T1 = K(round(w1.real)) + K(round(w1.imag)) * I; T2 = K(round(w2.real)) + K(round(w2.imag)) * I
                if T1 == 0 and T2 == 0: continue
                o1, o2 = a * T1 + b * T2, c * T1 + d * T2
                u, v, z = [kq(co[0]) * o1 ** 2 + kq(co[1]) * o1 * o2 + kq(co[2]) * o2 ** 2 for co in cl['X']]
                if z == 0: continue
                x = es[0] + d1 * (u / z) ** 2
                hs.append(float(x.global_height()))
            hs.sort()
            say('  %-3s %-10s boundary h_p(t) = H: h(x) min %.1f, median %.1f, max %.1f  (%d samples)' % (
                cv, cl['label'], hs[0], hs[len(hs) // 2], hs[-1], len(hs)))
    if WRITE:
        open(os.path.join(HERE, 'reach_H%d.txt' % int(H)), 'w').write('\n'.join(OUT) + '\n')

if __name__ == '__main__':
    main()
