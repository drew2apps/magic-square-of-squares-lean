#!/usr/bin/env python3
"""altmodels.py (conductor-t13): the other minimal models of the 2-covering quartics of covers.json at the primes of
multiplicative reduction, so that points on the non-identity components there (whose naive height on the first model
is inflated by up to a power of N(p)) become small on some model.  EXPERIMENT; written by an AI model (build agent
t13), not read by a person.

    python3 altmodels.py [--write]   -> altmodels.txt, altmodels.json (same keys as covers.json)

At a prime P of S, a move G -> pi^-2 G(pi t1 + r t2, t2) (r a root of G mod P of multiplicity >= 2, or of G/pi mod P
if pi exactly divides the content; also the move at infinity) that leaves an integral quartic keeps I and J: it gives
another minimal model of the same 2-covering, equivalent over K but in general not over O_P.  For each covering the
models reachable by at most DEPTH such moves at each chosen prime are listed (breadth-first, the start model
included), combined over the chosen primes, reduced as in covers.py, and de-duplicated (equal coefficient lists up to
sign).  Each model keeps the change of variable to the original parameter, so search.py and certify.py apply unchanged.
The control: on the torsion classes, torsion points whose preimage on the first model has h_p ~ 1e10 must appear with
small h_p on some model.
"""
from sage.all__sagemath_schemes import *
import sys, os, json, random, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', 'scripts', 'second_descent'))
import kctp_recheck as R
import covers as CV

K = R.K; I = R.I
WRITE = '--write' in sys.argv
OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(str(s))
def kq(p): return K(QQ(p[0])) + K(QQ(p[1])) * I

DEPTH = 2

def local_models(G, M, P):
    """models reachable from (G, M) by at most DEPTH invariant-preserving moves at P (the start included)."""
    pi = P.gens_reduced()[0]
    out = [(G, M)]; frontier = [(G, M)]; seen = {tuple(CV.coeffs(G))}
    for dep in range(DEPTH):
        new = []
        for H, MH in frontier:
            for r in CV.multiple_roots(H, P):
                H2, Mv = CV.move(H, P, r)
                if CV.cval(H2, P) >= 2:
                    H3 = H2 / pi ** 2
                    key = tuple(CV.coeffs(H3))
                    if key in seen: continue
                    seen.add(key); new.append((H3, CV.mmul(MH, Mv)))
        out += new; frontier = new
    return out

def main():
    cov = json.load(open(os.path.join(HERE, 'covers.json')))['curves']
    say('altmodels.py (conductor-t13).  Written by an AI model (build agent t13); not read by a person.')
    res = {'curves': {}}
    for cv in ('E', 'E3'):
        es = [kq(x) for x in cov[cv]['e']]
        D = 2 * (es[0] - es[1]) * (es[0] - es[2]) * (es[1] - es[2])
        primes = sorted(set(P for P, _ in K.ideal(D).factor()), key=lambda P: (-P.norm(), str(P.gens_reduced()[0])))
        odd = [P for P in primes if P.norm() % 2 == 1]       # the largest primes first (largest height inflation)
        out = []
        big = [P for P in odd if P.norm() >= 10 ** 4]
        for cl in cov[cv]['classes']:
            G = CV.from_coeffs([kq(x) for x in cl['G']])
            M = [[kq(x) for x in row] for row in cl['M']]
            # the base model; the other models at one large prime; and at the two largest primes together
            loc = {}
            for P in big:
                loc[P] = local_models(G, M, P)[1:]
            models = [(G, M)]
            for P in big: models += loc[P]
            if len(big) >= 2:
                for H, MH in loc[big[0]]:
                    models += local_models(H, MH, big[1])[1:]
            sizes = {str(P.gens_reduced()[0]): len(loc[P]) for P in big}
            # reduce and de-duplicate
            red = []; keys = set()
            for H, MH in models:
                H2, M2 = CV.reduce_quartic(H, MH)
                key = tuple(CV.coeffs(H2)); keyn = tuple(-c for c in key)
                if key in keys or keyn in keys: continue
                keys.add(key); red.append((H2, M2))
            say('  %s %s: models per prime %s -> %d models after reduction' % (cv, cl['label'], sizes, len(red)))
            for j, (H2, M2) in enumerate(red):
                zc, tc = CV.covariant(H2)
                ent = dict(cl); ent['label'] = '%s m%02d' % (cl['label'], j)
                ent['G'] = [CV.kstr(x) for x in CV.coeffs(H2)]
                ent['M'] = [[CV.kstr(x) for x in row] for row in M2]
                ent['covariant'] = [float('%.6f' % float(zc.real)), float('%.6f' % float(zc.imag)), float('%.6e' % float(tc))]
                out.append(ent)
        res['curves'][cv] = {'e': cov[cv]['e'], 'classes': out}
    if WRITE:
        json.dump(res, open(os.path.join(HERE, 'altmodels.json'), 'w'), indent=1, sort_keys=True)
        open(os.path.join(HERE, 'altmodels.txt'), 'w').write('\n'.join(OUT) + '\n')

if __name__ == '__main__':
    main()
